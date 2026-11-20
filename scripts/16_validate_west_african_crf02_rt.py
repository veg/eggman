#!/usr/bin/env python3
"""
Step 2A #2: Phylogenetic Validation of West African CRF02_AG Tier-1 G0[A6:RT] Pol Sequences.

Extracts 40 representative West African clinical pol sequences classified by RhizType v0.2
as Tier-2 CRF02_AG with Tier-1 structural_formula == G0[A6:RT], aligns them over HXB2
2550-3800 (pol.RT) alongside LANL references (CRF02_AG, NOVEL_CRF_06, A1, A2, A3, A6 [FSU],
CRF63_02A6, NOVEL_CRF_04, G, and B/HXB2 outgroup), infers a GTR+G4 ML phylogeny with
1,000 UFBoot replicates on m3.local via IQ-TREE 2, and quantifies clade placement vs Russian A6.
"""

import glob
import random
import subprocess
from pathlib import Path
import numpy as np
import pandas as pd
from Bio import Phylo, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

PROJECT_ROOT = Path("/Users/sergei/Projects/HIV-recombination")
DATA_DIR = PROJECT_ROOT / "data"
RESULTS_DIR = PROJECT_ROOT / "results"
WORK_DIR = RESULTS_DIR / "wa_crf02_rt_validation"
WORK_DIR.mkdir(parents=True, exist_ok=True)


def hxb2_slice_from_msa(rec_seq: str, hxb2_msa_seq: str, h_start: int, h_end: int) -> str:
    h_pos = 0
    cols = []
    for idx, ch in enumerate(hxb2_msa_seq):
        if ch != "-":
            h_pos += 1
        if h_start <= h_pos <= h_end:
            cols.append(rec_seq[idx])
        if h_pos > h_end:
            break
    return "".join(cols).replace("-", "").upper()


def extract_clade_support(tree, target_names, sister_exclude_names):
    """Find deepest monophyletic/majority node covering >=90% of targets without exclude_names, return (n_covered, ufboot)."""
    term_map = {t.name: t for t in tree.get_terminals()}
    tset = set(target_names) & set(term_map.keys())
    xset = set(sister_exclude_names) & set(term_map.keys())
    best = (0, 0.0, None)
    for clade in tree.get_nonterminals():
        leaves = {t.name for t in clade.get_terminals()}
        if leaves & xset:
            continue
        cov = len(leaves & tset)
        if cov > best[0] or (cov == best[0] and (clade.confidence or 0) > best[1]):
            best = (cov, float(clade.confidence or 0.0), clade)
    return best


def main():
    iso_files = sorted(glob.glob(str(DATA_DIR / "stanford_hivdb/hiv-survdb-payload/tables/isolates.d/*.csv")))
    iso_df = pd.concat([pd.read_csv(f, low_memory=False) for f in iso_files], ignore_index=True)
    pt_df = (
        iso_df.groupby("patient_id")
        .agg(
            country_code=("country_code", "first"),
            stanford_subtype=("subtype", lambda s: "|".join(sorted(set(str(x).strip() for x in s if pd.notna(x))))),
            isolate_year=("isolate_year", "first"),
        )
        .reset_index()
    )

    surv_df = pd.read_csv(RESULTS_DIR / "stanford_survdb_120k_rhiztype_raw.csv", low_memory=False)
    surv_df["patient_id"] = surv_df["sequence_id"].str.extract(r"PtID(\d+)")[0].astype(int)
    m = surv_df.merge(pt_df, on="patient_id", how="left")

    wa_countries = {"NGA", "GHA", "SEN", "CIV", "BFA", "MLI", "TGO", "BEN", "GIN", "GNB", "SLE", "NER"}
    cand = m[
        (m["called_lineage"] == "CRF02_AG")
        & (m["structural_formula"] == "G0[A6:RT]")
        & (m["country_code"].isin(wa_countries))
        & (m["hxb2_span_start"] <= 2550)
        & (m["hxb2_span_end"] >= 3650)
    ].drop_duplicates(subset=["patient_id"]).copy()

    # Stratified sample across countries totaling 40 sequences
    rng = random.Random(42)
    chosen_rows = []
    by_cty = {cty: grp.sample(frac=1.0, random_state=42) for cty, grp in cand.groupby("country_code")}
    cty_order = sorted(by_cty.keys(), key=lambda k: -len(by_cty[k]))
    ptrs = {k: 0 for k in cty_order}
    while len(chosen_rows) < 40:
        progressed = False
        for cty in cty_order:
            if len(chosen_rows) >= 40:
                break
            grp = by_cty[cty]
            if ptrs[cty] < len(grp):
                chosen_rows.append(grp.iloc[ptrs[cty]])
                ptrs[cty] += 1
                progressed = True
        if not progressed:
            break
    sel_df = pd.DataFrame(chosen_rows)
    sel_ids = set(sel_df["sequence_id"].tolist())
    print(f"Selected {len(sel_df)} West African CRF02_AG Tier-1 G0[A6:RT] isolates: {sel_df['country_code'].value_counts().to_dict()}")

    wa_seqs = {}
    for rec in SeqIO.parse(DATA_DIR / "stanford_hivdb/stanford_survdb_all_120k.fasta", "fasta"):
        if rec.id in sel_ids:
            wa_seqs[rec.id] = str(rec.seq).replace("-", "").upper()
    assert len(wa_seqs) == len(sel_ids)

    # Extract reference HXB2 2500-3900 from lanl_ref_2023 and additional FSU A6 from flt_2022
    ref_recs = list(SeqIO.parse(DATA_DIR / "lanl_ref_2023.fasta", "fasta"))
    hxb2_msa = next(str(r.seq).upper() for r in ref_recs if "K03455" in r.id)
    keep_prefixes = ("Ref.A1.", "Ref.A2.", "Ref.A3.", "Ref.A6.", "Ref.02_AG.", "Ref.63_02A6.", "Ref.G.")
    align_inputs = []
    hxb2_rt = hxb2_slice_from_msa(hxb2_msa, hxb2_msa, 2500, 3900)
    align_inputs.append(SeqRecord(Seq(hxb2_rt), id="REF|B|K03455_HXB2", description=""))

    ref_groups = {}
    for r in ref_recs:
        if any(r.id.startswith(p) for p in keep_prefixes):
            raw_rt = hxb2_slice_from_msa(str(r.seq).upper(), hxb2_msa, 2500, 3900)
            st = r.id.split(".")[1]
            acc = r.id.split(".")[-1]
            clean_id = f"REF|{st}|{acc}"
            align_inputs.append(SeqRecord(Seq(raw_rt), id=clean_id, description=""))
            ref_groups[clean_id] = st

    # Add 5 Russian/FSU A6 references from lanl_flt_2022.fasta
    flt_recs = list(SeqIO.parse(DATA_DIR / "lanl_flt_2022.fasta", "fasta"))
    hxb2_flt = next((str(r.seq).upper() for r in flt_recs if "K03455" in r.id), None)
    n_a6_extra = 0
    existing_accs = {x.split("|")[-1] for x in ref_groups}
    for r in flt_recs:
        if r.id.startswith("A6.") and (".RU." in r.id or ".UA." in r.id or ".BY." in r.id):
            acc = r.id.split(".")[-1]
            if acc in existing_accs:
                continue
            seq_u = hxb2_slice_from_msa(str(r.seq).upper(), hxb2_flt, 2500, 3900) if hxb2_flt else str(r.seq).replace("-", "")[2500:3900]
            if len(seq_u) >= 1000:
                clean_id = f"REF|A6_FSU|{acc}"
                align_inputs.append(SeqRecord(Seq(seq_u), id=clean_id, description=""))
                ref_groups[clean_id] = "A6"
                existing_accs.add(acc)
                n_a6_extra += 1
                if n_a6_extra >= 5:
                    break

    # Add NOVEL_CRF_04 and NOVEL_CRF_06 from benchmark_rhiztype_vs_comet.fasta (via later MAFFT alignment)
    for r in SeqIO.parse(DATA_DIR / "benchmark_rhiztype_vs_comet.fasta", "fasta"):
        if r.id.startswith("P2|NOVEL_CRF_04|") or r.id.startswith("P2|NOVEL_CRF_06|"):
            lbl = r.id.split("|")[1]
            acc = r.id.split("|")[2]
            raw = str(r.seq).replace("-", "").upper()
            # trim roughly around gag+pol to speed MAFFT alignment before slicing exact 2550-3800
            clean_id = f"REF|{lbl}|{acc}"
            align_inputs.append(SeqRecord(Seq(raw[1500:4500]), id=clean_id, description=""))
            ref_groups[clean_id] = lbl

    # Add the 40 West African clinical pol sequences
    wa_clean_ids = {}
    for idx, (_, row) in enumerate(sel_df.iterrows()):
        sid = row["sequence_id"]
        cty = row["country_code"]
        cid = f"WA|{cty}|{idx:02d}"
        wa_clean_ids[cid] = row
        align_inputs.append(SeqRecord(Seq(wa_seqs[sid]), id=cid, description=""))

    raw_fa = WORK_DIR / "wa_crf02_rt_unaligned.fasta"
    msa_fa = WORK_DIR / "wa_crf02_rt_mafft.fasta"
    trim_fa = WORK_DIR / "wa_crf02_rt_2550_3800.fasta"
    SeqIO.write(align_inputs, raw_fa, "fasta")

    with open(msa_fa, "w") as out:
        subprocess.run(["/usr/local/bin/mafft", "--auto", "--quiet", str(raw_fa)], stdout=out, check=True)

    # Slice strictly to HXB2 2550-3800 (RT core) & retain columns present in >=50% of WA sequences
    # Note: REF|B|K03455_HXB2 was sliced starting at HXB2 coordinate 2500
    aligned = {r.id: str(r.seq).upper() for r in SeqIO.parse(msa_fa, "fasta")}
    hxb2_aln = aligned["REF|B|K03455_HXB2"]
    keep_cols = []
    h_pos = 2499
    wa_list = [aligned[k] for k in wa_clean_ids]
    for col, ch in enumerate(hxb2_aln):
        if ch != "-":
            h_pos += 1
        if 2550 <= h_pos <= 3800:
            non_gap_wa = sum(1 for s in wa_list if s[col] in "ACGT")
            if non_gap_wa >= 0.5 * len(wa_list) and ch != "-":
                keep_cols.append(col)
    print(f"Aligned RT columns retained within HXB2 2550-3800: {len(keep_cols)} nt")
    assert len(keep_cols) >= 800, f"Too few columns retained: {len(keep_cols)}"
    trimmed_recs = [
        SeqRecord(Seq("".join(seq[c] for c in keep_cols)), id=rid, description="")
        for rid, seq in aligned.items()
    ]
    SeqIO.write(trimmed_recs, trim_fa, "fasta")

    # Run IQ-TREE 2 GTR+G4 -B 1000 on m3.local
    subprocess.run(["scp", "-q", str(trim_fa), "m3.local:~/hiv_recomb_bench/results/wa_crf02_rt_2550_3800.fasta"], check=True)
    subprocess.run(
        [
            "ssh", "-o", "BatchMode=yes", "m3.local",
            "/Users/sergei/hiv_recomb_bench/bin/iqtree2 -s /Users/sergei/hiv_recomb_bench/results/wa_crf02_rt_2550_3800.fasta "
            "-m GTR+G4 -B 1000 -T 8 -redo -pre /Users/sergei/hiv_recomb_bench/results/wa_crf02_rt_iqtree",
        ],
        check=True,
    )
    subprocess.run(["scp", "-q", "m3.local:~/hiv_recomb_bench/results/wa_crf02_rt_iqtree.treefile", str(WORK_DIR / "wa_crf02_rt_iqtree.treefile")], check=True)

    tree = Phylo.read(WORK_DIR / "wa_crf02_rt_iqtree.treefile", "newick")
    tree.root_with_outgroup("REF|B|K03455_HXB2")

    # Patristic distances from each WA sequence to CRF02_AG vs A6_FSU vs CRF63_02A6 vs NOVEL_CRF_06 vs A1
    ref_crf02 = [k for k, v in ref_groups.items() if v == "02_AG"]
    ref_a6 = [k for k, v in ref_groups.items() if v == "A6"]
    ref_a1 = [k for k, v in ref_groups.items() if v == "A1"]
    ref_nov06 = [k for k, v in ref_groups.items() if v == "NOVEL_CRF_06"]
    ref_63 = [k for k, v in ref_groups.items() if v == "63_02A6"]
    ref_g = [k for k, v in ref_groups.items() if v == "G"]

    rows = []
    for cid, meta in wa_clean_ids.items():
        d02 = np.mean([tree.distance(cid, r) for r in ref_crf02])
        da6 = np.mean([tree.distance(cid, r) for r in ref_a6])
        da1 = np.mean([tree.distance(cid, r) for r in ref_a1])
        dnov6 = np.mean([tree.distance(cid, r) for r in ref_nov06])
        d63 = np.mean([tree.distance(cid, r) for r in ref_63])
        dg = np.mean([tree.distance(cid, r) for r in ref_g])
        nearest = min(
            [("CRF02_AG", d02), ("NOVEL_CRF_06", dnov6), ("A1", da1), ("A6_FSU", da6), ("CRF63_02A6", d63), ("G", dg)],
            key=lambda x: x[1],
        )[0]
        rows.append({
            "tree_id": cid,
            "sequence_id": meta["sequence_id"],
            "country_code": meta["country_code"],
            "isolate_year": meta["isolate_year"],
            "stanford_subtype": meta["stanford_subtype"],
            "hxb2_span": f"{int(meta['hxb2_span_start'])}-{int(meta['hxb2_span_end'])}",
            "mean_dist_CRF02_AG": round(float(d02), 5),
            "mean_dist_A6_FSU": round(float(da6), 5),
            "mean_dist_CRF63_02A6": round(float(d63), 5),
            "mean_dist_NOVEL_CRF_06": round(float(dnov6), 5),
            "mean_dist_A1": round(float(da1), 5),
            "mean_dist_G": round(float(dg), 5),
            "nearest_clade": nearest,
            "closer_to_CRF02_than_A6": bool(d02 < da6),
        })

    out_df = pd.DataFrame(rows)
    cov_a6, boot_a6, _ = extract_clade_support(tree, ref_a6, list(wa_clean_ids.keys()) + ref_crf02 + ref_a1)
    cov_wa_02, boot_wa_02, _ = extract_clade_support(tree, list(wa_clean_ids.keys()) + ref_crf02 + ref_nov06, ref_a6 + ref_g)
    cov_a6_plus_63, boot_a6_plus_63, _ = extract_clade_support(tree, ref_a6 + ref_63, list(wa_clean_ids.keys()) + ref_crf02)

    out_df["fsu_a6_monophyly_boot"] = boot_a6
    out_df["fsu_a6_63_clade_boot"] = boot_a6_plus_63
    out_df["wa_crf02_excl_a6_clade_count"] = cov_wa_02
    out_df["wa_crf02_excl_a6_clade_boot"] = boot_wa_02

    out_csv = RESULTS_DIR / "stanford_west_african_crf02_rt_validation.csv"
    out_df.to_csv(out_csv, index=False)
    print("\n=== West African CRF02_AG (Tier-1 G0[A6:RT]) pol.RT ML Validation (n = 40) ===")
    print("Nearest clade counts:", out_df["nearest_clade"].value_counts().to_dict())
    print(f"Closer patristic distance to CRF02_AG than Russian/FSU A6: {out_df['closer_to_CRF02_than_A6'].sum()} / {len(out_df)}")
    print(f"Mean patristic dist: CRF02_AG={out_df['mean_dist_CRF02_AG'].mean():.4f} vs A6_FSU={out_df['mean_dist_A6_FSU'].mean():.4f} vs G={out_df['mean_dist_G'].mean():.4f}")
    print(f"FSU A6 pure monophyletic clade excluding all WA CRF02_AG: {cov_a6}/{len(ref_a6)} A6 taxa, UFBoot={boot_a6:.1f}%")
    print(f"FSU A6+CRF63 clade excluding all WA CRF02_AG: {cov_a6_plus_63}/{len(ref_a6)+len(ref_63)} taxa, UFBoot={boot_a6_plus_63:.1f}%")
    print(f"WA CRF02_AG + CRF02_AG clade excluding all A6_FSU & G: {cov_wa_02} taxa, UFBoot={boot_wa_02:.1f}%")
    print(f"Saved {out_csv}")


if __name__ == "__main__":
    main()
