#!/usr/bin/env python3
"""
Script 14: Head-to-Head Benchmark of EGGMAN vs. COMET v2.4 (Struck et al., NAR 2014).

Evaluates both tools on identical unaligned FASTA sequences across four stratified panels:
  1. LANL 2023 Full-Length Reference Compendium (N = 517 genomes: 58 pure-subtype refs + 459 CRF refs)
  2. Audited Member Genomes of the 15 Legacy Full-Genome Candidate Clusters (N = 61: 28 post-2023 official
     LANL CRFs + 3 CRF63_02A6 + 6 RU_G0B_gagMA + 24 polyphyletic convergent URFs)
  3. Representative Isolates of Partial pol & env Lineages (N = 112 sequences)
  4. Stratified Clinical Patient Validation Cohort from LANL 2022 Web Alignment (N = 400 non-reference genomes)
"""

import argparse
import gzip
import json
import os
import re
import subprocess
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path
import numpy as np
import pandas as pd

PROJECT_ROOT = Path("/Users/sergei/Projects/HIV-recombination")
DATA_DIR = PROJECT_ROOT / "data"
RESULTS_DIR = PROJECT_ROOT / "results"
RHIZTYPE_BIN = PROJECT_ROOT / "rhiztype" / "target" / "release" / "eggman"

# Canonical accession-level disposition of the 61 legacy Panel 2 sequences (NOVEL_CRF_01..15)
LEGACY_CRF04_SCRUB = {
    "JX500700": "CRF63_02A6",
    "MK984159": "CRF63_02A6",
    "PZ426393": "CRF63_02A6",
    "PQ585410": "CRF173_63A6",
    "PQ585412": "CRF173_63A6",
}


def read_fasta(path: Path):
    records = []
    header = None
    chunks = []
    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(chunks).upper()))
                header = line[1:].strip()
                chunks = []
            else:
                chunks.append(line)
        if header is not None:
            records.append((header, "".join(chunks).upper()))
    return records


def degap(seq: str) -> str:
    return re.sub(r"[^ACGTUacgtu]", "", seq).upper()


def submit_batch_to_comet(fasta_records, batch_idx: int) -> pd.DataFrame:
    """Submit a single batch of <=800 FASTA records to COMET v2.4 and return DataFrame."""
    tmp_fa = Path(f"/tmp/comet_batch_{batch_idx}.fasta")
    gz_path = Path(f"/tmp/comet_batch_{batch_idx}.fasta.gz")
    with open(tmp_fa, "w") as f:
        for sid, s_raw in fasta_records:
            f.write(f">{sid}\n{s_raw}\n")
    with open(tmp_fa, "rb") as f_in, gzip.open(gz_path, "wb") as f_out:
        f_out.writelines(f_in)

    cookie_file = f"/tmp/comet_bench_cookies_{batch_idx}.txt"
    if os.path.exists(cookie_file):
        os.unlink(cookie_file)

    html = subprocess.check_output(
        ["curl", "-s", "-c", cookie_file, "-b", cookie_file, "https://comet.lih.lu/index.php?cat=hiv1"],
        text=True,
    )
    m_csrf = re.search(r'name="csrf" value="([0-9a-f]+)"', html)
    if not m_csrf:
        raise RuntimeError("Could not parse CSRF token from comet.lih.lu")
    csrf = m_csrf.group(1)

    t0 = time.time()
    print(f"[*] Uploading batch {batch_idx} ({len(fasta_records)} seqs, {gz_path.stat().st_size / 1024:.1f} KB) to COMET v2.4...")
    post_headers = subprocess.check_output(
        [
            "curl", "-s", "-i", "-c", cookie_file, "-b", cookie_file,
            "-F", f"fastafile=@{gz_path}",
            "-F", "non_commercial=confirmed",
            "-F", f"csrf={csrf}",
            "-F", "submit=submit fasta file",
            "https://comet.lih.lu/index.php?cat=hiv1",
        ],
        text=True,
    )
    m_job = re.search(r"Location:\s*index\.php\?job=([0-9a-f]+)&cat=hiv1", post_headers)
    if not m_job:
        raise RuntimeError(f"Could not find job ID in COMET response:\n{post_headers[:800]}")
    job_id = m_job.group(1)
    poll_url = f"https://comet.lih.lu/index.php?job={job_id}&cat=hiv1"
    csv_url = f"https://comet.lih.lu/csv.php?job={job_id}"

    for attempt in range(90):
        time.sleep(2.0)
        page = subprocess.check_output(
            ["curl", "-s", "-c", cookie_file, "-b", cookie_file, poll_url],
            text=True,
        )
        if f"/csv.php?job={job_id}" in page:
            elapsed = time.time() - t0
            print(f"[+] COMET v2.4 batch {batch_idx} finished in {elapsed:.2f}s!")
            tsv_text = subprocess.check_output(
                ["curl", "-s", "-c", cookie_file, "-b", cookie_file, csv_url],
                text=True,
            )
            tmp_tsv = Path(f"/tmp/comet_batch_{batch_idx}.tsv")
            with open(tmp_tsv, "w") as f:
                f.write(tsv_text)
            return pd.read_csv(tmp_tsv, sep="\t")

    raise RuntimeError("COMET job timed out after 180 seconds")


def submit_to_comet(fasta_entries, out_tsv: Path):
    """Ensure all fasta_entries have COMET v2.4 predictions (accounting for COMET's 1,000-sequence upload cap)."""
    existing_df = pd.DataFrame()
    if out_tsv.exists() and out_tsv.stat().st_size > 1000:
        existing_df = pd.read_csv(out_tsv, sep="\t")
        done_ids = set(existing_df["name"].astype(str))
    else:
        done_ids = set()

    missing = [(sid, s_raw) for sid, s_raw in fasta_entries if sid not in done_ids]
    if not missing:
        print(f"[+] Using complete cached COMET v2.4 results ({len(existing_df):,} sequences) at {out_tsv}")
        return existing_df

    print(f"[*] Submitting {len(missing)} remaining sequences to COMET v2.4 (COMET web server enforces a 1,000-seq cap per upload)...")
    batch_dfs = [existing_df] if not existing_df.empty else []
    batch_size = 800
    for b_idx, i in enumerate(range(0, len(missing), batch_size), start=1):
        chunk = missing[i:i + batch_size]
        bdf = submit_batch_to_comet(chunk, b_idx)
        batch_dfs.append(bdf)

    full_df = pd.concat(batch_dfs, ignore_index=True).drop_duplicates(subset=["name"])
    full_df.to_csv(out_tsv, sep="\t", index=False)
    return full_df


def load_panel2_audit_map():
    """Return {acc: (audit_id, disposition, canonical_lineage, formula)} from 31k audit TSV."""
    aud_path = RESULTS_DIR / "novel_candidate_crfs_31k_audit.tsv"
    cand_df = pd.read_csv(RESULTS_DIR / "genbank_novel_circulating_recombinants.csv")
    cand_forms = dict(zip(cand_df.get("canonical_id", pd.Series(dtype=str)), cand_df.get("rhizaeon_formula", pd.Series(dtype=str))))
    acc_map = {}
    if aud_path.exists():
        aud_df = pd.read_csv(aud_path, sep="\t")
        for _, r in aud_df.iterrows():
            aid = str(r["audit_id"])
            disp = str(r["disposition"])
            canon = str(r["canonical_lineage"])
            accs = [a.strip() for a in str(r["legacy_accessions"]).split(",") if a.strip()]
            for a in accs:
                cl = LEGACY_CRF04_SCRUB.get(a, canon) if aid == "NOVEL_CRF_04" else canon
                d_eff = "Promoted_Official_LANL_CRF" if aid == "NOVEL_CRF_04" else disp
                form = cand_forms.get(aid, cand_forms.get(cl, ""))
                acc_map[a] = (aid, d_eff, cl, form)
    return acc_map


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--rebuild-fasta", action="store_true", help="Rebuild data/benchmark_rhiztype_vs_comet.fasta from source FASTAs (default: preserve existing de-overlapped 1,093-seq FASTA)")
    args = parser.parse_args()

    print("=" * 96)
    print("STEP 14: HEAD-TO-HEAD BENCHMARK OF EGGMAN vs. COMET v2.4 (STRUCK ET AL. 2014)")
    print("=" * 96)

    bench_fasta = DATA_DIR / "benchmark_rhiztype_vs_comet.fasta"
    p2_audit = load_panel2_audit_map()
    pure_set = {"A1", "A2", "A3", "A4", "A6", "A7", "A8", "B", "C", "D", "F1", "F2", "G", "H", "J", "K", "L"}

    if bench_fasta.exists() and not args.rebuild_fasta:
        print(f"[+] Preserving existing de-overlapped benchmark FASTA at {bench_fasta}")
        fasta_entries = read_fasta(bench_fasta)
        meta_rows = []
        for sid, raw in fasta_entries:
            parts = sid.split("|")
            if parts[0] == "P1":
                hdr = parts[1]
                cl = hdr.split(".")[1]
                if cl in pure_set:
                    sub_cat = "Panel1A_LANL_Ref_PureSubtype"
                else:
                    crf_num = int(cl.split("_")[0]) if cl.split("_")[0].isdigit() else 999
                    sub_cat = (
                        "Panel1B_LANL_Ref_CRF_01_to_118 (In COMET v2.4 Training)"
                        if crf_num <= 118
                        else "Panel1C_LANL_Ref_CRF_119_to_157 (Post-COMET v2.4 Cutoff)"
                    )
                meta_rows.append({"seq_id": sid, "panel": sub_cat, "ground_truth_label": cl, "original_id": hdr})
            elif parts[0] == "P2":
                lin_code, acc = parts[1], parts[2]
                _, disp, canon, form = p2_audit.get(acc, (lin_code, "Polyphyletic_Convergent_URF", "URF", ""))
                gt_lbl = f"{canon}:{form}" if form else canon
                meta_rows.append({
                    "seq_id": sid,
                    "panel": "Panel2_Novel_FullGenome_CRFs (Audited 61)",
                    "ground_truth_label": gt_lbl,
                    "original_id": acc,
                    "p2_disposition": disp,
                    "p2_canonical_lineage": canon,
                })
            elif parts[0] == "P3":
                lin_code, acc = parts[1], parts[2]
                p_cat = "Panel3A_Novel_Partial_pol_Lineages" if lin_code.startswith("POL_") else "Panel3B_Novel_Partial_env_Lineages"
                meta_rows.append({"seq_id": sid, "panel": p_cat, "ground_truth_label": lin_code, "original_id": acc})
            elif parts[0] == "P4":
                cat_key, hdr = parts[1], parts[2]
                p_map = {
                    "Pure": "Panel4A_Clinical_PureSubtype",
                    "CRF": "Panel4B_Clinical_NamedCRF",
                    "URF": "Panel4C_Clinical_URF",
                }
                meta_rows.append({
                    "seq_id": sid,
                    "panel": p_map.get(cat_key, "Panel4C_Clinical_URF"),
                    "ground_truth_label": hdr.split(".")[0],
                    "original_id": hdr,
                })
    else:
        meta_rows = []
        fasta_entries = []

        # 1. Panel 1: LANL 2023 Reference Compendium (N=517)
        ref_recs = read_fasta(DATA_DIR / "lanl_ref_2023.fasta")
        for hdr, seq in ref_recs:
            if hdr.startswith("Ref.HXB2"):
                continue
            raw = degap(seq)
            cl = hdr.split(".")[1]
            if cl in pure_set:
                sub_cat = "Panel1A_LANL_Ref_PureSubtype"
            else:
                crf_num = int(cl.split("_")[0]) if cl.split("_")[0].isdigit() else 999
                if crf_num <= 118:
                    sub_cat = "Panel1B_LANL_Ref_CRF_01_to_118 (In COMET v2.4 Training)"
                else:
                    sub_cat = "Panel1C_LANL_Ref_CRF_119_to_157 (Post-COMET v2.4 Cutoff)"
            seq_id = f"P1|{hdr}"
            fasta_entries.append((seq_id, raw))
            meta_rows.append({
                "seq_id": seq_id,
                "panel": sub_cat,
                "ground_truth_label": cl,
                "original_id": hdr,
            })

        # 2. Panel 2: Member Genomes of the 15 Legacy Candidate Clusters (Audited)
        nov_full_df = pd.read_csv(RESULTS_DIR / "genbank_novel_circulating_recombinants.csv")
        nov15_df = nov_full_df[nov_full_df["lineage_status"].str.startswith("Homologous_Circulating_Lineage")].copy()
        acc_to_nov_info = {}
        for idx, r in nov15_df.reset_index(drop=True).iterrows():
            lin_code = f"NOVEL_CRF_{idx+1:02d}"
            form = str(r["rhizaeon_formula"])
            for a in str(r["member_accessions"]).split(","):
                a = a.strip()
                if a:
                    acc_to_nov_info[a] = (lin_code, form)

        gb_31k_path = Path("/Users/sergei/Projects/TOGA_MEME/recombination/pan_group_m_nomenclature/data/ncbi_hiv1_all_complete_31k.fasta")
        if gb_31k_path.exists():
            cur_acc = None
            cur_chunks = []
            with open(gb_31k_path, "r") as f:
                for line in f:
                    if line.startswith(">"):
                        if cur_acc in acc_to_nov_info:
                            lin_code, form = acc_to_nov_info[cur_acc]
                            seq_id = f"P2|{lin_code}|{cur_acc}"
                            fasta_entries.append((seq_id, degap("".join(cur_chunks))))
                            meta_rows.append({
                                "seq_id": seq_id,
                                "panel": "Panel2_Novel_FullGenome_CRFs (Audited 61)",
                                "ground_truth_label": f"{lin_code}:{form}",
                                "original_id": cur_acc,
                            })
                        cur_acc = line[1:].strip().split()[0].split(".")[0]
                        cur_chunks = []
                    else:
                        cur_chunks.append(line.strip())
                if cur_acc in acc_to_nov_info:
                    lin_code, form = acc_to_nov_info[cur_acc]
                    seq_id = f"P2|{lin_code}|{cur_acc}"
                    fasta_entries.append((seq_id, degap("".join(cur_chunks))))
                    meta_rows.append({
                        "seq_id": seq_id,
                        "panel": "Panel2_Novel_FullGenome_CRFs (Audited 61)",
                        "ground_truth_label": f"{lin_code}:{form}",
                        "original_id": cur_acc,
                    })

        # 3. Panel 3: Representative Isolates of Partial pol & env Lineages
        nov_pol_df = pd.read_csv(RESULTS_DIR / "scale_partial_pol_novel_circulating_lineages.csv")
        nov_env_df = pd.read_csv(RESULTS_DIR / "scale_partial_env_novel_circulating_lineages.csv")

        pol_target_accs = {}
        for idx, r in nov_pol_df.reset_index(drop=True).iterrows():
            l_id = f"POL_NOV_{idx+1:02d}"
            form = str(r["window_formula"])
            for a in str(r["example_accessions"]).split(","):
                a = a.strip()
                if a:
                    pol_target_accs[a] = (l_id, form)

        env_target_accs = {}
        for idx, r in nov_env_df.reset_index(drop=True).iterrows():
            l_id = "ENV_ICEBERG" if idx == 0 else f"ENV_NOV_{idx:02d}"
            form = str(r["window_formula"])
            for a in str(r["example_accessions"]).split(","):
                a = a.strip()
                if a:
                    env_target_accs[a] = (l_id, form)

        for fasta_file, target_map, panel_label in [
            (DATA_DIR / "genbank_partial_pol_dedup.fasta", pol_target_accs, "Panel3A_Novel_Partial_pol_Lineages"),
            (DATA_DIR / "lanl_all_2022_pol.fasta", pol_target_accs, "Panel3A_Novel_Partial_pol_Lineages"),
            (DATA_DIR / "genbank_partial_env_dedup.fasta", env_target_accs, "Panel3B_Novel_Partial_env_Lineages"),
            (DATA_DIR / "lanl_all_2022_env.fasta", env_target_accs, "Panel3B_Novel_Partial_env_Lineages"),
        ]:
            if not fasta_file.exists():
                continue
            seen_in_panel = set(m["original_id"] for m in meta_rows)
            for hdr, seq in read_fasta(fasta_file):
                acc = hdr.split(".")[0].split()[0]
                if "." in hdr and len(hdr.split(".")) >= 5:
                    acc = hdr.split(".")[-1].split()[0]
                if acc in target_map and acc not in seen_in_panel:
                    l_id, form = target_map[acc]
                    seq_id = f"P3|{l_id}|{acc}"
                    fasta_entries.append((seq_id, degap(seq)))
                    meta_rows.append({
                        "seq_id": seq_id,
                        "panel": panel_label,
                        "ground_truth_label": f"{l_id}:{form}",
                        "original_id": acc,
                    })
                    seen_in_panel.add(acc)

        # 4. Panel 4: Stratified Clinical Patient Cohort from LANL 2022 Filtered Web Alignment
        flt_recs = read_fasta(DATA_DIR / "lanl_flt_2022.fasta")
        rng = np.random.default_rng(20261027)
        by_cat = defaultdict(list)
        for hdr, seq in flt_recs:
            tag = hdr.split(".")[0]
            if any(tag.startswith(x) for x in ("O", "N", "P", "CPZ", "GOR")):
                continue
            if tag in pure_set:
                by_cat["Pure"].append((hdr, seq, tag))
            elif "_" in tag and tag.split("_")[0].isdigit():
                by_cat["CRF"].append((hdr, seq, tag))
            else:
                by_cat["URF"].append((hdr, seq, tag))

        for cat_key, n_take, p_name in [
            ("Pure", 150, "Panel4A_Clinical_PureSubtype"),
            ("CRF", 150, "Panel4B_Clinical_NamedCRF"),
            ("URF", 100, "Panel4C_Clinical_URF"),
        ]:
            pool_list = by_cat[cat_key]
            idxs = rng.choice(len(pool_list), size=min(n_take, len(pool_list)), replace=False)
            for i in idxs:
                hdr, seq, tag = pool_list[i]
                seq_id = f"P4|{cat_key}|{hdr}"
                fasta_entries.append((seq_id, degap(seq)))
                meta_rows.append({
                    "seq_id": seq_id,
                    "panel": p_name,
                    "ground_truth_label": tag,
                    "original_id": hdr,
                })

        with open(bench_fasta, "w") as f:
            for sid, s_raw in fasta_entries:
                f.write(f">{sid}\n{s_raw}\n")

    meta_df = pd.DataFrame(meta_rows)
    print(f"[+] Benchmark FASTA ready: {len(meta_df):,} sequences ({bench_fasta.stat().st_size / (1024*1024):.2f} MB)")
    for p_lbl, cnt in meta_df["panel"].value_counts().items():
        print(f"    * {p_lbl:<56s}: {cnt:4d} sequences")

    # Run EGGMAN locally and measure wall-clock runtime
    rhiz_csv = RESULTS_DIR / "benchmark_rhiztype_output.csv"
    t0_rhiz = time.time()
    subprocess.run(
        [str(RHIZTYPE_BIN), str(bench_fasta), "--quiet", "--out", str(rhiz_csv)],
        check=True,
    )
    rhiz_sec = time.time() - t0_rhiz
    print(f"[+] EGGMAN native CLI screened {len(meta_df):,} unaligned sequences in {rhiz_sec:.3f}s ({len(meta_df)/rhiz_sec:.1f} seq/s)!")

    # Sync with COMET v2.4 TSV (reuses cached predictions when all seq_ids are present)
    comet_tsv = RESULTS_DIR / "benchmark_comet_v24_output.tsv"
    comet_df = submit_to_comet(fasta_entries, comet_tsv)
    rhiz_df = pd.read_csv(rhiz_csv)
    rhiz_df["rhiztype_is_clonal"] = rhiz_df["n_breakpoints"].fillna(0).astype(int) == 0

    # Merge all
    merged = meta_df.merge(
        comet_df.rename(columns={"name": "seq_id", "subtype": "comet_subtype", "bootstrap support": "comet_support"}),
        on="seq_id",
        how="left",
    ).merge(
        rhiz_df.rename(columns={
            "sequence_id": "seq_id",
            "structural_formula": "rhiztype_formula",
            "descriptive_name": "rhiztype_descriptive",
            "backbone_trunk": "rhiztype_backbone",
            "n_breakpoints": "rhiztype_n_breakpoints",
            "hxb2_barcode": "rhiztype_barcode",
            "compendium_matches": "rhiztype_top_match",
        }),
        on="seq_id",
        how="left",
    )

    merged_out = RESULTS_DIR / "benchmark_rhiztype_vs_comet_merged.csv"
    merged.to_csv(merged_out, index=False)
    print(f"[+] Saved full head-to-head comparison table ({len(merged):,} rows) to {merged_out}")


if __name__ == "__main__":
    main()
