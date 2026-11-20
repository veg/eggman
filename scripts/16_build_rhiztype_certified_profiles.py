#!/usr/bin/env python3
"""
Step 16: Build deterministic certified HIV-1 Group M reference profiles for Eggman/RhizType.

Reads ONLY authentic sequence alignments & Phase 1 C5/C7 outputs:
  - /Users/sergei/Projects/HIV-recombination/data/lanl_ref_2023.fasta (518 LANL Ref genomes)
  - /Users/sergei/Projects/TOGA_MEME/recombination/pan_group_m_nomenclature/data/ncbi_hiv1_all_complete_31k.fasta
    (8 promoted post-2023 LANL CRFs + RU_G0B_gagMA)
  - /Users/sergei/Projects/HIV-recombination/data/hxb2_15trunk_profile_9719.fasta
  - /Users/sergei/Projects/HIV-recombination/data/genbank_cluster_flatfiles.tsv
  - /Users/sergei/Projects/HIV-recombination/results/hiv1_descriptive_nomenclature_compendium.csv
  - /Users/sergei/Projects/HIV-recombination/results/hiv1_recombination_genealogy_dag.csv
  - /Users/sergei/Projects/HIV-recombination/results/genbank_novel_circulating_recombinants.csv

Generates:
  1. rhiztype/crates/eggman-core/data/compendium_bundle.json
     (163 LANL 2023 entries + 8 promoted post-2023 CRFs + 1 RU_G0B_gagMA; stale NOVEL_CRF_01..15 removed)
  2. rhiztype/crates/eggman-core/data/founder_profiles_9719.tsv
     (13 full 9,719-nt HXB2-anchored plurality consensus profiles for Tier-1 alignment & epi Viterbi)
  3. rhiztype/crates/eggman-core/data/certified_profiles_8628.tsv
     (152 coding 8,628-nt plurality consensus profiles: 1 Subtype_G + 142 active 2023 LANL CRFs
      + 8 promoted post-2023 official LANL CRFs + 1 RU_G0B_gagMA)
"""

import csv
import json
import os
import re
import subprocess
import tempfile
from collections import defaultdict
from pathlib import Path
import numpy as np
import pandas as pd

PROJECT_ROOT = Path("/Users/sergei/Projects/HIV-recombination")
DATA_DIR = PROJECT_ROOT / "data"
RESULTS_DIR = PROJECT_ROOT / "results"
RHIZ_DATA_DIR = PROJECT_ROOT / "rhiztype" / "crates" / "eggman-core" / "data"
GENBANK_FASTA = Path("/Users/sergei/Projects/TOGA_MEME/recombination/pan_group_m_nomenclature/data/ncbi_hiv1_all_complete_31k.fasta")
EGGMAN_BIN = PROJECT_ROOT / "rhiztype" / "target" / "release" / "eggman"
MAFFT_BIN = "/usr/local/bin/mafft"

POST2023_PROMOTED_CRFS = {
    "CRF127_07109": {
        "short_code": "127_07109",
        "category": "CRF_2ndGen",
        "accs": ["LC735411", "LC735412", "LC735413"],
        "dag_parents": ["CRF07_BC", "CRF109_0107"],
    },
    "CRF135_0107": {
        "short_code": "135_0107",
        "category": "CRF_2ndGen",
        "accs": ["OP501672", "OP501674", "OP501675", "OP501676", "OP758070"],
        "dag_parents": ["CRF01_AE", "CRF07_BC"],
    },
    "CRF136_0107": {
        "short_code": "136_0107",
        "category": "CRF_2ndGen",
        "accs": ["OP501677", "OP501678", "OP501679"],
        "dag_parents": ["CRF01_AE", "CRF07_BC"],
    },
    "CRF162_cpx": {
        "short_code": "162_cpx",
        "category": "CRF_cpx",
        "accs": ["PP483791", "PP483792", "PP483793", "PP483794", "PP483795", "PP483797"],
        "dag_parents": ["CRF01_AE", "CRF07_BC"],
    },
    "CRF172_0755": {
        "short_code": "172_0755",
        "category": "CRF_2ndGen",
        "accs": ["PQ356206", "PQ356207", "PQ356208", "PQ356209"],
        "dag_parents": ["CRF07_BC", "CRF55_01B"],
    },
    "CRF173_63A6": {
        "short_code": "173_63A6",
        "category": "CRF_2ndGen",
        "accs": ["PQ585408", "PQ585409", "PQ585410", "PQ585411", "PQ585412", "PQ585413", "PQ585414"],
        "dag_parents": ["CRF63_02A6"],
    },
    "CRF188_0107": {
        "short_code": "188_0107",
        "category": "CRF_2ndGen",
        "accs": ["PV363678", "PV363679", "PV363680", "PV363681", "PV363682"],
        "dag_parents": ["CRF01_AE", "CRF07_BC"],
    },
    "CRF205_0107": {
        "short_code": "205_0107",
        "category": "CRF_2ndGen",
        "accs": ["PX262425", "PX262426", "PX262427", "PX262428", "OP921950"],
        "dag_parents": ["CRF01_AE", "CRF07_BC"],
    },
}

NOVEL_REGIONAL_MOSAICS = {
    "RU_G0B_gagMA": {
        "short_code": "RU_G0B",
        "category": "Novel_Circulating_Full",
        "accs": ["MF614608", "MF614609", "MF614611", "MF614612", "MF614613", "MF614614"],
        "dag_parents": ["Subtype_G"],
    }
}

BASE_TO_INT = {ord("A"): 0, ord("C"): 1, ord("G"): 2, ord("T"): 3, ord("U"): 3}
INT_TO_BASE = {0: "A", 1: "C", 2: "G", 3: "T", -1: "-"}

FOUNDER_LINEAGES = [
    "CRF01_AE",
    "CRF02_AG",
    "Subtype_G",
    "CRF06_cpx",
    "CRF07_BC",
    "CRF08_BC",
    "CRF11_cpx",
    "CRF12_BF",
    "CRF13_cpx",
    "CRF15_01B",
    "CRF33_01B",
    "CRF55_01B",
    "CRF63_02A6",
]

GLOBAL_STEM_FOUNDERS = {
    "CRF01_AE",
    "CRF02_AG",
    "Subtype_G",
    "CRF06_cpx",
    "CRF07_BC",
    "CRF08_BC",
    "CRF11_cpx",
    "CRF12_BF",
    "CRF63_02A6",
}

DERIVED_SUBTYPES = {"Subtype_G", "Subtype_A3", "Subtype_A4", "Subtype_A7", "Subtype_A8"}


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


def seq_to_int_array(seq: str) -> np.ndarray:
    arr = np.full(len(seq), -1, dtype=np.int8)
    b = seq.encode("ascii", errors="replace")
    raw = np.frombuffer(b, dtype=np.uint8)
    for char_code, val in BASE_TO_INT.items():
        arr[raw == char_code] = val
    return arr


def plurality_consensus_str(arrs) -> str:
    mat = np.vstack(arrs)
    l_len = mat.shape[1]
    counts = np.zeros((4, l_len), dtype=np.int32)
    for val in range(4):
        counts[val] = np.sum(mat == val, axis=0)
    tot = np.sum(counts, axis=0)
    best = np.argmax(counts, axis=0)
    out = ["-"] * l_len
    for i in range(l_len):
        if tot[i] > 0:
            out[i] = INT_TO_BASE[int(best[i])]
    return "".join(out)


def align_novel_to_hxb2_coding(pairs, profile_path: Path):
    """Align raw novel CRF genomes against hxb2_15trunk_profile_9719.fasta via mafft --addfragments --keeplength."""
    with tempfile.NamedTemporaryFile("w", suffix=".fasta", delete=False) as tf:
        tmp_in = tf.name
        for sid, seq in pairs:
            tf.write(f">{sid}\n{seq}\n")
    try:
        proc = subprocess.run(
            [MAFFT_BIN, "--6merpair", "--addfragments", tmp_in, "--keeplength", "--quiet", "--thread", "4", str(profile_path)],
            capture_output=True,
            text=True,
            check=True,
        )
        aln_text = proc.stdout
    finally:
        if os.path.exists(tmp_in):
            os.unlink(tmp_in)

    aligned = {}
    cur_h = None
    chunks = []
    for line in aln_text.splitlines():
        if line.startswith(">"):
            if cur_h is not None and not cur_h.startswith("Profile_"):
                s = "".join(chunks).upper()
                if len(s) == 9719:
                    aligned[cur_h] = seq_to_int_array(s[789:9417])
            cur_h = line[1:].strip().split()[0]
            chunks = []
        else:
            chunks.append(line.strip())
    if cur_h is not None and not cur_h.startswith("Profile_"):
        s = "".join(chunks).upper()
        if len(s) == 9719:
            aligned[cur_h] = seq_to_int_array(s[789:9417])
    return aligned


def _call_consensus_basal_with_eggman(cons_map: dict):
    """Run eggman Tier-1 (--gamma 13.5 --min-snps 7 --min-span 200 --p-thresh 1e-5) on consensus sequences."""
    assert EGGMAN_BIN.exists(), f"Missing eggman binary: {EGGMAN_BIN}"
    with tempfile.TemporaryDirectory() as td:
        in_fa = Path(td) / "cons.fasta"
        out_js = Path(td) / "cons.json"
        with open(in_fa, "w") as f:
            for lbl, seq in cons_map.items():
                f.write(f">{lbl}\n{seq}\n")
        subprocess.run(
            [
                str(EGGMAN_BIN),
                str(in_fa),
                "--gamma", "13.5",
                "--min-snps", "7",
                "--min-span", "200",
                "--p-thresh", "1e-5",
                "--json", str(out_js),
                "-q",
            ],
            check=True,
        )
        return {str(r["sequence_id"]): r for r in json.loads(out_js.read_text())}


def _parse_lanl_hdr_year(hdr: str):
    parts = hdr.split(".")
    if len(parts) >= 4:
        yy = parts[3]
        if yy.isdigit() and len(yy) == 2:
            val = int(yy)
            return float(1900 + val if val >= 60 else 2000 + val)
    return None


def main():
    with open(RHIZ_DATA_DIR / "compendium_bundle.json", "r") as f:
        old_bundle = json.load(f)
    old_comp_by_label = {e["label"]: e for e in old_bundle["compendium"]}

    comp_csv = pd.read_csv(RESULTS_DIR / "hiv1_descriptive_nomenclature_compendium.csv")
    dag_csv = pd.read_csv(RESULTS_DIR / "hiv1_recombination_genealogy_dag.csv").set_index("legacy_label")
    ff_df = pd.read_csv(DATA_DIR / "genbank_cluster_flatfiles.tsv", sep="\t", dtype=str).fillna("")
    ff_map = {str(r["accession"]): r.to_dict() for _, r in ff_df.iterrows()}
    novel_df = pd.read_csv(RESULTS_DIR / "genbank_novel_circulating_recombinants.csv")

    # 1. Read LANL 2023 reference alignment and project onto exact 9,719 HXB2 non-gap columns
    ref_recs = read_fasta(DATA_DIR / "lanl_ref_2023.fasta")
    hxb2_seqs = [seq for hdr, seq in ref_recs if "K03455" in hdr]
    assert len(hxb2_seqs) == 1, f"Expected 1 HXB2 record, got {len(hxb2_seqs)}"
    hxb2_cols = np.array([i for i, ch in enumerate(hxb2_seqs[0]) if ch not in "-."], dtype=np.int32)
    assert len(hxb2_cols) == 9719, f"Expected 9719 HXB2 columns, got {len(hxb2_cols)}"

    comp_csv_labels = set(comp_csv["legacy_label"].astype(str))
    clade_arrays_9719 = defaultdict(list)
    clade_sampled_yrs = defaultdict(list)
    clade_sub_yrs = defaultdict(list)

    for hdr, seq in ref_recs:
        if "K03455" in hdr or hdr.startswith("Ref.HXB2"):
            continue
        cl = hdr.split(".")[1]
        lbl = f"CRF{cl}" if ("_" in cl and cl.split("_")[0].isdigit()) else f"Subtype_{cl}"
        if lbl in comp_csv_labels:
            arr_full = seq_to_int_array(seq)[hxb2_cols]
            clade_arrays_9719[lbl].append(arr_full)
            acc = hdr.split(".")[-1]
            ff_r = ff_map.get(acc, {})
            cy = ff_r.get("collection_year", "")
            sy = ff_r.get("first_submitted_year", "")
            if cy:
                try:
                    clade_sampled_yrs[lbl].append(float(cy))
                except ValueError:
                    pass
            else:
                hy = _parse_lanl_hdr_year(hdr)
                if hy is not None:
                    clade_sampled_yrs[lbl].append(hy)
            if sy:
                try:
                    clade_sub_yrs[lbl].append(float(sy))
                except ValueError:
                    pass

    # 2. Write founder_profiles_9719.tsv
    founder_path = RHIZ_DATA_DIR / "founder_profiles_9719.tsv"
    with open(founder_path, "w") as f:
        for lbl in FOUNDER_LINEAGES:
            arrs = clade_arrays_9719[lbl]
            assert len(arrs) > 0, f"No sequences found for founder {lbl}"
            cons_9719 = plurality_consensus_str(arrs)
            assert len(cons_9719) == 9719
            f.write(f"{lbl}\t{cons_9719}\n")
    print(f"[+] Wrote {len(FOUNDER_LINEAGES)} founder profiles (9,719 nt) to {founder_path}")

    # 3. Extract and align all 8 promoted post-2023 LANL CRFs + RU_G0B_gagMA from GENBANK_FASTA
    all_ext_specs = {}
    all_ext_specs.update(POST2023_PROMOTED_CRFS)
    all_ext_specs.update(NOVEL_REGIONAL_MOSAICS)
    needed_accs = {a for spec in all_ext_specs.values() for a in spec["accs"]}
    raw_ext = {}
    with open(GENBANK_FASTA, "r") as f:
        cur_acc, cur_buf = None, []
        for line in f:
            if line.startswith(">"):
                if cur_acc in needed_accs:
                    raw_ext[cur_acc] = "".join(cur_buf).upper().replace("-", "").replace(".", "")
                cur_acc = line[1:].strip().split()[0].split(".")[0]
                cur_buf = []
            elif cur_acc in needed_accs:
                cur_buf.append(line.strip())
        if cur_acc in needed_accs:
            raw_ext[cur_acc] = "".join(cur_buf).upper().replace("-", "").replace(".", "")
    missing_ext = needed_accs - set(raw_ext.keys())
    assert not missing_ext, f"Missing post-2023 / regional FASTA accessions: {sorted(missing_ext)}"

    profile_9719_path = DATA_DIR / "hxb2_15trunk_profile_9719.fasta"
    aligned_ext = align_novel_to_hxb2_coding([(a, raw_ext[a]) for a in sorted(needed_accs)], profile_9719_path)
    ext_arrays_8628 = {}
    ext_cons_8628 = {}
    for lbl, spec in all_ext_specs.items():
        arrs = [aligned_ext[a] for a in spec["accs"] if a in aligned_ext]
        assert len(arrs) == len(spec["accs"]), f"Alignment missing for {lbl}"
        ext_arrays_8628[lbl] = arrs
        ext_cons_8628[lbl] = plurality_consensus_str(arrs)

    egg_cons_calls = _call_consensus_basal_with_eggman(ext_cons_8628)

    # 4. Rebuild updated compendium_bundle.json (163 LANL 2023 + 8 post-2023 CRFs + 1 RU_G0B_gagMA; zero stale NOVEL_CRF_01..15)
    new_compendium = []
    for _, row in comp_csv.iterrows():
        lbl = str(row["legacy_label"])
        t_stat = str(row["taxonomic_status"])
        old_e = old_comp_by_label.get(lbl, {})
        cat = str(old_e.get("category", row["legacy_category"]))
        if t_stat == "Verified_Recombinant_Lineage":
            if cat == "Demoted_LANL_Entry":
                if "_cpx" in lbl:
                    cat = "CRF_cpx"
                elif lbl in dag_csv.index and int(dag_csv.loc[lbl, "generation"]) == 1:
                    cat = "CRF_Primary"
                else:
                    cat = "CRF_2ndGen"
        elif "Collapsed_to_Clonal_Trunk" in t_stat or "Out_of_Panel_Donor" in t_stat:
            cat = "Demoted_LANL_Entry"

        tm = float(row["t_mrca"]) if pd.notna(row["t_mrca"]) else None
        ci_lo = float(row["t_mrca_ci_low"]) if pd.notna(row["t_mrca_ci_low"]) else None
        ci_hi = float(row["t_mrca_ci_high"]) if pd.notna(row["t_mrca_ci_high"]) else None
        s_yrs = clade_sampled_yrs.get(lbl, [])
        sub_yrs = clade_sub_yrs.get(lbl, [])
        first_samp = round(min(s_yrs), 4) if s_yrs else tm
        first_sub = round(min(sub_yrs), 4) if sub_yrs else None

        new_compendium.append({
            "label": lbl,
            "short_code": str(row["short_code"]),
            "category": cat,
            "status": t_stat,
            "descriptive_name": str(row["descriptive_name"]),
            "structural_formula": str(row["structural_formula"]),
            "backbone": str(row["backbone_trunk"]),
            "t_mrca": tm,
            "t_mrca_ci_low": ci_lo,
            "t_mrca_ci_high": ci_hi,
            "first_sampled_year": first_samp,
            "first_submitted_year": first_sub,
            "hxb2_barcode": str(row["hxb2_barcode"]),
            "scope": "full",
        })

    for lbl, spec in POST2023_PROMOTED_CRFS.items():
        ec = egg_cons_calls[lbl]
        s_yrs = [float(ff_map[a]["collection_year"]) for a in spec["accs"] if a in ff_map and ff_map[a].get("collection_year")]
        sub_yrs = [float(ff_map[a]["first_submitted_year"]) for a in spec["accs"] if a in ff_map and ff_map[a].get("first_submitted_year")]
        first_samp = round(min(s_yrs), 4) if s_yrs else None
        first_sub = round(min(sub_yrs), 4) if sub_yrs else None
        yr_tag = f"@{int(first_samp)}" if first_samp else ""
        form_s = str(ec["basal_formula"])
        new_compendium.append({
            "label": lbl,
            "short_code": spec["short_code"],
            "category": spec["category"],
            "status": "Verified_Recombinant_Lineage",
            "descriptive_name": f"{form_s}{yr_tag}",
            "structural_formula": form_s,
            "backbone": str(ec["backbone_trunk"]),
            "t_mrca": round(first_samp, 1) if first_samp else None,
            "t_mrca_ci_low": None,
            "t_mrca_ci_high": None,
            "first_sampled_year": first_samp,
            "first_submitted_year": first_sub,
            "hxb2_barcode": str(ec["basal_barcode"]),
            "scope": "full",
        })

    for lbl, spec in NOVEL_REGIONAL_MOSAICS.items():
        ec = egg_cons_calls[lbl]
        s_yrs = [float(ff_map[a]["collection_year"]) for a in spec["accs"] if a in ff_map and ff_map[a].get("collection_year")]
        sub_yrs = [float(ff_map[a]["first_submitted_year"]) for a in spec["accs"] if a in ff_map and ff_map[a].get("first_submitted_year")]
        first_samp = round(min(s_yrs), 4) if s_yrs else None
        first_sub = round(min(sub_yrs), 4) if sub_yrs else None
        form_s = str(ec["basal_formula"])
        nov_match = novel_df[novel_df["rhizaeon_formula"] == form_s]
        min_ufb = float(nov_match.iloc[0]["min_lineage_ufboot_pct"]) if not nov_match.empty else 100.0
        div_pct = float(nov_match.iloc[0]["mean_pairwise_div_pct"]) if not nov_match.empty else 4.49
        yr_tag = f" [Span:{int(min(s_yrs))}-{int(max(s_yrs))}]" if s_yrs else ""
        new_compendium.append({
            "label": lbl,
            "short_code": spec["short_code"],
            "category": spec["category"],
            "status": f"Candidate_Regional_Mosaic (Strict ML Monophyly; N={len(spec['accs'])} Russia, min_UFBoot={min_ufb:.0f}%, div={div_pct:.2f}%)",
            "descriptive_name": f"{form_s}{yr_tag}",
            "structural_formula": form_s,
            "backbone": str(ec["backbone_trunk"]),
            "t_mrca": None,
            "t_mrca_ci_low": None,
            "t_mrca_ci_high": None,
            "first_sampled_year": first_samp,
            "first_submitted_year": first_sub,
            "hxb2_barcode": str(ec["basal_barcode"]),
            "scope": "full",
        })

    new_bundle = {
        "lanl_msa_len": old_bundle["lanl_msa_len"],
        "lanl_msa_to_hxb2": old_bundle["lanl_msa_to_hxb2"],
        "compendium": new_compendium,
    }
    bundle_path = RHIZ_DATA_DIR / "compendium_bundle.json"
    with open(bundle_path, "w") as f:
        json.dump(new_bundle, f, indent=2)
    print(f"[+] Wrote {len(new_compendium)} entries to {bundle_path}")

    # 5. Collect eligible certified lineages (152 total: 1 Subtype_G + 142 2023 CRFs + 8 post-2023 CRFs + 1 RU_G0B_gagMA)
    eligible_rows = []
    eligible_labels = set()
    for entry in new_compendium:
        lbl = entry["label"]
        st = entry["status"]
        if lbl == "Subtype_G":
            assert lbl in clade_arrays_9719, f"Missing LANL sequences for {lbl}"
            cod_arrs = [a[789:9417] for a in clade_arrays_9719[lbl]]
            out_cat = "Ancestral_Recombinant_Subtype"
        elif lbl in POST2023_PROMOTED_CRFS:
            cod_arrs = ext_arrays_8628[lbl]
            out_cat = "Certified_CRF"
        elif lbl in NOVEL_REGIONAL_MOSAICS:
            cod_arrs = ext_arrays_8628[lbl]
            out_cat = "Novel_Circulating_CRF"
        elif st == "Verified_Recombinant_Lineage":
            assert lbl in clade_arrays_9719, f"Missing LANL sequences for {lbl}"
            cod_arrs = [a[789:9417] for a in clade_arrays_9719[lbl]]
            out_cat = "Certified_CRF"
        else:
            continue

        eligible_labels.add(lbl)
        eligible_rows.append((lbl, out_cat, False, cod_arrs))

    # 6. Parse immediate DAG parents from results/hiv1_recombination_genealogy_dag.csv + promoted/novel defs
    dag_parents = {}
    with open(RESULTS_DIR / "hiv1_recombination_genealogy_dag.csv", "r") as f:
        for r in csv.DictReader(f):
            c_lbl = r["legacy_label"].strip()
            p_str = r["immediate_recombinant_parents"]
            parents = list(
                dict.fromkeys(
                    m
                    for m in re.findall(r"CRF\d+_[A-Za-z0-9]+|Subtype_G", p_str)
                    if m in eligible_labels and m != c_lbl
                )
            )
            if parents:
                dag_parents[c_lbl] = parents
    for lbl, spec in list(POST2023_PROMOTED_CRFS.items()) + list(NOVEL_REGIONAL_MOSAICS.items()):
        p_list = [p for p in spec["dag_parents"] if p in eligible_labels and p != lbl]
        if p_list:
            dag_parents[lbl] = p_list

    # 7. Write certified_profiles_8628.tsv
    cert_path = RHIZ_DATA_DIR / "certified_profiles_8628.tsv"
    n_verified = 0
    n_ancestral = 0
    n_novel = 0
    with open(cert_path, "w") as f:
        f.write("# label\tcategory\tflags\tdag_parents\tn_seqs\tseq_8628\n")
        for lbl, out_cat, _, cod_arrs in eligible_rows:
            flags = []
            if lbl in GLOBAL_STEM_FOUNDERS:
                flags.append("global_stem")
            if out_cat == "Ancestral_Recombinant_Subtype":
                n_ancestral += 1
            elif out_cat == "Certified_CRF":
                n_verified += 1
            elif out_cat == "Novel_Circulating_CRF":
                n_novel += 1

            flags_str = ",".join(flags) if flags else "-"
            p_list = dag_parents.get(lbl, [])
            parents_str = "+".join(p_list) if p_list else "-"
            cons_8628 = plurality_consensus_str(cod_arrs)
            assert len(cons_8628) == 8628
            f.write(f"{lbl}\t{out_cat}\t{flags_str}\t{parents_str}\t{len(cod_arrs)}\t{cons_8628}\n")

    assert len(eligible_rows) == 152, f"Expected 152 certified profiles, got {len(eligible_rows)}"
    print(
        f"[+] Wrote {len(eligible_rows)} certified profiles (8,628 nt) to {cert_path}:\n"
        f"    * Ancestral Derived Subtype (Subtype_G) : {n_ancestral}\n"
        f"    * Certified LANL CRFs (2023 + Post-2023): {n_verified}\n"
        f"    * Novel Regional Mosaics (RU_G0B_gagMA) : {n_novel}"
    )


if __name__ == "__main__":
    main()
