#!/usr/bin/env python3
"""
06_validate_all_lanl_full_genomes.py
====================================
Exhaustive RhizAeon validation of all full-length HIV-1 Group M genomes in LANL:
  1. Official LANL Subtype & CRF Reference Panel (`lanl_ref_2023.fasta`, N = 517)
  2. Official LANL Filtered Web Alignment (`lanl_flt_2022.fasta`, N = 5,088)
     Total = 5,605 full-length HIV-1 genomes across 279 LANL assigned labels.

Evaluates for every genome:
  - Its RhizAeon Descriptive Spatio-Temporal Formula and HXB2 Barcode
  - Concordance vs. discordance with its LANL-assigned subtype, CRF, or URF label
  - Identification of unrecognized recombinants (cryptic URFs/CRFs) inside nominal "Pure" subtypes
  - Validation of declared parental composition across 987 LANL URFs and 1,142 LANL CRFs
"""

import json
import sys
import time
from pathlib import Path
from collections import defaultdict, Counter
from multiprocessing import Pool, cpu_count
import numpy as np
import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))
from importlib import import_module

s3 = import_module("03_audit_and_reclassify_all_crfs")

PROJECT_ROOT = Path("/Users/sergei/Projects/HIV-recombination")
DATA_DIR = PROJECT_ROOT / "data"
RESULTS_DIR = PROJECT_ROOT / "results"

# Global variables shared across worker processes
G_TRUNK_MAT = None
G_TRUNK_NAMES = None
G_HXB2_COORDS = None
G_CM240_ARR = None


def get_hxb2_coding_cols(seq_str: str) -> np.ndarray:
    cols = []
    h_pos = 0
    for i, ch in enumerate(seq_str):
        if ch not in "-.":
            h_pos += 1
            if 790 <= h_pos <= 9417:
                cols.append(i)
    return np.array(cols, dtype=np.int32)


def init_worker(trunk_mat, trunk_names, hxb2_coords, cm240_arr):
    global G_TRUNK_MAT, G_TRUNK_NAMES, G_HXB2_COORDS, G_CM240_ARR
    G_TRUNK_MAT = trunk_mat
    G_TRUNK_NAMES = trunk_names
    G_HXB2_COORDS = hxb2_coords
    G_CM240_ARR = cm240_arr


def classify_single_task(task):
    dataset_name, hdr, lanl_tag, country, year, q_arr = task
    segs, bps = s3.rhizaeon_paff_viterbi(
        q_arr, G_TRUNK_MAT, G_TRUNK_NAMES, G_HXB2_COORDS, G_CM240_ARR,
        gamma=13.5, mismatch_cost=1.0, ghost_cost=0.58,
        min_informative_snps=6, min_hxb2_span_nt=180, p_fisher_thresh=1e-4
    )
    _, backbone, struct_formula = s3.format_descriptive_name(segs, None)
    barcode = s3.format_hxb2_barcode(segs)
    constituents = sorted(set(s["trunk"] for s in segs))
    return {
        "dataset": dataset_name,
        "sequence_id": hdr,
        "lanl_assigned_tag": lanl_tag,
        "country": country,
        "sampling_year": year,
        "rhizaeon_formula": struct_formula,
        "backbone_trunk": backbone,
        "constituent_trunks": "+".join(constituents),
        "n_segments": len(segs),
        "n_breakpoints": len(bps),
        "hxb2_barcode": barcode,
    }


def parse_lanl_tag_category(tag: str, crf_short_tags: set, pure_subtypes: set):
    if tag in pure_subtypes or tag == "A":
        return "Assigned_Pure_Subtype"
    if tag in crf_short_tags:
        return "Assigned_Named_CRF"
    if tag in ("U", "?", "CPX"):
        return "Assigned_Unclassified_U"
    return "Assigned_URF"


def evaluate_concordance(row, comp_constituents_map, pure_subtypes):
    tag = str(row["lanl_assigned_tag"])
    cat = row["lanl_category"]
    inferred_set = set(row["constituent_trunks"].split("+"))
    n_bps = int(row["n_breakpoints"])
    backbone = row["backbone_trunk"]

    # Map A -> A1/A6/A2, G -> G0, F -> F1/F2
    def normalize_expected(t):
        if t in ("A", "A1", "A3", "A4", "A8"):
            return "A1"
        if t == "G":
            return "G0"
        if t == "F":
            return "F1"
        return t

    if cat == "Assigned_Pure_Subtype":
        exp_trunk = normalize_expected(tag)
        if n_bps == 0 and backbone == exp_trunk:
            return "Concordant_Clonal_Subtype"
        elif tag == "G" and backbone == "G0" and inferred_set.issubset({"G0", "A1", "J"}):
            return "Concordant_G_Lineage (G0 or G0+A1)"
        elif tag == "A7" and "A1" in inferred_set:
            return "Confirmed_Mosaic_SubSubtype (A7)"
        elif n_bps > 0 and exp_trunk in inferred_set:
            return "Cryptic_Recombinant (Unrecognized Mosaic in 'Pure' Subtype)"
        elif n_bps == 0 and backbone != exp_trunk:
            return f"Mislabeled_Clonal_Subtype (Assigned {tag}, Actual {backbone})"
        else:
            return "Discordant_Complex_Mosaic"

    if cat == "Assigned_Named_CRF":
        exp_set = comp_constituents_map.get(tag, set())
        if inferred_set == exp_set:
            return "Exact_Constituent_Match"
        elif inferred_set.issubset(exp_set) and len(inferred_set) >= 1:
            return "Partial_or_Backcrossed_CRF_Match"
        elif exp_set and exp_set.issubset(inferred_set):
            return "Secondary_Recombinant_on_CRF_Backbone (CRF + Additional Donor)"
        else:
            return "Discordant_CRF_Assignment"

    if cat == "Assigned_URF":
        # Extract expected letters from URF tag (e.g. "A1D" -> A1, D; "01B" -> A1, ?E, B; "BF1" -> B, F1)
        exp_urf = set()
        if "01" in tag or "AE" in tag:
            exp_urf.update(["A1", "?E"])
        if "02" in tag or "AG" in tag:
            exp_urf.update(["A1", "G0"])
        if "07" in tag or "08" in tag or "BC" in tag:
            exp_urf.update(["B", "C"])
        if "06" in tag:
            exp_urf.update(["G0", "J", "K", "A1"])
        if "A1" in tag or ("A" in tag and "A2" not in tag and "A6" not in tag):
            exp_urf.add("A1")
        if "A2" in tag:
            exp_urf.add("A2")
        if "A6" in tag:
            exp_urf.add("A6")
        if "B" in tag:
            exp_urf.add("B")
        if "C" in tag:
            exp_urf.add("C")
        if "D" in tag:
            exp_urf.add("D")
        if "F1" in tag or ("F" in tag and "F2" not in tag):
            exp_urf.add("F1")
        if "F2" in tag:
            exp_urf.add("F2")
        if "G" in tag:
            exp_urf.add("G0")
        if "H" in tag:
            exp_urf.add("H")
        if "J" in tag:
            exp_urf.add("J")
        if "K" in tag:
            exp_urf.add("K")

        if n_bps == 0:
            return f"Collapsed_to_Clonal_Trunk (Mislabeled as URF {tag}, Actual {backbone})"
        elif exp_urf and inferred_set == exp_urf:
            return "Exact_URF_Parental_Match"
        elif exp_urf and len(inferred_set & exp_urf) >= 1:
            return "Overlapping_URF_Parental_Match (Refined Donors)"
        else:
            return "Discordant_URF_Parental_Label"

    return "Resolved_Unclassified_Genome"


def main():
    t0 = time.time()
    print("=" * 95)
    print("STEP 6: EXHAUSTIVE RHIZAEON VALIDATION OF ALL 5,605 FULL-LENGTH LANL HIV-1 GENOMES")
    print("=" * 95)

    ref_recs = s3.read_fasta(DATA_DIR / "lanl_ref_2023.fasta")
    flt_recs = s3.read_fasta(DATA_DIR / "lanl_flt_2022.fasta")

    ref_hxb2_seq = [seq for hdr, seq in ref_recs if "K03455" in hdr][0]
    flt_hxb2_seq = [seq for hdr, seq in flt_recs if "K03455" in hdr][0]
    ref_cols = get_hxb2_coding_cols(ref_hxb2_seq)
    flt_cols = get_hxb2_coding_cols(flt_hxb2_seq)
    trunk_mat, trunk_names, hxb2_coords, cm240_arr = s3.load_profiles_9719_coding()

    loop_mask_bool = np.zeros(len(hxb2_coords), dtype=bool)
    for h_lo, h_hi in s3.HYPERVARIABLE_MASK_RANGES:
        loop_mask_bool |= (hxb2_coords >= h_lo) & (hxb2_coords <= h_hi)

    ref_meta = pd.read_csv(DATA_DIR / "lanl_ref_2023_metadata.csv").set_index("seq_id")
    flt_meta = pd.read_csv(DATA_DIR / "lanl_flt_2022_metadata.csv").set_index("seq_id")

    tasks = []
    for hdr, seq in ref_recs:
        if hdr.startswith("Ref.HXB2"):
            continue
        m = ref_meta.loc[hdr] if hdr in ref_meta.index else {}
        tag = m.get("lanl_tag", hdr.split(".")[1])
        country = m.get("country", "")
        year = m.get("year", np.nan)
        q_a = s3.seq_to_int_array(seq)[ref_cols]
        q_a[loop_mask_bool] = -1
        tasks.append(("LANL_Ref_2023", hdr, str(tag), str(country), year, q_a))

    for hdr, seq in flt_recs:
        if "K03455" in hdr:
            continue
        m = flt_meta.loc[hdr] if hdr in flt_meta.index else {}
        tag = m.get("lanl_tag", hdr.split(".")[0])
        country = m.get("country", "")
        year = m.get("year", np.nan)
        q_a = s3.seq_to_int_array(seq)[flt_cols]
        q_a[loop_mask_bool] = -1
        tasks.append(("LANL_Flt_2022", hdr, str(tag), str(country), year, q_a))

    n_workers = min(8, cpu_count())
    print(f"[+] Dispatching {len(tasks)} full-length HIV-1 genomes across {n_workers} parallel workers...")
    with Pool(processes=n_workers, initializer=init_worker, initargs=(trunk_mat, trunk_names, hxb2_coords, cm240_arr)) as pool:
        results = pool.map(classify_single_task, tasks, chunksize=60)

    df = pd.DataFrame(results)

    # Annotate concordance against LANL assigned labels
    comp_df = pd.read_csv(RESULTS_DIR / "hiv1_descriptive_nomenclature_compendium.csv")
    pure_subtypes = {"A1", "A2", "A3", "A4", "A6", "A7", "A8", "B", "C", "D", "F1", "F2", "G", "H", "J", "K", "L"}
    crf_short_tags = set(comp_df[comp_df["legacy_category"] != "Nominal_Pure_Subtype"]["short_code"])
    with open(DATA_DIR / "lanl_official_crf_catalog.json") as f:
        crf_cat = json.load(f)
    for k, v in crf_cat.items():
        crf_short_tags.add(v["short_tag"])

    comp_constituents_map = {}
    for _, r in comp_df.iterrows():
        donors = [d for d in str(r["donor_trunks"]).split(",") if d and d != "None"]
        comp_constituents_map[r["short_code"]] = set([r["backbone_trunk"]] + donors)

    df["lanl_category"] = df["lanl_assigned_tag"].map(lambda t: parse_lanl_tag_category(t, crf_short_tags, pure_subtypes))
    df["validation_verdict"] = [evaluate_concordance(row, comp_constituents_map, pure_subtypes) for _, row in df.iterrows()]

    out_path = RESULTS_DIR / "lanl_all_full_genomes_rhizaeon_validation.csv"
    df.to_csv(out_path, index=False)
    elapsed = time.time() - t0
    print(f"[+] Completed full-genome validation of {len(df)} sequences in {elapsed:.1f}s -> {out_path.name}")

    # Summary statistics for LANL_Flt_2022 (N=5,088)
    flt_df = df[df["dataset"] == "LANL_Flt_2022"]
    print("\n" + "=" * 95)
    print(f"SUMMARY OF LANL FULL-LENGTH FILTERED WEB ALIGNMENT (N = {len(flt_df)} GENOMES)")
    print("=" * 95)
    for cat, grp in flt_df.groupby("lanl_category"):
        print(f"\n### {cat} (N = {len(grp)}) ###")
        for verd, cnt in grp["validation_verdict"].value_counts().items():
            print(f"  * {verd:<72s}: {cnt:4d} ({100.0*cnt/len(grp):5.1f}%)")

    # Breakdown by specific Pure Subtypes in LANL_Flt_2022
    print("\n--- BREAKDOWN OF ASSIGNED 'PURE SUBTYPE' GENOMES IN LANL WEB ALIGNMENT ---")
    pure_grp = flt_df[flt_df["lanl_category"] == "Assigned_Pure_Subtype"]
    summary_pure = []
    for tag, grp in pure_grp.groupby("lanl_assigned_tag"):
        n_tot = len(grp)
        n_clonal = int(np.sum(grp["n_breakpoints"] == 0))
        n_mosaic = int(np.sum(grp["n_breakpoints"] > 0))
        top_mosaics = grp[grp["n_breakpoints"] > 0]["rhizaeon_formula"].value_counts().head(3).to_dict()
        summary_pure.append({
            "lanl_subtype": tag,
            "n_genomes": n_tot,
            "n_truly_clonal": n_clonal,
            "pct_clonal": round(100.0 * n_clonal / n_tot, 1),
            "n_cryptic_recombinant": n_mosaic,
            "pct_cryptic_recombinant": round(100.0 * n_mosaic / n_tot, 1),
            "top_cryptic_mosaic_formulas": "; ".join(f"{k} (n={v})" for k, v in top_mosaics.items()),
        })
    pure_sum_df = pd.DataFrame(summary_pure).sort_values("n_genomes", ascending=False)
    pure_sum_df.to_csv(RESULTS_DIR / "lanl_pure_subtype_genome_validation_summary.csv", index=False)
    print(pure_sum_df.to_string(index=False))


if __name__ == "__main__":
    main()
