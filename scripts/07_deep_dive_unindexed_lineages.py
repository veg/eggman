#!/usr/bin/env python3
"""
07_deep_dive_unindexed_lineages.py
==================================
Forensic deep dive into the 3 Unindexed Circulating Recombinant Lineages discovered
inside LANL's URF bin during our full-genome RhizAeon audit:
  1. East African D[A1:env] / D[A1:gp120..gp41] (N = 16, 1993–2015, UG/KE/TZ/GB/SE)
  2. East African C[A1:env] & C[A1:gp120..nef] (N = 9, 1995–2019, TZ/KE/RW/IN/SE)
  3. Russian G0[B:gag.MA] (N = 6, 2014–2018, Russia)

For each lineage, computes:
  - Exact breakpoint coordinates, informative-site CIs, and Fisher 2x2 p-values
  - Pairwise genetic divergence (p-distance / TN93) inside the backbone vs. donor cassette
    compared to background pure subtypes (testing Single Clonal Origin vs. Recurrent Hotspot Convergence)
  - Nearest-neighbor phylogenetic affinity in backbone and donor partitions
  - ChronAeon molecular clock t_MRCA and 95% Fieller confidence intervals
"""

import sys
import json
from pathlib import Path
from collections import defaultdict
import numpy as np
import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))
from importlib import import_module

s3 = import_module("03_audit_and_reclassify_all_crfs")
from chronaeon.dating_models import run_ols_dating

PROJECT_ROOT = Path("/Users/sergei/Projects/HIV-recombination")
DATA_DIR = PROJECT_ROOT / "data"
RESULTS_DIR = PROJECT_ROOT / "results"


def get_hxb2_coding_cols(seq_str: str) -> np.ndarray:
    cols = []
    h_pos = 0
    for i, ch in enumerate(seq_str):
        if ch not in "-.":
            h_pos += 1
            if 790 <= h_pos <= 9417:
                cols.append(i)
    return np.array(cols, dtype=np.int32)


def pairwise_pdist(arr1: np.ndarray, arr2: np.ndarray, mask: np.ndarray) -> float:
    a1 = arr1[mask]
    a2 = arr2[mask]
    valid = (a1 >= 0) & (a2 >= 0)
    if np.sum(valid) < 50:
        return np.nan
    return float(np.sum(a1[valid] != a2[valid])) / float(np.sum(valid))


def mean_pairwise_div(seq_list, mask):
    dists = []
    for i in range(len(seq_list)):
        for j in range(i + 1, len(seq_list)):
            d = pairwise_pdist(seq_list[i], seq_list[j], mask)
            if not np.isnan(d):
                dists.append(d)
    return np.mean(dists) if dists else np.nan, np.min(dists) if dists else np.nan, np.max(dists) if dists else np.nan


def main():
    ref_recs = s3.read_fasta(DATA_DIR / "lanl_ref_2023.fasta")
    flt_recs = s3.read_fasta(DATA_DIR / "lanl_flt_2022.fasta")

    ref_cols = get_hxb2_coding_cols([seq for hdr, seq in ref_recs if "K03455" in hdr][0])
    flt_cols = get_hxb2_coding_cols([seq for hdr, seq in flt_recs if "K03455" in hdr][0])
    hxb2_coords = np.arange(790, 9418, dtype=np.int32)

    clade_seqs = defaultdict(list)
    seq_by_name = {}
    for hdr, seq in ref_recs:
        if not hdr.startswith("Ref.HXB2"):
            arr = s3.seq_to_int_array(seq)[ref_cols]
            seq_by_name[hdr] = arr
            clade_seqs[hdr.split(".")[1]].append((hdr, arr))

    contam_df = pd.read_csv(RESULTS_DIR / "phase1_contaminated_pure_references.csv")
    contaminated_masks = defaultdict(list)
    for _, row in contam_df.iterrows():
        if bool(row["is_ghost"]):
            continue
        iso = row["recombinant_id"]
        c0 = int(np.searchsorted(hxb2_coords, int(row["hxb2_start"])))
        c1 = int(np.searchsorted(hxb2_coords, int(row["hxb2_end"])))
        contaminated_masks[iso].append((c0, c1))

    basal_trunk_clades = [
        ("A1", "A1"), ("A2", "A2"), ("A6", "A6"), ("B", "B"), ("C", "C"),
        ("D", "D"), ("F1", "F1"), ("F2", "F2"), ("G0", "G"), ("H", "H"),
        ("J", "J"), ("K", "K"), ("L", "L"),
    ]
    trunk_names, trunk_rows = [], []
    for t_name, src_clade in basal_trunk_clades:
        arrs = [arr for _, arr in clade_seqs[src_clade]]
        masks = [contaminated_masks.get(hdr, []) for hdr, _ in clade_seqs[src_clade]]
        trunk_names.append(t_name)
        trunk_rows.append(s3.build_plurality_consensus(arrs, masks))

    e_cons = s3.build_plurality_consensus([arr for _, arr in clade_seqs["01_AE"]])
    e_cons[hxb2_coords < 5097] = -2
    trunk_names.append("?E")
    trunk_rows.append(e_cons)
    trunk_mat = np.vstack(trunk_rows)
    cm240_arr = seq_by_name["Ref.01_AE.TH.90.CM240.U54771"]

    flt_arr_map = {}
    flt_by_tag = defaultdict(list)
    for hdr, seq in flt_recs:
        if "K03455" in hdr:
            continue
        arr = s3.seq_to_int_array(seq)[flt_cols]
        flt_arr_map[hdr] = arr
        flt_by_tag[hdr.split(".")[0]].append((hdr, arr))

    val_df = pd.read_csv(RESULTS_DIR / "lanl_all_full_genomes_rhizaeon_validation.csv")
    flt_val = val_df[(val_df["dataset"] == "LANL_Flt_2022") & (val_df["lanl_category"] == "Assigned_URF")]

    lineages = [
        ("Lineage_1_EastAfrica_D_A1_env", flt_val[flt_val["rhizaeon_formula"].isin(["D[A1:env]", "D[A1:gp120..gp41]"])],
         "D", "A1", (hxb2_coords >= 790) & (hxb2_coords <= 6000), (hxb2_coords >= 6350) & (hxb2_coords <= 8200)),
        ("Lineage_2_EastAfrica_C_A1_env", flt_val[flt_val["rhizaeon_formula"].isin(["C[A1:env]", "C[A1:gp120..nef]"])],
         "C", "A1", (hxb2_coords >= 790) & (hxb2_coords <= 6000), (hxb2_coords >= 7100) & (hxb2_coords <= 8250)),
        ("Lineage_3_Russia_G0_B_gagMA", flt_val[flt_val["rhizaeon_formula"] == "G0[B:gag.MA]"],
         "G", "B", (hxb2_coords >= 1200) & (hxb2_coords <= 9417), (hxb2_coords >= 790) & (hxb2_coords <= 1150)),
    ]

    detailed_rows = []
    for lin_name, sub_df, bb_tag, don_tag, bb_mask, don_mask in lineages:
        print("\n" + "=" * 96)
        print(f"DEEP DIVE: {lin_name} (N = {len(sub_df)} full-length genomes)")
        print("=" * 96)

        seq_ids = list(sub_df["sequence_id"])
        seq_arrs = [flt_arr_map[sid] for sid in seq_ids]
        years = np.array(sub_df["sampling_year"], dtype=np.float64)

        # 1. Exact breakpoint details for each member
        print("\n--- A. Member Isolates, Exact Breakpoints, Confidence Intervals & Fisher p-values ---")
        for sid in seq_ids:
            q_arr = flt_arr_map[sid]
            segs, bps = s3.rhizaeon_paff_viterbi(
                q_arr, trunk_mat, trunk_names, hxb2_coords, cm240_arr,
                gamma=13.5, mismatch_cost=1.0, ghost_cost=0.58,
                min_informative_snps=6, min_hxb2_span_nt=180, p_fisher_thresh=1e-4
            )
            bp_summary = "; ".join(
                f"{b['left_trunk']}->{b['right_trunk']}@{b['bp_hxb2']}[CI:{b['ci_hxb2_low']}-{b['ci_hxb2_high']}](n11={b['n11']},n22={b['n22']},p={b['p_fisher']:.1e})"
                for b in bps
            )
            r_meta = sub_df[sub_df["sequence_id"] == sid].iloc[0]
            print(f"  * {sid:<34s} | {r_meta['country']} ({int(r_meta['sampling_year'])}) | {r_meta['rhizaeon_formula']:<18s} | {bp_summary}")
            detailed_rows.append({
                "lineage_group": lin_name,
                "sequence_id": sid,
                "country": r_meta["country"],
                "sampling_year": r_meta["sampling_year"],
                "formula": r_meta["rhizaeon_formula"],
                "barcode": r_meta["hxb2_barcode"],
                "breakpoints_detail": bp_summary,
            })

        # 2. Pairwise divergence within cluster vs background pure subtype
        bg_bb_arrs = [arr for _, arr in flt_by_tag[bb_tag][:80]]
        bg_don_arrs = [arr for _, arr in flt_by_tag[don_tag][:80]]

        mean_bb, min_bb, max_bb = mean_pairwise_div(seq_arrs, bb_mask)
        mean_don, min_don, max_don = mean_pairwise_div(seq_arrs, don_mask)
        bg_mean_bb, _, _ = mean_pairwise_div(bg_bb_arrs, bb_mask)
        bg_mean_don, _, _ = mean_pairwise_div(bg_don_arrs, don_mask)

        print("\n--- B. Within-Lineage vs. Background Subtype Genetic Divergence (p-distance) ---")
        print(f"  * Backbone ({bb_tag}) Partition : Within-Cluster Mean = {mean_bb:.4f} (min={min_bb:.4f}, max={max_bb:.4f}) vs. Background {bb_tag} Mean = {bg_mean_bb:.4f}")
        print(f"  * Donor ({don_tag}) Partition    : Within-Cluster Mean = {mean_don:.4f} (min={min_don:.4f}, max={max_don:.4f}) vs. Background {don_tag} Mean = {bg_mean_don:.4f}")

        # 3. Identify tight transmission sub-clusters (pairwise distance < 0.08 or identical breakpoint)
        print("\n--- C. Nearest-Neighbor Pairwise Affinity Among Members ---")
        for i in range(len(seq_ids)):
            best_j = -1
            best_d = 999.0
            for j in range(len(seq_ids)):
                if i == j:
                    continue
                d = pairwise_pdist(seq_arrs[i], seq_arrs[j], bb_mask | don_mask)
                if d < best_d:
                    best_d = d
                    best_j = j
            print(f"  * {seq_ids[i]:<34s} -> Closest cluster peer: {seq_ids[best_j]:<34s} (whole-genome dist = {best_d:.4f})")

        # 4. ChronAeon / Root-to-Tip Dating against earliest member / ancestral consensus
        valid_yr = ~np.isnan(years)
        if np.sum(valid_yr) >= 4 and (np.max(years[valid_yr]) - np.min(years[valid_yr])) >= 3.0:
            earliest_idx = int(np.argmin(years[valid_yr]))
            anc_cons = seq_arrs[earliest_idx]
            divs = np.array([pairwise_pdist(arr, anc_cons, bb_mask | don_mask) for arr in seq_arrs], dtype=np.float64)
            res = run_ols_dating(times=years[valid_yr], dists=divs[valid_yr], ci_method="fieller", seq_len=8000)
            ci = res.get("ci_mrca", [np.nan, np.nan])
            print("\n--- D. ChronAeon Molecular Clock Calibration ---")
            print(f"  * t_MRCA = {res['t_mrca']:.1f} [95% Fieller CI: {ci[0]:.1f} – {ci[1]:.1f}], mu = {res['mu']:.2e} subs/site/yr, r = {res['r']:.3f}")

    pd.DataFrame(detailed_rows).to_csv(RESULTS_DIR / "unindexed_circulating_lineages_deep_dive.csv", index=False)


if __name__ == "__main__":
    main()
