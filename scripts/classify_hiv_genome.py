#!/usr/bin/env python3
"""
classify_hiv_genome.py
======================
Standalone CLI classifier for HIV-1 Group M genomes under the RhizAeon + ChronAeon
Descriptive Spatio-Temporal Recombination Nomenclature.

Usage:
    python3 scripts/classify_hiv_genome.py --fasta <input.fasta> [--limit 10] [--out <out.csv>]
"""

import argparse
import json
import sys
from pathlib import Path
from collections import defaultdict, Counter
import numpy as np
import pandas as pd

# Import core RhizAeon PA-FF v2.1 engine from 03_audit_and_reclassify_all_crfs.py
SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))
from importlib import import_module

step3 = import_module("03_audit_and_reclassify_all_crfs")

PROJECT_ROOT = Path("/Users/sergei/Projects/HIV-recombination")
DATA_DIR = PROJECT_ROOT / "data"
RESULTS_DIR = PROJECT_ROOT / "results"


def load_clean_trunk_panel():
    records = step3.read_fasta(DATA_DIR / "lanl_ref_2023.fasta")
    with open(DATA_DIR / "lanl_ref_2023_hxb2_coord_map.json") as f:
        cmap = json.load(f)

    m2h = {int(k) - 1: int(v) for k, v in cmap["msa_to_hxb2"].items()}
    L_full = len(records[0][1])
    col_to_hxb2_full = np.array([m2h.get(i, 0) for i in range(L_full)], dtype=np.int32)
    coding_mask = (col_to_hxb2_full >= 790) & (col_to_hxb2_full <= 9417)
    coding_cols = np.where(coding_mask)[0]
    col_to_hxb2 = col_to_hxb2_full[coding_cols]

    clade_seqs = defaultdict(list)
    seq_by_name = {}
    for hdr, seq in records:
        if hdr.startswith("Ref.HXB2"):
            continue
        arr = step3.seq_to_int_array(seq)[coding_cols]
        seq_by_name[hdr] = arr
        parts = hdr.split(".")
        cl = parts[1] if len(parts) > 1 else "UNK"
        clade_seqs[cl].append((hdr, arr))

    contam_df = pd.read_csv(RESULTS_DIR / "phase1_contaminated_pure_references.csv")
    contaminated_masks = defaultdict(list)
    for _, row in contam_df.iterrows():
        if bool(row["is_ghost"]):
            continue
        iso = row["recombinant_id"]
        h0, h1 = int(row["hxb2_start"]), int(row["hxb2_end"])
        c0 = int(np.searchsorted(col_to_hxb2, h0))
        c1 = int(np.searchsorted(col_to_hxb2, h1))
        contaminated_masks[iso].append((c0, c1))

    basal_trunk_clades = [
        ("A1", "A1"), ("A2", "A2"), ("A6", "A6"), ("B", "B"), ("C", "C"),
        ("D", "D"), ("F1", "F1"), ("F2", "F2"), ("G0", "G"), ("H", "H"),
        ("J", "J"), ("K", "K"), ("L", "L"),
    ]
    trunk_names = []
    trunk_rows = []
    for t_name, src_clade in basal_trunk_clades:
        arrs = [arr for _, arr in clade_seqs[src_clade]]
        masks = [contaminated_masks.get(hdr, []) for hdr, _ in clade_seqs[src_clade]]
        cons = step3.build_plurality_consensus(arrs, masks)
        trunk_names.append(t_name)
        trunk_rows.append(cons)

    e_cons = step3.build_plurality_consensus([arr for _, arr in clade_seqs["01_AE"]])
    e_cons[col_to_hxb2 < 5097] = -2
    trunk_names.append("?E")
    trunk_rows.append(e_cons)

    trunk_mat = np.vstack(trunk_rows)
    cm240_arr = seq_by_name["Ref.01_AE.TH.90.CM240.U54771"]
    comp_df = pd.read_csv(RESULTS_DIR / "hiv1_descriptive_nomenclature_compendium.csv")
    return trunk_mat, trunk_names, coding_cols, col_to_hxb2, cm240_arr, comp_df, L_full


def main():
    parser = argparse.ArgumentParser(description="Classify HIV-1 Group M genomes using RhizAeon + ChronAeon nomenclature.")
    parser.add_argument("--fasta", type=str, required=True, help="Path to aligned HIV-1 FASTA file")
    parser.add_argument("--limit", type=int, default=10, help="Max sequences to classify (0 for all)")
    parser.add_argument("--out", type=str, default="", help="Optional output CSV path")
    args = parser.parse_args()

    trunk_mat, trunk_names, coding_cols, col_to_hxb2, cm240_arr, comp_df, L_full = load_clean_trunk_panel()
    formula_to_legacy = defaultdict(list)
    for _, r in comp_df.iterrows():
        formula_to_legacy[r["structural_formula"]].append(f"{r['legacy_label']} ({r['descriptive_name']})")

    records = step3.read_fasta(Path(args.fasta))
    if args.limit > 0:
        records = records[:args.limit]

    out_rows = []
    print("=" * 96)
    print(f"RHIZAEON + CHRONAEON SPATIO-TEMPORAL HIV-1 CLASSIFIER (N={len(records)} sequences)")
    print("=" * 96)

    for hdr, seq in records:
        if len(seq) == L_full:
            q_arr = step3.seq_to_int_array(seq)[coding_cols]
        elif len(seq) == len(coding_cols):
            q_arr = step3.seq_to_int_array(seq)
        else:
            raise ValueError(f"Sequence {hdr} length ({len(seq)}) does not match reference alignment ({L_full} or {len(coding_cols)}).")

        segs, bps = step3.rhizaeon_paff_viterbi(
            q_arr, trunk_mat, trunk_names, col_to_hxb2, cm240_arr,
            gamma=10.5, mismatch_cost=1.0, ghost_cost=0.58,
            min_informative_snps=5, min_hxb2_span_nt=120, p_fisher_thresh=1e-4
        )
        desc_name, backbone, struct_formula = step3.format_descriptive_name(segs, None)
        barcode = step3.format_hxb2_barcode(segs)
        ci_barcode = "|".join(
            f"{s['trunk']}:{s['hxb2_start']}-{s['hxb2_end']}[CI:{s['ci_hxb2_start']}-{s['ci_hxb2_end']}]"
            for s in segs
        )
        matches = formula_to_legacy.get(struct_formula, ["Novel_URF_or_Unindexed_Mosaic"])

        print(f"\n> {hdr}")
        print(f"  Descriptive Formula : {struct_formula}")
        print(f"  HXB2 Barcode        : {barcode}")
        print(f"  CI-Bracketed Map    : {ci_barcode}")
        print(f"  Compendium Match    : {', '.join(matches[:4])}")
        if bps:
            bp_strs = [f"HXB2 {b['bp_hxb2']} ({b['left_trunk']}->{b['right_trunk']}, CI {b['ci_hxb2_low']}-{b['ci_hxb2_high']}, p={b['p_fisher']:.2e})" for b in bps]
            print(f"  Verified Junctions  : {'; '.join(bp_strs)}")

        out_rows.append({
            "sequence_id": hdr,
            "structural_formula": struct_formula,
            "backbone_trunk": backbone,
            "n_segments": len(segs),
            "n_breakpoints": len(bps),
            "hxb2_barcode": barcode,
            "ci_bracketed_barcode": ci_barcode,
            "compendium_matches": ";".join(matches),
        })

    if args.out:
        pd.DataFrame(out_rows).to_csv(args.out, index=False)
        print(f"\n[+] Saved classification results to {args.out}")


if __name__ == "__main__":
    main()
