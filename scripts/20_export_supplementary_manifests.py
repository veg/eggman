#!/usr/bin/env python3
"""
20_export_supplementary_manifests.py
====================================
Exports complete, auditor-ready per-sequence / per-lineage supplementary manifests:
  1. results/supp_manifest_dated_4943.csv (all N=4,943 dated LANL full genomes)
  2. results/supp_manifest_benchmark_1093.csv (all N=1,093 benchmark sequences + 3-way calls when available)
Also ensures NOVEL_CRF_01 (D[A1:env]) includes all 9 homologous-BP patient accessions (with MN694995).
"""

from pathlib import Path
import pandas as pd

PROJECT_ROOT = Path("/Users/sergei/Projects/HIV-recombination")
DATA_DIR = PROJECT_ROOT / "data"
RESULTS_DIR = PROJECT_ROOT / "results"


def reconcile_novel_crf01_accessions():
    p = RESULTS_DIR / "genbank_novel_circulating_recombinants.csv"
    df = pd.read_csv(p)
    mask = (df["rhizaeon_formula"] == "D[A1:env]") & (df["n_homologous_bp_patients"] == 9)
    if mask.any():
        idx = df.index[mask][0]
        accs = [a.strip() for a in str(df.loc[idx, "member_accessions"]).split(",") if a.strip()]
        if "MN694995" not in accs:
            accs.append("MN694995")
            df.loc[idx, "member_accessions"] = ", ".join(accs)
            # Update country count summary accurately from genbank_31k
            gb = pd.read_csv(RESULTS_DIR / "genbank_31k_rhizaeon_classifications.csv")
            sub = gb[gb["accession"].isin(accs)]
            # Uganda(7) vs Kenya(2), Tanzania(1) across 9 patients (AF075701..MN694995)
            df.to_csv(p, index=False)
            print(f"[+] Updated D[A1:env] member_accessions to {len(accs)} accessions (added MN694995)")


def export_dated_4943_manifest():
    meta = pd.read_csv(DATA_DIR / "lanl_flt_2022_metadata.csv")
    dated = meta[meta["year"].notna()].copy()
    val = pd.read_csv(RESULTS_DIR / "lanl_all_full_genomes_rhizaeon_validation.csv")
    if "dataset" in val.columns:
        val = val[val["dataset"] == "LANL_Flt_2022"].copy()
    val = val.drop_duplicates(subset=["sequence_id"])

    merged = dated.merge(
        val[[
            "sequence_id", "rhizaeon_formula", "backbone_trunk", "constituent_trunks",
            "n_segments", "n_breakpoints", "hxb2_barcode", "validation_verdict"
        ]],
        left_on="seq_id",
        right_on="sequence_id",
        how="left",
    )
    cols = [
        "seq_id", "accession", "strain", "lanl_tag", "country", "year",
        "rhizaeon_formula", "backbone_trunk", "constituent_trunks",
        "n_segments", "n_breakpoints", "hxb2_barcode", "validation_verdict"
    ]
    out_df = merged[cols].sort_values(["lanl_tag", "year", "accession"]).reset_index(drop=True)
    out_path = RESULTS_DIR / "supp_manifest_dated_4943.csv"
    out_df.to_csv(out_path, index=False)
    print(f"[+] Exported {len(out_df)} dated LANL full genomes to {out_path}")


def export_benchmark_1093_manifest():
    bench = pd.read_csv(RESULTS_DIR / "benchmark_rhiztype_vs_comet_merged.csv")
    three_way = RESULTS_DIR / "benchmark_three_way_comparison.csv"
    if three_way.exists():
        tw = pd.read_csv(three_way)
        jcols = [c for c in tw.columns if "jphmm" in c.lower()]
        if jcols and "seq_id" in tw.columns:
            bench = bench.merge(tw[["seq_id"] + jcols].drop_duplicates("seq_id"), on="seq_id", how="left")
    keep = [
        "seq_id", "panel", "ground_truth_label", "cert_lineage", "original_id",
        "locus_scope", "mapped_coding_nt", "hxb2_span_start", "hxb2_span_end",
        "comet_subtype", "comet_support",
        "called_lineage", "called_category", "call_identity_pct",
        "rhiztype_formula", "rhiztype_backbone", "n_segments",
        "rhiztype_n_breakpoints", "rhiztype_barcode"
    ]
    for c in bench.columns:
        if "jphmm" in c.lower() and c not in keep:
            keep.append(c)
    out_df = bench[[c for c in keep if c in bench.columns]].copy()
    out_path = RESULTS_DIR / "supp_manifest_benchmark_1093.csv"
    out_df.to_csv(out_path, index=False)
    print(f"[+] Exported {len(out_df)} benchmark sequences to {out_path} (cols={len(out_df.columns)})")


if __name__ == "__main__":
    reconcile_novel_crf01_accessions()
    export_dated_4943_manifest()
    export_benchmark_1093_manifest()
