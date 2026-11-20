#!/usr/bin/env python3
"""
Script 15: Comprehensive Subtype & Recombination Audit of the Stanford HIV Drug Resistance Database (HIVDB).

Audits two official Stanford HIVDB resources using RhizType:
  1. Stanford HIVDB Genotyper Reference Panel (N = 720 reference pol sequences, `genotype-references.fasta`)
  2. Stanford HIV Drug Resistance Surveillance Database (`hiv-survdb-payload`, N = 120,468 clinical sequences
     across 517 published studies linked to 343,564 isolate rows in `tables/isolates.d/*.csv`).
"""

import re
import subprocess
import time
from pathlib import Path
import numpy as np
import pandas as pd

PROJECT_ROOT = Path("/Users/sergei/Projects/HIV-recombination")
DATA_DIR = PROJECT_ROOT / "data" / "stanford_hivdb"
SURV_DIR = DATA_DIR / "hiv-survdb-payload"
RESULTS_DIR = PROJECT_ROOT / "results"
RHIZTYPE_BIN = PROJECT_ROOT / "rhiztype" / "target" / "release" / "eggman"

PURE_STANFORD_MAP = {
    "A": {"A1", "A2", "A6"},
    "A1": {"A1"},
    "A2": {"A2"},
    "A3": {"A1"},
    "A4": {"A1"},
    "A6": {"A6"},
    "A7": {"A1"},
    "A8": {"A1"},
    "B": {"B"},
    "C": {"C"},
    "D": {"D"},
    "F": {"F1", "F2"},
    "F1": {"F1"},
    "F2": {"F2"},
    "G": {"G0", "G"},
    "H": {"H"},
    "J": {"J"},
    "K": {"K"},
    "L": {"L"},
}


def _crf_num(label: str) -> str:
    m = re.match(r"^CRF0*(\d+)", str(label).strip(), re.IGNORECASE)
    return m.group(1) if m else ""


def classify_stanford_vs_rhiztype(
    stanford_subtypes_set: str,
    n_distinct_stanford: int,
    rhiz_formula: str,
    rhiz_backbone: str,
    rhiz_matches: str,
    n_bp: int,
    called_lineage: str = "",
    called_category: str = "",
) -> str:
    """Classify the relationship between Stanford HIVDB's subtype(s) and RhizType v0.2 Tier-1/Tier-2 calls."""
    st = str(stanford_subtypes_set).strip()
    rf = str(rhiz_formula).strip()
    rb = str(rhiz_backbone).strip()
    rm = str(rhiz_matches).strip()
    cl = str(called_lineage).strip()
    cc = str(called_category).strip()

    if st == "Unlinked_No_Isolate_Row":
        return "Stanford_Unlinked_No_Isolate_Row"

    # 1. Intra-patient gene-row contradiction in Stanford HIVDB (split locus-expected recombinants vs true clonal contradictions)
    if n_distinct_stanford > 1:
        if n_bp >= 1 or cc in ("Certified_CRF", "Novel_Circulating_CRF", "URF"):
            return "Stanford_MultiEnzyme_LocusExpected_Recombinant"
        if cc in ("Unclassifiable", "OutOfPanel") and n_bp == 0:
            return "Unclassifiable_or_Short_Fragment"
        return "Stanford_IntraPatient_True_Clonal_Contradiction"

    if st.startswith("Group ") or st.startswith("HIV2") or st in ("O", "N", "P", "CPZ"):
        return "Non_Group_M"

    if cc in ("Unclassifiable", "OutOfPanel") and n_bp == 0:
        return "Unclassifiable_or_Short_Fragment"

    # 2. Ambiguous / Unknown / URF label in Stanford HIVDB ('Unknown', 'B + C', 'BC', 'B + CRF01_AE', etc.)
    if st in ("Unknown", "URF", "BC", "BF", "AB", "AD", "CD", "AG", "AE") or "+" in st:
        return "Stanford_Ambiguous_or_Unknown_Resolved"

    # 3. Single Pure Subtype in Stanford HIVDB (accept Pure_Subtype & G Ancestral_Recombinant_Subtype)
    if st in PURE_STANFORD_MAP:
        allowed_trunks = PURE_STANFORD_MAP[st]
        if cc in ("Pure_Subtype", "Ancestral_Recombinant_Subtype") and (cl in allowed_trunks or rb in allowed_trunks):
            return "Concordant_Pure"
        if cc in ("Pure_Subtype", "Ancestral_Recombinant_Subtype") and n_bp == 0:
            return "Stanford_Pure_vs_Different_Pure_Trunk"
        if n_bp == 0 and cc not in ("Certified_CRF", "Novel_Circulating_CRF", "URF"):
            if rb in allowed_trunks:
                return "Concordant_Pure"
            return "Stanford_Pure_vs_Different_Pure_Trunk"
        return "Stanford_Pure_vs_RhizType_Recombinant"

    # 4. Named CRF in Stanford HIVDB — check Tier-2 called_lineage first (e.g., CRF01_AE, CRF02_AG, CRF07_BC)
    st_num = _crf_num(st)
    cl_num = _crf_num(cl)
    if cl == st or (st_num and cl_num and st_num == cl_num):
        return "Concordant_CRF_or_Recomb"

    crf_short = st.replace("CRF", "")
    if n_bp == 0 and cc in ("Pure_Subtype", "Ancestral_Recombinant_Subtype"):
        if st == "CRF01_AE" and (cl in ("A1", "A6") or rb == "A1"):
            return "Stanford_CRF01_AE_Pol_Is_Pure_A1_Trunk"
        return "Stanford_CRF_vs_RhizType_Pure_Trunk"

    if n_bp == 0 and cl == "":
        if st == "CRF01_AE" and rb == "A1":
            return "Stanford_CRF01_AE_Pol_Is_Pure_A1_Trunk"
        return "Stanford_CRF_vs_RhizType_Pure_Trunk"

    if (crf_short in rm or st in rm) and cc != "Pure_Subtype":
        return "Concordant_CRF_or_Recomb"
    if st == "CRF02_AG" and ("G0" in rf and ("A1" in rf or "A6" in rf)) and cl in ("CRF02_AG", "URF", ""):
        return "Concordant_CRF_or_Recomb"
    if st in ("CRF07_BC", "CRF08_BC", "CRF57_BC", "CRF64_BC", "CRF85_BC", "CRF86_BC") and ("B" in rf and "C" in rf) and cc != "Pure_Subtype":
        return "Concordant_CRF_or_Recomb"
    if st == "CRF03_A6B" and ("A6" in rf and "B" in rf) and cc != "Pure_Subtype":
        return "Concordant_CRF_or_Recomb"
    if st in ("CRF12_BF", "CRF17_BF", "CRF28_BF", "CRF29_BF", "CRF38_BF", "CRF39_BF", "CRF40_BF", "CRF46_BF") and ("B" in rf and "F1" in rf) and cc != "Pure_Subtype":
        return "Concordant_CRF_or_Recomb"
    if st in ("CRF06_cpx",) and ("G0" in rf or "J" in rf or "K" in rf) and cc != "Pure_Subtype":
        return "Concordant_CRF_or_Recomb"
    return "Stanford_CRF_vs_Different_RhizType_Architecture"


def main():
    print("=" * 96)
    print("STEP 15: RESCANNING STANFORD HIV DRUG RESISTANCE DATABASE (HIVDB) WITH RHIZTYPE")
    print("=" * 96)

    # -------------------------------------------------------------------------
    # PART 1: Audit Stanford HIVDB Genotyper Reference Panel (N = 720)
    # -------------------------------------------------------------------------
    ref_fasta = DATA_DIR / "genotype-references.fasta"
    ref_out_csv = RESULTS_DIR / "stanford_hivdb_reference_panel_rhiztype.csv"
    t0 = time.time()
    subprocess.run(
        [str(RHIZTYPE_BIN), str(ref_fasta), "--quiet", "--out", str(ref_out_csv)],
        check=True,
    )
    ref_sec = time.time() - t0
    ref_df = pd.read_csv(ref_out_csv)

    parts = ref_df["sequence_id"].str.split("|", expand=True)
    ref_df["accession"] = parts[0]
    ref_df["stanford_ref_subtype"] = parts[1]
    ref_df["stanford_span"] = parts[2]
    ref_df["Audit_Category"] = [
        classify_stanford_vs_rhiztype(st, 1, rf, rb, rm, nbp, cl, cc)
        for st, rf, rb, rm, nbp, cl, cc in zip(
            ref_df["stanford_ref_subtype"],
            ref_df["structural_formula"],
            ref_df["backbone_trunk"],
            ref_df["compendium_matches"],
            ref_df["n_breakpoints"],
            ref_df["called_lineage"],
            ref_df["called_category"],
        )
    ]
    ref_df.to_csv(ref_out_csv, index=False)

    print(f"\n[Part 1] Screened Stanford HIVDB Reference Panel ({len(ref_df):,} sequences) in {ref_sec:.3f}s:")
    for cat, cnt in ref_df["Audit_Category"].value_counts().items():
        print(f"  * {cat:<58s}: {cnt:4d} ({100*cnt/len(ref_df):5.1f}%)")

    # -------------------------------------------------------------------------
    # PART 2: Audit All 120,468 Clinical Surveillance Sequences in Stanford HIV-SurvDB
    # -------------------------------------------------------------------------
    surv_fasta = DATA_DIR / "stanford_survdb_all_120k.fasta"
    surv_rhiz_csv = RESULTS_DIR / "stanford_survdb_120k_rhiztype_raw.csv"
    t1 = time.time()
    subprocess.run(
        [str(RHIZTYPE_BIN), str(surv_fasta), "--quiet", "--out", str(surv_rhiz_csv)],
        check=True,
    )
    print(f"\n[Part 2] Screened {surv_fasta.name} with EGGMAN in {time.time() - t1:.2f}s.")
    surv_df = pd.read_csv(surv_rhiz_csv, low_memory=False)
    surv_df["patient_id"] = surv_df["sequence_id"].str.extract(r"PtID(\d+)")[0].astype(int)
    surv_df["study_file"] = surv_df["sequence_id"].str.split("|").str[0]

    csv_files = sorted((SURV_DIR / "tables" / "isolates.d").glob("*.csv"))
    iso_all = pd.concat([pd.read_csv(f, low_memory=False) for f in csv_files], ignore_index=True)

    pt_grouped = (
        iso_all.groupby("patient_id")
        .agg(
            stanford_subtypes_set=("subtype", lambda s: "|".join(sorted(set(str(x).strip() for x in s if pd.notna(x))))),
            n_distinct_stanford=("subtype", lambda s: len(set(str(x).strip() for x in s if pd.notna(x)))),
            genes_present=("gene", lambda s: "+".join(sorted(set(str(x).strip() for x in s if pd.notna(x))))),
            country_code=("country_code", "first"),
            isolate_year=("isolate_year", "first"),
            dataset_name=("dataset_name", "first"),
            genbank_accn=("genbank_accn", "first"),
        )
        .reset_index()
    )

    # Exact contradiction denominators across isolates.d-linked patients vs unlinked FASTA sequences
    n_iso_patients_total = len(pt_grouped)
    n_iso_patients_contra = int((pt_grouped["n_distinct_stanford"] > 1).sum())
    merged = surv_df.merge(pt_grouped, on="patient_id", how="left")
    merged["in_isolates_table"] = merged["stanford_subtypes_set"].notna()
    merged["stanford_subtypes_set"] = merged["stanford_subtypes_set"].fillna("Unlinked_No_Isolate_Row")
    merged["n_distinct_stanford"] = merged["n_distinct_stanford"].fillna(0).astype(int)
    merged["country_code"] = merged["country_code"].fillna("Unlinked")
    merged["genes_present"] = merged["genes_present"].fillna("Unlinked")

    linked_df = merged[merged["in_isolates_table"]]
    unlinked_df = merged[~merged["in_isolates_table"]]
    n_linked_patients = int(linked_df["patient_id"].nunique())
    n_fasta_patients_contra = int(linked_df.loc[linked_df["n_distinct_stanford"] > 1, "patient_id"].nunique())
    n_fasta_seqs_contra = int((linked_df["n_distinct_stanford"] > 1).sum())
    print(
        f"\n[Denominator Reconciliation] isolates.d patients: {n_iso_patients_contra:,} / {n_iso_patients_total:,} "
        f"({100*n_iso_patients_contra/n_iso_patients_total:.2f}%) have >1 Stanford label; "
        f"Linked FASTA cohort: {len(linked_df):,} seqs across {n_linked_patients:,} patients "
        f"({n_fasta_patients_contra:,} patients / {n_fasta_seqs_contra:,} seqs >1 label); "
        f"Unlinked FASTA absent from isolates.d: {len(unlinked_df):,} seqs ({unlinked_df['patient_id'].nunique():,} patients)."
    )

    merged["Audit_Category"] = [
        classify_stanford_vs_rhiztype(st, nd, rf, rb, rm, nbp, cl, cc)
        for st, nd, rf, rb, rm, nbp, cl, cc in zip(
            merged["stanford_subtypes_set"],
            merged["n_distinct_stanford"],
            merged["structural_formula"],
            merged["backbone_trunk"],
            merged["compendium_matches"],
            merged["n_breakpoints"],
            merged["called_lineage"],
            merged["called_category"],
        )
    ]

    print(f"\n[Part 2A] Stanford HIV-SurvDB Audit Summary across isolates.d-linked Cohort (N = {len(merged[merged['in_isolates_table']]):,}):")
    linked_merged = merged[merged["in_isolates_table"]]
    linked_counts = linked_merged["Audit_Category"].value_counts()
    for cat, cnt in linked_counts.items():
        print(f"  * {cat:<58s}: {cnt:7,d} ({100*cnt/len(linked_merged):5.2f}%)")

    print(f"\n[Part 2B] All 120,468 Clinical pol Sequences (Including {len(unlinked_df):,} Unlinked):")
    cat_counts = merged["Audit_Category"].value_counts()
    for cat, cnt in cat_counts.items():
        print(f"  * {cat:<58s}: {cnt:7,d} ({100*cnt/len(merged):5.2f}%)")

    cat_summary = pd.DataFrame({
        "Audit_Category": cat_counts.index,
        "Sequence_Count_All_120k": cat_counts.values,
        "Percentage_All_120k": np.round(100.0 * cat_counts.values / len(merged), 2),
        "Sequence_Count_Linked_111905": [int(linked_counts.get(c, 0)) for c in cat_counts.index],
        "Percentage_Linked_111905": [round(100.0 * linked_counts.get(c, 0) / len(linked_merged), 2) for c in cat_counts.index],
    })
    cat_summary.to_csv(RESULTS_DIR / "stanford_hivdb_audit_category_summary.csv", index=False)

    # Separate breakdown of the 8,563 unlinked FASTA sequences
    unlinked_merged = merged[~merged["in_isolates_table"]]
    unlinked_summary = (
        unlinked_merged.groupby(["called_category", "called_lineage", "study_file"])
        .size()
        .reset_index(name="n_sequences")
        .sort_values("n_sequences", ascending=False)
    )
    unlinked_summary.to_csv(RESULTS_DIR / "stanford_hivdb_unlinked_8k_breakdown.csv", index=False)

    # Tier-1 vs Tier-2 Breakdown for CRF01_AE and CRF02_AG (including G0[A6:RT])
    t12_rows = []
    for st_label in ("CRF01_AE", "CRF02_AG", "CRF07_BC", "CRF06_cpx", "C", "B"):
        sub = merged[(merged["n_distinct_stanford"] == 1) & (merged["stanford_subtypes_set"] == st_label)]
        top_t2 = sub["called_lineage"].value_counts().head(5).to_dict()
        top_t1 = sub["structural_formula"].value_counts().head(5).to_dict()
        t12_rows.append({
            "stanford_label": st_label,
            "n_seqs": len(sub),
            "top_tier2_called_lineage": str(top_t2),
            "top_tier1_structural_formula": str(top_t1),
        })
        print(f"  [Tier1 vs Tier2 | Stanford={st_label}, n={len(sub):,}] Tier-2: {top_t2} | Tier-1: {top_t1}")
    pd.DataFrame(t12_rows).to_csv(RESULTS_DIR / "stanford_hivdb_tier1_vs_tier2_breakdown.csv", index=False)

    # 1. Intra-patient PR vs RT vs IN contradictions in Stanford HIVDB (split locus-expected vs true clonal contradictions)
    contra_df = merged[merged["Audit_Category"].isin([
        "Stanford_MultiEnzyme_LocusExpected_Recombinant",
        "Stanford_IntraPatient_True_Clonal_Contradiction",
    ])]
    print(f"\n=== 1. Intra-Patient Gene-Row Contradictions in Stanford HIVDB (N = {len(contra_df):,}) ===")
    top_contra = (
        contra_df.groupby(["Audit_Category", "stanford_subtypes_set", "genes_present", "structural_formula", "country_code"])
        .size()
        .reset_index(name="n_sequences")
        .sort_values("n_sequences", ascending=False)
    )
    print(top_contra.head(20).to_string(index=False))
    top_contra.to_csv(RESULTS_DIR / "stanford_hivdb_intrapatient_contradictions.csv", index=False)

    # 2. Consistent Stanford Pure Subtype -> RhizType Intersubtype Recombinant
    pure_to_recomb = merged[merged["Audit_Category"] == "Stanford_Pure_vs_RhizType_Recombinant"]
    print(f"\n=== 2. Consistent Stanford 'Pure Subtype' Calls that are Recombinants in RhizType (N = {len(pure_to_recomb):,}) ===")
    top_pure_recomb = (
        pure_to_recomb.groupby(["stanford_subtypes_set", "structural_formula", "country_code"])
        .size()
        .reset_index(name="n_sequences")
        .sort_values("n_sequences", ascending=False)
    )
    print(top_pure_recomb.head(20).to_string(index=False))
    top_pure_recomb.to_csv(RESULTS_DIR / "stanford_hivdb_pure_to_recombinant_mismatches.csv", index=False)

    # 3. Consistent Stanford Pure Subtype -> Different Pure Trunk
    pure_to_diff = merged[merged["Audit_Category"] == "Stanford_Pure_vs_Different_Pure_Trunk"]
    print(f"\n=== 3. Consistent Stanford 'Pure Subtype' Calls with Different Pure Trunk (N = {len(pure_to_diff):,}) ===")
    top_pure_diff = (
        pure_to_diff.groupby(["stanford_subtypes_set", "structural_formula", "country_code"])
        .size()
        .reset_index(name="n_sequences")
        .sort_values("n_sequences", ascending=False)
    )
    print(top_pure_diff.head(20).to_string(index=False))
    top_pure_diff.to_csv(RESULTS_DIR / "stanford_hivdb_pure_trunk_mismatches.csv", index=False)

    # 4. Consistent Stanford Named CRF -> Different RhizType Architecture
    crf_diff = merged[merged["Audit_Category"] == "Stanford_CRF_vs_Different_RhizType_Architecture"]
    print(f"\n=== 4. Consistent Stanford Named CRF Calls with Different RhizType Architecture (N = {len(crf_diff):,}) ===")
    top_crf_diff = (
        crf_diff.groupby(["stanford_subtypes_set", "structural_formula", "country_code"])
        .size()
        .reset_index(name="n_sequences")
        .sort_values("n_sequences", ascending=False)
    )
    print(top_crf_diff.head(20).to_string(index=False))
    top_crf_diff.to_csv(RESULTS_DIR / "stanford_hivdb_crf_architecture_mismatches.csv", index=False)

    # 5. Stanford Ambiguous / Unknown / 'B + C' Resolved by RhizType
    ambig_resolved = merged[merged["Audit_Category"] == "Stanford_Ambiguous_or_Unknown_Resolved"]
    print(f"\n=== 5. Stanford 'Unknown' / 'B + C' / Ambiguous Calls Resolved by RhizType (N = {len(ambig_resolved):,}) ===")
    top_ambig = (
        ambig_resolved.groupby(["stanford_subtypes_set", "structural_formula", "country_code"])
        .size()
        .reset_index(name="n_sequences")
        .sort_values("n_sequences", ascending=False)
    )
    print(top_ambig.head(20).to_string(index=False))
    top_ambig.to_csv(RESULTS_DIR / "stanford_hivdb_ambiguous_resolved.csv", index=False)


if __name__ == "__main__":
    main()
