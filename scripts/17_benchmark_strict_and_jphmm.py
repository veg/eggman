#!/usr/bin/env python3
"""
Step 17: Strict Out-of-Sample Benchmarking & 3-Way Comparison (RhizType v0.2 vs COMET v2.4 vs jpHMM).

Addresses BARS Review Tier A1/A4 & Tier B4:
  1. --sync-comet: Syncs COMET v2.4 calls for any new/replaced sequences in
     data/benchmark_rhiztype_vs_comet.fasta (e.g., the 48 de-overlapped Panel 4 sequences),
     benchmarks RhizType v0.2 wall-clock speed, and updates
     results/benchmark_comet_v24_output.tsv, results/benchmark_rhiztype_vs_comet_merged.csv,
     and results/benchmark_rhiztype_v02_vs_comet_summary.csv.
  2. --loso: Runs genuine Leave-One-Sequence-Out (LOSO) CV across all 579 sequences in Panels 1 & 2:
     - 41 basal subtype genomes (A1, A2, A6, B, C, D, F1, F2, H, J, K, L): basal trunk N-1 consensus check.
     - 538 derived/CRF/novel genomes in certified_profiles_8628.tsv:
       * Standard LOSO (recomputing 8,628-nt consensus from N-1 remaining members; stratified N=2 vs N>=3).
       * Simulated Unseen/Singleton Drop-Out (omitting the lineage completely from --profiles TSV).
     Saves results/benchmark_loso_cv_summary.csv & results/benchmark_loso_cv_details.csv.
  3. --ood: Runs strict Out-of-Dictionary Holdout filtering certified_profiles_8628.tsv to Pure + CRF01-118
     (matching COMET v2.4 cutoff), comparing RhizType v0.2 vs COMET v2.4 on Panel 1C (CRF119-157, N=94)
     and Panel 2 (NOVEL_CRF_01-15, N=61). Saves results/benchmark_out_of_dictionary_holdout.csv.
  4. --parse-jphmm: Parses jpHMM output from m3.local (mapped via alignment_to_msa.txt + HXB2 MSA row),
     computes 3-way accuracy, parental composition, median breakpoint error on major founder CRFs,
     and CPU/wall speed across all 1,093 de-overlapped benchmark sequences.
     Saves results/benchmark_three_way_comparison.csv.
"""

import argparse
import concurrent.futures
import csv
import importlib.util
import json
import os
import re
import subprocess
import sys
import tempfile
import time
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

PROJECT_ROOT = Path("/Users/sergei/Projects/HIV-recombination")
DATA_DIR = PROJECT_ROOT / "data"
RESULTS_DIR = PROJECT_ROOT / "results"
RHIZ_DATA_DIR = PROJECT_ROOT / "rhiztype" / "crates" / "eggman-core" / "data"
RHIZTYPE_BIN = PROJECT_ROOT / "rhiztype" / "target" / "release" / "eggman"
BENCH_FASTA = DATA_DIR / "benchmark_rhiztype_vs_comet.fasta"
CERT_PROFILES_TSV = RHIZ_DATA_DIR / "certified_profiles_8628.tsv"
BASAL_PROFILES_TSV = RHIZ_DATA_DIR / "profiles_9719.tsv"
FOUNDER_PROFILES_TSV = RHIZ_DATA_DIR / "founder_profiles_9719.tsv"
COMPENDIUM_BUNDLE_JSON = RHIZ_DATA_DIR / "compendium_bundle.json"

PURE_TRUNKS_ALL = {"A1", "A2", "A3", "A4", "A6", "A7", "A8", "B", "C", "D", "F1", "F2", "G", "H", "J", "K", "L"}
BASAL_PANEL1A = {"A1", "A2", "A6", "B", "C", "D", "F1", "F2", "H", "J", "K", "L"}
DERIVED_PANEL1A = {"A3", "A4", "A7", "A8", "G"}

OUT_OF_PANEL_OR_COLLAPSED_CRFS = {"26_A5U", "92_C2U", "46_BF1"}
LEGACY_CRF04_SCRUB = {
    "JX500700": "CRF63_02A6",
    "MK984159": "CRF63_02A6",
    "PZ426393": "CRF63_02A6",
    "PQ585410": "CRF173_63A6",
    "PQ585412": "CRF173_63A6",
}

MAJOR_FOUNDER_CRFS = {
    "01_AE", "02_AG", "06_cpx", "07_BC", "08_BC", "11_cpx", "12_BF",
    "13_cpx", "15_01B", "33_01B", "55_01B", "63_02A6",
}


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, str(path))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


s14 = load_module("s14", PROJECT_ROOT / "scripts" / "14_benchmark_rhiztype_vs_comet.py")
s16 = load_module("s16", PROJECT_ROOT / "scripts" / "16_build_rhiztype_certified_profiles.py")


def norm_pure(lbl: str) -> str:
    s = str(lbl).strip()
    if s.startswith("Subtype_"):
        s = s[len("Subtype_"):]
    if s == "G0":
        s = "G"
    if s.startswith("A (check for A1)") or s == "A (A1)":
        s = "A1"
    return s


def extract_crf_num(lbl: str):
    s = str(lbl).strip()
    if s.startswith("CRF"):
        s = s[3:]
    m = re.match(r"^(\d+)_", s)
    return int(m.group(1)) if m else None


def norm_crf(lbl: str) -> str:
    s = str(lbl).strip()
    if s.startswith("CRF"):
        s = s[3:]
    return s


def crf_matches(gt_crf: str, pred_lbl: str) -> bool:
    g = norm_crf(gt_crf)
    p = norm_crf(pred_lbl)
    if not g or not p:
        return False
    if g == p:
        return True
    gn = extract_crf_num(g)
    pn = extract_crf_num(p)
    return gn is not None and pn is not None and gn == pn


def parse_fasta_metadata(fasta_path: Path):
    recs = s14.read_fasta(fasta_path)
    p2_audit = s14.load_panel2_audit_map()
    comp_df = pd.read_csv(RESULTS_DIR / "hiv1_descriptive_nomenclature_compendium.csv")
    crf_forms = {}
    for _, r in comp_df.iterrows():
        crf_forms[str(r["legacy_label"]).strip()] = str(r.get("rhizaeon_formula", "")).strip()

    meta = []
    for sid, raw in recs:
        parts = sid.split("|")
        expected_cat = "Pure_Subtype"
        is_active_mosaic_crf = False
        is_post2023_crf = False
        if parts[0] == "P1":
            hdr = parts[1]
            cl = hdr.split(".")[1]
            if cl in PURE_TRUNKS_ALL:
                panel = "Panel1A_LANL_Ref_PureSubtype"
                cert_lbl = f"Subtype_{norm_pure(cl)}"
                expected_cat = "Ancestral_Recombinant_Subtype" if cl == "G" else "Pure_Subtype"
            else:
                cnum = int(cl.split("_")[0]) if cl.split("_")[0].isdigit() else 999
                panel = (
                    "Panel1B_LANL_Ref_CRF_01_to_118 (In COMET v2.4 Training)"
                    if cnum <= 118
                    else "Panel1C_LANL_Ref_CRF_119_to_157 (Post-COMET v2.4 Cutoff)"
                )
                cert_lbl = f"CRF{cl}"
                expected_cat = "Certified_CRF"
                is_active_mosaic_crf = norm_crf(cl) not in OUT_OF_PANEL_OR_COLLAPSED_CRFS
            gt = cl
            orig = hdr
        elif parts[0] == "P2":
            lin_code, acc = parts[1], parts[2]
            _, disp, canon, form = p2_audit.get(acc, (lin_code, "Polyphyletic_Convergent_URF", "URF", ""))
            if not form and canon in crf_forms:
                form = crf_forms[canon]
            panel = "Panel2_Novel_FullGenome_CRFs (Audited 61)"
            if canon == "URF" or "Polyphyletic" in disp:
                expected_cat = "URF"
                cert_lbl = "URF"
                gt = f"URF:{form}" if form else "URF"
            elif canon == "RU_G0B_gagMA":
                expected_cat = "Novel_Circulating_CRF"
                cert_lbl = "RU_G0B_gagMA"
                gt = f"RU_G0B_gagMA:{form}" if form else "RU_G0B_gagMA"
            else:
                expected_cat = "Certified_CRF"
                cert_lbl = canon
                gt = f"{canon}:{form}" if form else canon
                cnum = extract_crf_num(canon)
                is_post2023_crf = cnum is not None and cnum > 118
            orig = acc
        elif parts[0] == "P3":
            lin_code, acc = parts[1], parts[2]
            if lin_code.startswith("POL_"):
                panel = "Panel3A_Novel_Partial_pol_Lineages"
            else:
                panel = "Panel3B_Novel_Partial_env_Lineages"
            gt = lin_code
            cert_lbl = lin_code
            orig = acc
            expected_cat = "URF"
        elif parts[0] == "P4":
            cat_key, hdr = parts[1], parts[2]
            if cat_key == "Pure":
                panel = "Panel4A_Clinical_PureSubtype"
                expected_cat = "Ancestral_Recombinant_Subtype" if hdr.split(".")[0] == "G" else "Pure_Subtype"
            elif cat_key == "CRF":
                panel = "Panel4B_Clinical_NamedCRF"
                expected_cat = "Certified_CRF"
            else:
                panel = "Panel4C_Clinical_URF"
                expected_cat = "URF"
            gt = hdr.split(".")[0]
            cert_lbl = gt
            orig = hdr
        else:
            raise ValueError(f"Unrecognized seq_id: {sid}")
        meta.append({
            "seq_id": sid,
            "raw_seq": s14.degap(raw),
            "panel": panel,
            "ground_truth_label": gt,
            "cert_lineage": cert_lbl,
            "expected_category": expected_cat,
            "is_active_mosaic_crf": is_active_mosaic_crf,
            "is_post2023_crf": is_post2023_crf,
            "original_id": orig,
        })
    return pd.DataFrame(meta)


def classify_comet_call(cs: str):
    """Return (category_tag, extracted_label) for a COMET v2.4 subtype string."""
    s = str(cs).strip()
    sl = s.lower()
    if not s or s == "nan" or "unassigned" in sl or "truncate" in sl or re.search(r"recombinan", sl):
        return "unassigned", s
    m_chk = re.search(r"\(check for ([^)]+)\)", s)
    if m_chk:
        target = m_chk.group(1).strip()
        if "_" in target and target.split("_")[0].isdigit():
            return "named_crf", target
        return "pure", norm_pure(target)
    if s in PURE_TRUNKS_ALL or s == "A":
        return "pure", norm_pure(s)
    if re.match(r"^\d+_", s) and "-" not in s and "," not in s and " " not in s:
        return "named_crf", s
    # Multi-subtype / mosaic string without 'unassigned'
    if "-" in s or "," in s:
        return "unassigned", s
    return "pure", norm_pure(s)


def extract_parental_trunks_from_formula(text: str):
    """Extract canonical HIV-1 basal subtypes present in a formula or CRF/COMET string."""
    s = str(text)
    if not s or s == "nan":
        return set()
    # Map common CRF aliases inside strings when evaluating parental composition
    CRF_COMPONENTS = {
        "01_AE": {"A1", "?E"}, "CRF01_AE": {"A1", "?E"},
        "02_AG": {"A1", "A6", "G"}, "CRF02_AG": {"A1", "A6", "G"},
        "06_cpx": {"A1", "G", "J", "K"}, "CRF06_cpx": {"A1", "G", "J", "K"},
        "07_BC": {"B", "C"}, "CRF07_BC": {"B", "C"},
        "08_BC": {"B", "C"}, "CRF08_BC": {"B", "C"},
        "11_cpx": {"A1", "G", "J", "H"}, "CRF11_cpx": {"A1", "G", "J", "H"},
        "12_BF": {"B", "F1"}, "CRF12_BF": {"B", "F1"},
        "15_01B": {"A1", "?E", "B"}, "CRF15_01B": {"A1", "?E", "B"},
        "33_01B": {"A1", "?E", "B"}, "CRF33_01B": {"A1", "?E", "B"},
        "55_01B": {"A1", "?E", "B"}, "CRF55_01B": {"A1", "?E", "B"},
        "63_02A6": {"A1", "A6", "G"}, "CRF63_02A6": {"A1", "A6", "G"},
    }
    tokens = set()
    for crf_tok, comps in CRF_COMPONENTS.items():
        if crf_tok in s:
            tokens.update(comps)
    # Parse CRF<num>_<suffix> tokens (e.g. 16_A2D, 72_BF1, 115_01C, 114_0155)
    for suf in re.findall(r"(?:CRF)?\d+_([A-Za-z0-9]+)", s):
        if suf == "cpx":
            continue
        if "01" in suf or "AE" in suf:
            tokens.update({"A1", "?E"})
        if "02" in suf or "AG" in suf:
            tokens.update({"A1", "G"})
        if "07" in suf or "08" in suf:
            tokens.update({"B", "C"})
        if "55" in suf:
            tokens.update({"A1", "?E", "B"})
        for sub_m in re.findall(r"(A[1-8]|F[12]|[ABCDGHJKL])", suf):
            if sub_m in ("A3", "A4", "A7", "A8"):
                tokens.add("A1")
            tokens.add(sub_m)
    # Match standalone/bracketed subtype tokens like A1..A8, F1/F2, G0, ?E, B, C, D, G, H, J, K, L
    # Preserve '?E' distinct from 'A1' so parental set equality & Jaccard reflect true ?E ancestry
    for m in re.findall(r"(?<![A-Za-z0-9])(?:Subtype_)?(\?E|A[1-8]|F[12]|G0|[ABCDGHJKL])(?![A-Za-z0-9_])", s):
        tok = "G" if m == "G0" else m
        if tok in ("A3", "A4", "A7", "A8"):
            tokens.add("A1")
        tokens.add(tok)
    return tokens


def extract_jphmm_subtypes_set(jp_subtypes):
    """Normalize jpHMM subtype tokens, keeping ?E distinct via 01_AE -> {'A1', '?E'}."""
    out = set()
    for x in jp_subtypes:
        sx = str(x).strip()
        if sx in ("01_AE", "AE", "CRF01_AE"):
            out.update({"A1", "?E"})
        elif sx == "?E":
            out.add("?E")
        else:
            out.add(norm_pure(sx))
    return out


def compute_parent_set_stats(gt_tr, call_tr):
    if not gt_tr or not call_tr:
        return 0, 0, 0.0
    ov = int(len(gt_tr & call_tr) >= 1)
    eq = int(gt_tr == call_tr)
    jac = len(gt_tr & call_tr) / len(gt_tr | call_tr)
    return ov, eq, jac


def evaluate_panel_summary(merged_df: pd.DataFrame) -> pd.DataFrame:
    panel_slices = [
        ("Panel1A_LANL_Ref_PureSubtype", merged_df[merged_df["panel"] == "Panel1A_LANL_Ref_PureSubtype"]),
        ("Panel1B_LANL_Ref_CRF_01_to_118 (In COMET v2.4 Training)", merged_df[merged_df["panel"] == "Panel1B_LANL_Ref_CRF_01_to_118 (In COMET v2.4 Training)"]),
        ("Panel1C_LANL_Ref_CRF_119_to_157 (Post-COMET v2.4 Cutoff)", merged_df[merged_df["panel"] == "Panel1C_LANL_Ref_CRF_119_to_157 (Post-COMET v2.4 Cutoff)"]),
        ("Panel1BC_All_145_LANL_CRFs", merged_df[merged_df["panel"].str.startswith(("Panel1B", "Panel1C"))]),
        ("Panel1BC_Active_142_Mosaic_CRFs", merged_df[merged_df["panel"].str.startswith(("Panel1B", "Panel1C")) & (merged_df["is_active_mosaic_crf"] == True)]),
        ("Panel2_Novel_FullGenome_CRFs (Audited 61)", merged_df[merged_df["panel"] == "Panel2_Novel_FullGenome_CRFs (Audited 61)"]),
        ("Panel2_Promoted_Official_CRFs (31: Post2023+CRF63)", merged_df[(merged_df["panel"] == "Panel2_Novel_FullGenome_CRFs (Audited 61)") & (merged_df["expected_category"] == "Certified_CRF")]),
        ("Panel2_Novel_Circulating_RU_G0B_gagMA (6)", merged_df[(merged_df["panel"] == "Panel2_Novel_FullGenome_CRFs (Audited 61)") & (merged_df["expected_category"] == "Novel_Circulating_CRF")]),
        ("Panel2_Polyphyletic_Convergent_URF (24)", merged_df[(merged_df["panel"] == "Panel2_Novel_FullGenome_CRFs (Audited 61)") & (merged_df["expected_category"] == "URF")]),
        ("Panel3A_Novel_Partial_pol_Lineages", merged_df[merged_df["panel"] == "Panel3A_Novel_Partial_pol_Lineages"]),
        ("Panel3B_Novel_Partial_env_Lineages", merged_df[merged_df["panel"] == "Panel3B_Novel_Partial_env_Lineages"]),
        ("Panel4A_Clinical_PureSubtype", merged_df[merged_df["panel"] == "Panel4A_Clinical_PureSubtype"]),
        ("Panel4B_Clinical_NamedCRF", merged_df[merged_df["panel"] == "Panel4B_Clinical_NamedCRF"]),
        ("Panel4C_Clinical_URF", merged_df[merged_df["panel"] == "Panel4C_Clinical_URF"]),
    ]
    rows = []
    for p, sub in panel_slices:
        n = len(sub)
        if n == 0:
            continue
        c_exact = c_pure = c_crf = c_unas = 0
        r_exact = r_pure = r_cert_crf = r_novel = r_urf = r_wrong = r_comp = r_recomb = r_clonal = 0
        r_par_eq = 0
        r_par_jacs = []
        for _, row in sub.iterrows():
            gt = str(row["ground_truth_label"])
            exp_cat = str(row.get("expected_category", ""))
            c_cat, c_lbl = classify_comet_call(row.get("comet_subtype", ""))
            if c_cat == "pure":
                c_pure += 1
            elif c_cat == "named_crf":
                c_crf += 1
            else:
                c_unas += 1

            r_cat = str(row.get("called_category", ""))
            r_lin = str(row.get("called_lineage", ""))
            r_bp = int(row.get("rhiztype_n_breakpoints", row.get("n_breakpoints", 0)))
            r_bas_form = str(row.get("basal_formula", ""))
            r_comps = str(row.get("compendium_matches", ""))

            if r_cat in ("Pure_Subtype", "Ancestral_Recombinant_Subtype"):
                r_pure += 1
            elif r_cat == "Certified_CRF":
                r_cert_crf += 1
            elif r_cat == "Novel_Circulating_CRF":
                r_novel += 1
            else:
                r_urf += 1

            if "(Clonal)" in r_bas_form or (r_bp == 0 and r_cat in ("Pure_Subtype", "Ancestral_Recombinant_Subtype")):
                r_clonal += 1
            else:
                r_recomb += 1

            gt_base = gt.split(":")[0]
            gt_form = gt.split(":", 1)[1] if ":" in gt else ""
            call_tr = extract_parental_trunks_from_formula(row.get("rhiztype_formula", row.get("structural_formula", "")))
            if p.startswith(("Panel1A", "Panel4A")):
                gt_tr = {norm_pure(gt_base)}
                gt_p = norm_pure(gt_base)
                if c_cat == "pure" and norm_pure(c_lbl) == gt_p:
                    c_exact += 1
                if r_cat in ("Pure_Subtype", "Ancestral_Recombinant_Subtype") and norm_pure(r_lin) == gt_p:
                    r_exact += 1
                elif r_cat in ("Certified_CRF", "Novel_Circulating_CRF"):
                    r_wrong += 1
                if (r_cat in ("Pure_Subtype", "Ancestral_Recombinant_Subtype") and norm_pure(r_lin) == gt_p) or f"Subtype_{gt_p}" in r_comps:
                    r_comp += 1
            elif p.startswith(("Panel1B", "Panel1C", "Panel4B")):
                gt_tr = extract_parental_trunks_from_formula(row.get("basal_formula", "")) or extract_parental_trunks_from_formula(gt_base)
                is_oob = p.startswith("Panel1C")
                if not is_oob and c_cat == "named_crf" and crf_matches(gt_base, c_lbl):
                    c_exact += 1
                if r_cat == "Certified_CRF" and crf_matches(gt_base, r_lin):
                    r_exact += 1
                elif r_cat in ("Pure_Subtype", "Certified_CRF", "Novel_Circulating_CRF"):
                    r_wrong += 1
                gt_tag = f"CRF{norm_crf(gt_base)}"
                gt_num = extract_crf_num(gt_base)
                in_comp = gt_tag in r_comps or (gt_num is not None and re.search(rf"CRF0*{gt_num}_", r_comps) is not None)
                if (r_cat == "Certified_CRF" and crf_matches(gt_base, r_lin)) or in_comp:
                    r_comp += 1
            elif p.startswith("Panel2"):
                gt_tr = extract_parental_trunks_from_formula(gt_form or r_bas_form or gt_base)
                if exp_cat == "Certified_CRF":
                    if r_cat == "Certified_CRF" and crf_matches(gt_base, r_lin):
                        r_exact += 1
                    elif r_cat in ("Pure_Subtype", "Certified_CRF", "Novel_Circulating_CRF"):
                        r_wrong += 1
                    gt_num = extract_crf_num(gt_base)
                    if (r_cat == "Certified_CRF" and crf_matches(gt_base, r_lin)) or (gt_num is not None and re.search(rf"CRF0*{gt_num}_", r_comps) is not None):
                        r_comp += 1
                    if c_cat == "named_crf" and crf_matches(gt_base, c_lbl):
                        c_exact += 1
                elif exp_cat == "Novel_Circulating_CRF":
                    if r_cat == "Novel_Circulating_CRF" and r_lin == gt_base:
                        r_exact += 1
                    elif r_cat in ("Pure_Subtype", "Certified_CRF", "Novel_Circulating_CRF"):
                        r_wrong += 1
                    if (r_cat == "Novel_Circulating_CRF" and r_lin == gt_base) or gt_base in r_comps:
                        r_comp += 1
                else:
                    # Audited polyphyletic URF: should remain URF (not forced into CRF)
                    if r_cat == "URF":
                        r_exact += 1
                        r_comp += 1
                    elif r_cat in ("Pure_Subtype", "Certified_CRF", "Novel_Circulating_CRF"):
                        r_wrong += 1
                    if c_cat == "unassigned":
                        c_exact += 1
            elif p.startswith("Panel3"):
                gt_tr = extract_parental_trunks_from_formula(r_bas_form or gt_base)
                if r_cat in ("Pure_Subtype", "Certified_CRF"):
                    r_wrong += 1
                if gt_tr and call_tr and (gt_tr & call_tr):
                    r_comp += 1
            else:
                # Panel4C_Clinical_URF
                gt_tr = extract_parental_trunks_from_formula(gt_base)
                if r_cat in ("Pure_Subtype", "Certified_CRF", "Novel_Circulating_CRF"):
                    r_wrong += 1
                if (r_cat == "URF" and len(call_tr) >= 2) or (gt_tr and call_tr and len(gt_tr & call_tr) >= 1):
                    r_comp += 1

            _, eq, jac = compute_parent_set_stats(gt_tr, call_tr)
            r_par_eq += eq
            r_par_jacs.append(jac)

        rows.append({
            "panel": p,
            "n_sequences": n,
            "comet_exact_match": c_exact,
            "comet_exact_pct": round(100.0 * c_exact / n, 1),
            "comet_pure_calls": c_pure,
            "comet_named_crf_calls": c_crf,
            "comet_unassigned_calls": c_unas,
            "comet_unassigned_pct": round(100.0 * c_unas / n, 1),
            "rhiztype_exact_lineage_call": r_exact,
            "rhiztype_exact_call_pct": round(100.0 * r_exact / n, 1),
            "rhiztype_pure_calls": r_pure,
            "rhiztype_certified_crf_calls": r_cert_crf,
            "rhiztype_novel_crf_calls": r_novel,
            "rhiztype_urf_calls": r_urf,
            "rhiztype_urf_pct": round(100.0 * r_urf / n, 1),
            "rhiztype_wrong_crf_calls": r_wrong,
            "rhiztype_compendium_or_call_match": r_comp,
            "rhiztype_compendium_or_call_pct": round(100.0 * r_comp / n, 1),
            "rhiztype_parent_exact_set_pct": round(100.0 * r_par_eq / n, 1),
            "rhiztype_parent_mean_jaccard": round(float(np.mean(r_par_jacs)), 3) if r_par_jacs else 0.0,
            "rhiztype_basal_recomb_resolved": r_recomb,
            "rhiztype_basal_recomb_pct": round(100.0 * r_recomb / n, 1),
            "rhiztype_basal_clonal_calls": r_clonal,
        })
    return pd.DataFrame(rows)


def cmd_sync_comet():
    meta_df = parse_fasta_metadata(BENCH_FASTA)
    current_ids = set(meta_df["seq_id"])
    comet_tsv = RESULTS_DIR / "benchmark_comet_v24_output.tsv"
    if comet_tsv.exists():
        comet_df = pd.read_csv(comet_tsv, sep="\t")
        comet_df = comet_df[comet_df["name"].astype(str).isin(current_ids)].drop_duplicates(subset=["name"])
    else:
        comet_df = pd.DataFrame(columns=["name", "virus", "subtype", "bootstrap support"])

    done_ids = set(comet_df["name"].astype(str))
    missing = [(r["seq_id"], r["raw_seq"]) for _, r in meta_df.iterrows() if r["seq_id"] not in done_ids]
    print(f"[*] COMET sync: {len(done_ids)} kept, {len(missing)} missing to submit.")
    if missing:
        new_batches = [comet_df]
        for b_start in range(0, len(missing), 600):
            batch = missing[b_start:b_start + 600]
            b_df = s14.submit_batch_to_comet(batch, batch_idx=90 + b_start // 600)
            new_batches.append(b_df)
        comet_df = pd.concat(new_batches, ignore_index=True).drop_duplicates(subset=["name"])
    comet_df = comet_df[comet_df["name"].astype(str).isin(current_ids)].reset_index(drop=True)
    comet_df.to_csv(comet_tsv, sep="\t", index=False)
    print(f"[+] Updated {comet_tsv} ({len(comet_df)} rows; stale purged).")

    # Benchmark RhizType v0.2 locally
    rhiz_csv = RESULTS_DIR / "benchmark_rhiztype_v02_output.csv"
    t0 = time.time()
    subprocess.run([str(RHIZTYPE_BIN), str(BENCH_FASTA), "--quiet", "--out", str(rhiz_csv)], check=True)
    local_wall = time.time() - t0
    print(f"[+] RhizType v0.2 ran on {len(meta_df)} unaligned seqs locally in {local_wall:.3f}s ({len(meta_df)/local_wall:.1f} seq/s).")

    # Check m3.local timing as well
    m3_cmd = (
        "ssh -o BatchMode=yes m3.local 'python3 -c \"import time, subprocess; "
        "t=time.time(); subprocess.run([\\\"/Users/sergei/hiv_recomb_bench/bin/eggman\\\", "
        "\\\"/Users/sergei/hiv_recomb_bench/data/benchmark_rhiztype_vs_comet.fasta\\\", \\\"--quiet\\\", "
        "\\\"--out\\\", \\\"/tmp/rt_m3_bench.csv\\\"], check=True); "
        "print(round(time.time()-t, 4))\"'"
    )
    m3_res = subprocess.run(m3_cmd, shell=True, capture_output=True, text=True)
    m3_wall = float(m3_res.stdout.strip().splitlines()[-1]) if m3_res.returncode == 0 and m3_res.stdout.strip() else None

    rhiz_df = pd.read_csv(rhiz_csv)
    merged = meta_df.drop(columns=["raw_seq"]).merge(
        comet_df.rename(columns={"name": "seq_id", "subtype": "comet_subtype", "bootstrap support": "comet_support"}),
        on="seq_id", how="left",
    ).merge(
        rhiz_df.rename(columns={
            "sequence_id": "seq_id",
            "structural_formula": "rhiztype_formula",
            "descriptive_name": "rhiztype_descriptive",
            "backbone_trunk": "rhiztype_backbone",
            "n_breakpoints": "rhiztype_n_breakpoints",
            "hxb2_barcode": "rhiztype_barcode",
        }),
        on="seq_id", how="left",
    )
    merged.to_csv(RESULTS_DIR / "benchmark_rhiztype_vs_comet_merged.csv", index=False)

    sum_df = evaluate_panel_summary(merged)
    sum_df.to_csv(RESULTS_DIR / "benchmark_rhiztype_v02_vs_comet_summary.csv", index=False)
    print(sum_df[["panel", "n_sequences", "comet_exact_pct", "rhiztype_exact_call_pct", "rhiztype_urf_pct", "rhiztype_wrong_crf_calls"]].to_string(index=False))

    rt_meta = {
        "n_sequences": int(len(meta_df)),
        "rhiztype_local_wall_sec": round(local_wall, 4),
        "rhiztype_local_seq_per_sec": round(len(meta_df) / local_wall, 1),
        "rhiztype_m3_wall_sec": m3_wall,
        "rhiztype_m3_seq_per_sec": round(len(meta_df) / m3_wall, 1) if m3_wall else None,
    }
    with open(RESULTS_DIR / "benchmark_runtime_comparison.json", "w") as f:
        json.dump(rt_meta, f, indent=2)
    print("[+] Runtime meta:", rt_meta)


def load_certified_tsv_rows(path: Path):
    rows = []
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            rows.append(parts)
    return rows


def build_p1_p2_aligned_arrays(meta_df: pd.DataFrame):
    """Return {seq_id: int8_array_8628} for all 579 Panel 1+2 sequences + 9719 basal array for Panel 1A."""
    lanl_recs = s14.read_fasta(DATA_DIR / "lanl_ref_2023.fasta")
    hxb2_seq = [s for h, s in lanl_recs if "HXB2" in h][0]
    cols_9719 = [i for i, ch in enumerate(hxb2_seq) if ch != "-"]
    cols_8628 = cols_9719[789:9417]
    assert len(cols_9719) == 9719 and len(cols_8628) == 8628

    lanl_map = {h: s for h, s in lanl_recs}
    arrays_8628 = {}
    arrays_9719 = {}
    novel_pairs = []

    for _, r in meta_df.iterrows():
        sid = r["seq_id"]
        if sid.startswith("P1|"):
            orig = r["original_id"]
            aln_s = lanl_map[orig]
            s8628 = "".join(aln_s[i] for i in cols_8628)
            s9719 = "".join(aln_s[i] for i in cols_9719)
            arrays_8628[sid] = s16.seq_to_int_array(s8628)
            arrays_9719[sid] = s16.seq_to_int_array(s9719)
        elif sid.startswith("P2|"):
            novel_pairs.append((sid, r["raw_seq"]))

    print(f"[*] Aligning {len(novel_pairs)} Panel-2 novel CRF genomes to 8,628-nt HXB2 coding frame via mafft --addfragments...")
    novel_aligned = s16.align_novel_to_hxb2_coding(novel_pairs, DATA_DIR / "hxb2_15trunk_profile_9719.fasta")
    arrays_8628.update(novel_aligned)
    return arrays_8628, arrays_9719


def run_single_rhiztype(seq_id: str, raw_seq: str, prof_tsv_text: str):
    with tempfile.NamedTemporaryFile("w", suffix=".fasta", delete=False) as tf_fa, \
         tempfile.NamedTemporaryFile("w", suffix=".tsv", delete=False) as tf_pr, \
         tempfile.NamedTemporaryFile("w", suffix=".csv", delete=False) as tf_out:
        fa_p, pr_p, out_p = tf_fa.name, tf_pr.name, tf_out.name
        tf_fa.write(f">{seq_id}\n{raw_seq}\n")
        tf_pr.write(prof_tsv_text)
    try:
        subprocess.run(
            [str(RHIZTYPE_BIN), fa_p, "--profiles", pr_p, "--quiet", "--out", out_p],
            check=True,
        )
        df = pd.read_csv(out_p)
        return df.iloc[0].to_dict()
    finally:
        for p in (fa_p, pr_p, out_p):
            if os.path.exists(p):
                os.unlink(p)


def write_tsv_text(rows) -> str:
    lines = ["# label\tcategory\tflags\tdag_parents\tn_seqs\tseq_8628"]
    for parts in rows:
        lines.append("\t".join(parts))
    return "\n".join(lines) + "\n"


def eval_panel1a_basal_loo(p1a_basal_rows, arrays_9719):
    """Real eggman --pack LOO CV across the 12 pure basal trunks in profiles_9719.tsv."""
    import shutil
    base_rows = []
    with open(BASAL_PROFILES_TSV) as f:
        for line in f:
            s = line.rstrip("\n")
            if not s or s.startswith("#"):
                base_rows.append(s)
            else:
                parts = s.split("\t")
                base_rows.append(parts)

    cert_rows = load_certified_tsv_rows(CERT_PROFILES_TSV)
    by_gt = defaultdict(list)
    for _, r in p1a_basal_rows.iterrows():
        gt = str(r["ground_truth_label"])
        by_gt[gt].append(r["seq_id"])

    with tempfile.NamedTemporaryFile("w", suffix=".json", delete=False) as tf_m:
        manifest_path = Path(tf_m.name)
    subprocess.run([str(RHIZTYPE_BIN), "--export-default-pack", str(manifest_path)], check=True, stdout=subprocess.DEVNULL)
    manifest_txt = manifest_path.read_text()
    manifest_path.unlink(missing_ok=True)
    founder_txt = FOUNDER_PROFILES_TSV.read_text()
    comp_txt = COMPENDIUM_BUNDLE_JSON.read_text()

    def _run_basal_item(row_dict):
        sid = row_dict["seq_id"]
        gt = str(row_dict["ground_truth_label"])
        raw_seq = row_dict["raw_seq"]
        rem = [arrays_9719[x] for x in by_gt[gt] if x != sid]
        loo_cons_str = s16.plurality_consensus_str(rem)
        gt_norm = norm_pure(gt)

        with tempfile.TemporaryDirectory() as td:
            pdir = Path(td)
            (pdir / "manifest.json").write_text(manifest_txt)
            (pdir / "founder_profiles.tsv").write_text(founder_txt)
            (pdir / "compendium.json").write_text(comp_txt)
            filtered_cert = [cr for cr in cert_rows if cr[0] != f"Subtype_{gt_norm}"]
            (pdir / "certified_profiles.tsv").write_text(write_tsv_text(filtered_cert))

            prof_lines = []
            for item in base_rows:
                if isinstance(item, str):
                    prof_lines.append(item)
                elif item[0] == gt:
                    prof_lines.append("\t".join([item[0], loo_cons_str] + item[2:]))
                else:
                    prof_lines.append("\t".join(item))
            (pdir / "profiles.tsv").write_text("\n".join(prof_lines) + "\n")

            fa_p = pdir / "q.fa"
            out_p = pdir / "out.csv"
            fa_p.write_text(f">{sid}\n{raw_seq}\n")
            subprocess.run([str(RHIZTYPE_BIN), "--pack", str(pdir), str(fa_p), "--quiet", "--out", str(out_p)], check=True)
            res = pd.read_csv(out_p).iloc[0].to_dict()

        call_lin = str(res["called_lineage"])
        call_cat = str(res["called_category"])
        struct_f = str(res["structural_formula"])
        exact = int(call_cat in ("Pure_Subtype", "Ancestral_Recombinant_Subtype") and norm_pure(call_lin) == gt_norm)
        call_tr = extract_parental_trunks_from_formula(struct_f) | {norm_pure(call_lin)}
        return {
            "seq_id": sid,
            "panel": row_dict["panel"],
            "ground_truth_lineage": gt,
            "mode": "Basal_LOO",
            "n_ref_members": len(by_gt[gt]),
            "n_remaining_ref": len(rem),
            "called_lineage": call_lin,
            "called_category": call_cat,
            "exact_match": exact,
            "structural_formula": struct_f,
            "parent_recall_ok": int(gt_norm in call_tr),
            "nearest_pdist": np.nan,
            "margin_pdist": np.nan,
        }

    records = p1a_basal_rows.to_dict("records")
    with concurrent.futures.ThreadPoolExecutor(max_workers=8) as ex:
        results = list(ex.map(_run_basal_item, records))
    return results


def cmd_loso():
    meta_df = parse_fasta_metadata(BENCH_FASTA)
    p12_meta = meta_df[meta_df["seq_id"].str.startswith(("P1|", "P2|"))].copy()
    print(f"[*] Starting LOSO CV across Panels 1A-1C + Panel 2 ({len(p12_meta)} sequences)...")

    arrays_8628, arrays_9719 = build_p1_p2_aligned_arrays(p12_meta)
    cert_rows = load_certified_tsv_rows(CERT_PROFILES_TSV)
    cert_by_lbl = {r[0]: list(r) for r in cert_rows}

    p1a_basal = p12_meta[(p12_meta["panel"] == "Panel1A_LANL_Ref_PureSubtype") & (p12_meta["ground_truth_label"].isin(BASAL_PANEL1A))].copy()
    cert_eval_df = p12_meta[p12_meta["cert_lineage"].isin(cert_by_lbl)].copy()
    print(f"    * Basal 1A LOO seqs: {len(p1a_basal)} | Certified 8628 LOO seqs: {len(cert_eval_df)} (sum={len(p1a_basal)+len(cert_eval_df)})")

    # Group arrays by certified lineage
    members_by_lin = defaultdict(list)
    for _, r in cert_eval_df.iterrows():
        members_by_lin[r["cert_lineage"]].append(r["seq_id"])

    # Prepare tasks for both Standard LOSO (N-1 recomputed) & Unseen/Dropped-Lineage (N=1 mode)
    tasks = []
    for _, r in cert_eval_df.iterrows():
        sid = r["seq_id"]
        lin = r["cert_lineage"]
        mems = members_by_lin[lin]
        rem = [arrays_8628[m] for m in mems if m != sid]

        # (1) Standard LOSO: recompute consensus from N-1 remaining members when >= 1 remains
        if rem:
            new_cons = s16.plurality_consensus_str(rem)
            loo_rows = []
            for row in cert_rows:
                if row[0] == lin:
                    copy_r = list(row)
                    copy_r[4] = str(len(rem))
                    copy_r[5] = new_cons
                    loo_rows.append(copy_r)
                else:
                    loo_rows.append(row)
            tasks.append(("LOSO_N_minus_1", len(mems), len(rem), sid, r["raw_seq"], r["panel"], r["ground_truth_label"], lin, r["expected_category"], write_tsv_text(loo_rows)))

        # (2) Simulated Unseen/Singleton Drop-Out: omit lin completely
        drop_rows = [row for row in cert_rows if row[0] != lin]
        tasks.append(("Singleton_Omitted", len(mems), 0, sid, r["raw_seq"], r["panel"], r["ground_truth_label"], lin, r["expected_category"], write_tsv_text(drop_rows)))

    detail_rows = []
    # Run Real eggman --pack Basal 1A LOO
    basal_res = eval_panel1a_basal_loo(p1a_basal, arrays_9719)
    detail_rows.extend(basal_res)

    full_rt = pd.read_csv(RESULTS_DIR / "benchmark_rhiztype_v02_output.csv")
    full_basal_map = dict(zip(full_rt["sequence_id"], full_rt["basal_formula"]))

    print(f"[*] Dispatching {len(tasks)} EGGMAN --profiles evaluations across 8 workers...")
    t0 = time.time()
    def _worker(t):
        mode, n_ref, n_rem, sid, raw_seq, panel, gt, lin, exp_cat, tsv_txt = t
        res = run_single_rhiztype(sid, raw_seq, tsv_txt)
        call_lin = str(res["called_lineage"])
        call_cat = str(res["called_category"])
        struct_f = str(res["structural_formula"])
        if panel == "Panel1A_LANL_Ref_PureSubtype":
            exact = int(call_cat in ("Pure_Subtype", "Ancestral_Recombinant_Subtype") and norm_pure(call_lin) == norm_pure(gt))
        elif exp_cat == "Novel_Circulating_CRF":
            exact = int(call_cat == "Novel_Circulating_CRF" and call_lin == lin)
        else:
            exact = int(call_cat == "Certified_CRF" and crf_matches(gt.split(":")[0], call_lin))

        gt_form = gt.split(":", 1)[1] if ":" in gt else lin
        gt_tr = extract_parental_trunks_from_formula(full_basal_map.get(sid, "") or gt_form)
        call_tr = (
            extract_parental_trunks_from_formula(struct_f)
            | extract_parental_trunks_from_formula(str(res.get("basal_formula", "")))
            | extract_parental_trunks_from_formula(call_lin)
        )
        parent_ok = int(len(gt_tr & call_tr) >= 1) if gt_tr and call_tr else int(exact == 1)

        return {
            "seq_id": sid,
            "panel": panel,
            "ground_truth_lineage": lin,
            "mode": mode,
            "n_ref_members": n_ref,
            "n_remaining_ref": n_rem,
            "called_lineage": call_lin,
            "called_category": call_cat,
            "exact_match": exact,
            "structural_formula": struct_f,
            "parent_recall_ok": parent_ok,
            "nearest_pdist": np.nan,
            "margin_pdist": np.nan,
        }

    with concurrent.futures.ThreadPoolExecutor(max_workers=8) as ex:
        for out in ex.map(_worker, tasks):
            detail_rows.append(out)
    print(f"[+] Completed {len(tasks)} --profiles evaluations in {time.time()-t0:.2f}s.")

    det_df = pd.DataFrame(detail_rows)
    det_df.to_csv(RESULTS_DIR / "benchmark_loso_cv_details.csv", index=False)

    # Summarize cleanly
    sum_rows = []
    # Basal 1A
    b_sub = det_df[det_df["mode"] == "Basal_LOO"]
    sum_rows.append({
        "evaluation_stratum": "Panel1A_Basal_Subtypes_LOO (12 trunks, eggman --pack)",
        "mode": "Basal_LOO",
        "n_sequences": len(b_sub),
        "exact_or_basal_match": int(b_sub["exact_match"].sum()),
        "exact_or_basal_pct": round(100.0 * b_sub["exact_match"].mean(), 2),
        "urf_flagged": int((b_sub["called_category"] == "URF").sum()),
        "urf_flagged_pct": round(100.0 * (b_sub["called_category"] == "URF").mean(), 2),
        "false_legacy_crf_calls": int(((b_sub["exact_match"] == 0) & (b_sub["called_category"] != "URF")).sum()),
        "parent_composition_ok": int(b_sub["parent_recall_ok"].sum()),
        "parent_composition_pct": round(100.0 * b_sub["parent_recall_ok"].mean(), 2),
    })

    # Certified LOSO overall and stratified
    loso_df = det_df[det_df["mode"] == "LOSO_N_minus_1"]
    for label, sub in [
        ("Certified_LOSO_All (Derived+CRF+Novel)", loso_df),
        ("Certified_LOSO_N_eq_2 (N-1=1 single remaining template)", loso_df[loso_df["n_ref_members"] == 2]),
        ("Certified_LOSO_N_ge_3 (N-1>=2 consensus)", loso_df[loso_df["n_ref_members"] >= 3]),
        ("Certified_LOSO_Panel1A_Derived (G)", loso_df[loso_df["panel"] == "Panel1A_LANL_Ref_PureSubtype"]),
        ("Certified_LOSO_Panel1B (CRF01-118)", loso_df[loso_df["panel"].str.startswith("Panel1B")]),
        ("Certified_LOSO_Panel1C (CRF119-157)", loso_df[loso_df["panel"].str.startswith("Panel1C")]),
        ("Certified_LOSO_Panel2_Certified_And_Novel (37)", loso_df[loso_df["panel"].str.startswith("Panel2")]),
    ]:
        n = len(sub)
        if n == 0:
            continue
        urf_n = int((sub["called_category"] == "URF").sum())
        wrong_n = int(((sub["exact_match"] == 0) & (sub["called_category"] != "URF")).sum())
        sum_rows.append({
            "evaluation_stratum": label,
            "mode": "LOSO_N_minus_1",
            "n_sequences": n,
            "exact_or_basal_match": int(sub["exact_match"].sum()),
            "exact_or_basal_pct": round(100.0 * sub["exact_match"].mean(), 2),
            "urf_flagged": urf_n,
            "urf_flagged_pct": round(100.0 * urf_n / n, 2),
            "false_legacy_crf_calls": wrong_n,
            "parent_composition_ok": int(sub["parent_recall_ok"].sum()),
            "parent_composition_pct": round(100.0 * sub["parent_recall_ok"].mean(), 2),
        })

    # Combined 579-seq LOSO
    comb_579 = pd.concat([b_sub, loso_df], ignore_index=True)
    sum_rows.append({
        "evaluation_stratum": "Combined_579_LOSO (Panels 1A-1C & 2)",
        "mode": "Basal+LOSO",
        "n_sequences": len(comb_579),
        "exact_or_basal_match": int(comb_579["exact_match"].sum()),
        "exact_or_basal_pct": round(100.0 * comb_579["exact_match"].mean(), 2),
        "urf_flagged": int((comb_579["called_category"] == "URF").sum()),
        "urf_flagged_pct": round(100.0 * (comb_579["called_category"] == "URF").mean(), 2),
        "false_legacy_crf_calls": int(((comb_579["exact_match"] == 0) & (comb_579["called_category"] != "URF")).sum()),
        "parent_composition_ok": int(comb_579["parent_recall_ok"].sum()),
        "parent_composition_pct": round(100.0 * comb_579["parent_recall_ok"].mean(), 2),
    })

    # Unseen / Singleton drop-out mode
    drop_df = det_df[det_df["mode"] == "Singleton_Omitted"]
    for label, sub in [
        ("Singleton_Omitted_All (Lineage Dropped from Tier-2)", drop_df),
        ("Singleton_Omitted_Panel1C (Post-118 Unseen)", drop_df[drop_df["panel"].str.startswith("Panel1C")]),
        ("Singleton_Omitted_Panel2 (NOVEL_CRF_01-15 Unseen)", drop_df[drop_df["panel"].str.startswith("Panel2")]),
    ]:
        n = len(sub)
        if n == 0:
            continue
        urf_n = int((sub["called_category"] == "URF").sum())
        wrong_n = int((sub["called_category"] != "URF").sum())
        sum_rows.append({
            "evaluation_stratum": label,
            "mode": "Singleton_Omitted",
            "n_sequences": n,
            "exact_or_basal_match": urf_n,
            "exact_or_basal_pct": round(100.0 * urf_n / n, 2),
            "urf_flagged": urf_n,
            "urf_flagged_pct": round(100.0 * urf_n / n, 2),
            "false_legacy_crf_calls": wrong_n,
            "parent_composition_ok": int(sub["parent_recall_ok"].sum()),
            "parent_composition_pct": round(100.0 * sub["parent_recall_ok"].mean(), 2),
        })

    sum_df = pd.DataFrame(sum_rows)
    sum_df.to_csv(RESULTS_DIR / "benchmark_loso_cv_summary.csv", index=False)
    print(sum_df.to_string(index=False))


def cmd_ood():
    meta_df = parse_fasta_metadata(BENCH_FASTA)
    merged_comet = pd.read_csv(RESULTS_DIR / "benchmark_rhiztype_vs_comet_merged.csv")
    cert_rows = load_certified_tsv_rows(CERT_PROFILES_TSV)

    # Filter strictly to Pure/Ancestral + CRF01..CRF118 (including scrubbed CRF63_02A6)
    filtered_rows = []
    for row in cert_rows:
        lbl, cat = row[0], row[1]
        if cat in ("Pure_Subtype", "Ancestral_Recombinant_Subtype"):
            filtered_rows.append(row)
        elif cat == "Certified_CRF":
            cnum = extract_crf_num(lbl)
            if cnum is not None and cnum <= 118:
                filtered_rows.append(row)

    print(f"[*] OOD Holdout dictionary: {len(filtered_rows)}/{len(cert_rows)} profiles kept (Pure + CRF01-118 only).")
    with tempfile.NamedTemporaryFile("w", suffix=".tsv", delete=False) as tf_pr, \
         tempfile.NamedTemporaryFile("w", suffix=".csv", delete=False) as tf_out:
        pr_p, out_p = tf_pr.name, tf_out.name
        tf_pr.write(write_tsv_text(filtered_rows))
    try:
        subprocess.run([str(RHIZTYPE_BIN), str(BENCH_FASTA), "--profiles", pr_p, "--quiet", "--out", out_p], check=True)
        ood_df = pd.read_csv(out_p)
    finally:
        for p in (pr_p, out_p):
            if os.path.exists(p):
                os.unlink(p)

    # Build true basal parental set for CRF119-157 from DAG/compendium or full RhizType calls
    full_rt = pd.read_csv(RESULTS_DIR / "benchmark_rhiztype_v02_output.csv")
    full_basal_map = dict(zip(full_rt["sequence_id"], full_rt["basal_formula"]))

    joined = meta_df.merge(
        merged_comet[["seq_id", "comet_subtype"]], on="seq_id", how="left"
    ).merge(
        ood_df.rename(columns={"sequence_id": "seq_id"}), on="seq_id", how="left"
    )

    # Summarize out-of-dictionary holdout (strictly exclude in-dict CRF63_02A6 from Panel2 OOD)
    p1c = joined[joined["panel"] == "Panel1C_LANL_Ref_CRF_119_to_157 (Post-COMET v2.4 Cutoff)"]
    p2_ood = joined[(joined["panel"] == "Panel2_Novel_FullGenome_CRFs (Audited 61)") & (joined["cert_lineage"] != "CRF63_02A6")]
    p2_post23 = joined[(joined["panel"] == "Panel2_Novel_FullGenome_CRFs (Audited 61)") & (joined["is_post2023_crf"] == True)]
    p2_novel = joined[(joined["panel"] == "Panel2_Novel_FullGenome_CRFs (Audited 61)") & (joined["expected_category"] == "Novel_Circulating_CRF")]
    p2_poly = joined[(joined["panel"] == "Panel2_Novel_FullGenome_CRFs (Audited 61)") & (joined["expected_category"] == "URF")]
    p1a = joined[joined["panel"] == "Panel1A_LANL_Ref_PureSubtype"]
    p1b = joined[joined["panel"] == "Panel1B_LANL_Ref_CRF_01_to_118 (In COMET v2.4 Training)"]

    out_rows = []
    for p_label, sub, is_ood in [
        ("Panel1C_LANL_Ref_CRF_119_to_157 (Post-COMET v2.4 Cutoff)", p1c, True),
        ("Panel2_Post2023_And_Novel_OOD (34: Post2023 CRFs+RU_G0B_gagMA)", pd.concat([p2_post23, p2_novel], ignore_index=True), True),
        ("Panel2_Polyphyletic_Convergent_URF (24)", p2_poly, True),
        ("Panel2_All_OOD_58 (Excl 3 CRF63_02A6)", p2_ood, True),
        ("Panel1A_LANL_Ref_PureSubtype", p1a, False),
        ("Panel1B_LANL_Ref_CRF_01_to_118 (In COMET v2.4 Training)", p1b, False),
    ]:
        n = len(sub)
        if n == 0:
            continue

        c_false_legacy = c_chimera_legacy = c_urf = c_parent_ok = c_exact = 0
        r_false_legacy = r_chimera_legacy = r_urf = r_parent_ok = r_exact = 0

        for _, row in sub.iterrows():
            gt = str(row["ground_truth_label"])
            gt_base = gt.split(":")[0]
            if p_label.startswith("Panel2"):
                gt_tr = extract_parental_trunks_from_formula(full_basal_map.get(row["seq_id"], "") or (gt.split(":", 1)[1] if ":" in gt else ""))
            elif p_label.startswith(("Panel1B", "Panel1C")):
                gt_tr = extract_parental_trunks_from_formula(full_basal_map.get(row["seq_id"], "")) | extract_parental_trunks_from_formula(gt_base)
            else:
                gt_tr = extract_parental_trunks_from_formula(gt_base)

            # COMET
            c_raw = str(row["comet_subtype"])
            c_cat, c_lbl = classify_comet_call(c_raw)
            if is_ood:
                if c_cat in ("pure", "named_crf"):
                    c_false_legacy += 1
                    c_chimera_legacy += 1
                elif "unassigned_2" in c_raw or re.search(r"\d+_[A-Za-z0-9]+", c_raw):
                    c_chimera_legacy += 1
                if c_cat == "unassigned":
                    c_urf += 1
            else:
                if p_label.startswith("Panel1A") and c_cat == "pure" and norm_pure(c_lbl) == norm_pure(gt_base):
                    c_exact += 1
                elif p_label.startswith("Panel1B") and c_cat == "named_crf" and crf_matches(gt_base, c_lbl):
                    c_exact += 1
            c_tr = extract_parental_trunks_from_formula(c_raw)
            if (p_label.startswith("Panel1B") and c_cat == "named_crf" and crf_matches(gt_base, c_lbl)) or (gt_tr and c_tr and len(gt_tr & c_tr) >= 1):
                c_parent_ok += 1

            # EGGMAN (CRF01-118 dictionary only)
            r_cat = str(row["called_category"])
            r_lin = str(row["called_lineage"])
            if is_ood:
                if r_cat in ("Pure_Subtype", "Certified_CRF", "Ancestral_Recombinant_Subtype"):
                    r_false_legacy += 1
                    r_chimera_legacy += 1
                if r_cat in ("URF", "Novel_Circulating_CRF"):
                    r_urf += 1
            else:
                if p_label.startswith("Panel1A") and r_cat in ("Pure_Subtype", "Ancestral_Recombinant_Subtype") and norm_pure(r_lin) == norm_pure(gt_base):
                    r_exact += 1
                elif p_label.startswith("Panel1B") and r_cat == "Certified_CRF" and crf_matches(gt_base, r_lin):
                    r_exact += 1
            r_tr = (
                extract_parental_trunks_from_formula(row["structural_formula"])
                | extract_parental_trunks_from_formula(row.get("basal_formula", ""))
                | extract_parental_trunks_from_formula(r_lin)
            )
            if (p_label.startswith("Panel1B") and r_cat == "Certified_CRF" and crf_matches(gt_base, r_lin)) or (gt_tr and r_tr and len(gt_tr & r_tr) >= 1):
                r_parent_ok += 1

        for method, fl, cl, urf_c, p_ok, ex_c in [
            ("COMET_v2.4", c_false_legacy, c_chimera_legacy, c_urf, c_parent_ok, c_exact),
            ("EGGMAN_CRF01_118_Holdout", r_false_legacy, r_chimera_legacy, r_urf, r_parent_ok, r_exact),
        ]:
            out_rows.append({
                "panel": p_label,
                "is_out_of_dictionary": is_ood,
                "method": method,
                "n_sequences": n,
                "false_legacy_misassignment": fl if is_ood else 0,
                "false_legacy_misassignment_pct": round(100.0 * fl / n, 2) if is_ood else 0.0,
                "legacy_crf_chimera_or_misassigned": cl if is_ood else 0,
                "legacy_crf_chimera_or_misassigned_pct": round(100.0 * cl / n, 2) if is_ood else 0.0,
                "novel_urf_flagged": urf_c if is_ood else 0,
                "novel_urf_sensitivity_pct": round(100.0 * urf_c / n, 2) if is_ood else 0.0,
                "in_dict_exact_accuracy_pct": round(100.0 * ex_c / n, 2) if not is_ood else np.nan,
                "parent_subtype_composition_ok": p_ok,
                "parent_subtype_composition_pct": round(100.0 * p_ok / n, 2),
            })

    # Combined OOD row (Panel 1C + Panel 2 OOD 58 = 152)
    ood_sub = pd.DataFrame([r for r in out_rows if r["panel"] in ("Panel1C_LANL_Ref_CRF_119_to_157 (Post-COMET v2.4 Cutoff)", "Panel2_All_OOD_58 (Excl 3 CRF63_02A6)")])
    for method in ["COMET_v2.4", "EGGMAN_CRF01_118_Holdout"]:
        m_sub = ood_sub[ood_sub["method"] == method]
        n_tot = int(m_sub["n_sequences"].sum())
        fl_tot = int(m_sub["false_legacy_misassignment"].sum())
        cl_tot = int(m_sub["legacy_crf_chimera_or_misassigned"].sum())
        urf_tot = int(m_sub["novel_urf_flagged"].sum())
        p_tot = int(m_sub["parent_subtype_composition_ok"].sum())
        out_rows.append({
            "panel": "Combined_OOD_Holdout (Panel1C+Panel2_OOD58)",
            "is_out_of_dictionary": True,
            "method": method,
            "n_sequences": n_tot,
            "false_legacy_misassignment": fl_tot,
            "false_legacy_misassignment_pct": round(100.0 * fl_tot / n_tot, 2),
            "legacy_crf_chimera_or_misassigned": cl_tot,
            "legacy_crf_chimera_or_misassigned_pct": round(100.0 * cl_tot / n_tot, 2),
            "novel_urf_flagged": urf_tot,
            "novel_urf_sensitivity_pct": round(100.0 * urf_tot / n_tot, 2),
            "in_dict_exact_accuracy_pct": np.nan,
            "parent_subtype_composition_ok": p_tot,
            "parent_subtype_composition_pct": round(100.0 * p_tot / n_tot, 2),
        })

    res_df = pd.DataFrame(out_rows)
    res_df.to_csv(RESULTS_DIR / "benchmark_out_of_dictionary_holdout.csv", index=False)
    print(res_df.to_string(index=False))


def sanitize_jphmm_id(seq_id: str) -> str:
    return re.sub(r"[^A-Za-z0-9_]", "_", seq_id)


def parse_hxb2_msa_map(hiv_aln_fas: Path):
    """Build 1-based MSA col -> 1-based HXB2 coordinate mapping using HIV_alignment.fas."""
    hxb2_seq = None
    in_hxb2 = False
    chunks = []
    with open(hiv_aln_fas) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if in_hxb2:
                    break
                in_hxb2 = "HXB2" in line
                chunks = []
            elif in_hxb2:
                chunks.append(line)
    hxb2_seq = "".join(chunks)
    if not hxb2_seq:
        raise RuntimeError("Could not find HXB2 in jpHMM HIV_alignment.fas")
    # 1-based MSA column index -> nearest HXB2 1-based position
    msa_to_hxb2 = [0] * (len(hxb2_seq) + 1)
    h_pos = 0
    for idx, ch in enumerate(hxb2_seq, start=1):
        if ch != "-":
            h_pos += 1
        msa_to_hxb2[idx] = max(1, h_pos)
    return msa_to_hxb2


def parse_jphmm_dir(jphmm_local_dir: Path, msa_to_hxb2, id_lookup):
    """Parse all chunk_XX/{recombination.txt, alignment_to_msa.txt} outputs."""
    msa_paths = {}
    for aln_file in sorted(jphmm_local_dir.glob("chunk_*/alignment_to_msa.txt")):
        cur_id = None
        cols = []
        with open(aln_file) as f:
            for line in f:
                s = line.strip()
                if not s or s.startswith("#"):
                    continue
                if s.startswith(">"):
                    if cur_id is not None:
                        msa_paths[cur_id] = cols
                    cur_id = s[1:].split()[0]
                    cols = []
                else:
                    cols.extend(int(x) for x in s.split())
            if cur_id is not None:
                msa_paths[cur_id] = cols

    calls = {}
    for rec_file in sorted(jphmm_local_dir.glob("chunk_*/recombination.txt")):
        cur_id = None
        segs = []
        with open(rec_file) as f:
            for line in f:
                s = line.strip()
                if not s or s.startswith("#"):
                    continue
                if s.startswith(">"):
                    if cur_id is not None:
                        calls[cur_id] = segs
                    cur_id = s[1:].split()[0]
                    segs = []
                else:
                    parts = s.split("\t")
                    if len(parts) >= 3:
                        segs.append((int(parts[0]), int(parts[1]), parts[2].strip()))
            if cur_id is not None:
                calls[cur_id] = segs

    parsed = {}
    max_msa = len(msa_to_hxb2) - 1
    for san_id, segs in calls.items():
        orig_sid = id_lookup.get(san_id)
        if not orig_sid:
            continue
        q_cols = msa_paths.get(san_id, [])
        hxb2_segs = []
        bps = []
        subtypes = set()
        for idx_s, (p_start, p_end, st) in enumerate(segs):
            if "Insertion" in st or st == "N/A":
                continue
            subtypes.add(st)
            if q_cols and 1 <= p_start <= len(q_cols) and 1 <= p_end <= len(q_cols):
                m_s = min(max(1, q_cols[p_start - 1]), max_msa)
                m_e = min(max(1, q_cols[p_end - 1]), max_msa)
                hs, he = msa_to_hxb2[m_s], msa_to_hxb2[m_e]
            else:
                hs, he = p_start, p_end
            if hxb2_segs and hxb2_segs[-1][2] != st:
                bps.append(hs)
            hxb2_segs.append((hs, he, st))
        parsed[orig_sid] = {
            "jphmm_segments": hxb2_segs,
            "jphmm_subtypes": subtypes,
            "jphmm_breakpoints": bps,
            "jphmm_is_pure": len(subtypes) == 1,
            "jphmm_formula": "+".join(f"{st}:{a}-{b}" for a, b, st in hxb2_segs),
        }
    return parsed


def extract_breakpoints_from_barcode(barcode: str):
    """Return sorted list of internal transition coordinates from a HXB2 barcode string like A1:790-2100|B:2100-4500."""
    intervals = []
    for chunk in str(barcode).split("|"):
        for m in re.finditer(r"(\d+)-(\d+)", chunk):
            intervals.append((int(m.group(1)), int(m.group(2))))
    intervals.sort()
    bps = []
    for i in range(len(intervals) - 1):
        bps.append(intervals[i][1])
    return sorted(set(bps))


def nearest_breakpoint_errors(pred_bps, true_bps):
    if not pred_bps or not true_bps:
        return []
    errs = []
    for tb in true_bps:
        errs.append(min(abs(pb - tb) for pb in pred_bps))
    return errs


def load_official_lanl_crf_breakpoints():
    """Extract independent published LANL CRF breakpoint coordinates from lanl_official_crf_segments.csv."""
    seg_path = DATA_DIR / "lanl_official_crf_segments.csv"
    seg_df = pd.read_csv(seg_path)
    out = {}
    for tag, grp in seg_df.groupby("short_tag"):
        rows_s = grp.sort_values(["hxb2_start", "hxb2_end"]).to_dict("records")
        bps = []
        # Identify uncertainty transition intervals ('/') vs hard internal switches inside 790..9417
        for i, r in enumerate(rows_s):
            lbl = str(r["parent_label"]).strip()
            s, e = int(r["hxb2_start"]), int(r["hxb2_end"])
            if "/" in lbl:
                mid = int(round(0.5 * (s + e)))
                if 790 < mid < 9417:
                    bps.append(mid)
            elif i + 1 < len(rows_s):
                nxt = rows_s[i + 1]
                nlbl = str(nxt["parent_label"]).strip()
                if "/" not in nlbl and nlbl != lbl:
                    bp = int(round(0.5 * (e + int(nxt["hxb2_start"]))))
                    if 790 < bp < 9417:
                        bps.append(bp)
        out[norm_crf(tag)] = sorted(set(bps))
    return out


def score_breakpoint_set(pred_bps, true_bps):
    """Return nearest absolute error per true bp, recall within 100/250 nt, and unmatched pred count."""
    errs = [min(abs(pb - tb) for pb in pred_bps) for tb in true_bps] if (pred_bps and true_bps) else []
    rec_100 = sum(1 for e in errs if e <= 100)
    rec_250 = sum(1 for e in errs if e <= 250)
    false_bps = sum(1 for pb in pred_bps if not true_bps or min(abs(pb - tb) for tb in true_bps) > 250)
    return errs, rec_100, rec_250, false_bps


def cmd_gen_synth500():
    """Generate 500 GTR+G4 synthetic HIV-1 mosaics, run eggman locally, and launch jpHMM on m3.local."""
    s18 = load_module("s18", PROJECT_ROOT / "scripts" / "18_simulate_tier1_null_fpr.py")
    trunk_strs = s18.load_basal_coding_trunks()
    pi, Q, trunk_arrays = s18.estimate_empirical_hiv_gtr(trunk_strs)

    hxb2_full = ""
    with open(BASAL_PROFILES_TSV) as f:
        for ln in f:
            parts = ln.rstrip("\n").split("\t")
            if len(parts) >= 2 and parts[0] == "HXB2":
                hxb2_full = parts[1][789:9417].upper()
                break
    hxb2_arr = np.array([s18.BASE_TO_IDX.get(c, 0) for c in hxb2_full], dtype=np.int8)

    tr_list = sorted(trunk_strs.keys())
    rng = np.random.default_rng(20261118)
    gamma_rates = s18.discrete_gamma_rates(alpha=0.38, k=4)
    site_cats = rng.integers(0, 4, size=8628)
    base_mult = np.ones(8628, dtype=float)

    n_bp_plan = [1] * 200 + [2] * 200 + [3] * 100
    rng.shuffle(n_bp_plan)

    synth_fa = RESULTS_DIR / "benchmark_synth500_gtr_mosaics.fasta"
    truth_csv = RESULTS_DIR / "benchmark_synth500_gtr_truth.csv"
    egg_csv = RESULTS_DIR / "benchmark_synth500_eggman.csv"

    truth_rows = []
    fasta_recs = []
    for idx, n_bp in enumerate(n_bp_plan):
        # Draw n_bp junction coordinates in HXB2 1100..9100 separated by >= 500 nt
        bps = []
        while len(bps) < n_bp:
            cand = sorted(rng.integers(1100, 9100, size=n_bp).tolist())
            if all(cand[i + 1] - cand[i] >= 500 for i in range(n_bp - 1)):
                bps = cand
        # Pick n_bp + 1 distinct adjacent trunks
        parents = [tr_list[int(rng.integers(0, len(tr_list)))]]
        for _ in range(n_bp):
            choices = [t for t in tr_list if t != parents[-1]]
            parents.append(choices[int(rng.integers(0, len(choices)))])

        anc = np.empty(8628, dtype=np.int8)
        bounds_h = [790] + bps + [9418]
        for s_i, tk in enumerate(parents):
            lo_idx = bounds_h[s_i] - 790
            hi_idx = bounds_h[s_i + 1] - 790
            tpl = trunk_arrays[tk][lo_idx:hi_idx]
            anc[lo_idx:hi_idx] = np.where(tpl >= 0, tpl, hxb2_arr[lo_idx:hi_idx])

        total_div = float(rng.uniform(0.030, 0.095))
        sim_arr = s18.evolve_gtr_gamma4(anc, total_div, Q, site_cats, gamma_rates, base_mult, rng)
        seq_str = "".join(s18.IDX_TO_BASE[sim_arr].tolist())
        sid = f"SYNTH500_{idx:03d}_bp{n_bp}"
        fasta_recs.append((sid, seq_str))
        truth_rows.append({
            "seq_id": sid,
            "n_bp": n_bp,
            "true_trunks": "|".join(parents),
            "true_trunks_norm": "|".join(sorted({("G" if p == "G0" else p) for p in parents})),
            "true_bps_hxb2": ";".join(str(b) for b in bps),
            "total_div": round(total_div, 5),
        })

    with open(synth_fa, "w") as f:
        for sid, s_str in fasta_recs:
            f.write(f">{sid}\n{s_str}\n")
    pd.DataFrame(truth_rows).to_csv(truth_csv, index=False)

    t0_egg = time.time()
    subprocess.run([str(RHIZTYPE_BIN), str(synth_fa), "--quiet", "--out", str(egg_csv)], check=True)
    dt_egg = time.time() - t0_egg
    print(f"[+] Generated 500 GTR+G4 mosaics -> {synth_fa}; EGGMAN finished in {dt_egg:.3f}s ({500/dt_egg:.1f} seq/s).")

    # Prepare 100 5-seq chunks & launch 28-worker jpHMM on m3.local via nohup
    local_chunks = DATA_DIR / "synth500_jphmm_chunks"
    local_chunks.mkdir(parents=True, exist_ok=True)
    for ci in range(100):
        cdir = local_chunks / f"chunk_{ci:02d}"
        cdir.mkdir(exist_ok=True)
        with open(cdir / "query.fas", "w") as cf:
            for sid, s_str in fasta_recs[ci * 5:(ci + 1) * 5]:
                cf.write(f">{sid}\n{s_str}\n")

    runner_py = local_chunks / "run_synth500_jphmm.py"
    runner_py.write_text(
        "import concurrent.futures, json, os, pathlib, subprocess, time\n"
        "JPHMM_DIR = pathlib.Path('/Users/sergei/hiv_recomb_bench/jphmm_build/rootfs/usr/local/jpHMM')\n"
        "IN_ROOT = pathlib.Path('/Users/sergei/hiv_recomb_bench/synth500_jphmm')\n"
        "OUT_ROOT = pathlib.Path('/Users/sergei/hiv_recomb_bench/results/jphmm_synth500_out')\n"
        "OUT_ROOT.mkdir(parents=True, exist_ok=True)\n"
        "def run_chunk(ci):\n"
        "    tag = f'chunk_{ci:02d}'\n"
        "    qfa = IN_ROOT / tag / 'query.fas'\n"
        "    odir = OUT_ROOT / tag\n"
        "    odir.mkdir(parents=True, exist_ok=True)\n"
        "    rec = odir / 'recombination.txt'\n"
        "    if rec.exists() and rec.stat().st_size > 0:\n"
        "        return 0.0\n"
        "    t0 = time.time()\n"
        "    subprocess.run([\n"
        "        str(JPHMM_DIR / 'src' / 'jpHMM'),\n"
        "        '-s', str(qfa), '-v', 'HIV', '-P', './priors/', '-I', './input/', '-o', str(odir) + '/'\n"
        "    ], cwd=str(JPHMM_DIR), check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)\n"
        "    return time.time() - t0\n"
        "t_wall0 = time.time()\n"
        "with concurrent.futures.ThreadPoolExecutor(max_workers=28) as ex:\n"
        "    secs = list(ex.map(run_chunk, range(100)))\n"
        "wall = time.time() - t_wall0\n"
        "meta = {'total_seqs': 500, 'n_chunks': 100, 'workers': 28, 'wall_sec': round(wall, 2), 'sum_worker_sec': round(sum(secs), 2)}\n"
        "(OUT_ROOT / 'jphmm_synth500_timing.json').write_text(json.dumps(meta, indent=2))\n"
        "print('DONE synth500 jpHMM:', meta)\n"
    )
    subprocess.run(["ssh", "-o", "BatchMode=yes", "m3.local", "mkdir -p ~/hiv_recomb_bench/synth500_jphmm ~/hiv_recomb_bench/results/jphmm_synth500_out"], check=True)
    subprocess.run(["rsync", "-az", "--delete", str(local_chunks) + "/", "m3.local:/Users/sergei/hiv_recomb_bench/synth500_jphmm/"], check=True)
    subprocess.run([
        "ssh", "-o", "BatchMode=yes", "m3.local",
        "nohup python3 ~/hiv_recomb_bench/synth500_jphmm/run_synth500_jphmm.py > ~/hiv_recomb_bench/synth500_jphmm.log 2>&1 < /dev/null &"
    ], check=True)
    print("[+] Launched 28-worker synth500 jpHMM background job on m3.local!")


def cmd_parse_jphmm():
    local_jphmm_dir = DATA_DIR / "jphmm_bench_out"
    local_jphmm_dir.mkdir(parents=True, exist_ok=True)
    local_ready = (
        (local_jphmm_dir / "jphmm_timing.json").exists()
        and len(list(local_jphmm_dir.glob("chunk_*/recombination.txt"))) == 64
    )
    if not local_ready:
        chk = subprocess.run(
            "ssh -o BatchMode=yes m3.local 'cat ~/hiv_recomb_bench/results/jphmm_bench_out/jphmm_timing.json 2>/dev/null || true'",
            shell=True, capture_output=True, text=True,
        )
        if not chk.stdout.strip():
            print("[!] jpHMM has not yet written jphmm_timing.json on m3.local.")
            sys.exit(2)

        print("[*] Syncing jpHMM outputs & HXB2 MSA reference from m3.local...")
        subprocess.run([
            "rsync", "-az",
            "m3.local:/Users/sergei/hiv_recomb_bench/results/jphmm_bench_out/",
            str(local_jphmm_dir) + "/",
        ], check=True)
    else:
        print("[*] Using locally cached jpHMM outputs (64 chunks + jphmm_timing.json).")
    hiv_aln_local = DATA_DIR / "jphmm_HIV_alignment.fas"
    if not hiv_aln_local.exists():
        subprocess.run([
            "rsync", "-az",
            "m3.local:/Users/sergei/hiv_recomb_bench/jphmm_build/rootfs/usr/local/jpHMM/input/HIV_alignment.fas",
            str(hiv_aln_local),
        ], check=True)

    with open(local_jphmm_dir / "jphmm_timing.json") as f:
        jp_timing = json.load(f)

    rt_timing = {}
    rt_meta_path = RESULTS_DIR / "benchmark_runtime_comparison.json"
    if rt_meta_path.exists():
        with open(rt_meta_path) as f:
            rt_timing = json.load(f)

    merged = pd.read_csv(RESULTS_DIR / "benchmark_rhiztype_vs_comet_merged.csv")
    extra_cols = [
        c for c in ["called_lineage", "called_category", "basal_formula", "basal_barcode", "compendium_matches"]
        if c not in merged.columns
    ]
    if extra_cols:
        rhiz_full = pd.read_csv(RESULTS_DIR / "benchmark_rhiztype_v02_output.csv")
        merged = merged.merge(
            rhiz_full[["sequence_id"] + extra_cols],
            left_on="seq_id", right_on="sequence_id", how="left",
        )

    id_lookup = {sanitize_jphmm_id(sid): sid for sid in merged["seq_id"]}
    msa_to_hxb2 = parse_hxb2_msa_map(hiv_aln_local)
    jp_parsed = parse_jphmm_dir(local_jphmm_dir, msa_to_hxb2, id_lookup)
    print(f"[+] Parsed jpHMM predictions for {len(jp_parsed)}/{len(merged)} sequences.")

    # Secondary empirical breakpoint benchmark against official published LANL CRF map midpoints
    lanl_official_bps = load_official_lanl_crf_breakpoints()

    # Independent published LANL CRF -> canonical parent-set dictionary (A-collapsed for fair cross-tool parity)
    _seg_df = pd.read_csv(DATA_DIR / "lanl_official_crf_segments.csv")
    _num_to_tag = {}
    _raw_lanl_mp = defaultdict(set)
    for _tag, _grp in _seg_df.groupby("short_tag"):
        _nk = norm_crf(str(_tag))
        _m_n = re.match(r"^(?:CRF)?(\d+)_", str(_tag))
        if _m_n:
            _num_to_tag[str(int(_m_n.group(1)))] = _nk
        for _lbl in _grp["parent_label"].astype(str):
            if "/" not in _lbl:
                _raw_lanl_mp[_nk].add(_lbl.strip())

    def _canon_sub_fair(tok):
        t = str(tok).strip().upper()
        if t in ("", "U", "UNDETERMINED", "UNATTRIBUTED", "OUTOFPANEL", "N/A", "CPX"):
            return None
        if t in ("G0", "G", "SUBTYPE_G"):
            return "G"
        if t in ("?E", "E", "AE"):
            return "E"
        if re.match(r"^A[1-8]?$", t):
            return "A"
        if t in ("B", "C", "D", "F1", "F2", "F", "H", "J", "K", "L"):
            return "F1" if t == "F" else t
        return None

    def _expand_lanl_tok(tok, depth=0):
        if depth > 5:
            return set()
        tu = str(tok).strip().upper()
        if tu in ("U", "UNDETERMINED", ""):
            return set()
        if tu in ("01_AE", "CRF01_AE", "CRF01", "01", "AE"):
            return {"A", "E"}
        if tu in ("02_AG", "CRF02_AG", "CRF02", "02"):
            return {"A", "G"}
        if tu in ("07_BC", "CRF07_BC", "CRF07", "07", "08_BC", "CRF08_BC", "CRF08", "08"):
            return {"B", "C"}
        if tu in ("12_BF", "CRF12_BF", "CRF12", "12"):
            return {"B", "F1"}
        m_n = re.match(r"^(?:CRF)?(\d+)(?:_[A-Za-z0-9]+)?$", str(tok).strip(), re.I)
        if m_n and str(int(m_n.group(1))) in _num_to_tag:
            nk = _num_to_tag[str(int(m_n.group(1)))]
            out = set()
            for st in _raw_lanl_mp.get(nk, ()):
                out |= _expand_lanl_tok(st, depth + 1)
            if out:
                return out
        c = _canon_sub_fair(tok)
        return {c} if c else set()

    _lanl_parents_fair = {
        k: set().union(*(_expand_lanl_tok(t) for t in vs)) if vs else set()
        for k, vs in _raw_lanl_mp.items()
    }

    def _independent_crf_parents(name):
        nk = norm_crf(name)
        if _lanl_parents_fair.get(nk):
            return set(_lanl_parents_fair[nk])
        m = re.match(r"^(?:CRF)?(\d+)_([A-Za-z0-9]+)$", str(name).strip(), re.I)
        if not m:
            return set()
        nk2 = _num_to_tag.get(str(int(m.group(1))))
        if nk2 and _lanl_parents_fair.get(nk2):
            return set(_lanl_parents_fair[nk2])
        suf = m.group(2)
        out = set()
        for code in ("109", "163", "63", "56", "55", "01", "02", "07", "08"):
            if code in suf:
                out |= _expand_lanl_tok(code)
                suf = suf.replace(code, "")
        for tok in re.findall(r"(A[1-8]|F[12]|[ABCDGHJKL])", suf):
            c = _canon_sub_fair(tok)
            if c:
                out.add(c)
        return out

    def _eg_parents_fair(row, use_t2=True):
        cat = str(row.get("called_category", ""))
        lin = str(row.get("called_lineage", ""))
        if use_t2 and cat == "Certified_CRF":
            p_c = _independent_crf_parents(lin)
            if p_c:
                return p_c
        if use_t2 and lin == "Subtype_G":
            return {"G"}
        bc = str(row.get("basal_barcode", "") or "")
        out = set()
        for part in bc.split("|"):
            if ":" in part:
                c = _canon_sub_fair(part.split(":")[0].strip())
                if c:
                    out.add(c)
        return out

    def _cm_parents_fair(comet_raw):
        c_cat, c_lbl = classify_comet_call(str(comet_raw))
        if c_cat == "pure":
            c = _canon_sub_fair(c_lbl)
            return {c} if c else set()
        if c_cat == "named_crf":
            p_c = _independent_crf_parents(c_lbl)
            if p_c:
                return p_c
        out = set()
        s = str(comet_raw)
        for crf_m in re.findall(r"\b(\d+_[A-Za-z0-9]+)\b", s):
            out |= _independent_crf_parents(crf_m)
        s_rem = re.sub(r"\b\d+_[A-Za-z0-9]+\b", " ", s)
        for tok in re.findall(r"(?<![A-Za-z0-9])(A[1-8]|F[12]|[ABCDGHJKL])(?![A-Za-z0-9_])", s_rem):
            c = _canon_sub_fair(tok)
            if c:
                out.add(c)
        return out

    def _jp_parents_fair(jp_subtypes):
        out = set()
        for sx in jp_subtypes:
            sx = str(sx).strip()
            if sx in ("01_AE", "AE", "CRF01_AE"):
                out |= {"A", "E"}
            else:
                c = _canon_sub_fair(sx)
                if c:
                    out.add(c)
        return out

    # Compare 3 ways across panels
    panel_slices = [
        ("Panel1A_LANL_Ref_PureSubtype", merged[merged["panel"] == "Panel1A_LANL_Ref_PureSubtype"]),
        ("Panel1B_LANL_Ref_CRF_01_to_118 (In COMET v2.4 Training)", merged[merged["panel"] == "Panel1B_LANL_Ref_CRF_01_to_118 (In COMET v2.4 Training)"]),
        ("Panel1C_LANL_Ref_CRF_119_to_157 (Post-COMET v2.4 Cutoff)", merged[merged["panel"] == "Panel1C_LANL_Ref_CRF_119_to_157 (Post-COMET v2.4 Cutoff)"]),
        ("Panel1BC_All_145_LANL_CRFs", merged[merged["panel"].str.startswith(("Panel1B", "Panel1C"))]),
        ("Panel1BC_Active_142_Mosaic_CRFs", merged[merged["panel"].str.startswith(("Panel1B", "Panel1C")) & (merged.get("is_active_mosaic_crf", True) == True)]),
        ("Panel2_Novel_FullGenome_CRFs (Audited 61)", merged[merged["panel"].str.startswith("Panel2")]),
        ("Panel3A_Novel_Partial_pol_Lineages", merged[merged["panel"] == "Panel3A_Novel_Partial_pol_Lineages"]),
        ("Panel3B_Novel_Partial_env_Lineages", merged[merged["panel"] == "Panel3B_Novel_Partial_env_Lineages"]),
        ("Panel4A_Clinical_PureSubtype", merged[merged["panel"] == "Panel4A_Clinical_PureSubtype"]),
        ("Panel4B_Clinical_NamedCRF", merged[merged["panel"] == "Panel4B_Clinical_NamedCRF"]),
        ("Panel4C_Clinical_URF", merged[merged["panel"] == "Panel4C_Clinical_URF"]),
    ]
    rows = []
    rt_bp_errs_founder = []
    jp_bp_errs_founder = []

    for p, sub in panel_slices:
        n = len(sub)
        if n == 0:
            continue
        rt_exact = rt_mos_ok = rt_par_ok = rt_par_eq = rt_t1_par_eq = 0
        cm_exact = cm_mos_ok = cm_par_ok = cm_par_eq = 0
        jp_pure_exact = jp_mos_ok = jp_par_ok = jp_par_eq = jp_cov = 0
        par_eval_n = jp_par_eval_n = 0
        rt_jacs, cm_jacs, jp_jacs = [], [], []
        is_indep_par_panel = p.startswith(("Panel1A", "Panel1B", "Panel1C", "Panel4A", "Panel4B"))

        for _, row in sub.iterrows():
            sid = row["seq_id"]
            gt = str(row["ground_truth_label"])
            exp_cat = str(row.get("expected_category", ""))
            gt_base = gt.split(":")[0]
            jp = jp_parsed.get(sid)
            if jp:
                jp_cov += 1

            if p.startswith(("Panel1A", "Panel4A")):
                c_gt = _canon_sub_fair(gt_base)
                gt_tr = {c_gt} if c_gt else set()
            elif p.startswith(("Panel1B", "Panel1C", "Panel4B")):
                gt_tr = _independent_crf_parents(gt_base)
            else:
                gt_tr = set()  # Exclude circular Panel2/3 and unannotated Panel4C from parent-set scoring

            # EGGMAN
            r_cat = str(row["called_category"])
            r_lin = str(row["called_lineage"])
            r_tr = _eg_parents_fair(row, use_t2=True)
            r_tr_t1 = _eg_parents_fair(row, use_t2=False)
            if p.startswith(("Panel1A", "Panel4A")):
                if r_cat in ("Pure_Subtype", "Ancestral_Recombinant_Subtype") and norm_pure(r_lin) == norm_pure(gt_base):
                    rt_exact += 1
                    rt_mos_ok += 1
            elif p.startswith(("Panel1B", "Panel1C", "Panel4B")):
                if r_cat == "Certified_CRF" and crf_matches(gt_base, r_lin):
                    rt_exact += 1
                if r_cat in ("Certified_CRF", "URF", "Novel_Circulating_CRF") and int(row.get("rhiztype_n_breakpoints", 0)) >= 1:
                    rt_mos_ok += 1
            elif p.startswith("Panel2"):
                if (exp_cat == "Certified_CRF" and r_cat == "Certified_CRF" and crf_matches(gt_base, r_lin)) or (exp_cat == "Novel_Circulating_CRF" and r_cat == "Novel_Circulating_CRF" and r_lin == gt_base) or (exp_cat == "URF" and r_cat == "URF"):
                    rt_exact += 1
                if r_cat in ("Certified_CRF", "Novel_Circulating_CRF", "URF"):
                    rt_mos_ok += 1
            else:
                if r_cat == "URF":
                    rt_mos_ok += 1

            # COMET
            c_cat, c_lbl = classify_comet_call(row["comet_subtype"])
            c_tr = _cm_parents_fair(row["comet_subtype"])
            if p.startswith(("Panel1A", "Panel4A")):
                if c_cat == "pure" and norm_pure(c_lbl) == norm_pure(gt_base):
                    cm_exact += 1
                    cm_mos_ok += 1
            elif p.startswith(("Panel1B", "Panel4B")) and c_cat == "named_crf" and crf_matches(gt_base, c_lbl):
                cm_exact += 1
                cm_mos_ok += 1
            elif p.startswith(("Panel1C", "Panel2", "Panel3", "Panel4C")) and c_cat == "unassigned":
                cm_mos_ok += 1

            # jpHMM
            jp_tr = set()
            if jp:
                jp_tr = _jp_parents_fair(jp["jphmm_subtypes"])
                if p.startswith(("Panel1A", "Panel4A")):
                    if jp["jphmm_is_pure"] and _canon_sub_fair(gt_base) in jp_tr:
                        jp_pure_exact += 1
                        jp_mos_ok += 1
                else:
                    if len(jp_tr) >= 2:
                        jp_mos_ok += 1

            if is_indep_par_panel and gt_tr:
                par_eval_n += 1
                ov_r, eq_r, jac_r = compute_parent_set_stats(gt_tr, r_tr)
                _, eq_rt1, _ = compute_parent_set_stats(gt_tr, r_tr_t1)
                rt_par_ok += ov_r
                rt_par_eq += eq_r
                rt_t1_par_eq += eq_rt1
                rt_jacs.append(jac_r)

                ov_c, eq_c, jac_c = compute_parent_set_stats(gt_tr, c_tr)
                cm_par_ok += ov_c
                cm_par_eq += eq_c
                cm_jacs.append(jac_c)

                if jp:
                    jp_par_eval_n += 1
                    ov_j, eq_j, jac_j = compute_parent_set_stats(gt_tr, jp_tr)
                    jp_par_ok += ov_j
                    jp_par_eq += eq_j
                    jp_jacs.append(jac_j)

            if jp and p in ("Panel1B_LANL_Ref_CRF_01_to_118 (In COMET v2.4 Training)", "Panel4B_Clinical_NamedCRF") and norm_crf(gt_base) in MAJOR_FOUNDER_CRFS:
                ref_bps = lanl_official_bps.get(norm_crf(gt_base), [])
                if ref_bps:
                    rt_bps = extract_breakpoints_from_barcode(row.get("basal_barcode", ""))
                    rt_bp_errs_founder.extend(nearest_breakpoint_errors(rt_bps, ref_bps))
                    jp_bp_errs_founder.extend(nearest_breakpoint_errors(jp["jphmm_breakpoints"], ref_bps))

        rows.append({
            "panel": p,
            "n_sequences": n,
            "parent_eval_n": par_eval_n,
            "rhiztype_exact_lineage_pct": round(100.0 * rt_exact / n, 1),
            "comet_exact_lineage_pct": round(100.0 * cm_exact / n, 1),
            "jphmm_pure_exact_pct": round(100.0 * jp_pure_exact / n, 1) if p.startswith(("Panel1A", "Panel4A")) else np.nan,
            "rhiztype_mosaic_or_pure_state_pct": round(100.0 * rt_mos_ok / n, 1),
            "comet_mosaic_or_pure_state_pct": round(100.0 * cm_mos_ok / n, 1),
            "jphmm_mosaic_or_pure_state_pct": round(100.0 * jp_mos_ok / n, 1),
            "rhiztype_parent_comp_pct": round(100.0 * rt_par_ok / par_eval_n, 1) if par_eval_n else np.nan,
            "comet_parent_comp_pct": round(100.0 * cm_par_ok / par_eval_n, 1) if par_eval_n else np.nan,
            "jphmm_parent_comp_pct": round(100.0 * jp_par_ok / jp_par_eval_n, 1) if jp_par_eval_n else np.nan,
            "rhiztype_parent_exact_set_pct": round(100.0 * rt_par_eq / par_eval_n, 1) if par_eval_n else np.nan,
            "rhiztype_t1_parent_exact_set_pct": round(100.0 * rt_t1_par_eq / par_eval_n, 1) if par_eval_n else np.nan,
            "comet_parent_exact_set_pct": round(100.0 * cm_par_eq / par_eval_n, 1) if par_eval_n else np.nan,
            "jphmm_parent_exact_set_pct": round(100.0 * jp_par_eq / jp_par_eval_n, 1) if jp_par_eval_n else np.nan,
            "rhiztype_parent_mean_jaccard": round(float(np.mean(rt_jacs)), 3) if rt_jacs else np.nan,
            "comet_parent_mean_jaccard": round(float(np.mean(cm_jacs)), 3) if cm_jacs else np.nan,
            "jphmm_parent_mean_jaccard": round(float(np.mean(jp_jacs)), 3) if jp_jacs else np.nan,
            "jphmm_parsed_count": jp_cov,
        })

    out_df = pd.DataFrame(rows)
    out_df.to_csv(RESULTS_DIR / "benchmark_three_way_comparison.csv", index=False)

    # Check if 500-mosaic GTR+G4 jpHMM output is ready & score unbiased synthetic benchmark
    synth_summary = {}
    truth_csv = RESULTS_DIR / "benchmark_synth500_gtr_truth.csv"
    egg_csv = RESULTS_DIR / "benchmark_synth500_eggman.csv"
    local_s500_jp = DATA_DIR / "jphmm_synth500_out"
    if truth_csv.exists() and egg_csv.exists():
        local_s500_ready = (local_s500_jp / "jphmm_synth500_timing.json").exists() and len(list(local_s500_jp.glob("chunk_*/recombination.txt"))) == 100
        if not local_s500_ready:
            chk_s = subprocess.run("ssh -o BatchMode=yes m3.local 'cat ~/hiv_recomb_bench/results/jphmm_synth500_out/jphmm_synth500_timing.json 2>/dev/null || true'", shell=True, capture_output=True, text=True)
            if chk_s.stdout.strip():
                local_s500_jp.mkdir(parents=True, exist_ok=True)
                subprocess.run(["rsync", "-az", "m3.local:/Users/sergei/hiv_recomb_bench/results/jphmm_synth500_out/", str(local_s500_jp) + "/"], check=True)
                local_s500_ready = True
        if local_s500_ready:
            tr_df = pd.read_csv(truth_csv)
            eg_df = pd.read_csv(egg_csv).set_index("sequence_id")
            s500_lookup = {sanitize_jphmm_id(s): s for s in tr_df["seq_id"]}
            jp_s500 = parse_jphmm_dir(local_s500_jp, msa_to_hxb2, s500_lookup)
            s500_timing = json.loads((local_s500_jp / "jphmm_synth500_timing.json").read_text())

            rt_errs, jp_errs = [], []
            rt_r100 = rt_r250 = rt_fp = rt_tot_bp = rt_par_eq = 0
            jp_r100 = jp_r250 = jp_fp = jp_tot_bp = jp_par_eq = 0
            rt_jacs_s, jp_jacs_s = [], []
            for _, tr in tr_df.iterrows():
                sid = tr["seq_id"]
                t_bps = [int(x) for x in str(tr["true_bps_hxb2"]).split(";") if x]
                gt_set = set(str(tr["true_trunks_norm"]).split("|"))
                # For jpHMM fair comparison on synth500, map L/A6->A1 if evaluating jpHMM's panel, AND also score raw strict set
                er = eg_df.loc[sid]
                r_bps = extract_breakpoints_from_barcode(str(er["basal_barcode"]))
                e_r, r100, r250, fp_r = score_breakpoint_set(r_bps, t_bps)
                rt_errs.extend(e_r); rt_r100 += r100; rt_r250 += r250; rt_fp += fp_r; rt_tot_bp += len(r_bps)
                r_set = extract_parental_trunks_from_formula(str(er["basal_formula"]))
                _, eq_r, j_r = compute_parent_set_stats(gt_set, r_set)
                rt_par_eq += eq_r; rt_jacs_s.append(j_r)

                jp = jp_s500.get(sid)
                if jp:
                    j_bps = jp["jphmm_breakpoints"]
                    e_j, j100, j250, fp_j = score_breakpoint_set(j_bps, t_bps)
                    jp_errs.extend(e_j); jp_r100 += j100; jp_r250 += j250; jp_fp += fp_j; jp_tot_bp += len(j_bps)
                    j_set = extract_jphmm_subtypes_set(jp["jphmm_subtypes"])
                    _, eq_j, j_j = compute_parent_set_stats(gt_set, j_set)
                    jp_par_eq += eq_j; jp_jacs_s.append(j_j)

            n_true_bps = int(tr_df["n_bp"].sum())
            synth_rows = [
                {
                    "tool": "EGGMAN",
                    "n_sequences": len(tr_df),
                    "n_true_breakpoints": n_true_bps,
                    "n_predicted_breakpoints": rt_tot_bp,
                    "median_abs_bp_err_nt": round(float(np.median(rt_errs)), 1) if rt_errs else None,
                    "iqr_bp_err_nt": f"[{np.percentile(rt_errs, 25):.1f}, {np.percentile(rt_errs, 75):.1f}]" if rt_errs else "",
                    "recall_within_100nt_pct": round(100.0 * rt_r100 / n_true_bps, 2),
                    "recall_within_250nt_pct": round(100.0 * rt_r250 / n_true_bps, 2),
                    "false_junction_rate_pct": round(100.0 * rt_fp / max(1, rt_tot_bp), 2),
                    "strict_parent_set_eq_pct": round(100.0 * rt_par_eq / len(tr_df), 2),
                    "mean_parent_jaccard": round(float(np.mean(rt_jacs_s)), 3),
                },
                {
                    "tool": "jpHMM",
                    "n_sequences": len(jp_s500),
                    "n_true_breakpoints": n_true_bps,
                    "n_predicted_breakpoints": jp_tot_bp,
                    "median_abs_bp_err_nt": round(float(np.median(jp_errs)), 1) if jp_errs else None,
                    "iqr_bp_err_nt": f"[{np.percentile(jp_errs, 25):.1f}, {np.percentile(jp_errs, 75):.1f}]" if jp_errs else "",
                    "recall_within_100nt_pct": round(100.0 * jp_r100 / n_true_bps, 2),
                    "recall_within_250nt_pct": round(100.0 * jp_r250 / n_true_bps, 2),
                    "false_junction_rate_pct": round(100.0 * jp_fp / max(1, jp_tot_bp), 2),
                    "strict_parent_set_eq_pct": round(100.0 * jp_par_eq / max(1, len(jp_s500)), 2),
                    "mean_parent_jaccard": round(float(np.mean(jp_jacs_s)), 3) if jp_jacs_s else None,
                },
            ]
            pd.DataFrame(synth_rows).to_csv(RESULTS_DIR / "benchmark_synth500_breakpoint_comparison.csv", index=False)
            synth_summary = {"jphmm_synth500_timing": s500_timing, "comparison": synth_rows}
            print("[+] Saved 500-mosaic GTR+G4 benchmark -> results/benchmark_synth500_breakpoint_comparison.csv")

    meta_cmp = {
        "jphmm_timing": jp_timing,
        "rhiztype_timing": rt_timing,
        "jphmm_sec_per_seq_wall": round(jp_timing["wall_sec"] / jp_timing["total_seqs"], 3),
        "jphmm_sec_per_seq_cpu": round(jp_timing["sum_worker_sec"] / jp_timing["total_seqs"], 3),
        "founder_crf_ref_source": "lanl_official_crf_segments.csv",
        "founder_crf_n_bp_evaluated": len(jp_bp_errs_founder),
        "rhiztype_founder_median_bp_err_nt": round(float(np.median(rt_bp_errs_founder)), 1) if rt_bp_errs_founder else None,
        "rhiztype_founder_iqr_bp_err_nt": [round(float(np.percentile(rt_bp_errs_founder, 25)), 1), round(float(np.percentile(rt_bp_errs_founder, 75)), 1)] if rt_bp_errs_founder else None,
        "jphmm_founder_median_bp_err_nt": round(float(np.median(jp_bp_errs_founder)), 1) if jp_bp_errs_founder else None,
        "jphmm_founder_iqr_bp_err_nt": [round(float(np.percentile(jp_bp_errs_founder, 25)), 1), round(float(np.percentile(jp_bp_errs_founder, 75)), 1)] if jp_bp_errs_founder else None,
        "synth500_gtr_benchmark": synth_summary,
    }
    with open(RESULTS_DIR / "benchmark_three_way_meta.json", "w") as f:
        json.dump(meta_cmp, f, indent=2)
    print(out_df.to_string(index=False))
    print("[+] Three-way meta:", json.dumps(meta_cmp, indent=2))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sync-comet", action="store_true")
    parser.add_argument("--loso", action="store_true")
    parser.add_argument("--ood", action="store_true")
    parser.add_argument("--gen-synth500", action="store_true")
    parser.add_argument("--parse-jphmm", action="store_true")
    parser.add_argument("--all-local", action="store_true", help="Run sync-comet, loso, and ood sequentially")
    args = parser.parse_args()

    if args.all_local or not any([args.sync_comet, args.loso, args.ood, args.gen_synth500, args.parse_jphmm]):
        cmd_sync_comet()
        cmd_loso()
        cmd_ood()
    else:
        if args.sync_comet:
            cmd_sync_comet()
        if args.loso:
            cmd_loso()
        if args.ood:
            cmd_ood()
        if args.gen_synth500:
            cmd_gen_synth500()
        if args.parse_jphmm:
            cmd_parse_jphmm()


if __name__ == "__main__":
    main()
