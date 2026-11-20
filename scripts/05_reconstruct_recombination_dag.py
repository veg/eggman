#!/usr/bin/env python3
"""
05_reconstruct_recombination_dag.py
===================================
Reconstructs the multi-generational Recombination Directed Acyclic Graph (DAG)
and breakpoint stratigraphy of HIV-1 Group M across all 145 sequence-represented CRFs.

Separates:
  1. Constituent Core Genome Formula (Classification View):
     Deconvolves the genome strictly into elementary basal trunks (A1..L, ?E).
  2. Historical Event Sequence / Genealogical Formula (Evolutionary History View):
     Partitions every verified breakpoint into:
       - Inherited Ancestral Junctions (structural synapomorphies inherited from an older
         Gen-1 or Gen-2 circulating recombinant with t_MRCA(Ancestor) < t_MRCA(Child)), vs.
       - De Novo Junctions created at the birth of the recombinant lineage.
"""

from pathlib import Path
from collections import defaultdict
import numpy as np
import pandas as pd

PROJECT_ROOT = Path("/Users/sergei/Projects/HIV-recombination")
RESULTS_DIR = PROJECT_ROOT / "results"


def main():
    comp = pd.read_csv(RESULTS_DIR / "hiv1_descriptive_nomenclature_compendium.csv")
    bps = pd.read_csv(RESULTS_DIR / "hiv1_deconvolved_breakpoints.csv")

    tmrca_map = dict(zip(comp["legacy_label"], comp["t_mrca"]))
    formula_map = dict(zip(comp["legacy_label"], comp["descriptive_name"]))
    barcode_map = dict(zip(comp["legacy_label"], comp["hxb2_barcode"]))
    status_map = dict(zip(comp["legacy_label"], comp["taxonomic_status"]))

    # Primary Gen-1 and major Gen-2 founder recombinants ordered chronologically by compendium t_MRCA/origin
    founder_defs = [
        ("CRF02_AG",  1),
        ("CRF01_AE",  1),
        ("CRF11_cpx", 1),
        ("CRF06_cpx", 1),
        ("CRF12_BF",  1),
        ("CRF07_BC",  1),
        ("CRF13_cpx", 1),
        ("CRF08_BC",  1),
        ("CRF15_01B", 2),
        ("CRF03_A6B", 1),
        ("CRF33_01B", 2),
        ("CRF55_01B", 2),
        ("CRF63_02A6", 2),
    ]
    founders = sorted(
        [(f_lbl, float(tmrca_map[f_lbl]), f_gen) for f_lbl, f_gen in founder_defs],
        key=lambda x: x[1],
    )

    founder_bps = {}
    for f_lbl, f_t, f_gen in founders:
        sub = bps[bps["legacy_label"] == f_lbl]
        founder_bps[f_lbl] = [
            (r["left_trunk"], r["right_trunk"], int(r["bp_hxb2"]), int(r["ci_hxb2_low"]), int(r["ci_hxb2_high"]), f_t, f_gen)
            for _, r in sub.iterrows()
        ]

    # Note: CRF55_01B itself inherits A1->?E@5096 from CRF01_AE (Gen 1 -> Gen 2)
    dag_rows = []
    for _, row in comp.iterrows():
        lbl = row["legacy_label"]
        if row["legacy_category"] == "Nominal_Pure_Subtype":
            continue

        t_child = float(row["t_mrca"]) if pd.notna(row["t_mrca"]) else 2020.0
        sub_bps = bps[bps["legacy_label"] == lbl]

        if len(sub_bps) == 0:
            tax_st = str(status_map.get(lbl, ""))
            zero_reason = "None (Out-of-Panel Donor)" if "Out_of_Panel_Donor" in tax_st else "None (Clonal Trunk)"
            dag_rows.append({
                "legacy_label": lbl,
                "t_mrca": t_child,
                "generation": 0,
                "constituent_formula": row["descriptive_name"],
                "immediate_recombinant_parents": zero_reason,
                "n_total_bps": 0,
                "n_inherited_bps": 0,
                "n_denovo_bps": 0,
                "inherited_junctions": "None",
                "denovo_junctions": "None",
            })
            continue

        inherited_list = []
        denovo_list = []
        rec_parents = set()
        max_parent_gen = 0

        for _, b in sub_bps.iterrows():
            L, R = b["left_trunk"], b["right_trunk"]
            pos = int(b["bp_hxb2"])
            lo, hi = int(b["ci_hxb2_low"]), int(b["ci_hxb2_high"])

            # Find the most specific (latest-generation, still older than child) founder that already carries this exact junction
            matched_founder = None
            for f_lbl, f_t, f_gen in founders:
                if f_lbl == lbl:
                    continue
                if f_t >= t_child and f_lbl not in ("CRF01_AE", "CRF02_AG"):
                    continue
                for fL, fR, fpos, flo, fhi, _, _ in founder_bps[f_lbl]:
                    if L == fL and R == fR and (abs(pos - fpos) <= 50 or (max(lo, flo) <= min(hi, fhi))):
                        # Attribute the canonical CRF01_AE A1->?E@~5097-5151 junction to Gen-1 CRF01_AE
                        is_01ae_env_junc = (L == "A1" and R == "?E" and 5050 <= pos <= 5180)
                        if matched_founder is None or (f_gen > matched_founder[2] and not is_01ae_env_junc):
                            matched_founder = (f_lbl, fpos, f_gen, f_t)

            if matched_founder is not None:
                f_lbl, fpos, f_gen, f_t = matched_founder
                inherited_list.append(f"{L}->{R}@{pos}(from {f_lbl}@{int(round(f_t))})")
                rec_parents.add(f_lbl)
                if f_gen > max_parent_gen:
                    max_parent_gen = f_gen
            else:
                denovo_list.append(f"{L}->{R}@{pos}")

        # If a lineage inherited specific CRF55_01B RT junctions AND CRF01_AE's 5096 junction,
        # the 5096 junction came *via* CRF55_01B (3rd-gen cascade: CRF01_AE -> CRF55_01B -> Child)
        if "CRF55_01B" in rec_parents and "CRF01_AE" in rec_parents:
            rec_parents.remove("CRF01_AE")
            rec_parents.add("CRF55_01B(<-CRF01_AE)")

        gen = 1 if max_parent_gen == 0 else (max_parent_gen + 1)
        dag_rows.append({
            "legacy_label": lbl,
            "t_mrca": t_child,
            "generation": gen,
            "constituent_formula": row["descriptive_name"],
            "immediate_recombinant_parents": " + ".join(sorted(rec_parents)) if rec_parents else "None (De Novo Primary Mosaic)",
            "n_total_bps": len(sub_bps),
            "n_inherited_bps": len(inherited_list),
            "n_denovo_bps": len(denovo_list),
            "inherited_junctions": "; ".join(inherited_list) if inherited_list else "None",
            "denovo_junctions": "; ".join(denovo_list) if denovo_list else "None",
        })

    dag_df = pd.DataFrame(dag_rows)
    dag_df.to_csv(RESULTS_DIR / "hiv1_recombination_genealogy_dag.csv", index=False)

    print("=" * 95)
    print(f"STEP 5: MULTI-GENERATIONAL RECOMBINATION DAG & EVENT STRATIGRAPHY ACROSS {len(dag_df)} CRFs")
    print("=" * 95)
    for g, cnt in dag_df["generation"].value_counts().sort_index().items():
        if g == 0:
            sub0 = dag_df[dag_df["generation"] == 0]
            n_oop = int((sub0["immediate_recombinant_parents"] == "None (Out-of-Panel Donor)").sum())
            n_cln = int((sub0["immediate_recombinant_parents"] == "None (Clonal Trunk)").sum())
            print(f"  * Gen 0a: Out-of-Panel Donor (Unrepresented A5/U; 0 breakpoints)       : {n_oop:3d} ({100.0*n_oop/len(dag_df):5.1f}%)")
            print(f"  * Gen 0b: True Collapsed to Clonal Trunk (0 breakpoints)               : {n_cln:3d} ({100.0*n_cln/len(dag_df):5.1f}%)")
        else:
            g_desc = {
                1: "Gen 1: Primary De Novo Recombinants (between pure basal trunks)",
                2: "Gen 2: Secondary Recombinants (inheriting junctions from a Gen-1 CRF)",
                3: "Gen 3: Tertiary Serial Recombinants (inheriting junctions from a Gen-2 CRF)",
            }.get(g, f"Gen {g}")
            print(f"  * {g_desc:<72s}: {cnt:3d} ({100.0*cnt/len(dag_df):5.1f}%)")

    mos_dag = dag_df[dag_df["n_total_bps"] >= 1]
    print(f"  -> Active mosaic CRFs (n_total_bps >= 1): {len(mos_dag)} / {len(dag_df)}")

    print("\n--- SHOWCASE OF MULTI-GENERATIONAL SERIAL RECOMBINATION CASCADES ---")
    showcase_crfs = [
        "CRF02_AG", "CRF63_02A6", "CRF95_02B",
        "CRF01_AE", "CRF15_01B", "CRF33_01B", "CRF55_01B",
        "CRF07_BC", "CRF102_0107", "CRF117_0107",
        "CRF126_0755", "CRF155_0755", "CRF153_55B",
        "CRF06_cpx", "CRF32_06A6",
        "CRF12_BF", "CRF17_BF1", "CRF89_BF1"
    ]
    sub = dag_df[dag_df["legacy_label"].isin(showcase_crfs)]
    for _, r in sub.iterrows():
        print(f"\n  [{r['legacy_label']}] (Gen {r['generation']}, t={r['t_mrca']}) -> Constituent: {r['constituent_formula']}")
        print(f"      Recombinant Parents : {r['immediate_recombinant_parents']}")
        print(f"      Inherited Junctions : {r['inherited_junctions']}")
        print(f"      De Novo Junctions   : {r['denovo_junctions']}")


if __name__ == "__main__":
    main()
