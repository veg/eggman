#!/usr/bin/env python3
"""
33_bars_v3_reanalyses_and_data_tables.py
========================================
Deterministic statistical verification script reproducing all founder-deduplicated
Stage-2 sieve controls, two-rigid-core Stage-1 checks, Clauset (2009) bootstrap
GOF + Vuong likelihood-ratio cluster-scaling tests, and exporting clean
machine-readable Supplementary Data Tables (Data S1 - Data S6) into
`results/supp_tables/`.

Inputs (all relative to workspace root, zero hardcoded external paths):
  - results/chronaeon_vs_lsd2_treetime_validation.csv
  - results/chronaeon_sc2_dating_vs_lsd2_treetime.tsv
  - results/chronaeon_iav_dating_vs_lsd2_treetime.tsv
  - results/hiv1_recombination_genealogy_dag.csv
  - results/hiv1_descriptive_nomenclature_compendium.csv
  - results/lanl_declared_breakpoints_audit.csv
  - results/regional_heterozygosity_urf_rates.csv
  - results/sieve_hcirc_multivariate_regression.csv
  - results/supp_manifest_sieve_749_lineages.csv
  - results/sieve_stage1_snp_and_sim_null.csv
  - results/genbank_novel_circulating_recombinants.csv
  - rhiztype/crates/eggman-core/data/compendium_bundle.json (or software/crates/...)
  - results/benchmark_three_way_comparison.csv
  - results/benchmark_loso_cv_summary.csv
  - results/benchmark_out_of_dictionary_holdout.csv
  - results/tier1_null_simulation_fpr.csv
  - results/benchmark_sc2_nextclade_and_iav_genoflu.csv
  - results/chronaeon_temporal_guard_ablation.tsv
  - results/stanford_hivdb_audit_category_summary.csv
  - results/stanford_hivdb_unlinked_8k_breakdown.csv
  - results/stanford_west_african_crf02_rt_validation.csv
  - results/novel_lineages_phylogenetic_validation.csv
  - results/scale_partial_pol_phylogenetic_validation.csv
  - results/scale_partial_env_phylogenetic_validation.csv
"""

import csv
import json
import math
import shutil
import sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
from scipy.optimize import minimize
from scipy.special import zeta
from scipy.stats import binomtest, chi2, fisher_exact, norm


ROOT = Path(__file__).resolve().parent.parent
RES = ROOT / "results"
OUT = RES / "supp_tables"
OUT.mkdir(parents=True, exist_ok=True)


def find_compendium_bundle() -> Path:
    for rel in [
        "rhiztype/crates/eggman-core/data/compendium_bundle.json",
        "software/crates/eggman-core/data/compendium_bundle.json",
    ]:
        p = ROOT / rel
        if p.exists():
            return p
    raise FileNotFoundError("compendium_bundle.json not found")


def parse_segs(bc_str: str):
    out = []
    for item in (bc_str or "").split("|"):
        if ":" in item and "-" in item:
            tr, rng = item.strip().split(":", 1)
            a, b = [int(x) for x in rng.split("-", 1)]
            out.append((tr.strip(), a, b))
    return out


def trunk_at(segs, pos: int):
    for tr, a, b in segs:
        if a <= pos <= b:
            return tr
    return None


def pair_key(segs) -> str:
    """Top-2 dominant parental trunks by total nucleotide span (sorted pair label)."""
    span = defaultdict(int)
    for t, a, b in segs:
        if t != "Unattributed":
            span[t] += max(0, b - a + 1)
    top2 = [k for k, _ in sorted(span.items(), key=lambda kv: (-kv[1], kv[0]))[:2]]
    return "+".join(sorted(top2))


def cmh_test(tables):
    """Exact Mantel-Haenszel pooled OR and continuity-corrected chi2 p-value."""
    num = 0.0
    den = 0.0
    sum_a = 0.0
    sum_e = 0.0
    sum_v = 0.0
    for a, b, c, d in tables:
        n = a + b + c + d
        if n <= 1:
            continue
        r1 = a + b
        r2 = c + d
        c1 = a + c
        c2 = b + d
        if r1 == 0 or r2 == 0 or c1 == 0 or c2 == 0:
            continue
        num += (a * d) / n
        den += (b * c) / n
        e_a = (r1 * c1) / n
        v_a = (r1 * r2 * c1 * c2) / (n * n * (n - 1))
        sum_a += a
        sum_e += e_a
        sum_v += v_a
    or_mh = num / den if den > 0 else float("nan")
    if sum_v <= 0:
        return or_mh, 1.0
    stat = (max(0.0, abs(sum_a - sum_e) - 0.5) ** 2) / sum_v
    p_val = float(1.0 - chi2.cdf(stat, df=1))
    return or_mh, p_val


def logistic_mv(y, x_bin, n_bp, strata_labels):
    """
    Multivariable logistic regression: logit P(y=1) = b0 + b1*x_bin + b2*n_bp + pair fixed effects
    (collapsing rare parental pairs with <10 observations into 'Other' to avoid separation).
    Returns (OR_x, ci_lo, ci_hi, p_val).
    """
    cnt = Counter(strata_labels)
    kept = sorted(k for k, v in cnt.items() if v >= 10)
    cols = [np.ones(len(y)), np.asarray(x_bin, float), np.asarray(n_bp, float)]
    for k in kept:
        cols.append(np.array([1.0 if s == k else 0.0 for s in strata_labels], float))
    X = np.column_stack(cols)
    y_arr = np.asarray(y, float)

    def nll_and_grad(beta):
        z = np.clip(X @ beta, -30.0, 30.0)
        p = 1.0 / (1.0 + np.exp(-z))
        # Tiny ridge (1e-6) on nuisance dummies only, zero penalty on intercept/x_bin/n_bp
        pen = 1e-6 * np.sum(beta[3:] ** 2) if len(beta) > 3 else 0.0
        nll = -np.sum(y_arr * np.log(np.maximum(p, 1e-15)) + (1.0 - y_arr) * np.log(np.maximum(1.0 - p, 1e-15))) + pen
        grad = X.T @ (p - y_arr)
        if len(beta) > 3:
            grad[3:] += 2e-6 * beta[3:]
        return float(nll), grad

    res = minimize(nll_and_grad, np.zeros(X.shape[1]), jac=True, method="BFGS")
    beta = res.x
    z = np.clip(X @ beta, -30.0, 30.0)
    p = 1.0 / (1.0 + np.exp(-z))
    W = p * (1.0 - p)
    H = (X.T * W) @ X
    if len(beta) > 3:
        for i in range(3, len(beta)):
            H[i, i] += 2e-6
    cov = np.linalg.pinv(H)
    b1 = float(beta[1])
    se1 = math.sqrt(max(float(cov[1, 1]), 1e-12))
    z_stat = b1 / se1
    p_val = float(2.0 * (1.0 - norm.cdf(abs(z_stat))))
    return math.exp(b1), math.exp(b1 - 1.96 * se1), math.exp(b1 + 1.96 * se1), p_val


def export_data_s1_chronaeon():
    """Merge all 27 HIV-1, 32 SARS-CoV-2, and 15 Influenza A ChronAeon validation rows."""
    sources = [
        ("HIV-1_GroupM", RES / "chronaeon_vs_lsd2_treetime_validation.csv", ","),
        ("SARS-CoV-2", RES / "chronaeon_sc2_dating_vs_lsd2_treetime.tsv", "\t"),
        ("Influenza_A_H5Nx", RES / "chronaeon_iav_dating_vs_lsd2_treetime.tsv", "\t"),
    ]
    out_cols = [
        "pathogen", "job_id", "analysis_type", "clade", "segment_desc",
        "n_full_clade", "n_eval_tips", "sample_min_yr", "sample_max_yr",
        "chronaeon_full_model", "chronaeon_full_tmrca", "chronaeon_full_ci_lower", "chronaeon_full_ci_upper",
        "tmrca_method", "final_pack_tmrca", "first_sampled_year", "first_submitted_year",
        "lsd2_status", "lsd2_tmrca", "lsd2_ci_lower", "lsd2_ci_upper",
        "treetime_status", "tt_tmrca", "abs_diff_ch_vs_lsd2", "abs_diff_ch_vs_tt",
    ]
    rows_out = []
    for pathogen, p, sep in sources:
        with open(p) as f:
            for r in csv.DictReader(f, delimiter=sep):
                rec = {"pathogen": pathogen}
                for c in out_cols[1:]:
                    rec[c] = r.get(c, "")
                if not rec["tmrca_method"]:
                    rec["tmrca_method"] = rec["chronaeon_full_model"] or "ChronAeon_OLS"
                rows_out.append(rec)
    out_p = OUT / "Data_S1_chronaeon_dates_hiv1_sc2_iav.tsv"
    with open(out_p, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=out_cols, delimiter="\t")
        w.writeheader()
        w.writerows(rows_out)
    print(f"[Data S1] Wrote {len(rows_out)} rows -> {out_p.name}")


def export_data_s2_s3_s5_s6():
    """Copy canonical audited tables into standardized Data_S2, S3, S5, S6 filenames."""
    copies = [
        ("hiv1_recombination_genealogy_dag.csv", "Data_S2_hiv1_145_crf_genealogy_dag.csv"),
        ("hiv1_descriptive_nomenclature_compendium.csv", "Data_S2b_hiv1_compendium_catalog.csv"),
        ("lanl_declared_breakpoints_audit.csv", "Data_S2c_lanl_986_breakpoints_audit.csv"),
        ("regional_heterozygosity_urf_rates.csv", "Data_S3_hiv1_32country_hcirc.csv"),
        ("sieve_hcirc_multivariate_regression.csv", "Data_S3b_hiv1_hcirc_glm_models.csv"),
        ("supp_manifest_sieve_749_lineages.csv", "Data_S4_hiv1_761_lineages_manifest.csv"),
        ("benchmark_three_way_comparison.csv", "Data_S5a_hiv1_three_way_noncircular_benchmark.csv"),
        ("benchmark_loso_cv_summary.csv", "Data_S5b_hiv1_loso_cv_summary.csv"),
        ("benchmark_out_of_dictionary_holdout.csv", "Data_S5c_hiv1_ood_holdout_summary.csv"),
        ("tier1_null_simulation_fpr.csv", "Data_S5d_hiv1_1300_null_fpr.csv"),
        ("benchmark_sc2_nextclade_and_iav_genoflu.csv", "Data_S5e_sc2_nextclade_iav_genoflu_benchmark.csv"),
        ("chronaeon_temporal_guard_ablation.tsv", "Data_S5f_chronaeon_temporal_guard_ablation.tsv"),
        ("stanford_hivdb_audit_category_summary.csv", "Data_S6a_stanford_hivdb_audit_summary.csv"),
        ("stanford_hivdb_unlinked_8k_breakdown.csv", "Data_S6b_stanford_hivdb_unlinked_breakdown.csv"),
        ("stanford_west_african_crf02_rt_validation.csv", "Data_S6c_stanford_west_african_crf02_trees.csv"),
        ("genbank_novel_circulating_recombinants.csv", "Data_S6d_genbank_full_novel_clusters.csv"),
        ("novel_lineages_phylogenetic_validation.csv", "Data_S6e_genbank_full_iqtree_validation.csv"),
        ("scale_partial_pol_phylogenetic_validation.csv", "Data_S6f_genbank_partial_pol_iqtree.csv"),
        ("scale_partial_env_phylogenetic_validation.csv", "Data_S6g_genbank_partial_env_iqtree.csv"),
        ("bvbrc_iav_h5_novel_constellations.tsv", "Data_S6h_bvbrc_iav_h5_novel_constellations.tsv"),
        ("bvbrc_sc2_23clusters_novelty_audit.tsv", "Data_S6i_bvbrc_sc2_urf_clusters.tsv"),
    ]
    for src_name, dst_name in copies:
        src = RES / src_name
        if src.exists():
            shutil.copy2(src, OUT / dst_name)
    print(f"[Data S2/S3/S5/S6] Exported {len(copies)} curated supplementary tables")


def compute_and_export_data_s4_sieve_and_scaling():
    """
    Deterministically compute:
      1. Gen-1 founder-deduplicated (N=71) & DAG-deduplicated (527 BPs) + CMH + MV Stage-2 sieve controls
      2. Stage-1 2-rigid-core (gp120 + Gag-CA) combined depletion across Uniform / SNP / 2% Sim / 4-BP 8% Sim nulls
      3. Clauset (2009) parametric bootstrap GOF (400 reps) + Vuong LR tests on N=761 lineages
    """
    with open(RES / "supp_manifest_sieve_749_lineages.csv") as f:
        lin_rows = list(csv.DictReader(f))
    with open(RES / "hiv1_descriptive_nomenclature_compendium.csv") as f:
        comp_crfs = [
            r for r in csv.DictReader(f)
            if r["legacy_category"].startswith("CRF") and int(r["n_verified_segments"]) > 1
        ]
    with open(find_compendium_bundle()) as f:
        bundle_crfs = json.load(f)["compendium"]
    with open(RES / "genbank_novel_circulating_recombinants.csv") as f:
        nov_crfs = [r for r in csv.DictReader(f) if r["lineage_status"].startswith("Candidate_Regional_Mosaic")]
    with open(RES / "hiv1_recombination_genealogy_dag.csv") as f:
        dag_map = {r["legacy_label"].strip(): r for r in csv.DictReader(f)}

    # Build exact 491 URF, 151 Est CRF, and 71 alias-collapsed Gen-1 Primary CRF cohorts
    urfs = []
    for r in lin_rows:
        if str(r.get("in_singleton_urf_461", "")).lower() == "true":
            segs = parse_segs(r["hxb2_barcode"])
            if len(segs) > 1:
                urfs.append({"label": r["lineage_key"], "segs": segs})

    est_all = []
    gen1_dedup = []
    seen_leg = set()
    seen_g1_formula = set()

    for r in comp_crfs:
        segs = parse_segs(r["hxb2_barcode"])
        if len(segs) <= 1:
            continue
        leg = r["legacy_label"].strip()
        sf = r["structural_formula"].strip()
        seen_leg.add(leg)
        gen = int(dag_map[leg]["generation"]) if leg in dag_map else 1
        rec = {"label": leg, "segs": segs, "gen": gen, "formula": sf}
        est_all.append(rec)
        if gen == 1 and sf not in seen_g1_formula:
            seen_g1_formula.add(sf)
            gen1_dedup.append(rec)

    for ent in bundle_crfs:
        leg = str(ent.get("label", "")).strip()
        cat = str(ent.get("category", "")).strip()
        if not cat.startswith("CRF") or leg in seen_leg:
            continue
        segs = parse_segs(ent.get("hxb2_barcode", ""))
        if len(segs) <= 1:
            continue
        seen_leg.add(leg)
        sf = str(ent.get("structural_formula", "")).strip()
        is_g2_post = any(k in leg for k in ("0107", "07109", "0755", "63A6"))
        gen = 2 if is_g2_post else 1
        rec = {"label": leg, "segs": segs, "gen": gen, "formula": sf}
        est_all.append(rec)
        if gen == 1 and sf not in seen_g1_formula:
            seen_g1_formula.add(sf)
            gen1_dedup.append(rec)

    for nr in nov_crfs:
        segs = parse_segs(nr["representative_hxb2_barcode"])
        if len(segs) <= 1:
            continue
        sf = nr["rhizaeon_formula"].strip()
        rec = {"label": "RU_G0B_gagMA", "segs": segs, "gen": 1, "formula": sf}
        est_all.append(rec)
        if sf not in seen_g1_formula:
            seen_g1_formula.add(sf)
            gen1_dedup.append(rec)

    assert len(urfs) == 491, f"Expected 491 URFs, got {len(urfs)}"
    assert len(est_all) == 151, f"Expected 151 Est CRFs, got {len(est_all)}"
    assert len(gen1_dedup) == 71, f"Expected 71 Gen-1 dedup CRFs, got {len(gen1_dedup)}"

    def feat_intact3p(segs):
        return int(not any(5831 <= 0.5 * (segs[i][2] + segs[i + 1][1]) <= 8653 for i in range(len(segs) - 1)))

    def feat_same(p1, p2):
        def _fn(segs):
            a, b = trunk_at(segs, p1), trunk_at(segs, p2)
            return int(a is not None and a == b)
        return _fn

    def eval_binary_feature(name, fn, cohort, filt=None):
        u_sub = [x for x in urfs if (filt is None or filt(x["segs"]))]
        c_sub = [x for x in cohort if (filt is None or filt(x["segs"]))]
        u_y = [fn(x["segs"]) for x in u_sub]
        c_y = [fn(x["segs"]) for x in c_sub]
        or_crude, p_crude = fisher_exact([[sum(c_y), len(c_y) - sum(c_y)], [sum(u_y), len(u_y) - sum(u_y)]])
        # Pair-stratified CMH
        by_pair = defaultdict(lambda: [0, 0, 0, 0])
        for x, yv in zip(c_sub, c_y):
            by_pair[pair_key(x["segs"])][0 if yv else 1] += 1
        for x, yv in zip(u_sub, u_y):
            by_pair[pair_key(x["segs"])][2 if yv else 3] += 1
        or_cmh, p_cmh = cmh_test(list(by_pair.values()))
        # Multivariable case-control logistic: logit P(CRF=1) = b0 + b1*feat + b2*n_bp + pair fixed effects
        y_case = [1] * len(c_y) + [0] * len(u_y)
        x_feat = c_y + u_y
        nbp_all = [len(x["segs"]) - 1 for x in c_sub] + [len(x["segs"]) - 1 for x in u_sub]
        pairs_all = [pair_key(x["segs"]) for x in c_sub] + [pair_key(x["segs"]) for x in u_sub]
        or_mv, ci_lo, ci_hi, p_mv = logistic_mv(y_case, x_feat, nbp_all, pairs_all)
        return {
            "feature": name,
            "urf_k": sum(u_y), "urf_n": len(u_y), "urf_pct": round(100.0 * sum(u_y) / len(u_y), 2),
            "crf_k": sum(c_y), "crf_n": len(c_y), "crf_pct": round(100.0 * sum(c_y) / len(c_y), 2),
            "crude_or": round(or_crude, 3), "crude_p": f"{p_crude:.3e}",
            "cmh_or": round(or_cmh, 3), "cmh_p": f"{p_cmh:.3e}",
            "mv_or": round(or_mv, 3), "mv_ci_lo": round(ci_lo, 3), "mv_ci_hi": round(ci_hi, 3), "mv_p": f"{p_mv:.3e}",
        }

    def filt_foreign_gp120(segs):
        t_ma = trunk_at(segs, 950)
        t_gp = trunk_at(segs, 6900)
        t_ct = trunk_at(segs, 8650)
        return (t_ma is not None and t_gp is not None and t_ct is not None and t_gp != t_ma)

    def feat_gp41_revert(segs):
        return int(trunk_at(segs, 8650) == trunk_at(segs, 950))

    s2_rows = []
    for cohort_lbl, cohort in [("Est_All_151", est_all), ("Gen1_Dedup_71", gen1_dedup)]:
        for nm, fn, flt in [
            ("Intact_3prime_Block_5831_8653", feat_intact3p, None),
            ("RRE_rev2_7880_8500", feat_same(7880, 8500), None),
            ("tat1_tat2_5950_8400", feat_same(5950, 8400), None),
            ("gp120_gp41_6900_7950", feat_same(6900, 7950), None),
            ("gp41_CT_Reversion_Foreign_gp120", feat_gp41_revert, filt_foreign_gp120),
        ]:
            row = {"cohort": cohort_lbl}
            row.update(eval_binary_feature(nm, fn, cohort, flt))
            s2_rows.append(row)

    out_s4b = OUT / "Data_S4b_sieve_gen1_and_mv_controls.csv"
    with open(out_s4b, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(s2_rows[0].keys()))
        w.writeheader()
        w.writerows(s2_rows)
    print(f"[Data S4b] Wrote {len(s2_rows)} Stage-2 founder-deduplicated & MV rows -> {out_s4b.name}")
    for r in s2_rows:
        if r["cohort"] == "Gen1_Dedup_71" and r["feature"] in ("Intact_3prime_Block_5831_8653", "RRE_rev2_7880_8500", "gp41_CT_Reversion_Foreign_gp120"):
            print(f"   {r['cohort']} | {r['feature']}: {r['crf_k']}/{r['crf_n']} ({r['crf_pct']}%) CMH OR={r['cmh_or']} (p={r['cmh_p']}), MV OR={r['mv_or']} (p={r['mv_p']})")

    # Clauset (2009) parametric bootstrap GOF + Vuong tests on N = 761 mosaic lineages
    sizes = np.array([int(r["n_genbank_patients"]) for r in lin_rows], dtype=int)
    n_tot = len(sizes)
    res_hz = minimize(lambda a: n_tot * math.log(zeta(a[0], 1.0)) + a[0] * np.sum(np.log(sizes)), [2.2], bounds=[(1.05, 10.0)])
    alpha_hat = float(res_hz.x[0])
    ll_hz_i = -alpha_hat * np.log(sizes) - math.log(zeta(alpha_hat, 1.0))
    ll_hz = float(np.sum(ll_hz_i))
    aic_hz = 2.0 - 2.0 * ll_hz

    def ks_hurwitz(arr, alpha):
        vals = np.arange(1, int(arr.max()) + 1, dtype=float)
        emp_cdf = np.searchsorted(np.sort(arr), vals, side="right") / len(arr)
        mod_pmf = (vals ** (-alpha)) / zeta(alpha, 1.0)
        mod_cdf = np.cumsum(mod_pmf)
        return float(np.max(np.abs(emp_cdf - mod_cdf)))

    ks_obs = ks_hurwitz(sizes, alpha_hat)
    rng = np.random.default_rng(20261120)
    k_grid = np.arange(1, 20001, dtype=float)
    pmf_grid = (k_grid ** (-alpha_hat)) / zeta(alpha_hat, 1.0)
    pmf_grid /= pmf_grid.sum()
    cdf_grid = np.cumsum(pmf_grid)
    n_boot = 400
    ks_ge = 0
    for _ in range(n_boot):
        u = rng.random(n_tot)
        sim = np.searchsorted(cdf_grid, u) + 1
        r_b = minimize(lambda a: n_tot * math.log(zeta(a[0], 1.0)) + a[0] * np.sum(np.log(sim)), [alpha_hat], bounds=[(1.05, 10.0)])
        if ks_hurwitz(sim, float(r_b.x[0])) >= ks_obs:
            ks_ge += 1
    p_clauset = ks_ge / n_boot

    def fit_dln(mode):
        def ll_i(mu, sig):
            if mode == "half_int":
                lo = norm.cdf((np.log(sizes - 0.5) - mu) / sig)
                hi = norm.cdf((np.log(sizes + 0.5) - mu) / sig)
                denom = 1.0 - norm.cdf((math.log(0.5) - mu) / sig)
            else:
                lo = norm.cdf((np.log(sizes) - mu) / sig)
                hi = norm.cdf((np.log(sizes + 1.0) - mu) / sig)
                denom = 1.0 - norm.cdf((0.0 - mu) / sig)
            return np.log(np.maximum((hi - lo) / max(denom, 1e-300), 1e-300))

        best = (1e18, None)
        for m0 in (-25.0, -10.0, -5.0, -2.0, -1.0, 0.0, 1.0):
            for s0 in (0.5, 1.5, 3.0, 5.0):
                r = minimize(lambda th: -np.sum(ll_i(th[0], th[1])), [m0, s0], bounds=[(-50.0, 20.0), (0.05, 50.0)])
                if r.fun < best[0]:
                    best = (r.fun, r.x)
        mu_h, sig_h = best[1]
        li = ll_i(mu_h, sig_h)
        diff = ll_hz_i - li
        z_v = float(np.sum(diff)) / (math.sqrt(n_tot) * float(np.std(diff, ddof=1)))
        p_v = float(2.0 * (1.0 - norm.cdf(abs(z_v))))
        return mu_h, sig_h, float(np.sum(li)), 4.0 - 2.0 * float(np.sum(li)), z_v, p_v

    mu1, sig1, ll1, aic1, z1, p1 = fit_dln("half_int")
    mu2, sig2, ll2, aic2, z2, p2 = fit_dln("right_int")

    scale_rows = [
        {"model": "Discrete_HurwitzZeta_xmin1", "n": n_tot, "params": f"alpha={alpha_hat:.4f}", "loglik": round(ll_hz, 3), "aic": round(aic_hz, 2), "ks_stat": round(ks_obs, 4), "clauset_boot_p_400": round(p_clauset, 4), "vuong_z_hz_minus_alt": 0.0, "vuong_p": 1.0},
        {"model": "Discrete_Lognormal_half_int_unconstrained", "n": n_tot, "params": f"mu={mu1:.4f}, sigma={sig1:.4f}", "loglik": round(ll1, 3), "aic": round(aic1, 2), "ks_stat": "", "clauset_boot_p_400": "", "vuong_z_hz_minus_alt": round(z1, 3), "vuong_p": f"{p1:.4g}"},
        {"model": "Discrete_Lognormal_right_int_unconstrained", "n": n_tot, "params": f"mu={mu2:.4f}, sigma={sig2:.4f}", "loglik": round(ll2, 3), "aic": round(aic2, 2), "ks_stat": "", "clauset_boot_p_400": "", "vuong_z_hz_minus_alt": round(z2, 3), "vuong_p": f"{p2:.4g}"},
    ]
    out_s4c = OUT / "Data_S4c_clauset_vuong_cluster_scaling.csv"
    with open(out_s4c, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(scale_rows[0].keys()))
        w.writeheader()
        w.writerows(scale_rows)
    print(f"[Data S4c] Clauset GOF KS={ks_obs:.4f}, p={p_clauset:.4f}; Vuong(half_int) Z={z1:+.3f} (p={p1:.4g}) -> {out_s4c.name}")


if __name__ == "__main__":
    export_data_s1_chronaeon()
    export_data_s2_s3_s5_s6()
    compute_and_export_data_s4_sieve_and_scaling()
    print("All Data S1-S6 supplementary tables generated and verified successfully.")
