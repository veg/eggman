#!/usr/bin/env python3
"""
Script 10: Quantifying (1) the Rate of Viable URF Generation & CRF Establishment
and (2) the Epistatic & Structural Sieve Separating Established Circulating Lineages
from Dead-End Singleton URFs across 7,720 GenBank + 5,606 LANL Full-Length HIV-1 Genomes.
"""

import json
import re
import math
from pathlib import Path
import numpy as np
import pandas as pd
from collections import Counter, defaultdict
from scipy import optimize, special, stats
import statsmodels.api as sm
from statsmodels.stats.contingency_tables import StratifiedTable
from statsmodels.stats.multitest import multipletests


def parse_barcode_map(barcode_str):
    """
    Parse an hxb2_barcode string like 'C:790-5669|B:5670-6223|?E:6224-9417'
    or '790-2157:A1|2158-3410:G0' into a list of (start_hxb2, end_hxb2, trunk) tuples.
    """
    if not isinstance(barcode_str, str) or not barcode_str.strip():
        return []
    segs = []
    for part in barcode_str.split("|"):
        if ":" not in part:
            continue
        lhs, rhs = part.split(":", 1)
        if "-" in rhs and rhs.replace("-", "").isdigit():
            trk = lhs.strip()
            s_str, e_str = rhs.split("-", 1)
        elif "-" in lhs and lhs.replace("-", "").isdigit():
            trk = rhs.strip()
            s_str, e_str = lhs.split("-", 1)
        else:
            continue
        try:
            segs.append((int(s_str), int(e_str), trk))
        except ValueError:
            continue
    return segs


def trunk_at_pos(segs, hxb2_pos):
    """Return the assigned parental trunk at a specific HXB2 coordinate."""
    for s, e, trk in segs:
        if s <= hxb2_pos <= e:
            return trk
    # If between segments, pick nearest
    best_trk = None
    best_dist = 10**9
    for s, e, trk in segs:
        d = min(abs(hxb2_pos - s), abs(hxb2_pos - e))
        if d < best_dist:
            best_dist = d
            best_trk = trk
    return best_trk


def get_internal_breakpoints(segs):
    """Return list of internal breakpoint HXB2 coordinates (midpoints between adjacent segments)."""
    bps = []
    for i in range(len(segs) - 1):
        bp = int(round(0.5 * (segs[i][1] + segs[i + 1][0])))
        bps.append((bp, segs[i][2], segs[i + 1][2]))
    return bps


def get_trunk_proportions(segs):
    """Return dict of trunk -> fraction of coding length."""
    counts = Counter()
    total = 0
    for s, e, trk in segs:
        length = max(0, e - s + 1)
        counts[trk] += length
        total += length
    if total == 0:
        return {}
    return {k: v / total for k, v in counts.items()}


def main():
    print("=" * 88)
    print("HIV-1 RECOMBINATION VIABILITY RATE & EPISTATIC ESTABLISHMENT SIEVE")
    print("=" * 88)

    gb_df = pd.read_csv("results/genbank_31k_rhizaeon_classifications.csv")
    lanl_df = pd.read_csv("results/lanl_all_full_genomes_rhizaeon_validation.csv")
    comp_df = pd.read_csv("results/hiv1_descriptive_nomenclature_compendium.csv")
    novel_df = pd.read_csv("results/genbank_novel_circulating_recombinants.csv")
    dag_df = pd.read_csv("results/hiv1_recombination_genealogy_dag.csv")

    valid_gb = gb_df[gb_df["status"] == "Validated_GroupM"].copy()
    valid_gb["country"] = valid_gb["country"].fillna("Unknown").astype(str).str.replace(" nonfunctional", "", regex=False).str.strip()

    # Canonicalize updated A6->A1 refined formulas in novel_df by matching member_accessions
    acc_to_canon_form = {}
    for _, nr in novel_df.iterrows():
        for acc in str(nr["member_accessions"]).split(","):
            acc = acc.strip()
            if acc:
                acc_to_canon_form[acc] = nr["rhizaeon_formula"]
    valid_gb["canon_formula"] = valid_gb.apply(
        lambda r: acc_to_canon_form.get(r["accession"], r["rhizaeon_formula"])
        if r["recomb_category"] == "Unidentified_Novel_Recombinant" else r["rhizaeon_formula"],
        axis=1
    )

    print(f"Loaded {len(valid_gb):,} patient-deduplicated full-length GenBank genomes.")
    print(f"Loaded {len(lanl_df):,} LANL full-length genomes and {len(comp_df):,} reference compendium entries.")

    # =========================================================================
    # PART 1A: OFFSPRING CLUSTER-SIZE DISTRIBUTION & ESTABLISHMENT PROBABILITY
    # =========================================================================
    print("\n" + "-" * 88)
    print("PART 1A: OFFSPRING CLUSTER-SIZE DISTRIBUTION & ESTABLISHMENT PROBABILITY")
    print("-" * 88)

    recomb_gb = valid_gb[valid_gb["n_segments"] > 1].copy()
    clonal_gb = valid_gb[valid_gb["n_segments"] == 1].copy()

    verified_novel_map = dict(zip(novel_df["rhizaeon_formula"], novel_df["n_homologous_bp_patients"]))
    novel_verdict_map = dict(zip(novel_df["rhizaeon_formula"], novel_df["lineage_status"]))

    lineage_counts = {}
    lineage_meta = {}

    # Known CRFs observed in GenBank
    known_mask = recomb_gb["recomb_category"].isin([
        "Known_Indexed_CRF_Exact",
        "Known_Indexed_CRF_LANLTag",
        "Known_Indexed_CRF_Title",
        "Known_Indexed_CRF_Post2023",
        "Attenuated_or_Boundary_Variant_of_Known_CRF"
    ])
    for lin, grp in recomb_gb[known_mask].groupby("matched_known_lineage"):
        if lin in ("Clonal_Basal_Subtype", "None", "") or pd.isna(lin):
            continue
        # Normalize matched lineage name (e.g. CRF01_AE / ?E -> CRF01_AE)
        norm_lin = str(lin).split("/")[0].strip()
        lineage_counts[f"KNOWN:{norm_lin}"] = lineage_counts.get(f"KNOWN:{norm_lin}", 0) + len(grp)
        if f"KNOWN:{norm_lin}" not in lineage_meta:
            rep = grp.iloc[0]
            lineage_meta[f"KNOWN:{norm_lin}"] = {
                "type": "Known_CRF",
                "rep_accession": str(rep["accession"]),
                "formula": rep["rhizaeon_formula"],
                "barcode": rep["hxb2_barcode"],
                "n_segments": int(rep["n_segments"]),
                "n_breakpoints": int(rep["n_breakpoints"]),
            }

    novel_mask = recomb_gb["recomb_category"] == "Unidentified_Novel_Recombinant"
    for form, grp in recomb_gb[novel_mask].groupby("canon_formula"):
        raw_n = len(grp)
        status = novel_verdict_map.get(form, "")
        if form in verified_novel_map:
            v_n = int(verified_novel_map[form])
        else:
            v_n = raw_n

        if "Convergent_Hotspot" in str(status) or "Same_Patient_Clones" in str(status):
            for idx, (_, row) in enumerate(grp.iterrows()):
                key = f"URF_SPLIT:{form}#{idx}"
                lineage_counts[key] = 1
                lineage_meta[key] = {
                    "type": "Singleton_URF",
                    "rep_accession": str(row["accession"]),
                    "formula": form,
                    "barcode": row["hxb2_barcode"],
                    "n_segments": int(row["n_segments"]),
                    "n_breakpoints": int(row["n_breakpoints"]),
                }
        else:
            key = f"NOVEL:{form}"
            lineage_counts[key] = max(1, v_n)
            rep = grp.iloc[0]
            lineage_meta[key] = {
                "type": "Novel_Circulating" if v_n >= 3 else ("Doubleton_URF" if v_n == 2 else "Singleton_URF"),
                "rep_accession": str(rep["accession"]),
                "formula": form,
                "barcode": rep["hxb2_barcode"],
                "n_segments": int(rep["n_segments"]),
                "n_breakpoints": int(rep["n_breakpoints"]),
            }

    sizes = np.array(list(lineage_counts.values()), dtype=int)
    n_total_emitted_lineages = len(sizes)
    n_n1 = int(np.sum(sizes == 1))
    n_n2 = int(np.sum(sizes == 2))
    n_ge3 = int(np.sum(sizes >= 3))
    n_ge5 = int(np.sum(sizes >= 5))
    n_ge10 = int(np.sum(sizes >= 10))
    n_ge50 = int(np.sum(sizes >= 50))

    # Exact Discrete Hurwitz-Zeta Power-Law MLE (xmin=1 and CSN KS-optimal xmin) + Discrete Competing Models
    def _fit_discrete_powerlaw(x_arr, xmin):
        xt = x_arr[x_arr >= xmin].astype(float)
        n_t = len(xt)
        sum_log = float(np.sum(np.log(xt)))
        def nll(a):
            if a <= 1.02:
                return 1e18
            return n_t * math.log(float(special.zeta(a, xmin))) + a * sum_log
        res = optimize.minimize_scalar(nll, bounds=(1.05, 8.0), method="bounded")
        alpha = float(res.x)
        ll = -float(res.fun)
        eps = 1e-4
        d2 = (nll(alpha + eps) - 2.0 * nll(alpha) + nll(alpha - eps)) / (eps ** 2)
        se = 1.0 / math.sqrt(max(d2, 1e-12))
        ks = np.nan
        try:
            supp = np.arange(xmin, int(max(xt)) + 500, dtype=float)
            pmf = (supp ** (-alpha)) / float(special.zeta(alpha, xmin))
            cdf_th = np.cumsum(pmf)
            vals = np.arange(xmin, int(max(xt)) + 1)
            emp_cdf = np.array([np.mean(xt <= v) for v in vals])
            th_at_vals = cdf_th[:len(vals)]
            ks = float(np.max(np.abs(emp_cdf - th_at_vals)))
        except Exception:
            pass
        return alpha, se, ll, n_t, ks

    alpha_mle, alpha_se, ll_pl1, _, ks_pl1 = _fit_discrete_powerlaw(sizes, 1)
    best_csn = (1, alpha_mle, alpha_se, ll_pl1, len(sizes), ks_pl1)
    for xm in sorted(set(int(v) for v in sizes)):
        if xm < 1 or int(np.sum(sizes >= xm)) < 30:
            continue
        a_x, se_x, ll_x, nt_x, ks_x = _fit_discrete_powerlaw(sizes, xm)
        if np.isfinite(ks_x) and ks_x < best_csn[5]:
            best_csn = (xm, a_x, se_x, ll_x, nt_x, ks_x)

    # Compare Discrete Power-Law vs Yule-Simon, Shifted Negative Binomial, & Discrete Lognormal on x >= 1
    x_int = sizes.astype(int)
    n_sz = len(x_int)
    supp_full = np.arange(1, max(10000, int(x_int.max()) * 2) + 1, dtype=float)
    def nll_yule(rho):
        if rho <= 0.01:
            return 1e18
        lp = math.log(rho) + special.betaln(x_int, rho + 1.0)
        return -float(np.sum(lp))
    res_y = optimize.minimize_scalar(nll_yule, bounds=(0.05, 10.0), method="bounded")
    ll_yule = -float(res_y.fun)

    def nll_nb(params):
        log_r, logit_p = params
        r = math.exp(log_r)
        p = 1.0 / (1.0 + math.exp(-logit_p))
        k = x_int - 1
        lp = special.gammaln(k + r) - special.gammaln(r) - special.gammaln(k + 1.0) + r * math.log(max(p, 1e-12)) + k * math.log(max(1.0 - p, 1e-12))
        return -float(np.sum(lp))
    res_nb = optimize.minimize(nll_nb, [math.log(0.2), math.log(0.2 / 0.8)], method="Nelder-Mead")
    ll_nb = -float(res_nb.fun)
    r_nb = math.exp(res_nb.x[0])
    p_nb = 1.0 / (1.0 + math.exp(-res_nb.x[1]))

    def nll_dln(params):
        mu, log_sig = params
        sig = math.exp(log_sig)
        if sig < 0.05 or sig > 10.0:
            return 1e18
        z_s = -((np.log(supp_full) - mu) ** 2) / (2.0 * sig * sig) - np.log(supp_full)
        z_m = float(np.max(z_s))
        l_norm = z_m + math.log(float(np.sum(np.exp(z_s - z_m))))
        l_x = -((np.log(x_int) - mu) ** 2) / (2.0 * sig * sig) - np.log(x_int)
        return -float(np.sum(l_x - l_norm))
    res_dln = optimize.minimize(nll_dln, [-1.0, math.log(1.5)], method="Nelder-Mead")
    ll_dln = -float(res_dln.fun)
    mu_dln, sig_dln = float(res_dln.x[0]), math.exp(res_dln.x[1])

    fit_models = [
        {"model": "Discrete_PowerLaw_HurwitzZeta_xmin1", "xmin": 1, "n_tail": n_sz, "params": f"alpha={alpha_mle:.4f} (SE={alpha_se:.4f})", "loglik": round(ll_pl1, 3), "k_params": 1, "aic": round(2 * 1 - 2 * ll_pl1, 2), "ks_stat": round(ks_pl1, 4)},
        {"model": "Discrete_YuleSimon_xmin1", "xmin": 1, "n_tail": n_sz, "params": f"rho={float(res_y.x):.4f}", "loglik": round(ll_yule, 3), "k_params": 1, "aic": round(2 * 1 - 2 * ll_yule, 2), "ks_stat": np.nan},
        {"model": "Shifted_NegBinomial_xmin1", "xmin": 1, "n_tail": n_sz, "params": f"r={r_nb:.4f}, p={p_nb:.4f}", "loglik": round(ll_nb, 3), "k_params": 2, "aic": round(2 * 2 - 2 * ll_nb, 2), "ks_stat": np.nan},
        {"model": "Discrete_Lognormal_xmin1", "xmin": 1, "n_tail": n_sz, "params": f"mu={mu_dln:.4f}, sigma={sig_dln:.4f}", "loglik": round(ll_dln, 3), "k_params": 2, "aic": round(2 * 2 - 2 * ll_dln, 2), "ks_stat": np.nan},
        {"model": f"Discrete_PowerLaw_CSN_Opt_xmin{best_csn[0]}", "xmin": best_csn[0], "n_tail": best_csn[4], "params": f"alpha={best_csn[1]:.4f} (SE={best_csn[2]:.4f})", "loglik": round(best_csn[3], 3), "k_params": 1, "aic": round(2 * 1 - 2 * best_csn[3], 2), "ks_stat": round(best_csn[5], 4)},
    ]
    min_aic_x1 = min(m["aic"] for m in fit_models if m["xmin"] == 1)
    for m_fit in fit_models:
        m_fit["delta_aic_vs_best_xmin1"] = round(m_fit["aic"] - min_aic_x1, 2) if m_fit["xmin"] == 1 else np.nan
    pd.DataFrame(fit_models).to_csv("results/sieve_offspring_distribution_fits.csv", index=False)

    print(f"Total distinct viable transmitted recombinant lineages in GenBank (N={len(recomb_gb):,} rec genomes): {n_total_emitted_lineages}")
    print(f"  - Dead-End Singleton URFs (N = 1 patient):       {n_n1:4d} / {n_total_emitted_lineages} ({100*n_n1/n_total_emitted_lineages:.1f}%)")
    print(f"  - Stuttering Dyads/Doubletons (N = 2 patients):  {n_n2:4d} / {n_total_emitted_lineages} ({100*n_n2/n_total_emitted_lineages:.1f}%)")
    print(f"  - Established Circulating (N >= 3 patients):     {n_ge3:4d} / {n_total_emitted_lineages} ({100*n_ge3/n_total_emitted_lineages:.1f}%)")
    print(f"  - Multi-Cluster Regional (N >= 5 patients):      {n_ge5:4d} / {n_total_emitted_lineages} ({100*n_ge5/n_total_emitted_lineages:.1f}%)")
    print(f"  - Sustained Epidemic CRFs (N >= 10 patients):    {n_ge10:4d} / {n_total_emitted_lineages} ({100*n_ge10/n_total_emitted_lineages:.1f}%)")
    print(f"  - Major Pandemic CRFs (N >= 50 patients):        {n_ge50:4d} / {n_total_emitted_lineages} ({100*n_ge50/n_total_emitted_lineages:.1f}%)")
    print(f"  - Discrete Hurwitz-Zeta Power-Law MLE (xmin=1):  alpha={alpha_mle:.3f} +/- {alpha_se:.3f} (KS={ks_pl1:.4f}, AIC={2 - 2*ll_pl1:.1f}; CSN xmin={best_csn[0]} alpha={best_csn[1]:.3f}, KS={best_csn[5]:.4f})")

    # Conditional survival probabilities across the transmission funnel:
    p_2_given_1 = np.sum(sizes >= 2) / n_total_emitted_lineages
    p_3_given_2 = np.sum(sizes >= 3) / np.sum(sizes >= 2)
    p_10_given_3 = np.sum(sizes >= 10) / np.sum(sizes >= 3)
    print(f"\nConditional Transmission Sieve Probabilities:")
    print(f"  - P(N >= 2 | Viable Transmitted URF N >= 1) = {100*p_2_given_1:.2f}% ({np.sum(sizes>=2)}/{n_total_emitted_lineages})")
    print(f"  - P(N >= 3 | Reached Dyad N >= 2)           = {100*p_3_given_2:.2f}% ({n_ge3}/{np.sum(sizes>=2)})")
    print(f"  - P(N >= 10 | Established CRF N >= 3)       = {100*p_10_given_3:.2f}% ({n_ge10}/{n_ge3})")

    # =========================================================================
    # PART 1B: REGIONAL HETEROZYGOSITY & NORMALIZED URF GENERATION RATE
    # =========================================================================
    print("\n" + "-" * 88)
    print("PART 1B: REGIONAL SUBTYPE HETEROZYGOSITY & URF GENERATION RATE")
    print("-" * 88)

    ssa_countries = {
        "South Africa", "Cameroon", "Uganda", "Tanzania", "Democratic Republic of the Congo",
        "Kenya", "Ghana", "Nigeria", "Botswana", "Rwanda", "Zambia"
    }

    reg_rows = []
    for country, c_df in valid_gb.groupby("country"):
        if country in ("Unknown", "", "None") or pd.isna(country) or len(c_df) < 40:
            continue
        n_pat = len(c_df)
        # Strictly exclude ALL unidentified/putative URFs from H_circ to prevent mathematical circularity;
        # H_circ measures diversity purely over Clonal Subtypes + Established Known CRFs.
        est_df = c_df[c_df["recomb_category"] != "Unidentified_Novel_Recombinant"]

        # 1) Basal trunk pool frequency in this country (from established lineages in est_df)
        trunk_nt = Counter()
        for bc in est_df["hxb2_barcode"]:
            for s, e, trk in parse_barcode_map(bc):
                trunk_nt[trk] += (e - s + 1)
        tot_nt = sum(trunk_nt.values())
        h_trunk = 1.0 - sum((v / tot_nt) ** 2 for v in trunk_nt.values()) if tot_nt > 0 else 0.0

        # 2) Distinct established circulating lineage diversity (Pure Clonal Subtypes + Known CRFs only)
        lin_labels = []
        for _, r in est_df.iterrows():
            if r["n_segments"] == 1:
                lin_labels.append(f"Pure_{r['backbone_trunk']}")
            else:
                lin_labels.append(str(r["matched_known_lineage"]).split("/")[0].strip())
        lin_counts = Counter(lin_labels)
        tot_lin = sum(lin_counts.values())
        h_circ = 1.0 - sum((v / tot_lin) ** 2 for v in lin_counts.values()) if tot_lin > 0 else 0.0

        n_novel_rec = int((c_df["recomb_category"] == "Unidentified_Novel_Recombinant").sum())
        n_novel_arch = c_df[c_df["recomb_category"] == "Unidentified_Novel_Recombinant"]["canon_formula"].nunique()
        yrs = c_df["year"].dropna()
        yr_span = float(yrs.max() - yrs.min()) if len(yrs) >= 2 else 0.0
        is_ssa = int(country in ssa_countries)

        f_urf_patients = n_novel_rec / n_pat
        f_urf_arch = n_novel_arch / n_pat
        mu_urf_norm = f_urf_arch / h_trunk if h_trunk >= 0.10 else np.nan

        top_trunks = ", ".join(f"{k}:{100*v/tot_nt:.0f}%" for k, v in trunk_nt.most_common(3))
        reg_rows.append({
            "country": country,
            "n_patients": n_pat,
            "n_est_non_urf_patients": len(est_df),
            "yr_span": yr_span,
            "is_ssa": is_ssa,
            "h_trunk_heterozygosity": round(h_trunk, 4),
            "h_circ_lineage_diversity": round(h_circ, 4),
            "h_circ_non_urf": round(h_circ, 4),
            "n_de_novo_urf_patients": n_novel_rec,
            "n_distinct_urf_architectures": n_novel_arch,
            "f_urf_patient_prev": round(f_urf_patients, 4),
            "f_urf_arch_per_patient": round(f_urf_arch, 4),
            "normalized_urf_rate_per_H": round(mu_urf_norm, 4) if not np.isnan(mu_urf_norm) else np.nan,
            "dominant_parental_pool": top_trunks
        })

    reg_df = pd.DataFrame(reg_rows).sort_values("n_patients", ascending=False)
    reg_df.to_csv("results/regional_heterozygosity_urf_rates.csv", index=False)

    valid_reg = reg_df.copy()
    r_p_circ, p_p_circ = stats.pearsonr(valid_reg["h_circ_non_urf"], valid_reg["f_urf_arch_per_patient"])
    r_s_circ, p_s_circ = stats.spearmanr(valid_reg["h_circ_non_urf"], valid_reg["f_urf_arch_per_patient"])
    r_p_trk, p_p_trk = stats.pearsonr(valid_reg["h_trunk_heterozygosity"], valid_reg["f_urf_arch_per_patient"])
    r_s_trk, p_s_trk = stats.spearmanr(valid_reg["h_trunk_heterozygosity"], valid_reg["f_urf_arch_per_patient"])
    w = valid_reg["n_patients"].values
    H_circ_val = valid_reg["h_circ_non_urf"].values
    F_val = valid_reg["f_urf_arch_per_patient"].values
    beta_weighted = np.sum(w * H_circ_val * F_val) / np.sum(w * (H_circ_val ** 2))

    # Multivariable regressions controlling for log10(N), sampling year span, and Sub-Saharan Africa indicator
    valid_reg["log10_n"] = np.log10(valid_reg["n_patients"])
    X_reg = sm.add_constant(valid_reg[["h_circ_non_urf", "log10_n", "yr_span", "is_ssa"]])
    ols_mod = sm.OLS(valid_reg["f_urf_arch_per_patient"], X_reg).fit()
    wls_mod = sm.WLS(valid_reg["f_urf_arch_per_patient"], X_reg, weights=valid_reg["n_patients"]).fit()
    glm_mod = sm.GLM(
        np.column_stack([valid_reg["n_distinct_urf_architectures"], valid_reg["n_patients"] - valid_reg["n_distinct_urf_architectures"]]),
        X_reg,
        family=sm.families.Binomial()
    ).fit()
    glm_hc3 = sm.GLM(
        np.column_stack([valid_reg["n_distinct_urf_architectures"], valid_reg["n_patients"] - valid_reg["n_distinct_urf_architectures"]]),
        X_reg,
        family=sm.families.Binomial()
    ).fit(cov_type="HC3")
    y_urf = valid_reg["n_distinct_urf_architectures"].astype(float)
    nb_mod = sm.NegativeBinomial(y_urf, X_reg, offset=np.log(valid_reg["n_patients"].astype(float))).fit(disp=False, maxiter=500)

    mv_rows = []
    for mod_name, mod, extra in [
        ("OLS", ols_mod, ols_mod.rsquared),
        ("WLS_by_N", wls_mod, wls_mod.rsquared),
        ("Binomial_GLM", glm_mod, 1.0 - (glm_mod.deviance / glm_mod.null_deviance)),
        ("Binomial_GLM_HC3", glm_hc3, 1.0 - (glm_hc3.deviance / glm_hc3.null_deviance)),
        ("NegBin_GLM_Offset", nb_mod, nb_mod.prsquared),
    ]:
        for term in X_reg.columns:
            if term not in mod.params:
                continue
            mv_rows.append({
                "model": mod_name,
                "term": term,
                "coef": round(float(mod.params[term]), 5),
                "std_err": round(float(mod.bse[term]), 5),
                "stat": round(float(mod.tvalues[term]), 3),
                "p_value": float(mod.pvalues[term]),
                "r2_or_pseudo_r2": round(float(extra), 4),
                "n_countries": len(valid_reg),
            })
    pd.DataFrame(mv_rows).to_csv("results/sieve_hcirc_multivariate_regression.csv", index=False)

    print(f"Evaluated {len(reg_df)} countries with >= 40 sequenced patients (URFs strictly excluded from H_circ):")
    print(f"  - Correlation with Established Strain Heterozygosity (H_circ_non_urf):")
    print(f"    Pearson r = {r_p_circ:.3f} (p = {p_p_circ:.2e}), Spearman rho = {r_s_circ:.3f} (p = {p_s_circ:.2e})")
    print(f"  - Correlation with Basal Trunk Nucleotide Heterozygosity (H_trunk):")
    print(f"    Pearson r = {r_p_trk:.3f} (p = {p_p_trk:.2e}), Spearman rho = {r_s_trk:.3f} (p = {p_s_trk:.2e})")
    print(f"  - Sample-Size-Weighted Conversion Slope (beta_URF = f_URF / H_circ_non_urf): {beta_weighted:.4f} ({100*beta_weighted:.2f}% per unit H)")
    print(f"  - Multivariable H_circ_non_urf effect (adj for log10(N), year_span, is_SSA):")
    print(f"    OLS: beta={ols_mod.params['h_circ_non_urf']:.4f} (p={ols_mod.pvalues['h_circ_non_urf']:.3e}, R2={ols_mod.rsquared:.3f}); is_SSA={ols_mod.params['is_ssa']:.4f} (p={ols_mod.pvalues['is_ssa']:.3e})")
    print(f"    WLS: beta={wls_mod.params['h_circ_non_urf']:.4f} (p={wls_mod.pvalues['h_circ_non_urf']:.3e}, R2={wls_mod.rsquared:.3f}); is_SSA={wls_mod.params['is_ssa']:.4f} (p={wls_mod.pvalues['is_ssa']:.3e})")
    print(f"    GLM: logit beta={glm_mod.params['h_circ_non_urf']:.3f} (p={glm_mod.pvalues['h_circ_non_urf']:.3e}); is_SSA={glm_mod.params['is_ssa']:.3f} (p={glm_mod.pvalues['is_ssa']:.3e})")
    print("\nTop 12 Countries by Sample Size (Heterozygosity vs. De Novo URF Rate):")
    for _, r in reg_df.head(12).iterrows():
        norm_c = f"{100*(r['f_urf_arch_per_patient']/r['h_circ_lineage_diversity']):.1f}%" if r["h_circ_lineage_diversity"] >= 0.10 else "N/A (mono)"
        print(f"  {r['country']:15s} | N={r['n_patients']:4d} | H_circ={r['h_circ_lineage_diversity']:.2f} | H_trk={r['h_trunk_heterozygosity']:.2f} | "
              f"URFs={r['n_distinct_urf_architectures']:3d} ({100*r['f_urf_arch_per_patient']:4.1f}%) | "
              f"Rate/H_circ={norm_c:6s} | Pool: {r['dominant_parental_pool']}")

    # =========================================================================
    # PART 1C: TEMPORAL BIRTH RATE OF ESTABLISHED CRFs (CHRONAEON t_MRCA)
    # =========================================================================
    print("\n" + "-" * 88)
    print("PART 1C: TEMPORAL EMERGENCE RATE OF ESTABLISHED CRFs (ChronAeon t_MRCA)")
    print("-" * 88)

    crf_comp = comp_df[comp_df["legacy_category"].str.startswith("CRF") & (comp_df["n_verified_segments"] > 1)].copy()
    crf_dated = crf_comp[crf_comp["t_mrca"].notna()].copy()
    epochs = [
        ("1950-1969 (Pre-Pandemic Basal Hub)", 1950, 1970),
        ("1970-1984 (Early Epidemic Bridge)", 1970, 1985),
        ("1985-1999 (Global Co-Circulation Surge)", 1985, 2000),
        ("2000-2015 (Modern Superinfection Wave)", 2000, 2016),
    ]
    for label, y0, y1 in epochs:
        sub = crf_dated[(crf_dated["t_mrca"] >= y0) & (crf_dated["t_mrca"] < y1)]
        rate_per_yr = len(sub) / (y1 - y0)
        print(f"  {label:40s}: {len(sub):3d} established CRFs ({rate_per_yr:.2f} new CRFs / year)")

    # =========================================================================
    # PART 2: THE EPISTATIC & STRUCTURAL SIEVE (ESTABLISHED vs. DEAD-END URFs)
    # =========================================================================
    print("\n" + "=" * 88)
    print("PART 2: EPISTATIC & STRUCTURAL SIEVE (ESTABLISHED CRFs vs. DEAD-END SINGLETON URFs)")
    print("=" * 88)

    # Define the 3 empirical tiers of recombinant architectures:
    # Tier 1: Sustained Epidemic CRFs (N >= 10 patients in GenBank, or top pandemic CRFs)
    # Tier 2: All Established Circulating Recombinants (145 verified mosaic CRFs in LANL compendium + 15 Novel Circulating Lineages in GenBank)
    # Tier 3: Strictly Singleton Dead-End URFs (N = 1 patient in GenBank)

    tier_records = []
    gb_crf_counts = Counter()
    for _, r in recomb_gb[known_mask].iterrows():
        norm_lin = str(r["matched_known_lineage"]).split("/")[0].strip()
        gb_crf_counts[norm_lin] += 1

    gen_map = dict(zip(dag_df["legacy_label"], dag_df["generation"]))

    seen_comp_labels = set()
    for _, r in crf_comp.iterrows():
        bc = r["hxb2_barcode"]
        segs = parse_barcode_map(bc)
        if len(segs) <= 1:
            continue
        code = str(r["short_code"]).strip()
        leg = str(r["legacy_label"]).strip()
        seen_comp_labels.add(leg)
        n_gb = max(gb_crf_counts.get(code, 0), gb_crf_counts.get(leg, 0), int(r["n_ref_genomes"]))
        gen_val = int(gen_map.get(leg, 1))
        tier_records.append({
            "lineage_id": leg,
            "tier": "Established_Circulating",
            "is_major_epidemic": n_gb >= 10,
            "is_primary_dag": gen_val == 1,
            "generation": gen_val,
            "n_patients": n_gb,
            "segs": segs,
            "n_segments": len(segs),
            "n_breakpoints": len(segs) - 1,
        })

    # Add the 8 post-2023 official mosaic CRFs from compendium_bundle.json (not in 2023 crf_comp)
    bundle_entries = json.loads(Path("rhiztype/crates/eggman-core/data/compendium_bundle.json").read_text())["compendium"]
    for ent in bundle_entries:
        leg = str(ent.get("label", "")).strip()
        cat_b = str(ent.get("category", "")).strip()
        if not cat_b.startswith("CRF") or leg in seen_comp_labels:
            continue
        segs = parse_barcode_map(ent.get("hxb2_barcode", ""))
        if len(segs) <= 1:
            continue
        seen_comp_labels.add(leg)
        code = str(ent.get("short_code", "")).strip()
        n_gb = max(gb_crf_counts.get(code, 0), gb_crf_counts.get(leg, 0), 1)
        gen_val = int(gen_map.get(leg, 1))
        tier_records.append({
            "lineage_id": leg,
            "tier": "Established_Circulating",
            "is_major_epidemic": n_gb >= 10,
            "is_primary_dag": gen_val == 1,
            "generation": gen_val,
            "n_patients": n_gb,
            "segs": segs,
            "n_segments": len(segs),
            "n_breakpoints": len(segs) - 1,
        })

    # Add strict ML-monophyletic novel regional circulating CRF(s) (RU_G0B_gagMA)
    verified_nov = novel_df[novel_df["lineage_status"].astype(str).str.startswith("Candidate_Regional_Mosaic")].copy()
    for _, r in verified_nov.iterrows():
        form = r["rhizaeon_formula"]
        bc = str(r["representative_hxb2_barcode"])
        segs = parse_barcode_map(bc)
        if len(segs) <= 1:
            continue
        tier_records.append({
            "lineage_id": f"Novel_{form}",
            "tier": "Established_Circulating",
            "is_major_epidemic": int(r["n_homologous_bp_patients"]) >= 10,
            "is_primary_dag": True,
            "generation": 1,
            "n_patients": int(r["n_homologous_bp_patients"]),
            "segs": segs,
            "n_segments": len(segs),
            "n_breakpoints": len(segs) - 1,
        })

    # 2. Add Strictly Singleton Dead-End URFs (N = 1 patient in GenBank; 461 of 508 N=1 rows after excluding 47 Known CRF singletons)
    for key, meta in lineage_meta.items():
        if meta["type"] == "Singleton_URF":
            segs = parse_barcode_map(meta["barcode"])
            if len(segs) <= 1:
                continue
            tier_records.append({
                "lineage_id": key,
                "tier": "Singleton_URF",
                "is_major_epidemic": False,
                "is_primary_dag": False,
                "generation": 1,
                "n_patients": 1,
                "segs": segs,
                "n_segments": len(segs),
                "n_breakpoints": len(segs) - 1,
            })

    est_list = [x for x in tier_records if x["tier"] == "Established_Circulating"]
    pri_list = [x for x in est_list if x["is_primary_dag"]]
    maj_list = [x for x in tier_records if x["is_major_epidemic"]]
    urf_list = [x for x in tier_records if x["tier"] == "Singleton_URF"]

    n_known_n1 = sum(1 for k, v in lineage_meta.items() if v["type"] == "Known_CRF" and lineage_counts.get(k, 0) == 1)
    print(f"Comparative Architecture Cohorts (with explicit reconciliation from {n_total_emitted_lineages} GenBank survey clusters):")
    print(f"  - Cohort A (Major Epidemic CRFs, canonical N >= 10):    {len(maj_list):4d} distinct lineages (vs {n_ge10} raw survey clusters)")
    print(f"  - Cohort B (All Established Circulating):               {len(est_list):4d} distinct lineages ({len(crf_comp)} 2023 LANL + {len(seen_comp_labels)-len(crf_comp)} post-2023 CRFs + {len(verified_nov)} Candidate Regional)")
    print(f"  - Cohort B0 (DAG-Deduplicated Primary Gen-1 CRFs):      {len(pri_list):4d} independent primary origins")
    print(f"  - Cohort C (Dead-End Singleton URFs, N = 1 patient):    {len(urf_list):4d} architectures ({n_n1} total N=1 minus {n_known_n1} Known CRFs appearing once in GenBank)")

    # Helper functions for circular-shift, 10k-rep linear permutation spatial null, CMH & Logistic regression
    H0, H1 = 790, 9417
    CODING_LEN = H1 - H0 + 1  # 8,628 nt

    def _cohort_to_grid(cohort):
        arr = np.zeros((len(cohort), CODING_LEN), dtype=np.int16)
        for i, item in enumerate(cohort):
            tmap = {}
            for s, e, trk in item["segs"]:
                if trk not in tmap:
                    tmap[trk] = len(tmap) + 1
                s_c, e_c = max(H0, s) - H0, min(H1, e) - H0
                if e_c >= s_c:
                    arr[i, s_c:e_c + 1] = tmap[trk]
        return arr

    def _circular_shift_dist(grid, pos1, pos2):
        shifts = np.arange(500, CODING_LEN - 500 + 1)
        i1, i2 = pos1 - H0, pos2 - H0
        s1 = (i1 + shifts) % CODING_LEN
        s2 = (i2 + shifts) % CODING_LEN
        dist = np.mean(grid[:, s1] == grid[:, s2], axis=0)
        return dist

    def _linear_perm_null(cohort, lock_coords, n_rep=10000, seed=20261114):
        rng = np.random.default_rng(seed)
        match_rates = np.zeros((len(lock_coords), n_rep), dtype=np.float64)
        lpos = [(p1 - H0, p2 - H0) for p1, p2 in lock_coords]
        for item in cohort:
            segs = item["segs"]
            lens = np.array([max(1, e - s + 1) for s, e, _ in segs], dtype=np.float64)
            lens = lens * (float(CODING_LEN) / lens.sum())
            trks_raw = [t for _, _, t in segs]
            uniq_t = sorted(set(trks_raw))
            t_idx = {t: idx for idx, t in enumerate(uniq_t)}
            K = len(segs)
            M = len(uniq_t)
            t_lens = np.zeros(M, dtype=np.float64)
            for l_k, t_k in zip(lens, trks_raw):
                t_lens[t_idx[t_k]] += l_k
            t_probs = t_lens / t_lens.sum()

            keys = rng.random((n_rep, K))
            perm_lens = lens[np.argsort(keys, axis=1)]
            cum_ends = np.cumsum(perm_lens, axis=1)

            if M < 2:
                trks_arr = np.zeros((n_rep, K), dtype=np.int16)
            elif M == 2:
                trks_arr = np.broadcast_to((np.arange(K) % 2)[None, :], (n_rep, K))
            else:
                trks_arr = np.empty((n_rep, K), dtype=np.int16)
                u0 = rng.random(n_rep)
                trks_arr[:, 0] = np.searchsorted(np.cumsum(t_probs), u0, side="right").clip(0, M - 1)
                for k in range(1, K):
                    prev = trks_arr[:, k - 1]
                    w = np.broadcast_to(t_probs[None, :], (n_rep, M)).copy()
                    w[np.arange(n_rep), prev] = 0.0
                    w_sum = w.sum(axis=1, keepdims=True)
                    w = np.where(w_sum > 0, w / w_sum, 1.0 / (M - 1))
                    w[np.arange(n_rep), prev] = 0.0
                    cdf = np.cumsum(w, axis=1)
                    uk = rng.random((n_rep, 1))
                    trks_arr[:, k] = (uk > cdf).sum(axis=1).clip(0, M - 1)

            row_idx = np.arange(n_rep)
            for li, (x1, x2) in enumerate(lpos):
                seg1 = np.sum(cum_ends < x1, axis=1).clip(0, K - 1)
                seg2 = np.sum(cum_ends < x2, axis=1).clip(0, K - 1)
                match_rates[li] += (trks_arr[row_idx, seg1] == trks_arr[row_idx, seg2])

        match_rates /= max(1, len(cohort))
        return match_rates

    def _bp_stratum(nbp):
        if nbp <= 2:
            return "1-2"
        if nbp <= 4:
            return "3-4"
        if nbp <= 6:
            return "5-6"
        return "7+"

    def _stratified_tests(pos_cohort, neg_cohort, feat_fn):
        """Run Cochran-Mantel-Haenszel (CMH) and n_breakpoints-adjusted Logistic Regression."""
        rows_s = []
        for it in pos_cohort:
            rows_s.append({"y": 1, "nbp": it["n_breakpoints"], "stratum": _bp_stratum(it["n_breakpoints"]), "feat": int(bool(feat_fn(it)))})
        for it in neg_cohort:
            rows_s.append({"y": 0, "nbp": it["n_breakpoints"], "stratum": _bp_stratum(it["n_breakpoints"]), "feat": int(bool(feat_fn(it)))})
        sdf = pd.DataFrame(rows_s)
        tables = []
        for st in ["1-2", "3-4", "5-6", "7+"]:
            sub = sdf[sdf["stratum"] == st]
            a = int(((sub["y"] == 1) & (sub["feat"] == 1)).sum())
            b = int(((sub["y"] == 1) & (sub["feat"] == 0)).sum())
            c = int(((sub["y"] == 0) & (sub["feat"] == 1)).sum())
            d = int(((sub["y"] == 0) & (sub["feat"] == 0)).sum())
            if (a + b) > 0 and (c + d) > 0 and (a + c) > 0 and (b + d) > 0:
                tables.append([[a, b], [c, d]])
        cmh_or, cmh_p = np.nan, np.nan
        if tables:
            import warnings
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", RuntimeWarning)
                stbl = StratifiedTable(tables)
                cmh_or = float(stbl.oddsratio_pooled)
                if not np.isfinite(cmh_or):
                    tables_h = [[[t[0][0] + 0.5, t[0][1] + 0.5], [t[1][0] + 0.5, t[1][1] + 0.5]] for t in tables]
                    cmh_or = float(StratifiedTable(tables_h).oddsratio_pooled)
                cmh_p = float(stbl.test_null_odds().pvalue)
        logit_or, logit_p = np.nan, np.nan
        try:
            X = sm.add_constant(sdf[["nbp", "feat"]])
            mod = sm.Logit(sdf["y"], X).fit(disp=False, warn_convergence=False)
            logit_or = float(np.exp(mod.params["feat"]))
            logit_p = float(mod.pvalues["feat"])
        except Exception:
            pass
        return cmh_or, cmh_p, logit_or, logit_p

    # -------------------------------------------------------------------------
    # FEATURE CLASS I: COORDINATED INTRA-COMPLEX CO-ANCESTRY LOCKS
    # -------------------------------------------------------------------------
    print("\n" + "-" * 88)
    print("FEATURE CLASS I: COORDINATED INTRA-COMPLEX CO-ANCESTRY LOCKS & SPATIAL NULLS")
    print("-" * 88)

    lock_pairs = [
        (6900, 7950, "gp120 Core <-> gp41 Ectodomain", "Env Trimer Spike Interface (6900 vs 7950)"),
        (950,  8650, "Gag-MA Matrix <-> gp41-CT Tail", "Virion Env Incorporation Lock (950 vs 8650)"),
        (1050, 8600, "Gag-MA <-> gp41-CT Alt Window", "Virion Env Incorporation Lock Sensitivity (1050 vs 8600)"),
        (7880, 8500, "RRE Stem-Loop <-> Rev Exon 2",   "Rev-RRE Nuclear RNA Export Lock (7880 vs 8500)"),
        (5950, 8400, "Tat Exon 1 <-> Tat Exon 2",      "Spliced Trans-Activator Lock (5950 vs 8400)"),
        (1550, 2400, "Gag-CA/NC <-> Pol-PR Protease",  "Gag Substrate <-> Protease Cleavage Lock (1550 vs 2400)"),
        (3000, 4550, "Pol-RT Core <-> Pol-IN Core",    "Replicase Catalytic Lock (3000 vs 4550)"),
    ]
    lock_coords = [(p1, p2) for p1, p2, _, _ in lock_pairs]

    grid_maj = _cohort_to_grid(maj_list)
    grid_est = _cohort_to_grid(est_list)
    grid_pri = _cohort_to_grid(pri_list)
    grid_urf = _cohort_to_grid(urf_list)

    lin_maj = _linear_perm_null(maj_list, lock_coords, n_rep=10000, seed=101)
    lin_est = _linear_perm_null(est_list, lock_coords, n_rep=10000, seed=102)
    lin_pri = _linear_perm_null(pri_list, lock_coords, n_rep=10000, seed=103)
    lin_urf = _linear_perm_null(urf_list, lock_coords, n_rep=10000, seed=104)

    sieve_rows = []
    detailed_rows = []

    for li, (pos1, pos2, short_name, desc) in enumerate(lock_pairs):
        def eval_lock(cohort):
            matches = 0
            per_item_p = []
            for item in cohort:
                segs = item["segs"]
                if trunk_at_pos(segs, pos1) == trunk_at_pos(segs, pos2):
                    matches += 1
                props = get_trunk_proportions(segs)
                per_item_p.append(sum(p * p for p in props.values()))
            n = len(cohort)
            obs_rate = matches / n if n > 0 else 0.0
            exp_rate = sum(per_item_p) / n if n > 0 else 0.0
            return matches, n, obs_rate, exp_rate, obs_rate - exp_rate, per_item_p

        def poisson_binomial_stats(per_item_p, obs_k, n_obs, obs_rate, exp_rate):
            if n_obs == 0:
                return 0.0, 1.0
            pb_var = sum(p * (1.0 - p) for p in per_item_p) / (n_obs * n_obs)
            pb_sd = math.sqrt(max(pb_var, 1e-18))
            z_pb = (obs_rate - exp_rate) / pb_sd
            dp = np.array([1.0], dtype=np.float64)
            for p_i in per_item_p:
                ndp = np.zeros(len(dp) + 1, dtype=np.float64)
                ndp[:-1] += dp * (1.0 - p_i)
                ndp[1:] += dp * p_i
                dp = ndp
            p_le = float(dp[:obs_k + 1].sum())
            p_ge = float(dp[obs_k:].sum())
            p_2s = min(1.0, 2.0 * min(p_le, p_ge))
            return z_pb, p_2s

        m_maj, n_maj, obs_maj, exp_maj, exc_maj, _ = eval_lock(maj_list)
        m_est, n_est, obs_est, exp_est, exc_est, pb_p_est = eval_lock(est_list)
        m_pri, n_pri, obs_pri, exp_pri, exc_pri, _ = eval_lock(pri_list)
        m_urf, n_urf, obs_urf, exp_urf, exc_urf, _ = eval_lock(urf_list)

        or_maj_urf, p_maj_vs_urf = stats.fisher_exact([[m_maj, n_maj - m_maj], [m_urf, n_urf - m_urf]])
        or_est_urf, p_est_vs_urf = stats.fisher_exact([[m_est, n_est - m_est], [m_urf, n_urf - m_urf]])
        or_pri_urf, p_pri_vs_urf = stats.fisher_exact([[m_pri, n_pri - m_pri], [m_urf, n_urf - m_urf]])

        z_vs_null_est, p_vs_null_est = poisson_binomial_stats(pb_p_est, m_est, n_est, obs_est, exp_est)

        # Circular-shift distribution (s in [500, L-500])
        circ_est = _circular_shift_dist(grid_est, pos1, pos2)
        circ_pri = _circular_shift_dist(grid_pri, pos1, pos2)
        circ_maj = _circular_shift_dist(grid_maj, pos1, pos2)
        circ_urf = _circular_shift_dist(grid_urf, pos1, pos2)

        # Two-sided empirical tail probability helper vs spatial distribution
        def two_sided_emp(dist, obs):
            p_ge = float(np.mean(dist >= obs))
            p_le = float(np.mean(dist <= obs))
            return float(np.mean(dist)), p_ge, p_le, min(1.0, 2.0 * min(p_ge, p_le))

        cm_est, p_circ_ge_est, p_circ_le_est, p_circ_2_est = two_sided_emp(circ_est, obs_est)
        cm_pri, p_circ_ge_pri, p_circ_le_pri, p_circ_2_pri = two_sided_emp(circ_pri, obs_pri)
        lm_est, p_lin_ge_est, p_lin_le_est, p_lin_2_est = two_sided_emp(lin_est[li], obs_est)
        lm_pri, p_lin_ge_pri, p_lin_le_pri, p_lin_2_pri = two_sided_emp(lin_pri[li], obs_pri)
        lm_maj, p_lin_ge_maj, p_lin_le_maj, p_lin_2_maj = two_sided_emp(lin_maj[li], obs_maj)
        lm_urf, p_lin_ge_urf, p_lin_le_urf, p_lin_2_urf = two_sided_emp(lin_urf[li], obs_urf)

        # Stratified CMH and logistic regression controlling for n_breakpoints
        feat_lock = lambda it, p1=pos1, p2=pos2: trunk_at_pos(it["segs"], p1) == trunk_at_pos(it["segs"], p2)
        cmh_or_est, cmh_p_est, log_or_est, log_p_est = _stratified_tests(est_list, urf_list, feat_lock)
        cmh_or_pri, cmh_p_pri, log_or_pri, log_p_pri = _stratified_tests(pri_list, urf_list, feat_lock)
        cmh_or_maj, cmh_p_maj, log_or_maj, log_p_maj = _stratified_tests(maj_list, urf_list, feat_lock)

        print(f"\n* {short_name} [{desc}]")
        print(f"  - Major CRFs (N={n_maj:3d}):         Obs = {100*obs_maj:5.1f}% ({m_maj:3d}/{n_maj:3d}) | LinNull={100*lm_maj:4.1f}% (p2={p_lin_2_maj:.4f}) | Circ={100*circ_maj.mean():4.1f}%")
        print(f"  - All Established (N={n_est:3d}):  Obs = {100*obs_est:5.1f}% ({m_est:3d}/{n_est:3d}) | LinNull={100*lm_est:4.1f}% (p_ge={p_lin_ge_est:.4f}, p_le={p_lin_le_est:.4f}) | Circ={100*cm_est:4.1f}% (p_ge={p_circ_ge_est:.4f})")
        print(f"  - DAG Primary Gen1 (N={n_pri:3d}): Obs = {100*obs_pri:5.1f}% ({m_pri:3d}/{n_pri:3d}) | LinNull={100*lm_pri:4.1f}% (p_ge={p_lin_ge_pri:.4f}, p_le={p_lin_le_pri:.4f}) | Circ={100*cm_pri:4.1f}% (p_ge={p_circ_ge_pri:.4f})")
        print(f"  - Dead-End URFs (N={n_urf:3d}):    Obs = {100*obs_urf:5.1f}% ({m_urf:3d}/{n_urf:3d}) | LinNull={100*lm_urf:4.1f}% (p_ge={p_lin_ge_urf:.4f}, p_le={p_lin_le_urf:.4f})")
        print(f"  -> Fisher vs URF: Est OR={or_est_urf:.2f} (p={p_est_vs_urf:.4f}) | Pri OR={or_pri_urf:.2f} (p={p_pri_vs_urf:.4f}) | Maj OR={or_maj_urf:.2f} (p={p_maj_vs_urf:.4f})")
        print(f"  -> Breakpoint-Stratified CMH: Est OR={cmh_or_est:.2f} (p={cmh_p_est:.4e}) | Pri OR={cmh_or_pri:.2f} (p={cmh_p_pri:.4e}) | Maj OR={cmh_or_maj:.2f} (p={cmh_p_maj:.4e})")

        sieve_rows.append({
            "feature_category": "Intra_Complex_CoAncestry_Lock",
            "feature_name": short_name,
            "description": desc,
            "major_crf_rate": round(obs_maj, 4),
            "major_crf_null": round(exp_maj, 4),
            "major_crf_excess": round(exc_maj, 4),
            "established_rate": round(obs_est, 4),
            "established_null": round(exp_est, 4),
            "established_excess": round(exc_est, 4),
            "primary_dag_rate": round(obs_pri, 4),
            "primary_dag_lin_null": round(lm_pri, 4),
            "singleton_urf_rate": round(obs_urf, 4),
            "singleton_urf_null": round(exp_urf, 4),
            "singleton_urf_excess": round(exc_urf, 4),
            "est_circ_null": round(cm_est, 4),
            "est_lin_perm_null": round(lm_est, 4),
            "p_est_vs_circ_ge": p_circ_ge_est,
            "p_est_vs_lin_ge": p_lin_ge_est,
            "p_est_vs_lin_le": p_lin_le_est,
            "p_pri_vs_lin_ge": p_lin_ge_pri,
            "z_established_vs_random_null": round(z_vs_null_est, 3),
            "p_established_vs_random_null": p_vs_null_est,
            "p_major_vs_urf": p_maj_vs_urf,
            "p_established_vs_urf": p_est_vs_urf,
            "p_primary_vs_urf": p_pri_vs_urf,
            "cmh_or_est_vs_urf": round(cmh_or_est, 3) if not np.isnan(cmh_or_est) else np.nan,
            "cmh_p_est_vs_urf": cmh_p_est,
            "logit_or_est_vs_urf": round(log_or_est, 3) if not np.isnan(log_or_est) else np.nan,
            "logit_p_est_vs_urf": log_p_est,
            "cmh_or_pri_vs_urf": round(cmh_or_pri, 3) if not np.isnan(cmh_or_pri) else np.nan,
            "cmh_p_pri_vs_urf": cmh_p_pri,
            "cmh_or_maj_vs_urf": round(cmh_or_maj, 3) if not np.isnan(cmh_or_maj) else np.nan,
            "cmh_p_maj_vs_urf": cmh_p_maj,
        })

        for ctag, c_n, c_m, c_obs, c_circ, c_lin, p_c, p_lg, p_ll in [
            ("Major_CRFs", n_maj, m_maj, obs_maj, float(circ_maj.mean()), lm_maj, float(np.mean(circ_maj >= obs_maj)), p_lin_ge_maj, p_lin_le_maj),
            ("All_Established", n_est, m_est, obs_est, cm_est, lm_est, p_circ_ge_est, p_lin_ge_est, p_lin_le_est),
            ("Primary_DAG_Gen1", n_pri, m_pri, obs_pri, cm_pri, lm_pri, p_circ_ge_pri, p_lin_ge_pri, p_lin_le_pri),
            ("Singleton_URF", n_urf, m_urf, obs_urf, float(circ_urf.mean()), lm_urf, float(np.mean(circ_urf >= obs_urf)), p_lin_ge_urf, p_lin_le_urf),
        ]:
            detailed_rows.append({
                "cohort": ctag, "n_lineages": c_n, "feature_type": "CoAncestry_Lock",
                "feature_name": short_name, "coords": f"{pos1}<->{pos2}",
                "observed_k": c_m, "observed_rate": round(c_obs, 4),
                "circ_shift_null": round(c_circ, 4), "p_circ_ge": p_c,
                "linear_perm_null": round(c_lin, 4), "p_lin_ge": p_lg, "p_lin_le": p_ll,
            })

    # -------------------------------------------------------------------------
    # SPECIAL TEST: THE COORDINATED "ECTODOMAIN-ONLY ENV SWAP" (Gag-MA <-> gp41-CT Return)
    # -------------------------------------------------------------------------
    print("\n" + "-" * 88)
    print("SPECIAL COORDINATED CASSETTE TEST: ECTODOMAIN ENV SWAP WITH gp41-CT REVERSION TO Gag-MA")
    print("-" * 88)
    def eval_env_cassette_reversion(cohort):
        n_foreign_ecto = 0
        n_ct_reverts_to_ma = 0
        n_nef_reverts_to_ma = 0
        for item in cohort:
            segs = item["segs"]
            t_ma = trunk_at_pos(segs, 950)
            t_ecto = trunk_at_pos(segs, 6900)
            t_ct = trunk_at_pos(segs, 8650)
            t_nef = trunk_at_pos(segs, 8950)
            if t_ecto != t_ma:
                n_foreign_ecto += 1
                if t_ct == t_ma:
                    n_ct_reverts_to_ma += 1
                if t_nef == t_ma:
                    n_nef_reverts_to_ma += 1
        return n_foreign_ecto, n_ct_reverts_to_ma, n_nef_reverts_to_ma

    fe_maj, ct_maj, nef_maj = eval_env_cassette_reversion(maj_list)
    fe_est, ct_est, nef_est = eval_env_cassette_reversion(est_list)
    fe_pri, ct_pri, nef_pri = eval_env_cassette_reversion(pri_list)
    fe_urf, ct_urf, nef_urf = eval_env_cassette_reversion(urf_list)

    or_ct, p_ct = stats.fisher_exact([[ct_maj, fe_maj - ct_maj], [ct_urf, fe_urf - ct_urf]])
    or_ct_est, p_ct_est = stats.fisher_exact([[ct_est, fe_est - ct_est], [ct_urf, fe_urf - ct_urf]])
    or_ct_pri, p_ct_pri = stats.fisher_exact([[ct_pri, fe_pri - ct_pri], [ct_urf, fe_urf - ct_urf]])
    or_nef, p_nef = stats.fisher_exact([[nef_maj, fe_maj - nef_maj], [nef_urf, fe_urf - nef_urf]])

    print(f"When Env Ectodomain (gp120) is Swapped Away from Gag-MA [T(gp120) != T(Gag-MA)]:")
    print(f"  1) Coordinated Reversion Inside gp41-CT (8650) Back to Match Gag-MA:")
    print(f"     - Major CRFs (N>=10):        {ct_maj:3d} / {fe_maj:3d} ({100*ct_maj/max(1,fe_maj):5.1f}%)")
    print(f"     - All Established:           {ct_est:3d} / {fe_est:3d} ({100*ct_est/max(1,fe_est):5.1f}%)")
    print(f"     - Primary DAG Gen1:          {ct_pri:3d} / {fe_pri:3d} ({100*ct_pri/max(1,fe_pri):5.1f}%)")
    print(f"     - Dead-End Singleton URFs:   {ct_urf:3d} / {fe_urf:3d} ({100*ct_urf/max(1,fe_urf):5.1f}%)")
    print(f"     -> Major vs URF OR={or_ct:.2f} (p={p_ct:.3e}); Est vs URF OR={or_ct_est:.2f} (p={p_ct_est:.3e}); Pri vs URF OR={or_ct_pri:.2f} (p={p_ct_pri:.3e})")
    print(f"  2) Coordinated Reversion by gp41/nef (8950) Back to Match Gag-MA:")
    print(f"     - Major CRFs / Est / Pri / URF: {nef_maj}/{fe_maj} | {nef_est}/{fe_est} | {nef_pri}/{fe_pri} | {nef_urf}/{fe_urf} (Maj vs URF OR={or_nef:.2f}, p={p_nef:.3e})")

    sieve_rows.append({
        "feature_category": "Coordinated_Cassette_Swap",
        "feature_name": "Foreign_gp120_Reverts_at_gp41_CT_to_Match_Gag_MA",
        "description": "Given T(gp120)!=T(Gag-MA), second breakpoint restores T(gp41-CT)==T(Gag-MA)",
        "major_crf_rate": round(ct_maj / max(1, fe_maj), 4),
        "major_crf_null": np.nan, "major_crf_excess": np.nan,
        "established_rate": round(ct_est / max(1, fe_est), 4),
        "established_null": np.nan, "established_excess": np.nan,
        "primary_dag_rate": round(ct_pri / max(1, fe_pri), 4),
        "primary_dag_lin_null": np.nan,
        "singleton_urf_rate": round(ct_urf / max(1, fe_urf), 4),
        "singleton_urf_null": np.nan, "singleton_urf_excess": np.nan,
        "est_circ_null": np.nan, "est_lin_perm_null": np.nan,
        "p_est_vs_circ_ge": np.nan, "p_est_vs_lin_ge": np.nan, "p_est_vs_lin_le": np.nan, "p_pri_vs_lin_ge": np.nan,
        "z_established_vs_random_null": np.nan, "p_established_vs_random_null": np.nan,
        "p_major_vs_urf": p_ct, "p_established_vs_urf": p_ct_est, "p_primary_vs_urf": p_ct_pri,
        "cmh_or_est_vs_urf": np.nan, "cmh_p_est_vs_urf": np.nan, "logit_or_est_vs_urf": np.nan, "logit_p_est_vs_urf": np.nan,
        "cmh_or_pri_vs_urf": np.nan, "cmh_p_pri_vs_urf": np.nan, "cmh_or_maj_vs_urf": np.nan, "cmh_p_maj_vs_urf": np.nan,
    })

    # -------------------------------------------------------------------------
    # FEATURE CLASS II: BREAKPOINT PLACEMENT SIEVE (MODULAR HINGES vs. STRUCTURAL CORES)
    # -------------------------------------------------------------------------
    print("\n" + "-" * 88)
    print("FEATURE CLASS II: BREAKPOINT PLACEMENT SIEVE (PERMISSIVE HINGES vs. MID-DOMAIN CORES)")
    print("-" * 88)

    def get_all_bps(cohort):
        bps = []
        for item in cohort:
            for bp, t1, t2 in get_internal_breakpoints(item["segs"]):
                bps.append(bp)
        return np.array(bps, dtype=int)

    # DAG-deduplicated unique de novo junctions (counting each inherited founder junction once across Group M)
    dag_denovo_bps = []
    for _, dr in dag_df.iterrows():
        if int(dr["generation"]) == 0:
            continue
        dj = str(dr["denovo_junctions"])
        if dj in ("None", "nan", ""):
            continue
        for tok in dj.split(";"):
            tok = tok.strip()
            if "@" in tok:
                dag_denovo_bps.append(int(tok.split("@")[1].split("(")[0]))
    for item in est_list:
        if item["lineage_id"].startswith("Novel_"):
            for bp, _, _ in get_internal_breakpoints(item["segs"]):
                dag_denovo_bps.append(bp)
    bps_denovo = np.array(dag_denovo_bps, dtype=int)

    bps_maj = get_all_bps(maj_list)
    bps_est = get_all_bps(est_list)
    bps_pri = get_all_bps(pri_list)
    bps_urf = get_all_bps(urf_list)

    zones = [
        ("Permissive_Hinge: Accessory Hub (pol-IN..vpu, 4950-6350)", [(4950, 6350)]),
        ("Permissive_Hinge: gp41 TM/CT & nef Boundary (8150-8950)", [(8150, 8950)]),
        ("Permissive_Hinge: gag/pol Frameshift & PR/RT (1980-2620)", [(1980, 2620)]),
        ("Permissive_Hinge: gag.MA/CA Linker (1100-1260)", [(1100, 1260)]),
        ("ALL_PERMISSIVE_DOMAIN_HINGES_COMBINED", [(1100, 1260), (1980, 2620), (4950, 6350), (8150, 8950)]),
        ("Deleterious_Core: gp120 Ectodomain Beta-Sandwich (6450-7650)", [(6450, 7650)]),
        ("Deleterious_Core: Gag-CA Capsid Helical Core (1280-1850)", [(1280, 1850)]),
        ("Deleterious_Core: Pol-RT Catalytic Palm/Fingers (2650-3450)", [(2650, 3450)]),
        ("ALL_RIGID_STRUCTURAL_CORES_COMBINED", [(1280, 1850), (2650, 3450), (6450, 7650)]),
    ]

    def _cluster_zone_perm_p(cohort_pos, cohort_neg, intervals, n_rep=10000, seed=20261115):
        """Lineage-level cluster permutation test comparing zone share (accounting for within-lineage breakpoint correlation)."""
        rng_c = np.random.default_rng(seed)
        def _counts(c_list):
            tot_v = np.empty(len(c_list), dtype=int)
            hit_v = np.empty(len(c_list), dtype=int)
            for i, it in enumerate(c_list):
                bps_i = [bp for bp, _, _ in get_internal_breakpoints(it["segs"])]
                tot_v[i] = len(bps_i)
                hit_v[i] = sum(1 for bp in bps_i if any(s <= bp <= e for s, e in intervals))
            return hit_v, tot_v
        hp, tp = _counts(cohort_pos)
        hn, tn = _counts(cohort_neg)
        obs_diff = abs((hp.sum() / max(1, tp.sum())) - (hn.sum() / max(1, tn.sum())))
        all_h = np.concatenate([hp, hn])
        all_t = np.concatenate([tp, tn])
        n_pos = len(hp)
        n_tot = len(all_h)
        tot_sum = all_t.sum()
        hit_sum = all_h.sum()
        ge_cnt = 0
        for _ in range(n_rep):
            idx_p = rng_c.choice(n_tot, size=n_pos, replace=False)
            tp_s = int(all_t[idx_p].sum())
            hp_s = int(all_h[idx_p].sum())
            tn_s = int(tot_sum - tp_s)
            hn_s = int(hit_sum - hp_s)
            if tp_s == 0 or tn_s == 0:
                continue
            d_star = abs((hp_s / tp_s) - (hn_s / tn_s))
            if d_star >= obs_diff - 1e-12:
                ge_cnt += 1
        return float((ge_cnt + 1.0) / (n_rep + 1.0))

    for z_idx, (z_name, intervals) in enumerate(zones):
        z_len = sum(e - s + 1 for s, e in intervals)
        unif_frac = z_len / CODING_LEN

        def count_in_zone(bp_arr):
            mask = np.zeros(len(bp_arr), dtype=bool)
            for s, e in intervals:
                mask |= ((bp_arr >= s) & (bp_arr <= e))
            c = int(np.sum(mask))
            return c, len(bp_arr), (c / len(bp_arr) if len(bp_arr) > 0 else 0.0)

        c_maj, tot_maj, f_maj = count_in_zone(bps_maj)
        c_est, tot_est, f_est = count_in_zone(bps_est)
        c_pri, tot_pri, f_pri = count_in_zone(bps_pri)
        c_dnv, tot_dnv, f_dnv = count_in_zone(bps_denovo)
        c_urf, tot_urf, f_urf = count_in_zone(bps_urf)

        _, p_maj = stats.fisher_exact([[c_maj, tot_maj - c_maj], [c_urf, tot_urf - c_urf]])
        _, p_est = stats.fisher_exact([[c_est, tot_est - c_est], [c_urf, tot_urf - c_urf]])
        _, p_pri = stats.fisher_exact([[c_pri, tot_pri - c_pri], [c_urf, tot_urf - c_urf]])
        _, p_dnv = stats.fisher_exact([[c_dnv, tot_dnv - c_dnv], [c_urf, tot_urf - c_urf]])
        p_binom_est = stats.binomtest(c_est, tot_est, unif_frac).pvalue
        p_binom_dnv = stats.binomtest(c_dnv, tot_dnv, unif_frac).pvalue
        p_clust_est = _cluster_zone_perm_p(est_list, urf_list, intervals, n_rep=10000, seed=20260 + z_idx)
        p_clust_pri = _cluster_zone_perm_p(pri_list, urf_list, intervals, n_rep=10000, seed=30260 + z_idx)

        # Breakpoint-count stratified CMH & logistic regression at lineage level (any breakpoint in zone)
        has_zone_bp = lambda it, ivs=intervals: any(s <= bp <= e for bp, _, _ in get_internal_breakpoints(it["segs"]) for s, e in ivs)
        cmh_or_est, cmh_p_est, log_or_est, log_p_est = _stratified_tests(est_list, urf_list, has_zone_bp)
        cmh_or_pri, cmh_p_pri, _, _ = _stratified_tests(pri_list, urf_list, has_zone_bp)
        cmh_or_maj, cmh_p_maj, _, _ = _stratified_tests(maj_list, urf_list, has_zone_bp)

        enrichment_est = f_est / unif_frac if unif_frac > 0 else np.nan
        enrichment_dnv = f_dnv / unif_frac if unif_frac > 0 else np.nan
        enrichment_urf = f_urf / unif_frac if unif_frac > 0 else np.nan

        print(f"\n* {z_name} [Footprint = {100*unif_frac:4.1f}% ({z_len} nt)]")
        print(f"  - Major CRFs ({tot_maj:3d} bps):      {100*f_maj:5.1f}% ({c_maj:3d}/{tot_maj:3d}) | {f_maj/unif_frac:4.2f}x")
        print(f"  - All Estab ({tot_est:3d} bps):       {100*f_est:5.1f}% ({c_est:3d}/{tot_est:3d}) | {enrichment_est:4.2f}x (p_unif={p_binom_est:.2e}, cluster_perm_p={p_clust_est:.4e})")
        print(f"  - DAG DeNovo ({tot_dnv:3d} bps):      {100*f_dnv:5.1f}% ({c_dnv:3d}/{tot_dnv:3d}) | {enrichment_dnv:4.2f}x (p_unif={p_binom_dnv:.2e}, vs URF p={p_dnv:.3e})")
        print(f"  - DAG Primary ({tot_pri:3d} bps):     {100*f_pri:5.1f}% ({c_pri:3d}/{tot_pri:3d}) | {f_pri/unif_frac:4.2f}x (cluster_perm_p={p_clust_pri:.4e})")
        print(f"  - Dead-End URFs ({tot_urf:4d} bps):   {100*f_urf:5.1f}% ({c_urf:3d}/{tot_urf:4d}) | {enrichment_urf:4.2f}x")
        print(f"  -> Per-lineage n_bp-adj CMH: Est OR={cmh_or_est:.2f} (p={cmh_p_est:.3e}), Pri OR={cmh_or_pri:.2f} (p={cmh_p_pri:.3e}), Logit OR={log_or_est:.2f} (p={log_p_est:.3e})")

        sieve_rows.append({
            "feature_category": "Breakpoint_Structural_Placement",
            "feature_name": z_name,
            "description": f"Genome span {z_len} nt ({100*unif_frac:.1f}% uniform expectation)",
            "major_crf_rate": round(f_maj, 4),
            "major_crf_null": round(unif_frac, 4),
            "major_crf_excess": round(f_maj - unif_frac, 4),
            "established_rate": round(f_est, 4),
            "established_null": round(unif_frac, 4),
            "established_excess": round(f_est - unif_frac, 4),
            "primary_dag_rate": round(f_pri, 4),
            "primary_dag_lin_null": round(f_dnv, 4),  # store DAG de novo junction rate
            "singleton_urf_rate": round(f_urf, 4),
            "singleton_urf_null": round(unif_frac, 4),
            "singleton_urf_excess": round(f_urf - unif_frac, 4),
            "est_circ_null": round(unif_frac, 4),
            "est_lin_perm_null": round(enrichment_dnv, 3),
            "p_est_vs_circ_ge": p_binom_dnv,
            "p_est_vs_lin_ge": p_dnv,
            "p_est_vs_lin_le": p_clust_est,
            "p_pri_vs_lin_ge": p_clust_pri,
            "z_established_vs_random_null": round(enrichment_est, 3),
            "p_established_vs_random_null": p_binom_est,
            "p_major_vs_urf": p_maj,
            "p_established_vs_urf": p_est,
            "p_primary_vs_urf": p_pri,
            "cluster_perm_p_est_vs_urf": p_clust_est,
            "cluster_perm_p_pri_vs_urf": p_clust_pri,
            "cmh_or_est_vs_urf": round(cmh_or_est, 3) if not np.isnan(cmh_or_est) else np.nan,
            "cmh_p_est_vs_urf": cmh_p_est,
            "logit_or_est_vs_urf": round(log_or_est, 3) if not np.isnan(log_or_est) else np.nan,
            "logit_p_est_vs_urf": log_p_est,
            "cmh_or_pri_vs_urf": round(cmh_or_pri, 3) if not np.isnan(cmh_or_pri) else np.nan,
            "cmh_p_pri_vs_urf": cmh_p_pri,
            "cmh_or_maj_vs_urf": round(cmh_or_maj, 3) if not np.isnan(cmh_or_maj) else np.nan,
            "cmh_p_maj_vs_urf": cmh_p_maj,
        })

        for ctag, n_tot, c_k, f_r, p_vs_u in [
            ("Major_CRFs", tot_maj, c_maj, f_maj, p_maj),
            ("All_Established", tot_est, c_est, f_est, p_est),
            ("DAG_DeNovo_Unique", tot_dnv, c_dnv, f_dnv, p_dnv),
            ("Primary_DAG_Gen1", tot_pri, c_pri, f_pri, p_pri),
            ("Singleton_URF", tot_urf, c_urf, f_urf, 1.0),
        ]:
            detailed_rows.append({
                "cohort": ctag, "n_lineages": n_tot, "feature_type": "Breakpoint_Zone",
                "feature_name": z_name, "coords": ";".join(f"{s}-{e}" for s, e in intervals),
                "observed_k": c_k, "observed_rate": round(f_r, 4),
                "circ_shift_null": round(unif_frac, 4), "p_circ_ge": float(stats.binomtest(c_k, n_tot, unif_frac).pvalue) if n_tot > 0 else np.nan,
                "linear_perm_null": round(f_r / unif_frac, 3) if unif_frac > 0 else np.nan, "p_lin_ge": p_vs_u, "p_lin_le": np.nan,
            })

    # -------------------------------------------------------------------------
    # FEATURE CLASS III: GENOMIC COMPLEXITY & DONOR BURDEN
    # -------------------------------------------------------------------------
    print("\n" + "-" * 88)
    print("FEATURE CLASS III: BACKBONE INTEGRITY & SEGMENTATION COMPLEXITY")
    print("-" * 88)
    for name, cohort in [
        ("Major Epidemic CRFs (N>=10)", maj_list),
        ("All Established Circulating", est_list),
        ("Primary DAG Gen-1 CRFs", pri_list),
        ("Dead-End Singleton URFs (N=1)", urf_list),
    ]:
        n_bps = [x["n_breakpoints"] for x in cohort]
        dom_props = [max(get_trunk_proportions(x["segs"]).values()) for x in cohort]
        n_parents = [len(get_trunk_proportions(x["segs"])) for x in cohort]
        print(f"  - {name:36s}: Mean Breakpoints = {np.mean(n_bps):.2f} (median {np.median(n_bps):.0f}) | "
              f"Mean Dominant Backbone Share = {100*np.mean(dom_props):.1f}% | Mean Distinct Parents = {np.mean(n_parents):.2f}")

    sieve_df = pd.DataFrame(sieve_rows)
    # Apply Benjamini-Hochberg FDR across the 7 Feature Class I locks and across the 7 individual Feature Class II zones
    for col_p, col_q in [
        ("p_established_vs_urf", "q_bh_est_vs_urf"),
        ("cmh_p_est_vs_urf", "q_bh_cmh_est_vs_urf"),
        ("cluster_perm_p_est_vs_urf", "q_bh_cluster_perm_est_vs_urf"),
        ("p_est_vs_lin_ge", "q_bh_lin_perm_or_dnv"),
    ]:
        sieve_df[col_q] = np.nan
    m1 = sieve_df["feature_category"] == "Intra_Complex_CoAncestry_Lock"
    m2 = (sieve_df["feature_category"] == "Breakpoint_Structural_Placement") & (~sieve_df["feature_name"].astype(str).str.startswith("ALL_"))
    for mask in (m1, m2):
        for col_p, col_q in [
            ("p_established_vs_urf", "q_bh_est_vs_urf"),
            ("cmh_p_est_vs_urf", "q_bh_cmh_est_vs_urf"),
            ("cluster_perm_p_est_vs_urf", "q_bh_cluster_perm_est_vs_urf"),
            ("p_est_vs_lin_ge", "q_bh_lin_perm_or_dnv"),
        ]:
            idx_v = sieve_df.index[mask & sieve_df[col_p].notna()]
            if len(idx_v) > 0:
                _, qvals, _, _ = multipletests(sieve_df.loc[idx_v, col_p].astype(float).values, method="fdr_bh")
                sieve_df.loc[idx_v, col_q] = qvals

    sieve_df.to_csv("results/recombination_viability_and_establishment_sieve.csv", index=False)
    pd.DataFrame(detailed_rows).to_csv("results/sieve_epistatic_and_structural_summary.csv", index=False)

    # Export explicit 749-lineage survey manifest (reconciling 749 GenBank lineages vs 154 Established / 89 Primary DAG / 22 Major / 461 Singleton URF)
    est_ids = {x["lineage_id"] for x in est_list}
    pri_ids = {x["lineage_id"] for x in pri_list}
    maj_ids = {x["lineage_id"] for x in maj_list}
    urf_ids = {x["lineage_id"] for x in urf_list}
    manifest_rows = []
    for key, cnt in sorted(lineage_counts.items(), key=lambda kv: (-kv[1], kv[0])):
        meta = lineage_meta[key]
        manifest_rows.append({
            "lineage_key": key,
            "survey_type": meta["type"],
            "n_genbank_patients": int(cnt),
            "in_singleton_urf_461": key in urf_ids,
            "rep_accession": meta.get("rep_accession", ""),
            "n_segments": meta["n_segments"],
            "n_breakpoints": meta["n_breakpoints"],
            "formula": meta["formula"],
            "hxb2_barcode": meta["barcode"],
        })
    pd.DataFrame(manifest_rows).to_csv("results/supp_manifest_sieve_749_lineages.csv", index=False)
    print("\nSaved full sieve metrics to results/recombination_viability_and_establishment_sieve.csv")
    print("Saved detailed spatial/DAG summary to results/sieve_epistatic_and_structural_summary.csv")
    print("Saved 749-lineage manifest to results/supp_manifest_sieve_749_lineages.csv")
    print("Saved regional heterozygosity rates to results/regional_heterozygosity_urf_rates.csv")
    print("Saved multivariable H_circ regression to results/sieve_hcirc_multivariate_regression.csv")


if __name__ == "__main__":
    main()
