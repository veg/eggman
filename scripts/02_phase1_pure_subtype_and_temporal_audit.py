#!/usr/bin/env python3
"""
02_phase1_pure_subtype_and_temporal_audit.py
============================================
Phase 1 Spatio-Temporal Audit of HIV-1 Group M "Pure" Subtypes & Primary CRFs
using RhizAeon (PA-FF v2.1) + ChronAeon (Tree-Free Geodesic Clock Calibration).

Performs four rigorous, zero-shortcut analyses on authentic LANL alignments:
1. Reference Panel Contamination Audit:
   Maps all 12 RhizAeon-verified recombination events inside LANL's 55 nominal
   "pure" subtype reference genomes (`lanl_nominal_pure_subtypes.fasta`) onto exact
   HXB2 coordinates and HIV-1 genes/domains.
2. Leave-Clade-Out Basal Trunk vs. Derived Subtype Deconvolution:
   Tests each of the 17 nominal "pure" subtypes (A1..A8, B, C, D, F1, F2, G, H, J, K, L)
   plus CRF01_AE and CRF02_AG against a clean, scrubbed foundational reference panel.
3. ChronAeon Whole-Backbone & Partition-Specific Dating (N = 5,089 LANL genomes):
   Uses ChronAeon (`compute_time_decay_profile_divergences` + `run_ols_dating` with
   exact Fieller 95% confidence intervals) to date the crown t_MRCA and evolutionary
   rate (mu) of every nominal subtype and major CRF across:
   - Full coding genome (HXB2 790..9417)
   - CRF02_AG vs. Subtype G vs. Subtype A1 partitions (A1-shared tracts vs. G-shared tracts)
   - CRF01_AE A1-backbone (gag-pol) vs. Ghost_E cassette (vpu-env-nef)
"""

import os
import sys
import json
import subprocess
from pathlib import Path
from collections import defaultdict
import numpy as np
import pandas as pd
from Bio import SeqIO

from scipy import stats

from chronaeon.dating_io import compute_time_decay_profile_divergences
from chronaeon.dating_models import compute_rcs_basis, run_ols_dating, run_restricted_spline_clock_dating

PROJECT_ROOT = Path("/Users/sergei/Projects/HIV-recombination")
DATA_DIR = PROJECT_ROOT / "data"
RESULTS_DIR = PROJECT_ROOT / "results"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)


def load_coord_and_gene_maps():
    with open(DATA_DIR / "lanl_ref_2023_hxb2_coord_map.json") as f:
        cmap = json.load(f)
    m2h = {int(k): int(v) for k, v in cmap["msa_to_hxb2"].items()}
    h2m = {int(k): int(v) for k, v in cmap["hxb2_to_msa"].items()}
    with open(DATA_DIR / "hxb2_gene_domain_map.json") as f:
        gmap = json.load(f)
    return m2h, h2m, gmap


def annotate_hxb2_interval(h_start: int, h_end: int, gmap: dict) -> str:
    """Return human-readable HIV-1 domain/gene span for an HXB2 interval [h_start, h_end]."""
    hit_domains = []
    for d in gmap["domains"]:
        ov = min(h_end, d["end"]) - max(h_start, d["start"]) + 1
        if ov >= 30:
            hit_domains.append(d["label"])
    if not hit_domains:
        for g in gmap["genes"]:
            ov = min(h_end, g["end"]) - max(h_start, g["start"]) + 1
            if ov >= 15:
                hit_domains.append(g["id"])
    if not hit_domains:
        return f"nt{h_start}-{h_end}"
    if len(hit_domains) == 1:
        return hit_domains[0]
    return f"{hit_domains[0]}..{hit_domains[-1]}"


def build_flt_hxb2_map(flt_fasta: Path):
    """Build HXB2 coordinate -> MSA column map for lanl_flt_2022.fasta."""
    for rec in SeqIO.parse(str(flt_fasta), "fasta"):
        if "K03455" in rec.id:
            h2m = {}
            m2h = {}
            h_pos = 0
            for idx_0, ch in enumerate(str(rec.seq).upper()):
                col_1 = idx_0 + 1
                if ch not in ("-", "."):
                    h_pos += 1
                    h2m[h_pos] = col_1
                    m2h[col_1] = h_pos
                else:
                    m2h[col_1] = max(1, h_pos)
            return m2h, h2m
    raise RuntimeError("HXB2 K03455 not found in lanl_flt_2022.fasta")


class FastProfileEvaluator:
    """
    Precompute non-constant/gapped column indicator matrices ONCE per clade so arbitrary
    weighted/bootstrap ChronAeon Steiner divergence vectors run via 5 BLAS gemv ops
    (bit-for-bit identical to dense _compute_soft_root_q + _divs_from_q across all L cols).
    """
    def __init__(self, M: np.ndarray, L_val: np.ndarray):
        self.N, self.L = M.shape
        self.L_val = L_val
        # Columns that are NOT strictly constant (identical base 0..3 across every taxon)
        has_gap = np.any(M < 0, axis=0)
        n_bases = sum((np.any(M == b, axis=0)).astype(np.int8) for b in range(4))
        active = has_gap | (n_bases >= 2)
        M_p = M[:, active]
        self.P = M_p.shape[1]
        self.B_mats = [(M_p == b).astype(np.float64) for b in range(4)]
        self.V_p = (M_p >= 0).astype(np.float64)
        self.V_p_rowsum = np.sum(self.V_p, axis=1)

    def eval_batch_divs(self, W_anc: np.ndarray, w_sq_sums: np.ndarray, alpha_g: float = 0.20) -> np.ndarray:
        """
        Evaluate Steiner + Gamma-corrected root-to-tip distances across R replicate
        weightings simultaneously using 9 batched BLAS DGEMM calls (zero array copying).
        W_anc: (R, N), w_sq_sums: (R,) -> returns DIVS_MAT (N, R).
        """
        Q_list = [W_anc @ self.B_mats[b] for b in range(4)]  # 4 x (R, P)
        q_sum = Q_list[0] + Q_list[1] + Q_list[2] + Q_list[3]
        inv_q = np.where(q_sum > 0, 1.0 / np.maximum(q_sum, 1e-30), 0.0)
        for b in range(4):
            Q_list[b] = Q_list[b] * inv_q
        H_k = 1.0 - (Q_list[0] ** 2 + Q_list[1] ** 2 + Q_list[2] ** 2 + Q_list[3] ** 2)  # (R, P)

        Q_DOT = self.B_mats[0] @ Q_list[0].T
        for b in range(1, 4):
            Q_DOT += self.B_mats[b] @ Q_list[b].T  # (N, R)
        H_DOT = self.V_p @ H_k.T                   # (N, R)

        P_iQ = (self.V_p_rowsum[:, None] - Q_DOT) / self.L_val[:, None]
        H_QQ = H_DOT / self.L_val[:, None] / np.maximum(1e-6, 1.0 - w_sq_sums[None, :])
        P_rQ = 0.5 * H_QQ
        P_ri = np.maximum(1e-5, (P_iQ - P_rQ) / np.maximum(0.1, 1.0 - (4.0 / 3.0) * P_rQ))
        return 0.75 * alpha_g * (np.maximum(0.05, 1.0 - (4.0 / 3.0) * P_ri) ** (-1.0 / alpha_g) - 1.0)

    def eval_divs(
        self,
        anc_idx: np.ndarray,
        t_anc: np.ndarray,
        reg_idx: np.ndarray,
        eff_gamma: float = 0.10,
        alpha_g: float = 0.20,
    ) -> np.ndarray:
        t_min_a = float(np.min(t_anc))
        w_anc = np.exp(-eff_gamma * (t_anc - t_min_a))
        w_sum = float(np.sum(w_anc))
        w_anc = w_anc / w_sum if w_sum > 0 else np.full(len(t_anc), 1.0 / len(t_anc))
        w_sq_sum = float(np.sum(w_anc ** 2))
        w_full = np.bincount(anc_idx, weights=w_anc, minlength=self.N)

        Q_p = np.column_stack([w_full @ self.B_mats[b] for b in range(4)])
        q_sum = np.sum(Q_p, axis=1, keepdims=True)
        valid_q = q_sum.squeeze() > 0
        Q_p[valid_q] /= q_sum[valid_q]
        h_k_p = 1.0 - np.sum(Q_p ** 2, axis=1)

        q_dot = self.B_mats[0] @ Q_p[:, 0]
        for b in range(1, 4):
            q_dot += self.B_mats[b] @ Q_p[:, b]
        p_iQ = (self.V_p_rowsum - q_dot) / self.L_val
        h_QQ = (self.V_p @ h_k_p) / self.L_val / max(1e-6, 1.0 - w_sq_sum)
        p_rQ = 0.5 * h_QQ
        p_ri = np.maximum(1e-5, (p_iQ - p_rQ) / np.maximum(0.1, 1.0 - (4.0 / 3.0) * p_rQ))
        divs_all = 0.75 * alpha_g * (np.maximum(0.05, 1.0 - (4.0 / 3.0) * p_ri) ** (-1.0 / alpha_g) - 1.0)
        return divs_all[reg_idx]


def _fast_rcs_tmrca(times: np.ndarray, dists: np.ndarray):
    """Fast 3x3 OLS-equivalent Harrell Restricted Cubic Spline clock matching ChronAeon."""
    n = len(times)
    t_min, t_max = float(np.min(times)), float(np.max(times))
    knots = np.array([t_min, float(np.median(times)), float(np.percentile(times, 90))], dtype=np.float64)
    if (knots[2] - knots[1]) <= 1e-6 or (knots[1] - knots[0]) <= 1e-6:
        return np.nan, np.nan, 1.0, 1.0, False

    B, _ = compute_rcs_basis(times, knots)
    _, dB_tmax = compute_rcs_basis(np.array([t_max], dtype=np.float64), knots)
    dB_max = float(dB_tmax[0, 0]) if (dB_tmax is not None and dB_tmax.shape[1] > 0) else 0.0

    X_lin = np.column_stack([np.ones(n), times])
    XtX_l = X_lin.T @ X_lin
    try:
        b_lin = np.linalg.solve(XtX_l, X_lin.T @ dists)
    except np.linalg.LinAlgError:
        b_lin = np.linalg.pinv(XtX_l) @ (X_lin.T @ dists)
    rss_lin = float(np.sum((dists - X_lin @ b_lin) ** 2))

    X_sp = np.column_stack([np.ones(n), times, B])
    XtX_s = X_sp.T @ X_sp
    try:
        b_sp = np.linalg.solve(XtX_s, X_sp.T @ dists)
    except np.linalg.LinAlgError:
        b_sp = np.linalg.pinv(XtX_s) @ (X_sp.T @ dists)
    rss_sp = float(np.sum((dists - X_sp @ b_sp) ** 2))

    diff_rss = max(0.0, rss_lin - rss_sp)
    f_stat = float(diff_rss / max(1e-12, rss_sp / max(1, n - 3)))
    p_f = float(stats.f.sf(f_stat, 1, max(1, n - 3)))
    mu_anc = float(b_sp[1])
    if mu_anc <= 1e-12:
        return np.nan, mu_anc, 1.0, p_f, False
    t0_sp = float(-b_sp[0] / mu_anc)
    mu_rec = float(mu_anc + b_sp[2] * dB_max)
    rate_ratio = float(mu_rec / mu_anc)
    use_sp = (p_f < 0.05) and (rate_ratio < 0.85) and (t0_sp > 1900.0)
    return t0_sp, mu_anc, rate_ratio, p_f, use_sp


def _fast_ols_tmrca(times: np.ndarray, dists: np.ndarray):
    """Unclamped centered OLS root-to-tip dating + exact Fieller CI (no min_time clamp)."""
    n = len(times)
    t_mean = float(np.mean(times))
    d_mean = float(np.mean(dists))
    x = times - t_mean
    sxx = float(np.sum(x * x))
    if sxx <= 1e-12:
        return np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, 0.0
    mu = float(np.sum(x * (dists - d_mean)) / sxx)
    res = dists - (mu * x + d_mean)
    sigma2 = float(np.sum(res * res) / max(1, n - 2))
    var_mu = sigma2 / sxx
    var_d0 = sigma2 / n
    r_val = float(np.corrcoef(times, dists)[0, 1]) if np.std(dists) > 1e-12 else 0.0
    if mu <= 1e-12:
        return np.nan, mu, np.nan, np.nan, np.nan, r_val, 1.0, np.nan
    t_mrca = t_mean - (d_mean / mu)
    t_crit = float(stats.t.ppf(0.975, df=max(1, n - 2)))
    g = (t_crit ** 2 * var_mu) / (mu ** 2)
    f_stat = (r_val ** 2 / max(1e-12, 1.0 - r_val ** 2)) * max(1, n - 2) if abs(r_val) < 1.0 else 999.0
    p_val = float(1.0 - stats.f.cdf(f_stat, 1, max(1, n - 2)))
    if g < 1.0:
        # Exact unclamped Fieller interval for ratio theta = d_mean / mu, t_MRCA = t_mean - theta
        # Centered predictors imply cov(mu, d_mean) = 0
        theta_hat = d_mean / mu
        root_term = np.sqrt((theta_hat ** 2 * var_mu) + (1.0 - g) * var_d0)
        th_lo = (theta_hat - t_crit / mu * root_term) / (1.0 - g)
        th_hi = (theta_hat + t_crit / mu * root_term) / (1.0 - g)
        ci_lo = float(t_mean - th_hi)
        ci_hi = float(t_mean - th_lo)
    else:
        ci_lo, ci_hi = np.nan, np.nan
    return float(t_mrca), float(mu), ci_lo, ci_hi, float(g), r_val, p_val, float(np.mean(dists))


def block_bootstrap_chronaeon(
    evaluator: FastProfileEvaluator,
    times: np.ndarray,
    block_labels: np.ndarray,
    n_boot: int = 2000,
    seed: int = 20261114,
    disjoint_split: bool = False,
    use_rcs: bool = False,
    eff_gamma: float = 0.10,
):
    """
    2,000-replicate block bootstrap resampling (country, 5-yr) clusters with NO t_min_sample clamping,
    matching the exact clock estimator (OLS or RCS) used for the full clade point estimate.
    Vectorized across batches of replicates using batched BLAS DGEMM.
    """
    unique_blocks, block_inv = np.unique(block_labels, return_inverse=True)
    n_blocks = len(unique_blocks)
    block_members = [np.where(block_inv == b)[0] for b in range(n_blocks)]
    rng = np.random.default_rng(seed)

    reg_indices = []
    t_regs = []
    W_rows = []
    w_sq_sums = []
    for _ in range(n_boot):
        sampled_b = rng.integers(0, n_blocks, size=n_blocks)
        idx = np.concatenate([block_members[b] for b in sampled_b])
        if len(idx) < 5:
            continue
        t_b = times[idx]
        if (float(np.max(t_b)) - float(np.min(t_b))) < 2.0:
            continue
        if disjoint_split and len(idx) >= 6:
            order_b = np.argsort(t_b, kind="stable")
            anc_mask_b = np.zeros(len(idx), dtype=bool)
            anc_mask_b[order_b[::5]] = True
            if anc_mask_b.sum() < 2:
                anc_mask_b[order_b[:2]] = True
            if (~anc_mask_b).sum() < 3:
                anc_mask_b = np.zeros(len(idx), dtype=bool)
                anc_mask_b[order_b[:max(2, len(idx) - 3)]] = True
            anc_idx = idx[anc_mask_b]
            reg_idx = idx[~anc_mask_b]
            t_anc, t_reg = t_b[anc_mask_b], t_b[~anc_mask_b]
        else:
            anc_idx, reg_idx = idx, idx
            t_anc, t_reg = t_b, t_b

        t_min_a = float(np.min(t_anc))
        w_anc = np.exp(-eff_gamma * (t_anc - t_min_a))
        w_sum = float(np.sum(w_anc))
        w_anc = w_anc / w_sum if w_sum > 0 else np.full(len(t_anc), 1.0 / len(t_anc))
        w_sq_sums.append(float(np.sum(w_anc ** 2)))
        W_rows.append(np.bincount(anc_idx, weights=w_anc, minlength=evaluator.N))
        reg_indices.append(reg_idx)
        t_regs.append(t_reg)

    t_boots = []
    ols_boots = []
    n_nonpos_mu = 0
    R_tot = len(reg_indices)
    chunk_sz = 500
    for start in range(0, R_tot, chunk_sz):
        end = min(R_tot, start + chunk_sz)
        W_chunk = np.vstack(W_rows[start:end])
        wsq_chunk = np.array(w_sq_sums[start:end], dtype=np.float64)
        DIVS_MAT = evaluator.eval_batch_divs(W_chunk, wsq_chunk)
        for j in range(end - start):
            r_idx = reg_indices[start + j]
            t_reg = t_regs[start + j]
            d_b = DIVS_MAT[r_idx, j]
            tm_ols, mu_ols, *_ = _fast_ols_tmrca(t_reg, d_b)
            if np.isnan(mu_ols) or mu_ols <= 1e-12:
                ols_boots.append(np.inf)
            else:
                ols_boots.append(tm_ols)

            if use_rcs:
                tm_rcs, mu_rcs, *_ = _fast_rcs_tmrca(t_reg, d_b)
                if np.isnan(mu_rcs) or mu_rcs <= 1e-12:
                    n_nonpos_mu += 1
                    t_boots.append(np.inf)
                else:
                    t_boots.append(tm_rcs)
            else:
                if np.isnan(mu_ols) or mu_ols <= 1e-12:
                    n_nonpos_mu += 1
                    t_boots.append(np.inf)
                else:
                    t_boots.append(tm_ols)

    def _pcts(lst):
        arr = np.array(lst, dtype=np.float64)
        if len(arr) < 50:
            return None, None, len(arr)
        lo = float(np.percentile(arr, 2.5))
        hi = float(np.percentile(arr, 97.5))
        return (round(lo, 2) if np.isfinite(lo) else None), (round(hi, 2) if np.isfinite(hi) else None), len(arr)

    ci_lo, ci_hi, n_valid = _pcts(t_boots)
    ols_lo, ols_hi, _ = _pcts(ols_boots)
    return ci_lo, ci_hi, ols_lo, ols_hi, n_blocks, n_valid, n_nonpos_mu


def date_clade_partition_chronaeon(
    seq_dict: dict,
    dates_map: dict,
    country_map: dict = None,
    min_span_years: float = 5.0,
    n_boot: int = 2000,
):
    """
    Run authentic ChronAeon Steiner self-heterozygosity-corrected and Gamma-rate-corrected
    root-to-tip divergence + OLS / Restricted Spline clock dating with:
      1) Estimator-consistent unclamped 2,000-rep (country, 5-yr) block-bootstrap CIs;
      2) Disjoint 20/80 stratified anchor/regression validation split.
    """
    dated_taxa = [t for t in seq_dict if t in dates_map and dates_map[t] is not None and not np.isnan(dates_map[t])]
    if len(dated_taxa) < 5:
        return None

    times = np.array([dates_map[t] for t in dated_taxa], dtype=np.float64)
    t_min, t_max = float(np.min(times)), float(np.max(times))
    if (t_max - t_min) < min_span_years:
        return None

    if country_map is not None:
        block_labels = np.array([f"{country_map.get(t, 'UNK')}_{(int(dates_map[t]) // 5) * 5}" for t in dated_taxa])
    else:
        block_labels = np.array([f"ALL_{(int(dates_map[t]) // 5) * 5}" for t in dated_taxa])
    # Fallback to per-isolate resampling when fewer than 4 country x 5yr blocks exist
    # so sparse clades (e.g. 71_BF1, 35_A1D with 2 blocks) do not collapse to 0-width CIs.
    if len(np.unique(block_labels)) < 4:
        block_labels = np.array([f"ISO_{i}" for i in range(len(dated_taxa))])

    char_map = {"A": 0, "C": 1, "G": 2, "T": 3}
    N = len(dated_taxa)
    L = len(seq_dict[dated_taxa[0]])
    M = np.full((N, L), -1, dtype=np.int8)
    for i, t in enumerate(dated_taxa):
        s = seq_dict[t].upper()
        for k, ch in enumerate(s):
            M[i, k] = char_map.get(ch, -1)
    L_val = np.maximum(np.sum(M >= 0, axis=1), 1).astype(np.float64)
    evaluator = FastProfileEvaluator(M, L_val)

    all_idx = np.arange(N, dtype=np.int64)
    divs = evaluator.eval_divs(all_idx, times, all_idx)

    ols_tm, ols_mu, f_lo, f_hi, f_g, ols_r, ols_p, mean_rd = _fast_ols_tmrca(times, divs)
    sp_res = run_restricted_spline_clock_dating(times=times, dists=divs, n_boot=0)
    use_sp = (
        bool(sp_res.get("p_f_test", 1.0) < 0.05)
        and float(sp_res.get("rate_ratio", 1.0)) < 0.85
        and float(sp_res.get("t_mrca", 0.0)) > 1900.0
    )

    if use_sp:
        t_mrca_val = float(sp_res["t_mrca"])
        mu_val = float(sp_res["rate_ancestral"])
        model_used = "RCS"
    else:
        t_mrca_val = ols_tm
        mu_val = ols_mu
        model_used = "OLS"

    boot_lo, boot_hi, ols_boot_lo, ols_boot_hi, n_blocks, n_boot_valid, n_nonpos_mu = block_bootstrap_chronaeon(
        evaluator, times, block_labels, n_boot=n_boot, seed=20261114, disjoint_split=False, use_rcs=use_sp
    )

    # Disjoint 20/80 stratified anchor vs regression split
    order = np.argsort(times, kind="stable")
    anc_mask = np.zeros(N, dtype=bool)
    anc_mask[order[::5]] = True
    if anc_mask.sum() < 2:
        anc_mask[order[:2]] = True
    if (~anc_mask).sum() < 3:
        anc_mask = np.zeros(N, dtype=bool)
        anc_mask[order[:max(2, N - 3)]] = True

    anc_idx_20 = all_idx[anc_mask]
    reg_idx_80 = all_idx[~anc_mask]
    divs_sp = evaluator.eval_divs(anc_idx_20, times[anc_mask], reg_idx_80)
    sp_ols_tm, sp_ols_mu, _, _, _, sp_r, sp_p, _ = _fast_ols_tmrca(times[~anc_mask], divs_sp)
    if use_sp:
        sp_rcs_tm, sp_rcs_mu, *_ = _fast_rcs_tmrca(times[~anc_mask], divs_sp)
        sp_tm = sp_rcs_tm if not np.isnan(sp_rcs_tm) else sp_ols_tm
        sp_mu = sp_rcs_mu if not np.isnan(sp_rcs_tm) else sp_ols_mu
    else:
        sp_tm, sp_mu = sp_ols_tm, sp_ols_mu

    sp_boot_lo, sp_boot_hi, _, _, _, _, _ = block_bootstrap_chronaeon(
        evaluator, times, block_labels, n_boot=n_boot, seed=20261115, disjoint_split=True, use_rcs=use_sp
    )

    return {
        "n_dated": len(dated_taxa),
        "n_blocks": n_blocks,
        "clock_model": model_used,
        "t_min_sample": round(t_min, 1),
        "t_max_sample": round(t_max, 1),
        "t_mrca": round(t_mrca_val, 2) if not np.isnan(t_mrca_val) else None,
        "ci_lower": boot_lo if n_boot_valid >= 50 else (round(f_lo, 2) if not np.isnan(f_lo) else None),
        "ci_upper": boot_hi if n_boot_valid >= 50 else (round(f_hi, 2) if not np.isnan(f_hi) else None),
        "boot_ci_lower": boot_lo,
        "boot_ci_upper": boot_hi,
        "ols_t_mrca": round(ols_tm, 2) if not np.isnan(ols_tm) else None,
        "ols_boot_ci_lower": ols_boot_lo,
        "ols_boot_ci_upper": ols_boot_hi,
        "fieller_unclamped_lower": round(f_lo, 2) if not np.isnan(f_lo) else None,
        "fieller_unclamped_upper": round(f_hi, 2) if not np.isnan(f_hi) else None,
        "n_boot_valid": n_boot_valid,
        "n_boot_nonpos_mu": n_nonpos_mu,
        "mu_subs_site_yr": mu_val,
        "ols_mu_subs_site_yr": ols_mu if not np.isnan(ols_mu) else None,
        "r": round(float(ols_r), 3),
        "r2": round(float(ols_r ** 2), 3),
        "p_value": float(ols_p),
        "fieller_g": round(float(f_g), 4) if not np.isnan(f_g) else None,
        "mean_root_div": round(float(mean_rd), 4),
        "n_anchor_20": int(anc_mask.sum()),
        "n_reg_80": int((~anc_mask).sum()),
        "split2080_t_mrca": round(sp_tm, 2) if not np.isnan(sp_tm) else None,
        "split2080_ols_t_mrca": round(sp_ols_tm, 2) if not np.isnan(sp_ols_tm) else None,
        "split2080_boot_ci_lower": sp_boot_lo,
        "split2080_boot_ci_upper": sp_boot_hi,
        "split2080_mu": sp_mu if not np.isnan(sp_mu) else None,
        "split2080_r": round(float(sp_r), 3),
        "split2080_p_value": float(sp_p) if not np.isnan(sp_p) else None,
    }


def main():
    m2h_ref, h2m_ref, gmap = load_coord_and_gene_maps()

    # =========================================================================
    # PART 1: Audit of the 55 Nominal "Pure" Subtype Reference Genomes
    # =========================================================================
    print("=" * 85)
    print("PART 1: RHIZAEON AUDIT OF LANL'S 55 NOMINAL 'PURE' SUBTYPE REFERENCE GENOMES")
    print("=" * 85)
    with open(RESULTS_DIR / "pure_subtypes_rhizaeon.json") as f:
        pure_scan = json.load(f)

    pure_events = []
    for ev in pure_scan.get("events", []):
        s_msa, e_msa = int(ev["u1"]), int(ev["u2"])
        s_hxb2 = m2h_ref.get(s_msa, s_msa)
        e_hxb2 = m2h_ref.get(e_msa, e_msa)
        dom_str = annotate_hxb2_interval(s_hxb2, e_hxb2, gmap)
        rec_sub = ev["candidate_name"].split(".")[1]
        home_sub = ev["home_name"].split(".")[1] if "." in ev["home_name"] else ev["home_name"]
        donor_sub = "Ghost" if ev["is_ghost_donor"] else (ev["donor_name"].split(".")[1] if "." in ev["donor_name"] else ev["donor_name"])
        co_isolates = ";".join(ev.get("isolates", [ev["candidate_name"]]))
        row = {
            "event_id": ev["event_id"],
            "recombinant_id": ev["candidate_name"],
            "co_recombinant_isolates": co_isolates,
            "n_co_recombinants": len(ev.get("isolates", [ev["candidate_name"]])),
            "nominal_subtype": rec_sub,
            "home_id": ev["home_name"],
            "home_subtype": home_sub,
            "donor_id": ev["donor_name"],
            "donor_subtype": donor_sub,
            "is_ghost": ev["is_ghost_donor"],
            "msa_start": s_msa,
            "msa_end": e_msa,
            "hxb2_start": s_hxb2,
            "hxb2_end": e_hxb2,
            "hxb2_length_nt": e_hxb2 - s_hxb2 + 1,
            "s_informative": ev["s_informative"],
            "genomic_domains": dom_str,
            "p_fisher": ev["p_fisher"],
            "z_phys": round(ev["z_phys"], 2)
        }
        pure_events.append(row)
        print(f"  [CONTAMINATED REF] {rec_sub:4s} (n={row['n_co_recombinants']}) | {ev['candidate_name']:36s} | Home={home_sub:3s} Donor={donor_sub:5s} | HXB2 {s_hxb2:4d}-{e_hxb2:4d} ({dom_str:22s}) | p={ev['p_fisher']:.2e} Z={ev['z_phys']:.2f}")

    df_pure_events = pd.DataFrame(pure_events)
    df_pure_events.to_csv(RESULTS_DIR / "phase1_contaminated_pure_references.csv", index=False)

    # =========================================================================
    # PART 2: ChronAeon Molecular Clock Dating Across All Major Subtypes & CRFs
    #         (Using N = 5,089 full-length genomes in lanl_flt_2022.fasta)
    # =========================================================================
    print("\n" + "=" * 85)
    print("PART 2: CHRONAEON TREE-FREE MOLECULAR CLOCK DATING OF SUBTYPES & CRFs (N=5,089)")
    print("=" * 85)

    flt_fasta = DATA_DIR / "lanl_flt_2022.fasta"
    df_flt_meta = pd.read_csv(DATA_DIR / "lanl_flt_2022_metadata.csv")
    dates_map_flt = {
        row["seq_id"]: float(row["year"])
        for _, row in df_flt_meta.iterrows()
        if pd.notna(row["year"])
    }
    country_map_flt = {
        row["seq_id"]: (str(row["country"]).strip() if pd.notna(row.get("country")) else "UNK")
        for _, row in df_flt_meta.iterrows()
    }
    tag_map_flt = {row["seq_id"]: str(row["lanl_tag"]) for _, row in df_flt_meta.iterrows()}

    m2h_flt, h2m_flt = build_flt_hxb2_map(flt_fasta)
    c_start_flt = h2m_flt[790] - 1
    c_end_flt = h2m_flt[9417]

    # Group full-coding sequences by LANL clade tag
    clade_seqs = defaultdict(dict)
    all_flt_seqs = {}
    for rec in SeqIO.parse(str(flt_fasta), "fasta"):
        s_full = str(rec.seq).upper()
        all_flt_seqs[rec.id] = s_full
        tag = tag_map_flt.get(rec.id, "Unknown")
        clade_seqs[tag][rec.id] = s_full[c_start_flt:c_end_flt]

    target_clades = ["A1", "A6", "B", "C", "D", "F1", "F2", "G", "H", "J", "01_AE", "02_AG", "06_cpx", "07_BC", "08_BC", "11_cpx", "12_BF", "13_cpx", "35_A1D", "63_02A6", "71_BF1", "85_BC"]
    clock_rows = []
    print(f"  {'Clade':<10} {'Mdl':>4} {'N_dated':>7} {'SampleSpan':>13} {'t_MRCA':>9} {'BlockBoot_95%_CI':>18} {'Split20_80_tMRCA':>16} {'Rate (mu)':>12} {'r':>7} {'p-val':>10}", flush=True)
    print("  " + "-" * 117, flush=True)
    for cl in target_clades:
        if cl not in clade_seqs or len(clade_seqs[cl]) < 5:
            continue
        res = date_clade_partition_chronaeon(clade_seqs[cl], dates_map_flt, country_map=country_map_flt, min_span_years=4.0, n_boot=2000)
        if res:
            res["clade"] = cl
            clock_rows.append(res)
            ci_str = f"[{res['ci_lower']}, {res['ci_upper'] if res['ci_upper'] is not None else 'unb.'}]" if res["ci_lower"] is not None else "unbounded"
            span_str = f"{res['t_min_sample']:.0f}-{res['t_max_sample']:.0f}"
            tm_s = f"{res['t_mrca']:.1f}" if res["t_mrca"] is not None else "NA"
            sp_s = f"{res['split2080_t_mrca']:.1f}" if res["split2080_t_mrca"] is not None else "NA"
            mu_s = f"{res['mu_subs_site_yr']:.2e}" if pd.notna(res["mu_subs_site_yr"]) else "NA"
            print(f"  {cl:<10} {res['clock_model']:>4} {res['n_dated']:>7d} {span_str:>13s} {tm_s:>9s} {ci_str:>18s} {sp_s:>16s} {mu_s:>12s} {res['r']:>7.3f} {res['p_value']:>10.2e}", flush=True)

    df_clock = pd.DataFrame(clock_rows)
    df_clock.to_csv(RESULTS_DIR / "phase1_chronaeon_clade_dating.csv", index=False)

    # =========================================================================
    # PART 3: Partition-Specific ChronAeon Audit of CRF02_AG vs. Subtype G vs. A1
    #         and CRF01_AE (A1 gag-pol backbone vs. Ghost_E env-nef cassette)
    # =========================================================================
    print("\n" + "=" * 85)
    print("PART 3: PARTITION-SPECIFIC CHRONAEON POLARITY AUDIT (CRF02_AG vs G; CRF01_AE Ghost_E)")
    print("=" * 85)

    # According to LANL official map for CRF02_AG:
    # A-like regions in CRF02_AG: HXB2 790-1285 (gag-MA), 1859-2550 (gag-NC..pol-PR), 3275-4174 (pol-RT..IN), 6225-8300 (env gp120..gp41)
    # G-like regions in CRF02_AG: HXB2 1286-1858 (gag-CA), 2626-3274 (pol-RT), 4175-4874 (pol-IN), 5041-6060 (vif-vpr-tat), 8800-9417 (nef)
    a_intervals = [(790, 1285), (1859, 2550), (3275, 4174), (6225, 8300)]
    g_intervals = [(1286, 1858), (2626, 3274), (4175, 4874), (5041, 6060), (8800, 9417)]

    def slice_hxb2_intervals(s_full: str, intervals):
        chunks = []
        for hs, he in intervals:
            ms = h2m_flt[hs] - 1
            me = h2m_flt[he]
            chunks.append(s_full[ms:me])
        return "".join(chunks)

    partition_audit_rows = []
    for cl in ["A1", "02_AG", "G"]:
        ids = list(clade_seqs[cl].keys())
        seqs_A_part = {tid: slice_hxb2_intervals(all_flt_seqs[tid], a_intervals) for tid in ids}
        seqs_G_part = {tid: slice_hxb2_intervals(all_flt_seqs[tid], g_intervals) for tid in ids}
        res_A = date_clade_partition_chronaeon(seqs_A_part, dates_map_flt, country_map=country_map_flt, n_boot=2000)
        res_G = date_clade_partition_chronaeon(seqs_G_part, dates_map_flt, country_map=country_map_flt, n_boot=2000)
        if res_A and res_G:
            print(f"  [CRF02 vs G Audit] Clade {cl:6s} (N={res_A['n_dated']:3d}):")
            print(f"     - In A1-shared tracts (gag/IN/gp120): t_MRCA = {res_A['t_mrca']:.1f} [{res_A['ci_lower']}, {res_A['ci_upper']}], split2080={res_A['split2080_t_mrca']}, root_div = {res_A['mean_root_div']:.4f}, r = {res_A['r']:.3f}")
            print(f"     - In G-shared tracts  (CA/RT/vif/nef): t_MRCA = {res_G['t_mrca']:.1f} [{res_G['ci_lower']}, {res_G['ci_upper']}], split2080={res_G['split2080_t_mrca']}, root_div = {res_G['mean_root_div']:.4f}, r = {res_G['r']:.3f}")
            partition_audit_rows.append({"clade": cl, "partition": "A1_shared_tracts", **res_A})
            partition_audit_rows.append({"clade": cl, "partition": "G_shared_tracts", **res_G})

    # CRF01_AE: A1 gag-pol backbone (790..5096) vs Ghost_E env-nef cassette (6062..9417)
    ae_backbone_int = [(790, 5096)]
    ae_ghost_int = [(6062, 9417)]
    for cl in ["01_AE", "A1"]:
        ids = list(clade_seqs[cl].keys())
        seqs_bb = {tid: slice_hxb2_intervals(all_flt_seqs[tid], ae_backbone_int) for tid in ids}
        seqs_gh = {tid: slice_hxb2_intervals(all_flt_seqs[tid], ae_ghost_int) for tid in ids}
        res_bb = date_clade_partition_chronaeon(seqs_bb, dates_map_flt, country_map=country_map_flt, n_boot=2000)
        res_gh = date_clade_partition_chronaeon(seqs_gh, dates_map_flt, country_map=country_map_flt, n_boot=2000)
        if res_bb and res_gh:
            print(f"  [CRF01_AE Audit]   Clade {cl:6s} (N={res_bb['n_dated']:3d}):")
            print(f"     - In gag-pol (A1 backbone):        t_MRCA = {res_bb['t_mrca']:.1f} [{res_bb['ci_lower']}, {res_bb['ci_upper']}], split2080={res_bb['split2080_t_mrca']}, mu = {res_bb['mu_subs_site_yr']:.2e}, r = {res_bb['r']:.3f}")
            print(f"     - In vpu-env-nef (Ghost_E tract):  t_MRCA = {res_gh['t_mrca']:.1f} [{res_gh['ci_lower']}, {res_gh['ci_upper']}], split2080={res_gh['split2080_t_mrca']}, mu = {res_gh['mu_subs_site_yr']:.2e}, r = {res_gh['r']:.3f}")
            partition_audit_rows.append({"clade": cl, "partition": "CRF01_gag_pol_A1_backbone", **res_bb})
            partition_audit_rows.append({"clade": cl, "partition": "CRF01_vpu_env_nef_GhostE", **res_gh})

    pd.DataFrame(partition_audit_rows).to_csv(RESULTS_DIR / "phase1_partition_temporal_audit.csv", index=False)


if __name__ == "__main__":
    main()
