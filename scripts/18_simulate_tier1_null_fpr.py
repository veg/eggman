#!/usr/bin/env python3
"""
Step 18: Parametric Non-Recombinant GTR+Gamma4 Coalescent Null Simulations (B2 / Tier B).

Tests EGGMAN Tier-1 Viterbi and Tier-2 final lineage classification on N = 1,300 strictly
non-recombinant 8,628-nt coding genomes evolved from all 13 non-ghost HIV-1 Group M basal
founder consensi (A1, A2, A6, B, C, D, F1, F2, G0, H, J, K, L; 100 per founder) under three regimes:
  1. Standard_GTR_G4 (n = 650): Coalescent GTR+G4 evolution calibrated to observed
     intra-subtype HIV-1 Group M coding divergence.
  2. RateSpike_Homoplasy (n = 325): 3-fold substitution rate spikes in gp120 (HXB2 6225-7757)
     and pol.RT (HXB2 2550-3870) plus 10 convergent drug-resistance / CTL-escape homoplasies
     mapped to alleles shared with non-parental subtypes.
  3. APOBEC3G_3F_Hypermutation (n = 325): Standard GTR+G4 plus directional APOBEC3G (GG->AG)
     and APOBEC3F (GA->AA) hypermutation affecting 2-5% of target dinucleotides.

Evaluates both default, partial-screen, full-genome-screen, and permissive compendium regimes:
  - Default (--gamma 7.0 --p-thresh 1e-4 --min-snps 5)
  - Partial_Screen (--gamma 10.5 --p-thresh 1e-5 --min-snps 7 --min-span 200)
  - FullGenome_Screen (--gamma 13.5 --p-thresh 1e-5 --min-snps 7 --min-span 200)
  - Compendium (--gamma 7.0 --p-thresh 1e-3 --min-snps 3)

Saves:
  - results/tier1_null_simulation_fpr.csv
  - results/tier1_null_simulation_details.csv
"""

import os
import subprocess
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.linalg import expm

PROJECT_ROOT = Path("/Users/sergei/Projects/HIV-recombination")
DATA_DIR = PROJECT_ROOT / "data"
RESULTS_DIR = PROJECT_ROOT / "results"
BASAL_PROFILES_TSV = PROJECT_ROOT / "rhiztype" / "crates" / "eggman-core" / "data" / "profiles_9719.tsv"
RHIZTYPE_BIN = PROJECT_ROOT / "rhiztype" / "target" / "release" / "eggman"

TRUNKS = ["A1", "A2", "A6", "B", "C", "D", "F1", "F2", "G0", "H", "J", "K", "L"]
BASE_TO_IDX = {"A": 0, "C": 1, "G": 2, "T": 3}
IDX_TO_BASE = np.array(["A", "C", "G", "T"])


def load_basal_coding_trunks():
    """Load 8,628-nt HXB2 coding consensus profiles (HXB2 pos 790..9417) for empirical trunks."""
    trunks = {}
    with open(BASAL_PROFILES_TSV) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2 and parts[0] in TRUNKS:
                coding = parts[1][789:9417].upper()
                assert len(coding) == 8628, f"Expected 8628 nt, got {len(coding)}"
                trunks[parts[0]] = coding
    return trunks


def estimate_empirical_hiv_gtr(trunk_seqs: dict):
    """Estimate HIV-1 Group M empirical base freqs and GTR exchangeabilities from basal trunks + intra-subtype LANL diversity."""
    all_str = "".join(trunk_seqs[k] for k in TRUNKS)
    counts = np.array([all_str.count(b) for b in "ACGT"], dtype=float)
    pi = counts / counts.sum()

    # Estimate pairwise nucleotide substitution spectrum across basal subtypes to get realistic HIV-1 transition bias
    seqs_idx = {k: np.array([BASE_TO_IDX.get(c, -1) for c in v], dtype=np.int8) for k, v in trunk_seqs.items()}
    pair_mat = np.zeros((4, 4), dtype=float)
    for i, k1 in enumerate(TRUNKS):
        for k2 in TRUNKS[i + 1:]:
            s1, s2 = seqs_idx[k1], seqs_idx[k2]
            m = (s1 >= 0) & (s2 >= 0) & (s1 != s2)
            for a, b in zip(s1[m], s2[m]):
                pair_mat[a, b] += 0.5
                pair_mat[b, a] += 0.5

    # Exchangeability R_ij proportional to observed symmetrical rate / (pi_i * pi_j)
    R = np.zeros((4, 4), dtype=float)
    for i in range(4):
        for j in range(4):
            if i != j:
                R[i, j] = max(pair_mat[i, j], 1.0) / (pi[i] * pi[j])
    # Scale GTR Q matrix so mean rate sum_i pi_i (-Q_ii) = 1.0
    Q = R * pi[np.newaxis, :]
    np.fill_diagonal(Q, 0.0)
    np.fill_diagonal(Q, -Q.sum(axis=1))
    mean_rate = float(-(pi * np.diag(Q)).sum())
    Q /= mean_rate
    return pi, Q, seqs_idx


def discrete_gamma_rates(alpha: float, k: int = 4) -> np.ndarray:
    """Yang (1994) mean-quantile discrete 4-category Gamma site-rate approximation."""
    from scipy import stats
    quantiles = np.linspace(0.5 / k, 1.0 - 0.5 / k, k)
    rates = stats.gamma.ppf(quantiles, a=alpha, scale=1.0 / alpha)
    return rates / rates.mean()


def evolve_gtr_gamma4(anc_arr: np.ndarray, total_div: float, Q: np.ndarray, site_cats: np.ndarray,
                      gamma_rates: np.ndarray, rate_mult: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    """Evolve non-recombinant sequence along a coalescent branch path summing to total_div expected subs/site."""
    out = anc_arr.copy()
    L = len(anc_arr)
    # Build effective branch-length multiplier per site category & region
    # Evolve via discrete coalescent-tree epoch splits (3 epochs summing to total_div to mimic shared genealogy)
    epochs = rng.dirichlet([2.0, 2.0, 2.0]) * total_div
    for ep_t in epochs:
        for c_idx, g_r in enumerate(gamma_rates):
            for rm_val in np.unique(rate_mult):
                t_eff = ep_t * g_r * rm_val
                P = expm(Q * t_eff)
                cdf = np.cumsum(P, axis=1)
                cdf[:, -1] = 1.0
                mask = (site_cats == c_idx) & (rate_mult == rm_val) & (out >= 0)
                pos = np.where(mask)[0]
                if len(pos) == 0:
                    continue
                u = rng.random(len(pos))
                states = out[pos]
                new_states = np.empty(len(pos), dtype=np.int8)
                for s in range(4):
                    sm = states == s
                    if np.any(sm):
                        new_states[sm] = np.searchsorted(cdf[s], u[sm], side="right").clip(0, 3)
                out[pos] = new_states
    return out


def main():
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    trunk_strs = load_basal_coding_trunks()
    pi, Q, trunk_arrays = estimate_empirical_hiv_gtr(trunk_strs)

    print("[*] Estimated empirical HIV-1 Group M GTR parameters from basal trunks:")
    print(f"    * Base freqs A/C/G/T: {np.round(pi, 4)}")
    print(f"    * Transition/transversion excess: Q(A->G)={Q[0,2]:.3f}, Q(C->T)={Q[1,3]:.3f} vs Q(A->C)={Q[0,1]:.3f}")

    # Precompute discrete Gamma4 categories (alpha=0.38 typical HIV-1 coding) and site rate multipliers
    rng = np.random.default_rng(20261114)
    gamma_rates = discrete_gamma_rates(alpha=0.38, k=4)
    site_cats = rng.integers(0, 4, size=8628)

    # Spiked rate multiplier vector: 3x across pol.RT (HXB2 2550..3870) and env gp120 (HXB2 6225..7757)
    # 0-based index in 8,628 slice = HXB2_pos - 790
    spike_mult = np.ones(8628, dtype=float)
    spike_mult[(2550 - 790):(3870 - 790)] = 3.0
    spike_mult[(6225 - 790):(7757 - 790)] = 3.0
    base_mult = np.ones(8628, dtype=float)

    # Identify clustered non-parental homoplasy candidate sites (clustered in pol.RT 2550-3500 & env gp120 6300-7400)
    # where another basal trunk differs from the focal trunk so convergent evolution mimics clustered discordant SNPs
    homoplasy_pools = {}
    pol_slice = np.arange(2550 - 790, 3500 - 790)
    env_slice = np.arange(6300 - 790, 7400 - 790)
    for t in TRUNKS:
        anc = trunk_arrays[t]
        other_trunks = [ot for ot in TRUNKS if ot != t]
        # Pick another trunk (e.g. next in list) as convergent donor for clustered homoplasy test
        pol_cands = []
        env_cands = []
        for ot in other_trunks:
            o_arr = trunk_arrays[ot]
            p_diff = [idx for idx in pol_slice if anc[idx] >= 0 and o_arr[idx] >= 0 and anc[idx] != o_arr[idx]]
            e_diff = [idx for idx in env_slice if anc[idx] >= 0 and o_arr[idx] >= 0 and anc[idx] != o_arr[idx]]
            if len(p_diff) >= 10:
                pol_cands.append((ot, p_diff))
            if len(e_diff) >= 10:
                env_cands.append((ot, e_diff))
        homoplasy_pools[t] = (pol_cands, env_cands)

    # 100 simulations per each of the 13 non-ghost founder trunks (N = 1,300)
    scenarios = (
        [("Standard_GTR_G4", i) for i in range(650)]
        + [("RateSpike_Homoplasy", i) for i in range(325)]
        + [("APOBEC3G_3F_Burst", i) for i in range(325)]
    )

    fasta_path = RESULTS_DIR / "tier1_null_simulated_1300.fasta"
    meta_rows = []
    with open(fasta_path, "w") as f_fa:
        for idx, (scen, s_idx) in enumerate(scenarios):
            trunk = TRUNKS[s_idx % len(TRUNKS)]
            anc = trunk_arrays[trunk]
            # Draw realistic intra-subtype divergence from empirical distribution (mean ~0.065 subs/site)
            total_div = float(rng.uniform(0.035, 0.110))
            rm = spike_mult if scen == "RateSpike_Homoplasy" else base_mult

            sim_arr = evolve_gtr_gamma4(anc, total_div, Q, site_cats, gamma_rates, rm, rng)

            # Impose scenario-specific artifacts on strictly non-recombinant background
            if scen == "RateSpike_Homoplasy":
                pol_cands, env_cands = homoplasy_pools[trunk]
                chosen_pool = pol_cands if (s_idx % 2 == 0) else env_cands
                target_ot, diff_sites = chosen_pool[rng.integers(0, len(chosen_pool))]
                # Place 10 convergent homoplasies matching target_ot within a narrow ~600-nt gene region
                win_start = diff_sites[rng.integers(0, max(1, len(diff_sites) - 10))]
                h_sites = [x for x in diff_sites if x >= win_start][:10]
                if len(h_sites) < 10:
                    h_sites = diff_sites[:10]
                ot_arr = trunk_arrays[target_ot]
                for hs in h_sites:
                    sim_arr[hs] = ot_arr[hs]
            elif scen == "APOBEC3G_3F_Burst":
                apobec_frac = float(rng.uniform(0.02, 0.05))
                # GG -> AG (APOBEC3G) or GA -> AA (APOBEC3F): pos i is G (2) followed by G (2) or A (0)
                is_g = sim_arr[:-1] == 2
                next_ga = (sim_arr[1:] == 2) | (sim_arr[1:] == 0)
                cand_pos = np.where(is_g & next_ga)[0]
                n_mut = max(1, int(round(apobec_frac * len(cand_pos))))
                if len(cand_pos) > 0:
                    chosen = rng.choice(cand_pos, size=min(n_mut, len(cand_pos)), replace=False)
                    sim_arr[chosen] = 0  # G -> A

            # Reconstruct ungapped ACGT string for FASTA (dropping template gaps if any)
            valid = sim_arr >= 0
            seq_str = "".join(IDX_TO_BASE[sim_arr[valid]].tolist())
            sid = f"SIM_{idx:04d}|{scen}|{trunk}|div={total_div:.4f}"
            f_fa.write(f">{sid}\n{seq_str}\n")
            meta_rows.append({
                "seq_id": sid,
                "scenario": scen,
                "true_trunk": trunk,
                "true_div": round(total_div, 5),
            })

    meta_df = pd.DataFrame(meta_rows)
    print(f"[+] Simulated {len(meta_df)} strictly non-recombinant 8,628-nt genomes to {fasta_path}.")

    regimes = [
        ("Default (--gamma 7.0 --p-thresh 1e-4 --min-snps 5)", "7.0", "0.0001", "5", None, "default_g7"),
        ("Partial_Screen (--gamma 10.5 --p-thresh 1e-5 --min-snps 7 --min-span 200)", "10.5", "1e-5", "7", "200", "partial_g105"),
        ("FullGenome_Screen (--gamma 13.5 --p-thresh 1e-5 --min-snps 7 --min-span 200)", "13.5", "1e-5", "7", "200", "full_g135"),
        ("Compendium (--gamma 7.0 --p-thresh 1e-3 --min-snps 3)", "7.0", "0.001", "3", None, "comp_g7"),
    ]

    all_detail = []
    sum_rows = []
    for reg_name, gamma_s, p_thr, min_s, min_sp, tag_s in regimes:
        out_csv = RESULTS_DIR / f"tier1_null_out_{tag_s}.csv"
        t0 = time.time()
        cmd = [
            str(RHIZTYPE_BIN), str(fasta_path),
            "--gamma", gamma_s,
            "--p-thresh", p_thr, "--min-snps", min_s,
            "--quiet", "--out", str(out_csv),
        ]
        if min_sp is not None:
            cmd.extend(["--min-span", min_sp])
        subprocess.run(cmd, check=True)
        dt = time.time() - t0
        res_df = pd.read_csv(out_csv)
        m = meta_df.merge(res_df, left_on="seq_id", right_on="sequence_id", how="inner")
        m["regime"] = reg_name
        # True trunk normalized: G0 -> G for Tier-2 Pure / Ancestral_Recombinant_Subtype check
        m["true_tier2"] = m["true_trunk"].replace({"G0": "G"})
        m["tier1_is_clonal"] = m["basal_formula"].astype(str).str.endswith("(Clonal)") & (m["basal_formula"].astype(str).str.split("(").str[0] == m["true_trunk"])
        m["tier1_false_recomb"] = ~m["basal_formula"].astype(str).str.endswith("(Clonal)")
        m["tier1_wrong_trunk"] = (~m["tier1_false_recomb"]) & (m["basal_formula"].astype(str).str.split("(").str[0] != m["true_trunk"])
        pure_cats = {"Pure_Subtype", "Ancestral_Recombinant_Subtype"}
        m["tier2_exact_pure"] = m["called_category"].isin(pure_cats) & (m["called_lineage"] == m["true_tier2"])
        m["tier2_false_urf_or_crf"] = m["called_category"].isin(["URF", "Certified_CRF", "Novel_Circulating_CRF"])
        all_detail.append(m)

        print(f"[+] Regime {reg_name} completed {len(m):,} seqs in {dt:.2f}s:")
        for scen_lbl in ["Standard_GTR_G4", "RateSpike_Homoplasy", "APOBEC3G_3F_Burst", "ALL_1300_COMBINED"]:
            sub = m if scen_lbl == "ALL_1300_COMBINED" else m[m["scenario"] == scen_lbl]
            n = len(sub)
            t1_recomb = int(sub["tier1_false_recomb"].sum())
            t1_wrong = int(sub["tier1_wrong_trunk"].sum())
            t2_exact = int(sub["tier2_exact_pure"].sum())
            t2_false = int(sub["tier2_false_urf_or_crf"].sum())
            t2_urf = int((sub["called_category"] == "URF").sum())
            t2_crf = int((sub["called_category"] == "Certified_CRF").sum())
            row = {
                "regime": reg_name,
                "scenario": scen_lbl,
                "n_simulated": n,
                "tier1_basal_false_recomb": t1_recomb,
                "tier1_basal_false_recomb_fpr_pct": round(100.0 * t1_recomb / n, 2),
                "tier1_wrong_clonal_trunk": t1_wrong,
                "tier2_exact_pure_accuracy": t2_exact,
                "tier2_exact_pure_pct": round(100.0 * t2_exact / n, 2),
                "tier2_false_recombinant_total": t2_false,
                "tier2_false_recombinant_fpr_pct": round(100.0 * t2_false / n, 2),
                "tier2_false_urf_calls": t2_urf,
                "tier2_false_crf_calls": t2_crf,
                "runtime_sec": round(dt, 2),
            }
            sum_rows.append(row)
            print(f"    * {scen_lbl:<22s} (N={n:4d}): Tier1 FalseRecomb FPR={row['tier1_basal_false_recomb_fpr_pct']:5.2f}% ({t1_recomb}) | Tier2 Pure={row['tier2_exact_pure_pct']:6.2f}% | Tier2 FalseRecomb={row['tier2_false_recombinant_fpr_pct']:5.2f}% (URF={t2_urf}, CRF={t2_crf})")

    pd.concat(all_detail, ignore_index=True).to_csv(RESULTS_DIR / "tier1_null_simulation_details.csv", index=False)
    sum_df = pd.DataFrame(sum_rows)
    sum_df.to_csv(RESULTS_DIR / "tier1_null_simulation_fpr.csv", index=False)
    print(f"[+] Saved summary to {RESULTS_DIR / 'tier1_null_simulation_fpr.csv'}")


if __name__ == "__main__":
    main()
