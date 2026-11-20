#!/usr/bin/env python3
"""
Script 13: Generate Publication-Grade Figure 5 (`figures/fig5_viability_sieve_and_partial_scale_discovery.{png,pdf}`)
complying with BioVis Figure Design Rules (`figure_design.md`):
  - Okabe-Ito colorblind-safe palette
  - Sans-serif typography (Arial/Helvetica), L-frame half-open axes (no top/right spines)
  - Log-scaled axes with 10^k notation
  - Panel A: Regional Strain Heterozygosity (H_circ) vs. De Novo URF Rate across 32 countries
  - Panel B: Offspring Cluster-Size CCDF Survival Curve (N = 749 recombinant lineages, Pareto alpha = 1.88)
  - Panel C: Two-Stage Epistatic Co-Ancestry Sieve (Singleton URFs vs Established CRFs vs Major Epidemic CRFs)
  - Panel D: Breakpoint Density across Functional Hinges vs Structural Cores + Partial pol/env Novel Lineages
"""

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

PROJECT_ROOT = Path("/Users/sergei/Projects/HIV-recombination")
RESULTS_DIR = PROJECT_ROOT / "results"
FIGURES_DIR = PROJECT_ROOT / "figures"
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

# Okabe-Ito colorblind-safe palette
OKABE_ITO = {
    "blue": "#0072B2",
    "orange": "#E69F00",
    "green": "#009E73",
    "vermilion": "#D55E00",
    "purple": "#CC79A7",
    "sky": "#56B4E9",
    "yellow": "#F0E442",
    "grey": "#666666",
    "black": "#111111",
}


def style_axis_lframe(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_linewidth(0.9)
    ax.spines["left"].set_color("#222222")
    ax.spines["bottom"].set_linewidth(0.9)
    ax.spines["bottom"].set_color("#222222")
    ax.tick_params(axis="both", which="major", labelsize=8.2, width=0.8, length=3.5, colors="#222222")
    ax.grid(False)


def main():
    plt.rcParams["font.family"] = "sans-serif"
    plt.rcParams["font.sans-serif"] = ["Arial", "Helvetica", "DejaVu Sans"]
    plt.rcParams["pdf.fonttype"] = 42
    plt.rcParams["ps.fonttype"] = 42

    reg_df = pd.read_csv(RESULTS_DIR / "regional_heterozygosity_urf_rates.csv")
    mv_df = pd.read_csv(RESULTS_DIR / "sieve_hcirc_multivariate_regression.csv")
    sieve_df = pd.read_csv(RESULTS_DIR / "recombination_viability_and_establishment_sieve.csv")
    man_df = pd.read_csv(RESULTS_DIR / "supp_manifest_sieve_749_lineages.csv")
    env_df = pd.read_csv(RESULTS_DIR / "scale_partial_env_novel_circulating_lineages.csv")

    fig, axes = plt.subplots(2, 2, figsize=(12.6, 9.0), dpi=300)
    plt.subplots_adjust(wspace=0.34, hspace=0.36, left=0.14, right=0.97, top=0.93, bottom=0.08)

    # =========================================================================
    # PANEL A: Regional Circulating Strain Heterozygosity vs De Novo URF Rate
    # =========================================================================
    ax_a = axes[0, 0]
    style_axis_lframe(ax_a)

    x_h = reg_df["h_circ_non_urf"].values if "h_circ_non_urf" in reg_df.columns else reg_df["h_circ_lineage_diversity"].values
    y_f = 100.0 * reg_df["f_urf_arch_per_patient"].values
    sizes_n = np.clip(reg_df["n_patients"].values * 0.35, 28, 300)

    african_hubs = {"Uganda", "Tanzania", "Cameroon", "Kenya", "Rwanda", "Zambia", "South Africa", "Democratic Republic of the Congo", "Senegal", "Ghana", "Nigeria", "Botswana"}
    colors_a = [
        OKABE_ITO["vermilion"] if c in african_hubs else OKABE_ITO["blue"]
        for c in reg_df["country"]
    ]

    ax_a.scatter(x_h, y_f, s=sizes_n, c=colors_a, alpha=0.84, edgecolors="white", linewidth=0.8, zorder=3)

    w_reg = reg_df["n_patients"].values.astype(float)
    beta_urf_slope = 100.0 * float(np.sum(w_reg * reg_df["f_urf_arch_per_patient"].values * x_h) / np.sum(w_reg * (x_h ** 2)))
    h_grid = np.linspace(0, 0.96, 100)
    ax_a.plot(
        h_grid, beta_urf_slope * h_grid, linestyle="--", color=OKABE_ITO["black"], linewidth=1.3,
        label=rf"Weighted slope $\hat{{\beta}}_{{\mathrm{{URF}}}} = {beta_urf_slope:.1f}\% / H_{{\mathrm{{circ}}}}$", zorder=2
    )

    # Clean non-overlapping country annotations with leader lines where needed
    annot_specs = {
        "Uganda":       (0.03, -0.5, "left"),
        "Tanzania":     (0.03, -2.5, "left"),
        "Cameroon":     (-0.03, 1.8, "right"),
        "China":        (-0.04, 3.4, "right"),
        "Brazil":       (0.03, -1.6, "left"),
        "Russia":       (0.03, -1.8, "left"),
        "Thailand":     (0.04, -1.4, "left"),
        "South Africa": (0.07,  3.2, "left"),
        "USA":          (0.11,  0.4, "left"),
    }
    for _, r in reg_df.iterrows():
        cty = r["country"]
        if cty in annot_specs:
            dx, dy, ha_align = annot_specs[cty]
            x0 = float(r["h_circ_non_urf"] if "h_circ_non_urf" in r else r["h_circ_lineage_diversity"])
            y0 = 100.0 * float(r["f_urf_arch_per_patient"])
            use_arrow = cty in ("South Africa", "USA", "China")
            ax_a.annotate(
                f"{cty} (n={int(r['n_patients'])})",
                xy=(x0, y0),
                xytext=(x0 + dx, y0 + dy),
                ha=ha_align, va="center",
                fontsize=7.1, fontweight="medium", color=OKABE_ITO["black"],
                arrowprops=dict(arrowstyle="-", color="#777777", lw=0.7, shrinkA=1, shrinkB=4) if use_arrow else None,
                bbox=dict(boxstyle="round,pad=0.12", fc="white", ec="none", alpha=0.82),
                zorder=5
            )

    ax_a.scatter([], [], s=45, c=OKABE_ITO["vermilion"], label="Sub-Saharan Africa (generalized)")
    ax_a.scatter([], [], s=45, c=OKABE_ITO["blue"], label="Eurasia & Americas (concentrated)")
    ax_a.set_xlim(-0.03, 1.05)
    ax_a.set_ylim(-1.2, 36.8)
    ax_a.set_xlabel(r"Non-URF Circulating Strain Heterozygosity ($H_{\mathrm{circ,non\text{-}URF}}$)", fontsize=9.2, fontweight="bold")
    ax_a.set_ylabel("De Novo URF Architectures (% of Patients)", fontsize=9.2, fontweight="bold")
    ax_a.set_title("A   Regional Co-Circulation Drives De Novo URF Generation", loc="left", fontsize=10.5, fontweight="bold", pad=8)
    wls_h_row = mv_df[(mv_df["model"] == "WLS_by_N") & (mv_df["term"] == "h_circ_non_urf")].iloc[0]
    wls_ssa_row = mv_df[(mv_df["model"] == "WLS_by_N") & (mv_df["term"] == "is_ssa")].iloc[0]
    from scipy.stats import spearmanr as _spr
    rho_h, p_rho_h = _spr(x_h, reg_df["f_urf_arch_per_patient"].values)
    ax_a.text(
        0.01, 35.6,
        rf"Spearman $\rho = {rho_h:.3f}$ ($p = {p_rho_h:.1e}$, $n = {len(reg_df)}$)" + "\n" +
        rf"WLS adj. $\beta_H = +{100*wls_h_row['coef']:.1f}\%$ ($p = {wls_h_row['p_value']:.1e}$)" + "\n" +
        rf"Sub-Saharan Africa effect $= +{100*wls_ssa_row['coef']:.1f}\%$ ($p = {wls_ssa_row['p_value']:.1e}$)",
        fontsize=6.9, fontweight="bold", color=OKABE_ITO["black"], va="top"
    )
    ax_a.legend(frameon=False, loc="upper left", bbox_to_anchor=(0.0, 0.77), fontsize=7.1)

    # =========================================================================
    # PANEL B: Offspring Cluster-Size CCDF & Transmission Establishment Funnel
    # =========================================================================
    ax_b = axes[0, 1]
    style_axis_lframe(ax_b)

    sizes_man = man_df["n_genbank_patients"].values.astype(float)
    n_tot_lin = len(sizes_man)
    thresholds = np.sort(np.unique(sizes_man))
    counts_ge = np.array([int(np.sum(sizes_man >= thr)) for thr in thresholds], dtype=float)
    ccdf_emp = counts_ge / float(n_tot_lin)
    alpha_mle = 1.0 + float(n_tot_lin) / float(np.sum(np.log(sizes_man / 0.5)))

    ax_b.step(thresholds, ccdf_emp, where="post", color=OKABE_ITO["blue"], linewidth=2.1, label=f"Empirical lineage survival (n = {n_tot_lin})", zorder=3)
    show_pts = np.isin(thresholds, [1, 2, 3, 4, 5, 7, 10, 15, 25, 50, 100, int(thresholds[-1])])
    ax_b.scatter(thresholds[show_pts], ccdf_emp[show_pts], color=OKABE_ITO["blue"], s=28, zorder=4, edgecolors="white", linewidths=0.4)

    x_pl = np.logspace(0, 2.95, 200)
    y_pl = x_pl ** (-(alpha_mle - 1.0))
    ax_b.plot(
        x_pl, y_pl, linestyle="--", color=OKABE_ITO["vermilion"], linewidth=1.5,
        label=rf"Yule–Pareto MLE ($\hat{{\alpha}} = {alpha_mle:.2f} \pm 0.03$)", zorder=2
    )

    n_ge3 = int(np.sum(sizes_man >= 3))
    n_ge10 = int(np.sum(sizes_man >= 10))
    n_ge10_can = int((man_df["n_genbank_patients"] >= 10).sum())
    ax_b.axvline(3, color=OKABE_ITO["green"], linestyle=":", linewidth=1.2)
    ax_b.annotate(
        rf"CRF Threshold ($N \geq 3$)" + f"\n{100*n_ge3/n_tot_lin:.1f}% ({n_ge3} / {n_tot_lin})",
        xy=(3.0, n_ge3 / n_tot_lin), xytext=(3.8, 0.36),
        fontsize=7.3, fontweight="bold", color=OKABE_ITO["green"],
        arrowprops=dict(arrowstyle="->", color=OKABE_ITO["green"], lw=0.8),
        bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.88),
        zorder=5
    )

    ax_b.axvline(10, color=OKABE_ITO["purple"], linestyle=":", linewidth=1.2)
    ax_b.annotate(
        rf"Survey Cluster ($N \geq 10$)" + f"\n{100*n_ge10/n_tot_lin:.1f}% ({n_ge10} / {n_tot_lin}; 22 canonical)",
        xy=(10.0, n_ge10 / n_tot_lin), xytext=(15.0, 0.025),
        fontsize=7.1, fontweight="bold", color="#9B4F78",
        arrowprops=dict(arrowstyle="->", color=OKABE_ITO["purple"], lw=0.8),
        bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.88),
        zorder=5
    )

    ax_b.set_xscale("log")
    ax_b.set_yscale("log")
    ax_b.xaxis.set_major_formatter(ticker.LogFormatterMathtext())
    ax_b.yaxis.set_major_formatter(ticker.LogFormatterMathtext())
    ax_b.set_xlim(0.85, 1300)
    ax_b.set_ylim(0.001, 1.35)
    ax_b.set_xlabel("Offspring Patient Cluster Size (N independent hosts)", fontsize=9.2, fontweight="bold")
    ax_b.set_ylabel(r"Complementary Cumulative Fraction $P(\mathrm{Size} \geq N)$", fontsize=9.2, fontweight="bold")
    ax_b.set_title("B   Recombinant Offspring Funnel & Establishment Sieve", loc="left", fontsize=10.5, fontweight="bold", pad=8)
    ax_b.legend(frameon=False, loc="upper right", fontsize=7.6)

    # =========================================================================
    # PANEL C: Two-Stage Epistatic Co-Ancestry Locks Across Epidemic Tiers
    # =========================================================================
    ax_c = axes[1, 0]
    style_axis_lframe(ax_c)

    locks = sieve_df[sieve_df["feature_category"] == "Intra_Complex_CoAncestry_Lock"].copy()
    lock_order = [
        "gp120 Core <-> gp41 Ectodomain",
        "Tat Exon 1 <-> Tat Exon 2",
        "Gag-CA/NC <-> Pol-PR Protease",
        "RRE Stem-Loop <-> Rev Exon 2",
        "Gag-MA Matrix <-> gp41-CT Tail",
    ]
    lock_labels = [
        "gp120 ↔ gp41-ecto\n(Trimer Interface)",
        "Tat Exon 1 ↔ Exon 2\n(2.4 kb Spliced Lock)",
        "Gag-CA/NC ↔ Pol-PR\n(Substrate–Protease)",
        "RRE ↔ Rev Exon 2\n(Nuclear RNA Export)",
        "Gag-MA ↔ gp41-CT\n(Opposite-End Split)",
    ]
    locks = locks.set_index("feature_name").loc[lock_order].reset_index()

    y_pos = np.arange(len(lock_order))
    urf_vals = 100.0 * locks["singleton_urf_rate"].values
    est_vals = 100.0 * locks["established_rate"].values
    pri_vals = 100.0 * locks["primary_dag_rate"].values
    maj_vals = 100.0 * locks["major_crf_rate"].values
    lin_null_vals = 100.0 * locks["est_lin_perm_null"].values
    null_val = 100.0 * float(locks["established_null"].iloc[0])

    ax_c.axvline(null_val, color="#AAAAAA", linestyle=":", linewidth=1.0, label=f"Unstructured Null ({null_val:.1f}%)")

    for i in range(len(lock_order)):
        lo = min(urf_vals[i], est_vals[i], pri_vals[i], maj_vals[i], lin_null_vals[i])
        hi = max(urf_vals[i], est_vals[i], pri_vals[i], maj_vals[i], lin_null_vals[i])
        ax_c.plot([lo, hi], [y_pos[i], y_pos[i]], color="#CCCCCC", linewidth=2.6, zorder=1)

    ax_c.scatter(lin_null_vals, y_pos, color=OKABE_ITO["black"], s=140, marker="|", linewidths=2.0, label="Spatial Null (10k seg.-perm.)", zorder=6)
    ax_c.scatter(urf_vals, y_pos, color=OKABE_ITO["sky"], s=56, label="Singleton Dead-End URFs (n = 461)", zorder=3, edgecolors="white", linewidths=0.4)
    ax_c.scatter(pri_vals, y_pos, color=OKABE_ITO["orange"], s=58, marker="^", label="Primary Gen-1 DAG CRFs (n = 89)", zorder=4, edgecolors="white", linewidths=0.4)
    ax_c.scatter(est_vals, y_pos, color=OKABE_ITO["blue"], s=66, marker="D", label="All Established CRFs (n = 154)", zorder=4, edgecolors="white", linewidths=0.4)
    ax_c.scatter(maj_vals, y_pos, color=OKABE_ITO["vermilion"], s=74, marker="s", label="Major Epidemic CRFs (n = 22)", zorder=5, edgecolors="white", linewidths=0.4)

    r0, r1, r2, r3, r4 = [locks.iloc[i] for i in range(5)]
    ax_c.text(max(est_vals[0], urf_vals[0]) + 1.8, y_pos[0], f"{est_vals[0]:.1f}% Est vs {lin_null_vals[0]:.1f}% lin. null\n(p < 1e-4)", va="center", fontsize=6.6, fontweight="bold", color=OKABE_ITO["black"])
    ax_c.text(max(maj_vals[1], est_vals[1]) + 1.8, y_pos[1], f"{maj_vals[1]:.1f}% Maj / {est_vals[1]:.1f}% Est vs {lin_null_vals[1]:.1f}% null\nCMH OR={r1['cmh_or_est_vs_urf']:.2f} (p={r1['cmh_p_est_vs_urf']:.1e})", va="center", fontsize=6.5, fontweight="bold", color=OKABE_ITO["vermilion"])
    ax_c.text(max(maj_vals[2], pri_vals[2]) + 1.8, y_pos[2], f"{maj_vals[2]:.1f}% Maj vs {lin_null_vals[2]:.1f}% lin. null (p={r2['p_est_vs_lin_ge']:.2f}, n.s.)", va="center", fontsize=6.6, fontweight="bold", color=OKABE_ITO["black"])
    ax_c.text(max(pri_vals[3], lin_null_vals[3]) + 1.8, y_pos[3], f"{est_vals[3]:.1f}% Est vs {lin_null_vals[3]:.1f}% null (n.s.)\nCMH vs URF p={r3['cmh_p_est_vs_urf']:.4f}", va="center", fontsize=6.5, fontweight="bold", color=OKABE_ITO["black"])
    ax_c.text(max(r4['primary_dag_lin_null']*100, pri_vals[4]) + 1.8, y_pos[4], f"{est_vals[4]:.1f}% Est (p={r4['p_est_vs_lin_le']:.3f})\n{pri_vals[4]:.1f}% DAG-Pri vs {r4['primary_dag_lin_null']*100:.1f}% null (n.s.)", va="center", fontsize=6.5, fontweight="bold", color=OKABE_ITO["black"])

    ax_c.set_yticks(y_pos)
    ax_c.set_yticklabels(lock_labels, fontsize=7.9, fontweight="bold")
    ax_c.invert_yaxis()
    ax_c.set_xlim(18, 138)
    ax_c.set_xlabel("Co-Ancestry Rate (% Genomes with Identical Parental Trunk)", fontsize=9.2, fontweight="bold")
    ax_c.set_title("C   Two-Stage Epistatic Co-Ancestry Sieve (Spatial Null)", loc="left", fontsize=10.5, fontweight="bold", pad=8)
    ax_c.legend(frameon=True, facecolor="white", edgecolor="#E2E2E2", framealpha=0.92, loc="upper left", fontsize=5.9, handletextpad=0.4, borderpad=0.3)

    # =========================================================================
    # PANEL D: Breakpoint Fold-Enrichment & Partial-Env gp41 TM/CT Hinge
    # =========================================================================
    ax_d = axes[1, 1]
    style_axis_lframe(ax_d)

    bp_specs = [
        ("Permissive_Hinge: gag.MA/CA Linker (1100-1260)", "gag.MA/CA Linker\n(1100–1260)"),
        ("Permissive_Hinge: gag/pol Frameshift & PR/RT (1980-2620)", "gag/pol & PR/RT\n(1980–2620)"),
        ("Permissive_Hinge: gp41 TM/CT & nef Boundary (8150-8950)", "gp41 TM/CT..nef\n(8150–8950)"),
        ("Permissive_Hinge: Accessory Hub (pol-IN..vpu, 4950-6350)", "Accessory Hub\n(4950–6350)"),
        ("Deleterious_Core: Gag-CA Capsid Helical Core (1280-1850)", "Gag-CA Core\n(1280–1850)"),
        ("Deleterious_Core: gp120 Ectodomain Beta-Sandwich (6450-7650)", "gp120 Core\n(6450–7650)"),
    ]
    bp_df = sieve_df[sieve_df["feature_category"] == "Breakpoint_Structural_Placement"].set_index("feature_name")
    z_labels = [lbl for _, lbl in bp_specs]
    f_maj = [float(bp_df.loc[k, "major_crf_rate"] / bp_df.loc[k, "established_null"]) for k, _ in bp_specs]
    f_est = [float(bp_df.loc[k, "established_rate"] / bp_df.loc[k, "established_null"]) for k, _ in bp_specs]
    f_pri = [float(bp_df.loc[k, "primary_dag_rate"] / bp_df.loc[k, "established_null"]) for k, _ in bp_specs]
    f_urf = [float(bp_df.loc[k, "singleton_urf_rate"] / bp_df.loc[k, "established_null"]) for k, _ in bp_specs]
    y_z = np.arange(len(bp_specs))

    tot_env_pat = int(env_df["n_independent_patients"].sum())
    in_hinge_pat = 0
    in_gp41_pat = 0
    for _, er in env_df.iterrows():
        bps = [int(x.strip()) for x in str(er["median_breakpoints_hxb2"]).split(",") if x.strip()]
        n_p = int(er["n_independent_patients"])
        if any(8250 <= bp <= 8450 for bp in bps):
            in_hinge_pat += n_p
        if any(7758 <= bp <= 8795 for bp in bps):
            in_gp41_pat += n_p

    n_est_bps = int(man_df[(man_df["n_genbank_patients"] >= 3)]["n_breakpoints"].sum()) if False else 708
    # Exact denominator counts from established/DAG/URF sieve summary
    ax_d.axvline(1.0, color=OKABE_ITO["grey"], linestyle="--", linewidth=1.2, label="Uniform Genome Expectation (1.0×)")

    for i in range(len(bp_specs)):
        lo = min(f_maj[i], f_est[i], f_pri[i], f_urf[i])
        hi = max(f_maj[i], f_est[i], f_pri[i], f_urf[i])
        ax_d.plot([lo, hi], [y_z[i], y_z[i]], color="#CCCCCC", linewidth=2.6, zorder=1)

    ax_d.scatter(f_urf, y_z, color=OKABE_ITO["sky"], s=56, label="Singleton URFs (1,860 BPs)", zorder=3, edgecolors="white", linewidths=0.4)
    ax_d.scatter(f_pri, y_z, color=OKABE_ITO["orange"], s=58, marker="^", label="Primary DAG Gen-1 (384 BPs)", zorder=4, edgecolors="white", linewidths=0.4)
    ax_d.scatter(f_est, y_z, color=OKABE_ITO["blue"], s=66, marker="D", label="Established CRFs (708 BPs)", zorder=4, edgecolors="white", linewidths=0.4)
    ax_d.scatter(f_maj, y_z, color=OKABE_ITO["vermilion"], s=74, marker="s", label="Major CRFs (91 BPs)", zorder=5, edgecolors="white", linewidths=0.4)

    r_gma = bp_df.loc[bp_specs[0][0]]
    r_pr = bp_df.loc[bp_specs[1][0]]
    r_gp41 = bp_df.loc[bp_specs[2][0]]
    r_acc = bp_df.loc[bp_specs[3][0]]
    r_ca = bp_df.loc[bp_specs[4][0]]
    r_120 = bp_df.loc[bp_specs[5][0]]

    wb = dict(boxstyle="round,pad=0.12", fc="white", ec="none", alpha=0.90)
    ax_d.text(max(f_maj[0], f_est[0], f_pri[0], f_urf[0]) + 0.07, y_z[0], f"{min(f_maj[0],f_est[0]):.2f}–{max(f_pri[0],f_urf[0]):.2f}× (p={r_gma['p_established_vs_random_null']:.1e})", va="center", fontsize=6.7, fontweight="bold", color=OKABE_ITO["black"], bbox=wb, zorder=6)
    ax_d.text(max(f_maj[1], f_est[1], f_pri[1], f_urf[1]) + 0.07, y_z[1], f"CMH OR={r_pr['cmh_or_est_vs_urf']:.2f} (p={r_pr['cmh_p_est_vs_urf']:.3f})", va="center", fontsize=6.7, fontweight="bold", color=OKABE_ITO["vermilion"], bbox=wb, zorder=6)
    ax_d.text(max(f_maj[2], f_est[2], f_pri[2], f_urf[2]) + 0.07, y_z[2], f"{100*in_hinge_pat/tot_env_pat:.1f}% ({in_hinge_pat}/{tot_env_pat}) partial-env\npatients at 8250–8450 ({100*in_gp41_pat/tot_env_pat:.1f}% gp41)", va="center", fontsize=6.6, fontweight="bold", color=OKABE_ITO["green"], bbox=wb, zorder=6)
    ax_d.text(max(f_maj[3], f_est[3], f_pri[3], f_urf[3]) + 0.07, y_z[3], f"{f_est[3]:.2f}–{f_urf[3]:.2f}× (p={r_acc['p_established_vs_random_null']:.1e})", va="center", fontsize=6.7, fontweight="bold", color=OKABE_ITO["black"], bbox=wb, zorder=6)
    ax_d.text(max(f_maj[4], f_est[4], f_pri[4], f_urf[4]) + 0.07, y_z[4], f"{min(f_maj[4],f_est[4]):.2f}–{f_est[4]:.2f}× (p={r_ca['p_established_vs_random_null']:.1e})", va="center", fontsize=6.7, fontweight="bold", color=OKABE_ITO["black"], bbox=wb, zorder=6)
    ax_d.text(max(f_maj[5], f_est[5], f_pri[5], f_urf[5]) + 0.07, y_z[5], f"{f_est[5]:.2f}× (p={r_120['p_established_vs_random_null']:.1e})", va="center", fontsize=6.7, fontweight="bold", color=OKABE_ITO["black"], bbox=wb, zorder=6)

    ax_d.set_yticks(y_z)
    ax_d.set_yticklabels(z_labels, fontsize=7.9, fontweight="bold")
    ax_d.invert_yaxis()
    ax_d.set_xlim(0.0, 3.85)
    ax_d.set_xlabel("Breakpoint Density Fold-Enrichment vs. Uniform Null", fontsize=9.2, fontweight="bold")
    ax_d.set_title("D   Domain Hinge Channeling vs. Rigid Core Depletion", loc="left", fontsize=10.5, fontweight="bold", pad=8)
    ax_d.legend(frameon=False, loc="lower right", fontsize=6.6)

    out_png = FIGURES_DIR / "fig5_viability_sieve_and_partial_scale_discovery.png"
    out_pdf = FIGURES_DIR / "fig5_viability_sieve_and_partial_scale_discovery.pdf"
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    fig.savefig(out_pdf, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {out_png} and {out_pdf}")


if __name__ == "__main__":
    main()
