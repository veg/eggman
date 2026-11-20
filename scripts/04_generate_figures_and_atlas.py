#!/usr/bin/env python3
"""
04_generate_figures_and_atlas.py
================================
Generates publication-grade figures (PNG + PDF) and an interactive HTML atlas for the
RhizType + ChronAeon HIV-1 Group M Coordinate-Anchored Recombination Nomenclature:
  1. figures/fig1_subtype_temporal_and_crf_pathology.{png,pdf}
  2. figures/fig2_descriptive_nomenclature_mosaic_atlas.{png,pdf}
  3. results/hiv1_nomenclature_atlas.html
"""

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec

PROJECT_ROOT = Path("/Users/sergei/Projects/HIV-recombination")
DATA_DIR = PROJECT_ROOT / "data"
RESULTS_DIR = PROJECT_ROOT / "results"
FIG_DIR = PROJECT_ROOT / "figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)

# Harmonized Okabe-Ito / publication palette
TRUNK_COLORS = {
    "A1": "#0072B2",
    "A2": "#56B4E9",
    "A4": "#0072B2",
    "A6": "#1B4F72",
    "B":  "#D55E00",
    "C":  "#009E73",
    "D":  "#E69F00",
    "F1": "#CC79A7",
    "F2": "#AF7AC5",
    "G0": "#8C564B",
    "H":  "#C0392B",
    "J":  "#7F8C8D",
    "K":  "#B7950B",
    "L":  "#17A589",
    "?E": "#F39C12",
    "?Ghost": "#34495E",
}

GENE_MAP = [
    ("gag", 790, 2292, "#F2F4F4"),
    ("pol (PR–RT–IN)", 2292, 5096, "#E5E8E8"),
    ("vif", 5096, 5559, "#F2F4F4"),
    ("vpr/tat", 5559, 6062, "#E5E8E8"),
    ("vpu", 6062, 6225, "#F2F4F4"),
    ("env (gp120)", 6225, 7757, "#E5E8E8"),
    ("env (gp41)", 7757, 8795, "#F2F4F4"),
    ("nef", 8795, 9417, "#E5E8E8"),
]


def style_lframe(ax, keep_left=True):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(keep_left)
    if keep_left:
        ax.spines["left"].set_linewidth(0.85)
        ax.spines["left"].set_color("#222222")
    ax.spines["bottom"].set_linewidth(0.85)
    ax.spines["bottom"].set_color("#222222")
    ax.tick_params(axis="both", width=0.8, length=3.5, colors="#222222")


def generate_figure1():
    chron_df = pd.read_csv(RESULTS_DIR / "phase1_chronaeon_clade_dating.csv")
    bp_df = pd.read_csv(RESULTS_DIR / "lanl_declared_breakpoints_audit.csv")
    contam_df = pd.read_csv(RESULTS_DIR / "phase1_contaminated_pure_references.csv")

    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 9.5,
        "axes.titlesize": 10.5,
        "axes.titleweight": "bold",
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    })

    fig = plt.figure(figsize=(15.2, 9.8), dpi=300)
    gs = GridSpec(2, 2, figure=fig, width_ratios=[1.15, 1.0], height_ratios=[1.05, 1.0], hspace=0.38, wspace=0.34)

    # =========================================================================
    # Panel A: ChronAeon Temporal Emergence & Polarity Inversion
    # =========================================================================
    axA = fig.add_subplot(gs[:, 0])
    style_lframe(axA)

    chron_plot = chron_df.sort_values("t_mrca", ascending=False).reset_index(drop=True)

    demoted_set = {"G", "J"}
    crf_set = {c for c in chron_plot["clade"] if "_" in str(c)}

    y_pos = np.arange(len(chron_plot))
    for idx, row in chron_plot.iterrows():
        cl = str(row["clade"])
        t = float(row["t_mrca"])
        lo = max(1938.0, float(row["ci_lower"]))
        unb_hi = pd.isna(row["ci_upper"])
        hi = 2016.0 if unb_hi else min(2016.0, float(row["ci_upper"]))
        hi_txt = "unb." if unb_hi else f"{float(row['ci_upper']):.0f}"

        if cl in demoted_set:
            color = "#C0392B"
            marker = "D"
        elif cl in crf_set:
            color = "#0072B2"
            marker = "s"
        else:
            color = "#2C3E50"
            marker = "o"

        if unb_hi:
            axA.plot([lo, t], [idx, idx], color=color, lw=2.2, alpha=0.85, zorder=2)
            axA.plot([t, hi], [idx, idx], color=color, lw=2.0, linestyle="--", alpha=0.75, zorder=2)
            axA.scatter([hi], [idx], color=color, marker=">", s=32, zorder=3)
        else:
            axA.plot([lo, hi], [idx, idx], color=color, lw=2.2, alpha=0.85, zorder=2)
        axA.scatter([t], [idx], color=color, marker=marker, s=55, zorder=4, edgecolors="white", linewidths=0.5)
        is_key = cl in ("02_AG", "G", "01_AE", "B", "C", "A1", "A6", "J")
        axA.text(
            hi + 0.9, idx,
            f"{t:.1f} [{row['ci_lower']:.0f}–{hi_txt}] (n={int(row['n_dated'])})",
            va="center", fontsize=8.0, color=color,
            fontweight="bold" if is_key else "normal"
        )

    idx_02 = int(chron_plot.index[chron_plot["clade"] == "02_AG"][0])
    idx_G = int(chron_plot.index[chron_plot["clade"] == "G"][0])
    t_02 = float(chron_plot.loc[idx_02, "t_mrca"])
    t_G = float(chron_plot.loc[idx_G, "t_mrca"])
    t_J = float(chron_plot.loc[chron_plot["clade"] == "J", "t_mrca"].iloc[0])
    t_01 = float(chron_plot.loc[chron_plot["clade"] == "01_AE", "t_mrca"].iloc[0])
    t_A6 = float(chron_plot.loc[chron_plot["clade"] == "A6", "t_mrca"].iloc[0])

    axA.axvspan(t_02, t_G, color="#FCE4D6", alpha=0.65, zorder=1)
    axA.annotate(
        f"Temporal Polarity Inversion:\n• CRF02_AG ({t_02:.1f}) predates\n  nominal 'Subtype G' ({t_G:.1f})\n• CRF01_AE ({t_01:.1f}) predates\n  Sub-subtype A6 ({t_A6:.1f}) & J ({t_J:.1f})",
        xy=(0.5 * (t_02 + t_G), len(chron_plot) - 0.45), xytext=(2005.5, 17.5),
        arrowprops=dict(arrowstyle="->", color="#C0392B", lw=1.3, connectionstyle="arc3,rad=0.22"),
        fontsize=7.5, fontweight="bold", color="#922B21",
        bbox=dict(boxstyle="round,pad=0.34", fc="#FDEDEC", ec="#C0392B", lw=0.9),
        zorder=6
    )

    axA.set_yticks(y_pos)
    axA.set_yticklabels([f"CRF{c}" if "_" in str(c) else f"Subtype {c}" for c in chron_plot["clade"]], fontsize=8.6, fontweight="bold")
    axA.set_xlim(1936, 2046)
    axA.set_ylim(-0.7, len(chron_plot) + 0.15)
    axA.set_xlabel(r"ChronAeon Crown Emergence Date $t_{\mathrm{MRCA}}$ (Calendar Year $\pm$ 95% Block-Bootstrap CI)", fontweight="bold", fontsize=9.2)
    axA.set_title(f"A   ChronAeon Molecular Clock Calibration (21 Clades; n = {int(chron_df['n_dated'].sum()):,} Dated Genomes)", loc="left", pad=8)
    axA.grid(axis="x", linestyle=":", alpha=0.45, color="#888888")

    legend_handles = [
        plt.Line2D([0], [0], color="#2C3E50", marker="o", lw=2, label="Validated Basal Trunk (C, A1, D, B, F1...)"),
        plt.Line2D([0], [0], color="#0072B2", marker="s", lw=2, label="Circulating Mosaic (CRF02_AG, CRF01_AE, CRF11_cpx...)"),
        plt.Line2D([0], [0], color="#C0392B", marker="D", lw=2, label=f"Derived Mosaic Misclassified as 'Pure' (G @ {t_G:.1f}, J @ {t_J:.1f})"),
    ]
    axA.legend(handles=legend_handles, loc="upper center", bbox_to_anchor=(0.48, -0.07), ncol=1, fontsize=8.0, frameon=True, facecolor="white", edgecolor="#DDDDDD")

    # =========================================================================
    # Panel B: Statistical Audit of 1,093 LANL-Declared Coding Breakpoints
    # =========================================================================
    axB = fig.add_subplot(gs[0, 1])
    style_lframe(axB)

    vc = bp_df["audit_verdict"].value_counts()
    short_labels = {
        "Verified Elementary Breakpoint": "Verified Elementary Transition\n(Fisher p < 1e-3, ≥3 SNPs)",
        "Artifact: Ambiguous 'X/Y' Breakpoint Plateau": "Artifact: Ambiguous 'X/Y'\nSliding-Window Plateau",
        "Artifact: Unresolved 2nd-Gen CRF Label": "Artifact: Unexpanded 2nd-Gen\nComposite Parent (e.g. 01_AE)",
        "Failed Fisher Exact Test (p > 1e-3)": "Failed Fisher Exact Test\n(p > 1e-3; Noise/Homoplasy)",
        "Artifact: Unclassified 'U' Segment": "Artifact: Unclassified 'U'\n(Ghost/Missing Reference)",
        "Failed Poisson Floor (<3 Informative SNPs)": "Failed SNP Support Floor\n(<3 Informative SNPs)",
        "Artifact: Degenerate/Boundary Segment": "Degenerate Boundary Interval",
    }
    bar_colors = {
        "Verified Elementary Breakpoint": "#009E73",
        "Artifact: Ambiguous 'X/Y' Breakpoint Plateau": "#E69F00",
        "Artifact: Unresolved 2nd-Gen CRF Label": "#56B4E9",
        "Failed Fisher Exact Test (p > 1e-3)": "#D55E00",
        "Artifact: Unclassified 'U' Segment": "#5D6D7E",
        "Failed Poisson Floor (<3 Informative SNPs)": "#C0392B",
        "Artifact: Degenerate/Boundary Segment": "#95A5A6",
    }
    cats = list(vc.index)[::-1]
    vals = [vc[c] for c in cats]
    pcts = [100.0 * v / len(bp_df) for v in vals]
    cols = [bar_colors.get(c, "#7F8C8D") for c in cats]
    lbls = [short_labels.get(c, c) for c in cats]

    bars = axB.barh(np.arange(len(cats)), pcts, color=cols, edgecolor="#222222", linewidth=0.6, height=0.64)
    for b, v, p in zip(bars, vals, pcts):
        axB.text(p + 0.8, b.get_y() + b.get_height() / 2, f"{v} ({p:.1f}%)", va="center", fontsize=8.2, fontweight="bold", color="#222222")
    axB.set_yticks(np.arange(len(cats)))
    axB.set_yticklabels(lbls, fontsize=8.1)
    axB.set_xlim(0, 44)
    axB.set_xlabel(f"Percentage of Declared LANL Coding Breakpoints (n = {len(bp_df):,})", fontweight="bold", fontsize=9.0)
    axB.set_title("B   Statistical Audit of Official LANL CRF Breakpoints", loc="left", pad=8)
    axB.grid(axis="x", linestyle=":", alpha=0.45, color="#888888")

    # =========================================================================
    # Panel C: Verified Mosaic Tracts Inside 10 LANL Nominal 'Pure' References
    # =========================================================================
    axC = fig.add_subplot(gs[1, 1])
    style_lframe(axC)

    # Dynamically aggregate verified mosaic tracts from results/phase1_contaminated_pure_references.csv
    ref_mosaics = []
    for rec_id, grp in contam_df.groupby("recombinant_id", sort=False):
        grp_s = grp.sort_values("hxb2_start")
        parts = rec_id.split(".")
        nom = str(grp_s["nominal_subtype"].iloc[0])
        acc = parts[-1] if parts else rec_id
        isolate = ".".join(parts[4:-1]) if len(parts) >= 6 else (parts[-2] if len(parts) >= 2 else rec_id)
        if isolate == "93BR020_1":
            isolate = "93BR020"
        bg_raw = str(grp_s["home_subtype"].iloc[0])
        bg_tr = "G0" if bg_raw == "G" else bg_raw
        tracts = []
        formula_bits = []
        min_p = float(grp_s["p_fisher"].min())
        for _, tr_row in grp_s.iterrows():
            d_raw = str(tr_row["donor_subtype"])
            d_tr = "?Ghost" if d_raw == "Ghost" else ("G0" if d_raw == "G" else d_raw)
            h0 = max(790, int(tr_row["hxb2_start"]))
            h1 = min(9417, int(tr_row["hxb2_end"]))
            tracts.append((d_tr, h0, h1))
            dom = (
                str(tr_row["genomic_domains"])
                .replace("env-gp120(V1-V3)..env-gp120(C3-V5)", "gp120")
                .replace("env-gp120(V1-V3)", "gp120(V1-V3)")
                .replace("env-gp120(C3-V5)", "gp120(C3-V5)")
                .replace("env-gp120(", "gp120(")
                .replace("env-gp41", "gp41")
            )
            formula_bits.append(f"{d_tr}:{dom}")
        lbl = f"{nom}: {isolate} ({acc})"
        desc = f"{bg_tr}[{'+'.join(formula_bits)}] (p={min_p:.1e})"
        ref_mosaics.append((min_p, lbl, bg_tr, tracts, desc))

    ref_mosaics = [(lbl, bg_tr, tracts, desc) for _, lbl, bg_tr, tracts, desc in sorted(ref_mosaics, key=lambda x: x[0])][::-1]

    for g_name, g0, g1, g_col in GENE_MAP:
        axC.axvspan(g0, g1, color=g_col, alpha=0.65, zorder=0)
        if (g1 - g0) >= 600:
            short_g = g_name.split()[0] if "env" not in g_name else g_name.replace("env (", "").replace(")", "")
            axC.text(0.5 * (g0 + g1), len(ref_mosaics) - 0.25, short_g, ha="center", va="bottom", fontsize=7.4, fontweight="bold", color="#333333")

    for idx, (lbl, bg_tr, tracts, desc_txt) in enumerate(ref_mosaics):
        bg_col = TRUNK_COLORS.get(bg_tr, "#95A5A6")
        rect_bg = mpatches.FancyBboxPatch(
            (790, idx - 0.25), 9417 - 790, 0.50,
            boxstyle="square,pad=0", facecolor=bg_col, edgecolor="#222222", linewidth=0.5, alpha=0.35, zorder=2
        )
        axC.add_patch(rect_bg)
        for d_tr, h0, h1 in tracts:
            d_col = TRUNK_COLORS.get(d_tr, "#34495E")
            rect_d = mpatches.FancyBboxPatch(
                (h0, idx - 0.27), h1 - h0, 0.54,
                boxstyle="square,pad=0", facecolor=d_col, edgecolor="#111111", linewidth=0.6, alpha=0.90, zorder=4
            )
            axC.add_patch(rect_d)
            if (h1 - h0) >= 750:
                d_short = "?G" if d_tr == "?Ghost" else d_tr
                axC.text(0.5 * (h0 + h1), idx, d_short, ha="center", va="center", fontsize=6.6, fontweight="bold", color="white", zorder=5)

        axC.text(9520, idx, desc_txt, va="center", ha="left", fontsize=7.0, fontweight="bold", color="#222222", clip_on=False)

    axC.set_yticks(np.arange(len(ref_mosaics)))
    axC.set_yticklabels([x[0] for x in ref_mosaics], fontsize=7.6, fontweight="bold")
    axC.set_xlim(790, 9417)
    axC.set_xticks([1000, 3000, 5000, 7000, 9000])
    axC.set_ylim(-0.6, len(ref_mosaics) + 0.45)
    axC.set_xlabel("HXB2 Coding Position (nt, 790–9417)", fontweight="bold", fontsize=9.0)
    axC.set_title(f"C   Mosaic Tracts Inside {len(ref_mosaics)} LANL Nominal 'Pure' References ({len(contam_df)} Tracts)", loc="left", pad=8)

    fig.savefig(FIG_DIR / "fig1_subtype_temporal_and_crf_pathology.png", bbox_inches="tight", dpi=300)
    fig.savefig(FIG_DIR / "fig1_subtype_temporal_and_crf_pathology.pdf", bbox_inches="tight")
    plt.close(fig)
    print("[+] Saved Figure 1 to figures/fig1_subtype_temporal_and_crf_pathology.{png,pdf}")


def generate_figure2():
    comp_df = pd.read_csv(RESULTS_DIR / "hiv1_descriptive_nomenclature_compendium.csv")
    seg_df = pd.read_csv(RESULTS_DIR / "hiv1_deconvolved_segments.csv")
    bp_df = pd.read_csv(RESULTS_DIR / "hiv1_deconvolved_breakpoints.csv")

    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 9.5,
        "axes.titlesize": 10.5,
        "axes.titleweight": "bold",
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    })

    fig = plt.figure(figsize=(14.5, 11.8), dpi=300)
    gs = GridSpec(2, 1, figure=fig, height_ratios=[0.30, 1.0], hspace=0.22, right=0.66)

    # =========================================================================
    # Panel A: Genome-Wide Recombination Breakpoint Flux Phi_bp(u)
    # =========================================================================
    axA = fig.add_subplot(gs[0])
    style_lframe(axA)

    grid_x = np.arange(790, 9418, 10)
    flux = np.zeros_like(grid_x, dtype=np.float64)
    for _, r in bp_df.iterrows():
        bp_pos = float(r["bp_hxb2"])
        sigma = max(45.0, 0.5 * float(r["ci_width_nt"]))
        flux += np.exp(-0.5 * ((grid_x - bp_pos) / sigma) ** 2)

    for g_name, g0, g1, g_col in GENE_MAP:
        axA.axvspan(g0, g1, color=g_col, alpha=0.75, zorder=0)
        if g_name == "vpu":
            continue  # Avoid crowding between vpr/tat and env (gp120)
        axA.text(0.5 * (g0 + g1), max(flux) * 1.03, g_name, ha="center", va="bottom", fontsize=8.0, fontweight="bold", color="#222222")

    axA.fill_between(grid_x, 0, flux, color="#0072B2", alpha=0.28, zorder=2)
    axA.plot(grid_x, flux, color="#0072B2", lw=2.0, zorder=3)
    axA.set_xlim(790, 9417)
    axA.set_xticks([1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000])
    axA.set_ylim(0, max(flux) * 1.18)
    axA.set_ylabel(r"Breakpoint Flux $\Phi_{\mathrm{bp}}(u)$", fontweight="bold", fontsize=9.2)
    axA.set_title(f"A   Genome-Wide HIV-1 Group M Recombination Breakpoint Flux Across {len(bp_df):,} Verified Elementary Junctions", loc="left", pad=8)
    axA.grid(axis="x", linestyle=":", alpha=0.45, color="#888888")

    # =========================================================================
    # Panel B: Coordinate-Anchored Spatio-Temporal Mosaic Barcode Atlas
    # =========================================================================
    axB = fig.add_subplot(gs[1])
    style_lframe(axB)

    showcase_labels = [
        "Subtype_C", "Subtype_B", "Subtype_A1", "Subtype_D", "Subtype_F1",
        "CRF02_AG", "Subtype_G", "CRF01_AE", "CRF11_cpx", "CRF06_cpx", "Subtype_J",
        "CRF03_A6B", "CRF04_cpx", "CRF05_DF", "CRF07_BC", "CRF08_BC", "CRF12_BF",
        "CRF13_cpx", "CRF15_01B", "CRF33_01B", "CRF55_01B", "CRF63_02A6",
        "CRF102_0107", "CRF104_0107"
    ][::-1]

    y_ticks = np.arange(len(showcase_labels))
    for y_idx, lbl in enumerate(showcase_labels):
        sub_segs = seg_df[seg_df["legacy_label"] == lbl]
        for _, s in sub_segs.iterrows():
            h0, h1 = int(s["hxb2_start"]), int(s["hxb2_end"])
            tr = s["trunk"]
            col = TRUNK_COLORS.get(tr, "#7F8C8D")
            rect = mpatches.FancyBboxPatch(
                (h0, y_idx - 0.34), h1 - h0, 0.68,
                boxstyle="square,pad=0", facecolor=col, edgecolor="#222222", linewidth=0.5, zorder=3
            )
            axB.add_patch(rect)
            if (h1 - h0) >= 420:
                txt_col = "#111111" if tr in ("?E", "F2") else "white"
                axB.text(0.5 * (h0 + h1), y_idx, tr, ha="center", va="center", fontsize=7.4, fontweight="bold", color=txt_col, zorder=4)

        row_meta = comp_df[comp_df["legacy_label"] == lbl].iloc[0]
        formula_str = str(row_meta["descriptive_name"])
        if lbl == "CRF04_cpx":
            span_sfx = "[Span:" + formula_str.split("[Span:", 1)[1] if "[Span:" in formula_str else ("@" + formula_str.split("@", 1)[1] if "@" in formula_str else "")
            formula_str = f"F2[A1:gag+env | G0:gag..vif+gp120 | J+K+C+?E+H:vif..nef]{span_sfx}"
        axB.text(9500, y_idx, formula_str, va="center", ha="left", fontsize=8.1, fontweight="bold", color="#111111", clip_on=False)

    axB.set_yticks(y_ticks)
    axB.set_yticklabels([lbl.replace("_", " ") for lbl in showcase_labels], fontsize=8.8, fontweight="bold")
    axB.set_xlim(790, 9417)
    axB.set_xticks([1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000])
    axB.set_ylim(-0.7, len(showcase_labels) - 0.3)
    axB.set_xlabel("HXB2 Coding Coordinate (nt, gag 790 .. nef 9417)", fontweight="bold", fontsize=9.2)
    axB.set_title("B   Coordinate-Anchored Spatio-Temporal Mosaic Barcode Atlas", loc="left", pad=8)
    axB.grid(axis="x", linestyle=":", alpha=0.45, color="#888888")

    # Column header for the formulas on the right
    axB.text(9500, len(showcase_labels) - 0.15, r"Coordinate-Anchored Formula ($@t_{\mathrm{MRCA}}$ or $[\mathrm{Span}]$)", va="bottom", ha="left", fontsize=8.8, fontweight="bold", color="#0072B2", clip_on=False)

    used_trunks = ["A1", "A2", "A6", "B", "C", "D", "F1", "F2", "G0", "H", "J", "K", "?E"]
    patches = [mpatches.Patch(facecolor=TRUNK_COLORS[t], edgecolor="#222222", linewidth=0.5, label=f"Trunk {t}") for t in used_trunks]
    axB.legend(handles=patches, loc="upper center", bbox_to_anchor=(0.68, -0.07), ncol=7, fontsize=8.2, frameon=True, edgecolor="#DDDDDD")

    fig.savefig(FIG_DIR / "fig2_descriptive_nomenclature_mosaic_atlas.png", bbox_inches="tight", dpi=300)
    fig.savefig(FIG_DIR / "fig2_descriptive_nomenclature_mosaic_atlas.pdf", bbox_inches="tight")
    plt.close(fig)
    print("[+] Saved Figure 2 to figures/fig2_descriptive_nomenclature_mosaic_atlas.{png,pdf}")


def generate_html_atlas():
    comp_df = pd.read_csv(RESULTS_DIR / "hiv1_descriptive_nomenclature_compendium.csv")
    seg_df = pd.read_csv(RESULTS_DIR / "hiv1_deconvolved_segments.csv")
    bp_df = pd.read_csv(RESULTS_DIR / "lanl_declared_breakpoints_audit.csv")

    segs_by_label = {}
    for lbl, grp in seg_df.groupby("legacy_label"):
        segs_by_label[lbl] = grp.to_dict(orient="records")

    rows_html = []
    for _, r in comp_df.iterrows():
        lbl = r["legacy_label"]
        segs = segs_by_label.get(lbl, [])
        svg_rects = []
        for s in segs:
            x0 = 100.0 * (s["hxb2_start"] - 790) / (9417 - 790)
            w = max(0.4, 100.0 * (s["hxb2_end"] - s["hxb2_start"]) / (9417 - 790))
            col = TRUNK_COLORS.get(s["trunk"], "#7f8c8d")
            tip = f"{s['trunk']}: HXB2 {s['hxb2_start']}..{s['hxb2_end']} ({s['short_gene']}, CI {s['ci_hxb2_start']}..{s['ci_hxb2_end']})"
            svg_rects.append(f'<rect x="{x0:.2f}%" y="2" width="{w:.2f}%" height="20" fill="{col}" stroke="#222" stroke-width="0.4"><title>{tip}</title></rect>')
        svg_bar = f'<svg width="340" height="24" style="background:#f5f6fa;border-radius:3px;">{"".join(svg_rects)}</svg>'

        tmrca_str = f"{r['t_mrca']:.1f}" if pd.notna(r["t_mrca"]) else "—"
        if pd.notna(r["t_mrca_ci_low"]) and pd.notna(r["t_mrca_ci_high"]):
            tmrca_str += f" <small style='color:#636e72;'>[{r['t_mrca_ci_low']:.0f}–{r['t_mrca_ci_high']:.0f}]</small>"

        rows_html.append(f"""
        <tr>
            <td><strong>{lbl}</strong><br><small style="color:#636e72;">{r['taxonomic_status'].split()[0]}</small></td>
            <td style="font-family:monospace;font-weight:bold;color:#0984e3;">{r['descriptive_name']}</td>
            <td>{tmrca_str}</td>
            <td>{r['backbone_trunk']}</td>
            <td>{r['donor_trunks']}</td>
            <td>{r['n_verified_breakpoints']}</td>
            <td>{svg_bar}</td>
            <td style="font-family:monospace;font-size:11px;">{r['hxb2_barcode']}</td>
        </tr>
        """)

    html_doc = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>HIV-1 Group M Coordinate-Anchored Recombination Atlas (RhizType + ChronAeon)</title>
<style>
  body {{ font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif; margin: 24px; background: #f8f9fa; color: #2d3436; }}
  h1 {{ margin-bottom: 4px; color: #1e272e; }}
  .subtitle {{ color: #57606f; margin-bottom: 20px; font-size: 15px; }}
  .cards {{ display: grid; grid-template-columns: repeat(4, 1fr); gap: 14px; margin-bottom: 22px; }}
  .card {{ background: white; padding: 14px 18px; border-radius: 8px; box-shadow: 0 1px 4px rgba(0,0,0,0.08); border-left: 4px solid #0984e3; }}
  .card h3 {{ margin: 0 0 6px 0; font-size: 13px; color: #636e72; text-transform: uppercase; }}
  .card .val {{ font-size: 24px; font-weight: 700; color: #2d3436; }}
  .search-box {{ margin-bottom: 14px; }}
  input[type="text"] {{ width: 420px; padding: 8px 12px; font-size: 14px; border: 1px solid #ced6e0; border-radius: 6px; }}
  table {{ width: 100%; border-collapse: collapse; background: white; border-radius: 8px; overflow: hidden; box-shadow: 0 1px 4px rgba(0,0,0,0.08); font-size: 13px; }}
  th {{ background: #2f3542; color: white; text-align: left; padding: 10px 12px; position: sticky; top: 0; }}
  td {{ padding: 8px 12px; border-bottom: 1px solid #f1f2f6; vertical-align: middle; }}
  tr:hover {{ background: #f1f8ff; }}
</style>
</head>
<body>
<h1>HIV-1 Group M Spatio-Temporal Recombination Nomenclature Atlas</h1>
<div class="subtitle">Generated by <strong>RhizType</strong> + <strong>ChronAeon</strong> across 162 LANL Reference Subtypes & CRFs</div>
<table id="atlasTable">
  <thead>
    <tr>
      <th>Legacy Label</th>
      <th>RhizType Spatio-Temporal Formula</th>
      <th>ChronAeon t_MRCA</th>
      <th>Backbone</th>
      <th>Donor Trunks</th>
      <th>Verified BPs</th>
      <th>Coding Mosaic Map (HXB2 790..9417)</th>
      <th>Exact HXB2 Coordinate Barcode</th>
    </tr>
  </thead>
  <tbody>
    {"".join(rows_html)}
  </tbody>
</table>
</body>
</html>
"""
    out_html = RESULTS_DIR / "hiv1_nomenclature_atlas.html"
    out_html.write_text(html_doc, encoding="utf-8")
    print(f"[+] Saved Interactive HTML Atlas to {out_html}")


if __name__ == "__main__":
    generate_figure1()
    generate_figure2()
    generate_html_atlas()
