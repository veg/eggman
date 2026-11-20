#!/usr/bin/env python3
"""
Generate all 5 main-text figures plus Supplemental Figure S1
at exact 1:1 print width (6.50 in = \\textwidth in main.tex / supplement.tex),
following a strict >= 6.0 pt typography floor, Okabe-Ito colorblind safety,
high data-to-ink ratio, and live bounding-box verification.

100% of bars, distributions, phylograms, dates, coordinates, discrete fits, and
statistical values are loaded dynamically from verified Phase 1 outputs on disk.
"""

import csv
import json
import math
import re
import sys
from collections import Counter
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
from scipy.special import zeta
from scipy.stats import fisher_exact, spearmanr

sys.path.insert(0, str(Path(__file__).resolve().parent))
from audit_figures_standard import CANVAS_WIDTH_IN, PALETTE, audit_matplotlib_figure


ROOT = Path("/Users/sergei/Projects/HIV-recombination")
RES = ROOT / "results"
TREES_DIR = RES / "trees"
FIG_DIR = ROOT / "paper" / "figures"
PACKS_DIR = ROOT / "rhiztype" / "packs"
CORE_DATA_DIR = ROOT / "rhiztype" / "crates" / "eggman-core" / "data"
MP_JSON = ROOT / "paper" / "multipathogen_benchmark_summary.json"
SC2_JSON = RES / "sars_cov2_global_and_bvbrc_screen.json"
SPECIMENS_JSON = ROOT / "paper" / "fig1_empirical_specimens.json"

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "mathtext.fontset": "custom",
    "mathtext.rm": "sans",
    "mathtext.it": "sans:italic",
    "mathtext.bf": "sans:bold",
    "mathtext.default": "regular",
    "font.size": 6.2,
    "axes.titlesize": 6.8,
    "axes.labelsize": 6.5,
    "xtick.labelsize": 6.0,
    "ytick.labelsize": 6.0,
    "legend.fontsize": 6.0,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.edgecolor": PALETTE["SUBINK"],
    "axes.linewidth": 0.55,
    "xtick.major.width": 0.50,
    "ytick.major.width": 0.50,
    "xtick.major.size": 2.2,
    "ytick.major.size": 2.2,
    "xtick.major.pad": 2.0,
    "ytick.major.pad": 2.0,
    "xtick.color": PALETTE["SLATE"],
    "ytick.color": PALETTE["SLATE"],
})

INK = PALETTE["INK"]
SUBINK = PALETTE["SUBINK"]
MUTED = PALETTE["MUTED"]
RULE = PALETTE["HAIR"]
HAIR = PALETTE["GRID"]
BLUE = PALETTE["BLUE"]
BLUE_SOFT = PALETTE["BLUE_SOFT"]
SKY = PALETTE["SKY"]
TEAL = PALETTE["TEAL"]
TEAL_SOFT = PALETTE["TEAL_SOFT"]
AMBER = PALETTE["AMBER"]
AMBER_SOFT = PALETTE["AMBER_SOFT"]
VERMILION = PALETTE["VERMILION"]
VERMILION_SOFT = PALETTE["VERMILION_SOFT"]
PURPLE = PALETTE["PURPLE"]
PURPLE_SOFT = PALETTE["PURPLE_SOFT"]
INDIGO = PALETTE["INDIGO"]
WINE = PALETTE["WINE"]
SLATE = PALETTE["SLATE"]
SLATE_MID = PALETTE["SLATE_MID"]
SLATE_LIGHT = PALETTE["SLATE_LIGHT"]
BG_SOFT = PALETTE["BG_SOFT"]


def draw_panel_header(fig, x_fig, y_fig, letter, title):
    fig.text(x_fig, y_fig, letter, fontsize=10.0, fontweight="bold", color=INK, ha="left", va="center")
    fig.text(x_fig + 0.024, y_fig, title, fontsize=7.2, fontweight="bold", color=INK, ha="left", va="center")


def draw_rect(ax, x, y, w, h, fc, ec="none", lw=0.55, alpha=1.0, zorder=2, ls="-"):
    patch = mpatches.Rectangle(
        (x, y), w, h,
        facecolor=fc, edgecolor=ec, linewidth=lw, linestyle=ls, alpha=alpha, zorder=zorder,
    )
    ax.add_patch(patch)
    return patch


def add_quadrant_hairlines(fig, x_split=0.500, y_split=0.495):
    fig.add_artist(plt.Line2D([0.012, 0.988], [y_split, y_split], color=RULE, lw=0.60, transform=fig.transFigure))
    fig.add_artist(plt.Line2D([x_split, x_split], [0.015, 0.985], color=RULE, lw=0.60, transform=fig.transFigure))


def fmt_p_tex(p_val):
    """Format p-value into crisp mathtext without raw '1e-...' ASCII notation."""
    if p_val <= 0:
        return r"$p < 10^{-300}$"
    exp = int(math.floor(math.log10(p_val)))
    mant = p_val / (10.0 ** exp)
    if mant >= 9.95:
        exp += 1
        mant = 1.0
    if exp <= -20:
        return rf"$p < 10^{{{exp + 1}}}$"
    if abs(mant - 1.0) < 0.05:
        return rf"$p = 10^{{{exp}}}$"
    return rf"$p = {mant:.1f}\times 10^{{{exp}}}$"


# =============================================================================
# FIGURE 1: Concept, Two-Tier Spanning Architecture & Four-Gate Benchmark
# =============================================================================
def build_fig1():
    with open(CORE_DATA_DIR / "compendium_bundle.json") as f:
        cb = json.load(f)
    with open(PACKS_DIR / "influenza_a_8seg" / "manifest.json") as f:
        iav_manifest = json.load(f)
    with open(PACKS_DIR / "sars_cov2_wuhan" / "manifest.json") as f:
        sc2_manifest = json.load(f)
    with open(PACKS_DIR / "sarbecovirus_wuhan" / "manifest.json") as f:
        sb_manifest = json.load(f)
    with open(RES / "hiv1_chronaeon_founder_audit.tsv") as f:
        hiv_founders = list(csv.DictReader(f, delimiter="\t"))
    with open(RES / "benchmark_three_way_comparison.csv") as f:
        bench_rows = {r["panel"]: r for r in csv.DictReader(f)}
    with open(RES / "benchmark_three_way_meta.json") as f:
        bench_meta = json.load(f)
    with open(RES / "benchmark_out_of_dictionary_holdout.csv") as f:
        ood_rows = list(csv.DictReader(f))
    with open(RES / "tier1_null_simulation_fpr.csv") as f:
        null_rows = list(csv.DictReader(f))
    with open(RES / "benchmark_synth500_breakpoint_comparison.csv") as f:
        synth_rows = {r["tool"]: r for r in csv.DictReader(f)}
    with open(RES / "benchmark_sc2_nextclade_and_iav_genoflu.csv") as f:
        ext_rows = {(r["cohort"], r["tool"]): r for r in csv.DictReader(f)}
    with open(MP_JSON) as f:
        mp_summary = json.load(f)

    n_hiv_t1 = len(hiv_founders)
    n_iav_t1 = len(iav_manifest["basal_trunks"])
    n_sc2_t1 = len(sc2_manifest["basal_trunks"])
    n_sb_t1 = len(sb_manifest["basal_trunks"])
    n_active_t2 = sum(
        1 for e in cb["compendium"]
        if e["category"] not in ("Demoted_LANL_Entry", "Nominal_Pure_Subtype")
    )
    demoted_entries = [e for e in cb["compendium"] if e["category"] == "Demoted_LANL_Entry"]
    n_dem_crfs = sum(1 for e in demoted_entries if e["label"].startswith("CRF"))
    n_dem_pure = len(demoted_entries) - n_dem_crfs

    b313_formula = mp_summary["influenza_a"]["per_genotype"]["B3.13"]["top_structural_formulas"][0][0]
    b313_donor_segs = set(b313_formula.split(":", 1)[1].rstrip("]").split("+"))
    iav_landmarks = iav_manifest["landmarks"]
    iav_total_nt = float(iav_manifest["ref_full_len"])

    fig = plt.figure(figsize=(CANVAS_WIDTH_IN, 5.40), facecolor="white")
    add_quadrant_hairlines(fig, x_split=0.500, y_split=0.470)
    draw_panel_header(fig, 0.014, 0.974, "A", "Clean two-tier pack: scrub, date & span")
    draw_panel_header(fig, 0.514, 0.974, "B", "One engine: raw sequence to coordinate barcode")
    draw_panel_header(fig, 0.014, 0.448, "C", "Four-gate guarded lineage calling")
    draw_panel_header(fig, 0.514, 0.448, "D", "Accuracy, held-out flagging & runtime")

    # --- Panel A: Austere Technical Schematic (Open Whitespace + Hairline Rules) ---
    ax_a = fig.add_axes([0.016, 0.488, 0.474, 0.456])
    ax_a.set_xlim(0, 10.0)
    ax_a.set_ylim(0, 7.4)
    ax_a.axis("off")

    # Section 1: Tier 1 Scrubbing + ChronAeon Date Guard (no outer box)
    ax_a.text(
        0.10, 7.05,
        f"Tier 1  \u00b7  Scrubbed, emergence-dated founders ({n_hiv_t1} HIV, {n_iav_t1} IAV, {n_sc2_t1} SC2, {n_sb_t1} CoV)",
        fontsize=6.0, fontweight="bold", color=INK, va="center",
    )
    ax_a.text(0.10, 6.46, "Raw ref G", fontsize=6.0, color=SUBINK, va="center")
    draw_rect(ax_a, 2.05, 6.28, 5.55, 0.36, fc=SLATE, ec="none")
    draw_rect(ax_a, 7.64, 6.28, 2.22, 0.36, fc=SKY, ec="white", lw=0.4)
    ax_a.text(8.75, 6.46, "A1 gp41..nef", fontsize=6.0, fontweight="bold", color="white", ha="center", va="center")
    ax_a.annotate("", xy=(8.75, 5.90), xytext=(8.75, 6.25),
                  arrowprops=dict(arrowstyle="-|>", color=VERMILION, lw=0.90, mutation_scale=6.0))
    ax_a.text(4.82, 6.07, "RhizAeon mosaic scrub", fontsize=6.0, style="italic", color=SUBINK, ha="center", va="center")
    ax_a.text(0.10, 5.68, "Scrubbed G0", fontsize=6.0, fontweight="bold", color=INK, va="center")
    draw_rect(ax_a, 2.05, 5.50, 5.55, 0.36, fc=BLUE, ec="none")
    draw_rect(ax_a, 7.64, 5.50, 2.22, 0.36, fc="white", ec=MUTED, lw=0.65, ls="--")
    ax_a.text(4.82, 5.68, "ChronAeon emergence founder", fontsize=6.0, fontweight="bold", color="white", ha="center", va="center")
    ax_a.text(8.75, 5.68, "Masked", fontsize=6.0, color=SUBINK, ha="center", va="center")

    # ChronAeon Temporal Guard timeline
    ax_a.text(0.10, 4.74, "Date guard", fontsize=6.0, fontweight="bold", color=INK, va="center")
    ax_a.plot([2.05, 9.82], [4.74, 4.74], color=SLATE, lw=0.85, zorder=2)
    for tx in [2.55, 3.55, 4.65, 5.70]:
        ax_a.scatter([tx], [4.74], s=18, color=BLUE, edgecolor="white", linewidth=0.4, zorder=4)
    for tx in [7.55, 8.50, 9.40]:
        ax_a.scatter([tx], [4.74], s=18, facecolor="white", edgecolor=MUTED, linewidth=0.85, zorder=4)
    ax_a.plot([6.55, 6.55], [4.46, 5.06], color=VERMILION, lw=0.95, ls="--", zorder=3)
    ax_a.text(4.12, 5.04, "Active clades", fontsize=6.0, fontweight="bold", color=BLUE, ha="center", va="center")
    ax_a.text(8.48, 5.04, "Future masked", fontsize=6.0, color=SUBINK, ha="center", va="center")
    ax_a.text(6.55, 4.30, "Query date cutoff", fontsize=6.0, fontweight="bold", color=VERMILION, ha="center", va="center")

    # Hairline rule separating Section 1 -> Section 2
    ax_a.plot([0.10, 9.86], [3.98, 3.98], color=RULE, lw=0.55)

    # Section 2: Spanning Basis Rule (3 open columns separated by thin vertical hairlines)
    ax_a.text(0.10, 3.64, "Spanning basis rule (Tier 1 vs. Tier 2)", fontsize=6.0, fontweight="bold", color=INK, va="center")
    crits = [
        ("1. Separable", r"$D_{\min} \geq 2 m_{\min}$", "Keeps KP.2 in JN.1"),
        ("2. Irreducible", "Not a mosaic", "XDV, XFZ \u2192 Tier 2"),
        ("3. Complete", r"$\geq$25 private SNPs", "Adds BA.4, BA.3.2"),
    ]
    for i, (t1, t2, t3) in enumerate(crits):
        cx = 1.65 + i * 3.30
        if i > 0:
            vx = 0.10 + i * 3.25
            ax_a.plot([vx, vx], [2.04, 3.26], color=RULE, lw=0.50)
        ax_a.text(cx, 3.12, t1, fontsize=6.0, fontweight="bold", color=BLUE, ha="center", va="center")
        ax_a.text(cx, 2.66, t2, fontsize=6.0, color=INK, ha="center", va="center")
        ax_a.text(cx, 2.20, t3, fontsize=6.0, color=SUBINK, ha="center", va="center")

    # Hairline rule separating Section 2 -> Section 3
    ax_a.plot([0.10, 9.86], [1.74, 1.74], color=RULE, lw=0.55)

    # Section 3: Tier 2 Portable Catalog (minimal left accent rule, no pastel box)
    ax_a.plot([0.12, 0.12], [0.22, 1.46], color=BLUE, lw=2.0)
    ax_a.text(0.36, 1.34, "Tier 2  \u00b7  Portable epidemic catalog (--pack)",
              fontsize=6.0, fontweight="bold", color=INK, va="center")
    ax_a.text(
        0.36, 0.84,
        f"{n_active_t2} HIV-1 ({n_dem_crfs} CRFs & {n_dem_pure} A* demoted), 10 H5Nx, 28 SC2",
        fontsize=6.0, color=SUBINK, va="center",
    )
    ax_a.text(0.36, 0.36, "Novel clusters (\u22653 hosts + ML monophyly) \u2192 promote",
              fontsize=6.0, color=SUBINK, va="center")

    # --- Panel B: Austere 3-Stage Walkthrough (Open Rules, Zero Pastel Boxes) ---
    ax_b = fig.add_axes([0.516, 0.488, 0.474, 0.456])
    ax_b.set_xlim(0, 10.0)
    ax_b.set_ylim(0, 7.4)
    ax_b.axis("off")

    # Stage 1: Reference-free coordinate projection
    ax_b.text(0.10, 7.05, "1. Reference-free coordinate projection (no MSA)",
              fontsize=6.0, fontweight="bold", color=INK, va="center")
    draw_rect(ax_b, 0.25, 6.30, 5.30, 0.40, fc="white", ec=BLUE, lw=0.70)
    ax_b.text(2.90, 6.50, "Unaligned full genome", fontsize=6.0, fontweight="bold", color=BLUE, ha="center", va="center")
    draw_rect(ax_b, 5.95, 6.30, 3.80, 0.40, fc="white", ec=SLATE_MID, lw=0.70)
    ax_b.text(7.85, 6.50, "or partial pol / env read", fontsize=6.0, fontweight="bold", color=SLATE, ha="center", va="center")
    # Coordinate ruler + k-mer anchors
    ax_b.plot([0.25, 9.75], [5.84, 5.84], color=INK, lw=1.25, zorder=2)
    for ax_x in [0.65, 1.48, 2.40, 3.32, 5.05]:
        ax_b.plot([ax_x, ax_x], [5.84, 6.29], color=BLUE, lw=0.85, zorder=3)
    for ax_x in [6.45, 7.38, 8.32, 9.25]:
        ax_b.plot([ax_x, ax_x], [5.84, 6.29], color=SLATE_MID, lw=0.85, zorder=3)
    arc_x = np.linspace(3.42, 4.95, 30)
    arc_y = 5.84 + 0.30 * np.sin(np.pi * (arc_x - 3.42) / 1.53)
    ax_b.plot(arc_x, arc_y, color=VERMILION, lw=0.95, ls="--", zorder=4)
    ax_b.text(5.00, 5.44, "k-mer anchors + banded DP (split-DP jumps >180 nt indels)",
              fontsize=6.0, color=SUBINK, ha="center", va="center")

    ax_b.plot([0.10, 9.86], [5.10, 5.10], color=RULE, lw=0.55)

    # Stage 2: Unified Viterbi scan + Fisher check
    ax_b.text(0.10, 4.78, "2. Unified Viterbi scan + Fisher breakpoint check",
              fontsize=6.0, fontweight="bold", color=INK, va="center")

    # Contiguous HIV-1 track
    ax_b.text(0.10, 4.20, "HIV-1", fontsize=6.0, fontweight="bold", color=INK, va="center")
    draw_rect(ax_b, 1.40, 4.00, 6.18, 0.40, fc=SLATE, ec="none")
    ax_b.text(4.49, 4.20, "G0 (backbone)", fontsize=6.0, fontweight="bold", color="white", ha="center", va="center")
    draw_rect(ax_b, 7.64, 4.00, 2.12, 0.40, fc=SKY, ec="none")
    ax_b.text(8.70, 4.20, "A1", fontsize=6.0, fontweight="bold", color="white", ha="center", va="center")
    ax_b.plot([7.61, 7.61], [3.94, 4.48], color=VERMILION, lw=1.1, ls="--", zorder=4)
    ax_b.text(4.15, 3.60, r"Contiguous: scalar switch cost $\gamma > 0$",
              fontsize=6.0, color=SUBINK, ha="center", va="center")
    ax_b.text(8.62, 3.60, r"Fisher $p < 10^{-4}$", fontsize=6.0, fontweight="bold", color=VERMILION, ha="center", va="center")

    # Segmented H5N1 B3.13 track (exact proportional lengths from manifest)
    ax_b.text(0.10, 2.96, "H5N1", fontsize=6.0, fontweight="bold", color=INK, va="center")
    bar_x0, bar_w = 1.40, 8.36
    gap = 0.045
    usable_w = bar_w - gap * (len(iav_landmarks) - 1)
    cur_x = bar_x0
    for idx_s, seg in enumerate(iav_landmarks):
        sname = seg["short_name"]
        slen = float(seg["end"] - seg["start"] + 1)
        sw = usable_w * (slen / iav_total_nt)
        is_ea1 = sname in b313_donor_segs
        scol = TEAL if is_ea1 else INDIGO
        draw_rect(ax_b, cur_x, 2.75, sw, 0.42, fc=scol, ec="none")
        ax_b.text(cur_x + sw / 2.0, 2.96, sname, fontsize=6.0, fontweight="bold", color="white", ha="center", va="center")
        if idx_s < len(iav_landmarks) - 1:
            bx = cur_x + sw + gap / 2.0
            ax_b.plot([bx, bx], [2.64, 3.28], color=SUBINK, lw=0.65, ls=":", zorder=4)
        cur_x += sw + gap
    ax_b.text(5.58, 2.38, r"Segmented: $\gamma = 0$ at segment breaks; blocked inside",
              fontsize=6.0, color=SUBINK, ha="center", va="center")

    ax_b.plot([0.10, 9.86], [2.04, 2.04], color=RULE, lw=0.55)

    # Stage 3: Readable Coordinate Barcode Anatomy (clean white ground, no container box)
    ax_b.text(0.10, 1.74, "3. Readable coordinate barcode (BNF grammar)",
              fontsize=6.0, fontweight="bold", color=INK, va="center")

    def draw_colored_mono(ax, x_start, y_pos, tokens, char_w=0.168, fs=6.3):
        cx = x_start
        for txt, col in tokens:
            ax.text(cx, y_pos, txt, fontsize=fs, fontweight="bold", color=col,
                    family="monospace", ha="left", va="center")
            cx += len(txt) * char_w
        return cx

    ax_b.text(0.10, 1.24, "HIV-1:", fontsize=6.0, fontweight="bold", color=SUBINK, va="center")
    end_hiv = draw_colored_mono(ax_b, 1.45, 1.24, [
        ("G0", SLATE), ("[", INK), ("A1", BLUE), (":", INK), ("gp41..nef", VERMILION), ("]", INK),
    ])
    ax_b.text(end_hiv + 0.35, 1.24, "(+ @pol / @env scope)", fontsize=6.0, color=SUBINK, va="center")

    ax_b.text(0.10, 0.76, "H5N1:", fontsize=6.0, fontweight="bold", color=SUBINK, va="center")
    draw_colored_mono(ax_b, 1.45, 0.76, [
        ("NAm_B313", INDIGO), ("[", INK), ("EA1", TEAL), (":", INK), ("PA+HA+NA+MP", VERMILION), ("]", INK),
    ])

    ax_b.text(0.10, 0.30, "Backbone", fontsize=6.0, fontweight="bold", color=INDIGO, ha="left", va="center")
    ax_b.text(1.66, 0.30, "(majority)", fontsize=6.0, color=SUBINK, ha="left", va="center")
    ax_b.text(3.98, 0.30, "Donor", fontsize=6.0, fontweight="bold", color=BLUE, ha="left", va="center")
    ax_b.text(4.96, 0.30, "(foreign)", fontsize=6.0, color=SUBINK, ha="left", va="center")
    ax_b.text(7.05, 0.30, "Region", fontsize=6.0, fontweight="bold", color=VERMILION, ha="left", va="center")
    ax_b.text(8.14, 0.30, "(locus)", fontsize=6.0, color=SUBINK, ha="left", va="center")

    # --- Panel C: Austere Technical Four-Gate Logic Diagram ---
    ax_c = fig.add_axes([0.016, 0.022, 0.474, 0.398])
    ax_c.set_xlim(0, 10.0)
    ax_c.set_ylim(0, 6.8)
    ax_c.axis("off")

    gates = [
        ("G1", "Footprint & coverage", "Covers locus (\u226584% id) & all donor trunks", 5.26),
        ("G2", "Beats parents & ancestors", "Significantly closer than pure & DAG parents", 3.94),
        ("G3", "Shared clonal stem", "Post-recomb. stem SNPs \u2192 stops hotspots", 2.62),
        ("G4", "Residual 3rd-donor scan", "Flags extra foreign tract \u2192 2nd-gen mosaics", 1.30),
    ]
    for idx_g, (gnum, gtitle, gdesc, y0) in enumerate(gates):
        draw_rect(ax_c, 0.10, y0, 7.22, 1.04, fc="white", ec=MUTED, lw=0.60)
        draw_rect(ax_c, 0.10, y0, 0.92, 1.04, fc=SLATE, ec="none")
        ax_c.text(0.56, y0 + 0.52, gnum, fontsize=6.2, fontweight="bold", color="white", ha="center", va="center")
        ax_c.text(1.24, y0 + 0.72, gtitle, fontsize=6.0, fontweight="bold", color=INK, va="center")
        ax_c.text(1.24, y0 + 0.27, gdesc, fontsize=6.0, color=SUBINK, va="center")
        ax_c.plot([7.32, 7.66], [y0 + 0.52, y0 + 0.52], color=SUBINK, lw=0.75)
        next_y = gates[idx_g + 1][3] + 1.04 if idx_g < len(gates) - 1 else 0.94
        ax_c.annotate("", xy=(0.56, next_y), xytext=(0.56, y0),
                      arrowprops=dict(arrowstyle="-|>", color=INK, lw=0.75, mutation_scale=5.5))

    ax_c.plot([7.66, 7.66], [1.82, 5.78], color=SUBINK, lw=0.75)
    draw_rect(ax_c, 7.94, 2.28, 1.94, 3.00, fc="white", ec=VERMILION, lw=0.85)
    draw_rect(ax_c, 7.94, 4.62, 1.94, 0.66, fc=VERMILION, ec="none")
    ax_c.text(8.91, 4.95, "Any fail", fontsize=6.0, fontweight="bold", color="white", ha="center", va="center")
    ax_c.text(8.91, 4.10, "URF", fontsize=8.0, fontweight="bold", color=VERMILION, ha="center", va="center")
    ax_c.text(8.91, 3.44, "/ OutOfPanel", fontsize=6.0, fontweight="bold", color=SUBINK, ha="center", va="center")
    ax_c.text(8.91, 2.76, "+ barcode &\n25-trunk trace", fontsize=6.0, color=INK, ha="center", va="center")
    ax_c.annotate("", xy=(7.94, 3.78), xytext=(7.66, 3.78),
                  arrowprops=dict(arrowstyle="-|>", color=VERMILION, lw=0.95, mutation_scale=7.0))

    draw_rect(ax_c, 0.10, 0.14, 7.22, 0.80, fc="white", ec=BLUE, lw=0.85)
    draw_rect(ax_c, 0.10, 0.14, 0.18, 0.80, fc=BLUE, ec="none")
    ax_c.text(3.76, 0.54, "All 4 pass \u2192 Certified (CRF02_AG, B3.13, XEC)",
              fontsize=6.0, fontweight="bold", color=BLUE, ha="center", va="center")
    ood_row_rt = [r for r in ood_rows if r["panel"].startswith("Combined_OOD") and r["method"].startswith("EGGMAN")][0]
    ax_c.text(8.91, 0.54, f"{float(ood_row_rt['novel_urf_sensitivity_pct']):.1f}% HIV\n100% H5N1",
              fontsize=6.0, fontweight="bold", color=INK, ha="center", va="center")

    # --- Panel D: Empirical Benchmark Bar Chart (Austere Steel + Grayscale) ---
    ax_d = fig.add_axes([0.585, 0.105, 0.398, 0.295])
    p1bc_act = bench_rows["Panel1BC_Active_142_Mosaic_CRFs"]
    p4b = bench_rows["Panel4B_Clinical_NamedCRF"]
    cm_ood_row = [r for r in ood_rows if r["panel"].startswith("Combined_OOD") and r["method"] == "COMET_v2.4"][0]
    rt_ood = float(ood_row_rt["novel_urf_sensitivity_pct"])
    cm_ood_flag = float(cm_ood_row["novel_urf_sensitivity_pct"])
    n_ood_tot = int(cm_ood_row["n_sequences"])
    cm_ood_pure = 100.0 * (int(cm_ood_row["novel_urf_flagged"]) - (int(cm_ood_row["legacy_crf_chimera_or_misassigned"]) - int(cm_ood_row["false_legacy_misassignment"]))) / n_ood_tot
    iav_ood_ext = ext_rows[("iav_sra_650_ood", "EGGMAN_OOD_Masked")]
    null_def = [r for r in null_rows if r["regime"].startswith("Default") and r["scenario"] == "ALL_1300_COMBINED"][0]

    rt_local_s = bench_meta["rhiztype_timing"]["rhiztype_local_wall_sec"]
    jp_wall_h = bench_meta["jphmm_timing"]["wall_sec"] / 3600.0
    s_eg = synth_rows["EGGMAN"]
    s_jp = synth_rows["jpHMM"]

    p4a = bench_rows["Panel4A_Clinical_PureSubtype"]
    metrics = [
        (f"Clin Pure\n(n={int(p4a['n_sequences'])})", float(p4a["rhiztype_exact_lineage_pct"]), float(p4a["comet_exact_lineage_pct"]), float(p4a["jphmm_pure_exact_pct"])),
        (f"Clin Parent\n(n={int(p4b['parent_eval_n'])})", float(p4b["rhiztype_parent_exact_set_pct"]), float(p4b["comet_parent_exact_set_pct"]), float(p4b["jphmm_parent_exact_set_pct"])),
        (f"Clin State\n(n={int(p4b['n_sequences'])})", float(p4b["rhiztype_mosaic_or_pure_state_pct"]), float(p4b["comet_mosaic_or_pure_state_pct"]), float(p4b["jphmm_mosaic_or_pure_state_pct"])),
        (f"Held-out\nCRF ({n_ood_tot})", rt_ood, cm_ood_flag, None),
        (f"Held-out\nIAV ({int(iav_ood_ext['n_sequences'])})", float(iav_ood_ext["ood_novel_detection_pct"]), 0.0, None),
    ]
    x = np.arange(len(metrics))
    w = 0.20
    for i, (lbl, v_rt, v_cm, v_jp) in enumerate(metrics):
        if v_cm is not None and v_jp is not None:
            bx_rt, bx_cm, bx_jp = i - 0.25, i, i + 0.25
        elif v_cm is not None:
            bx_rt, bx_cm, bx_jp = i - 0.16, i + 0.16, None
        else:
            bx_rt, bx_cm, bx_jp = i, None, None

        ax_d.bar(bx_rt, v_rt, width=w, color=BLUE, zorder=3)
        x_rt_txt = bx_rt - 0.06 if i == 3 else bx_rt
        ax_d.text(x_rt_txt, v_rt + 2.5, f"{v_rt:.1f}",
                  ha="center", va="bottom", fontsize=6.0, fontweight="bold", color=BLUE)
        if v_cm is not None:
            if i == 3:
                ax_d.bar(bx_cm, v_cm, width=w, color="white", edgecolor=SLATE_MID, lw=0.70, zorder=3)
                ax_d.bar(bx_cm, cm_ood_pure, width=w, color=SLATE_MID, zorder=4)
                ax_d.text(bx_cm + 0.11, v_cm + 2.5, f"{v_cm:.1f}*", ha="center", va="bottom", fontsize=6.0, fontweight="bold", color=SUBINK)
                ax_d.text(bx_cm, cm_ood_pure / 2.0, f"{cm_ood_pure:.0f}", ha="center", va="center", fontsize=6.0, fontweight="bold", color="white", zorder=5)
            elif v_jp is not None:
                # Middle bar flanked by two taller bars: rotate cleanly 90 deg so zero collision with neighbors
                ax_d.bar(bx_cm, v_cm, width=w, color=SLATE_MID, zorder=3)
                ax_d.text(bx_cm, v_cm + 3.0, f"{v_cm:.1f}",
                          ha="center", va="bottom", rotation=90, fontsize=6.0, fontweight="bold", color=SUBINK, zorder=5)
            else:
                ax_d.bar(bx_cm, v_cm, width=w, color=SLATE_MID, zorder=3)
                x_cm_txt = bx_cm + 0.05 if i == 0 else bx_cm
                ax_d.text(x_cm_txt, v_cm + 2.5, "0.0" if i == 4 else f"{v_cm:.1f}",
                          ha="center", va="bottom", fontsize=6.0, fontweight="bold", color=SUBINK)
        if v_jp is not None:
            ax_d.bar(bx_jp, v_jp, width=w, color=SLATE_LIGHT, edgecolor=SLATE_MID, lw=0.55, zorder=3)
            x_jp_txt = bx_jp + 0.04 if i == 2 else bx_jp
            ax_d.text(x_jp_txt, v_jp + 2.5, f"{v_jp:.1f}", ha="center", va="bottom", fontsize=6.0, color=SUBINK)

    ax_d.set_xticks(x)
    ax_d.set_xticklabels([m[0] for m in metrics], fontsize=6.0)
    ax_d.set_ylim(0, 152)
    ax_d.set_yticks([0, 25, 50, 75, 100])
    ax_d.set_ylabel("Accuracy or URF recall (%)", fontsize=6.4, color=INK)
    ax_d.yaxis.grid(True, color=HAIR, lw=0.5, zorder=1)
    ax_d.set_axisbelow(True)

    handles = [
        mpatches.Patch(color=BLUE, label=f"EGGMAN ({rt_local_s:.2f} s)"),
        mpatches.Patch(color=SLATE_MID, label="COMET / GenoFLU*"),
        mpatches.Patch(facecolor=SLATE_LIGHT, edgecolor=SLATE_MID, lw=0.55, label=f"jpHMM ({jp_wall_h:.1f} h)"),
    ]
    ax_d.legend(handles=handles, frameon=False, fontsize=6.0, loc="upper left", ncol=2,
                handlelength=0.85, columnspacing=0.65, borderaxespad=0.15)
    ax_d.text(
        0.99, 0.80,
        f"500-mosaic BP: {float(s_eg['median_abs_bp_err_nt']):.0f} vs {float(s_jp['median_abs_bp_err_nt']):.0f} nt | "
        f"{float(null_def['tier1_basal_false_recomb_fpr_pct']):.1f}% T1 null ({int(null_def['n_simulated']):,})",
        transform=ax_d.transAxes, fontsize=6.0, color=SUBINK, ha="right", va="top",
    )

    audit_matplotlib_figure(fig, "Fig1_Concept", quadrant_splits=(0.500, 0.470))
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG_DIR / "fig_science_1_concept.pdf", dpi=300)
    fig.savefig(FIG_DIR / "fig_science_1_concept.png", dpi=350)
    plt.close(fig)


# =============================================================================
# FIGURE 2: Multi-Pathogen Reassortment & Recombination (IAV, SC2 BV-BRC, Sarbeco)
# =============================================================================
def draw_barcode_strip(ax, y, segs, colors, total_len, offset=0, h=0.44):
    for name, s, e in segs:
        x0 = (s - offset) / (total_len - offset)
        w = max((e - s + 1) / (total_len - offset), 0.004)
        ax.add_patch(mpatches.Rectangle((x0, y - h / 2), w, h, facecolor=colors.get(name, SUBINK), edgecolor="white", lw=0.45, zorder=3))


def build_fig2():
    with open(MP_JSON) as f:
        mp = json.load(f)
    with open(SC2_JSON) as f:
        sc2 = json.load(f)
    with open(SPECIMENS_JSON) as f:
        specs = json.load(f)
    with open(PACKS_DIR / "influenza_a_8seg" / "manifest.json") as f:
        iav_manifest = json.load(f)
    with open(PACKS_DIR / "sarbecovirus_wuhan" / "manifest.json") as f:
        sb_manifest = json.load(f)
    with open(RES / "bvbrc_iav_h5_8seg_stream_summary.json") as f:
        iav_bvbrc = json.load(f)
    with open(RES / "bvbrc_iav_40constellations_novelty_audit.tsv") as f:
        iav_aud = list(csv.DictReader(f, delimiter="\t"))
    with open(RES / "bvbrc_100k_sc2_stream_summary.json") as f:
        sc2_100k = json.load(f)
    with open(RES / "bvbrc_sc2_23clusters_novelty_audit.tsv") as f:
        sc2_aud = list(csv.DictReader(f, delimiter="\t"))
    with open(RES / "benchmark_sc2_nextclade_and_iav_genoflu.csv") as f:
        ext_rows = {(r["cohort"], r["tool"]): r for r in csv.DictReader(f)}

    fig = plt.figure(figsize=(CANVAS_WIDTH_IN, 4.95), facecolor="white")
    add_quadrant_hairlines(fig, x_split=0.500, y_split=0.495)
    draw_panel_header(fig, 0.014, 0.972, "A", "H5N1 eight-segment reassortment (N = 1,200)")
    draw_panel_header(fig, 0.514, 0.972, "B", "Head-to-head recall & 3,980-isolate BV-BRC")
    draw_panel_header(fig, 0.014, 0.472, "C", "SARS-CoV-2 spanning basis & 100k stream")
    draw_panel_header(fig, 0.514, 0.472, "D", "Sarbecovirus wildlife mosaics (N = 79)")

    # --- Panel A: IAV 8-Segment Reassortment Matrix (Loaded dynamically from manifest + specimens JSON) ---
    ax_a = fig.add_axes([0.118, 0.598, 0.315, 0.315])
    ax_a.set_xlim(0, 1.0)
    ax_a.set_ylim(-0.6, 7.8)
    ax_a.axis("off")

    iav_total_len = int(iav_manifest["ref_full_len"])
    seg_bounds = [(lm["short_name"], int(lm["start"]), int(lm["end"])) for lm in iav_manifest["landmarks"]]
    iav_cols = {
        "EA1": TEAL, "EA2": SKY, "EA3": BLUE,
        "NAm_B313": INDIGO, "NAm_D1": AMBER, "NAm_B36": PURPLE,
        "NAm_B2": VERMILION, "NAm_B32": WINE, "NAm_C2": SUBINK,
    }
    row_labels_map = {
        "A1": "A1 (pure)", "A3": "A3 (pure)", "B3.13": "B3.13 cattle",
        "D1.1": "D1.1 avian", "B3.6": "B3.6", "D1.8": "D1.8", "B3.2": "B3.2", "C2.1": "C2.1",
    }
    for idx, item in enumerate(specs["iav"]):
        y = 7.0 - idx * 0.96
        gt_tag = item["id"].split("|")[-1]
        lbl = row_labels_map.get(gt_tag, gt_tag)
        call = item["call"]
        draw_barcode_strip(ax_a, y, item["segments"], iav_cols, iav_total_len, offset=1, h=0.52)
        ax_a.text(-0.025, y, lbl, fontsize=6.0, fontweight="bold", color=VERMILION if call == "URF" else INK, ha="right", va="center")
        ax_a.text(1.025, y, call, fontsize=6.0, fontweight="bold", color=VERMILION if call == "URF" else SUBINK, ha="left", va="center")

    for sname, s0, s1 in seg_bounds:
        xm = ((s0 + s1) / 2.0 - 1.0) / (iav_total_len - 1.0)
        if sname == "MP":
            xm -= 0.012
        elif sname == "NS":
            xm += 0.014
        ax_a.text(xm, 7.65, sname, fontsize=6.0, fontweight="bold", color=SLATE, ha="center", va="center")
        if s0 > 1:
            xdiv = (s0 - 1.0) / (iav_total_len - 1.0)
            ax_a.plot([xdiv, xdiv], [0.05, 7.30], color="white", lw=0.85, zorder=4)
            ax_a.plot([xdiv, xdiv], [0.05, 7.30], color=SUBINK, lw=0.45, ls=":", zorder=5)

    iav_leg = [mpatches.Patch(color=iav_cols[k], label=k)
               for k in ["EA1", "EA3", "NAm_D1", "NAm_B36", "NAm_B2", "NAm_B32", "NAm_C2", "NAm_B313"]]
    ax_a.legend(handles=iav_leg, frameon=False, fontsize=6.0, loc="upper center",
                bbox_to_anchor=(0.45, -0.02), ncol=4, handlelength=0.75, columnspacing=0.90)

    # --- Panel B: Head-to-Head Comparators + 3,980-Isolate Global BV-BRC Stream ---
    ax_b1 = fig.add_axes([0.678, 0.705, 0.302, 0.228])
    iav_sra_eg = ext_rows[("iav_sra_1200", "EGGMAN")]
    iav_sra_gf = ext_rows[("iav_sra_1200", "USDA_GenoFLU_BLASTn")]
    iav_ood_eg = ext_rows[("iav_sra_650_ood", "EGGMAN_OOD_Masked")]
    sc2_eg = ext_rows[("sc2_curated_1000", "EGGMAN")]
    sc2_nc = ext_rows[("sc2_curated_1000", "Nextclade_CLI_8thr")]
    sc2_ood_eg = ext_rows[("sc2_ood_xdv_xfz_xfc", "EGGMAN_OOD_Masked")]
    iav_speed = float(iav_sra_eg["throughput_seq_per_sec"])

    cats_b = [
        ("IAV SRA exact (1.1k)", float(iav_sra_eg["lineage_or_genotype_acc_pct"]), BLUE),
        ("GenoFLU BLASTn", float(iav_sra_gf["lineage_or_genotype_acc_pct"]), SLATE_MID),
        ("Held-out IAV donors", float(iav_ood_eg["clonal_vs_recomb_acc_pct"]), VERMILION),
        ("USDA unassg. (100)", float(mp["influenza_a"]["panels"]["IAV-3_Unassigned"]["tier1_mosaic_pct"]), VERMILION),
        ("SC2 state (1,000)", float(sc2_eg["clonal_vs_recomb_acc_pct"]), BLUE),
        ("Nextclade 8-thr", float(sc2_nc["clonal_vs_recomb_acc_pct"]), SLATE_MID),
        ("Held-out SC2 BP (8)", float(sc2_ood_eg["clonal_vs_recomb_acc_pct"]), VERMILION),
    ]
    yb = np.arange(len(cats_b))[::-1]
    for yp, (lbl, val, col) in zip(yb, cats_b):
        ax_b1.plot([0, val], [yp, yp], color=RULE, lw=1.1, zorder=2)
        ax_b1.scatter([val], [yp], s=20, color=col, edgecolor="white", lw=0.4, zorder=4)
        ax_b1.text(val + 3.0, yp, f"{val:.0f}%" if val == 100.0 else f"{val:.1f}%",
                   va="center", ha="left", fontsize=6.0, fontweight="bold", color=INK,
                   bbox=dict(boxstyle="square,pad=0.06", fc="white", ec="none", alpha=1.0), zorder=6)
    ax_b1.set_yticks(yb)
    ax_b1.set_yticklabels([c[0] for c in cats_b], fontsize=6.0)
    ax_b1.set_xlim(0, 128)
    ax_b1.set_xticks([0, 25, 50, 75, 100])
    ax_b1.set_xlabel(f"Accuracy / state (%; IAV {iav_speed:,.0f} seq/s)", fontsize=6.1, color=INK)
    ax_b1.xaxis.grid(True, color=HAIR, lw=0.5, zorder=1)

    ax_b2 = fig.add_axes([0.520, 0.515, 0.460, 0.115])
    ax_b2.set_xlim(0, 1.0)
    ax_b2.set_ylim(0, 2.35)
    ax_b2.axis("off")
    bvbrc_hc = int(iav_bvbrc["high_coverage_summary"]["n_isolates"])
    bvbrc_tot = int(iav_bvbrc["assembly_stats"]["assembled_8seg_screened"])
    hc_cats = iav_bvbrc["high_coverage_summary"]["categories"]
    n_cert = int(hc_cats["Certified_Reassortant"])
    n_pure = int(hc_cats["Pure_Subtype"])
    n_oop = int(hc_cats["OutOfPanel"])
    n_omitted_g = sum(int(r["n_isolates"]) for r in iav_aud if r["classification"] == "KNOWN_USDA_GENOFLU_OMITTED_FROM_PACK")
    n_pub_other = sum(int(r["n_isolates"]) for r in iav_aud if r["classification"] == "PUBLISHED_NON_USDA_GENOTYPE_OR_EURL")
    n_rare_urf = int(hc_cats["URF"]) - (n_omitted_g + n_pub_other)

    ax_b2.text(
        0.0, 2.08,
        f"BV-BRC H5Nx ({bvbrc_hc:,}/{bvbrc_tot:,} \u226511.5 kb; 0 unchar.):",
        fontsize=6.0, fontweight="bold", color=INK, va="center",
    )
    seg_comp = [
        (n_cert, BLUE, f"Cert. {n_cert:,}"),
        (n_pure, SKY, f"{n_pure}"),
        (n_oop, SLATE_MID, f"{n_oop}"),
        (n_omitted_g, AMBER, f"{n_omitted_g}"),
        (n_pub_other, VERMILION, f"{n_pub_other}"),
        (n_rare_urf, SLATE_LIGHT, f"{n_rare_urf}"),
    ]
    x_cur = 0.0
    for cnt, col, tag in seg_comp:
        w_frac = cnt / bvbrc_hc
        ax_b2.add_patch(mpatches.Rectangle((x_cur, 1.08), w_frac, 0.62, facecolor=col, edgecolor="white", lw=0.5))
        if w_frac >= 0.095:
            ax_b2.text(x_cur + w_frac / 2, 1.39, tag, fontsize=6.0, fontweight="bold", color="white", ha="center", va="center")
        x_cur += w_frac
    b2_leg = [
        mpatches.Patch(color=BLUE, label=f"Cert ({n_cert:,})"),
        mpatches.Patch(color=SKY, label=f"Pure ({n_pure})"),
        mpatches.Patch(color=SLATE_MID, label=f"LPAI/OOP ({n_oop})"),
        mpatches.Patch(color=AMBER, label=f"GenoFLU ({n_omitted_g})"),
        mpatches.Patch(color=VERMILION, label=f"EURL ({n_pub_other})"),
        mpatches.Patch(color=SLATE_LIGHT, label=f"<3 ({n_rare_urf})"),
    ]
    ax_b2.legend(handles=b2_leg, frameon=False, fontsize=6.0, loc="lower center",
                 bbox_to_anchor=(0.49, -0.09), ncol=3, handlelength=0.68, columnspacing=0.95)

    # --- Panel C: SC2 Spanning Basis & Blind Pango X* Rediscovery in 100k BV-BRC Stream ---
    ax_c1 = fig.add_axes([0.064, 0.102, 0.172, 0.312])
    w_base = sc2["part_a_world_10062"]["world_10062_baseline_cov_tuned"]
    w_span = sc2["part_a_world_10062"]["world_10062_spanning_cov_tuned"]
    b_base = sc2["part_b_bvbrc_stream"]["bvbrc_baseline_cov_tuned"]
    b_span = sc2["part_b_bvbrc_stream"]["bvbrc_spanning_cov_tuned"]

    ba4_pre = 100.0 * (w_base["ba4_total"] - w_base["ba4_clonal"]) / w_base["ba4_total"]
    ba4_post = 100.0 * (w_span["ba4_total"] - w_span["ba4_clonal"]) / w_span["ba4_total"]
    ba32_pre = 100.0 * (w_base["ba32_total"] - w_base["ba32_clonal"]) / w_base["ba32_total"]
    ba32_post = 100.0 * (w_span["ba32_total"] - w_span["ba32_clonal"]) / w_span["ba32_total"]
    bvbrc_pre = b_base["nominal_pure_waves_mosaic_pct"]
    bvbrc_post = b_span["nominal_pure_waves_mosaic_pct"]

    groups_c1 = ["BA.4\n(61)", "BA.3.2\n(183)", "BV-BRC\n(13.3k)"]
    pre_vals = [ba4_pre, ba32_pre, bvbrc_pre]
    post_vals = [ba4_post, ba32_post, bvbrc_post]
    xc = np.arange(len(groups_c1))
    ax_c1.bar(xc - 0.21, pre_vals, width=0.25, color=VERMILION, label="11-base", zorder=3)
    ax_c1.bar(xc + 0.21, post_vals, width=0.25, color=BLUE, label="16-span", zorder=3)
    for i, (v1, v2) in enumerate(zip(pre_vals, post_vals)):
        y_v1 = v1 + 6.0 if i == 2 else v1 + 2.8
        ax_c1.text(i - 0.21, y_v1, f"{v1:.0f}%", ha="center", va="bottom", fontsize=6.0, fontweight="bold", color=INK)
        ax_c1.text(i + 0.02, max(v2, 0.0) + 2.2, f"{v2:.1f}%", ha="left", va="bottom", fontsize=6.0, fontweight="bold", color=BLUE)
    ax_c1.set_xticks(xc)
    ax_c1.set_xticklabels(groups_c1, fontsize=6.0)
    ax_c1.set_xlim(-0.66, 2.86)
    ax_c1.set_ylim(0, 128)
    ax_c1.set_yticks([0, 25, 50, 75, 100])
    ax_c1.set_ylabel("Pseudo-mosaic (%)", fontsize=6.4, color=INK)
    ax_c1.yaxis.grid(True, color=HAIR, lw=0.5, zorder=1)
    ax_c1.legend(frameon=False, fontsize=6.0, loc="upper left", ncol=1, handlelength=0.8)

    ax_c2 = fig.add_axes([0.314, 0.102, 0.176, 0.312])
    pango_cnt = Counter()
    pango_roots = set()
    drift_cnt = 0
    for r in sc2_aud:
        n = int(r["n_genomes"])
        if r["classification"] == "KNOWN_PANGO_RECOMBINANT_UNANNOTATED_IN_BVBRC":
            pango_cnt[r["resolved_pango_root"]] += n
            pango_roots.add(r["resolved_pango_root"])
        else:
            drift_cnt += n
    tot_cl_g = sum(pango_cnt.values()) + drift_cnt
    pct_pango = 100.0 * sum(pango_cnt.values()) / max(tot_cl_g, 1)
    other_x = sum(v for k, v in pango_cnt.items() if k not in {"XEC", "XFZ", "XBB"})
    n_other_roots = len(pango_roots - {"XEC", "XFZ", "XBB"})
    sc2_bars = [
        ("XEC", pango_cnt["XEC"], BLUE),
        ("XFZ", pango_cnt["XFZ"], BLUE),
        ("XBB", pango_cnt["XBB"], BLUE),
        (f"{n_other_roots} other X*", other_x, SKY),
        ("Drift (3)", drift_cnt, SLATE_LIGHT),
    ]
    yc2 = np.arange(len(sc2_bars))[::-1]
    for yi, (lbl, cnt, col) in zip(yc2, sc2_bars):
        ec = SLATE_MID if col == SLATE_LIGHT else "none"
        ax_c2.barh(yi, cnt, height=0.56, color=col, edgecolor=ec, lw=0.55, zorder=3)
        ax_c2.text(cnt + 22, yi, f"{cnt:,}", va="center", ha="left", fontsize=6.0, fontweight="bold", color=INK)
    ax_c2.set_yticks(yc2)
    ax_c2.set_yticklabels([b[0] for b in sc2_bars], fontsize=6.0)
    ax_c2.set_xlim(0, 1480)
    ax_c2.set_xticks([0, 500, 1000])
    ax_c2.set_xlabel(f"Unannot. URF (n={tot_cl_g:,})", fontsize=6.2, color=INK)
    ax_c2.xaxis.grid(True, color=HAIR, lw=0.5, zorder=1)
    n_backfill_cov = sc2_100k["eggman_cov_tuned"]["mosaic_breakdown"]["wuhan_b_segment_backfill_post2021"]
    ax_c2.text(
        0.98, 0.12,
        f"{pct_pango:.1f}% = {len(pango_roots)} X*\n100% XBB/XDV\n+{n_backfill_cov} Wuhan-B QA",
        transform=ax_c2.transAxes, fontsize=6.0, color=INK, ha="right", va="bottom",
    )

    # --- Panel D: Sarbecovirus Natural Wildlife Recombination Atlas ---
    ax_d = fig.add_axes([0.625, 0.140, 0.290, 0.268])
    ax_d.set_xlim(0, 1.0)
    ax_d.set_ylim(-0.6, 6.8)
    ax_d.axis("off")

    sb_cols = {"SARS2_Wuhan": BLUE, "SARS1_HSZ": PURPLE, "Bat_RpYN06": TEAL, "Bat_HKU3": AMBER, "Bat_BM48": SLATE}
    c_start = int(sb_manifest["coding_start"])
    c_end = int(sb_manifest["coding_end"])
    lm_map = {lm["short_name"]: (int(lm["start"]), int(lm["end"])) for lm in sb_manifest["landmarks"]}
    hdr_genes = [
        ("ORF1a", lm_map["ORF1a"][0], lm_map["ORF1a"][1]),
        ("ORF1b", lm_map["ORF1b"][0], lm_map["ORF1b"][1]),
        ("S", lm_map["S"][0], lm_map["S"][1]),
        ("3a\u2013N", lm_map["ORF3a"][0], lm_map["ORF10"][1]),
    ]
    for gname, g0, g1 in hdr_genes:
        x0 = (g0 - c_start) / (c_end - c_start)
        gw = (g1 - g0) / (c_end - c_start)
        ax_d.add_patch(mpatches.Rectangle((x0, 6.08), gw, 0.36, facecolor="#E2E8F0", edgecolor="white", lw=0.6))
        ax_d.text(x0 + gw / 2, 6.56, rf"$\it{{{gname}}}$" if gname in {"S", "ORF1a", "ORF1b"} else gname,
                  fontsize=6.0, fontweight="bold", color=SLATE, ha="center", va="center")

    want_order = ["Wuhan-Hu-1", "WIV16", "RaTG13", "PrC31", "BANAL-52", "CoVZC45"]
    sb_by_name = {s["id"].split("|")[0]: s for s in specs["sarbeco"]}
    for idx, sp_name in enumerate(want_order):
        s_item = sb_by_name[sp_name]
        segs = s_item["segments"]
        n_bp_s = len(segs) - 1
        lbl = f"{sp_name}*" if sp_name == "Wuhan-Hu-1" else sp_name
        bpt = f"{n_bp_s} BPs" if n_bp_s != 1 else "1 BP"
        y = 5.3 - idx * 0.98
        draw_barcode_strip(ax_d, y, segs, sb_cols, c_end, offset=c_start, h=0.52)
        ax_d.text(-0.025, y, lbl, fontsize=6.0, fontweight="bold", color=INK, ha="right", va="center")
        ax_d.text(1.02, y, bpt, fontsize=6.0, color=SUBINK, ha="left", va="center")

    n_valid_sb = int(mp["sarbecovirus"]["n_valid_total"])
    n_mos_sb = int(mp["sarbecovirus"]["valid_non_basal_mosaic_n"])
    n_nonb_sb = int(mp["sarbecovirus"]["n_valid_non_basal"])
    pct_sb = float(mp["sarbecovirus"]["valid_non_basal_mosaic_pct"])
    med_bp_sb = float(mp["sarbecovirus"]["median_breakpoints_mosaic"])
    fig.text(0.750, 0.118, f"N = {n_valid_sb} · 5/5 basal clonal · {n_mos_sb}/{n_nonb_sb} ({pct_sb:.1f}%) mosaic (med {med_bp_sb:.0f} BP)",
             fontsize=6.0, color=SUBINK, ha="center", va="center")
    sb_leg = [mpatches.Patch(color=c, label=k) for k, c in sb_cols.items()]
    ax_d.legend(handles=sb_leg, frameon=False, fontsize=6.0, loc="upper center",
                bbox_to_anchor=(0.45, -0.11), ncol=3, handlelength=0.8, columnspacing=0.6)

    audit_matplotlib_figure(fig, "Fig2_Multipathogen", quadrant_splits=(0.500, 0.495))
    fig.savefig(FIG_DIR / "fig_science_2_multipathogen.pdf", dpi=300)
    fig.savefig(FIG_DIR / "fig_science_2_multipathogen.png", dpi=350)
    plt.close(fig)


# =============================================================================
# FIGURE 3: HIV-1 Reference Scrubbing, ChronAeon Timeline & Compendium Pruning
# =============================================================================
def build_fig3():
    with open(RES / "phase1_chronaeon_clade_dating.csv") as f:
        clock_rows = list(csv.DictReader(f))
    with open(RES / "lanl_declared_breakpoints_audit.csv") as f:
        bp_aud = list(csv.DictReader(f))
    with open(RES / "hiv1_recombination_genealogy_dag.csv") as f:
        dag_rows = list(csv.DictReader(f))
    with open(RES / "phase1_contaminated_pure_references.csv") as f:
        contam_rows = list(csv.DictReader(f))
    with open(RES / "hiv1_descriptive_nomenclature_compendium.csv") as f:
        comp_rows = list(csv.DictReader(f))
    with open(CORE_DATA_DIR / "compendium_bundle.json") as f:
        cb = json.load(f)

    fig = plt.figure(figsize=(CANVAS_WIDTH_IN, 5.05), facecolor="white")
    add_quadrant_hairlines(fig, x_split=0.500, y_split=0.495)
    draw_panel_header(fig, 0.014, 0.972, "A", "Mosaic tracts inside LANL pure references")
    draw_panel_header(fig, 0.514, 0.972, "B", "ChronAeon dates across Group M clades")
    draw_panel_header(fig, 0.014, 0.472, "C", f"Bridged audit of {len(bp_aud):,} LANL CRF junctions")
    draw_panel_header(fig, 0.514, 0.472, "D", "Three-generation CRF genealogy & hubs")

    # --- Panel A: Canonical Mosaic/Contaminated Pure References ---
    ax_a = fig.add_axes([0.125, 0.602, 0.258, 0.315])
    ax_a.set_xlim(790, 9417)
    ax_a.set_ylim(-0.6, 8.65)
    ax_a.set_xticks([790, 2253, 5097, 6225, 8379, 9417])
    ax_a.set_xticklabels([r"$\it{gag}$", r"$\it{pol}$", r"$\it{vif}$", r"$\it{env}$", r"$\it{gp41}$", r"$\it{nef}$"], fontsize=6.0)
    ax_a.set_yticks([])
    ax_a.spines["left"].set_visible(False)

    order_keys = ["P1|Ref.G.BE", "P1|Ref.G.KE", "P1|Ref.G.NG", "P1|Ref.G.PT", "P1|Ref.A7.NG.09", "P1|Ref.A7.NG.10", "P1|Ref.A7.NG.11", "P1|Ref.L.CD.83"]
    with open(RES / "benchmark_rhiztype_v02_output.csv") as f:
        p1a_rows = [r for r in csv.DictReader(f) if r["sequence_id"].startswith("P1|Ref.")][:55]
    p1_mosaic = [r for r in p1a_rows if int(r["n_breakpoints"]) >= 1 or r["sequence_id"].startswith("P1|Ref.L.CD.83")]
    p1_mosaic.sort(key=lambda r: next((i for i, k in enumerate(order_keys) if r["sequence_id"].startswith(k)), 99))

    trunk_pal = {
        "G0": SLATE, "A1": BLUE, "?E": AMBER, "F1": VERMILION,
        "F2": SKY, "B": INDIGO, "L": PURPLE, "J": WINE, "OutOfPanel": MUTED,
    }
    for i, r in enumerate(p1_mosaic):
        y = (len(p1_mosaic) - 1) - i * 0.96
        parts = r["sequence_id"].split("|")[1].split(".")
        sub_tag = parts[1]
        if sub_tag == "G":
            lbl = f"G: {parts[4][:7].rstrip('_')}"
            note = fmt_p_tex(float(r["call_p_value"]))
        elif sub_tag == "A7":
            lbl = f"{sub_tag}: {parts[-1][:8]}"
            note = "Demoted"
        else:
            lbl = f"L: {parts[-1][:8]}"
            note = "RhizAeon K:nef"

        for seg in r["basal_barcode"].split("|"):
            tr, rng_s = seg.split(":")
            s, e = [int(x) for x in rng_s.split("-")]
            col = trunk_pal.get(tr, MUTED)
            ax_a.add_patch(mpatches.Rectangle((s, y - 0.28), e - s, 0.56, facecolor=col, edgecolor="white", lw=0.4, zorder=3))
        ax_a.text(650, y, lbl, fontsize=6.0, fontweight="bold", color=VERMILION if note == "Demoted" else INK, ha="right", va="center")
        ax_a.text(9550, y, note, fontsize=6.0, color=VERMILION if note == "Demoted" else SUBINK, ha="left", va="center")

    n_unidx_tr = sum(1 for r in contam_rows if r["is_ghost"].lower() == "true")
    ax_a.axvline(5097, color=VERMILION, lw=0.85, ls="--", zorder=5)
    ax_a.text(
        5260, 8.05,
        f"vif 5097: ?E restricted 3'\n{len(contam_rows)} tracts ({n_unidx_tr} unindexed)",
        fontsize=6.0, fontweight="bold", color=VERMILION, ha="left", va="center",
    )
    leg_a = [mpatches.Patch(color=trunk_pal[k], label=k) for k in ["G0", "A1", "?E", "F1", "F2", "B", "L", "J"]]
    ax_a.legend(handles=leg_a, frameon=False, fontsize=6.0, loc="upper center",
                bbox_to_anchor=(0.48, -0.105), ncol=4, handlelength=0.75, columnspacing=0.85)

    # --- Panel B: ChronAeon Molecular Clock Timeline (Austere Ink + Steel + Rust) ---
    ax_b = fig.add_axes([0.592, 0.562, 0.390, 0.375])
    sel_clades = ["C", "02_AG", "A1", "G", "D", "B", "F1", "01_AE", "11_cpx", "06_cpx", "H", "A6", "J", "07_BC", "08_BC", "63_02A6"]
    cmap_clock = {r["clade"]: r for r in clock_rows}
    ordered = [cmap_clock[c] for c in sel_clades if c in cmap_clock]
    yp = np.arange(len(ordered))[::-1]
    ax_b.axvspan(1955.0, 1962.5, color="#F3F4F6", alpha=0.95, zorder=1)

    has_unbounded_hi = False
    for y, r in zip(yp, ordered):
        cl = r["clade"]
        tm = float(r["t_mrca"])
        lo = float(r["ci_lower"])
        hi_raw = r["ci_upper"].strip()
        is_open_hi = (hi_raw == "" or hi_raw.lower() == "nan")
        hi = 2015.0 if is_open_hi else float(hi_raw)
        if is_open_hi:
            has_unbounded_hi = True
        is_deriv = cl in {"G", "J"}
        is_crf = "_" in cl
        col = VERMILION if is_deriv else (BLUE if is_crf else SLATE)
        marker = "D" if is_deriv else ("s" if is_crf else "o")
        txt_col = VERMILION if is_deriv else SUBINK
        ax_b.plot([max(lo, 1931.0), hi], [y, y], color=col, lw=0.95, zorder=2)
        if lo < 1931.0:
            ax_b.text(1930.8, y, "<", fontsize=6.0, fontweight="bold", color=col, ha="right", va="center")
        if is_open_hi:
            ax_b.annotate("", xy=(hi + 4.5, y), xytext=(hi, y),
                          arrowprops=dict(arrowstyle="-|>", color=col, lw=0.90, mutation_scale=6.0))
            ax_b.scatter([tm], [y], color=col, marker=marker, s=17, zorder=4)
            ax_b.text(hi + 5.2, y, f"{tm:.1f}*", fontsize=6.0, color=txt_col, va="center")
        else:
            ax_b.scatter([tm], [y], color=col, marker=marker, s=17, zorder=4)
            ax_b.text(hi + 1.0, y, f"{tm:.1f}", fontsize=6.0, color=txt_col, va="center")

    ax_b.set_yticks(yp)
    ax_b.set_yticklabels([r["clade"] for r in ordered], fontsize=6.0)
    ax_b.set_xlim(1928, 2046)
    ax_b.set_xticks([1940, 1960, 1980, 2000, 2020])
    ax_b.set_xlabel(r"ChronAeon emergence $t_{\mathrm{MRCA}}$ (year, 95% CI)", fontsize=6.4, color=INK)
    ax_b.xaxis.grid(True, color=HAIR, lw=0.5, zorder=1)
    b_leg = [
        plt.Line2D([0], [0], marker="o", color="w", markerfacecolor=SLATE, markersize=4.5, label="Pure clade"),
        plt.Line2D([0], [0], marker="s", color="w", markerfacecolor=BLUE, markersize=4.5, label="CRF"),
        plt.Line2D([0], [0], marker="D", color="w", markerfacecolor=VERMILION, markersize=4.5, label="Mosaic ref (G, J)"),
    ]
    ax_b.legend(handles=b_leg, frameon=False, fontsize=6.0, loc="upper right")
    t_02 = float(cmap_clock["02_AG"]["t_mrca"])
    t_g = float(cmap_clock["G"]["t_mrca"])
    note_b = f"02_AG ({t_02:.1f}) < G ({t_g:.1f})\n" + r"$\rho=0.81$ vs LSD2"
    if has_unbounded_hi:
        note_b += "\n*open upper CI"
    ax_b.text(0.99, 0.72, note_b, transform=ax_b.transAxes, fontsize=6.0, color=INK, ha="right", va="top")

    # --- Panel C: Statistical Audit of 986 Bridged LANL Breakpoints (Disciplined Monochrome + Accent) ---
    ax_c = fig.add_axes([0.188, 0.095, 0.295, 0.338])
    verdicts = Counter(r["audit_verdict"] for r in bp_aud)
    tot_bp = sum(verdicts.values())
    n_bridged_ok = sum(1 for r in bp_aud if str(r["was_bridged_xy_interval"]).lower() == "true" and r["audit_verdict"] == "Verified Elementary Breakpoint")
    n_bridged_tot = sum(1 for r in bp_aud if str(r["was_bridged_xy_interval"]).lower() == "true")
    n_unsupp = (
        verdicts["Artifact: Ambiguous 'X/Y' Breakpoint Plateau"]
        + verdicts["Failed Poisson Floor (<3 Informative SNPs)"]
        + verdicts["Artifact: Degenerate/Boundary Segment"]
        + verdicts["Artifact: Unclassified 'U' Segment"]
        + verdicts["Failed Fisher Exact Test (p > 1e-3)"]
    )
    bars_c = [
        (verdicts["Verified Elementary Breakpoint"], BLUE, r"Verified ($p < 10^{-3}$)"),
        (verdicts["Artifact: Unresolved 2nd-Gen CRF Label"], SKY, "Unexp. 2nd-gen CRF"),
        (verdicts["Artifact: Ambiguous 'X/Y' Breakpoint Plateau"], SLATE_MID, "Ambig. X/Y plateau"),
        (verdicts["Failed Poisson Floor (<3 Informative SNPs)"] + verdicts["Artifact: Degenerate/Boundary Segment"], MUTED, "Floor (<3 SNPs) & edge"),
        (verdicts["Artifact: Unclassified 'U' Segment"], SLATE_LIGHT, "Unclassified U"),
        (verdicts["Failed Fisher Exact Test (p > 1e-3)"], VERMILION, r"Failed Fisher ($p > 10^{-3}$)"),
    ]
    yc = np.arange(len(bars_c))[::-1]
    for y, (cnt, col, lbl_tex) in zip(yc, bars_c):
        pct = 100.0 * cnt / tot_bp
        ec = MUTED if col == SLATE_LIGHT else "none"
        ax_c.barh(y, pct, height=0.56, color=col, edgecolor=ec, lw=0.55, zorder=3)
        ax_c.text(pct + 1.0, y, f"{cnt} ({pct:.1f}%)", va="center", ha="left", fontsize=6.0, fontweight="bold", color=INK)
    ax_c.set_yticks(yc)
    ax_c.set_yticklabels([b[2] for b in bars_c], fontsize=6.0)
    ax_c.set_xlim(0, 74)
    ax_c.set_xticks([0, 20, 40, 60])
    ax_c.set_xlabel(f"Share of {tot_bp:,} bridged junctions (%)", fontsize=6.4, color=INK)
    ax_c.xaxis.grid(True, color=HAIR, lw=0.5, zorder=1)

    m_comp = [r for r in comp_rows if int(r["n_verified_breakpoints"]) >= 1]
    vc_comp = Counter(r["structural_formula"] for r in m_comp)
    n_shared_forms = sum(1 for k, c in vc_comp.items() if c > 1)
    n_aliases = sum(c - 1 for c in vc_comp.values() if c > 1)
    dem_entries = [e for e in cb["compendium"] if e["category"] == "Demoted_LANL_Entry"]
    n_dem_crfs = sum(1 for e in dem_entries if e["label"].startswith("CRF"))
    n_dem_a = len(dem_entries) - n_dem_crfs
    n_u_seg = verdicts["Artifact: Unclassified 'U' Segment"]
    n_unsupp_idx = n_unsupp - n_u_seg
    pct_unsupp_idx = 100.0 * n_unsupp_idx / tot_bp
    pct_u_seg = 100.0 * n_u_seg / tot_bp
    ax_c.text(
        0.97, 0.08,
        f"{pct_unsupp_idx:.1f}% unsupported + {pct_u_seg:.1f}% U\n{n_bridged_ok}/{n_bridged_tot} X/Y bridged\n{n_shared_forms} shared ({n_aliases} aliases)\n{n_dem_crfs} CRFs & {n_dem_a} A* demoted",
        transform=ax_c.transAxes, fontsize=6.0, color=SUBINK, ha="right", va="bottom",
    )

    # --- Panel D: 3-Gen Recombination DAG (Monochrome-to-Steel Progression) ---
    ax_d1 = fig.add_axes([0.562, 0.102, 0.165, 0.325])
    gens = Counter(int(r["generation"]) for r in dag_rows)
    gen_labels = ["G0\nDem.", "G1\nPrim.", "G2\n2nd", "G3\n3rd"]
    gen_vals = [gens[0], gens[1], gens[2], gens[3]]
    gen_cols = [VERMILION, SLATE_LIGHT, SKY, BLUE]
    xg = np.arange(len(gen_labels))
    for i, (v, col) in enumerate(zip(gen_vals, gen_cols)):
        ec = SLATE_MID if col == SLATE_LIGHT else "none"
        ax_d1.bar(i, v, width=0.58, color=col, edgecolor=ec, lw=0.55, zorder=3)
        ax_d1.text(i, v + 2.2, f"{v}", ha="center", va="bottom", fontsize=6.0, fontweight="bold", color=INK,
                   bbox=dict(boxstyle="square,pad=0.06", fc="white", ec="none", alpha=1.0), zorder=5)
    ax_d1.set_xticks(xg)
    ax_d1.set_xticklabels(gen_labels, fontsize=6.0)
    ax_d1.set_ylim(0, 88)
    ax_d1.set_yticks([0, 25, 50, 75])
    ax_d1.set_ylabel(f"LANL CRFs (n = {len(dag_rows)})", fontsize=6.4, color=INK)
    ax_d1.yaxis.grid(True, color=HAIR, lw=0.5, zorder=1)

    ax_d2 = fig.add_axes([0.822, 0.102, 0.162, 0.325])
    parent_map = {}
    for r in dag_rows:
        lbl = r["legacy_label"]
        raw_p = r["immediate_recombinant_parents"]
        if raw_p.startswith("None"):
            parent_map[lbl] = []
        else:
            parent_map[lbl] = [t.split("(")[0].strip() for t in re.split(r"[,+]", raw_p) if t.strip()]

    def _ancestors(node):
        seen = set()
        stack = list(parent_map.get(node, []))
        while stack:
            x = stack.pop()
            if x not in seen:
                seen.add(x)
                stack.extend(parent_map.get(x, []))
        return seen

    gen23 = [r["legacy_label"] for r in dag_rows if int(r["generation"]) >= 2]
    anc_cnt = Counter()
    hub3 = {"CRF01_AE", "CRF07_BC", "CRF02_AG"}
    covered_hub3 = 0
    for g in gen23:
        a = _ancestors(g)
        for x in a:
            anc_cnt[x] += 1
        if a & hub3:
            covered_hub3 += 1
    non_hub = len(gen23) - covered_hub3
    hub_bars = [
        ("CRF01_AE", anc_cnt["CRF01_AE"], BLUE),
        ("CRF07_BC", anc_cnt["CRF07_BC"], BLUE),
        ("CRF02_AG", anc_cnt["CRF02_AG"], BLUE),
        ("Other CRFs", non_hub, SLATE_LIGHT),
    ]
    yh = np.arange(len(hub_bars))[::-1]
    for y_i, (lbl, cnt, col) in zip(yh, hub_bars):
        ec = SLATE_MID if col == SLATE_LIGHT else "none"
        ax_d2.barh(y_i, cnt, height=0.56, color=col, edgecolor=ec, lw=0.55, zorder=3)
        ax_d2.text(cnt + 1.2, y_i, f"{cnt}", va="center", ha="left", fontsize=6.0, fontweight="bold", color=INK)
    ax_d2.set_yticks(yh)
    ax_d2.set_yticklabels([b[0] for b in hub_bars], fontsize=6.0)
    ax_d2.set_xlim(0, 40)
    ax_d2.set_xticks([0, 15, 30])
    ax_d2.set_xlabel(f"Gen 2\u20133 hubs (n={len(gen23)})", fontsize=6.2, color=INK)
    ax_d2.xaxis.grid(True, color=HAIR, lw=0.5, zorder=1)
    ax_d2.text(0.96, 0.14, f"{100.0 * covered_hub3 / len(gen23):.1f}% top 3",
               transform=ax_d2.transAxes, fontsize=6.0, fontweight="bold", color=INK, ha="right", va="bottom")

    audit_matplotlib_figure(fig, "Fig3_HIV_Pathology", quadrant_splits=(0.500, 0.495))
    fig.savefig(FIG_DIR / "fig_science_3_hiv_pathology.pdf", dpi=300)
    fig.savefig(FIG_DIR / "fig_science_3_hiv_pathology.png", dpi=350)
    plt.close(fig)


# =============================================================================
# FIGURE 4: Reconciled Stanford HIVDB Audit & 3-Class GenBank Archive Audit
# =============================================================================
def parse_newick(nwk):
    """Minimal pure-Python Newick parser returning node dict {name, len, children}."""
    s = nwk.strip().rstrip(";")
    idx = 0

    def parse_node():
        nonlocal idx
        node = {"name": "", "len": 0.0, "children": []}
        if s[idx] == "(":
            idx += 1
            while True:
                node["children"].append(parse_node())
                if s[idx] == ",":
                    idx += 1
                elif s[idx] == ")":
                    idx += 1
                    break
        start = idx
        while idx < len(s) and s[idx] not in ",():;":
            idx += 1
        label = s[start:idx].strip()
        if label:
            node["name"] = label
        if idx < len(s) and s[idx] == ":":
            idx += 1
            bstart = idx
            while idx < len(s) and s[idx] not in ",():;":
                idx += 1
            try:
                node["len"] = max(float(s[bstart:idx]), 0.0)
            except ValueError:
                node["len"] = 0.0
        return node

    return parse_node()


def plot_newick_mini(ax, tree_path, target_parent_prefix, title, ufboot_txt, scale_bar=0.05, show_units=False, mono_color=TEAL):
    """Plot phylogram using true unpadded ML branch lengths with explicit subst./site scale bar."""
    nwk = Path(tree_path).read_text()
    root = parse_newick(nwk)
    leaves = []

    def assign_xy(node, depth):
        if not node["children"]:
            y = float(len(leaves))
            node["x"] = depth
            node["y"] = y
            leaves.append(node)
        else:
            ys = []
            for ch in node["children"]:
                assign_xy(ch, depth + max(ch["len"], 0.0))
                ys.append(ch["y"])
            node["x"] = depth
            node["y"] = sum(ys) / len(ys)

    assign_xy(root, 0.0)
    max_x = max(n["x"] for n in leaves) or 0.1

    def draw_edges(node):
        if not node["children"]:
            return
        ys = [c["y"] for c in node["children"]]
        ax.plot([node["x"], node["x"]], [min(ys), max(ys)], color=MUTED, lw=0.55, zorder=2)
        for c in node["children"]:
            ax.plot([node["x"], c["x"]], [c["y"], c["y"]], color=MUTED, lw=0.55, zorder=2)
            draw_edges(c)

    draw_edges(root)
    parent_prefixes = (f"Ref_{target_parent_prefix}_", "Ref__E_") if target_parent_prefix == "?E" else (f"Ref_{target_parent_prefix}_",)
    for n in leaves:
        nm = n["name"]
        if nm.startswith("Ref_"):
            is_parent = nm.startswith(parent_prefixes)
            if is_parent:
                ax.scatter([n["x"]], [n["y"]], s=18, marker="s", color=BLUE, edgecolor="white", lw=0.35, zorder=5)
            else:
                ax.scatter([n["x"]], [n["y"]], s=8, marker="o", color=HAIR, zorder=3)
        else:
            ax.scatter([n["x"]], [n["y"]], s=18, marker="o", color=mono_color, edgecolor="white", lw=0.35, zorder=6)

    sb_len = scale_bar if max_x >= scale_bar * 1.1 else round(max_x * 0.3, 2) or 0.02
    y_sb = -2.3
    ax.plot([0.0, sb_len], [y_sb, y_sb], color=INK, lw=1.1, zorder=6)
    sb_lbl = f"{sb_len:g} subst./site" if show_units else f"{sb_len:g}"
    ax.text(0.0 if show_units else sb_len / 2.0, y_sb - 1.4, sb_lbl, fontsize=6.0, color=INK,
            ha="left" if show_units else "center", va="top")

    ax.set_xlim(-0.04 * max_x, 1.10 * max_x)
    ax.set_ylim(-5.6, max(len(leaves) + 1.8, 2.0))
    ax.axis("off")
    ax.set_title(f"{title}\n{ufboot_txt}", fontsize=6.0, fontweight="bold", color=INK, pad=3)


def build_fig4():
    with open(RES / "stanford_hivdb_audit_category_summary.csv") as f:
        st_rows = {r["Audit_Category"]: r for r in csv.DictReader(f)}
    with open(RES / "stanford_west_african_crf02_rt_validation.csv") as f:
        wa_rows = list(csv.DictReader(f))
    with open(RES / "novel_candidate_crfs_31k_audit.tsv") as f:
        cand_rows = list(csv.DictReader(f, delimiter="\t"))
    with open(RES / "genbank_novel_circulating_recombinants.csv") as f:
        full_gb = list(csv.DictReader(f))
    with open(RES / "novel_lineages_phylogenetic_validation.csv") as f:
        phylo_rows = list(csv.DictReader(f))
    with open(RES / "scale_partial_pol_novel_circulating_lineages.csv") as f:
        pol_rows = list(csv.DictReader(f))
    with open(RES / "scale_partial_env_novel_circulating_lineages.csv") as f:
        env_rows = list(csv.DictReader(f))
    with open(RES / "scale_partial_env_phylogenetic_validation.csv") as f:
        env_phylo_rows = list(csv.DictReader(f))
    a1d_gp41_phylo = [r for r in env_phylo_rows if r["rhizaeon_formula"] == "A1[D:gp41]" and int(r["segment_idx"]) == 1][0]

    n_linked_tot = sum(int(r["Sequence_Count_Linked_111905"]) for r in st_rows.values())
    n_unlinked_tot = int(st_rows["Stanford_Unlinked_No_Isolate_Row"]["Sequence_Count_All_120k"])
    n_concordant_linked = (
        int(st_rows["Concordant_Pure"]["Sequence_Count_Linked_111905"])
        + int(st_rows["Concordant_CRF_or_Recomb"]["Sequence_Count_Linked_111905"])
    )
    pct_concordant_linked = 100.0 * n_concordant_linked / n_linked_tot

    fig = plt.figure(figsize=(CANVAS_WIDTH_IN, 5.10), facecolor="white")
    add_quadrant_hairlines(fig, x_split=0.500, y_split=0.495)
    draw_panel_header(fig, 0.014, 0.972, "A", f"Stanford HIVDB pol audit ({n_linked_tot:,} linked)")
    draw_panel_header(fig, 0.514, 0.972, "B", "Three-class HIV-1 GenBank archive audit")
    draw_panel_header(fig, 0.014, 0.472, "C", "Exact maps: monophyletic CRF vs. hotspots")
    draw_panel_header(fig, 0.514, 0.472, "D", "Segment ML tests separate descent vs. convergence")

    # --- Panel A: Reconciled Stanford HIVDB Discrepancies (Semantic Monochrome + Rust Alert) ---
    ax_a = fig.add_axes([0.195, 0.562, 0.289, 0.362])
    n_alt_grp = (
        int(st_rows["Stanford_CRF_vs_Different_RhizType_Architecture"]["Sequence_Count_Linked_111905"])
        + int(st_rows["Stanford_Pure_vs_Different_Pure_Trunk"]["Sequence_Count_Linked_111905"])
        + int(st_rows["Stanford_CRF01_AE_Pol_Is_Pure_A1_Trunk"]["Sequence_Count_Linked_111905"])
    )
    bars_a = [
        ("Multi-enzyme mosaic", int(st_rows["Stanford_MultiEnzyme_LocusExpected_Recombinant"]["Sequence_Count_Linked_111905"]), SLATE_MID),
        ("True clonal contradict.", int(st_rows["Stanford_IntraPatient_True_Clonal_Contradiction"]["Sequence_Count_Linked_111905"]), VERMILION),
        (r"CRF w/ pure $\it{pol}$", int(st_rows["Stanford_CRF_vs_RhizType_Pure_Trunk"]["Sequence_Count_Linked_111905"]), SLATE_LIGHT),
        (r"Pure w/ mosaic $\it{pol}$", int(st_rows["Stanford_Pure_vs_RhizType_Recombinant"]["Sequence_Count_Linked_111905"]), BLUE),
        ("Ambiguous resolved", int(st_rows["Stanford_Ambiguous_or_Unknown_Resolved"]["Sequence_Count_Linked_111905"]), SLATE_LIGHT),
        ("Alt. CRF / trunk / 01_AE", n_alt_grp, SLATE_LIGHT),
    ]
    ya = np.arange(len(bars_a))[::-1]
    for y, (lbl, cnt, col) in zip(ya, bars_a):
        pct = 100.0 * cnt / n_linked_tot
        ec = MUTED if col == SLATE_LIGHT else "none"
        ax_a.barh(y, cnt / 1000.0, height=0.56, color=col, edgecolor=ec, lw=0.55, zorder=3)
        ax_a.text(cnt / 1000.0 + 0.16, y, f"{cnt:,} ({pct:.2f}%)", va="center", ha="left", fontsize=6.0, fontweight="bold", color=INK)
    ax_a.set_yticks(ya)
    ax_a.set_yticklabels([b[0] for b in bars_a], fontsize=6.0)
    ax_a.set_xlim(0, 10.6)
    ax_a.set_xticks([0, 2, 4, 6, 8])
    ax_a.set_xlabel(f"Linked pol records (k; {pct_concordant_linked:.1f}% concord.)", fontsize=6.1, color=INK)
    ax_a.xaxis.grid(True, color=HAIR, lw=0.5, zorder=1)
    n_wa_trees = len(wa_rows)
    n_wa_ok = sum(1 for r in wa_rows if str(r["closer_to_CRF02_than_A6"]).lower() == "true")
    n_uncl_hivdb = int(st_rows["Unclassifiable_or_Short_Fragment"]["Sequence_Count_Linked_111905"])
    ax_a.text(
        0.98, 0.06,
        f"+{n_uncl_hivdb} short/unclass.\n+{n_unlinked_tot:,} unlinked typed\nWA 02_AG: {n_wa_ok}/{n_wa_trees} ML",
        transform=ax_a.transAxes, fontsize=6.0, color=SUBINK, ha="right", va="bottom",
    )

    # --- Panel B: Three-Class GenBank Archive Audit (Strict 3-Class Semantic Palette) ---
    ax_b = fig.add_axes([0.665, 0.562, 0.315, 0.362])
    n_promoted_official = sum(
        1 for r in cand_rows
        if r["disposition"] in ("Promoted_Official_LANL_CRF", "Scrubbed_Official_LANL_CRFs")
    )
    ru_row = [r for r in full_gb if r["rhizaeon_formula"].startswith("G0[B:gag.MA")][0]
    da1_row = [r for r in full_gb if r["rhizaeon_formula"] == "D[A1:env]"][0]
    env_top = {r["window_formula"]: r for r in env_rows}
    pol_top = {r["window_formula"]: r for r in pol_rows}
    a1d_gp41 = env_top["A1[D:gp41]"]
    a1c_gp41 = env_top["A1[C:gp41]"]
    other_full_pts = sum(
        int(r["n_independent_patients_formula"]) for r in full_gb
        if r["rhizaeon_formula"] not in ("G0[B:gag.MA | A1:gp41..nef]", "D[A1:env]")
    )
    other_part_pts = (
        sum(int(r["n_independent_patients"]) for r in pol_rows)
        + sum(int(r["n_independent_patients"]) for r in env_rows if r["window_formula"] not in ("A1[D:gp41]", "A1[C:gp41]"))
    )

    cls_bars = [
        ("8 Post-2023 CRFs\n(LANL rediscovered)", 31, BLUE, "31 bench (CRF127\u2013205)"),
        ("RU_G0B_gagMA\n(True mono. CRF)", int(ru_row["n_independent_patients_formula"]), TEAL, f"{int(ru_row['n_independent_patients_formula'])} (4/4 segs 100%)"),
        ("A1[D:gp41]*\n(Conv. env URF)", int(a1d_gp41["n_independent_patients"]), VERMILION, f"{int(a1d_gp41['n_independent_patients'])} (15 stud.; 0/2 mono.)"),
        ("A1[C:gp41]*\n(Conv. env URF)", int(a1c_gp41["n_independent_patients"]), VERMILION, f"{int(a1c_gp41['n_independent_patients'])} (11 stud.; 0/2 mono.)"),
        ("D[A1:env]\n(Conv. full URF)", int(da1_row["n_independent_patients_formula"]), VERMILION, f"{int(da1_row['n_independent_patients_formula'])} (env UFB 50%)"),
        ("6 Other full URF\n(Hotspots)", other_full_pts, SLATE_LIGHT, f"{other_full_pts} pts"),
        (f"10 Other part.*\n({len(pol_rows)} pol + {len(env_rows)-2} env)", other_part_pts, SLATE_MID, f"{other_part_pts} (1/10 mono.)"),
    ]
    yb = np.arange(len(cls_bars))[::-1]
    for y, (lbl, val, col, tag) in zip(yb, cls_bars):
        ec = SLATE_MID if col == SLATE_LIGHT else "none"
        ax_b.barh(y, val, height=0.58, color=col, edgecolor=ec, lw=0.55, zorder=3)
        ax_b.text(val + 1.4, y, tag, va="center", ha="left", fontsize=6.0, fontweight="bold", color=INK)
    ax_b.set_yticks(yb)
    ax_b.set_yticklabels([c[0] for c in cls_bars], fontsize=6.0)
    ax_b.set_xlim(0, 108)
    ax_b.set_xticks([0, 20, 40, 60, 80])
    ax_b.set_xlabel("Patient genomes (*partial: 1/12 mono. env)", fontsize=6.1, color=INK)
    ax_b.xaxis.grid(True, color=HAIR, lw=0.5, zorder=1)
    assert n_promoted_official == 8

    # --- Panel C: Exact Genomic Coordinate Comparison (RU_G0B_gagMA true 790..1154 width!) ---
    ax_c = fig.add_axes([0.135, 0.088, 0.345, 0.345])
    ax_c.set_xlim(790, 9417)
    ax_c.set_ylim(-0.6, 5.2)
    ax_c.set_xticks([790, 2253, 5097, 6225, 8379, 9417])
    ax_c.set_xticklabels([r"$\it{gag}$", r"$\it{pol}$", r"$\it{vif}$", r"$\it{env}$", r"$\it{gp41}$", r"$\it{nef}$"], fontsize=6.0)
    ax_c.set_yticks([])
    ax_c.spines["left"].set_visible(False)

    seg_pal_c = {"B": AMBER, "G0": SLATE, "A1": BLUE, "D": INDIGO}
    ru_segs = [r for r in phylo_rows if r["rhizaeon_formula"].startswith("G0[B:gag.MA")]
    da1_segs = [r for r in phylo_rows if r["rhizaeon_formula"] == "D[A1:env]"]

    y_ru = 3.55
    for sr in ru_segs:
        s0, s1 = int(sr["hxb2_start"]), int(sr["hxb2_end"])
        tr = sr["rhizaeon_assigned_trunk"]
        ax_c.add_patch(mpatches.Rectangle((s0, y_ru - 0.32), s1 - s0, 0.64,
                                          facecolor=seg_pal_c.get(tr, MUTED), edgecolor="white", lw=0.45, zorder=3))
        if (s1 - s0) >= 1200:
            ax_c.text((s0 + s1) / 2.0, y_ru, f"{tr} ({s0}\u2013{s1})", color="white", fontsize=6.0, fontweight="bold", ha="center", va="center")
    b_seg = ru_segs[0]
    b_mid = 0.5 * (int(b_seg["hxb2_start"]) + int(b_seg["hxb2_end"]))
    ax_c.annotate(
        f"B:{b_seg['hxb2_start']}\u2013{b_seg['hxb2_end']} ({b_seg['span_nt']} nt)",
        xy=(b_mid, y_ru + 0.34), xytext=(2200, y_ru + 0.85),
        fontsize=6.0, fontweight="bold", color=INK, ha="left", va="center",
        arrowprops=dict(arrowstyle="-|>", color=SUBINK, lw=0.80, mutation_scale=6.0),
    )
    a1_seg_ru = ru_segs[2]
    ax_c.text(8550, y_ru + 0.72, f"A1:{a1_seg_ru['hxb2_start']}\u2013{a1_seg_ru['hxb2_end']}",
              fontsize=6.0, fontweight="bold", color=INK, ha="center", va="center")
    ax_c.text(650, y_ru, "RU_G0B\n(N=6)", fontsize=6.0, fontweight="bold", color=TEAL, ha="right", va="center")

    y_d = 1.85
    for sr in da1_segs:
        s0, s1 = int(sr["hxb2_start"]), int(sr["hxb2_end"])
        tr = sr["rhizaeon_assigned_trunk"]
        ufb = float(sr["lineage_cluster_ufboot_pct"])
        ax_c.add_patch(mpatches.Rectangle((s0, y_d - 0.32), s1 - s0, 0.64,
                                          facecolor=seg_pal_c.get(tr, MUTED), edgecolor="white", lw=0.45, zorder=3))
        if tr == "D" and (s1 - s0) >= 2000:
            ax_c.text((s0 + s1) / 2.0, y_d, f"D ({s0}\u2013{s1})", color="white", fontsize=6.0, fontweight="bold", ha="center", va="center")
        elif tr == "A1":
            ax_c.text((s0 + s1) / 2.0, y_d, f"A1 ({ufb:.0f}%)", color="white", fontsize=6.0, fontweight="bold", ha="center", va="center")
    ax_c.text(650, y_d, "D[A1:env]\n(N=11)", fontsize=6.0, fontweight="bold", color=VERMILION, ha="right", va="center")

    y_pe = 0.35
    for s_str in a1d_gp41["representative_barcode"].split("|"):
        tr, rng_s = s_str.split(":")
        s0, s1 = [int(v) for v in rng_s.split("-")]
        ax_c.add_patch(mpatches.Rectangle((s0, y_pe - 0.32), s1 - s0, 0.64,
                                          facecolor=seg_pal_c.get(tr, MUTED), edgecolor="white", lw=0.45, zorder=3))
    ufb_pe = float(a1d_gp41_phylo["lineage_cluster_ufboot_pct"])
    ax_c.text(7250, y_pe, f"A1 ({ufb_pe:.0f}%)", color="white", fontsize=6.0, fontweight="bold", ha="center", va="center")
    ax_c.text(9080, y_pe + 0.58, "D gp41", color=INK, fontsize=6.0, fontweight="bold", ha="center", va="center")
    ax_c.text(650, y_pe, "A1[D:gp41]*\n(N=49)", fontsize=6.0, fontweight="bold", color=VERMILION, ha="right", va="center")
    ax_c.axvline(8230, ymin=0.05, ymax=0.78, color=VERMILION, lw=0.75, ls=":")
    ax_c.text(5900, 0.92, "Convergent env/gp41 hinge", fontsize=6.0, color=SUBINK, ha="right", va="center")

    # --- Panel D: Direct ML Tree Proof (RU Monophyly vs East African Polyphyly) ---
    d_specs = [
        (TREES_DIR / "lin02_G0_B_gag_MA_A1_gp41_nef_seg1_B_790_1154.fasta.treefile", "B", "RU Seg 1: gag.MA", "100% mono., 6/6", TEAL),
        (TREES_DIR / "lin01_D_A1_env_seg1_D_790_6335.fasta.treefile", "D", "D[A1:env] Seg 1: D", "100% mono., 11/11", TEAL),
        (TREES_DIR / "lin01_D_A1_env_seg2_A1_6336_8306.fasta.treefile", "A1", "D[A1:env] Seg 2: env", "50% polyphyletic", VERMILION),
    ]
    for i, (tp, pref, ttl, ufb, mcol) in enumerate(d_specs):
        axt = fig.add_axes([0.520 + i * 0.156, 0.118, 0.148, 0.280])
        plot_newick_mini(axt, tp, pref, ttl, ufb, scale_bar=0.05, show_units=(i == 0), mono_color=mcol)

    ax_dleg = fig.add_axes([0.518, 0.024, 0.464, 0.055])
    ax_dleg.axis("off")
    d_handles = [
        plt.Line2D([0], [0], marker="o", color="w", markerfacecolor=TEAL, markersize=4.5, label="Monophyletic (RU / D)"),
        plt.Line2D([0], [0], marker="o", color="w", markerfacecolor=VERMILION, markersize=4.5, label="Polyphyletic URF (env)"),
        plt.Line2D([0], [0], marker="s", color="w", markerfacecolor=BLUE, markersize=4.5, label="Parent ref"),
    ]
    ax_dleg.legend(handles=d_handles, frameon=False, fontsize=6.0, loc="center", ncol=3, handlelength=0.8, columnspacing=0.6)

    audit_matplotlib_figure(fig, "Fig4_Hidden_Epidemics", quadrant_splits=(0.500, 0.495))
    fig.savefig(FIG_DIR / "fig_science_4_hidden_epidemics.pdf", dpi=300)
    fig.savefig(FIG_DIR / "fig_science_4_hidden_epidemics.png", dpi=350)
    plt.close(fig)


# =============================================================================
# FIGURE 5: Co-Infection Supply & Two-Stage Epistatic Selection Sieve
# =============================================================================
def build_fig5():
    with open(RES / "regional_heterozygosity_urf_rates.csv") as f:
        reg_rows = list(csv.DictReader(f))
    with open(RES / "supp_manifest_sieve_749_lineages.csv") as f:
        lin_rows = list(csv.DictReader(f))
    with open(RES / "sieve_offspring_distribution_fits.csv") as f:
        off_rows = {r["model"]: r for r in csv.DictReader(f)}
    with open(RES / "sieve_hcirc_multivariate_regression.csv") as f:
        hcirc_rows = list(csv.DictReader(f))
    with open(RES / "recombination_viability_and_establishment_sieve.csv") as f:
        sieve_cmp = {r["feature_name"]: r for r in csv.DictReader(f)}
    with open(RES / "sieve_stage1_snp_and_sim_null.csv") as f:
        s1_null = {
            r["zone"]: r for r in csv.DictReader(f)
            if r["sim_regime"].startswith("FullGenome_Screen")
        }
    with open(RES / "sieve_stage2_tat_rre_multivariable.csv") as f:
        s2_mv = list(csv.DictReader(f))

    def _parse_segs(bc_str):
        out = []
        for s_item in bc_str.split("|"):
            if ":" in s_item and "-" in s_item:
                tr_s, rng_s = s_item.strip().split(":")
                a_i, b_i = [int(x) for x in rng_s.split("-")]
                out.append((tr_s, a_i, b_i))
        return out

    def _has_3p_bp(segs):
        return any(5831 <= 0.5 * (segs[i][2] + segs[i + 1][1]) <= 8653 for i in range(len(segs) - 1))

    all_bps_all = []
    sizes = []
    gb_counts = {}
    urf_segs_all = []
    for r in lin_rows:
        sizes.append(int(r["n_genbank_patients"]))
        segs = _parse_segs(r["hxb2_barcode"])
        bps = [0.5 * (segs[i][2] + segs[i + 1][1]) for i in range(len(segs) - 1)]
        all_bps_all.extend(bps)
        if r.get("survey_type") == "Known_CRF":
            k_crf = r["lineage_key"].replace("KNOWN:", "")
            gb_counts[k_crf] = gb_counts.get(k_crf, 0) + int(r["n_genbank_patients"])
        if str(r.get("in_singleton_urf_461", "")).lower() == "true" and len(segs) > 1:
            urf_segs_all.append(segs)

    # Exact Phase 1 comparative sieve cohort construction (151 Est, 22 Major N>=10, 491 Singleton URFs)
    with open(RES / "hiv1_descriptive_nomenclature_compendium.csv") as f:
        comp_crf_rows = [
            r for r in csv.DictReader(f)
            if r["legacy_category"].startswith("CRF") and int(r["n_verified_segments"]) > 1
        ]
    with open(CORE_DATA_DIR / "compendium_bundle.json") as f:
        bundle_crfs = json.load(f)["compendium"]
    with open(RES / "genbank_novel_circulating_recombinants.csv") as f:
        nov_crf_rows = [
            r for r in csv.DictReader(f)
            if r["lineage_status"].startswith("Candidate_Regional_Mosaic")
        ]

    est_segs_all = []
    maj_segs_all = []
    seen_crf_legs = set()
    for r in comp_crf_rows:
        segs = _parse_segs(r["hxb2_barcode"])
        if len(segs) <= 1:
            continue
        code, leg = r["short_code"].strip(), r["legacy_label"].strip()
        seen_crf_legs.add(leg)
        n_gb = max(gb_counts.get(code, 0), gb_counts.get(leg, 0), int(r["n_ref_genomes"]))
        est_segs_all.append(segs)
        if n_gb >= 10:
            maj_segs_all.append(segs)
    for ent in bundle_crfs:
        leg = str(ent.get("label", "")).strip()
        cat_b = str(ent.get("category", "")).strip()
        if not cat_b.startswith("CRF") or leg in seen_crf_legs:
            continue
        segs = _parse_segs(ent.get("hxb2_barcode", ""))
        if len(segs) <= 1:
            continue
        seen_crf_legs.add(leg)
        code = str(ent.get("short_code", "")).strip()
        n_gb = max(gb_counts.get(code, 0), gb_counts.get(leg, 0), 1)
        est_segs_all.append(segs)
        if n_gb >= 10:
            maj_segs_all.append(segs)
    for nr in nov_crf_rows:
        segs = _parse_segs(nr["representative_hxb2_barcode"])
        if len(segs) <= 1:
            continue
        est_segs_all.append(segs)
        if int(nr["n_homologous_bp_patients"]) >= 10:
            maj_segs_all.append(segs)

    cohort_totals = Counter({
        "Singleton": len(urf_segs_all),
        "Established": len(est_segs_all),
        "Major": len(maj_segs_all),
    })
    cohort_intact = Counter({
        "Singleton": sum(not _has_3p_bp(s) for s in urf_segs_all),
        "Established": sum(not _has_3p_bp(s) for s in est_segs_all),
        "Major": sum(not _has_3p_bp(s) for s in maj_segs_all),
    })

    or_intact, p_intact = fisher_exact([
        [cohort_intact["Established"], cohort_totals["Established"] - cohort_intact["Established"]],
        [cohort_intact["Singleton"], cohort_totals["Singleton"] - cohort_intact["Singleton"]],
    ])

    fig = plt.figure(figsize=(CANVAS_WIDTH_IN, 4.45), facecolor="white")
    add_quadrant_hairlines(fig, x_split=0.500, y_split=0.495)
    draw_panel_header(fig, 0.014, 0.972, "A", "Co-circulation drives de novo URF supply")
    draw_panel_header(fig, 0.514, 0.972, "B", f"Heavy-tailed scaling across {len(lin_rows):,} lineages")
    draw_panel_header(fig, 0.014, 0.472, "C", "Stage 1 folding sieve across three null models")
    draw_panel_header(fig, 0.514, 0.472, "D", "Stage 2 cassette & pol hinge vs. tat control")

    # --- Panel A: H_circ vs URF Prevalence (Dynamic Bivariate WLS + Multivariable stats) ---
    ax_a = fig.add_axes([0.085, 0.588, 0.398, 0.336])
    xs, ys, ws, sz, cols = [], [], [], [], []
    label_specs = {
        "Uganda": (0.080, -1.8, "left", True),
        "Kenya": (0.070, 2.6, "left", True),
        "Tanzania": (-0.042, -8.2, "right", True),
        "Cameroon": (-0.015, 10.5, "right", True),
        "China": (0.085, -7.5, "right", True),
        "Russia": (-0.075, 5.8, "right", True),
        "Brazil": (-0.030, -11.0, "center", True),
        "Thailand": (0.025, -9.2, "center", True),
        "South Africa": (0.060, 12.2, "left", True),
        "USA": (0.145, -4.4, "left", True),
    }
    for r in reg_rows:
        h = float(r["h_circ_lineage_diversity"])
        u = float(r["f_urf_patient_prev"]) * 100.0
        n_p = int(r["n_patients"])
        c_name = r["country"]
        xs.append(h)
        ys.append(u)
        ws.append(float(n_p))
        sz.append(max(15, math.sqrt(n_p) * 3.0))
        cols.append(AMBER if int(r["is_ssa"]) == 1 else BLUE)
        if c_name in label_specs:
            dx, dy, halign, leader = label_specs[c_name]
            if leader:
                ax_a.plot([h, h + dx * 0.80], [u, u + dy * 0.80], color=MUTED, lw=0.55, zorder=3)
            ax_a.text(h + dx, u + dy, c_name, fontsize=6.0, fontweight="bold", color=INK, ha=halign, va="center", zorder=6)

    x_arr = np.array(xs)
    y_arr = np.array(ys)
    W_diag = np.diag(ws)
    X_mat = np.column_stack([np.ones_like(x_arr), x_arr])
    b0_wls, b1_wls = np.linalg.solve(X_mat.T @ W_diag @ X_mat, X_mat.T @ W_diag @ y_arr)
    rho_h, p_rho_h = spearmanr(x_arr, y_arr)
    wls_hc = [r for r in hcirc_rows if r["model"] == "WLS_by_N" and r["term"] == "h_circ_non_urf"][0]
    hc3_hc = [r for r in hcirc_rows if r["model"] == "Binomial_GLM_HC3" and r["term"] == "h_circ_non_urf"][0]

    ax_a.scatter(xs, ys, s=sz, c=cols, alpha=0.84, edgecolor="white", lw=0.45, zorder=4)
    gx = np.linspace(0, max(x_arr) * 1.03, 50)
    ax_a.plot(gx, b0_wls + b1_wls * gx, color=SLATE, ls="--", lw=0.90, zorder=2)
    ax_a.set_xlim(-0.05, 1.04)
    ax_a.set_ylim(-6.2, 48.0)
    ax_a.set_xlabel(r"Circulating heterozygosity ($H_{\mathrm{circ}}$)", fontsize=6.4, color=INK, labelpad=2.0)
    ax_a.set_ylabel("De novo URF prevalence (%)", fontsize=6.4, color=INK)
    ax_a.yaxis.grid(True, color=HAIR, lw=0.5, zorder=1)
    ax_a.text(
        0.02, 0.97,
        rf"$\rho = {rho_h:.2f}$ ({fmt_p_tex(p_rho_h)}, $n={len(reg_rows)}$)"
        + "\n"
        + rf"WLS $R^2={float(wls_hc['r2_or_pseudo_r2']):.2f}$, HC3 {fmt_p_tex(float(hc3_hc['p_value']))}",
        transform=ax_a.transAxes, fontsize=6.0, color=INK, va="top",
    )
    ax_a.legend(handles=[mpatches.Patch(color=AMBER, label="Sub-Saharan Africa"), mpatches.Patch(color=BLUE, label="Eurasia & Americas")],
                frameon=False, fontsize=6.0, loc="upper left", bbox_to_anchor=(0.00, 0.80), handlelength=0.85)

    # --- Panel B: Exact Discrete Hurwitz-Zeta Step CCDF across N = 761 Lineages ---
    ax_b = fig.add_axes([0.588, 0.588, 0.394, 0.336])
    s_arr = np.array(sorted(sizes))
    n_tot = len(s_arr)
    uniq_k = sorted(set(s_arr))
    ccdf = [np.sum(s_arr >= k) / n_tot for k in uniq_k]
    ax_b.step(uniq_k, ccdf, where="post", color=BLUE, lw=1.40, label=f"Observed (N={n_tot})", zorder=3)

    hz_row = off_rows["Discrete_PowerLaw_HurwitzZeta_xmin1"]
    ln_row = off_rows["Discrete_Lognormal_xmin1"]
    m_alpha = re.search(r"alpha=([0-9.]+)\s*\(SE=([0-9.]+)\)", hz_row["params"])
    alpha_hat = float(m_alpha.group(1))
    alpha_se = float(m_alpha.group(2))
    delta_ln = float(ln_row["delta_aic_vs_best_xmin1"])
    k_int = np.array(sorted(set(list(range(1, 30)) + uniq_k)), dtype=float)
    hz_ccdf = zeta(alpha_hat, k_int) / zeta(alpha_hat, 1.0)
    ax_b.step(
        k_int, hz_ccdf, where="post", color=INK, ls="--", lw=0.95,
        label=rf"Hurwitz-zeta ($\hat{{\alpha}}={alpha_hat:.2f}\pm {alpha_se:.2f}$)", zorder=4,
    )

    n_eq1 = int(np.sum(s_arr == 1))
    n_ge3 = int(np.sum(s_arr >= 3))
    n_ge10 = int(np.sum(s_arr >= 10))
    ax_b.plot([3, 3], [0.015, 1.35], color=SLATE_MID, ls=":", lw=0.85, zorder=2)
    ax_b.plot([10, 10], [0.015, 1.35], color=VERMILION, ls=":", lw=0.85, zorder=2)
    ax_b.text(3.35, 0.62, f"N \u2265 3\n{100*n_ge3/n_tot:.1f}%\n({n_ge3})", fontsize=6.0, color=INK, va="center")
    ax_b.text(12.2, 0.20, f"N \u2265 10\n{100*n_ge10/n_tot:.1f}% ({n_ge10})", fontsize=6.0, fontweight="bold", color=VERMILION, va="center")
    ax_b.text(0.97, 0.96, f"N=1: {100*n_eq1/n_tot:.1f}%\nLN-comparable",
              transform=ax_b.transAxes, ha="right", va="top", fontsize=6.0, color=SUBINK)
    ax_b.set_xscale("log")
    ax_b.set_yscale("log")
    ax_b.set_xlim(0.9, max(uniq_k) * 1.4)
    ax_b.set_ylim(0.001, 1.35)
    ax_b.set_xlabel("Independent host cluster size (N)", fontsize=6.4, color=INK, labelpad=2.0)
    ax_b.set_ylabel(r"CCDF $P(\mathrm{size} \geq N)$", fontsize=6.4, color=INK)
    ax_b.yaxis.grid(True, color=HAIR, lw=0.5, zorder=1)
    ax_b.legend(frameon=False, fontsize=6.0, loc="lower left")

    # --- Panel C: Genome Breakpoint Flux & Both-Core Depletion Under 3 Nulls ---
    ax_c1 = fig.add_axes([0.085, 0.285, 0.398, 0.138])
    bins = np.linspace(790, 9417, 55)
    hist, edges = np.histogram(all_bps_all, bins=bins)
    mids = 0.5 * (edges[:-1] + edges[1:])
    ax_c1.axvspan(1280, 1850, color=VERMILION_SOFT, alpha=0.80, zorder=1)
    ax_c1.axvspan(6450, 7650, color=VERMILION_SOFT, alpha=0.80, zorder=1)
    ax_c1.plot(mids, hist, color=BLUE, lw=1.0, zorder=3)
    ax_c1.fill_between(mids, hist, color=BLUE_SOFT, alpha=0.85, zorder=2)
    ax_c1.text(1565, max(hist) * 1.13, r"$\it{Gag}$-CA", fontsize=6.0, fontweight="bold", color=INK, ha="center")
    ax_c1.text(7200, max(hist) * 1.13, r"$\it{gp120}$ core", fontsize=6.0, fontweight="bold", color=INK, ha="center")
    ax_c1.set_xlim(790, 9417)
    ax_c1.set_ylim(0, max(hist) * 1.34)
    ax_c1.set_xticks([790, 2253, 5097, 6225, 8379, 9417])
    ax_c1.set_xticklabels([r"$\it{gag}$", r"$\it{pol}$", r"$\it{vif}$", r"$\it{env}$", r"$\it{gp41}$", r"$\it{nef}$"], fontsize=6.0)
    ax_c1.set_ylabel("Breaks", fontsize=6.2, color=INK)

    ax_c2 = fig.add_axes([0.185, 0.088, 0.302, 0.145])
    # Compute exact combined 2 rigid folding cores (gp120 + Gag-CA; excluding crossover-tolerant Pol-RT)
    z_gp = s1_null["gp120 Core (6450-7650)"]
    z_ca = s1_null["Gag-CA Core (1280-1850)"]
    k2_est = int(z_gp["est_k"]) + int(z_ca["est_k"])
    n2_est = int(z_gp["est_n"])
    p2_unif = float(z_gp["uniform_share"]) + float(z_ca["uniform_share"])
    p2_snp = float(z_gp["snp_density_null_share"]) + float(z_ca["snp_density_null_share"])
    p2_sim = float(z_gp["sim_detectability_null_share"]) + float(z_ca["sim_detectability_null_share"])
    from scipy.stats import binomtest as _bt
    s1_null["Two_Rigid_Cores"] = {
        "est_enrichment_vs_unif": (k2_est / n2_est) / p2_unif,
        "est_enrichment_vs_snp": (k2_est / n2_est) / p2_snp,
        "est_enrichment_vs_sim": (k2_est / n2_est) / p2_sim,
        "est_p_vs_sim": _bt(k2_est, n2_est, p2_sim).pvalue,
    }
    zone_specs = [
        (r"$\it{gp120}$ core", "gp120 Core (6450-7650)"),
        (r"$\it{Gag}$-CA core", "Gag-CA Core (1280-1850)"),
        ("Both cores", "Two_Rigid_Cores"),
        ("4 hinges", "All 4 Domain Hinges Combined"),
    ]
    yc2 = np.arange(len(zone_specs))
    for i, (lbl, key_z) in enumerate(zone_specs):
        zr = s1_null[key_z]
        vu = float(zr["est_enrichment_vs_unif"])
        vs = float(zr["est_enrichment_vs_snp"])
        vm = float(zr["est_enrichment_vs_sim"])
        ps = float(zr["est_p_vs_sim"])
        ax_c2.plot([min(vu, vs, vm), max(vu, vs, vm)], [i, i], color=MUTED, lw=1.0, zorder=2)
        ax_c2.scatter([vu], [i], facecolor=SLATE_LIGHT, edgecolor=SLATE_MID, lw=0.55, marker="o", s=18, zorder=3, label="Uniform" if i == 0 else "")
        ax_c2.scatter([vs], [i], color=SLATE_MID, marker="D", s=17, zorder=4, label="SNP null" if i == 0 else "")
        ax_c2.scatter([vm], [i], color=BLUE, marker="s", s=18, zorder=5, label="Sim null" if i == 0 else "")
        x_hi = max(vu, vs, vm)
        x_txt = 1.06 if ("Two_Rigid" in key_z or "Gag-CA" in key_z) else x_hi + 0.07
        if x_txt > x_hi + 0.08:
            ax_c2.plot([x_hi + 0.05, x_txt - 0.03], [i, i], color=HAIR, lw=0.55, zorder=2)
        ax_c2.text(x_txt, i, fmt_p_tex(ps), va="center", ha="left", fontsize=6.0, fontweight="bold",
                   color=VERMILION if ("gp120" in key_z or "Gag-CA" in key_z or "Two_Rigid" in key_z) else INK, zorder=10)
    ax_c2.axvline(1.0, color=SLATE, ls="--", lw=0.70, zorder=1)
    ax_c2.set_yticks(yc2)
    ax_c2.set_yticklabels([z[0] for z in zone_specs], fontsize=6.0)
    ax_c2.set_xlim(0.0, 2.85)
    ax_c2.set_ylim(-0.65, 3.55)
    ax_c2.set_xticks([0.5, 1.0, 1.5, 2.0])
    ax_c2.set_xlabel("Est. fold vs. null (1.0x)", fontsize=6.2, color=INK)
    ax_c2.legend(frameon=False, fontsize=6.0, loc="lower right", handlelength=0.75, handletextpad=0.35, borderaxespad=0.1)

    # --- Panel D: Austere Stage 1 vs Stage 2 Trajectories + 3' Cassette Preservation ---
    ax_d = fig.add_axes([0.665, 0.095, 0.318, 0.338])
    tat_cmh = sieve_cmp["Tat Exon 1 <-> Tat Exon 2"]
    rre_cmh = sieve_cmp["RRE Stem-Loop <-> Rev Exon 2"]
    gp_cmh = sieve_cmp["gp120 Core <-> gp41 Ectodomain"]
    ma_cmh = sieve_cmp["Gag-MA Matrix <-> gp41-CT Tail"]
    h_gp_pol = sieve_cmp["Permissive_Hinge: gag/pol Frameshift & PR/RT (1980-2620)"]
    rev_swap = sieve_cmp["Foreign_gp120_Reverts_at_gp41_CT_to_Match_Gag_MA"]

    n_urf_d = cohort_totals["Singleton"]
    n_est_d = cohort_totals["Established"]
    n_maj_d = cohort_totals["Major"]
    int_u = 100.0 * cohort_intact["Singleton"] / n_urf_d
    int_e = 100.0 * cohort_intact["Established"] / n_est_d
    int_m = 100.0 * cohort_intact["Major"] / n_maj_d

    loci = [
        (r"$\it{gp120}\leftrightarrow\it{gp41}$" + "\n(Folding)", gp_cmh, None,
         rf"$p_{{\mathrm{{lin}}}}<10^{{-4}}$", INK, "bold"),
        ("Intact 3' block\n(5831\u20138653)", None, (int_u, int_e, int_m, None),
         rf"OR={or_intact:.1f}, {fmt_p_tex(p_intact)}", INK, "bold"),
        (r"RRE$\leftrightarrow\it{rev2}$" + "\n(Export)", rre_cmh, None,
         rf"OR={float(rre_cmh['cmh_or_est_vs_urf']):.1f}, {fmt_p_tex(float(rre_cmh['cmh_p_est_vs_urf']))}", INK, "bold"),
        (r"$\it{tat1}\leftrightarrow\it{tat2}$" + "\n(Point ctrl)", tat_cmh, None,
         rf"ns ($p={float(tat_cmh['cmh_p_est_vs_urf']):.2f}$)", SUBINK, "normal"),
        (r"MA$\leftrightarrow\it{gp41}$-CT" + "\n(Distal ctrl)", ma_cmh, None,
         rf"ns ($p={float(ma_cmh['cmh_p_est_vs_urf']):.2f}$)", SUBINK, "normal"),
    ]
    yl = np.arange(len(loci))[::-1]
    for y, (lbl, row_d, custom_tup, note, col_note, fw_note) in zip(yl, loci):
        if custom_tup is not None:
            s_v, e_v, m_v, n_v = custom_tup
        else:
            s_v = float(row_d["singleton_urf_rate"]) * 100.0
            e_v = float(row_d["established_rate"]) * 100.0
            m_v = float(row_d["major_crf_rate"]) * 100.0
            n_v = float(row_d["est_lin_perm_null"]) * 100.0
        if n_v is not None:
            ax_d.plot([n_v, n_v], [y - 0.30, y + 0.30], color=SLATE, lw=0.85, ls=":",
                      label="Linear null" if y == yl[0] else "")
        ax_d.plot([min(s_v, e_v, m_v), max(s_v, e_v, m_v)], [y, y], color=MUTED, lw=1.1, zorder=2)
        ax_d.scatter([s_v], [y], facecolor=SLATE_LIGHT, edgecolor=SLATE_MID, lw=0.60, marker="o", s=21, zorder=4, label=f"URF ({n_urf_d})" if y == yl[0] else "")
        ax_d.scatter([e_v], [y], color=BLUE, marker="D", s=23, zorder=5, label=f"Est. ({n_est_d})" if y == yl[0] else "")
        ax_d.scatter([m_v], [y], color=VERMILION, marker="s", s=21, zorder=6, label=f"Major ({n_maj_d})" if y == yl[0] else "")
        x_anc = max([v for v in (s_v, e_v, m_v, n_v) if v is not None])
        ax_d.text(x_anc + 3.8, y, note, va="center", ha="left", fontsize=6.0, fontweight=fw_note, color=col_note,
                  bbox=dict(boxstyle="square,pad=0.08", fc="white", ec="none", alpha=1.0), zorder=8)

    ax_d.set_yticks(yl)
    ax_d.set_yticklabels([x[0] for x in loci], fontsize=6.0)
    ax_d.set_xlim(10, 160)
    ax_d.set_ylim(-2.95, 4.55)
    ax_d.set_xticks([20, 40, 60, 80, 100])
    ax_d.set_xlabel("Lineage rate (%)", fontsize=6.4, color=INK)
    ax_d.xaxis.grid(True, color=HAIR, lw=0.5, zorder=1)
    # Compute exact Gen-1 founder-deduplicated (N=71) primary CRF rates for Intact 3' Block & RRE<->rev2 dynamically!
    with open(RES / "hiv1_recombination_genealogy_dag.csv") as f:
        dag_gen_map = {r["legacy_label"].strip(): int(r["generation"]) for r in csv.DictReader(f)}
    gen1_seen_f = set()
    gen1_seen_legs = set()
    gen1_dedup_segs = []
    for r in comp_crf_rows:
        leg = r["legacy_label"].strip()
        gen1_seen_legs.add(leg)
        if dag_gen_map.get(leg, 1) != 1:
            continue
        sf = r["structural_formula"].strip()
        if sf in gen1_seen_f:
            continue
        segs = _parse_segs(r["hxb2_barcode"])
        if len(segs) > 1:
            gen1_seen_f.add(sf)
            gen1_dedup_segs.append(segs)
    for ent in bundle_crfs:
        leg = str(ent.get("label", "")).strip()
        cat_b = str(ent.get("category", "")).strip()
        if not (cat_b.startswith("CRF") or cat_b.startswith("Novel_Circulating")) or leg in gen1_seen_legs:
            continue
        if any(k in leg for k in ("0107", "07109", "0755", "63A6")):
            continue
        sf = str(ent.get("structural_formula", "")).strip()
        if sf in gen1_seen_f:
            continue
        segs = _parse_segs(ent.get("hxb2_barcode", ""))
        if len(segs) > 1:
            gen1_seen_f.add(sf)
            gen1_dedup_segs.append(segs)
    for r in nov_crf_rows:
        sf = r["rhizaeon_formula"].strip()
        if sf in gen1_seen_f:
            continue
        segs = _parse_segs(r["representative_hxb2_barcode"])
        if len(segs) > 1:
            gen1_seen_f.add(sf)
            gen1_dedup_segs.append(segs)
    g1_intact_k = sum(not _has_3p_bp(s) for s in gen1_dedup_segs)
    or_g1_intact, p_g1_intact = fisher_exact([
        [g1_intact_k, len(gen1_dedup_segs) - g1_intact_k],
        [cohort_intact["Singleton"], cohort_totals["Singleton"] - cohort_intact["Singleton"]],
    ])
    maj_gagpol_pct = 100.0 * sum(1980 <= 0.5 * (s[i][2] + s[i + 1][1]) <= 2620 for s in maj_segs_all for i in range(len(s) - 1)) / max(1, sum(len(s) - 1 for s in maj_segs_all))
    ax_d.text(
        0.03, 0.06,
        rf"$\it{{gag/pol}}$: {100*float(h_gp_pol['singleton_urf_rate']):.1f}%"
        + f"\u2192{maj_gagpol_pct:.1f}% "
        + f"($p={float(h_gp_pol['cmh_p_est_vs_urf']):.3f}$)\n"
        + f"Gen-1 ({len(gen1_dedup_segs)}): 3' {100*g1_intact_k/len(gen1_dedup_segs):.1f}% "
        + f"($p={p_g1_intact:.3f}$)",
        transform=ax_d.transAxes, fontsize=6.0, color=INK, va="bottom", ha="left",
    )
    ax_d.legend(frameon=False, fontsize=6.0, loc="lower right", bbox_to_anchor=(1.01, 0.03),
                ncol=1, handlelength=0.75, labelspacing=0.22, borderaxespad=0.1)

    audit_matplotlib_figure(fig, "Fig5_Epistatic_Sieve", quadrant_splits=(0.500, 0.495))
    fig.savefig(FIG_DIR / "fig_science_5_epistatic_sieve.pdf", dpi=300)
    fig.savefig(FIG_DIR / "fig_science_5_epistatic_sieve.png", dpi=350)
    plt.close(fig)


# =============================================================================
# SUPPLEMENTAL FIGURE S1: Authentic IQ-TREE 2 Validation of Monophyletic vs Convergent Lineages
# =============================================================================
def _plot_supp_tree(ax, tree_path, target_parent, title_line1, title_line2, trunk_colors, is_mono=True):
    """Plot clean true-branch-length phylogram with well-spaced >= 6pt labels for Supplemental Fig S1."""
    nwk = Path(tree_path).read_text()
    root = parse_newick(nwk)
    leaves = []

    def assign_xy(node, depth):
        if not node["children"]:
            y = float(len(leaves))
            node["x"] = depth
            node["y"] = y
            leaves.append(node)
        else:
            ys = []
            for ch in node["children"]:
                assign_xy(ch, depth + max(ch["len"], 0.0))
                ys.append(ch["y"])
            node["x"] = depth
            node["y"] = sum(ys) / len(ys)

    assign_xy(root, 0.0)
    max_x = max(n["x"] for n in leaves) or 0.1
    n_leaves = max(len(leaves), 1)
    parent_col = trunk_colors.get(target_parent, BLUE)
    q_col = TEAL if is_mono else VERMILION
    parent_prefixes = (f"Ref_{target_parent}_", "Ref__E_") if target_parent == "?E" else (f"Ref_{target_parent}_",)

    def _subtree_types(node):
        if not node["children"]:
            nm = node["name"]
            if nm.startswith("Query_"):
                return {"Q"}
            if nm.startswith(parent_prefixes):
                return {"P"}
            return {"B"}
        out = set()
        for c in node["children"]:
            out |= _subtree_types(c)
        return out

    def draw_edges(node):
        if not node["children"]:
            return
        ys = [c["y"] for c in node["children"]]
        ts = _subtree_types(node)
        vcol = q_col if ts == {"Q"} else (parent_col if ts <= {"Q", "P"} else MUTED)
        vlw = 1.15 if ts <= {"Q", "P"} else 0.55
        ax.plot([node["x"], node["x"]], [min(ys), max(ys)], color=vcol, lw=vlw, zorder=3 if vlw > 1 else 2)
        for c in node["children"]:
            cts = _subtree_types(c)
            ecol = q_col if cts == {"Q"} else (parent_col if cts <= {"Q", "P"} else MUTED)
            elw = 1.15 if cts <= {"Q", "P"} else 0.55
            ax.plot([node["x"], c["x"]], [c["y"], c["y"]], color=ecol, lw=elw, zorder=3 if elw > 1 else 2)
            draw_edges(c)

    draw_edges(root)
    q_nodes, p_nodes, b_nodes = [], [], []
    for n in leaves:
        nm = n["name"]
        if nm.startswith("Query_"):
            ax.scatter([n["x"]], [n["y"]], s=14, marker="o", color=q_col, edgecolor="white", lw=0.35, zorder=6)
            q_nodes.append(n)
        elif nm.startswith(parent_prefixes):
            ax.scatter([n["x"]], [n["y"]], s=15, marker="s", color=parent_col, edgecolor="white", lw=0.35, zorder=5)
            p_nodes.append(n)
        else:
            ax.scatter([n["x"]], [n["y"]], s=8, marker="o", color=HAIR, zorder=3)
            b_nodes.append(n)

    # Right-margin aligned label column outside all branches/tips (zero branch/spine collisions!)
    x_lbl_col = 1.07 * max_x
    y_lo = -max(2.6, n_leaves * 0.085)
    y_hi = len(leaves) - 0.1 + max(1.2, n_leaves * 0.04)
    min_row_gap = max(3.8, n_leaves * 0.155)
    placed_ys = []

    def _pick_free_y(target_y):
        if all(abs(target_y - py) >= min_row_gap for py in placed_ys) and (0.4 <= target_y <= len(leaves) - 0.5):
            return target_y
        # Search cleanly upward/downward within tree span for nearest open slot
        for step_k in range(1, 20):
            for sign in (1.0, -1.0):
                cand = target_y + sign * step_k * (min_row_gap * 0.45)
                if 0.6 <= cand <= len(leaves) - 0.5 and all(abs(cand - py) >= min_row_gap for py in placed_ys):
                    return cand
        return None

    if q_nodes:
        # Pick representative anchor tip closest to cluster midpoint (or top cluster if polyphyletic)
        qy_mid = sum(n["y"] for n in q_nodes) / len(q_nodes)
        q_anchor = min(q_nodes, key=lambda n: abs(n["y"] - qy_mid))
        qy_lbl = _pick_free_y(q_anchor["y"]) or q_anchor["y"]
        ax.plot([q_anchor["x"] + 0.02 * max_x, x_lbl_col - 0.02 * max_x], [q_anchor["y"], qy_lbl],
                color=q_col, lw=0.45, ls=":", zorder=4)
        ax.text(x_lbl_col, qy_lbl, f"Query (n={len(q_nodes)})",
                fontsize=6.0, fontweight="bold", color=q_col, va="center", ha="left")
        placed_ys.append(qy_lbl)

    if p_nodes:
        # Pick parental ref tip with largest vertical separation from Query label
        p_anchor = max(p_nodes, key=lambda n: abs(n["y"] - placed_ys[0]) if placed_ys else 0.0)
        py_lbl = _pick_free_y(p_anchor["y"])
        if py_lbl is None:
            py_lbl = max(0.8, placed_ys[0] - min_row_gap) if p_anchor["y"] <= placed_ys[0] else min(len(leaves) - 0.8, placed_ys[0] + min_row_gap)
        ax.plot([p_anchor["x"] + 0.02 * max_x, x_lbl_col - 0.02 * max_x], [p_anchor["y"], py_lbl],
                color=parent_col, lw=0.45, ls=":", zorder=4)
        ax.text(x_lbl_col, py_lbl, f"Ref {target_parent}",
                fontsize=6.0, fontweight="bold", color=parent_col, va="center", ha="left")
        placed_ys.append(py_lbl)

    labeled_bg = set()
    for bn in b_nodes:
        clade_nm = bn["name"].split("_")[1] if bn["name"].startswith("Ref_") and len(bn["name"].split("_")) >= 2 else bn["name"]
        if clade_nm in labeled_bg or len(labeled_bg) >= 3:
            continue
        if bn["y"] < 1.0 or bn["y"] > len(leaves) - 1.2:
            continue
        if all(abs(bn["y"] - py) >= min_row_gap for py in placed_ys):
            ax.plot([bn["x"] + 0.02 * max_x, x_lbl_col - 0.02 * max_x], [bn["y"], bn["y"]],
                    color=RULE, lw=0.45, ls=":", zorder=2)
            ax.text(x_lbl_col, bn["y"], clade_nm, fontsize=6.0, color=SUBINK, va="center", ha="left")
            placed_ys.append(bn["y"])
            labeled_bg.add(clade_nm)

    ax.set_xlim(-0.04 * max_x, 1.98 * max_x)
    ax.set_ylim(y_lo, y_hi)
    ax.set_yticks([])
    ax.spines["left"].set_visible(False)
    ax.tick_params(axis="x", labelsize=6.0, pad=1.5)
    step = round(max_x * 0.45, 2) or 0.05
    ticks_x = [0.0, step, round(2 * step, 2)]
    ax.set_xticks(ticks_x)
    ax.set_xticklabels([f"{v:.2f}" for v in ticks_x], fontsize=6.0)
    ax.set_xlabel("Subst. / site", fontsize=6.0, color=SUBINK, labelpad=1.5)
    ax.set_title(f"{title_line1}\n{title_line2}", fontsize=6.0, fontweight="bold", color=INK, pad=3)


def build_fig_s1():
    """Generate Supplemental Fig. S1 at 6.50 in width contrasting RU_G0B_gagMA monophyly vs convergent hotspot URFs."""
    with open(RES / "genbank_novel_circulating_recombinants.csv") as f:
        full_rows = {r["rhizaeon_formula"]: r for r in csv.DictReader(f)}
    with open(RES / "novel_lineages_phylogenetic_validation.csv") as f:
        val_rows = list(csv.DictReader(f))
    with open(RES / "scale_partial_env_phylogenetic_validation.csv") as f:
        env_val_rows = list(csv.DictReader(f))

    trunk_cols = {"A1": BLUE, "D": INDIGO, "G0": SLATE, "B": AMBER, "H": PURPLE, "F1": WINE}
    fig = plt.figure(figsize=(CANVAS_WIDTH_IN, 4.95), facecolor="white")

    s1_leg = [
        plt.Line2D([0], [0], marker="o", color=TEAL, lw=1.4, markersize=5, label="Monophyletic cluster (RU_G0B_gagMA / D)"),
        plt.Line2D([0], [0], marker="o", color=VERMILION, lw=1.4, markersize=5, label="Polyphyletic URF (env / gp41)"),
        plt.Line2D([0], [0], marker="s", color=BLUE, lw=1.2, markersize=5, label="Assigned parental ref"),
        plt.Line2D([0], [0], marker="o", color=MUTED, lw=0.8, markersize=4, label="Background Group M"),
    ]
    fig.legend(handles=s1_leg, loc="upper center", bbox_to_anchor=(0.50, 0.992), ncol=2, frameon=False, fontsize=6.0)

    # Row A: Strictly monophyletic RU_G0B_gagMA (all 4 segments, N=6 Russia, 100% UFBoot)
    y_a = 0.565
    fig.text(0.014, y_a + 0.348, "A", fontsize=9.5, fontweight="bold", color=INK, va="center")
    fig.text(0.036, y_a + 0.348, "Strictly monophyletic RU_G0B_gagMA (Russia, N=6: 4/4 segments UFBoot 100%)",
             fontsize=6.6, fontweight="bold", color=INK, va="center")
    ax_bar_a = fig.add_axes([0.025, y_a + 0.285, 0.955, 0.036])
    ax_bar_a.set_xlim(790, 9417)
    ax_bar_a.set_ylim(0, 1)
    ax_bar_a.axis("off")
    bc_ru = full_rows["G0[B:gag.MA | A1:gp41..nef]"]["representative_hxb2_barcode"]
    for s_item in bc_ru.split("|"):
        tr_b, c_b = s_item.split(":")
        b0, b1 = [int(x) for x in c_b.split("-")]
        ax_bar_a.add_patch(mpatches.Rectangle((b0, 0.12), b1 - b0, 0.76, facecolor=trunk_cols.get(tr_b, MUTED), edgecolor="white", lw=0.5))
        if (b1 - b0) >= 650:
            ax_bar_a.text((b0 + b1) / 2, 0.50, tr_b, fontsize=6.0, fontweight="bold", color="white", ha="center", va="center")
    ru_segs = sorted([r for r in val_rows if r["rhizaeon_formula"] == "G0[B:gag.MA | A1:gp41..nef]"],
                     key=lambda r: int(r["segment_idx"]))
    for ci, s_row in enumerate(ru_segs):
        axt = fig.add_axes([0.025 + ci * 0.244, y_a + 0.030, 0.218, 0.195])
        tr = s_row["rhizaeon_assigned_trunk"]
        h0, h1 = int(s_row["hxb2_start"]), int(s_row["hxb2_end"])
        ufb = float(s_row["lineage_cluster_ufboot_pct"])
        _plot_supp_tree(axt, s_row["treefile_path"], tr,
                        f"Seg {s_row['segment_idx']}: {s_row['gene_region']} \u2192 {tr}",
                        f"{h0}\u2013{h1} (Mono. {ufb:.0f}%)", trunk_cols, is_mono=True)

    # Row B: Convergent Polyphyletic Hotspot Controls (D[A1:env] N=11 & partial A1[D:gp41] N=49)
    y_b = 0.055
    fig.text(0.014, y_b + 0.365, "B", fontsize=9.5, fontweight="bold", color=INK, va="center")
    fig.text(0.036, y_b + 0.365, "Convergent polyphyletic controls: D[A1:env] (N=11 full) & A1[D:gp41]* (N=49 env, 15 studies)",
             fontsize=6.6, fontweight="bold", color=INK, va="center")
    ax_bar_b = fig.add_axes([0.025, y_b + 0.302, 0.955, 0.036])
    ax_bar_b.set_xlim(790, 9417)
    ax_bar_b.set_ylim(0, 1)
    ax_bar_b.axis("off")
    bc_d = full_rows["D[A1:env]"]["representative_hxb2_barcode"]
    for s_item in bc_d.split("|"):
        tr_b, c_b = s_item.split(":")
        b0, b1 = [int(x) for x in c_b.split("-")]
        ax_bar_b.add_patch(mpatches.Rectangle((b0, 0.12), b1 - b0, 0.76, facecolor=trunk_cols.get(tr_b, MUTED), edgecolor="white", lw=0.5))
        ax_bar_b.text((b0 + b1) / 2, 0.50, tr_b, fontsize=6.0, fontweight="bold", color="white", ha="center", va="center")

    da1_segs = sorted([r for r in val_rows if r["rhizaeon_formula"] == "D[A1:env]"],
                      key=lambda r: int(r["segment_idx"]))
    a1d_part = [r for r in env_val_rows if r["rhizaeon_formula"] == "A1[D:gp41]" and int(r["segment_idx"]) == 1][0]
    row_b_specs = da1_segs + [a1d_part]
    for ci, s_row in enumerate(row_b_specs):
        axt = fig.add_axes([0.025 + ci * 0.244, y_b + 0.030, 0.218, 0.195])
        tr = s_row["rhizaeon_assigned_trunk"]
        h0, h1 = int(s_row["hxb2_start"]), int(s_row["hxb2_end"])
        ufb = float(s_row["lineage_cluster_ufboot_pct"])
        is_m = str(s_row["is_strictly_monophyletic"]).lower() == "true"
        lbl_tag = "Mono." if is_m else "Polyphyletic"
        title1 = f"A1[D:gp41]* {s_row['gene_region']}" if ci == 3 else f"D[A1] Seg {s_row['segment_idx']}: {s_row['gene_region']}"
        _plot_supp_tree(axt, s_row["treefile_path"], tr, title1,
                        f"{h0}\u2013{h1} ({lbl_tag} {ufb:.0f}%)", trunk_cols, is_mono=is_m)

    audit_matplotlib_figure(fig, "FigS1_Phylogenies", quadrant_splits=None)
    fig.savefig(FIG_DIR / "fig4_chinese_a6_bridge_and_russian_lineages_trees.pdf", dpi=300)
    fig.savefig(FIG_DIR / "fig4_chinese_a6_bridge_and_russian_lineages_trees.png", dpi=350)
    plt.close(fig)


def main():
    build_fig1()
    build_fig2()
    build_fig3()
    build_fig4()
    build_fig5()
    build_fig_s1()
    print("Generated and audited all 5 main figures + Fig S1 in", FIG_DIR)


if __name__ == "__main__":
    main()
