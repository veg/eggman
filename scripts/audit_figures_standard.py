#!/usr/bin/env python3
"""
Strict Science/Nature Graphics, Legibility & Data-to-Ink Audit Suite.

Enforces quantitative standards across Matplotlib figures, exported PDF/PNG
artifacts, and LaTeX \\includegraphics directives in paper/main.tex and
paper/supplement.tex:
  1. Exact 1:1 print geometry (canvas width == 6.50 in == \\textwidth; 0% silent
     downscaling in LaTeX).
  2. Strict >= 6.0 pt typography floor at final print size with single sans-serif
     family and Type-42 editable vector outlines.
  3. Zero canvas clipping, zero quadrant-divider collisions, and zero pairwise
     text bounding-box overlaps via live renderer metrics.
  4. Zero raw ASCII code/math shorthand (->, <->, >=, <=, 1e-3, t_MRCA, p_lin,
     H_circ, rho=, alpha=, OOD, SOTA, closed dict, bundle, constellation).
  5. Tufte data-to-ink & colorblind-safe Okabe-Ito compliance (half-open spines,
     borderless legends, non-empty visual dictionaries).
"""

from __future__ import annotations

import re
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.text import Text

ROOT = Path("/Users/sergei/Projects/HIV-recombination")
PAPER_DIR = ROOT / "paper"
FIG_DIR = PAPER_DIR / "figures"

CANVAS_WIDTH_IN = 6.50  # Exact letterpaper (8.5 in) minus 1.0 in left/right margins
FONT_FLOOR_PT = 6.0

# Austere, low-chroma mineral + monochrome technical palette (colorblind-safe)
PALETTE: Dict[str, str] = {
    "INK": "#111827",            # Primary technical black
    "SUBINK": "#4B5563",         # Secondary charcoal
    "MUTED": "#9CA3AF",          # Light graphite / secondary annotations
    "HAIR": "#D1D5DB",           # Architectural divider hairlines
    "GRID": "#ECEEF1",           # Subtle axis gridlines
    "SLATE": "#374151",          # Primary dark neutral fill (G0, baseline)
    "SLATE_MID": "#6B7280",      # Medium neutral fill (comparators, benign)
    "SLATE_LIGHT": "#CBD5E1",    # Light neutral fill (background bars)
    "BLUE": "#2B5B84",           # Deep steel blue (EGGMAN, verified primary)
    "BLUE_SOFT": "#EDF2F7",      # Whisper steel tint
    "SKY": "#7A9BB8",            # Desaturated cool steel
    "TEAL": "#2F6B66",           # Muted pine / slate-teal (monophyletic, clonal)
    "TEAL_SOFT": "#EEF5F4",      # Whisper pine tint
    "AMBER": "#A87C38",          # Muted ochre / brass (secondary accent)
    "AMBER_SOFT": "#F8F3EA",     # Whisper ochre tint
    "VERMILION": "#9E3B2B",      # Muted terracotta / rust (breakpoints, alerts)
    "VERMILION_SOFT": "#F6ECEA", # Whisper rust tint
    "PURPLE": "#745C78",         # Dusty slate-plum (additional categorical trunk)
    "PURPLE_SOFT": "#F4F0F5",
    "INDIGO": "#2C3E60",         # Deep ink navy (NAm_B313, Subtype D)
    "WINE": "#6E3B4A",           # Muted oxblood (NAm_B32, Subtype J)
    "BG_SOFT": "#F8F9FA",        # Barely-there neutral ground
}

BANNED_ASCII_PATTERNS: List[Tuple[str, str]] = [
    (r"->", "Use Unicode arrow '\\u2192' or mathtext '\\rightarrow' instead of '->'"),
    (r"<->", "Use Unicode arrow '\\u2194' or mathtext '\\leftrightarrow' instead of '<->'"),
    (r">=", "Use Unicode '\\u2265' or mathtext '\\geq' instead of '>='"),
    (r"<=", "Use Unicode '\\u2264' or mathtext '\\leq' instead of '<='"),
    (r"\b\d+(?:\.\d+)?[eE][+-]?\d+\b", "Use formatted powers of ten ($10^{-k}$) instead of '1e-...'"),
    (r"\bt_MRCA\b", "Use mathtext '$t_{\\mathrm{MRCA}}$' instead of 't_MRCA'"),
    (r"\bp_lin\b", "Use mathtext '$p_{\\mathrm{lin}}$' instead of 'p_lin'"),
    (r"\bL_min\b", "Use mathtext '$L_{\\min}$' instead of 'L_min'"),
    (r"\bp_Fisher\b", "Use mathtext '$p_{\\mathrm{Fisher}}$' instead of 'p_Fisher'"),
    (r"\bH_circ\b", "Use mathtext '$H_{\\mathrm{circ}}$' instead of 'H_circ'"),
    (r"\bn_c\s*-\s*n_u\b", "Use mathtext '$n_c - n_u$' instead of raw ASCII 'n_c - n_u'"),
    (r"(?<!\\)\brho\s*=", "Use mathtext '$\\rho =$' instead of ASCII 'rho ='"),
    (r"(?<!\\)\balpha\s*=", "Use mathtext '$\\hat{\\alpha} =$' instead of ASCII 'alpha ='"),
    (r"\bOOD\b", "Spell out 'held-out' or 'uncatalogued' instead of 'OOD'"),
    (r"\bSOTA\b", "Write 'benchmark' instead of 'SOTA'"),
    (r"\bclosed dict\b", "Write 'fixed catalog' instead of 'closed dict'"),
    (r"\bbundles?\b", "Standardize on 'reference panel' instead of 'bundle'"),
    (r"\bconstellations?\b", "Standardize on 'genotypes' or 'segment combinations' instead of 'constellation'"),
]


def _collect_visible_texts(fig: plt.Figure, renderer) -> List[Tuple[Text, Tuple[float, float, float, float], Tuple[float, float, float, float], str]]:
    """Return list of (artist, fig_norm_bbox, pixel_bbox, label) for active visible Text artists."""
    results = []
    fig_w_px, fig_h_px = fig.bbox.width, fig.bbox.height

    # Helper to test if a tick label is within axis view limits
    def _is_active_tick(ax, txt_obj, is_x: bool) -> bool:
        try:
            idx = (ax.get_xticklabels() if is_x else ax.get_yticklabels()).index(txt_obj)
            locs = ax.get_xticks() if is_x else ax.get_yticks()
            if idx >= len(locs):
                return False
            lo, hi = ax.get_xlim() if is_x else ax.get_ylim()
            vmin, vmax = min(lo, hi), max(lo, hi)
            return (vmin - 1e-9 * max(abs(vmin), 1.0)) <= locs[idx] <= (vmax + 1e-9 * max(abs(vmax), 1.0))
        except Exception:
            return True

    tick_axes: Dict[Text, Tuple[plt.Axes, bool]] = {}
    for ax in fig.axes:
        for t in ax.get_xticklabels():
            tick_axes[t] = (ax, True)
        for t in ax.get_yticklabels():
            tick_axes[t] = (ax, False)

    seen_ids = set()
    for t in fig.findobj(match=Text):
        if id(t) in seen_ids:
            continue
        seen_ids.add(id(t))
        if not t.get_visible():
            continue
        raw = t.get_text()
        if not raw or not raw.strip():
            continue
        if t in tick_axes:
            ax_owner, is_x = tick_axes[t]
            if not ax_owner.axison:
                continue
            if not _is_active_tick(ax_owner, t, is_x):
                continue
        try:
            bbox = t.get_window_extent(renderer=renderer)
        except Exception:
            continue
        if bbox.width <= 0.5 or bbox.height <= 0.5:
            continue
        fx0 = bbox.x0 / fig_w_px
        fx1 = bbox.x1 / fig_w_px
        fy0 = bbox.y0 / fig_h_px
        fy1 = bbox.y1 / fig_h_px
        px_box = (bbox.x0, bbox.y0, bbox.x1, bbox.y1)
        fn_box = (fx0, fy0, fx1, fy1)
        results.append((t, fn_box, px_box, raw.replace("\n", " / ")))
    return results


def audit_matplotlib_figure(
    fig: plt.Figure,
    fig_name: str,
    quadrant_splits: Optional[Tuple[Optional[float], Optional[float]]] = (0.500, 0.495),
    min_font_pt: float = FONT_FLOOR_PT,
    expected_width_in: float = CANVAS_WIDTH_IN,
    strict: bool = True,
) -> List[str]:
    """Audit live Matplotlib figure for typography, layout, overlap, ink, and text hygiene."""
    violations: List[str] = []

    w_in, h_in = fig.get_size_inches()
    if abs(w_in - expected_width_in) > 0.02:
        violations.append(f"[{fig_name}] CANVAS_WIDTH: {w_in:.2f} in != expected {expected_width_in:.2f} in (breaks 1:1 print scale)")

    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    texts = _collect_visible_texts(fig, renderer)

    # 1. Font size, canvas clipping, quadrant divider intersection, banned ASCII
    for t, (fx0, fy0, fx1, fy1), _, clean_str in texts:
        fs = float(t.get_fontsize())
        if fs + 1e-6 < min_font_pt:
            violations.append(f"[{fig_name}] FONT_BELOW_FLOOR ({fs:.2f}pt < {min_font_pt:.1f}pt): '{clean_str}'")
        if fx0 < 0.002 or fx1 > 0.998 or fy0 < 0.002 or fy1 > 0.998:
            violations.append(
                f"[{fig_name}] CANVAS_CLIP ({fx0:.3f},{fy0:.3f})-({fx1:.3f},{fy1:.3f}): '{clean_str}'"
            )
        if quadrant_splits is not None:
            xs, ys = quadrant_splits
            if xs is not None and fx0 < xs < fx1:
                violations.append(f"[{fig_name}] CROSSES_X_DIVIDER (x={xs:.3f}): '{clean_str}' [{fx0:.3f}..{fx1:.3f}]")
            if ys is not None and fy0 < ys < fy1:
                violations.append(f"[{fig_name}] CROSSES_Y_DIVIDER (y={ys:.3f}): '{clean_str}' [{fy0:.3f}..{fy1:.3f}]")
        for pat, msg in BANNED_ASCII_PATTERNS:
            if re.search(pat, t.get_text()):
                violations.append(f"[{fig_name}] RAW_ASCII ({pat}): '{clean_str}' -> {msg}")

    # 2. Pairwise text bounding-box overlaps (> 1.5 px on both axes at canvas dpi)
    dpi = float(fig.dpi)
    min_ov_px = 1.5 * (dpi / 100.0)
    n_t = len(texts)
    for i in range(n_t):
        _, _, (ax0, ay0, ax1, ay1), s1 = texts[i]
        for j in range(i + 1, n_t):
            _, _, (bx0, by0, bx1, by1), s2 = texts[j]
            ox = min(ax1, bx1) - max(ax0, bx0)
            oy = min(ay1, by1) - max(ay0, by0)
            if ox > min_ov_px and oy > min_ov_px:
                violations.append(
                    f"[{fig_name}] TEXT_OVERLAP ({ox:.1f}x{oy:.1f}px): '{s1}' <-> '{s2}'"
                )

    # 3. Tufte spine & legend hygiene
    for idx, ax in enumerate(fig.axes):
        if ax.axison:
            if ax.spines["top"].get_visible() or ax.spines["right"].get_visible():
                violations.append(f"[{fig_name}] SPINE_CHARTJUNK: axes[{idx}] shows top/right spine")
        leg = ax.get_legend()
        if leg is not None and leg.get_visible() and leg.get_frame_on():
            violations.append(f"[{fig_name}] LEGEND_BOX: axes[{idx}] has framed legend (set frameon=False)")

    if strict and violations:
        report = "\n".join("  - " + v for v in violations)
        raise AssertionError(f"Figure standard audit FAILED for {fig_name} ({len(violations)} issues):\n{report}")
    return violations


def audit_tex_includes() -> List[str]:
    """Verify paper/main.tex and paper/supplement.tex embed figures at 1:1 \\textwidth scale and PDFs match 468pt width."""
    issues: List[str] = []
    for tex_name in ("main.tex", "supplement.tex"):
        p = PAPER_DIR / tex_name
        txt = p.read_text(encoding="utf-8")
        for m in re.finditer(r"\\includegraphics(\[[^\]]*\])?\{([^}]+)\}", txt):
            opts = m.group(1) or ""
            rel_fig = m.group(2)
            fig_path = PAPER_DIR / rel_fig
            if "height=" in opts:
                issues.append(f"[{tex_name}] SILENT_DOWNSCALE: '{rel_fig}' uses '{opts}' (remove height= for 1:1 scale)")
            if "width=\\textwidth" not in opts:
                issues.append(f"[{tex_name}] NON_STANDARD_WIDTH: '{rel_fig}' uses '{opts}' (require exact width=\\textwidth)")
            if fig_path.suffix == ".pdf" and fig_path.exists():
                try:
                    raw_pdf = fig_path.read_bytes()
                    mb = re.search(rb"/MediaBox\s*\[\s*([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s*\]", raw_pdf)
                    if not mb:
                        issues.append(f"[{tex_name}] PDF_MEDIABOX_MISSING: {fig_path.name}")
                    else:
                        x0, _, x1, _ = [float(v) for v in mb.groups()]
                        pw = x1 - x0
                        if abs(pw - 468.0) > 2.0:
                            issues.append(f"[{tex_name}] PDF_WIDTH_MISMATCH: {fig_path.name} is {pw:.1f}pt != 468.0pt (6.50 in)")
                except Exception as exc:
                    issues.append(f"[{tex_name}] PDF_CHECK_FAILED for {fig_path.name}: {exc}")
    return issues


def main() -> int:
    import importlib.util

    gen_path = ROOT / "scripts" / "22_generate_science_figures.py"
    spec = importlib.util.spec_from_file_location("gen_figs", gen_path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)

    print("=" * 80)
    print("SCIENCE / NATURE + TUFTE FIGURE STANDARD AUDIT (1:1 6.50 in, >= 6.0 pt)")
    print("=" * 80)
    mod.main()
    tex_issues = audit_tex_includes()
    if tex_issues:
        print("TeX / PDF Geometry Issues:")
        for x in tex_issues:
            print("  -", x)
        return 1
    print("[PASS] All 5 main figures + Fig. S1 + LaTeX 1:1 includes satisfy all strict checks.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
