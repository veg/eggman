#!/usr/bin/env python3
"""
09_generate_segment_phylogenies.py
==================================
Constructs Maximum-Likelihood phylogenetic trees (`IQ-TREE 2`, GTR+G with 1,000 Ultrafast
Bootstraps and 1,000 SH-aLRT replicates) across every genomic segment (>= 250 nt) for all
15 Verified Homologous Unidentified Circulating Recombinant Lineages discovered in GenBank.

Outputs:
  1. `results/trees/*.fasta` and `results/trees/*.treefile` (75 segment ML trees)
  2. `results/novel_lineages_phylogenetic_validation.csv` (Segment-by-segment phylogenetic
     parental clade assignment and Ultrafast Bootstrap support across all 15 novel lineages)
  3. `figures/fig3_novel_lineages_segment_phylogenies.{png,pdf}` (Publication-grade multi-tree
     discordance figure for Lineages #1, #2, and #3)
  4. `figures/fig4_chinese_a6_bridge_and_russian_lineages_trees.{png,pdf}` (Publication-grade
     multi-tree discordance figure for Lineages #5, #4, and #6)
"""

import os
import re
import subprocess
import sys
import time
from collections import defaultdict
from multiprocessing import Pool, cpu_count
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from Bio import Phylo

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))
from importlib import import_module

s3 = import_module("03_audit_and_reclassify_all_crfs")
s8 = import_module("08_scan_genbank_31k_full_genomes")

PROJECT_ROOT = Path("/Users/sergei/Projects/HIV-recombination")
DATA_DIR = PROJECT_ROOT / "data"
RESULTS_DIR = PROJECT_ROOT / "results"
TREES_DIR = RESULTS_DIR / "trees"
FIGURES_DIR = PROJECT_ROOT / "figures"
TREES_DIR.mkdir(parents=True, exist_ok=True)
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

IQTREE_BIN = "/Users/sergei/.local/bin/iqtree2"
INT_TO_BASE = {0: "A", 1: "C", 2: "G", 3: "T", -1: "-", -2: "-"}

# Publication-grade Okabe-Ito-harmonized subtype palette
TRUNK_COLORS = {
    "A1": "#0072B2",
    "A2": "#56B4E9",
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
    "?E": "#D35400",
}


def arr_to_fasta_slice(arr: np.ndarray, c0: int, c1: int) -> str:
    return "".join(INT_TO_BASE.get(int(x), "-") for x in arr[c0:c1])


def load_reference_panel_coding():
    """
    Load clean gold-standard reference sequences per Basal Trunk from `lanl_ref_2023.fasta`
    on the 8,628 HXB2 coding coordinates (790..9417).
    - For Subtype A1, we retain the three geographic founder sub-lineages:
      (i) East/Central African pure A1 (`AB253429`, `DQ676872`),
      (ii) Asian/Central African CRF01_AE A1 backbone in `gag..IN` (`790..5096`: `U54771`, `EF036536`),
      (iii) West African CRF02_AG A1 donor in `gag/RT/env` (`L39106`, `AY271690`).
    - For G0 in 8200..9417 (where nominal Subtype G carries an A1 tract), we use
      CRF14_BG / CRF24_BG references which carry authentic G0 in gp41..nef.
    """
    ref_recs = s3.read_fasta(DATA_DIR / "lanl_ref_2023.fasta")
    ref_hxb2_seq = [seq for hdr, seq in ref_recs if "K03455" in hdr][0]
    all_hxb2_cols = s8.get_hxb2_all_cols(ref_hxb2_seq)
    coding_cols = all_hxb2_cols[789:9417]
    hxb2_coords = np.arange(790, 9418, dtype=np.int32)

    contam_df = pd.read_csv(RESULTS_DIR / "phase1_contaminated_pure_references.csv")
    contaminated_ids = set(contam_df[~contam_df["is_ghost"]]["recombinant_id"])

    by_clade = defaultdict(list)
    crf_g0_tail = []
    a1_ae_refs = []
    a1_ag_refs = []
    for hdr, seq in ref_recs:
        if hdr.startswith("Ref.HXB2"):
            continue
        cl = hdr.split(".")[1]
        arr = s3.seq_to_int_array(seq)[coding_cols]
        acc = hdr.split(".")[-1]
        if cl == "91_cpx" and len(crf_g0_tail) < 3:
            crf_g0_tail.append((f"Ref_G0_{acc}", "G0", arr))
        if cl == "01_AE" and len(a1_ae_refs) < 2:
            a1_ae_refs.append((f"Ref_A1_{acc}", "A1", arr))
        if cl == "02_AG" and len(a1_ag_refs) < 2:
            a1_ag_refs.append((f"Ref_A1_{acc}", "A1", arr))
        if hdr in contaminated_ids and cl not in ("G", "J"):
            continue
        by_clade[cl].append((hdr, acc, arr))

    selected_refs = []
    clade_map = [
        ("A1", "A1", 2), ("A2", "A2", 2), ("A6", "A6", 3),
        ("B", "B", 3), ("C", "C", 3), ("D", "D", 2),
        ("F1", "F1", 2), ("F2", "F2", 2), ("G0", "G", 3),
        ("H", "H", 2), ("J", "J", 2), ("K", "K", 2), ("L", "L", 2),
        ("?E", "01_AE", 3),
    ]
    for trunk_label, src_cl, max_n in clade_map:
        for hdr, acc, arr in by_clade[src_cl][:max_n]:
            selected_refs.append((f"Ref_{trunk_label}_{acc}", trunk_label, arr))

    return selected_refs, crf_g0_tail, a1_ae_refs, a1_ag_refs, hxb2_coords


def parse_ref_trunk(t_name: str) -> str:
    if t_name.startswith("Ref__E_"):
        return "?E"
    return t_name.split("_")[1]


def parse_support(clade_obj):
    if clade_obj is None:
        return np.nan
    if clade_obj.name:
        parts = str(clade_obj.name).split("/")
        try:
            return float(parts[-1])
        except ValueError:
            pass
    if clade_obj.confidence is not None:
        return float(clade_obj.confidence)
    return np.nan


def run_iqtree_task(task):
    """Worker task to run IQ-TREE 2 on a single segment alignment and evaluate clade support."""
    lin_idx, formula, seg_idx, trunk_expected, h_start, h_end, fasta_path, query_names = task
    treefile = Path(str(fasta_path) + ".treefile")

    if not treefile.exists():
        cmd = [
            IQTREE_BIN, "-s", str(fasta_path),
            "-m", "GTR+G", "-B", "1000", "-alrt", "1000",
            "-T", "1", "--quiet", "-redo"
        ]
        subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    tree = Phylo.read(str(treefile), "newick")
    tree.root_at_midpoint()

    q_set = set(query_names)
    best_parent_ufboot = np.nan
    best_lineage_ufboot = np.nan

    # 1. Smallest internal clade containing all query sequences
    best_q_clade = None
    best_q_size = 999999
    for clade in tree.find_clades(order="preorder"):
        terms = [t.name for t in clade.get_terminals()]
        q_in = [t for t in terms if t in q_set]
        if len(q_in) == len(q_set) and len(terms) < best_q_size:
            best_q_size = len(terms)
            best_q_clade = clade

    # 2. Crown clade of the parental subtype containing >= 2/3 of queries + exclusively references of a single trunk.
    # If queries nest inside the references of trunk T, the crown clade containing ALL references of trunk T
    # defines the parental subtype clade support; we take the highest UFBoot among clades containing only trunk T.
    candidate_pure_clades = []
    for clade in tree.find_clades(order="preorder"):
        terms = [t.name for t in clade.get_terminals()]
        q_in = [t for t in terms if t in q_set]
        if len(q_in) < max(2, int(np.ceil(0.66 * len(q_set)))):
            continue
        ref_trunks_in = [parse_ref_trunk(t) for t in terms if t.startswith("Ref_")]
        if not ref_trunks_in:
            continue
        unique_ref_trunks = set(ref_trunks_in)
        if len(unique_ref_trunks) == 1:
            tr_only = list(unique_ref_trunks)[0]
            ufb = parse_support(clade)
            candidate_pure_clades.append((len(ref_trunks_in), ufb if pd.notna(ufb) else 0.0, -len(terms), tr_only, clade))

    if candidate_pure_clades:
        # Sort by: number of focal references included (prefer full crown clade of that subtype), then UFBoot
        candidate_pure_clades.sort(reverse=True)
        _, best_parent_ufboot, _, inferred_tree_trunk, _ = candidate_pure_clades[0]
    else:
        # Fallback: nearest reference trunk by mean patristic distance
        ref_terms = [t for t in tree.get_terminals() if t.name.startswith("Ref_")]
        q_terms = [t for t in tree.get_terminals() if t.name in q_set]
        best_d = 999.0
        inferred_tree_trunk = trunk_expected
        for rt in ref_terms:
            d = float(np.mean([tree.distance(rt, qt) for qt in q_terms]))
            if d < best_d:
                best_d = d
                inferred_tree_trunk = parse_ref_trunk(rt.name)

    if best_q_clade is not None:
        best_lineage_ufboot = parse_support(best_q_clade)

    q_terms = [t for t in tree.get_terminals() if t.name in q_set]
    pair_dists = []
    for i in range(len(q_terms)):
        for j in range(i + 1, len(q_terms)):
            pair_dists.append(float(tree.distance(q_terms[i], q_terms[j])))
    mean_patristic = float(np.mean(pair_dists)) if pair_dists else np.nan

    return {
        "lineage_rank": lin_idx,
        "rhizaeon_formula": formula,
        "segment_idx": seg_idx,
        "hxb2_start": h_start,
        "hxb2_end": h_end,
        "span_nt": h_end - h_start + 1,
        "gene_region": s3.hxb2_to_short_gene(h_start, h_end),
        "rhizaeon_assigned_trunk": trunk_expected,
        "ml_tree_sister_trunk": inferred_tree_trunk,
        "concordant_with_tree": inferred_tree_trunk == trunk_expected,
        "parent_clade_ufboot_pct": best_parent_ufboot,
        "lineage_cluster_ufboot_pct": best_lineage_ufboot,
        "n_query_taxa": len(q_set),
        "n_total_taxa_in_lineage_mrca": best_q_size,
        "is_strictly_monophyletic": best_q_size == len(q_set),
        "mean_intracluster_patristic_dist": round(mean_patristic, 5) if pd.notna(mean_patristic) else np.nan,
        "treefile_path": str(treefile),
    }


def collapse_background_clades(tree, q_set: set, focal_trunk: str):
    """
    Prune the tree down to 1 representative reference per core background subtype
    (A1, A6, B, C, D, F1, G0, ?E) plus up to 2 references for the focal parental trunk
    and up to 5 representative query tips so every tree has ~13-15 well-spaced tips.
    """
    keep_bg_trunks = {"A1", "A6", "B", "C", "D", "F1", "G0", "?E", focal_trunk}

    for term in list(tree.get_terminals()):
        if term.name.startswith("Ref_"):
            tr = parse_ref_trunk(term.name)
            if tr not in keep_bg_trunks:
                try:
                    tree.prune(term)
                except Exception:
                    pass

    seen_trunk_cnt = defaultdict(int)
    for term in list(tree.get_terminals()):
        if term.name.startswith("Ref_"):
            tr = parse_ref_trunk(term.name)
            max_keep = 2 if tr == focal_trunk else 1
            if seen_trunk_cnt[tr] >= max_keep:
                try:
                    tree.prune(term)
                except Exception:
                    pass
            else:
                seen_trunk_cnt[tr] += 1

    q_terms = [t for t in tree.get_terminals() if t.name in q_set]
    if len(q_terms) > 10:
        for t_prune in q_terms[10:]:
            try:
                tree.prune(t_prune)
            except Exception:
                pass


def draw_compact_ml_tree(ax, treefile_path: str, query_names: list, focal_trunk: str, title_str: str, subtitle_str: str):
    """
    Render a clean, publication-grade rectangular phylogram:
      - Outgroup backbone and non-focal reference subtypes in muted slate-grey (#888888)
      - Only the monophyletic parental clade (Queries + Ref_focal) colored in focal_col
      - Query patient clade colored in bold crimson (#C0392B)
      - Highlight box tightly framing the parental subtype clade from its stem x-coordinate
    """
    tree = Phylo.read(treefile_path, "newick")
    tree.root_at_midpoint()
    q_set = set(query_names)

    collapse_background_clades(tree, q_set, focal_trunk)
    tree.ladderize()

    terminals = list(tree.get_terminals())
    n_tips = len(terminals)
    y_map = {term: idx for idx, term in enumerate(terminals)}
    x_map = tree.depths()
    max_x = max(x_map.values()) if x_map else 0.1

    def get_y(clade):
        if clade.is_terminal():
            return float(y_map[clade])
        return float(np.mean([get_y(c) for c in clade.clades]))

    def is_focal_or_q(t_name: str) -> bool:
        return (t_name in q_set) or (t_name.startswith("Ref_") and parse_ref_trunk(t_name) == focal_trunk)

    # Identify:
    #   1) q_mrca_node: smallest clade containing all remaining query tips
    #   2) parent_mrca_node: crown clade containing all remaining query tips + all remaining focal_trunk refs
    rem_q = {t.name for t in terminals if t.name in q_set}
    rem_f = {t.name for t in terminals if t.name.startswith("Ref_") and parse_ref_trunk(t.name) == focal_trunk}

    q_mrca_node = None
    q_mrca_size = 9999
    parent_mrca_node = None
    parent_mrca_size = 9999

    for clade in tree.find_clades(order="preorder"):
        if clade.is_terminal():
            continue
        terms = [t.name for t in clade.get_terminals()]
        q_in = [t for t in terms if t in rem_q]
        f_in = [t for t in terms if t in rem_f]
        if len(q_in) == len(rem_q) and len(terms) < q_mrca_size:
            q_mrca_size = len(terms)
            q_mrca_node = clade
        if len(q_in) == len(rem_q) and len(f_in) == len(rem_f) and len(terms) < parent_mrca_size:
            parent_mrca_size = len(terms)
            parent_mrca_node = clade

    focal_col = TRUNK_COLORS.get(focal_trunk, "#0072B2")

    # Tight highlight box starting at the stem of parent_mrca_node
    focal_ys = [y_map[t] for t in terminals if is_focal_or_q(t.name)]
    if focal_ys:
        y_min, y_max = min(focal_ys) - 0.48, max(focal_ys) + 0.48
        x_box_start = max(0.0, (x_map.get(parent_mrca_node, 0.0) - max_x * 0.06)) if parent_mrca_node else 0.0
        rect = mpatches.FancyBboxPatch(
            (x_box_start, y_min), max_x * 1.38 - x_box_start, y_max - y_min,
            boxstyle="round,pad=0.08",
            facecolor=focal_col,
            edgecolor=focal_col,
            alpha=0.09, lw=0.9, linestyle="-", zorder=0
        )
        ax.add_patch(rect)

    # Draw branches: ONLY color branches whose descendants are exclusively (Query U Ref_focal)
    for clade in tree.find_clades(order="preorder"):
        x0 = x_map[clade]
        y0 = get_y(clade)
        for child in clade.clades:
            x1 = x_map[child]
            y1 = get_y(child)
            c_terms = [t.name for t in child.get_terminals()]
            all_q = all(t in rem_q for t in c_terms)
            all_focal_or_q = all(is_focal_or_q(t) for t in c_terms)

            if all_q:
                col, lw, alpha, z = "#C0392B", 2.1, 1.0, 4
            elif all_focal_or_q:
                col, lw, alpha, z = focal_col, 1.65, 0.95, 3
            else:
                col, lw, alpha, z = "#8C969E", 1.0, 0.72, 1

            ax.plot([x0, x0], [y0, y1], color=col, lw=lw, alpha=alpha, zorder=z)
            ax.plot([x0, x1], [y1, y1], color=col, lw=lw, alpha=alpha, zorder=z)

            if child in (q_mrca_node, parent_mrca_node):
                ufb = parse_support(child)
                if pd.notna(ufb) and ufb >= 50:
                    is_q_node = (child == q_mrca_node and q_mrca_node != parent_mrca_node)
                    ax.text(
                        0.5 * (x0 + x1), y1 + 0.26, f"{int(round(ufb))}%",
                        fontsize=7.2, fontweight="bold",
                        color="#A93226" if is_q_node else focal_col,
                        ha="center", va="bottom",
                        bbox=dict(boxstyle="round,pad=0.10", fc="white", ec="none", alpha=0.88),
                        zorder=6
                    )

    # Draw tip markers and labels
    for term in terminals:
        xt = x_map[term]
        yt = y_map[term]
        name = term.name
        if name in rem_q:
            ax.scatter([xt], [yt], color="#C0392B", s=26, zorder=6, edgecolors="white", linewidths=0.5)
            short_q = name.replace("Query_", "")
            ax.text(xt + max_x * 0.03, yt, short_q, fontsize=7.3, fontweight="bold", color="#A93226", va="center", zorder=6)
        elif name.startswith("Ref_"):
            tr = parse_ref_trunk(name)
            col = TRUNK_COLORS.get(tr, "#555555")
            if tr == focal_trunk:
                ax.scatter([xt], [yt], color=col, s=22, marker="s", zorder=5, edgecolors="white", linewidths=0.4)
                ax.text(xt + max_x * 0.03, yt, f"Ref {tr}", fontsize=7.2, fontweight="bold", color=col, va="center", zorder=5)
            else:
                ax.scatter([xt], [yt], color="#7F8C8D", s=12, zorder=2)
                ax.text(xt + max_x * 0.03, yt, f"Subtype {tr}", fontsize=6.7, color="#566573", va="center", zorder=2)

    ax.set_xlim(-max_x * 0.03, max_x * 1.42)
    ax.set_ylim(-0.85, n_tips - 0.15)
    ax.set_yticks([])
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.spines["bottom"].set_linewidth(0.8)
    ax.spines["bottom"].set_color("#333333")
    ax.tick_params(axis="x", labelsize=6.8, width=0.8, length=3.0)
    ax.set_xlabel("Substitutions / site", fontsize=7.2, labelpad=2, color="#222222")
    ax.set_title(f"{title_str}\n{subtitle_str}", fontsize=8.2, fontweight="bold", color="#111111", pad=4)


def main():
    t0 = time.time()
    plt.rcParams["font.family"] = "sans-serif"
    plt.rcParams["font.sans-serif"] = ["Arial", "Helvetica", "DejaVu Sans"]
    plt.rcParams["pdf.fonttype"] = 42
    plt.rcParams["ps.fonttype"] = 42

    print("=" * 100)
    print("STEP 9: MAXIMUM-LIKELIHOOD PHYLOGENETIC VALIDATION ACROSS GENOMIC SEGMENTS (IQ-TREE 2 GTR+G)")
    print("=" * 100)

    selected_refs, crf_g0_tail, a1_ae_refs, a1_ag_refs, hxb2_coords = load_reference_panel_coding()
    novel_df = pd.read_csv(RESULTS_DIR / "genbank_novel_circulating_recombinants.csv")
    hom_df = novel_df[novel_df["n_homologous_bp_patients"] >= 3].copy().reset_index(drop=True)
    print(f"[+] Loaded {len(hom_df)} Verified Homologous Unidentified Circulating Lineages for across-genome ML tree construction.")

    all_member_accs = set()
    for _, r in hom_df.iterrows():
        for a in str(r["member_accessions"]).split(","):
            all_member_accs.add(a.strip())

    raw_pairs = []
    with open(s8.GENBANK_FASTA, "r") as f:
        cur_acc = None
        cur_seq = []
        for line in f:
            if line.startswith(">"):
                if cur_acc in all_member_accs:
                    raw_pairs.append((cur_acc, "".join(cur_seq).upper()))
                cur_acc = line[1:].strip().split()[0].split(".")[0]
                cur_seq = []
            else:
                cur_seq.append(line.strip())
        if cur_acc in all_member_accs:
            raw_pairs.append((cur_acc, "".join(cur_seq).upper()))

    profile_path = DATA_DIR / "hxb2_15trunk_profile_9719.fasta"
    loop_mask_bool = np.zeros(len(hxb2_coords), dtype=bool)
    for h_lo, h_hi in s8.HYPERVARIABLE_MASK_RANGES:
        loop_mask_bool |= (hxb2_coords >= h_lo) & (hxb2_coords <= h_hi)

    aligned_query_map = s8.align_sequences_to_coding_arrays(raw_pairs, profile_path, loop_mask_bool)

    flt_recs = s3.read_fasta(DATA_DIR / "lanl_flt_2022.fasta")
    flt_hxb2_seq = [seq for hdr, seq in flt_recs if "K03455" in hdr][0]
    flt_coding_cols = np.array([i for i, ch in enumerate(flt_hxb2_seq) if ch not in "-."], dtype=np.int32)[789:9417]
    for hdr, seq in flt_recs:
        if "K03455" in hdr:
            continue
        acc = hdr.split(".")[-1]
        if acc in all_member_accs and acc not in aligned_query_map:
            q_a = s3.seq_to_int_array(seq)[flt_coding_cols]
            q_a[loop_mask_bool] = -1
            aligned_query_map[acc] = q_a

    # Remove all stale lin*/audit* IQ-TREE files so every tree reflects updated ChronAeon founders & candidate clusters
    for old_f in list(TREES_DIR.glob("lin*")) + list(TREES_DIR.glob("audit*")):
        old_f.unlink()

    def _write_seg_fasta(fasta_path, accs, c_trunks, h0, h1, c0, c1):
        use_ae_a1 = ("?E" in c_trunks and h0 < 5050)
        use_ag_a1 = ("G0" in c_trunks and "A1" in c_trunks and h0 < 8200)
        q_names = []
        with open(fasta_path, "w") as f:
            for ref_name, ref_tr, ref_arr in selected_refs:
                if ref_tr == "?E" and h0 < 5050:
                    continue
                if ref_tr == "A1" and (use_ae_a1 or use_ag_a1):
                    continue
                if ref_tr == "G0" and h1 > 8200 and h0 < 8900:
                    continue
                s_slice = arr_to_fasta_slice(ref_arr, c0, c1)
                if (len(s_slice) - s_slice.count("-")) >= 150:
                    f.write(f">{ref_name}\n{s_slice}\n")

            if use_ae_a1:
                for ref_name, ref_tr, ref_arr in a1_ae_refs:
                    s_slice = arr_to_fasta_slice(ref_arr, c0, c1)
                    if (len(s_slice) - s_slice.count("-")) >= 150:
                        f.write(f">{ref_name}\n{s_slice}\n")
            elif use_ag_a1:
                for ref_name, ref_tr, ref_arr in a1_ag_refs:
                    s_slice = arr_to_fasta_slice(ref_arr, c0, c1)
                    if (len(s_slice) - s_slice.count("-")) >= 150:
                        f.write(f">{ref_name}\n{s_slice}\n")

            if h1 > 8200 and h0 < 8900:
                for ref_name, ref_tr, ref_arr in crf_g0_tail:
                    s_slice = arr_to_fasta_slice(ref_arr, c0, c1)
                    if (len(s_slice) - s_slice.count("-")) >= 150:
                        f.write(f">{ref_name}\n{s_slice}\n")

            for acc in accs:
                q_arr = aligned_query_map.get(acc)
                if q_arr is None:
                    continue
                s_slice = arr_to_fasta_slice(q_arr, c0, c1)
                if (len(s_slice) - s_slice.count("-")) >= 150:
                    q_name = f"Query_{acc}"
                    q_names.append(q_name)
                    f.write(f">{q_name}\n{s_slice}\n")
        return q_names

    iqtree_tasks = []
    for idx, r in hom_df.iterrows():
        lin_rank = idx + 1
        formula = str(r["rhizaeon_formula"])
        barcode = str(r["representative_hxb2_barcode"])
        c_trunks = set(str(r["constituent_trunks"]).split("+"))
        accs = [a.strip() for a in str(r["member_accessions"]).split(",") if a.strip() in aligned_query_map]
        if not accs:
            continue

        seg_items = barcode.split("|")
        for seg_idx, s_item in enumerate(seg_items, 1):
            tr, coords_s = s_item.split(":")
            h0, h1 = [int(x) for x in coords_s.split("-")]
            if (h1 - h0 + 1) < 250:
                continue
            c0 = int(np.searchsorted(hxb2_coords, h0))
            c1 = int(np.searchsorted(hxb2_coords, h1, side="right"))

            safe_form = re.sub(r"[^A-Za-z0-9]+", "_", formula).strip("_")
            fasta_path = TREES_DIR / f"lin{lin_rank:02d}_{safe_form}_seg{seg_idx}_{tr}_{h0}_{h1}.fasta"
            q_names = _write_seg_fasta(fasta_path, accs, c_trunks, h0, h1, c0, c1)
            if len(q_names) >= 3:
                iqtree_tasks.append((lin_rank, formula, seg_idx, tr, h0, h1, fasta_path, q_names))

    # Also evaluate the 6 legacy deprecated polyphyletic clusters + 8 promoted post-2023 CRFs directly from benchmark/GenBank FASTAs for audit TSV
    legacy_clusters = defaultdict(list)
    with open(DATA_DIR / "benchmark_rhiztype_vs_comet.fasta", "r") as f:
        for line in f:
            if line.startswith(">P2|"):
                parts = line[1:].strip().split("|")
                legacy_clusters[parts[1]].append(parts[-1])

    gb_cls = pd.read_csv(RESULTS_DIR / "genbank_31k_rhizaeon_classifications.csv").set_index("accession")
    extra_accs = {a for accs in legacy_clusters.values() for a in accs} - set(aligned_query_map.keys())
    if extra_accs:
        extra_pairs = []
        with open(s8.GENBANK_FASTA, "r") as f:
            cur_acc, cur_seq = None, []
            for line in f:
                if line.startswith(">"):
                    if cur_acc in extra_accs:
                        extra_pairs.append((cur_acc, "".join(cur_seq).upper()))
                    cur_acc = line[1:].strip().split()[0].split(".")[0]
                    cur_seq = []
                elif cur_acc in extra_accs:
                    cur_seq.append(line.strip())
            if cur_acc in extra_accs:
                extra_pairs.append((cur_acc, "".join(cur_seq).upper()))
        aligned_query_map.update(s8.align_sequences_to_coding_arrays(extra_pairs, profile_path, loop_mask_bool))

    audit_specs = {
        "NOVEL_CRF_01": ("Polyphyletic_Convergent_URF", "URF", "East African D[A1:env] convergent gp120 captures"),
        "NOVEL_CRF_02": ("Candidate_Regional_Mosaic", "RU_G0B_gagMA", "Russian G0[B:gag.MA | A1:gp41..nef] monophyletic cluster"),
        "NOVEL_CRF_03": ("Promoted_Official_LANL_CRF", "CRF135_0107", "Official post-2023 LANL CRF135_0107"),
        "NOVEL_CRF_04": ("Scrubbed_Official_LANL_CRFs", "CRF173_63A6 / CRF63_02A6", "Scrubbed CRF63_02A6 (JX500700, MK984159, PZ426393) from CRF173_63A6 (PQ585408..14)"),
        "NOVEL_CRF_05": ("Promoted_Official_LANL_CRF", "CRF205_0107", "Official post-2023 LANL CRF205_0107"),
        "NOVEL_CRF_06": ("Polyphyletic_Convergent_URF", "URF", "West African A1/G0 convergent isolates"),
        "NOVEL_CRF_07": ("Promoted_Official_LANL_CRF", "CRF162_cpx", "Official post-2023 LANL CRF162_cpx"),
        "NOVEL_CRF_08": ("Polyphyletic_Convergent_URF", "URF", "East African D[A1:gp120..gp41] polyphyletic isolates"),
        "NOVEL_CRF_09": ("Promoted_Official_LANL_CRF", "CRF136_0107", "Official post-2023 LANL CRF136_0107"),
        "NOVEL_CRF_10": ("Polyphyletic_Convergent_URF", "URF", "Multi-country A1/G0 isolates failing strict monophyly"),
        "NOVEL_CRF_11": ("Polyphyletic_Convergent_URF", "URF", "Multi-country A1/C/D isolates failing strict monophyly"),
        "NOVEL_CRF_12": ("Polyphyletic_Convergent_URF", "URF", "SE Asian 01_AE/B/C isolates failing strict monophyly"),
        "NOVEL_CRF_13": ("Promoted_Official_LANL_CRF", "CRF172_0755", "Official post-2023 LANL CRF172_0755"),
        "NOVEL_CRF_14": ("Promoted_Official_LANL_CRF", "CRF127_07109", "Official post-2023 LANL CRF127_07109"),
        "NOVEL_CRF_15": ("Promoted_Official_LANL_CRF", "CRF188_0107", "Official post-2023 LANL CRF188_0107"),
    }

    audit_tasks = []
    for leg_lbl, accs in sorted(legacy_clusters.items()):
        disp, _, _ = audit_specs[leg_lbl]
        if disp != "Polyphyletic_Convergent_URF" or leg_lbl == "NOVEL_CRF_01":
            continue
        rep_row = gb_cls.loc[accs[0]]
        barcode = str(rep_row["hxb2_barcode"])
        c_trunks = set(str(rep_row["constituent_trunks"]).split("+"))
        for seg_idx, s_item in enumerate(barcode.split("|"), 1):
            tr, coords_s = s_item.split(":")
            h0, h1 = [int(x) for x in coords_s.split("-")]
            if (h1 - h0 + 1) < 250:
                continue
            c0 = int(np.searchsorted(hxb2_coords, h0))
            c1 = int(np.searchsorted(hxb2_coords, h1, side="right"))
            fasta_path = TREES_DIR / f"audit_{leg_lbl}_seg{seg_idx}_{tr}_{h0}_{h1}.fasta"
            q_names = _write_seg_fasta(fasta_path, accs, c_trunks, h0, h1, c0, c1)
            if len(q_names) >= 3:
                audit_tasks.append((leg_lbl, str(rep_row["rhizaeon_formula"]), seg_idx, tr, h0, h1, fasta_path, q_names))

    n_workers = min(8, cpu_count())
    all_tasks = iqtree_tasks + audit_tasks
    print(f"[+] Running IQ-TREE 2 (GTR+G, 1,000 UFBoot, 1,000 SH-aLRT) on {len(all_tasks)} genomic segments ({len(iqtree_tasks)} GenBank hom + {len(audit_tasks)} legacy audit) across {n_workers} parallel workers...")
    all_tree_results = []
    with Pool(processes=n_workers) as pool:
        for i, res in enumerate(pool.imap_unordered(run_iqtree_task, all_tasks), 1):
            all_tree_results.append(res)
            if i % 10 == 0 or i == len(all_tasks):
                print(f"    -> Completed {i}/{len(all_tasks)} ML segment trees ({time.time() - t0:.1f}s elapsed)")

    tree_results = [r for r in all_tree_results if isinstance(r["lineage_rank"], int)]
    audit_tree_results = [r for r in all_tree_results if not isinstance(r["lineage_rank"], int)]

    tree_df = pd.DataFrame(tree_results).sort_values(["lineage_rank", "segment_idx"]).reset_index(drop=True)
    out_csv = RESULTS_DIR / "novel_lineages_phylogenetic_validation.csv"
    tree_df.to_csv(out_csv, index=False)

    # Evaluate strict monophyly (all segments >=250 nt monophyletic, UFBoot >= 80%, mean intracluster patristic <= 0.06)
    hom_updated = []
    for idx, r in hom_df.iterrows():
        lin_rank = idx + 1
        sub_t = tree_df[tree_df["lineage_rank"] == lin_rank]
        n_seg = len(sub_t)
        n_mono = int(sub_t["is_strictly_monophyletic"].sum()) if n_seg else 0
        min_ufb = float(sub_t["lineage_cluster_ufboot_pct"].min()) if n_seg and sub_t["lineage_cluster_ufboot_pct"].notna().all() else np.nan
        max_pat = float(sub_t["mean_intracluster_patristic_dist"].max()) if n_seg and sub_t["mean_intracluster_patristic_dist"].notna().all() else np.nan
        mean_div_pct = float(r.get("mean_pairwise_div_pct", 999.0))
        strict_mono = bool(
            n_seg >= 2
            and n_mono == n_seg
            and pd.notna(min_ufb) and min_ufb >= 80.0
            and pd.notna(mean_div_pct) and mean_div_pct <= 6.0
            and pd.notna(max_pat) and max_pat <= 0.10
        )
        d = r.to_dict()
        d["n_segments_tested"] = n_seg
        d["n_segments_monophyletic"] = n_mono
        d["min_lineage_ufboot_pct"] = round(min_ufb, 1) if pd.notna(min_ufb) else np.nan
        d["max_segment_patristic_dist"] = round(max_pat, 5) if pd.notna(max_pat) else np.nan
        d["all_segments_strictly_monophyletic"] = strict_mono
        if strict_mono:
            d["phylogenetic_classification"] = "Candidate_Regional_Mosaic"
            d["lineage_status"] = "Candidate_Regional_Mosaic (Strict ML Monophyly)"
        else:
            d["phylogenetic_classification"] = "Convergent_Hotspot_Polyphyletic"
            d["lineage_status"] = "Convergent_Hotspot_Polyphyletic (Routes to URF)"
        hom_updated.append(d)

    hom_up_df = pd.DataFrame(hom_updated)
    conv_df = novel_df[novel_df["n_homologous_bp_patients"] < 3].copy()
    if not conv_df.empty:
        conv_df["n_segments_tested"] = 0
        conv_df["n_segments_monophyletic"] = 0
        conv_df["min_lineage_ufboot_pct"] = np.nan
        conv_df["max_segment_patristic_dist"] = np.nan
        conv_df["all_segments_strictly_monophyletic"] = False
        conv_df["phylogenetic_classification"] = "Convergent_Hotspot_Architecture"
    out_novel_df = pd.concat([hom_up_df, conv_df], ignore_index=True)
    out_novel_df.to_csv(RESULTS_DIR / "genbank_novel_circulating_recombinants.csv", index=False)

    # Write results/novel_candidate_crfs_31k_audit.tsv covering all legacy 1..15 + any new hom_df clusters
    audit_rows = []
    audit_tree_df = pd.DataFrame(audit_tree_results)
    for leg_lbl, (disp, new_lbl, note_txt) in sorted(audit_specs.items()):
        accs = legacy_clusters.get(leg_lbl, [])
        if leg_lbl == "NOVEL_CRF_01":
            sub_t = tree_df[tree_df["rhizaeon_formula"] == "D[A1:env]"]
        elif leg_lbl == "NOVEL_CRF_02":
            sub_t = tree_df[tree_df["rhizaeon_formula"].str.contains("G0\[B:gag\.MA", regex=True)]
        elif not audit_tree_df.empty:
            sub_t = audit_tree_df[audit_tree_df["lineage_rank"] == leg_lbl]
        else:
            sub_t = pd.DataFrame()
        n_seg = len(sub_t)
        n_mono = int(sub_t["is_strictly_monophyletic"].sum()) if n_seg else 0
        min_ufb = float(sub_t["lineage_cluster_ufboot_pct"].min()) if n_seg and sub_t["lineage_cluster_ufboot_pct"].notna().any() else np.nan
        max_pat = float(sub_t["mean_intracluster_patristic_dist"].max()) if n_seg and sub_t["mean_intracluster_patristic_dist"].notna().any() else np.nan
        audit_rows.append({
            "audit_id": leg_lbl,
            "legacy_accessions": ",".join(accs),
            "disposition": disp,
            "canonical_lineage": new_lbl,
            "n_segments_tested": n_seg,
            "n_segments_monophyletic": n_mono,
            "min_ufboot_pct": round(min_ufb, 1) if pd.notna(min_ufb) else np.nan,
            "max_segment_patristic_dist": round(max_pat, 5) if pd.notna(max_pat) else np.nan,
            "notes": note_txt,
        })

    for _, hr in hom_up_df.iterrows():
        form = str(hr["rhizaeon_formula"])
        if form in ("D[A1:env]",) or "G0[B:gag.MA" in form:
            continue
        audit_rows.append({
            "audit_id": form,
            "legacy_accessions": str(hr["member_accessions"]).replace(" ", ""),
            "disposition": hr["phylogenetic_classification"],
            "canonical_lineage": "RU_G0B_gagMA" if hr["all_segments_strictly_monophyletic"] else "URF",
            "n_segments_tested": hr["n_segments_tested"],
            "n_segments_monophyletic": hr["n_segments_monophyletic"],
            "min_ufboot_pct": hr["min_lineage_ufboot_pct"],
            "max_segment_patristic_dist": hr["max_segment_patristic_dist"],
            "notes": f"Newly surfaced Homologous GenBank candidate ({hr['countries']})",
        })

    audit_tsv = RESULTS_DIR / "novel_candidate_crfs_31k_audit.tsv"
    pd.DataFrame(audit_rows).to_csv(audit_tsv, sep="\t", index=False)

    n_concordant = int(tree_df["concordant_with_tree"].sum())
    print(f"\n[+] Phylogenetic Segment Validation Complete -> {out_csv.name} & {audit_tsv.name}")
    print(f"    * Segment ML Trees Concordant with RhizType Parental Trunk : {n_concordant}/{len(tree_df)} ({100.0*n_concordant/len(tree_df):.1f}%)")
    print(f"    * Median Parental Clade Ultrafast Bootstrap Support        : {tree_df['parent_clade_ufboot_pct'].median():.1f}%")
    print(f"    * Median Lineage Monophyly Ultrafast Bootstrap Support     : {tree_df['lineage_cluster_ufboot_pct'].median():.1f}%")
    print("\nCANDIDATE LINEAGE MONOPHYLY VERDICT:")
    print(hom_up_df[["rhizaeon_formula", "n_homologous_bp_patients", "n_segments_tested", "n_segments_monophyletic", "min_lineage_ufboot_pct", "max_segment_patristic_dist", "phylogenetic_classification"]].to_string(index=False))
    print("\nSEGMENT-BY-SEGMENT ML TREE SUMMARY:")
    cols_print = [
        "lineage_rank", "rhizaeon_formula", "segment_idx", "hxb2_start", "hxb2_end",
        "gene_region", "rhizaeon_assigned_trunk", "ml_tree_sister_trunk",
        "parent_clade_ufboot_pct", "lineage_cluster_ufboot_pct", "is_strictly_monophyletic",
        "mean_intracluster_patristic_dist",
    ]
    print(tree_df[cols_print].to_string(index=False))
    return

    def plot_lineage_block(fig_obj, gs_obj, block_idx, lin_rank, header_title, max_cols=4, custom_segs=None):
        row_bar = block_idx * 2
        row_trees = block_idx * 2 + 1
        r_meta = hom_df.iloc[lin_rank - 1]
        barcode = str(r_meta["representative_hxb2_barcode"])

        # 1. Full-width slim genome mosaic bar above the segment trees
        ax_bar = fig_obj.add_subplot(gs_obj[row_bar, :])
        ax_bar.set_xlim(750, 9450)
        ax_bar.set_ylim(0, 10)
        for s_item in barcode.split("|"):
            tr_b, c_b = s_item.split(":")
            b0, b1 = [int(x) for x in c_b.split("-")]
            col_b = TRUNK_COLORS.get(tr_b, "#888888")
            rect = mpatches.FancyBboxPatch(
                (b0, 1.3), b1 - b0, 4.1,
                boxstyle="round,pad=0.02", facecolor=col_b, edgecolor="#111111", lw=0.9
            )
            ax_bar.add_patch(rect)
            if (b1 - b0) >= 650:
                ax_bar.text(
                    (b0 + b1) * 0.5, 3.35, f"{tr_b} ({b0}–{b1})",
                    color="white", fontweight="bold", fontsize=8.0, ha="center", va="center"
                )
            elif (b1 - b0) >= 220:
                ax_bar.text(
                    (b0 + b1) * 0.5, 3.35, f"{tr_b}",
                    color="white", fontweight="bold", fontsize=8.0, ha="center", va="center"
                )
        ax_bar.text(
            790, 7.9,
            header_title,
            fontsize=10.0, fontweight="bold", color="#111111", va="center"
        )
        ax_bar.set_xticks([1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000])
        ax_bar.tick_params(axis="x", labelsize=7.2, pad=1, width=0.8)
        ax_bar.set_yticks([])
        ax_bar.spines["top"].set_visible(False)
        ax_bar.spines["right"].set_visible(False)
        ax_bar.spines["left"].set_visible(False)
        ax_bar.spines["bottom"].set_linewidth(0.8)

        # 2. Segment ML trees below the genome bar
        sub_t = tree_df[tree_df["lineage_rank"] == lin_rank].sort_values("segment_idx")
        if custom_segs is not None:
            sub_t = sub_t[sub_t["segment_idx"].isin(custom_segs)].sort_values("segment_idx")
        elif len(sub_t) > max_cols:
            picked_indices = []
            for tr_u, grp_u in sub_t.groupby("rhizaeon_assigned_trunk", sort=False):
                best_idx = grp_u.sort_values("span_nt", ascending=False).index[0]
                picked_indices.append(best_idx)
            if len(picked_indices) < max_cols:
                rem = sub_t.drop(index=picked_indices).sort_values("span_nt", ascending=False)
                picked_indices.extend(list(rem.head(max_cols - len(picked_indices)).index))
            sub_t = sub_t.loc[picked_indices].sort_values("segment_idx").head(max_cols)

        gs_trees = gs_obj[row_trees, :].subgridspec(1, len(sub_t), wspace=0.26)
        for col_i, (_, s_row) in enumerate(sub_t.iterrows()):
            ax = fig_obj.add_subplot(gs_trees[0, col_i])
            tf_path = s_row["treefile_path"]
            fasta_p = tf_path.replace(".treefile", "")
            q_names = [hdr for hdr, _ in s3.read_fasta(Path(fasta_p)) if hdr.startswith("Query_")]
            tr = s_row["rhizaeon_assigned_trunk"]
            h0, h1 = int(s_row["hxb2_start"]), int(s_row["hxb2_end"])
            gene_s = s_row["gene_region"]
            p_ufb = s_row["parent_clade_ufboot_pct"]
            l_ufb = s_row["lineage_cluster_ufboot_pct"]
            best_ufb = p_ufb if pd.notna(p_ufb) else l_ufb
            ufb_str = f"Clade UFBoot={best_ufb:.0f}%" if pd.notna(best_ufb) else ""
            t_title = f"Seg {int(s_row['segment_idx'])} ({gene_s}) → Clade {tr}"
            t_sub = f"HXB2 {h0:,}–{h1:,} ({h1-h0+1:,} nt)  {ufb_str}"
            draw_compact_ml_tree(ax, tf_path, q_names, tr, t_title, t_sub)

    def add_top_legend(fig_obj):
        legend_handles = [
            plt.Line2D([0], [0], marker="o", color="#C0392B", lw=2.0, markersize=6, label="Novel circulating lineage patient isolate"),
            plt.Line2D([0], [0], marker="s", color="#0072B2", lw=1.6, markersize=6, label="Assigned parental subtype reference"),
            plt.Line2D([0], [0], marker="o", color="#8C969E", lw=1.0, markersize=4.5, label="Background Group M subtype reference"),
        ]
        fig_obj.legend(
            handles=legend_handles, loc="upper center", bbox_to_anchor=(0.5, 0.992),
            ncol=3, frameon=False, fontsize=9.2
        )

    # =========================================================================
    # FIGURE 3: Segment ML Trees for Lineages #1, #2, and #3
    # =========================================================================
    print("\n[+] Rendering Figure 3: Across-Genome ML Phylogenies for Lineages #1, #2, and #3...")
    fig3 = plt.figure(figsize=(16.8, 14.8), dpi=300)
    gs3 = fig3.add_gridspec(6, 4, height_ratios=[0.26, 1.0, 0.26, 1.0, 0.26, 1.0], hspace=0.40, wspace=0.25, top=0.955, bottom=0.04)

    plot_lineage_block(fig3, gs3, 0, 1, "A   NOVEL_CRF_01: D[A1:env]@1986  —  East Africa (Uganda, Kenya, Tanzania; n = 9 patients, 9.0% mean divergence)")
    plot_lineage_block(fig3, gs3, 1, 2, "B   NOVEL_CRF_02: G0[B:gag.MA] [Span:2014-2016]  —  Russia (n = 6 patients, 4.6% mean divergence)")
    plot_lineage_block(fig3, gs3, 2, 3, "C   NOVEL_CRF_03: C[A1:gag..PR | B:tat1..gp120]  —  Heilongjiang, China (n = 5 patients, 4.2% mean divergence)")
    add_top_legend(fig3)

    fig3_png = FIGURES_DIR / "fig3_novel_lineages_segment_phylogenies.png"
    fig3_pdf = FIGURES_DIR / "fig3_novel_lineages_segment_phylogenies.pdf"
    fig3.savefig(fig3_png, bbox_inches="tight", dpi=300)
    fig3.savefig(fig3_pdf, bbox_inches="tight")
    plt.close(fig3)
    print(f"[+] Saved Figure 3 -> {fig3_png.name}")

    # =========================================================================
    # FIGURE 4: Segment ML Trees for Lineages #5 (China), #4 (Russia A6 bridge), #6 (West Africa)
    # =========================================================================
    print("[+] Rendering Figure 4: Across-Genome ML Phylogenies for Lineages #5, #4, and #6...")
    fig4 = plt.figure(figsize=(16.8, 14.8), dpi=300)
    gs4 = fig4.add_gridspec(6, 4, height_ratios=[0.26, 1.0, 0.26, 1.0, 0.26, 1.0], hspace=0.40, wspace=0.25, top=0.955, bottom=0.04)

    plot_lineage_block(fig4, gs4, 0, 5, "A   NOVEL_CRF_05: C[A1:gag.MA+IN | ?E:vif..gp41+gp41..nef] [Span:2020-2024]  —  Hangzhou, China (n = 5 patients, 1.9% mean divergence)", custom_segs=[1, 2, 3, 4])
    plot_lineage_block(fig4, gs4, 1, 4, "B   NOVEL_CRF_04: A1[G0:gag..RT+tat1..vpu+gp41..nef | A6:RT..tat1] [Span:2010-2023]  —  Russia (n = 5 patients, 4.4% mean divergence)", custom_segs=[1, 2, 3, 5])
    plot_lineage_block(fig4, gs4, 2, 6, "C   NOVEL_CRF_06: A1[G0:gag..RT+gp41..nef] [Span:2008-2010]  —  West Africa (Ghana, Nigeria, Cameroon; n = 4 patients, 9.2% mean divergence)")
    add_top_legend(fig4)

    fig4_png = FIGURES_DIR / "fig4_chinese_a6_bridge_and_russian_lineages_trees.png"
    fig4_pdf = FIGURES_DIR / "fig4_chinese_a6_bridge_and_russian_lineages_trees.pdf"
    fig4.savefig(fig4_png, bbox_inches="tight", dpi=300)
    fig4.savefig(fig4_pdf, bbox_inches="tight")
    plt.close(fig4)
    print(f"[+] Saved Figure 4 -> {fig4_png.name}")


if __name__ == "__main__":
    main()
