#!/usr/bin/env python3
"""
Step 3: Exhaustive RhizAeon + ChronAeon Audit & Re-classification of All LANL CRFs & Subtypes.

1. Audits all 1,054 LANL-declared coding breakpoints across the 145 sequence-represented CRFs
   using RhizAeon's local 2x2 Fisher exact test + Poisson floor (>=3 informative SNPs).
2. Deconvolves all 162 sequence-represented LANL Subtypes & CRFs (and 19 coordinate-only CRFs)
   directly against clean Basal Reference Trunk consensuses (scrubbed of the 12 verified
   contaminated tracts identified in Phase 1) + an unpolarized [GHOST] sink channel that
   resolves the unsampled `?E` (`CRF01_AE` vif..nef) donor lineage.
3. Combines RhizAeon's PA-FF v2.1 segments with ChronAeon's t_MRCA dates to generate the
   Descriptive Spatio-Temporal Nomenclature Compendium and Barcode Atlas.
"""

import json
import math
from pathlib import Path
from collections import defaultdict, Counter
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

PROJECT_ROOT = Path("/Users/sergei/Projects/HIV-recombination")
DATA_DIR = PROJECT_ROOT / "data"
RESULTS_DIR = PROJECT_ROOT / "results"
PROFILES_9719_TSV = PROJECT_ROOT / "rhiztype" / "crates" / "eggman-core" / "data" / "profiles_9719.tsv"

BASE_TO_INT = {ord('A'): 0, ord('C'): 1, ord('G'): 2, ord('T'): 3}
HYPERVARIABLE_MASK_RANGES = [
    (2130, 2185),  # gag p6 PTAP
    (6615, 6812),  # gp120 V1/V2
    (7110, 7217),  # gp120 V3
    (7377, 7478),  # gp120 V4
    (7596, 7637),  # gp120 V5
]


def load_profiles_9719_coding():
    """
    Load the scrubbed C1d-A ChronAeon basal trunk profiles + CM240 from
    eggman-core's profiles_9719.tsv on the exact 8,628-nt coding frame (790..9417).
    """
    hxb2_coords = np.arange(790, 9418, dtype=np.int32)
    raw_map = {}
    with open(PROFILES_9719_TSV, "r") as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            name, seq = line.split("\t", 1)
            raw_map[name] = seq[789:9417].upper()

    trunk_names = ["A1", "A2", "A6", "B", "C", "D", "F1", "F2", "G0", "H", "J", "K", "L", "?E"]
    trunk_rows = []
    for t_name in trunk_names:
        arr = seq_to_int_array(raw_map[t_name])
        if t_name == "?E":
            arr[hxb2_coords < 5097] = -2
        trunk_rows.append(arr)
    trunk_mat = np.vstack(trunk_rows)
    cm240_arr = seq_to_int_array(raw_map["CM240"])
    return trunk_mat, trunk_names, hxb2_coords, cm240_arr


def read_fasta(path: Path):
    records = []
    header = None
    chunks = []
    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(chunks).upper()))
                header = line[1:].strip()
                chunks = []
            else:
                chunks.append(line)
        if header is not None:
            records.append((header, "".join(chunks).upper()))
    return records


def seq_to_int_array(seq: str) -> np.ndarray:
    arr = np.full(len(seq), -1, dtype=np.int8)
    b = seq.encode("ascii")
    for char_code, val in BASE_TO_INT.items():
        arr[np.frombuffer(b, dtype=np.uint8) == char_code] = val
    return arr


def build_plurality_consensus(int_arrays, mask_intervals=None):
    """Build plurality consensus across (seq_arr, mask_list) entries."""
    L = len(int_arrays[0])
    mat = np.vstack(int_arrays).copy()
    if mask_intervals:
        for row_idx, intervals in enumerate(mask_intervals):
            for c0, c1 in intervals:
                mat[row_idx, c0:c1] = -1
    cons = np.full(L, -1, dtype=np.int8)
    counts = np.zeros((4, L), dtype=np.int32)
    for val in range(4):
        counts[val] = np.sum(mat == val, axis=0)
    max_c = np.max(counts, axis=0)
    valid = max_c > 0
    cons[valid] = np.argmax(counts[:, valid], axis=0)
    return cons


def hxb2_to_short_gene(h0: int, h1: int) -> str:
    mid = 0.5 * (h0 + h1)
    span = h1 - h0

    def pt_gene(pos):
        if pos < 790:
            return "5'LTR"
        if pos <= 2292:
            if pos <= 1185:
                return "gag.MA"
            if pos <= 1878:
                return "gag.CA"
            return "gag.NC"
        if pos <= 5096:
            if pos <= 2549:
                return "PR"
            if pos <= 4229:
                return "RT"
            return "IN"
        if pos <= 5559:
            return "vif"
        if pos <= 5850:
            return "vpr"
        if pos <= 6062:
            return "tat1"
        if pos <= 6225:
            return "vpu"
        if pos <= 8795:
            if pos <= 7757:
                return "gp120"
            return "gp41"
        if pos <= 9417:
            return "nef"
        return "3'LTR"

    g0 = pt_gene(h0 + 15)
    g1 = pt_gene(h1 - 15)
    if g0 == g1:
        if g0.startswith("gag.") and span > 650:
            return "gag"
        if g0 in ("PR", "RT", "IN") and span > 1800:
            return "pol"
        if g0 in ("gp120", "gp41") and span > 1800:
            return "env"
        return g0
    if g0.startswith("gag.") and g1.startswith("gag."):
        if span > 900:
            return "gag"
        return f"{g0}..{g1}"
    if g0 in ("PR", "RT", "IN") and g1 in ("PR", "RT", "IN") and g0 == "PR" and g1 == "IN":
        return "pol"
    if g0 == "gp120" and g1 == "gp41":
        if span > 1800:
            return "env"
        return "gp120..gp41"
    if g0.startswith("gag.") and not g1.startswith("gag."):
        g0 = "gag"
    if g1.startswith("gag.") and not g0.startswith("gag."):
        g1 = "gag"
    return f"{g0}..{g1}"


def test_boundary_significance(q_arr, trunk_mat, ghost_idx, col_start, col_bp, col_end, tL, tR,
                               min_snps=3, p_thresh=1e-3):
    """
    Evaluate 2x2 Fisher exact test and informative-site support across a boundary col_bp
    between left state tL and right state tR. Supports both standard trunk transitions
    and transitions involving the unpolarized [GHOST] sink channel.
    """
    valid_q_L = (q_arr[col_start:col_bp] >= 0)
    valid_q_R = (q_arr[col_bp:col_end] >= 0)

    if tL != ghost_idx and tR != ghost_idx:
        refL_L = trunk_mat[tL, col_start:col_bp]
        refR_L = trunk_mat[tR, col_start:col_bp]
        # If refR_L is domain-masked across >50% of the left segment (e.g. ?E outside vif..nef), compare refL_L against alternative trunks
        if np.mean(refR_L >= 0) < 0.5:
            other_match = np.any(
                (trunk_mat[:, col_start:col_bp] != refL_L[None, :]) & (trunk_mat[:, col_start:col_bp] >= 0),
                axis=0
            )
            info_L = valid_q_L & (refL_L >= 0) & other_match
            n11 = int(np.sum(info_L & (q_arr[col_start:col_bp] == refL_L)))
            n12 = int(np.sum(info_L & (q_arr[col_start:col_bp] != refL_L)))
        else:
            info_L = valid_q_L & (refL_L >= 0) & (refR_L >= 0) & (refL_L != refR_L)
            n11 = int(np.sum(info_L & (q_arr[col_start:col_bp] == refL_L)))
            n12 = int(np.sum(info_L & (q_arr[col_start:col_bp] == refR_L)))

        refL_R = trunk_mat[tL, col_bp:col_end]
        refR_R = trunk_mat[tR, col_bp:col_end]
        if np.mean(refL_R >= 0) < 0.5:
            other_match = np.any(
                (trunk_mat[:, col_bp:col_end] != refR_R[None, :]) & (trunk_mat[:, col_bp:col_end] >= 0),
                axis=0
            )
            info_R = valid_q_R & (refR_R >= 0) & other_match
            n22 = int(np.sum(info_R & (q_arr[col_bp:col_end] == refR_R)))
            n21 = int(np.sum(info_R & (q_arr[col_bp:col_end] != refR_R)))
        else:
            info_R = valid_q_R & (refL_R >= 0) & (refR_R >= 0) & (refL_R != refR_R)
            n21 = int(np.sum(info_R & (q_arr[col_bp:col_end] == refL_R)))
            n22 = int(np.sum(info_R & (q_arr[col_bp:col_end] == refR_R)))
    else:
        known_t = tR if tL == ghost_idx else tL
        refK_L = trunk_mat[known_t, col_start:col_bp]
        any_trunk_match_L = np.any(trunk_mat[:, col_start:col_bp] == q_arr[col_start:col_bp], axis=0)
        valid_L = valid_q_L & (refK_L >= 0)
        match_K_L = int(np.sum(valid_L & (q_arr[col_start:col_bp] == refK_L)))
        orphan_L = int(np.sum(valid_L & (~any_trunk_match_L)))

        refK_R = trunk_mat[known_t, col_bp:col_end]
        any_trunk_match_R = np.any(trunk_mat[:, col_bp:col_end] == q_arr[col_bp:col_end], axis=0)
        valid_R = valid_q_R & (refK_R >= 0)
        match_K_R = int(np.sum(valid_R & (q_arr[col_bp:col_end] == refK_R)))
        orphan_R = int(np.sum(valid_R & (~any_trunk_match_R)))

        if tL == ghost_idx:
            n11, n12 = orphan_L, match_K_L
            n21, n22 = orphan_R, match_K_R
        else:
            n11, n12 = match_K_L, orphan_L
            n21, n22 = match_K_R, orphan_R

    if n11 >= min_snps and n22 >= min_snps and (n11 * n22 > n12 * n21):
        _, pval = fisher_exact([[n11, n12], [n21, n22]], alternative="greater")
    else:
        pval = 1.0

    # Informative site brackets
    if tL != ghost_idx and tR != ghost_idx:
        refL = trunk_mat[tL]
        refR = trunk_mat[tR]
        info_cols = np.where(
            (q_arr >= 0) & (refL >= 0) & (refR >= 0) & (refL != refR) &
            ((q_arr == refL) | (q_arr == refR))
        )[0]
    else:
        known_t = tR if tL == ghost_idx else tL
        refK = trunk_mat[known_t]
        any_match = np.any(trunk_mat == q_arr[None, :], axis=0)
        info_cols = np.where((q_arr >= 0) & (refK >= 0) & ((q_arr == refK) | (~any_match)))[0]

    left_info = info_cols[(info_cols >= col_start) & (info_cols < col_bp)]
    right_info = info_cols[(info_cols >= col_bp) & (info_cols < col_end)]
    ci_low_col = int(left_info[-1]) if len(left_info) > 0 else int(col_bp)
    ci_high_col = int(right_info[0]) if len(right_info) > 0 else int(col_bp)

    return pval, n11, n12, n21, n22, ci_low_col, ci_high_col


def rhizaeon_paff_viterbi(
    q_arr: np.ndarray,
    trunk_mat: np.ndarray,
    trunk_names: list,
    col_to_hxb2: np.ndarray,
    cm240_arr: np.ndarray,
    gamma: float = 10.5,
    mismatch_cost: float = 1.0,
    ghost_cost: float = 0.58,
    min_informative_snps: int = 5,
    min_hxb2_span_nt: int = 120,
    p_fisher_thresh: float = 1e-4,
    e_penalty: float = 0.035,
):
    """
    RhizAeon PA-FF v2.1:
      1. K-state + [GHOST] Viterbi DP over coding alignment columns.
      2. Iterative statistical boundary pruning (Fisher exact test p < 1e-3 & >= 4 informative SNPs).
      3. Boundary polishing & exact informative-site CI calculation.
      4. Post-hoc identification of [GHOST] segments matching CRF01_AE's unsampled `?E` lineage.
    """
    K, L = trunk_mat.shape
    K_total = K + 1
    ghost_idx = K

    valid_q = q_arr >= 0
    valid_ref = trunk_mat >= 0
    forbidden_ref = trunk_mat == -2
    indel_gap_ref = trunk_mat == -1
    mismatch = (trunk_mat != q_arr[None, :]) & valid_q[None, :] & valid_ref
    emit = np.zeros((K_total, L), dtype=np.float64)
    emit[:K, :] = mismatch.astype(np.float64) * mismatch_cost
    # Ordinary alignment indel gaps (-1) incur neutral substitution cost (0.35);
    # Forbidden domain masks (-2, e.g. ?E in gag..pol < 5097) incur high penalty (2.0)
    emit[:K, :] += (indel_gap_ref & valid_q[None, :]).astype(np.float64) * 0.35
    emit[:K, :] += (forbidden_ref & valid_q[None, :]).astype(np.float64) * 2.00
    if e_penalty > 0.0 and "?E" in trunk_names:
        e_idx = trunk_names.index("?E")
        emit[e_idx, :] += (valid_q & valid_ref[e_idx, :]).astype(np.float64) * e_penalty
    emit[ghost_idx, :] = valid_q.astype(np.float64) * ghost_cost

    dp = np.zeros((K_total, L), dtype=np.float64)
    back = np.zeros((K_total, L), dtype=np.int16)
    dp[:, 0] = emit[:, 0]

    for u in range(1, L):
        prev = dp[:, u - 1]
        best_prev_idx = int(np.argmin(prev))
        best_prev_val = prev[best_prev_idx]
        stay_cost = prev
        switch_cost = best_prev_val + gamma
        switched = switch_cost < stay_cost
        dp[:, u] = np.where(switched, switch_cost, stay_cost) + emit[:, u]
        back[:, u] = np.where(switched, best_prev_idx, np.arange(K_total))

    path = np.zeros(L, dtype=np.int16)
    path[-1] = int(np.argmin(dp[:, -1]))
    for u in range(L - 1, 0, -1):
        path[u - 1] = back[path[u], u]

    segs = []
    cur_t = int(path[0])
    start_u = 0
    for u in range(1, L):
        if int(path[u]) != cur_t:
            segs.append([start_u, u, cur_t])
            cur_t = int(path[u])
            start_u = u
    segs.append([start_u, L, cur_t])

    # Iterative boundary pruning
    changed = True
    while changed and len(segs) > 1:
        changed = False
        best_merge_idx = -1
        worst_pval = -1.0

        for i in range(len(segs) - 1):
            s1, e1, t1 = segs[i]
            s2, e2, t2 = segs[i + 1]
            if t1 == t2:
                best_merge_idx = i
                worst_pval = 2.0
                break
            span1 = int(col_to_hxb2[e1 - 1] - col_to_hxb2[s1])
            span2 = int(col_to_hxb2[e2 - 1] - col_to_hxb2[s2])
            pval, n11, n12, n21, n22, _, _ = test_boundary_significance(
                q_arr, trunk_mat, ghost_idx, s1, e1, e2, t1, t2,
                min_snps=min_informative_snps, p_thresh=p_fisher_thresh
            )
            if (span1 < min_hxb2_span_nt or span2 < min_hxb2_span_nt or
                    pval > p_fisher_thresh or n11 < min_informative_snps or n22 < min_informative_snps):
                score = pval + (1.5 if (span1 < min_hxb2_span_nt or span2 < min_hxb2_span_nt) else 0.0)
                if score > worst_pval:
                    worst_pval = score
                    best_merge_idx = i

        if best_merge_idx >= 0:
            i = best_merge_idx
            s1, e1, t1 = segs[i]
            s2, e2, t2 = segs[i + 1]
            span1 = int(col_to_hxb2[e1 - 1] - col_to_hxb2[s1])
            span2 = int(col_to_hxb2[e2 - 1] - col_to_hxb2[s2])
            # Never allow a partially-masked reference (like ?E outside 5820..9417) to absorb a region where it is masked
            valid1 = np.mean(trunk_mat[t1, s1:e2] >= 0) if t1 < K else 1.0
            valid2 = np.mean(trunk_mat[t2, s1:e2] >= 0) if t2 < K else 1.0
            if valid1 < 0.7 and valid2 >= 0.7:
                keep_t = t2
            elif valid2 < 0.7 and valid1 >= 0.7:
                keep_t = t1
            else:
                keep_t = t1 if span1 >= span2 else t2
            segs[i] = [s1, e2, keep_t]
            segs.pop(i + 1)
            changed = True

    # Merge identical adjacent states
    merged = []
    for s in segs:
        if merged and merged[-1][2] == s[2]:
            merged[-1][1] = s[1]
        else:
            merged.append(s)
    segs = merged

    # Build final segment and breakpoint metadata
    breakpoints = []
    for i in range(len(segs) - 1):
        s1, e1, t1 = segs[i]
        s2, e2, t2 = segs[i + 1]
        pval, n11, n12, n21, n22, ci_low_col, ci_high_col = test_boundary_significance(
            q_arr, trunk_mat, ghost_idx, s1, e1, e2, t1, t2,
            min_snps=min_informative_snps, p_thresh=p_fisher_thresh
        )
        breakpoints.append({
            "bp_col": int(e1),
            "bp_hxb2": int(col_to_hxb2[e1]),
            "ci_hxb2_low": int(col_to_hxb2[ci_low_col]),
            "ci_hxb2_high": int(col_to_hxb2[ci_high_col]),
            "ci_width_nt": int(col_to_hxb2[ci_high_col] - col_to_hxb2[ci_low_col]),
            "left_trunk": trunk_names[t1] if t1 < K else "?Ghost",
            "right_trunk": trunk_names[t2] if t2 < K else "?Ghost",
            "p_fisher": float(pval),
            "n11": n11,
            "n12": n12,
            "n21": n21,
            "n22": n22,
        })

    out_segs = []
    for i, (s, e, t) in enumerate(segs):
        h0 = int(col_to_hxb2[s])
        h1 = int(col_to_hxb2[e - 1])
        ci_start = breakpoints[i - 1]["ci_hxb2_low"] if i > 0 else h0
        ci_end = breakpoints[i]["ci_hxb2_high"] if i < len(breakpoints) else h1
        if t < K:
            t_label = trunk_names[t]
        else:
            valid_sub = (q_arr[s:e] >= 0) & (cm240_arr[s:e] >= 0)
            if np.sum(valid_sub) > 50 and h0 >= 5000:
                match_e = float(np.sum(q_arr[s:e][valid_sub] == cm240_arr[s:e][valid_sub])) / float(np.sum(valid_sub))
                t_label = "?E" if match_e >= 0.85 else "?Ghost"
            else:
                t_label = "?Ghost"

        out_segs.append({
            "seg_idx": i + 1,
            "trunk": t_label,
            "col_start": int(s),
            "col_end": int(e),
            "hxb2_start": h0,
            "hxb2_end": h1,
            "span_nt": h1 - h0 + 1,
            "ci_hxb2_start": ci_start,
            "ci_hxb2_end": ci_end,
            "short_gene": hxb2_to_short_gene(h0, h1),
        })

    for i, bp in enumerate(breakpoints):
        bp["left_trunk"] = out_segs[i]["trunk"]
        bp["right_trunk"] = out_segs[i + 1]["trunk"]

    return out_segs, breakpoints


def format_descriptive_name(segs, t_mrca=None, origin_basis="Clock_tMRCA", span_years=None):
    trunk_spans = Counter()
    for s in segs:
        trunk_spans[s["trunk"]] += s["span_nt"]

    # Prefer full-length known basal trunks over incomplete/extinct ghost cassettes (?E, ?Ghost) as backbone
    known_spans = [(tr, sp) for tr, sp in trunk_spans.most_common() if not tr.startswith("?")]
    backbone = known_spans[0][0] if known_spans else trunk_spans.most_common(1)[0][0]
    if origin_basis == "Clock_tMRCA" and t_mrca is not None and not np.isnan(t_mrca):
        yr_suffix = f"@{int(round(t_mrca))}"
    elif span_years is not None and span_years[0] is not None and not np.isnan(span_years[0]):
        y0 = int(round(span_years[0]))
        y1 = int(round(span_years[1])) if (span_years[1] is not None and not np.isnan(span_years[1])) else y0
        yr_suffix = f"[Span:{y0}-{y1}]" if y0 != y1 else f"[Span:{y0}]"
    else:
        yr_suffix = ""

    if len(trunk_spans) == 1:
        return f"{backbone}(Clonal){yr_suffix}", backbone, f"{backbone}(Clonal)"

    donor_genes = defaultdict(list)
    for s in segs:
        if s["trunk"] != backbone:
            g = s["short_gene"]
            if not donor_genes[s["trunk"]] or donor_genes[s["trunk"]][-1] != g:
                donor_genes[s["trunk"]].append(g)

    donor_parts = []
    for donor, genes in donor_genes.items():
        donor_parts.append(f"{donor}:{'+'.join(genes)}")

    structural_formula = f"{backbone}[{' | '.join(donor_parts)}]"
    full_name = f"{structural_formula}{yr_suffix}"
    return full_name, backbone, structural_formula


def format_hxb2_barcode(segs):
    parts = []
    for s in segs:
        parts.append(f"{s['trunk']}:{s['hxb2_start']}-{s['hxb2_end']}")
    return "|".join(parts)


def main():
    print("=" * 90)
    print("STEP 3: EXHAUSTIVE RHIZAEON + CHRONAEON AUDIT & RE-CLASSIFICATION OF ALL LANL CRFs")
    print("=" * 90)

    records = read_fasta(DATA_DIR / "lanl_ref_2023.fasta")
    with open(DATA_DIR / "lanl_ref_2023_hxb2_coord_map.json") as f:
        cmap = json.load(f)

    # Map directly onto the exact 8,628 HXB2 coding positions (790..9417) matching profiles_9719.tsv
    h2m = {int(k): int(v) - 1 for k, v in cmap["hxb2_to_msa"].items()}
    coding_cols = np.array([h2m[h] for h in range(790, 9418)], dtype=np.int32)
    trunk_mat, trunk_names, col_to_hxb2, cm240_arr = load_profiles_9719_coding()
    trunk_name_to_idx = {name: i for i, name in enumerate(trunk_names)}

    loop_mask_bool = np.zeros(len(col_to_hxb2), dtype=bool)
    for h_lo, h_hi in HYPERVARIABLE_MASK_RANGES:
        loop_mask_bool |= (col_to_hxb2 >= h_lo) & (col_to_hxb2 <= h_hi)

    clade_seqs = defaultdict(list)
    for hdr, seq in records:
        if hdr.startswith("Ref.HXB2"):
            continue
        arr = seq_to_int_array(seq)[coding_cols]
        parts = hdr.split(".")
        cl = parts[1] if len(parts) > 1 else "UNK"
        clade_seqs[cl].append((hdr, arr))

    # Load ChronAeon dates
    chron_df = pd.read_csv(RESULTS_DIR / "phase1_chronaeon_clade_dating.csv")
    clade_tmrca = {}
    clade_ci = {}
    for _, r in chron_df.iterrows():
        clade_tmrca[r["clade"]] = float(r["t_mrca"])
        clade_ci[r["clade"]] = (float(r["ci_lower"]), float(r["ci_upper"]))

    flt_meta = pd.read_csv(DATA_DIR / "lanl_flt_2022_metadata.csv")
    ref_meta = pd.read_csv(DATA_DIR / "lanl_ref_2023_metadata.csv")
    earliest_year = {}
    latest_year = {}
    for df in (flt_meta, ref_meta):
        for cl, grp in df.groupby("lanl_tag"):
            yrs = grp["year"].dropna()
            if len(yrs) > 0:
                min_y = float(yrs.min())
                max_y = float(yrs.max())
                if cl not in earliest_year or min_y < earliest_year[cl]:
                    earliest_year[cl] = min_y
                if cl not in latest_year or max_y > latest_year[cl]:
                    latest_year[cl] = max_y

    with open(DATA_DIR / "lanl_official_crf_catalog.json") as f:
        crf_catalog = json.load(f)

    def map_lanl_sub(s):
        s = s.strip()
        if s in trunk_name_to_idx:
            return s, "elementary"
        if s == "A":
            return "A1", "elementary"
        if s == "F":
            return "F1", "elementary"
        if s == "G":
            return "G0", "elementary"
        if s in ("01_AE", "CRF01_AE", "E"):
            return "?E", "composite_01AE"
        if "/" in s:
            return s, "ambiguous_slash"
        if s == "U" or "?" in s:
            return "U", "unclassified_U"
        if "CRF" in s or "_" in s:
            return s, "composite_2nd_gen"
        return s, "other"

    def slash_tokens_match_flanks(slash_lbl, mL, mR):
        parts = [map_lanl_sub(p.strip())[0] for p in slash_lbl.split("/") if p.strip()]
        part_set = set(parts)
        return part_set.issubset({mL, mR}) and len(part_set) >= 1

    # 1. Audit all LANL-declared coding breakpoints across sequence-represented CRFs (bridging X -> X/Y -> Y plateaus)
    lanl_bp_audit_rows = []
    no_map_seq_crfs = []
    for crf_key, entry in sorted(crf_catalog.items(), key=lambda x: int(x[0].split("_")[0].replace("CRF", ""))):
        crf_code = entry["crf"]
        short_code = entry["short_tag"]
        if short_code not in clade_seqs:
            continue

        segs = entry.get("segments", [])
        coding_segs = [s for s in segs if s["hxb2_end"] >= 850 and s["hxb2_start"] <= 9400]
        if not coding_segs:
            no_map_seq_crfs.append(crf_code)
            continue

        q_cons = build_plurality_consensus([arr for _, arr in clade_seqs[short_code]])
        q_cons[loop_mask_bool] = -1

        bp_counter = 0
        i = 0
        while i < len(coding_segs) - 1:
            sL = coding_segs[i]
            subL = sL["parent_label"].strip()
            mL, catL = map_lanl_sub(subL)

            # Check for contiguous X/Y slash plateau between two non-slash elementary trunks X and Y
            j = i + 1
            while j < len(coding_segs) and "/" in coding_segs[j]["parent_label"]:
                j += 1

            bridged = False
            if j > i + 1 and j < len(coding_segs):
                sR_cand = coding_segs[j]
                subR_cand = sR_cand["parent_label"].strip()
                mR_cand, catR_cand = map_lanl_sub(subR_cand)
                if (
                    mL in trunk_name_to_idx
                    and mR_cand in trunk_name_to_idx
                    and mL != mR_cand
                    and all(slash_tokens_match_flanks(coding_segs[k]["parent_label"], mL, mR_cand) for k in range(i + 1, j))
                ):
                    bridged = True
                    sR = sR_cand
                    subR = subR_cand
                    mR, catR = mR_cand, catR_cand
                    slash_start = int(coding_segs[i + 1]["hxb2_start"])
                    slash_end = int(coding_segs[j - 1]["hxb2_end"])
                    bp_hxb2 = int(round(0.5 * (slash_start + slash_end)))
                    bridged_span = f"{slash_start}-{slash_end}"
                    next_i = j
            if not bridged:
                sR = coding_segs[i + 1]
                subR = sR["parent_label"].strip()
                mR, catR = map_lanl_sub(subR)
                bp_hxb2 = int(sL["hxb2_end"])
                bridged_span = ""
                next_i = i + 1

            bp_counter += 1
            col_start = int(np.searchsorted(col_to_hxb2, max(790, int(sL["hxb2_start"]))))
            col_bp = int(np.searchsorted(col_to_hxb2, bp_hxb2))
            col_end = int(np.searchsorted(col_to_hxb2, min(9417, int(sR["hxb2_end"]))))

            if mL in trunk_name_to_idx and mR in trunk_name_to_idx and mL != mR and col_start < col_bp < col_end:
                tL = trunk_name_to_idx[mL]
                tR = trunk_name_to_idx[mR]
                pval, n11, n12, n21, n22, _, _ = test_boundary_significance(
                    q_cons, trunk_mat, len(trunk_names), col_start, col_bp, col_end, tL, tR
                )
                if pval < 1e-3 and n11 >= 3 and n22 >= 3:
                    status = "Verified Elementary Breakpoint"
                elif n11 < 3 or n22 < 3:
                    status = "Failed Poisson Floor (<3 Informative SNPs)"
                else:
                    status = "Failed Fisher Exact Test (p > 1e-3)"
            else:
                n11 = n12 = n21 = n22 = 0
                pval = 1.0
                if "ambiguous_slash" in (catL, catR):
                    status = "Artifact: Ambiguous 'X/Y' Breakpoint Plateau"
                elif "unclassified_U" in (catL, catR):
                    status = "Artifact: Unclassified 'U' Segment"
                elif "composite" in catL or "composite" in catR:
                    status = "Artifact: Unresolved 2nd-Gen CRF Label"
                else:
                    status = "Artifact: Degenerate/Boundary Segment"

            lanl_bp_audit_rows.append({
                "crf_code": crf_code,
                "short_code": short_code,
                "bp_index": bp_counter,
                "lanl_bp_hxb2": bp_hxb2,
                "lanl_left_label": subL,
                "lanl_right_label": subR,
                "was_bridged_xy_interval": bridged,
                "bridged_slash_span": bridged_span,
                "left_span_nt": int(sL["hxb2_end"]) - int(sL["hxb2_start"]) + 1,
                "right_span_nt": int(sR["hxb2_end"]) - int(sR["hxb2_start"]) + 1,
                "n11_support_left": n11,
                "n12_conflict_left": n12,
                "n21_conflict_right": n21,
                "n22_support_right": n22,
                "fisher_p_value": pval,
                "audit_verdict": status,
            })
            i = next_i

    bp_audit_df = pd.DataFrame(lanl_bp_audit_rows)
    bp_audit_df.to_csv(RESULTS_DIR / "lanl_declared_breakpoints_audit.csv", index=False)
    n_mapped_crfs = bp_audit_df["crf_code"].nunique()
    n_bridged = int(bp_audit_df["was_bridged_xy_interval"].sum())
    print(f"\n--- AUDIT OF {len(bp_audit_df)} LANL-DECLARED CODING BREAKPOINTS ACROSS {n_mapped_crfs} MAPPED CRFs (+ {len(no_map_seq_crfs)} UNMAPPED: {', '.join(no_map_seq_crfs)} = {n_mapped_crfs + len(no_map_seq_crfs)} CRFs) ---")
    print(f"    (Bridged {n_bridged} contiguous X -> X/Y -> Y uncertain plateaus; verified {int((bp_audit_df.loc[bp_audit_df['was_bridged_xy_interval'], 'audit_verdict'] == 'Verified Elementary Breakpoint').sum())}/{n_bridged})")
    for verdict, cnt in bp_audit_df["audit_verdict"].value_counts().items():
        print(f"  * {verdict:<48s}: {cnt:4d} ({100.0*cnt/len(bp_audit_df):5.1f}%)")

    # 2. Deconvolve all 162 sequence-represented LANL Subtypes & CRFs (17 pure + 145 CRFs)
    compendium_rows = []
    all_deconvolved_segments = []
    all_deconvolved_bps = []

    pure_clades_set = {"A1", "A2", "A3", "A4", "A6", "A7", "A8", "B", "C", "D", "F1", "F2", "G", "H", "J", "K", "L"}
    out_of_panel_crfs = {"CRF26_A5U", "CRF92_C2U"}
    sorted_clades = sorted(
        clade_seqs.keys(),
        key=lambda c: (0 if c in pure_clades_set else 1,
                       int(c.split("_")[0]) if c[0].isdigit() and "_" in c else 999,
                       c)
    )

    for cl in sorted_clades:
        q_cons = build_plurality_consensus([arr for _, arr in clade_seqs[cl]])
        q_cons[loop_mask_bool] = -1
        cur_trunk_mat = trunk_mat.copy()
        if cl == "G":
            # For nominal Subtype G, mask G0 across the Phase-1 verified A1-shared gp41..nef tract (8262..8880)
            g0_idx = trunk_name_to_idx["G0"]
            cur_trunk_mat[g0_idx, (col_to_hxb2 >= 8262) & (col_to_hxb2 <= 8880)] = -2
        gamma_used = 10.5
        segs, bps = rhizaeon_paff_viterbi(
            q_cons, cur_trunk_mat, trunk_names, col_to_hxb2, cm240_arr,
            gamma=10.5, mismatch_cost=1.0, ghost_cost=0.58,
            min_informative_snps=5, min_hxb2_span_nt=120, p_fisher_thresh=1e-4
        )
        if cl not in pure_clades_set and len(segs) == 1:
            segs_g7, bps_g7 = rhizaeon_paff_viterbi(
                q_cons, cur_trunk_mat, trunk_names, col_to_hxb2, cm240_arr,
                gamma=7.0, mismatch_cost=1.0, ghost_cost=0.58,
                min_informative_snps=5, min_hxb2_span_nt=120, p_fisher_thresh=1e-4
            )
            if len(segs_g7) > 1:
                segs, bps = segs_g7, bps_g7
                gamma_used = 7.0

        if cl in clade_tmrca:
            t_mrca = clade_tmrca[cl]
            origin_basis = "Clock_tMRCA"
            t_src = "ChronAeon_RCS_or_OLS_BlockBootstrap"
            ci_low, ci_high = clade_ci[cl]
            span_tup = None
        else:
            t_mrca = earliest_year.get(cl, np.nan)
            origin_basis = "Sampling_Span_Only" if not np.isnan(t_mrca) else "Undated"
            t_src = "Sampling_Span_Only" if not np.isnan(t_mrca) else "Undated"
            ci_low, ci_high = (np.nan, np.nan)
            span_tup = (earliest_year.get(cl, np.nan), latest_year.get(cl, np.nan))

        desc_name, backbone, struct_formula = format_descriptive_name(
            segs, t_mrca=t_mrca, origin_basis=origin_basis, span_years=span_tup
        )
        barcode = format_hxb2_barcode(segs)

        if cl in pure_clades_set:
            legacy_name = f"Subtype_{cl}"
            category = "Nominal_Pure_Subtype"
        else:
            num = int(cl.split("_")[0])
            suffix = cl.split("_", 1)[1]
            legacy_name = f"CRF{num:02d}_{suffix}"
            category = "CRF_cpx" if suffix == "cpx" else ("CRF_2ndGen" if any(ch.isdigit() for ch in suffix) else "CRF_Primary")

        # Annotate taxonomic status under RhizAeon + ChronAeon rules
        if cl in ("A1", "A2", "B", "C", "D", "F1", "F2", "H", "K"):
            tax_status = "Validated_Basal_Trunk"
        elif cl == "G":
            tax_status = f"Demoted_Derived_Lineage (Postdates CRF02_AG: {clade_tmrca['G']:.1f} vs {clade_tmrca['02_AG']:.1f}; A1 gp41..nef + PT2695 A1 env)"
        elif cl == "J":
            tax_status = f"Demoted_Derived_Lineage (Postdates CRF11_cpx {clade_tmrca['11_cpx']:.1f} & CRF06_cpx {clade_tmrca['06_cpx']:.1f}: t_MRCA={clade_tmrca['J']:.1f})"
        elif cl == "L":
            tax_status = "Contaminated_Ref_Trunk (Ref 83CD003_Z3 carries K:nef 8976-9417)"
        elif cl in ("A3", "A4", "A6", "A7", "A8"):
            tax_status = "Sub-Subtype_Mosaic_or_Regional_Branch"
        elif legacy_name in out_of_panel_crfs and len(segs) == 1:
            tax_status = "Out_of_Panel_Donor (Unrepresented Basal Lineage A5/U)"
        elif len(segs) == 1:
            tax_status = "Collapsed_to_Clonal_Trunk (No Verified Mosaic Breakpoints)"
        else:
            tax_status = "Verified_Recombinant_Lineage"

        donors = sorted(set(s["trunk"] for s in segs if s["trunk"] != backbone))
        span_label = (
            f"{int(round(earliest_year[cl]))}-{int(round(latest_year[cl]))}"
            if cl in earliest_year and int(round(earliest_year[cl])) != int(round(latest_year[cl]))
            else (f"{int(round(earliest_year[cl]))}" if cl in earliest_year else "Undated")
        )
        compendium_rows.append({
            "legacy_label": legacy_name,
            "short_code": cl,
            "legacy_category": category,
            "taxonomic_status": tax_status,
            "n_ref_genomes": len(clade_seqs[cl]),
            "gamma_used": gamma_used,
            "descriptive_name": desc_name,
            "structural_formula": struct_formula,
            "backbone_trunk": backbone,
            "donor_trunks": ",".join(donors) if donors else "None",
            "n_verified_segments": len(segs),
            "n_verified_breakpoints": len(bps),
            "t_mrca": round(t_mrca, 1) if pd.notna(t_mrca) else None,
            "t_mrca_ci_low": round(ci_low, 1) if pd.notna(ci_low) else None,
            "t_mrca_ci_high": round(ci_high, 1) if pd.notna(ci_high) else None,
            "origin_basis": origin_basis,
            "sampling_span": span_label,
            "dating_method": t_src,
            "hxb2_barcode": barcode,
        })

        for s in segs:
            s_row = dict(s)
            s_row["legacy_label"] = legacy_name
            s_row["short_code"] = cl
            s_row["descriptive_name"] = desc_name
            all_deconvolved_segments.append(s_row)

        for b in bps:
            b_row = dict(b)
            b_row["legacy_label"] = legacy_name
            b_row["short_code"] = cl
            b_row["descriptive_name"] = desc_name
            all_deconvolved_bps.append(b_row)

    comp_df = pd.DataFrame(compendium_rows)
    comp_df.to_csv(RESULTS_DIR / "hiv1_descriptive_nomenclature_compendium.csv", index=False)
    pd.DataFrame(all_deconvolved_segments).to_csv(RESULTS_DIR / "hiv1_deconvolved_segments.csv", index=False)
    pd.DataFrame(all_deconvolved_bps).to_csv(RESULTS_DIR / "hiv1_deconvolved_breakpoints.csv", index=False)

    crf_only = comp_df[comp_df["legacy_category"] != "Nominal_Pure_Subtype"]
    print(f"\n[+] Generated Descriptive Spatio-Temporal Nomenclature for all {len(comp_df)} LANL Subtypes & CRFs ({len(crf_only)} CRFs):")
    for st_v, st_c in crf_only["taxonomic_status"].value_counts().items():
        print(f"    * {st_v:<62s}: {st_c:3d}")

    print("\n--- SAMPLE DESCRIPTIVE RE-CLASSIFICATIONS (Subtypes, Classic CRFs, and 2nd-Gen CRFs) ---")
    showcase = [
        "Subtype_A1", "Subtype_A2", "Subtype_A3", "Subtype_A6", "Subtype_A7", "Subtype_A8",
        "Subtype_B", "Subtype_C", "Subtype_D", "Subtype_F1", "Subtype_F2",
        "Subtype_G", "Subtype_H", "Subtype_J", "Subtype_K", "Subtype_L",
        "CRF01_AE", "CRF02_AG", "CRF03_A6B", "CRF04_cpx", "CRF05_DF", "CRF06_cpx",
        "CRF07_BC", "CRF08_BC", "CRF09_cpx", "CRF11_cpx", "CRF12_BF", "CRF13_cpx", "CRF15_01B",
        "CRF26_A5U", "CRF31_BC", "CRF33_01B", "CRF35_A1D", "CRF46_BF1", "CRF55_01B", "CRF63_02A6",
        "CRF92_C2U", "CRF146_BC", "CRF102_0107", "CRF104_0107"
    ]
    sub_df = comp_df[comp_df["legacy_label"].isin(showcase)]
    for _, r in sub_df.iterrows():
        print(f"  {r['legacy_label']:<14s} (t={str(r['t_mrca']):<6s}, g={r['gamma_used']}) -> {r['descriptive_name']}")
        print(f"      Status : {r['taxonomic_status']}")
        print(f"      Barcode: {r['hxb2_barcode']}")

    # Check redundant mosaic CRFs (n_verified_breakpoints >= 1) that collapse to identical structural formulas
    mosaic_crf_df = crf_only[crf_only["n_verified_breakpoints"] >= 1]
    formula_groups = mosaic_crf_df.groupby("structural_formula")["legacy_label"].apply(list)
    redundant = {k: v for k, v in formula_groups.items() if len(v) > 1}
    print(f"\n[+] Found {sum(len(v) for v in redundant.values())} mosaic CRFs (out of {len(mosaic_crf_df)} with >=1 bp) collapsing into {len(redundant)} identical gene-level structural formulas:")
    for form, crfs in sorted(redundant.items(), key=lambda x: -len(x[1])):
        print(f"    * {form:<48s} <- [{', '.join(crfs)}]")


if __name__ == "__main__":
    main()
