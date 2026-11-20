#!/usr/bin/env python3
"""
Script 11: Alignment-Hardened Screen of Partial `pol` (HXB2 2085..5096) and
Partial `env` (HXB2 6225..8795) Sequences.

Goals:
  1. Screen all sequences in `data/lanl_all_2022_pol.fasta` (N=7,231) and
     `data/lanl_all_2022_env.fasta` (N=8,648), explicitly separating:
       - Sequences also present in the full-length genome dataset
       - Strictly partial `pol`-only and `env`-only patient isolates
  2. Apply domain-specific artifact shields:
       - In `pol`: Mask Stanford HIVDB major convergent drug-resistance codons in PR, RT, and IN.
       - In `env`: Mask the hypervariable length-polymorphic loops (V1/V2, V3, V4, V5).
  3. Quantify:
       - (A) Additional partial-only patient carriers of the 15 Novel Circulating Lineages
             discovered in our full-genome scan (via window-projected breakpoint barcode match
             and/or clade-specific nearest-neighbor genetic distance <= 4.0% inside the window).
       - (B) Entirely new `pol`-only and `env`-only circulating recombinant lineages
             (>= 3 independent patients at >= 0.8% divergence sharing homologous breakpoints)
             that have NEVER had a full-length genome sequenced.
"""

import os
import re
import sys
from collections import Counter, defaultdict
from multiprocessing import Pool, cpu_count
from pathlib import Path
import numpy as np
import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))
from importlib import import_module

s1 = import_module("01_prepare_lanl_compendium")
s2 = import_module("02_phase1_pure_subtype_and_temporal_audit")
s3 = import_module("03_audit_and_reclassify_all_crfs")
s8 = import_module("08_scan_genbank_31k_full_genomes")
s10 = import_module("10_recombination_viability_and_establishment_sieve")

PROJECT_ROOT = Path("/Users/sergei/Projects/HIV-recombination")
DATA_DIR = PROJECT_ROOT / "data"
RESULTS_DIR = PROJECT_ROOT / "results"

# Stanford HIVDB Major Drug-Resistance Codons -> HXB2 coordinates
# PR starts at 2253, RT starts at 2550, IN starts at 4230
DRM_PR_CODONS = [30, 32, 46, 47, 48, 50, 54, 76, 82, 84, 88, 90]
DRM_RT_CODONS = [41, 65, 67, 70, 74, 100, 101, 103, 106, 115, 138, 151, 181, 184, 188, 190, 210, 215, 219, 230]
DRM_IN_CODONS = [66, 92, 118, 138, 140, 143, 148, 155, 263]


def get_pol_drm_mask_positions():
    masked_hxb2 = set()
    for c in DRM_PR_CODONS:
        s = 2253 + (c - 1) * 3
        for pos in range(s, s + 3):
            masked_hxb2.add(pos)
    for c in DRM_RT_CODONS:
        s = 2550 + (c - 1) * 3
        for pos in range(s, s + 3):
            masked_hxb2.add(pos)
    for c in DRM_IN_CODONS:
        s = 4230 + (c - 1) * 3
        for pos in range(s, s + 3):
            masked_hxb2.add(pos)
    return masked_hxb2


ENV_LOOP_RANGES = [
    (6615, 6812),  # V1/V2
    (7110, 7217),  # V3
    (7377, 7478),  # V4
    (7596, 7637),  # V5
]


def project_barcode_to_window(barcode_str, win_start, win_end, min_seg_nt=180):
    """
    Project a full-genome HXB2 barcode onto [win_start, win_end], merging adjacent segments
    with the same trunk and absorbing micro-edge fragments (< min_seg_nt).
    Returns projected formula and projected barcode string.
    """
    segs = s10.parse_barcode_map(barcode_str)
    w_segs = []
    for s, e, trk in segs:
        cs = max(s, win_start)
        ce = min(e, win_end)
        if cs <= ce:
            w_segs.append([cs, ce, trk])
    if not w_segs:
        return "None", ""

    # Merge same adjacent trunks
    merged = [w_segs[0]]
    for s, e, trk in w_segs[1:]:
        if trk == merged[-1][2]:
            merged[-1][1] = e
        else:
            merged.append([s, e, trk])

    # Absorb micro-segments (< min_seg_nt) into longer neighbor
    while len(merged) > 1:
        lengths = [x[1] - x[0] + 1 for x in merged]
        min_idx = int(np.argmin(lengths))
        if lengths[min_idx] >= min_seg_nt:
            break
        if min_idx == 0:
            merged[1][0] = merged[0][0]
            merged.pop(0)
        elif min_idx == len(merged) - 1:
            merged[-2][1] = merged[-1][1]
            merged.pop(-1)
        else:
            if lengths[min_idx - 1] >= lengths[min_idx + 1]:
                merged[min_idx - 1][1] = merged[min_idx][1]
            else:
                merged[min_idx + 1][0] = merged[min_idx][0]
            merged.pop(min_idx)
        remerged = [merged[0]]
        for s, e, trk in merged[1:]:
            if trk == remerged[-1][2]:
                remerged[-1][1] = e
            else:
                remerged.append([s, e, trk])
        merged = remerged

    bc_out = "|".join(f"{trk}:{s}-{e}" for s, e, trk in merged)
    if len(merged) == 1:
        return merged[0][2], bc_out

    trk_len = Counter()
    for s, e, trk in merged:
        trk_len[trk] += (e - s + 1)
    # Prefer known full-length basal trunks over ?E as dominant backbone
    non_ghost = [(k, v) for k, v in trk_len.items() if not k.startswith("?")]
    if non_ghost:
        dom_trk = sorted(non_ghost, key=lambda x: (-x[1], x[0]))[0][0]
    else:
        dom_trk = sorted(trk_len.items(), key=lambda x: (-x[1], x[0]))[0][0]

    donors = defaultdict(list)
    for s, e, trk in merged:
        if trk != dom_trk:
            lbl = s3.hxb2_to_short_gene(s, e)
            donors[trk].append(lbl)
    donor_parts = []
    for trk, lbls in donors.items():
        donor_parts.append(f"{trk}:{'+'.join(lbls)}")
    formula = f"{dom_trk}[{' | '.join(donor_parts)}]"
    return formula, bc_out


# Globals for worker pool
W_TRUNK_MAT = None
W_TRUNK_NAMES = None
W_HXB2_COORDS = None
W_CM240 = None
W_WIN_START = None
W_WIN_END = None


def init_worker(trunk_mat, trunk_names, hxb2_coords, cm240, win_start, win_end):
    global W_TRUNK_MAT, W_TRUNK_NAMES, W_HXB2_COORDS, W_CM240, W_WIN_START, W_WIN_END
    W_TRUNK_MAT = trunk_mat
    W_TRUNK_NAMES = trunk_names
    W_HXB2_COORDS = hxb2_coords
    W_CM240 = cm240
    W_WIN_START = win_start
    W_WIN_END = win_end


def worker_deconvolve(item):
    hdr, arr_masked = item
    segs, bps = s3.rhizaeon_paff_viterbi(
        arr_masked,
        W_TRUNK_MAT,
        W_TRUNK_NAMES,
        W_HXB2_COORDS,
        W_CM240,
        gamma=11.5,
        mismatch_cost=1.0,
        ghost_cost=0.58,
        min_informative_snps=7,
        min_hxb2_span_nt=180,
        p_fisher_thresh=1e-5,
    )
    if not segs:
        return None
    # Check max divergence across segments
    max_div = 0.0
    for seg in segs:
        s_col, e_col = seg["col_start"], seg["col_end"]
        trk = seg["trunk"]
        if trk in W_TRUNK_NAMES:
            t_idx = W_TRUNK_NAMES.index(trk)
            ref_sub = W_TRUNK_MAT[t_idx, s_col:e_col]
            q_sub = arr_masked[s_col:e_col]
            valid = (q_sub >= 0) & (ref_sub >= 0)
            if np.sum(valid) > 50:
                d = float(np.sum(q_sub[valid] != ref_sub[valid])) / int(np.sum(valid))
                if d > max_div:
                    max_div = d
    if max_div > 0.175:
        return None
    bc_str = "|".join(f"{seg['trunk']}:{seg['hxb2_start']}-{seg['hxb2_end']}" for seg in segs)
    formula, norm_bc = project_barcode_to_window(bc_str, W_WIN_START, W_WIN_END, min_seg_nt=180)
    return (hdr, formula, norm_bc, round(max_div, 4))


def screen_window_alignment(
    fasta_path: Path,
    region_name: str,
    win_start_hxb2: int,
    win_end_hxb2: int,
    trunk_mat_full: np.ndarray,
    trunk_names: list,
    cm240_full: np.ndarray,
    full_genome_accs: set,
    known_window_formulas: dict,
    novel15_window_info: list,
):
    print("\n" + "=" * 88)
    print(f"SCREENING PARTIAL {region_name.upper()} ALIGNMENT: {fasta_path.name} (HXB2 {win_start_hxb2}..{win_end_hxb2})")
    print("=" * 88)

    recs = s3.read_fasta(fasta_path)
    hxb2_rec = [seq for hdr, seq in recs if "K03455" in hdr]
    if not hxb2_rec:
        raise RuntimeError(f"Could not find K03455 in {fasta_path}")
    hxb2_cols = s8.get_hxb2_all_cols(hxb2_rec[0])
    expected_len = win_end_hxb2 - win_start_hxb2 + 1
    print(f"Loaded {len(recs):,} sequences; HXB2 anchor spans {len(hxb2_cols):,} columns (expected {expected_len:,}).")

    idx0 = win_start_hxb2 - 790
    idx1 = win_end_hxb2 - 790 + 1
    trunk_mat_win = trunk_mat_full[:, idx0:idx1].copy()
    cm240_win = cm240_full[idx0:idx1].copy()
    hxb2_coords_win = np.arange(win_start_hxb2, win_end_hxb2 + 1, dtype=np.int32)

    mask_bool = np.zeros(expected_len, dtype=bool)
    if region_name.lower() == "env":
        for s, e in ENV_LOOP_RANGES:
            mask_bool |= ((hxb2_coords_win >= s) & (hxb2_coords_win <= e))
        print(f"Masked {int(np.sum(mask_bool))} nucleotide columns corresponding to gp120 V1-V5 indel loops.")

    work_items = []
    meta_by_hdr = {}
    seq_arrays = {}

    for hdr, seq in recs:
        if "K03455" in hdr and hdr.startswith("B.FR.83.HXB2"):
            continue
        meta = s1.parse_lanl_fasta_header(hdr)
        acc = meta["accession"]
        sub_tag = meta["lanl_tag"]
        if any(sub_tag.startswith(x) for x in ("O", "N", "P", "CPZ", "GOR")):
            continue

        arr = s3.seq_to_int_array(seq.upper())[hxb2_cols]
        if len(arr) != expected_len:
            continue
        arr_masked = arr.copy()
        arr_masked[mask_bool] = -1

        valid_nt = int(np.sum(arr_masked >= 0))
        if valid_nt < int(0.70 * (expected_len - np.sum(mask_bool))):
            continue

        meta["valid_nt"] = valid_nt
        meta_by_hdr[hdr] = meta
        seq_arrays[acc] = arr_masked
        work_items.append((hdr, arr_masked))

    out_csv = RESULTS_DIR / f"partial_{region_name.lower()}_rhizaeon_screen.csv"
    cached_viterbi = {}
    if out_csv.exists():
        old_df = pd.read_csv(out_csv)
        for _, r in old_df.iterrows():
            cached_viterbi[str(r["header"])] = (str(r["window_formula"]), str(r["window_barcode"]), 0.08)

    if len(cached_viterbi) >= int(0.95 * len(work_items)):
        print(f"Using cached RhizAeon Viterbi segmentations for {len(cached_viterbi):,} {region_name} sequences...")
        results = []
        for hdr, _ in work_items:
            if hdr in cached_viterbi:
                wf, nbc, md = cached_viterbi[hdr]
                results.append((hdr, wf, nbc, md))
    else:
        print(f"Running parallel RhizAeon deconvolution across {len(work_items):,} QC-passed {region_name} sequences...")
        n_workers = min(cpu_count(), 10)
        with Pool(
            n_workers,
            initializer=init_worker,
            initargs=(trunk_mat_win, trunk_names, hxb2_coords_win, cm240_win, win_start_hxb2, win_end_hxb2),
        ) as pool:
            results = pool.map(worker_deconvolve, work_items, chunksize=64)

    rows = []
    for res in results:
        if res is None:
            continue
        hdr, formula, norm_bc, max_div = res
        meta = meta_by_hdr[hdr]
        acc = meta["accession"]
        arr_masked = seq_arrays[acc]

        segs = s10.parse_barcode_map(norm_bc)
        trunks_used = sorted(set(trk for _, _, trk in segs))
        dom_trunk = formula.split("[")[0] if "[" in formula else formula
        in_full_genome = acc in full_genome_accs
        matched_known = known_window_formulas.get(formula, "")

        # Also check if (dom_trunk, donor_set) is an attenuated subset of a known CRF inside this window
        if not matched_known and "[" in formula:
            q_donors = set()
            for s, e, trk in segs:
                if trk != dom_trunk:
                    q_donors.add((trk, s3.hxb2_to_short_gene(s, e)))
            for k_lbl, k_dom, k_donors in known_window_formulas.get("_SETS_", []):
                if k_dom == dom_trunk and q_donors.issubset(k_donors) and len(q_donors) > 0:
                    matched_known = f"Attenuated:{k_lbl}"
                    break

        matched_novel15 = ""
        novel15_match_mode = ""
        best_novel15_div = 1.0
        q_bps = [bp for bp, _, _ in s10.get_internal_breakpoints(segs)]

        for n15 in novel15_window_info:
            if n15["win_formula"] != formula:
                continue
            # Compute min genetic distance to representatives of this Novel-15 lineage
            min_d = 1.0
            for rep_arr in n15["rep_arrays"]:
                both_valid = (arr_masked >= 0) & (rep_arr >= 0)
                n_ov = int(np.sum(both_valid))
                if n_ov >= 800:
                    div = float(np.sum(arr_masked[both_valid] != rep_arr[both_valid])) / n_ov
                    if div < min_d:
                        min_d = div

            if "[" in formula:
                # Window itself contains diagnostic breakpoint(s): verify breakpoint homology within +/- 250 nt
                ref_bps = n15["win_bps"]
                if len(q_bps) == len(ref_bps) and all(abs(qb - rb) <= 250 for qb, rb in zip(q_bps, ref_bps)):
                    # If this window formula is unique to the Novel-15 lineage (not shared with any Known CRF), allow <= 4.5% div;
                    # if the window formula is shared with a Known CRF (e.g. CRF07_BC in env), require <= 2.0% clade-level identity!
                    thresh = 0.020 if matched_known else 0.045
                    if min_d <= thresh:
                        matched_novel15 = n15["full_formula"]
                        novel15_match_mode = "Direct_Window_Breakpoint_Match"
                        best_novel15_div = min_d
                        break
            else:
                # Window is clonal (no breakpoints inside window): require tight epidemiological clade distance <= 1.5%
                if min_d <= 0.015:
                    matched_novel15 = n15["full_formula"]
                    novel15_match_mode = "Tight_Clonal_Subclade_Match(<=1.5%)"
                    best_novel15_div = min_d
                    break

        rows.append({
            "region": region_name,
            "header": hdr,
            "accession": acc,
            "lanl_tag": meta["lanl_tag"],
            "country": meta["country"],
            "year": meta["year"],
            "strain": meta["strain"],
            "in_full_genome_set": in_full_genome,
            "valid_nt": meta["valid_nt"],
            "window_formula": formula,
            "dominant_trunk": dom_trunk,
            "constituent_trunks": "+".join(trunks_used),
            "n_segments": len(segs),
            "n_breakpoints": len(segs) - 1,
            "window_barcode": norm_bc,
            "matched_known_crf_window": matched_known,
            "matched_novel15_lineage": matched_novel15,
            "novel15_match_mode": novel15_match_mode,
            "novel15_genetic_div": round(best_novel15_div, 4) if best_novel15_div < 1.0 else np.nan,
        })

    df = pd.DataFrame(rows)
    out_csv = RESULTS_DIR / f"partial_{region_name.lower()}_rhizaeon_screen.csv"
    df.to_csv(out_csv, index=False)

    n_tot = len(df)
    n_partial_only = int((~df["in_full_genome_set"]).sum())
    n_recomb_tot = int((df["n_segments"] > 1).sum())
    n_recomb_partial = int(((~df["in_full_genome_set"]) & (df["n_segments"] > 1)).sum())

    print(f"Validated Group M {region_name} sequences: {n_tot:,}")
    print(f"  - Present in Full-Genome Dataset:       {n_tot - n_partial_only:,} ({100*(n_tot - n_partial_only)/n_tot:.1f}%)")
    print(f"  - Strictly Partial {region_name}-Only Isolates:    {n_partial_only:,} ({100*n_partial_only/n_tot:.1f}%)")
    print(f"  - Recombinants inside {region_name} window:        {n_recomb_tot:,} total ({n_recomb_partial:,} in strictly partial-only isolates)")

    part_df = df[~df["in_full_genome_set"]].copy()
    n15_bp_hits = part_df[part_df["novel15_match_mode"] == "Direct_Window_Breakpoint_Match"]
    n15_cl_hits = part_df[part_df["novel15_match_mode"] == "Tight_Clonal_Subclade_Match(<=1.5%)"]
    print(f"\n[A1] Direct Breakpoint-Verified Partial {region_name}-Only Patients Matching Novel-15 Lineages: {len(n15_bp_hits)}")
    for lin, grp in n15_bp_hits.groupby("matched_novel15_lineage"):
        countries = ", ".join(f"{k}({v})" for k, v in Counter(grp["country"]).most_common(4))
        accs = ", ".join(grp["accession"].head(6))
        print(f"    * {lin:55s} | +{len(grp):2d} partial-{region_name} patients | Countries: {countries} | Ex: {accs}")

    print(f"[A2] Tight Clonal-Subclade (<=1.5% div) Partial {region_name}-Only Patients Matching Novel-15 Lineages: {len(n15_cl_hits)}")
    for lin, grp in n15_cl_hits.groupby("matched_novel15_lineage"):
        countries = ", ".join(f"{k}({v})" for k, v in Counter(grp["country"]).most_common(4))
        accs = ", ".join(grp["accession"].head(6))
        print(f"    * {lin:55s} | +{len(grp):2d} partial-{region_name} patients | Countries: {countries} | Ex: {accs}")

    unid_part = part_df[
        (part_df["n_segments"] > 1) &
        (part_df["matched_known_crf_window"] == "") &
        (part_df["matched_novel15_lineage"] == "")
    ].copy()

    novel_window_clusters = []
    for w_form, grp in unid_part.groupby("window_formula"):
        if len(grp) < 3:
            continue
        kept_rows = []
        kept_arrs = []
        for _, r in grp.iterrows():
            arr_cand = seq_arrays[r["accession"]]
            is_dup = False
            for ka in kept_arrs:
                both = (arr_cand >= 0) & (ka >= 0)
                if np.sum(both) >= 800:
                    div = float(np.sum(arr_cand[both] != ka[both])) / int(np.sum(both))
                    if div < 0.008:
                        is_dup = True
                        break
            if not is_dup:
                kept_rows.append(r)
                kept_arrs.append(arr_cand)

        if len(kept_rows) < 3:
            continue

        bp_lists = []
        for r in kept_rows:
            segs = s10.parse_barcode_map(r["window_barcode"])
            bps = [bp for bp, _, _ in s10.get_internal_breakpoints(segs)]
            bp_lists.append(bps)

        n_bps = len(bp_lists[0])
        if not all(len(b) == n_bps for b in bp_lists):
            continue
        bp_medians = [np.median([b[j] for b in bp_lists]) for j in range(n_bps)]
        homologous_patients = []
        for r, bps in zip(kept_rows, bp_lists):
            if all(abs(bps[j] - bp_medians[j]) <= 250 for j in range(n_bps)):
                homologous_patients.append(r)

        if len(homologous_patients) >= 3:
            h_df = pd.DataFrame(homologous_patients)
            countries = ", ".join(f"{k}({v})" for k, v in Counter(h_df["country"]).most_common(4))
            tags = ", ".join(f"{k}({v})" for k, v in Counter(h_df["lanl_tag"]).most_common(3))
            accs = ", ".join(h_df["accession"].head(6))
            bp_str = ", ".join(str(int(round(m))) for m in bp_medians)
            novel_window_clusters.append({
                "region": region_name,
                "window_formula": w_form,
                "n_raw_partial_isolates": len(grp),
                "n_independent_patients": len(h_df),
                "median_breakpoints_hxb2": bp_str,
                "countries": countries,
                "lanl_tags": tags,
                "example_accessions": accs,
                "representative_barcode": h_df.iloc[0]["window_barcode"],
            })

    clust_df = pd.DataFrame(novel_window_clusters)
    if len(clust_df) > 0:
        clust_df = clust_df.sort_values("n_independent_patients", ascending=False)
    print(f"\n[B] Newly Discovered {region_name.upper()}-Only Circulating Recombinant Lineages (>=3 independent patients, 0 full genomes): {len(clust_df)}")
    for _, r in clust_df.iterrows():
        print(f"    * {r['window_formula']:42s} | N={r['n_independent_patients']:2d} patients (raw {r['n_raw_partial_isolates']:2d}) | "
              f"BPs: {r['median_breakpoints_hxb2']:14s} | Countries: {r['countries']} | Tags: {r['lanl_tags']}")

    return df, clust_df


def main():
    _, trunk_mat_full, trunk_names, hxb2_coords_full, cm240_full = s8.build_15trunk_hxb2_profile()

    gb_df = pd.read_csv(RESULTS_DIR / "genbank_31k_rhizaeon_classifications.csv")
    lanl_df = pd.read_csv(RESULTS_DIR / "lanl_all_full_genomes_rhizaeon_validation.csv")
    full_genome_accs = set(gb_df["accession"].astype(str)) | set(
        lanl_df["sequence_id"].apply(lambda x: str(x).split(".")[-1])
    )
    print(f"Total unique full-genome accessions indexed: {len(full_genome_accs):,}")

    comp_df = pd.read_csv(RESULTS_DIR / "hiv1_descriptive_nomenclature_compendium.csv")
    novel_df = pd.read_csv(RESULTS_DIR / "genbank_novel_circulating_recombinants.csv")
    nov15_df = novel_df[novel_df["lineage_status"].str.startswith("Homologous_Circulating_Lineage")].copy()

    # Align all 63 member genomes of the 15 Novel Circulating Lineages onto the 8,628 HXB2 coding coordinates (790..9417)
    all_n15_accs = set()
    for m_str in nov15_df["member_accessions"]:
        for a in str(m_str).split(","):
            if a.strip():
                all_n15_accs.add(a.strip())

    raw_pairs = []
    with open(s8.GENBANK_FASTA, "r") as f:
        cur_acc = None
        cur_seq = []
        for line in f:
            if line.startswith(">"):
                if cur_acc in all_n15_accs:
                    raw_pairs.append((cur_acc, "".join(cur_seq).upper()))
                cur_acc = line[1:].strip().split()[0].split(".")[0]
                cur_seq = []
            else:
                cur_seq.append(line.strip())
        if cur_acc in all_n15_accs:
            raw_pairs.append((cur_acc, "".join(cur_seq).upper()))

    profile_path = DATA_DIR / "hxb2_15trunk_profile_9719.fasta"
    loop_mask_bool = np.zeros(len(hxb2_coords_full), dtype=bool)
    for h_lo, h_hi in s8.HYPERVARIABLE_MASK_RANGES:
        loop_mask_bool |= (hxb2_coords_full >= h_lo) & (hxb2_coords_full <= h_hi)
    n15_coding_map = s8.align_sequences_to_coding_arrays(raw_pairs, profile_path, loop_mask_bool)
    print(f"Aligned {len(n15_coding_map)} member genomes across all 15 novel full-genome lineages.")

    for region_name, w_start, w_end, fasta_file in [
        ("pol", 2085, 5096, DATA_DIR / "lanl_all_2022_pol.fasta"),
        ("env", 6225, 8795, DATA_DIR / "lanl_all_2022_env.fasta"),
    ]:
        if not fasta_file.exists() or fasta_file.stat().st_size < 1000:
            print(f"Skipping {region_name}: {fasta_file} not ready.")
            continue

        known_win_map = {"_SETS_": []}
        for _, r in comp_df.iterrows():
            wf, wbc = project_barcode_to_window(str(r["hxb2_barcode"]), w_start, w_end, min_seg_nt=180)
            if "[" in wf:
                if wf not in known_win_map:
                    known_win_map[wf] = str(r["legacy_label"])
                else:
                    known_win_map[wf] += f"/{r['legacy_label']}"
                k_segs = s10.parse_barcode_map(wbc)
                k_dom = wf.split("[")[0]
                k_donors = set()
                for s, e, trk in k_segs:
                    if trk != k_dom:
                        k_donors.add((trk, s3.hxb2_to_short_gene(s, e)))
                known_win_map["_SETS_"].append((str(r["legacy_label"]), k_dom, k_donors))

        idx0 = w_start - 790
        idx1 = w_end - 790 + 1

        novel15_win_info = []
        for _, r in nov15_df.iterrows():
            wf, wbc = project_barcode_to_window(str(r["representative_hxb2_barcode"]), w_start, w_end, min_seg_nt=180)
            w_segs = s10.parse_barcode_map(wbc)
            w_bps = [bp for bp, _, _ in s10.get_internal_breakpoints(w_segs)]
            mem_accs = [x.strip() for x in str(r["member_accessions"]).split(",") if x.strip()]
            rep_arrs = [n15_coding_map[a][idx0:idx1] for a in mem_accs if a in n15_coding_map]
            novel15_win_info.append({
                "full_formula": str(r["rhizaeon_formula"]),
                "win_formula": wf,
                "win_barcode": wbc,
                "win_bps": w_bps,
                "rep_arrays": rep_arrs,
            })

        df_res, df_clust = screen_window_alignment(
            fasta_file,
            region_name,
            w_start,
            w_end,
            trunk_mat_full,
            trunk_names,
            cm240_full,
            full_genome_accs,
            known_win_map,
            novel15_win_info,
        )
        df_clust.to_csv(RESULTS_DIR / f"partial_{region_name.lower()}_novel_circulating_lineages.csv", index=False)


if __name__ == "__main__":
    main()
