#!/usr/bin/env python3
"""
19_chronaeon_lsd2_treetime_validation.py
========================================
Tree-based molecular clock validation of ChronAeon (OLS / Restricted Cubic Spline +
country x 5-yr block bootstrap) against IQ-TREE 2 (--date / LSD2 v2.4) and TreeTime
(v0.12) across:
  Part A: All 21 HIV-1 Group M subtypes & CRFs from `phase1_chronaeon_clade_dating.csv`
          plus longest non-recombinant segment sensitivity checks for 6 major CRFs
          (`01_AE`, `02_AG`, `06_cpx`, `07_BC`, `08_BC`, `12_BF`).
          Uses deterministic temporal-stratified subsampling (<= 150 tips) for clades
          with N > 150 and evaluates ChronAeon, LSD2, and TreeTime on identical taxa,
          sites, and sampling years.
  Part B: All 15 Homologous_Circulating_Lineage rows (`NOVEL_CRF_01..15`) using their
          longest non-recombinant segment (> 1,000 nt) and applying the strict
          qualification rule (`n_dated >= 4`, `sample_span >= 3.0 yr`, positive clock
          slope) to separate `@tMRCA` (`Clock_tMRCA`) from `[Span:YYYY-YYYY]`
          (`Sampling_Span_Only`) and `Undated`.

Heavy IQ-TREE 2 + LSD2 + TreeTime execution is dispatched to `m3.local` (32-core Mac
Studio) via SSH/tarball and results are saved locally under `results/`.
"""

import glob
import json
import math
import re
import subprocess
import sys
import tarfile
import tempfile
import time
from importlib import import_module
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import Phylo, SeqIO
from scipy import stats

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent
DATA_DIR = PROJECT_ROOT / "data"
RESULTS_DIR = PROJECT_ROOT / "results"
TREES_DIR = RESULTS_DIR / "trees"

sys.path.insert(0, str(SCRIPT_DIR))
s2 = import_module("02_phase1_pure_subtype_and_temporal_audit")

REMOTE_HOST = "m3.local"
REMOTE_DIR = "/Users/sergei/hiv_recomb_bench/results/clock_validation"
REMOTE_IQTREE = "/Users/sergei/hiv_recomb_bench/bin/iqtree2"
REMOTE_TREETIME = "/Users/sergei/epistasis_env/bin/treetime"


def load_fasta_dict(path: Path) -> dict:
    seqs = {}
    for rec in SeqIO.parse(str(path), "fasta"):
        seqs[rec.id] = str(rec.seq).upper()
    return seqs


def get_hxb2_cols(seqs: dict, hxb2_sub: str = "K03455") -> dict:
    """Map 1-indexed HXB2 nucleotide positions to 0-indexed MSA columns."""
    hxb2_id = next((k for k in seqs if hxb2_sub in k), None)
    if hxb2_id is None:
        raise RuntimeError(f"HXB2 ({hxb2_sub}) not found in FASTA.")
    hxb2_seq = seqs[hxb2_id]
    h_to_col = {}
    h_pos = 0
    for col_0, ch in enumerate(hxb2_seq):
        if ch not in "-.":
            h_pos += 1
            h_to_col[h_pos] = col_0
    return h_to_col


def extract_filtered_alignment(seq_ids: list, raw_seqs: dict, cols: list, max_gap_frac: float = 0.50) -> dict:
    """Slice MSA columns and retain positions with <= max_gap_frac gaps across seq_ids."""
    n_seq = len(seq_ids)
    sliced = {sid: "".join(raw_seqs[sid][c] for c in cols) for sid in seq_ids}
    arr = np.array([list(sliced[sid]) for sid in seq_ids], dtype="S1")
    is_gap = (arr == b"-") | (arr == b".")
    keep_mask = is_gap.mean(axis=0) <= max_gap_frac
    kept_idx = np.where(keep_mask)[0]
    arr_kept = arr[:, kept_idx]
    return {sid: arr_kept[i].tobytes().decode("ascii") for i, sid in enumerate(seq_ids)}


def temporal_stratified_subsample(df: pd.DataFrame, max_tips: int = 150, seed: int = 20261114) -> pd.DataFrame:
    """Deterministic temporal-stratified subsample preserving oldest/newest anchors and year spread."""
    if len(df) <= max_tips:
        return df.copy().reset_index(drop=True)
    rng = np.random.default_rng(seed)
    df_sorted = df.sort_values(["year", "seq_id"]).reset_index(drop=True)
    # Allocate counts proportionally across year strata while ensuring every year bin with >=1 sample gets >=1 if possible
    years = df_sorted["year"].astype(int).values
    unique_years, counts = np.unique(years, return_counts=True)
    if len(unique_years) >= max_tips:
        # Spread evenly across sorted index
        pick_idx = np.linspace(0, len(df_sorted) - 1, max_tips, dtype=int)
        return df_sorted.iloc[pick_idx].reset_index(drop=True)

    # First reserve 1 deterministic representative per unique year
    selected_indices = []
    remaining_indices = []
    for y in unique_years:
        y_idx = np.where(years == y)[0]
        # Deterministic pick centered in stratum
        chosen = y_idx[0]
        selected_indices.append(int(chosen))
        for idx in y_idx[1:]:
            remaining_indices.append(int(idx))

    n_need = max_tips - len(selected_indices)
    if n_need > 0 and remaining_indices:
        # Stratified draw weighted evenly across remaining sorted order
        rem = np.array(sorted(remaining_indices), dtype=int)
        if len(rem) <= n_need:
            selected_indices.extend(rem.tolist())
        else:
            step_picks = np.unique(np.linspace(0, len(rem) - 1, n_need, dtype=int))
            selected_indices.extend(rem[step_picks].tolist())
            # Fill any rounding deficit deterministically with fixed-seed shuffle of leftover
            if len(selected_indices) < max_tips:
                picked_set = set(selected_indices)
                leftover = [int(i) for i in rem if int(i) not in picked_set]
                rng.shuffle(leftover)
                selected_indices.extend(leftover[: max_tips - len(selected_indices)])

    selected_indices = sorted(selected_indices[:max_tips])
    return df_sorted.iloc[selected_indices].reset_index(drop=True)


def run_remote_clock_jobs(job_specs: list) -> dict:
    """
    Package all alignments & date files, push to m3.local, execute IQ-TREE2 + LSD2 + TreeTime
    in parallel on m3.local, and fetch back a single summary JSON + outputs tarball.
    """
    with tempfile.TemporaryDirectory(prefix="chronaeon_val_") as tmp_dir:
        tmp_path = Path(tmp_dir)
        jobs_dir = tmp_path / "jobs"
        jobs_dir.mkdir(parents=True, exist_ok=True)

        manifest = []
        for spec in job_specs:
            jid = spec["job_id"]
            jdir = jobs_dir / jid
            jdir.mkdir(parents=True, exist_ok=True)
            fa_path = jdir / "aln.fasta"
            lsd_dates = jdir / "dates.tsv"
            tt_dates = jdir / "dates_tt.csv"
            with open(fa_path, "w") as ff, open(lsd_dates, "w") as lf, open(tt_dates, "w") as tf:
                tf.write("name,date\n")
                for sid, seq_str in spec["seqs"].items():
                    yr = float(spec["dates"][sid])
                    ff.write(f">{sid}\n{seq_str}\n")
                    lf.write(f"{sid}\t{yr:.2f}\n")
                    tf.write(f"{sid},{yr:.2f}\n")
            manifest.append({"job_id": jid, "n_taxa": len(spec["seqs"])})

        with open(tmp_path / "manifest.json", "w") as f:
            json.dump(manifest, f, indent=2)

        remote_runner = tmp_path / "run_remote_clock.py"
        remote_runner.write_text(f"""#!/usr/bin/env python3
import json, re, subprocess, time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

ROOT = Path("{REMOTE_DIR}")
IQTREE = "{REMOTE_IQTREE}"
TREETIME = "{REMOTE_TREETIME}"
manifest = json.loads((ROOT / "manifest.json").read_text())

def run_one(item):
    jid = item["job_id"]
    jdir = ROOT / "jobs" / jid
    fa = jdir / "aln.fasta"
    dt_lsd = jdir / "dates.tsv"
    dt_tt = jdir / "dates_tt.csv"
    pref = jdir / "iq"
    tt_dir = jdir / "tt_out"

    t0 = time.time()
    out = {{"job_id": jid, "lsd2_status": "NOT_RUN", "treetime_status": "NOT_RUN"}}

    # 1. IQ-TREE 2 ML + LSD2 dating
    cmd_iq = [
        IQTREE, "-s", str(fa), "-m", "GTR+G4", "--fast",
        "--date", str(dt_lsd), "--date-ci", "100", "--date-options", "-l 0",
        "-T", "2", "-redo", "--prefix", str(pref)
    ]
    res_iq = subprocess.run(cmd_iq, capture_output=True, text=True)
    lsd_file = Path(str(pref) + ".timetree.lsd")
    treefile = Path(str(pref) + ".treefile")
    if res_iq.returncode == 0 and lsd_file.exists():
        txt = lsd_file.read_text()
        out["lsd2_raw_text"] = txt
        m_ci = re.search(
            r"rate\\s+([0-9.eE+-]+)\\s*\\[([0-9.eE+-]+);\\s*([0-9.eE+-]+)\\],\\s*tMRCA\\s+([0-9.eE+-]+)\\s*\\[([^;]+);\\s*([^\\]]+)\\]",
            txt
        )
        m_pt = re.search(r"-\\s*Dating results:\\s*\\n\\s*rate\\s+([0-9.eE+-]+),\\s*tMRCA\\s+([0-9.eE+-]+)", txt)
        if m_ci:
            out["lsd2_rate"] = float(m_ci.group(1))
            out["lsd2_rate_lo"] = float(m_ci.group(2))
            out["lsd2_rate_hi"] = float(m_ci.group(3))
            out["lsd2_tmrca"] = float(m_ci.group(4))
            try:
                out["lsd2_ci_lo"] = float(m_ci.group(5))
            except ValueError:
                out["lsd2_ci_lo"] = None
            try:
                out["lsd2_ci_hi"] = float(m_ci.group(6))
            except ValueError:
                out["lsd2_ci_hi"] = None
            out["lsd2_status"] = "OK"
        elif m_pt:
            out["lsd2_rate"] = float(m_pt.group(1))
            out["lsd2_tmrca"] = float(m_pt.group(2))
            out["lsd2_status"] = "OK_NO_CI"
        else:
            out["lsd2_status"] = "PARSE_ERR"
    else:
        out["lsd2_status"] = f"FAILED_RC_{{res_iq.returncode}}"
        out["lsd2_err"] = (res_iq.stderr or res_iq.stdout)[-400:]

    # 2. TreeTime on IQ-TREE ML tree
    if treefile.exists():
        cmd_tt = [
            TREETIME, "--tree", str(treefile), "--aln", str(fa),
            "--dates", str(dt_tt), "--clock-filter", "0", "--confidence",
            "--outdir", str(tt_dir)
        ]
        res_tt = subprocess.run(cmd_tt, capture_output=True, text=True)
        dates_out = tt_dir / "dates.tsv"
        clock_out = tt_dir / "molecular_clock.txt"
        tt_nex = tt_dir / "timetree.nexus"
        if res_tt.returncode == 0 and dates_out.exists():
            # Parse root node name from nexus tree
            root_name = None
            if tt_nex.exists():
                nex_txt = tt_nex.read_text()
                # Last internal node before trailing ';' in NEXUS tree is root
                m_root = re.search(r"\\)([A-Za-z0-9_]+)(?::[0-9.eE+-]+)?(?:\\[[^\\]]*\\])?;", nex_txt)
                if m_root:
                    root_name = m_root.group(1)
            node_dates = {{}}
            for line in dates_out.read_text().splitlines():
                if not line.strip() or line.startswith("#"):
                    continue
                parts = line.split("\\t")
                if len(parts) >= 3 and parts[0].startswith("NODE_"):
                    try:
                        node_dates[parts[0]] = float(parts[2])
                    except ValueError:
                        pass
            if root_name and root_name in node_dates:
                out["tt_tmrca"] = node_dates[root_name]
                out["tt_root_node"] = root_name
            elif node_dates:
                oldest = min(node_dates.items(), key=lambda kv: kv[1])
                out["tt_tmrca"] = oldest[1]
                out["tt_root_node"] = oldest[0]
            if clock_out.exists():
                clk_txt = clock_out.read_text()
                m_r = re.search(r"--rate:\\s*([0-9.eE+-]+)", clk_txt)
                m_r2 = re.search(r"--r\\^2:\\s*([0-9.eE+-]+)", clk_txt)
                if m_r: out["tt_rate"] = float(m_r.group(1))
                if m_r2: out["tt_r2"] = float(m_r2.group(1))
            out["treetime_status"] = "OK" if "tt_tmrca" in out else "PARSE_ERR"
        else:
            out["treetime_status"] = f"FAILED_RC_{{res_tt.returncode}}"
            out["tt_err"] = (res_tt.stderr or res_tt.stdout)[-400:]

    out["wall_sec"] = round(time.time() - t0, 2)
    return out

t_start = time.time()
results = {{}}
with ThreadPoolExecutor(max_workers=10) as ex:
    futs = {{ex.submit(run_one, it): it["job_id"] for it in manifest}}
    for fut in as_completed(futs):
        res = fut.result()
        jid = res["job_id"]
        results[jid] = res
        print(f"[m3.local] {{jid:24s}} | LSD2={{res.get('lsd2_tmrca')}} ({{res['lsd2_status']}}) | TT={{res.get('tt_tmrca')}} ({{res['treetime_status']}}) | {{res['wall_sec']}}s", flush=True)

(ROOT / "clock_results.json").write_text(json.dumps(results, indent=2))
print(f"ALL REMOTE CLOCK JOBS FINISHED in {{time.time() - t_start:.1f}}s")
""")

        tar_path = tmp_path / "clock_payload.tar.gz"
        with tarfile.open(tar_path, "w:gz") as tar:
            tar.add(jobs_dir, arcname="jobs")
            tar.add(tmp_path / "manifest.json", arcname="manifest.json")
            tar.add(remote_runner, arcname="run_remote_clock.py")

        subprocess.run(["ssh", "-o", "BatchMode=yes", REMOTE_HOST, f"rm -rf {REMOTE_DIR} && mkdir -p {REMOTE_DIR}"], check=True)
        subprocess.run(["scp", "-q", str(tar_path), f"{REMOTE_HOST}:{REMOTE_DIR}/clock_payload.tar.gz"], check=True)
        print(f"[+] Staged {len(job_specs)} alignments ({tar_path.stat().st_size / 1e6:.2f} MB) to {REMOTE_HOST}; executing IQ-TREE2 + LSD2 + TreeTime...", flush=True)

        run_cmd = f"cd {REMOTE_DIR} && tar -xzf clock_payload.tar.gz && python3 -u run_remote_clock.py"
        subprocess.run(["ssh", "-o", "BatchMode=yes", "-o", "ServerAliveInterval=60", REMOTE_HOST, run_cmd], check=True)

        local_json = tmp_path / "clock_results.json"
        subprocess.run(["scp", "-q", f"{REMOTE_HOST}:{REMOTE_DIR}/clock_results.json", str(local_json)], check=True)
        return json.loads(local_json.read_text())


def main():
    t0 = time.time()
    print("=" * 90)
    print("STEP 19: TREE-BASED CLOCK VALIDATION (ChronAeon vs IQ-TREE 2 LSD2 vs TreeTime)")
    print("=" * 90)

    p1_df = pd.read_csv(RESULTS_DIR / "phase1_chronaeon_clade_dating.csv")
    flt_meta = pd.read_csv(DATA_DIR / "lanl_flt_2022_metadata.csv")
    flt_seqs = load_fasta_dict(DATA_DIR / "lanl_flt_2022.fasta")
    h_to_col = get_hxb2_cols(flt_seqs, "K03455")

    coding_h_range = [h_to_col[h] for h in range(790, 9418) if h in h_to_col]
    print(f"[+] Loaded lanl_flt_2022 ({len(flt_seqs)} seqs, {len(coding_h_range)} HXB2 coding cols 790..9417)")

    seg_df = pd.read_csv(RESULTS_DIR / "hiv1_deconvolved_segments.csv")
    crf_longest = {}
    for sc in ["01_AE", "02_AG", "06_cpx", "07_BC", "08_BC", "12_BF"]:
        sub = seg_df[seg_df["short_code"] == sc].sort_values("span_nt", ascending=False)
        best = sub.iloc[0]
        crf_longest[sc] = {
            "trunk": str(best["trunk"]),
            "hxb2_start": int(best["hxb2_start"]),
            "hxb2_end": int(best["hxb2_end"]),
            "span_nt": int(best["span_nt"]),
            "short_gene": str(best["short_gene"]),
        }

    job_specs = []
    part_a_meta = []

    # --- Build Part A: 21 Clades (<= 150 tip stratified subsample) ---
    for _, p_row in p1_df.iterrows():
        cl = str(p_row["clade"])
        cl_df = flt_meta[(flt_meta["lanl_tag"] == cl) & flt_meta["year"].notna() & flt_meta["seq_id"].isin(flt_seqs)].copy()
        sub_df = temporal_stratified_subsample(cl_df, max_tips=150, seed=20261114)
        sids = list(sub_df["seq_id"])
        dates_map = {r["seq_id"]: float(r["year"]) for _, r in sub_df.iterrows()}
        country_map = {r["seq_id"]: (str(r["country"]) if pd.notna(r["country"]) else "UNK") for _, r in sub_df.iterrows()}

        # Filter alignment on full coding slice
        aln_sub = extract_filtered_alignment(sids, flt_seqs, coding_h_range, max_gap_frac=0.50)

        # Evaluate ChronAeon on the exact same subsample alignment & taxa
        ch_res = s2.date_clade_partition_chronaeon(aln_sub, dates_map, country_map=country_map, n_boot=500)

        jid = f"partA_clade_{cl}"
        job_specs.append({"job_id": jid, "seqs": aln_sub, "dates": dates_map})
        part_a_meta.append({
            "job_id": jid,
            "analysis_type": "clade_full_coding",
            "clade": cl,
            "segment_desc": "full_coding_790_9417",
            "n_full_clade": int(p_row["n_dated"]),
            "n_eval_tips": len(sids),
            "sample_min_yr": float(sub_df["year"].min()),
            "sample_max_yr": float(sub_df["year"].max()),
            "chronaeon_full_model": str(p_row["clock_model"]),
            "chronaeon_full_tmrca": float(p_row["t_mrca"]),
            "chronaeon_full_ci_lower": p_row["ci_lower"],
            "chronaeon_full_ci_upper": p_row["ci_upper"],
            "chronaeon_sub_model": ch_res["clock_model"] if ch_res else "NA",
            "chronaeon_sub_tmrca": ch_res["t_mrca"] if ch_res else np.nan,
            "chronaeon_sub_ci_lower": ch_res["ci_lower"] if ch_res else np.nan,
            "chronaeon_sub_ci_upper": ch_res["ci_upper"] if ch_res else np.nan,
            "chronaeon_sub_mu": ch_res["mu_subs_site_yr"] if ch_res else np.nan,
            "chronaeon_sub_r": ch_res["r"] if ch_res else np.nan,
        })

        # Sensitivity check on longest non-recombinant segment for 6 major CRFs
        if cl in crf_longest:
            seg_info = crf_longest[cl]
            h0, h1 = seg_info["hxb2_start"], seg_info["hxb2_end"]
            seg_cols = [h_to_col[h] for h in range(h0, h1 + 1) if h in h_to_col]
            seg_aln = extract_filtered_alignment(sids, flt_seqs, seg_cols, max_gap_frac=0.50)
            ch_seg = s2.date_clade_partition_chronaeon(seg_aln, dates_map, country_map=country_map, n_boot=500)
            jid_seg = f"partA_crfseg_{cl}"
            job_specs.append({"job_id": jid_seg, "seqs": seg_aln, "dates": dates_map})
            part_a_meta.append({
                "job_id": jid_seg,
                "analysis_type": "crf_longest_nonrecomb_segment",
                "clade": cl,
                "segment_desc": f"{seg_info['trunk']}:{h0}-{h1} ({seg_info['short_gene']}; {seg_info['span_nt']}nt)",
                "n_full_clade": int(p_row["n_dated"]),
                "n_eval_tips": len(sids),
                "sample_min_yr": float(sub_df["year"].min()),
                "sample_max_yr": float(sub_df["year"].max()),
                "chronaeon_full_model": str(p_row["clock_model"]),
                "chronaeon_full_tmrca": float(p_row["t_mrca"]),
                "chronaeon_full_ci_lower": p_row["ci_lower"],
                "chronaeon_full_ci_upper": p_row["ci_upper"],
                "chronaeon_sub_model": ch_seg["clock_model"] if ch_seg else "NA",
                "chronaeon_sub_tmrca": ch_seg["t_mrca"] if ch_seg else np.nan,
                "chronaeon_sub_ci_lower": ch_seg["ci_lower"] if ch_seg else np.nan,
                "chronaeon_sub_ci_upper": ch_seg["ci_upper"] if ch_seg else np.nan,
                "chronaeon_sub_mu": ch_seg["mu_subs_site_yr"] if ch_seg else np.nan,
                "chronaeon_sub_r": ch_seg["r"] if ch_seg else np.nan,
            })

    # --- Build Part B: All 15 Homologous Circulating Lineages (NOVEL_CRF_01..15) ---
    nov_df = pd.read_csv(RESULTS_DIR / "genbank_novel_circulating_recombinants.csv")
    gb31 = pd.read_csv(RESULTS_DIR / "genbank_31k_rhizaeon_classifications.csv")
    gb_yr = dict(zip(gb31["accession"].astype(str), pd.to_numeric(gb31["year"], errors="coerce")))
    flt_yr = dict(zip(flt_meta["accession"].astype(str), pd.to_numeric(flt_meta["year"], errors="coerce")))
    gb_ctry = dict(zip(gb31["accession"].astype(str), gb31["country"].fillna("UNK").astype(str)))
    flt_ctry = dict(zip(flt_meta["accession"].astype(str), flt_meta["country"].fillna("UNK").astype(str)))
    flt_acc_to_id = dict(zip(flt_meta["accession"].astype(str), flt_meta["seq_id"].astype(str)))

    hom_mask = nov_df["lineage_status"].astype(str).str.startswith("Homologous")
    hom_indices = nov_df.index[hom_mask].tolist()

    part_b_meta = []
    for lin_rank, row_idx in enumerate(hom_indices, 1):
        r = nov_df.loc[row_idx]
        lin_name = f"NOVEL_CRF_{lin_rank:02d}"
        formula = str(r["rhizaeon_formula"])
        barcode = str(r["representative_hxb2_barcode"])
        accs = [a.strip() for a in str(r["member_accessions"]).split(",") if a.strip()]

        # Find longest non-recombinant segment (> 1,000 nt)
        segs = []
        for s_idx, s_item in enumerate(barcode.split("|"), 1):
            tr, cs = s_item.split(":")
            h0, h1 = [int(x) for x in cs.split("-")]
            segs.append((s_idx, tr, h0, h1, h1 - h0 + 1))
        best_sidx, best_tr, best_h0, best_h1, best_len = max(segs, key=lambda x: x[4])

        # Gather dates & countries for all member accessions
        dated_accs = []
        all_years = []
        for a in accs:
            y = gb_yr.get(a, np.nan)
            if pd.isna(y):
                y = flt_yr.get(a, np.nan)
            if pd.notna(y):
                dated_accs.append((a, float(y)))
                all_years.append(float(y))
        n_dated = len(dated_accs)
        min_y = min(all_years) if all_years else np.nan
        max_y = max(all_years) if all_years else np.nan
        span_yr = (max_y - min_y) if len(all_years) >= 2 else 0.0

        # Extract longest-segment sequences for dated_accs when eligible for clock analysis
        jid_b = f"partB_{lin_name}"
        ch_nov = None
        eval_flag = (best_len > 1000 and n_dated >= 3 and span_yr >= 2.0)
        if eval_flag:
            seg_seqs = {}
            if lin_rank == 1:
                # Extract all 9 East African isolates cleanly from lanl_flt_2022.fasta along D:790-6335
                seg_cols = [h_to_col[h] for h in range(best_h0, best_h1 + 1) if h in h_to_col]
                raw_sub = {}
                for a, _ in dated_accs:
                    sid = flt_acc_to_id.get(a)
                    if sid and sid in flt_seqs:
                        raw_sub[a] = flt_seqs[sid]
                if len(raw_sub) == len(dated_accs):
                    seg_seqs = extract_filtered_alignment([a for a, _ in dated_accs], raw_sub, seg_cols, max_gap_frac=0.50)
            else:
                matches = glob.glob(str(TREES_DIR / f"lin{lin_rank:02d}_*_seg{best_sidx}_*.fasta"))
                if matches:
                    fa_recs = load_fasta_dict(Path(matches[0]))
                    raw_q = {}
                    for a, _ in dated_accs:
                        qk = f"Query_{a}"
                        if qk in fa_recs:
                            raw_q[a] = fa_recs[qk]
                    if len(raw_q) >= 3:
                        L_fa = len(next(iter(raw_q.values())))
                        seg_seqs = extract_filtered_alignment(list(raw_q.keys()), raw_q, list(range(L_fa)), max_gap_frac=0.50)

            if len(seg_seqs) >= 3:
                d_map = {a: y for a, y in dated_accs if a in seg_seqs}
                c_map = {a: gb_ctry.get(a, flt_ctry.get(a, "UNK")) for a in seg_seqs}
                # Fallback to isolate resampling since sample size <= 9
                ch_nov = s2.date_clade_partition_chronaeon(seg_seqs, d_map, country_map=None, min_span_years=2.0, n_boot=500)
                job_specs.append({"job_id": jid_b, "seqs": seg_seqs, "dates": d_map})
            else:
                eval_flag = False

        part_b_meta.append({
            "row_idx": row_idx,
            "lineage_id": lin_name,
            "rhizaeon_formula": formula,
            "n_patients": int(r["n_homologous_bp_patients"]),
            "n_dated": n_dated,
            "min_year": int(min_y) if pd.notna(min_y) else None,
            "max_year": int(max_y) if pd.notna(max_y) else None,
            "sample_span_yr": span_yr,
            "longest_seg_trunk": best_tr,
            "longest_seg_hxb2": f"{best_h0}-{best_h1}",
            "longest_seg_nt": best_len,
            "clock_evaluated": eval_flag,
            "job_id": jid_b if eval_flag else None,
            "chronaeon_tmrca": ch_nov["t_mrca"] if ch_nov else np.nan,
            "chronaeon_ci_lower": ch_nov["ci_lower"] if ch_nov else np.nan,
            "chronaeon_ci_upper": ch_nov["ci_upper"] if ch_nov else np.nan,
            "chronaeon_mu": ch_nov["mu_subs_site_yr"] if ch_nov else np.nan,
            "chronaeon_r": ch_nov["r"] if ch_nov else np.nan,
            "chronaeon_p": ch_nov["p_value"] if ch_nov else np.nan,
        })

    # Dispatch all jobs to m3.local
    remote_results = run_remote_clock_jobs(job_specs)

    # Assemble Part A DataFrame
    part_a_rows = []
    for meta in part_a_meta:
        rr = remote_results.get(meta["job_id"], {})
        ch_tm = meta["chronaeon_sub_tmrca"]
        lsd_tm = rr.get("lsd2_tmrca", np.nan)
        tt_tm = rr.get("tt_tmrca", np.nan)
        lsd_lo = rr.get("lsd2_ci_lo", np.nan)
        lsd_hi = rr.get("lsd2_ci_hi", np.nan)
        lsd_lo_note = "unbounded (<1800)" if (pd.notna(lsd_lo) and lsd_lo < 1800.0) else "bounded"
        part_a_rows.append({
            **meta,
            "lsd2_status": rr.get("lsd2_status", "NOT_RUN"),
            "lsd2_tmrca": round(lsd_tm, 2) if pd.notna(lsd_tm) else np.nan,
            "lsd2_ci_lower": round(lsd_lo, 2) if pd.notna(lsd_lo) else np.nan,
            "lsd2_ci_upper": round(lsd_hi, 2) if pd.notna(lsd_hi) else np.nan,
            "lsd2_ci_bound_note": lsd_lo_note if pd.notna(lsd_lo) else "NA",
            "lsd2_rate": rr.get("lsd2_rate", np.nan),
            "treetime_status": rr.get("treetime_status", "NOT_RUN"),
            "tt_tmrca": round(tt_tm, 2) if pd.notna(tt_tm) else np.nan,
            "tt_rate": rr.get("tt_rate", np.nan),
            "tt_r2": rr.get("tt_r2", np.nan),
            "abs_diff_ch_vs_lsd2": round(abs(ch_tm - lsd_tm), 2) if pd.notna(ch_tm) and pd.notna(lsd_tm) else np.nan,
            "abs_diff_ch_vs_tt": round(abs(ch_tm - tt_tm), 2) if pd.notna(ch_tm) and pd.notna(tt_tm) else np.nan,
            "abs_diff_lsd2_vs_tt": round(abs(lsd_tm - tt_tm), 2) if pd.notna(lsd_tm) and pd.notna(tt_tm) else np.nan,
        })

    df_a = pd.DataFrame(part_a_rows)
    out_a = RESULTS_DIR / "chronaeon_vs_lsd2_treetime_validation.csv"
    df_a.to_csv(out_a, index=False)

    print("\n--- PART A SUMMARY: 21 CLADES FULL CODING (IDENTICAL <=150-TIP SUBSAMPLE) ---")
    sub_21 = df_a[df_a["analysis_type"] == "clade_full_coding"].dropna(subset=["chronaeon_sub_tmrca", "lsd2_tmrca", "tt_tmrca"])
    for x_col, y_col, lbl in [
        ("chronaeon_sub_tmrca", "lsd2_tmrca", "ChronAeon_sub vs LSD2"),
        ("chronaeon_sub_tmrca", "tt_tmrca", "ChronAeon_sub vs TreeTime"),
        ("chronaeon_full_tmrca", "lsd2_tmrca", "ChronAeon_full vs LSD2"),
        ("chronaeon_full_tmrca", "tt_tmrca", "ChronAeon_full vs TreeTime"),
        ("lsd2_tmrca", "tt_tmrca", "LSD2 vs TreeTime"),
    ]:
        xv, yv = sub_21[x_col].values.astype(float), sub_21[y_col].values.astype(float)
        r_p, p_p = stats.pearsonr(xv, yv)
        r_s, p_s = stats.spearmanr(xv, yv)
        mad = float(np.median(np.abs(xv - yv)))
        print(f"  {lbl:28s} (n={len(sub_21)}): Pearson r={r_p:.3f} (p={p_p:.2e}) | Spearman rho={r_s:.3f} (p={p_s:.2e}) | MedAbsDiff={mad:.2f} yr")

    # Assemble Part B DataFrame and apply strict qualification rule
    part_b_rows = []
    for meta in part_b_meta:
        rr = remote_results.get(meta["job_id"], {}) if meta["job_id"] else {}
        ch_tm = meta["chronaeon_tmrca"]
        ch_r = meta["chronaeon_r"]
        lsd_tm = rr.get("lsd2_tmrca", np.nan)
        lsd_rate = rr.get("lsd2_rate", np.nan)
        tt_tm = rr.get("tt_tmrca", np.nan)
        tt_rate = rr.get("tt_rate", np.nan)
        tt_r2 = rr.get("tt_r2", np.nan)

        n_d = meta["n_dated"]
        sp = meta["sample_span_yr"]
        min_y, max_y = meta["min_year"], meta["max_year"]
        pos_slope = (
            pd.notna(ch_tm)
            and pd.notna(ch_r) and ch_r > 0
            and pd.notna(lsd_rate) and lsd_rate > 1e-5
            and pd.notna(lsd_tm) and lsd_tm <= min_y
        )
        if n_d == 0:
            origin_basis = "Undated"
            suffix = ""
            final_yr = "Undated"
        elif n_d >= 4 and sp >= 3.0 and pos_slope:
            origin_basis = "Clock_tMRCA"
            tm_int = int(round(ch_tm))
            suffix = f"@{tm_int}"
            final_yr = f"{ch_tm:.1f}"
        else:
            origin_basis = "Sampling_Span_Only"
            span_str = f"{min_y}-{max_y}" if min_y != max_y else f"{min_y}"
            suffix = f"[Span:{span_str}]"
            final_yr = "Undated (Insufficient Clock)"

        row_idx = meta["row_idx"]
        nov_df.loc[row_idx, "origin_basis"] = origin_basis
        nov_df.loc[row_idx, "formatted_suffix"] = suffix
        nov_df.loc[row_idx, "chronaeon_tmrca"] = final_yr

        part_b_rows.append({
            "lineage_id": meta["lineage_id"],
            "rhizaeon_formula": meta["rhizaeon_formula"],
            "n_patients": meta["n_patients"],
            "n_dated": n_d,
            "min_year": min_y,
            "max_year": max_y,
            "sample_span_yr": sp,
            "longest_seg_trunk": meta["longest_seg_trunk"],
            "longest_seg_hxb2": meta["longest_seg_hxb2"],
            "longest_seg_nt": meta["longest_seg_nt"],
            "origin_basis": origin_basis,
            "formatted_suffix": suffix,
            "chronaeon_tmrca": round(ch_tm, 2) if pd.notna(ch_tm) else np.nan,
            "chronaeon_ci_lower": round(meta["chronaeon_ci_lower"], 2) if pd.notna(meta["chronaeon_ci_lower"]) else np.nan,
            "chronaeon_ci_upper": round(meta["chronaeon_ci_upper"], 2) if pd.notna(meta["chronaeon_ci_upper"]) else np.nan,
            "chronaeon_r": round(ch_r, 3) if pd.notna(ch_r) else np.nan,
            "chronaeon_p": meta["chronaeon_p"],
            "lsd2_status": rr.get("lsd2_status", "NOT_ELIGIBLE"),
            "lsd2_tmrca": round(lsd_tm, 2) if pd.notna(lsd_tm) else np.nan,
            "lsd2_ci_lower": rr.get("lsd2_ci_lo", np.nan),
            "lsd2_ci_upper": rr.get("lsd2_ci_hi", np.nan),
            "lsd2_rate": lsd_rate,
            "treetime_status": rr.get("treetime_status", "NOT_ELIGIBLE"),
            "tt_tmrca": round(tt_tm, 2) if pd.notna(tt_tm) else np.nan,
            "tt_rate": tt_rate,
            "tt_r2": tt_r2,
        })

    # Populate origin_basis & formatted_suffix for remaining convergent rows in genbank_novel_circulating_recombinants.csv
    for idx in nov_df.index:
        if idx in hom_indices:
            continue
        sp_s = str(nov_df.loc[idx, "sampling_years"])
        if sp_s in ("Undated", "nan", ""):
            nov_df.loc[idx, "origin_basis"] = "Undated"
            nov_df.loc[idx, "formatted_suffix"] = ""
        else:
            nov_df.loc[idx, "origin_basis"] = "Sampling_Span_Only"
            pts = sp_s.split("-")
            s_fmt = pts[0] if (len(pts) == 2 and pts[0] == pts[1]) else sp_s
            nov_df.loc[idx, "formatted_suffix"] = f"[Span:{s_fmt}]"
            nov_df.loc[idx, "chronaeon_tmrca"] = "Undated (Convergent)"

    df_b = pd.DataFrame(part_b_rows)
    out_b = RESULTS_DIR / "novel_crfs_15_clock_tmrca.csv"
    df_b.to_csv(out_b, index=False)
    nov_df.to_csv(RESULTS_DIR / "genbank_novel_circulating_recombinants.csv", index=False)

    print(f"\n[+] Saved {out_a.name} ({len(df_a)} rows) and {out_b.name} ({len(df_b)} rows) in {time.time() - t0:.1f}s")
    print("\nPART B: NOVEL_CRF_01..15 ORIGIN QUALIFICATION:")
    print(df_b[["lineage_id", "n_patients", "n_dated", "sample_span_yr", "longest_seg_hxb2", "origin_basis", "formatted_suffix", "chronaeon_tmrca", "lsd2_tmrca", "tt_tmrca"]].to_string(index=False))


if __name__ == "__main__":
    main()
