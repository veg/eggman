#!/usr/bin/env python3
"""
26_unified_chronaeon_temporal_basis_sc2_iav.py
==============================================
Unifies RhizAeon + ChronAeon + EGGMAN across SARS-CoV-2, Influenza A (H5Nx),
and Sarbecoviruses:
  1. Constructs ChronAeon exponential time-decay consensus bases (t_{1/2}=0.25 yr
     for SC2 across all 16 spanning Tier-1 trunks + 16 Tier-2 recombinants;
     t_{1/2}=0.50 yr for H5Nx across all 10 8-segment donor trunks + certified
     reassortants), sweeping half-lives vs unweighted full-clade majority & the
     earlier first<=25 heuristic, quantifying coding SNP diffs & stripped
     late-wave homoplasies.
  2. Runs ChronAeon Steiner + Gamma4 (alpha=0.20) OLS / RCS tree-free molecular
     clock dating + block bootstrap alongside orthogonal IQ-TREE 2 (GTR+G4) +
     LSD2 + TreeTime across SC2 waves/recombinants and H5Nx genotypes +
     shared-segment chronological symmetry breaking (A1 vs B3.13; A3 vs D1.1).
  3. Upgrades results/sc2_pack_v3_chronaeon/, results/iav_pack_v2_chronaeon/,
     rhiztype/packs/{sars_cov2_wuhan,influenza_a_8seg,sarbecovirus_wuhan}/, and
     rhiztype/docs/examples/sc2_spanning_16.pack.json.
  4. Re-verifies the 3-criterion Spanning Basis Rule on 10,062 SC2 genomes and
     runs standard vs --temporal-guard ablations (collection_date vs
     submission/deposition date) across:
       - SC2: 1,000 curated benchmark, 10,062 World Union, 23,404 BV-BRC,
              and 100,000-genome BV-BRC stream.
       - IAV: 1,200 SRA cohort, 650 OOD holdout, and 3,980 BV-BRC 2020-2026
              8-segment H5Nx isolates.
"""

from __future__ import annotations

import argparse
import collections
import concurrent.futures
import csv
import datetime as dt_mod
import glob
import gzip
import importlib
import json
import lzma
import math
import os
import re
import shutil
import statistics
import subprocess
import sys
import tempfile
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Sequence, Tuple

import numpy as np
from scipy import stats

ROOT = Path("/Users/sergei/Projects/HIV-recombination")
SCRIPTS_DIR = ROOT / "scripts"
RESULTS_DIR = ROOT / "results"
PACKS_DIR = ROOT / "rhiztype" / "packs"
EXAMPLES_DIR = ROOT / "rhiztype" / "docs" / "examples"
EGGMAN_BIN = ROOT / "rhiztype" / "target" / "release" / "eggman"

sys.path.insert(0, str(SCRIPTS_DIR))
sys.path.insert(0, "/Users/sergei/Projects/TOGA_MEME/axomeme_repo/chronaeon/src")
sys.path.insert(0, "/Users/sergei/Projects/TOGA_MEME/axomeme_repo/aeon-core/src")

from chronaeon.dating_io import (  # noqa: E402
    generate_consensus_sequence,
    generate_time_decay_consensus_sequence,
    parse_date_to_decimal,
)
s2 = importlib.import_module("02_phase1_pure_subtype_and_temporal_audit")

IQTREE_BIN = "/Users/sergei/.local/bin/iqtree2"
TREETIME_BIN = "/Users/sergei/miniconda3/bin/treetime"
NEXTCLADE_BIN = "/Users/sergei/miniconda3/bin/nextclade"

WORLD_FASTA_XZ = Path("/Users/sergei/Projects/trees/data/alignments/viral_and_bwz_codec/sars_cov_2_10062_world_union.fasta.xz")
WORLD_META_GLOB = "/Users/sergei/Projects/trees/data/alignments/viral_and_bwz_codec/sc2_*.metadata.tsv.xz"
BVBRC_23K_DIR = ROOT / "data" / "bvbrc_sc2_stream"
SC2_HITS_FASTA_GZ = ROOT / "data" / "bvbrc_stream_hits" / "sc2_novel_hits.fasta.gz"
IAV_HITS_FASTA_GZ = ROOT / "data" / "bvbrc_stream_hits" / "iav_h5_novel_hits.fasta.gz"

SRA_FASTA_DIR = Path("/Users/sergei/Projects/stuff/bovine/data/andersen_avian_influenza/fasta")
SRA_GENOFLU_TSV = Path("/Users/sergei/Projects/stuff/bovine/data/andersen_metadata/genoflu_results.tsv")
SRA_TABLE_TSV = Path("/Users/sergei/Projects/stuff/bovine/data/andersen_metadata/SraRunTable_automated_normalized.tsv")

CACHE_DIR = Path("/Users/sergei/.gemini/antigravity-cli/brain/d7d3d768-93f4-475c-81ac-ccb936c7805d/scratch/chronaeon_cache")
CACHE_DIR.mkdir(parents=True, exist_ok=True)

SC2_V3_DIR = RESULTS_DIR / "sc2_pack_v3_chronaeon"
IAV_V2_DIR = RESULTS_DIR / "iav_pack_v2_chronaeon"

OUT_SUMMARY_JSON = RESULTS_DIR / "chronaeon_multipathogen_basis_and_dating.json"
OUT_SC2_CLOCK_TSV = RESULTS_DIR / "chronaeon_sc2_dating_vs_lsd2_treetime.tsv"
OUT_IAV_CLOCK_TSV = RESULTS_DIR / "chronaeon_iav_dating_vs_lsd2_treetime.tsv"
OUT_ABLATION_TSV = RESULTS_DIR / "chronaeon_temporal_guard_ablation.tsv"

BVBRC_GENOME_URL = "https://www.bv-brc.org/api/genome/"
BVBRC_SEQ_URL = "https://www.bv-brc.org/api/genome_sequence/"
HEADERS_JSON = {"Accept": "application/json", "User-Agent": "curl/8.7.1"}
HEADERS_FASTA_POST = {
    "Content-Type": "application/rqlquery+x-www-form-urlencoded",
    "Accept": "application/dna+fasta",
    "User-Agent": "curl/8.7.1",
}

SC2_LEN = 29903
SC2_COD_START = 266
SC2_COD_END = 29674
SC2_TIER1_16 = [
    "B", "B.1.1.7", "P.1", "B.1.617.2", "BA.1", "BA.2", "BA.4", "BA.5",
    "BJ.1", "BM.1.1.1", "JN.1", "KP.3.3", "KS.1.1", "BA.3.2", "LF.7", "LP.8.1",
]
SC2_TIER2_LABELS = [
    "XBB", "XBB.1.5", "XBB.1.16", "XD", "XE", "XF", "XBC", "XBF",
    "XAY", "XDK", "XEC", "XBK", "XDV", "XFZ", "XFC", "XFN",
]
SC2_BASE_11_LABELS = [
    "B", "B.1.1.7", "B.1.617.2", "BA.1", "BA.2", "BA.5",
    "BJ.1", "BM.1.1.1", "JN.1", "KP.3.3", "KS.1.1",
]
SC2_ORFS = [
    ("ORF1a", 266, 13468),
    ("ORF1b", 13468, 21555),
    ("S", 21563, 25384),
    ("ORF3a", 25393, 26220),
    ("E", 26245, 26472),
    ("M", 26523, 27191),
    ("ORF6", 27202, 27387),
    ("ORF7a", 27394, 27759),
    ("ORF7b", 27756, 27887),
    ("ORF8", 27894, 28259),
    ("N", 28274, 29533),
    ("ORF10", 29558, 29674),
]
CODON_TABLE = {
    "TTT":"F","TTC":"F","TTA":"L","TTG":"L","CTT":"L","CTC":"L","CTA":"L","CTG":"L",
    "ATT":"I","ATC":"I","ATA":"I","ATG":"M","GTT":"V","GTC":"V","GTA":"V","GTG":"V",
    "TCT":"S","TCC":"S","TCA":"S","TCG":"S","CCT":"P","CCC":"P","CCA":"P","CCG":"P",
    "ACT":"T","ACC":"T","ACA":"T","ACG":"T","GCT":"A","GCC":"A","GCA":"A","GCG":"A",
    "TAT":"Y","TAC":"Y","TAA":"*","TAG":"*","CAT":"H","CAC":"H","CAA":"Q","CAG":"Q",
    "AAT":"N","AAC":"N","AAA":"K","AAG":"K","GAT":"D","GAC":"D","GAA":"E","GAG":"E",
    "TGT":"C","TGC":"C","TGA":"*","TGG":"W","CGT":"R","CGC":"R","CGA":"R","CGG":"R",
    "AGT":"S","AGC":"S","AGA":"R","AGG":"R","GGT":"G","GGC":"G","GGA":"G","GGG":"G",
}

IAV_SEGS = [
    ("PB2", 1, 2280),
    ("PB1", 2281, 4554),
    ("PA", 4555, 6705),
    ("HA", 6706, 8409),
    ("NP", 8410, 9906),
    ("NA", 9907, 11316),
    ("MP", 11317, 12298),
    ("NS", 12299, 13136),
]
IAV_SEG_NUM = {"1":"PB2","2":"PB1","3":"PA","4":"HA","5":"NP","6":"NA","7":"MP","8":"NS"}
IAV_LEN = 13136
IAV_TRUNKS = [
    "EA1", "EA3", "EA4", "EA2", "NAm_B313", "NAm_D1",
    "NAm_B36", "NAm_C2", "NAm_B2", "NAm_B32",
]
IAV_DONOR_SEGMENTS = {
    "EA1": ["PB2","PB1","PA","HA","NP","NA","MP","NS"],
    "EA3": ["PB2","PB1","PA","HA","NP","NA","MP","NS"],
    "EA4": ["PB2","PB1","PA","HA","NP","NA","MP","NS"],
    "EA2": ["PB1","HA","NA"],
    "NAm_B313": ["PB2","PB1","NP","NS"],
    "NAm_D1": ["PB2","PA","NP","NA"],
    "NAm_B36": ["PB2","NP"],
    "NAm_C2": ["PB2","PB1","PA","NP"],
    "NAm_B2": ["PB2","PB1","NP"],
    "NAm_B32": ["PB2"],
}


def log(msg: str) -> None:
    ts = dt_mod.datetime.now().strftime("%H:%M:%S")
    print(f"[{ts}] {msg}", flush=True)


def parse_clean_date(raw: Any) -> Tuple[Optional[float], str]:
    """
    Parse date strictly requiring at least YYYY-MM for high-resolution dating/gating,
    or fall back to YYYY.5 for bare-year labels.
    Returns (decimal_year, precision_flag in {'ymd','ym','y','none'}).
    """
    if raw is None:
        return None, "none"
    s = str(raw).strip()
    if not s or s.lower() in ("na", "nan", "none", "null", "missing", "unknown"):
        return None, "none"
    s10 = s[:10]
    m_ymd = re.match(r"^(\d{4})-(\d{2})-(\d{2})$", s10)
    if m_ymd and "XX" not in s10.upper():
        try:
            y, m, d = int(m_ymd.group(1)), int(m_ymd.group(2)), int(m_ymd.group(3))
            if 1900 <= y <= 2100 and 1 <= m <= 12 and 1 <= d <= 31:
                val = float(parse_date_to_decimal(s10))
                if math.isfinite(val):
                    return val, "ymd"
        except Exception:
            pass
    m_ym = re.match(r"^(\d{4})-(\d{2})$", s[:7])
    if m_ym and "XX" not in s[:7].upper():
        try:
            y, m = int(m_ym.group(1)), int(m_ym.group(2))
            if 1900 <= y <= 2100 and 1 <= m <= 12:
                val = float(parse_date_to_decimal(s[:7]))
                if math.isfinite(val):
                    return val, "ym"
        except Exception:
            pass
    m_y = re.match(r"^(\d{4})(?:-XX(?:-XX)?)?$", s, re.IGNORECASE)
    if m_y:
        y = int(m_y.group(1))
        if 1900 <= y <= 2100:
            return float(y) + 0.5, "y"
    try:
        val = float(parse_date_to_decimal(s))
        if math.isfinite(val) and 1900.0 <= val <= 2100.0:
            return val, "ymd" if re.search(r"\d{4}-\d{2}-\d{2}", s) else "ym"
    except Exception:
        pass
    return None, "none"


def read_fasta(path: Path) -> Dict[str, str]:
    opener = lzma.open if str(path).endswith(".xz") else (gzip.open if str(path).endswith(".gz") else open)
    out = {}
    hdr = None
    chunks: List[str] = []
    with opener(path, "rt", errors="replace") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if hdr is not None:
                    out[hdr] = "".join(chunks).upper()
                hdr = line[1:].strip()
                chunks = []
            else:
                chunks.append(line)
        if hdr is not None:
            out[hdr] = "".join(chunks).upper()
    return out


def http_get_json(url: str, timeout: int = 45, retries: int = 5) -> Tuple[str, Any]:
    last_err = None
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers=HEADERS_JSON)
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                cr = resp.headers.get("Content-Range", "")
                body = resp.read().decode("utf-8", errors="replace")
                return cr, json.loads(body)
        except urllib.error.HTTPError as e:
            last_err = e
            wait = 1.5 * (2 ** attempt)
            if e.code in (429, 503):
                wait *= 2.0
            time.sleep(wait)
        except Exception as e:
            last_err = e
            time.sleep(1.2 * (2 ** attempt))
    raise RuntimeError(f"BV-BRC GET failed for {url}: {last_err}")


def bvbrc_fetch_pages(base_rql: str, select_fields: str, max_rows: Optional[int] = None, page_size: int = 10000) -> List[Dict[str, Any]]:
    collected: List[Dict[str, Any]] = []
    seen: set = set()
    last_gid = None
    while True:
        if max_rows is not None and len(collected) >= max_rows:
            break
        q = f"{base_rql}&sort(+genome_id)&select({select_fields})&limit({page_size})"
        if last_gid is not None:
            q = f"{base_rql}&gt(genome_id,{urllib.parse.quote(last_gid)})&sort(+genome_id)&select({select_fields})&limit({page_size})"
        _, data = http_get_json(f"{BVBRC_GENOME_URL}?{q}")
        if not data:
            break
        added = 0
        for row in data:
            gid = row.get("genome_id")
            if gid and gid not in seen:
                seen.add(gid)
                collected.append(row)
                added += 1
                if max_rows is not None and len(collected) >= max_rows:
                    break
        new_last = data[-1].get("genome_id")
        if len(data) < page_size or new_last == last_gid or added == 0:
            break
        last_gid = new_last
    return collected


def post_fasta_batch(genome_ids: List[str], retries: int = 5) -> str:
    body = f"in(genome_id,({','.join(genome_ids)}))&limit(2000)".encode("utf-8")
    last_err = None
    for attempt in range(retries):
        try:
            req = urllib.request.Request(BVBRC_SEQ_URL, data=body, headers=HEADERS_FASTA_POST, method="POST")
            with urllib.request.urlopen(req, timeout=60) as resp:
                return resp.read().decode("utf-8", errors="replace")
        except Exception as e:
            last_err = e
            time.sleep(1.5 * (2 ** attempt))
    raise RuntimeError(f"BV-BRC FASTA POST failed: {last_err}")


def parse_bvbrc_fasta_chunks(texts: Sequence[str], valid_gids: set, acc_to_gid: Dict[str, str]) -> Dict[str, str]:
    GID_RE = re.compile(r"\|\s*([0-9]+\.[0-9]+)\s*\]\s*$")
    ACCN_RE = re.compile(r"^>accn\|([^\s|]+)")
    seqs_parts: Dict[str, List[str]] = collections.defaultdict(list)
    for txt in texts:
        cur_gid = None
        cur_chunks: List[str] = []
        for line in txt.splitlines():
            if line.startswith(">"):
                if cur_gid and cur_chunks:
                    seqs_parts[cur_gid].append("".join(cur_chunks).upper())
                cur_chunks = []
                m = GID_RE.search(line)
                gid = m.group(1) if m else None
                if not gid or gid not in valid_gids:
                    am = ACCN_RE.match(line)
                    if am:
                        gid = acc_to_gid.get(am.group(1))
                cur_gid = gid if (gid in valid_gids) else None
            elif cur_gid is not None:
                s = line.strip()
                if s:
                    cur_chunks.append(s)
        if cur_gid and cur_chunks:
            seqs_parts[cur_gid].append("".join(cur_chunks).upper())
    return {gid: "".join(parts) for gid, parts in seqs_parts.items()}


def fetch_bvbrc_seqs_by_gids(rows: List[Dict[str, Any]], batch_size: int = 200, workers: int = 8) -> Dict[str, str]:
    valid_gids = {r["genome_id"] for r in rows if r.get("genome_id")}
    acc_to_gid = {}
    for r in rows:
        gid = r.get("genome_id")
        acc = str(r.get("genbank_accessions") or "").strip().split(",")[0].strip()
        if acc and gid:
            acc_to_gid[acc] = gid
    gids = sorted(valid_gids)
    batches = [gids[i:i + batch_size] for i in range(0, len(gids), batch_size)]
    texts: List[str] = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as ex:
        futs = [ex.submit(post_fasta_batch, b) for b in batches]
        for fut in concurrent.futures.as_completed(futs):
            texts.append(fut.result())
    return parse_bvbrc_fasta_chunks(texts, valid_gids, acc_to_gid)


def fetch_bvbrc_fields_for_gids(gids: List[str], select_fields: str, batch_size: int = 400, workers: int = 8) -> Dict[str, Dict[str, Any]]:
    batches = [gids[i:i + batch_size] for i in range(0, len(gids), batch_size)]
    out: Dict[str, Dict[str, Any]] = {}

    def _one(b: List[str]) -> List[Dict[str, Any]]:
        url = f"{BVBRC_GENOME_URL}?in(genome_id,({','.join(b)}))&select({select_fields})&limit({len(b) + 20})"
        _, data = http_get_json(url, timeout=45, retries=5)
        return data

    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as ex:
        for res in ex.map(_one, batches):
            for r in res:
                if r.get("genome_id"):
                    out[r["genome_id"]] = r
    return out


def align_sc2_to_wuhan_nextclade(seqs_dict: Dict[str, str], wuhan_ref_seq: str) -> Dict[str, str]:
    """Align raw SC2 genomes onto exact 29,903-nt MN908947 reference frame via Nextclade."""
    if not seqs_dict:
        return {}
    aligned: Dict[str, str] = {}
    need: Dict[str, str] = {}
    for k, s in seqs_dict.items():
        if len(s) == SC2_LEN and s[:15].count("-") < 15:
            aligned[k] = s
        else:
            need[k] = s
    if not need:
        return aligned
    with tempfile.TemporaryDirectory() as td:
        td_p = Path(td)
        ref_fa = td_p / "ref.fasta"
        q_fa = td_p / "q.fasta"
        out_fa = td_p / "aligned.fasta"
        ref_fa.write_text(f">Wuhan\n{wuhan_ref_seq}\n")
        with open(q_fa, "w") as f:
            for k, s in need.items():
                f.write(f">{k}\n{s}\n")
        subprocess.run(
            [NEXTCLADE_BIN, "run", "--input-ref", str(ref_fa), "--output-fasta", str(out_fa), "--silent", str(q_fa)],
            check=True,
        )
        got = read_fasta(out_fa)
        for k, s in got.items():
            if len(s) == SC2_LEN:
                aligned[k] = s
    return aligned


def seq_dict_to_matrix(seq_dict: Dict[str, str], taxa: List[str], start0: int = 0, end0: Optional[int] = None) -> np.ndarray:
    cmap = {ord("A"): 0, ord("C"): 1, ord("G"): 2, ord("T"): 3, ord("U"): 3}
    L0 = len(seq_dict[taxa[0]])
    e0 = L0 if end0 is None else end0
    L = e0 - start0
    M = np.full((len(taxa), L), -1, dtype=np.int8)
    for i, t in enumerate(taxa):
        sub = seq_dict[t][start0:e0].encode("ascii", "ignore")
        arr = np.frombuffer(sub, dtype=np.uint8)
        row = np.full(len(arr), -1, dtype=np.int8)
        for ch, val in cmap.items():
            row[arr == ch] = val
        M[i, :len(row)] = row
    return M


def chronaeon_clock_eval(
    seq_dict: Dict[str, str],
    dates_map: Dict[str, float],
    meta_map: Dict[str, Dict[str, Any]],
    taxa: List[str],
    half_life: float,
    slice_ranges: Optional[List[Tuple[int, int]]] = None,
    n_boot: int = 1000,
    seed: int = 42,
) -> Dict[str, Any]:
    """
    Evaluate ChronAeon Steiner + Gamma4 (alpha=0.20) tree-free clock dating
    with country x year-quarter (or year-month when <4 blocks) block bootstrap.
    """
    valid = [t for t in taxa if t in dates_map and math.isfinite(dates_map[t])]
    if len(valid) < 6:
        return {"status": "TOO_FEW", "n": len(valid)}

    if slice_ranges is None:
        seqs_eval = {t: seq_dict[t] for t in valid}
    else:
        seqs_eval = {
            t: "".join(seq_dict[t][a:b] for (a, b) in slice_ranges)
            for t in valid
        }

    M = seq_dict_to_matrix(seqs_eval, valid)
    L_val = np.maximum(np.sum(M >= 0, axis=1).astype(np.float64), 1.0)
    ev = s2.FastProfileEvaluator(M, L_val)

    times = np.array([dates_map[t] for t in valid], dtype=np.float64)
    t_min = float(np.min(times))
    t_max = float(np.max(times))
    eff_gamma = float(math.log(2.0) / half_life)
    w0 = np.exp(-eff_gamma * (times - t_min))
    w0 /= np.sum(w0)
    w_sq0 = np.array([float(np.sum(w0 * w0))], dtype=np.float64)
    dists0 = ev.eval_batch_divs(w0[None, :], w_sq0, alpha_g=0.20)[:, 0]

    # Triage clock-outlier artifacts via MAD of OLS residuals before final clock fit
    t_mrca_raw, mu_raw, ci_lo_raw, ci_hi_raw, g_raw, r_raw, p_raw, mean_d_raw = s2._fast_ols_tmrca(times, dists0)
    kept_mask = np.ones(len(valid), dtype=bool)
    if np.isfinite(mu_raw) and len(valid) >= 12:
        res0 = dists0 - (mu_raw * (times - np.mean(times)) + np.mean(dists0))
        med = float(np.median(res0))
        mad = float(np.median(np.abs(res0 - med)))
        scale = max(1e-5, 1.4826 * mad)
        if np.sum(np.abs(res0 - med) <= 4.5 * scale) >= max(6, int(0.8 * len(valid))):
            kept_mask = np.abs(res0 - med) <= 4.5 * scale

    if int(np.sum(kept_mask)) < len(valid):
        valid = [t for t, k in zip(valid, kept_mask) if k]
        M = M[kept_mask]
        L_val = L_val[kept_mask]
        ev = s2.FastProfileEvaluator(M, L_val)
        times = times[kept_mask]
        t_min = float(np.min(times))
        t_max = float(np.max(times))
        w0 = np.exp(-eff_gamma * (times - t_min))
        w0 /= np.sum(w0)
        w_sq0 = np.array([float(np.sum(w0 * w0))], dtype=np.float64)
        dists0 = ev.eval_batch_divs(w0[None, :], w_sq0, alpha_g=0.20)[:, 0]

    t_ols, mu_ols, f_lo, f_hi, fieller_g, r_val, p_val, mean_d = s2._fast_ols_tmrca(times, dists0)
    t_rcs, mu_rcs, rate_ratio, p_rcs, use_rcs = s2._fast_rcs_tmrca(times, dists0)

    # Check whether RCS triggers cleanly
    model_used = "RCS" if use_rcs else "OLS"
    point_est = float(t_rcs if use_rcs else t_ols)

    # Construct country x year-quarter (or year-month when <12 quarterly blocks or narrow span) blocks
    span = t_max - t_min
    use_monthly = span < 0.8
    for _pass in range(2):
        blocks_map: Dict[Tuple[str, str], List[int]] = collections.defaultdict(list)
        for i, t in enumerate(valid):
            ctry = str(meta_map.get(t, {}).get("country") or "UNK")
            y = dates_map[t]
            if use_monthly:
                bucket = f"{int(y)}M{min(12, int((y - int(y)) * 12) + 1):02d}"
            else:
                bucket = f"{int(y)}Q{min(4, int((y - int(y)) * 4) + 1)}"
            blocks_map[(ctry, bucket)].append(i)
        if len(blocks_map) >= 12 or use_monthly:
            break
        use_monthly = True
    if len(blocks_map) < 5:
        blocks_map = {(t, "single"): [i] for i, t in enumerate(valid)}

    block_lists = [np.array(idxs, dtype=np.int32) for idxs in blocks_map.values()]
    rng = np.random.default_rng(seed)
    n_b = len(block_lists)
    obs_mu = float(mu_rcs if use_rcs else mu_ols) if math.isfinite(mu_rcs if use_rcs else mu_ols) else 0.0

    # Resample blocks in vectorized batches of 250
    boot_tmrcas = []
    batch_sz = 250
    for b0 in range(0, n_boot, batch_sz):
        bs = min(batch_sz, n_boot - b0)
        W_batch = np.zeros((bs, len(valid)), dtype=np.float64)
        wsq_batch = np.zeros(bs, dtype=np.float64)
        chosen_sets = []
        for r_i in range(bs):
            picks = rng.integers(0, n_b, size=n_b)
            samp = np.concatenate([block_lists[p] for p in picks])
            chosen_sets.append(samp)
            t_s = times[samp]
            tm_s = float(np.min(t_s))
            counts = np.bincount(samp, minlength=len(valid)).astype(np.float64)
            raw_w = counts * np.exp(-eff_gamma * (times - tm_s))
            ws = float(np.sum(raw_w))
            if ws > 0:
                raw_w /= ws
            W_batch[r_i] = raw_w
            wsq_batch[r_i] = float(np.sum(raw_w * raw_w))
        D_batch = ev.eval_batch_divs(W_batch, wsq_batch, alpha_g=0.20)
        for r_i in range(bs):
            samp = chosen_sets[r_i]
            ts_b = times[samp]
            if span > 0 and (float(np.max(ts_b)) - float(np.min(ts_b))) < max(0.03, 0.25 * span):
                continue
            ds_b = D_batch[samp, r_i]
            if model_used == "RCS":
                val_b, mu_b, _, _, ok_b = s2._fast_rcs_tmrca(ts_b, ds_b)
                if not ok_b or not math.isfinite(val_b):
                    val_b, mu_b, *_ = s2._fast_ols_tmrca(ts_b, ds_b)
            else:
                val_b, mu_b, *_ = s2._fast_ols_tmrca(ts_b, ds_b)
            if obs_mu > 0 and (not math.isfinite(mu_b) or mu_b <= 0.15 * obs_mu):
                continue
            if val_b is not None and math.isfinite(val_b) and 1900.0 <= val_b <= t_max:
                boot_tmrcas.append(float(val_b))

    ci_lo = float(np.percentile(boot_tmrcas, 2.5)) if len(boot_tmrcas) >= 20 else f_lo
    ci_hi = float(np.percentile(boot_tmrcas, 97.5)) if len(boot_tmrcas) >= 20 else f_hi
    if (ci_lo is None or not math.isfinite(ci_lo)) and math.isfinite(point_est):
        ci_lo = f_lo
    if (ci_hi is None or not math.isfinite(ci_hi)) and math.isfinite(point_est):
        ci_hi = f_hi

    return {
        "status": "OK" if math.isfinite(point_est) else "NONPOS_SLOPE",
        "n": len(valid),
        "n_blocks": n_b,
        "t_min": round(t_min, 4),
        "t_max": round(t_max, 4),
        "model": model_used,
        "t_mrca": round(point_est, 4) if math.isfinite(point_est) else None,
        "ci_lower": round(ci_lo, 4) if (ci_lo is not None and math.isfinite(ci_lo)) else None,
        "ci_upper": round(ci_hi, 4) if (ci_hi is not None and math.isfinite(ci_hi)) else None,
        "ols_t_mrca": round(t_ols, 4) if math.isfinite(t_ols) else None,
        "fieller_ci_lower": round(f_lo, 4) if math.isfinite(f_lo) else None,
        "fieller_ci_upper": round(f_hi, 4) if math.isfinite(f_hi) else None,
        "mu": float(mu_rcs if use_rcs else mu_ols) if math.isfinite(mu_rcs if use_rcs else mu_ols) else None,
        "r": round(r_val, 4) if math.isfinite(r_val) else None,
        "r2": round(r_val ** 2, 4) if math.isfinite(r_val) else None,
        "p_value": float(p_val) if math.isfinite(p_val) else None,
        "mean_root_div": round(mean_d, 6) if math.isfinite(mean_d) else None,
        "n_boot_valid": len(boot_tmrcas),
        "kept_taxa": valid,
        "dists": dists0,
    }


def stratified_clock_subsample(
    taxa: List[str],
    seq_dict: Dict[str, str],
    dates_map: Dict[str, float],
    max_tips: int = 80,
    slice_ranges: Optional[List[Tuple[int, int]]] = None,
) -> List[str]:
    """Select up to max_tips non-identical sequences evenly across sampling time."""
    if slice_ranges is None:
        sub_seqs = {t: seq_dict[t] for t in taxa}
    else:
        sub_seqs = {t: "".join(seq_dict[t][a:b] for (a, b) in slice_ranges) for t in taxa}
    # Prefer dated taxa and deduplicate exact identical strings keeping earliest
    ordered = sorted([t for t in taxa if t in dates_map], key=lambda t: (dates_map[t], t))
    uniq: List[str] = []
    seen_seqs: set = set()
    for t in ordered:
        s = sub_seqs[t]
        if s not in seen_seqs:
            seen_seqs.add(s)
            uniq.append(t)
    if len(uniq) <= max_tips:
        return uniq
    # Stratified uniform index pick keeping endpoints (min and max date)
    idxs = sorted(set(np.round(np.linspace(0, len(uniq) - 1, max_tips)).astype(int).tolist()))
    return [uniq[i] for i in idxs]


def run_iqtree_lsd2_treetime_job(
    job_id: str,
    seqs_sub: Dict[str, str],
    dates_sub: Dict[str, float],
    min_plausible_year: float = 1990.0,
) -> Dict[str, Any]:
    """Run IQ-TREE 2 (GTR+G4 --fast) + LSD2 + TreeTime on a clean sanitized alignment."""
    out: Dict[str, Any] = {"job_id": job_id, "n_eval_tips": len(seqs_sub), "lsd2_status": "NOT_RUN", "treetime_status": "NOT_RUN"}
    if len(seqs_sub) < 6 or (max(dates_sub.values()) - min(dates_sub.values())) <= 1e-4:
        out["lsd2_status"] = "TOO_NARROW"
        out["treetime_status"] = "TOO_NARROW"
        return out

    with tempfile.TemporaryDirectory() as td:
        td_p = Path(td)
        fa = td_p / "aln.fasta"
        dt_lsd = td_p / "dates.tsv"
        dt_tt = td_p / "dates_tt.csv"
        pref = td_p / "iq"
        tt_dir = td_p / "tt_out"

        id_map = {orig: f"T{i+1:04d}" for i, orig in enumerate(sorted(seqs_sub))}
        with open(fa, "w") as ff, open(dt_lsd, "w") as lf, open(dt_tt, "w") as tf:
            tf.write("name,date\n")
            for orig, tid in id_map.items():
                yr = float(dates_sub[orig])
                ff.write(f">{tid}\n{seqs_sub[orig]}\n")
                lf.write(f"{tid}\t{yr:.5f}\n")
                tf.write(f"{tid},{yr:.5f}\n")

        cmd_iq = [
            IQTREE_BIN, "-s", str(fa), "-m", "GTR+G4", "--fast",
            "--date", str(dt_lsd), "--date-ci", "100", "--date-options", "-l 0",
            "-T", "2", "-redo", "--prefix", str(pref), "--quiet",
        ]
        res_iq = subprocess.run(cmd_iq, capture_output=True, text=True)
        lsd_file = Path(str(pref) + ".timetree.lsd")
        treefile = Path(str(pref) + ".treefile")
        if res_iq.returncode == 0 and lsd_file.exists():
            txt = lsd_file.read_text()
            m_ci = re.search(
                r"rate\s+([0-9.eE+-]+)\s*\[([0-9.eE+-]+);\s*([0-9.eE+-]+)\],\s*tMRCA\s+([0-9.eE+-]+)\s*\[([^;]+);\s*([^\]]+)\]",
                txt,
            )
            m_pt = re.search(r"-\s*Dating results:\s*\n\s*rate\s+([0-9.eE+-]+),\s*tMRCA\s+([0-9.eE+-]+)", txt)
            if m_ci:
                r_val = float(m_ci.group(1))
                t_val = round(float(m_ci.group(4)), 4)
                if not math.isfinite(t_val) or t_val < min_plausible_year or r_val <= 1e-7:
                    out["lsd2_status"] = "DEGENERATE_RATE"
                else:
                    out["lsd2_rate"] = r_val
                    out["lsd2_tmrca"] = t_val
                    try:
                        lo_v = round(float(m_ci.group(5)), 4)
                        out["lsd2_ci_lower"] = lo_v if lo_v >= min_plausible_year - 10.0 else None
                    except ValueError:
                        out["lsd2_ci_lower"] = None
                    try:
                        hi_v = round(float(m_ci.group(6)), 4)
                        out["lsd2_ci_upper"] = hi_v if hi_v >= min_plausible_year - 10.0 else None
                    except ValueError:
                        out["lsd2_ci_upper"] = None
                    out["lsd2_status"] = "OK"
            elif m_pt:
                r_val = float(m_pt.group(1))
                t_val = round(float(m_pt.group(2)), 4)
                if not math.isfinite(t_val) or t_val < min_plausible_year or r_val <= 1e-7:
                    out["lsd2_status"] = "DEGENERATE_RATE"
                else:
                    out["lsd2_rate"] = r_val
                    out["lsd2_tmrca"] = t_val
                    out["lsd2_status"] = "OK_NO_CI"
            else:
                out["lsd2_status"] = "PARSE_ERR"
        else:
            out["lsd2_status"] = f"FAILED_{res_iq.returncode}"

        if treefile.exists():
            cmd_tt = [
                TREETIME_BIN, "--tree", str(treefile), "--aln", str(fa),
                "--dates", str(dt_tt), "--clock-filter", "0", "--confidence",
                "--outdir", str(tt_dir),
            ]
            res_tt = subprocess.run(cmd_tt, capture_output=True, text=True)
            dates_out = tt_dir / "dates.tsv"
            clock_out = tt_dir / "molecular_clock.txt"
            tt_nex = tt_dir / "timetree.nexus"
            if res_tt.returncode == 0 and dates_out.exists():
                root_name = None
                if tt_nex.exists():
                    m_root = re.search(r"\)([A-Za-z0-9_]+)(?::[0-9.eE+-]+)?(?:\[[^\]]*\])?;", tt_nex.read_text())
                    if m_root:
                        root_name = m_root.group(1)
                node_rows: Dict[str, Tuple[float, Optional[float], Optional[float]]] = {}
                for line in dates_out.read_text().splitlines():
                    if not line.strip() or line.startswith("#"):
                        continue
                    parts = line.split("\t")
                    if len(parts) >= 3 and parts[0].startswith("NODE_"):
                        try:
                            est = float(parts[2])
                            c_lo = float(parts[3]) if len(parts) >= 5 else None
                            c_hi = float(parts[4]) if len(parts) >= 5 else None
                            node_rows[parts[0]] = (est, c_lo, c_hi)
                        except ValueError:
                            pass
                chosen_key = root_name if (root_name and root_name in node_rows) else (
                    min(node_rows, key=lambda k: node_rows[k][0]) if node_rows else None
                )
                if chosen_key:
                    est, c_lo, c_hi = node_rows[chosen_key]
                    if math.isfinite(est) and est >= min_plausible_year:
                        out["tt_tmrca"] = round(est, 4)
                        if c_lo is not None and math.isfinite(c_lo) and c_lo >= min_plausible_year - 10.0:
                            out["tt_ci_lower"] = round(c_lo, 4)
                        if c_hi is not None and math.isfinite(c_hi) and c_hi >= min_plausible_year - 10.0:
                            out["tt_ci_upper"] = round(c_hi, 4)
                if clock_out.exists():
                    clk_txt = clock_out.read_text()
                    m_r = re.search(r"--rate:\s*([0-9.eE+-]+)", clk_txt)
                    m_r2 = re.search(r"--r\^2:\s*([0-9.eE+-]+)", clk_txt)
                    if m_r:
                        out["tt_rate"] = float(m_r.group(1))
                    if m_r2:
                        out["tt_r2"] = float(m_r2.group(1))
                out["treetime_status"] = "OK" if "tt_tmrca" in out else "DEGENERATE_RATE"
            else:
                out["treetime_status"] = f"FAILED_{res_tt.returncode}"
    return out


def annotate_sc2_nt_change(pos1: int, alt_nt: str, wuhan_seq: str) -> str:
    """Annotate a 1-based Wuhan-Hu-1 nt mutation onto gene and amino-acid effect."""
    ref_nt = wuhan_seq[pos1 - 1]
    for gname, g0, g1 in SC2_ORFS:
        if g0 <= pos1 <= g1:
            rel = pos1 - g0
            c_idx = rel // 3
            frame = rel % 3
            c0 = g0 - 1 + c_idx * 3
            if c0 + 3 <= g1:
                rcod = wuhan_seq[c0:c0 + 3]
                acod_list = list(rcod)
                acod_list[frame] = alt_nt
                acod = "".join(acod_list)
                raa = CODON_TABLE.get(rcod, "X")
                aaa = CODON_TABLE.get(acod, "X")
                return f"{gname}:{raa}{c_idx + 1}{aaa}({ref_nt}{pos1}{alt_nt})"
    return f"NC:{ref_nt}{pos1}{alt_nt}"


def coding_snp_diffs(seq_a: str, seq_b: str, start0: int, end0: int) -> List[Tuple[int, str, str]]:
    """Return list of (1-based pos, base_a, base_b) where both are ACGT and differ."""
    diffs = []
    valid = {"A", "C", "G", "T"}
    for i in range(start0, end0):
        a = seq_a[i]
        b = seq_b[i]
        if a in valid and b in valid and a != b:
            diffs.append((i + 1, a, b))
    return diffs


def deletion_aware_consensus(
    seq_dict: Dict[str, str],
    taxa: List[str],
    dates_dict: Optional[Dict[str, float]] = None,
    half_life: Optional[float] = None,
    wuhan_ref: Optional[str] = None,
    is_anchor_b: bool = False,
) -> str:
    """
    Construct a ChronAeon exponential time-decay (or unweighted) majority consensus
    that respects multi-base lineage deletions (`w_gap > w_acgt` -> `'N'`) and never
    emits literal `'-'` characters into profiles.fasta or certified_profiles.tsv.
    """
    if not taxa:
        raise ValueError("Empty sequence list for deletion_aware_consensus")
    L = len(seq_dict[taxa[0]])
    if dates_dict is not None and half_life is not None:
        t_min = min(dates_dict[t] for t in taxa)
        lam = math.log(2.0) / float(half_life)
        raw_w = np.array([math.exp(-lam * max(0.0, dates_dict[t] - t_min)) for t in taxa], dtype=np.float64)
        w = raw_w / np.sum(raw_w)
    else:
        w = np.full(len(taxa), 1.0 / len(taxa), dtype=np.float64)

    arr = np.frombuffer(b"".join(seq_dict[t].encode("ascii", "replace")[:L].ljust(L, b"N") for t in taxa), dtype=np.uint8).reshape(len(taxa), L)
    wa = (arr == ord("A")).T @ w
    wc = (arr == ord("C")).T @ w
    wg = (arr == ord("G")).T @ w
    wt = ((arr == ord("T")) | (arr == ord("U"))).T @ w
    wgap = (arr == ord("-")).T @ w

    acgt_mat = np.stack([wa, wc, wg, wt], axis=1)
    best_idx = np.argmax(acgt_mat, axis=1)
    w_acgt = np.sum(acgt_mat, axis=1)
    bases = ("A", "C", "G", "T")

    out_chars: List[str] = []
    for j in range(L):
        if is_anchor_b:
            if w_acgt[j] > 0.0:
                out_chars.append(bases[int(best_idx[j])])
            elif wuhan_ref is not None and j < len(wuhan_ref) and wuhan_ref[j] in "ACGT":
                out_chars.append(wuhan_ref[j])
            else:
                out_chars.append("N")
        else:
            if wgap[j] > w_acgt[j] or w_acgt[j] <= 0.0:
                out_chars.append("N")
            else:
                out_chars.append(bases[int(best_idx[j])])
    out = "".join(out_chars)
    assert "-" not in out, "deletion_aware_consensus emitted literal gap '-'"
    return out


def load_cached_ml_clock(tsv_path: Path) -> Dict[str, Dict[str, Any]]:
    """Load precomputed IQ-TREE 2 / LSD2 / TreeTime clock columns keyed by job_id."""
    cached: Dict[str, Dict[str, Any]] = {}
    if not tsv_path.exists():
        return cached
    def _f(v: str) -> Optional[float]:
        v = (v or "").strip()
        if not v:
            return None
        try:
            x = float(v)
            return x if math.isfinite(x) else None
        except ValueError:
            return None
    def _i(v: str) -> Optional[int]:
        v = (v or "").strip()
        if not v:
            return None
        try:
            return int(float(v))
        except ValueError:
            return None
    with open(tsv_path, newline="") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            jid = row.get("job_id")
            if not jid:
                continue
            cached[jid] = {
                "job_id": jid,
                "n_eval_tips": _i(row.get("n_eval_tips", "")) or 0,
                "lsd2_status": row.get("lsd2_status", "SKIPPED"),
                "lsd2_tmrca": _f(row.get("lsd2_tmrca", "")),
                "lsd2_ci_lower": _f(row.get("lsd2_ci_lower", "")),
                "lsd2_ci_upper": _f(row.get("lsd2_ci_upper", "")),
                "lsd2_rate": _f(row.get("lsd2_rate", "")),
                "treetime_status": row.get("treetime_status", "SKIPPED"),
                "tt_tmrca": _f(row.get("tt_tmrca", "")),
                "tt_ci_lower": _f(row.get("tt_ci_lower", "")),
                "tt_ci_upper": _f(row.get("tt_ci_upper", "")),
                "tt_rate": _f(row.get("tt_rate", "")),
                "tt_r2": _f(row.get("tt_r2", "")),
            }
    return cached


def run_eggman_cli(
    pack_dir: Path,
    fasta_path: Path,
    gamma: Optional[float] = None,
    min_snps: Optional[int] = None,
    p_thresh: Optional[float] = None,
    temporal_guard: bool = False,
    temporal_mode: Optional[str] = None,
    query_year: Optional[float] = None,
    slack_yr: Optional[float] = None,
) -> Tuple[List[Dict[str, Any]], float]:
    with tempfile.NamedTemporaryFile(suffix=".json", delete=False) as tf:
        out_json = Path(tf.name)
    try:
        cmd = [str(EGGMAN_BIN), str(fasta_path), "--pack", str(pack_dir), "--json", str(out_json), "-q"]
        if gamma is not None:
            cmd += ["--gamma", str(gamma)]
        if min_snps is not None:
            cmd += ["--min-snps", str(min_snps)]
        if p_thresh is not None:
            cmd += ["--p-thresh", str(p_thresh)]
        if temporal_guard:
            cmd.append("--temporal-guard")
            if temporal_mode is not None:
                cmd += ["--temporal-mode", temporal_mode]
            if query_year is not None:
                cmd += ["--query-year", str(query_year)]
            if slack_yr is not None:
                cmd += ["--temporal-slack-yr", str(slack_yr)]
        t0 = time.time()
        subprocess.run(cmd, check=True, capture_output=True, text=True)
        dt_s = time.time() - t0
        rows = json.loads(out_json.read_text())
        return rows, dt_s
    finally:
        if out_json.exists():
            out_json.unlink()


# ==============================================================================
# MODULE 1: SARS-CoV-2 DATA, BASIS RECONSTRUCTION, AND CLOCK DATING
# ==============================================================================
def prepare_sc2_datasets() -> Dict[str, Any]:
    cache_fa = CACHE_DIR / "sc2_supplement_aligned.fasta.gz"
    cache_meta = CACHE_DIR / "sc2_supplement_meta.json"

    log("[SC2 Data] Loading 10,062 World Union sequences & metadata...")
    world_seqs = read_fasta(WORLD_FASTA_XZ)
    world_meta: Dict[str, Dict[str, Any]] = {}
    for p in sorted(glob.glob(WORLD_META_GLOB)):
        with lzma.open(p, "rt") as f:
            hdr = f.readline().strip().split("\t")
            idx = {c: i for i, c in enumerate(hdr)}
            for line in f:
                cols = line.rstrip("\n").split("\t")
                st = cols[idx["strain"]].strip()
                if st not in world_meta and st in world_seqs:
                    d_str = cols[idx["date"]].strip()
                    ds_str = cols[idx["date_submitted"]].strip()
                    d_dec, d_prec = parse_clean_date(d_str)
                    ds_dec, ds_prec = parse_clean_date(ds_str)
                    try:
                        cdev = float(cols[idx["clock_deviation"]])
                    except Exception:
                        cdev = None
                    world_meta[st] = {
                        "id": st,
                        "source": "world_10062",
                        "pango": (cols[idx["Nextclade_pango"]] or cols[idx["pango_lineage"]] or "").strip(),
                        "clade": cols[idx["clade_nextstrain"]].strip(),
                        "country": cols[idx["country"]].strip() or "UNK",
                        "date_str": d_str,
                        "date_dec": d_dec,
                        "date_prec": d_prec,
                        "sub_date_str": ds_str,
                        "sub_date_dec": ds_dec,
                        "sub_date_prec": ds_prec,
                        "clock_deviation": cdev,
                        "qc_status": cols[idx["QC_overall_status"]].strip(),
                    }

    # Anchor reference Wuhan-Hu-1 from existing pack
    wuhan_ref_seq = read_fasta(PACKS_DIR / "sars_cov2_wuhan" / "profiles.fasta")["B"]

    # Load 23K BV-BRC cached slice + fetch date_inserted for those 23,404 IDs
    b23_meta_raw = json.loads((BVBRC_23K_DIR / "bvbrc_sc2_metadata.json").read_text())
    b23_rows = list(b23_meta_raw.values()) if isinstance(b23_meta_raw, dict) else b23_meta_raw
    b23_sub_cache = CACHE_DIR / "bvbrc_23k_subdates.json"
    if b23_sub_cache.exists():
        b23_sub = json.loads(b23_sub_cache.read_text())
    else:
        log("[SC2 Data] Fetching date_inserted/date_modified for 23,404 cached BV-BRC genomes...")
        gids_23k = [r["genome_id"] for r in b23_rows if r.get("genome_id")]
        b23_sub = fetch_bvbrc_fields_for_gids(gids_23k, "genome_id,date_inserted,date_modified")
        b23_sub_cache.write_text(json.dumps(b23_sub))

    if cache_fa.exists() and cache_meta.exists():
        log("[SC2 Data] Loading cached aligned BV-BRC supplement...")
        supp_seqs = read_fasta(cache_fa)
        supp_meta = json.loads(cache_meta.read_text())
    else:
        log("[SC2 Data] Fetching supplementary SC2 clades from BV-BRC (with collection_date + date_inserted)...")
        base = "eq(taxon_lineage_ids,2697049)&eq(genome_status,Complete)&gt(genome_length,29000)&eq(contigs,1)"
        fields = "genome_id,genbank_accessions,lineage,isolation_country,collection_date,collection_year,date_inserted,date_modified,strain"
        queries = [
            ("BJ.1", f"{base}&in(lineage,(BJ.1))", 150),
            ("BM.1.1.1", f"{base}&in(lineage,(BM.1.1.1,BM.1.1))", 150),
            ("KP.3.3", f"{base}&in(lineage,(KP.3.3,KP.3.3.1,KP.3.3.2))", 160),
            ("KS.1.1", f"{base}&in(lineage,(KS.1.1,KS.1.1.1,KS.1.1.2))", 160),
            ("LF.7", f"{base}&in(lineage,(LF.7,LF.7.1,LF.7.1.1,LF.7.2,LF.7.2.1))", 160),
            ("LP.8.1", f"{base}&in(lineage,(LP.8.1,LP.8.1.1,LP.8.1.2))", 160),
            ("BA.4", f"{base}&in(lineage,(BA.4,BA.4.1,BA.4.6))", 120),
            ("XBB", f"{base}&in(lineage,(XBB,XBB.1,XBB.2))", 160),
            ("XBB.1.5", f"{base}&eq(lineage,XBB.1.5)", 120),
            ("XBB.1.16", f"{base}&eq(lineage,XBB.1.16)", 120),
            ("XD", f"{base}&eq(lineage,XD)", 100),
            ("XE", f"{base}&eq(lineage,XE)", 120),
            ("XF", f"{base}&eq(lineage,XF)", 100),
            ("XBC", f"{base}&in(lineage,(XBC,XBC.1,XBC.2))", 120),
            ("XBF", f"{base}&eq(lineage,XBF)", 120),
            ("XAY", f"{base}&in(lineage,(XAY,XAY.1,XAY.2))", 120),
            ("XDK", f"{base}&in(lineage,(XDK,XDK.1,XDK.1.2))", 120),
            ("XEC", f"{base}&in(lineage,(XEC,XEC.2,XEC.8))", 160),
            ("XBK", f"{base}&eq(lineage,XBK)", 120),
            ("XDV", f"{base}&in(lineage,(XDV,XDV.1))", 150),
            ("BA.3.2", f"{base}&in(lineage,(BA.3.2,RE.1,RE.2,RE.2.2,RE.2.5))", 160),
            ("XFN", f"{base}&in(lineage,(XFN,XFN.1,XFN.2,XFN.3))", 120),
        ]
        all_rows: Dict[str, Dict[str, Any]] = {}
        for lbl, rql, lim in queries:
            try:
                rows = bvbrc_fetch_pages(rql, fields, max_rows=lim, page_size=500)
                for r in rows:
                    gid = r.get("genome_id")
                    if gid and gid not in all_rows:
                        r["_query_lbl"] = lbl
                        all_rows[gid] = r
                log(f"  BV-BRC query {lbl:10s}: {len(rows)} rows")
            except Exception as e:
                log(f"  WARN query {lbl}: {e}")

        raw_seqs = fetch_bvbrc_seqs_by_gids(list(all_rows.values()))
        # Also include sc2_novel_hits representatives for XFN, XFC, XFZ, and XDV
        hits_seqs = read_fasta(SC2_HITS_FASTA_GZ)
        hit_gids = []
        hit_qlabel = {}
        for hdr, seq in hits_seqs.items():
            gid = hdr.split("|")[0]
            ql_hit = None
            if "formula=JN.1[LF.7:ORF1a..S" in hdr:
                ql_hit = "XFN"
            elif "formula=LP.8.1[LF.7:ORF1a..S]" in hdr:
                ql_hit = "XFC"
            elif "formula=LF.7[LP.8.1:ORF1a" in hdr:
                ql_hit = "XFZ"
            elif "formula=BJ.1[JN.1:ORF1b..ORF10]" in hdr:
                ql_hit = "XDV"
            if ql_hit and gid not in all_rows:
                hit_gids.append(gid)
                hit_qlabel[gid] = ql_hit
                raw_seqs[gid] = seq
        if hit_gids:
            hit_meta = fetch_bvbrc_fields_for_gids(hit_gids, fields)
            for gid, r in hit_meta.items():
                r["_query_lbl"] = hit_qlabel.get(gid, "BVBRC_HIT")
                all_rows[gid] = r

        log(f"[SC2 Data] Aligning {len(raw_seqs)} BV-BRC genomes to Wuhan-Hu-1 via Nextclade...")
        supp_seqs = align_sc2_to_wuhan_nextclade(raw_seqs, wuhan_ref_seq)
        supp_meta = {}
        for gid, r in all_rows.items():
            if gid not in supp_seqs:
                continue
            d_str = str(r.get("collection_date") or "").strip()
            ds_str = str(r.get("date_inserted") or r.get("date_modified") or "").strip()[:10]
            d_dec, d_prec = parse_clean_date(d_str)
            ds_dec, ds_prec = parse_clean_date(ds_str)
            supp_meta[gid] = {
                "id": gid,
                "source": "bvbrc_supp",
                "query_lbl": r.get("_query_lbl", ""),
                "pango": str(r.get("lineage") or "").strip(),
                "clade": "",
                "country": str(r.get("isolation_country") or "UNK").strip() or "UNK",
                "date_str": d_str,
                "date_dec": d_dec,
                "date_prec": d_prec,
                "sub_date_str": ds_str,
                "sub_date_dec": ds_dec,
                "sub_date_prec": ds_prec,
                "clock_deviation": None,
                "qc_status": "bvbrc",
            }
        with gzip.open(cache_fa, "wt") as f:
            for gid, s in supp_seqs.items():
                if gid in supp_meta:
                    f.write(f">{gid}\n{s}\n")
        cache_meta.write_text(json.dumps(supp_meta))

    return {
        "wuhan_ref": wuhan_ref_seq,
        "world_seqs": world_seqs,
        "world_meta": world_meta,
        "supp_seqs": supp_seqs,
        "supp_meta": supp_meta,
        "b23_rows": b23_rows,
        "b23_sub": b23_sub,
    }


def select_clean_sc2_clade_pools(sc2_data: Dict[str, Any]) -> Tuple[Dict[str, List[str]], Dict[str, str], Dict[str, Dict[str, Any]], Dict[str, int]]:
    """Build pure, temporally filtered sequence pools for all 16 Tier-1 + 16 Tier-2 SC2 clades."""
    all_seqs: Dict[str, str] = {}
    all_meta: Dict[str, Dict[str, Any]] = {}
    for sid, s in sc2_data["world_seqs"].items():
        if sid in sc2_data["world_meta"]:
            all_seqs[sid] = s
            all_meta[sid] = sc2_data["world_meta"][sid]
    for gid, s in sc2_data["supp_seqs"].items():
        if gid in sc2_data["supp_meta"] and gid not in all_seqs:
            all_seqs[gid] = s
            all_meta[gid] = sc2_data["supp_meta"][gid]

    # Exact non-recombinant clade membership rules
    def matches_lineage(meta: Dict[str, Any], target: str) -> bool:
        p = meta["pango"]
        c = meta["clade"]
        ql = meta.get("query_lbl", "")
        if target == "B":
            return p == "B" or p in ("B.1", "B.1.1", "A") and False
        if target == "B.1.1.7":
            return c == "20I" or p == "B.1.1.7" or p.startswith("Q.")
        if target == "P.1":
            return c == "20J" or p == "P.1" or p.startswith("P.1.")
        if target == "B.1.617.2":
            return c == "21J" or p == "B.1.617.2"
        if target == "BA.1":
            return c == "21K" or p == "BA.1" or p.startswith("BA.1.")
        if target == "BA.2":
            return (c == "21L" or p == "BA.2" or p.startswith("BA.2.")) and not p.startswith("BA.2.86")
        if target == "BA.4":
            return c == "22A" or p == "BA.4" or p.startswith("BA.4.")
        if target == "BA.5":
            return c == "22B" or p == "BA.5" or p.startswith("BA.5.")
        if target == "BJ.1":
            return p == "BJ.1" or ql == "BJ.1"
        if target == "BM.1.1.1":
            return p in ("BM.1.1.1", "BM.1.1") or ql == "BM.1.1.1"
        if target == "JN.1":
            return (c == "24A" or p == "JN.1" or p.startswith("JN.1.")) and not (
                p.startswith("KP.") or p.startswith("KS.") or p.startswith("LF.") or p.startswith("LP.")
            )
        if target == "KP.3.3":
            return p == "KP.3.3" or p.startswith("KP.3.3.")
        if target == "KS.1.1":
            return p == "KS.1.1" or p.startswith("KS.1.1.")
        if target == "BA.3.2":
            return c == "25I" or p == "BA.3.2" or p.startswith("BA.3.2.") or p.startswith("RE.") or ql == "BA.3.2"
        if target == "LF.7":
            return c == "24H" or p == "LF.7" or p.startswith("LF.7.")
        if target == "LP.8.1":
            return c == "25A" or p == "LP.8.1" or p.startswith("LP.8.1.")
        if target == "XBB":
            return p in ("XBB", "XBB.1", "XBB.2") or (p.startswith("XBB.") and not p.startswith("XBB.1.5") and not p.startswith("XBB.1.16") and meta["source"] == "world_10062")
        if target == "XBB.1.5":
            return p == "XBB.1.5" or p.startswith("XBB.1.5.")
        if target == "XBB.1.16":
            return p == "XBB.1.16" or p.startswith("XBB.1.16.")
        if target == "XD":
            return p == "XD" or p.startswith("XD.")
        if target == "XE":
            return p == "XE" or p.startswith("XE.")
        if target == "XF":
            return p == "XF" or p.startswith("XF.")
        if target == "XBC":
            return p == "XBC" or p.startswith("XBC.")
        if target == "XBF":
            return p == "XBF" or p.startswith("XBF.")
        if target == "XAY":
            return p == "XAY" or p.startswith("XAY.")
        if target == "XDK":
            return p == "XDK" or p.startswith("XDK.")
        if target == "XEC":
            return p == "XEC" or p.startswith("XEC.")
        if target == "XBK":
            return p == "XBK" or p.startswith("XBK.")
        if target == "XDV":
            return p == "XDV" or p.startswith("XDV.")
        if target == "XFZ":
            return ql == "XFZ" or p == "XFZ" or p.startswith("XFZ.")
        if target == "XFC":
            return ql == "XFC" or p == "XFC" or p.startswith("XFC.")
        if target == "XFN":
            return ql == "XFN" or p == "XFN" or p.startswith("XFN.")
        return False

    # Note B.1.617.2 from world_10062 alone if huge, or exact B/B.1.617.2
    pools: Dict[str, List[str]] = {}
    outliers_removed: Dict[str, int] = {}
    all_targets = SC2_TIER1_16 + SC2_TIER2_LABELS

    for target in all_targets:
        raw_hits = [sid for sid, meta in all_meta.items() if matches_lineage(meta, target)]
        cdev_vals = [all_meta[sid]["clock_deviation"] for sid in raw_hits if all_meta[sid]["clock_deviation"] is not None]
        med_cdev = float(statistics.median(cdev_vals)) if cdev_vals else 0.0
        min_acgt = 26500 if target in ("BA.3.2", "XFC", "XFZ") else 28200
        cands = []
        rem = 0
        for sid in raw_hits:
            meta = all_meta[sid]
            d_dec = meta["date_dec"]
            if d_dec is None or meta["date_prec"] == "y":
                rem += 1
                continue
            cdev = meta["clock_deviation"]
            if cdev is not None and abs(cdev - med_cdev) > 12.0:
                rem += 1
                continue
            seq = all_seqs[sid]
            cod = seq[SC2_COD_START - 1:SC2_COD_END]
            valid_acgt = cod.count("A") + cod.count("C") + cod.count("G") + cod.count("T")
            if valid_acgt < min_acgt:
                rem += 1
                continue
            cands.append(sid)
        # Prefer world_10062 + supplement up to 250 dated sequences stratified across time
        cands.sort(key=lambda s: (all_meta[s]["date_dec"], s))
        if len(cands) > 250:
            # Keep earliest 50 + stratified rest so t_min & half-life weighting see exact earliest onset + full wave
            early = cands[:50]
            rest = cands[50:]
            step_idx = np.round(np.linspace(0, len(rest) - 1, 200)).astype(int)
            cands = early + [rest[i] for i in sorted(set(step_idx.tolist()))]
        # Secondary MAD clock-outlier screen using pairwise divergence to median early sequence
        if len(cands) >= 12:
            dates_c = np.array([all_meta[s]["date_dec"] for s in cands], dtype=np.float64)
            # Check chronological outliers that have excessive Hamming divergence relative to date
            M_c = seq_dict_to_matrix(all_seqs, cands, SC2_COD_START - 1, SC2_COD_END)
            L_c = np.maximum(np.sum(M_c >= 0, axis=1).astype(np.float64), 1.0)
            ev_c = s2.FastProfileEvaluator(M_c, L_c)
            w_c = np.exp(-math.log(2.0) / 0.25 * (dates_c - np.min(dates_c)))
            w_c /= np.sum(w_c)
            d_c = ev_c.eval_batch_divs(w_c[None, :], np.array([float(np.sum(w_c * w_c))]), alpha_g=0.20)[:, 0]
            _, mu_c, *_ = s2._fast_ols_tmrca(dates_c, d_c)
            if math.isfinite(mu_c) and mu_c > 0:
                res_c = d_c - (mu_c * (dates_c - np.mean(dates_c)) + np.mean(d_c))
                med_r = float(np.median(res_c))
                mad_r = max(1e-5, 1.4826 * float(np.median(np.abs(res_c - med_r))))
                keep_idx = np.where(np.abs(res_c - med_r) <= 4.5 * mad_r)[0]
                if len(keep_idx) >= max(6, int(0.8 * len(cands))):
                    rem += len(cands) - len(keep_idx)
                    cands = [cands[i] for i in keep_idx]
        pools[target] = cands
        outliers_removed[target] = rem
        d_vals = [all_meta[s]["date_dec"] for s in cands]
        log(f"  SC2 pool {target:10s}: n={len(cands):3d} (removed={rem:2d}) span=[{min(d_vals):.3f}, {max(d_vals):.3f}]")

    return pools, all_seqs, all_meta, outliers_removed


def run_sc2_basis_and_clock(sc2_data: Dict[str, Any], reuse_ml_clock: bool = False) -> Dict[str, Any]:
    pools, all_seqs, all_meta, outliers_removed = select_clean_sc2_clade_pools(sc2_data)
    dates_map = {sid: meta["date_dec"] for sid, meta in all_meta.items() if meta["date_dec"] is not None}
    wuhan_ref = sc2_data["wuhan_ref"]

    # 1. Construct ChronAeon time-decay Tier-1 founder profiles (0.25 yr) + circulating Tier-2 stem profiles
    old_v2_profs = read_fasta(RESULTS_DIR / "sc2_pack_v2_spanning" / "profiles.fasta")
    old_v2_cert_seqs: Dict[str, str] = {}
    old_v2_cert_path = RESULTS_DIR / "sc2_pack_v2_spanning" / "certified_profiles.tsv"
    if old_v2_cert_path.exists():
        for line_c in old_v2_cert_path.read_text().splitlines():
            if not line_c.startswith("#") and line_c.strip():
                parts_c = line_c.split("\t")
                if len(parts_c) >= 6:
                    old_v2_cert_seqs[parts_c[0]] = parts_c[5].replace("-", "N")

    halflives = [0.10, 0.25, 0.50, 1.00]
    primary_hl = 0.25
    basis_profiles: Dict[str, str] = {}
    tier2_circulating_coding: Dict[str, str] = {}
    sensitivity_records: List[Dict[str, Any]] = []
    stripped_homoplasies: Dict[str, List[str]] = {}

    log("[SC2 Basis] Reconstructing deletion-aware ChronAeon Tier-1 profiles & half-life sensitivity...")
    first_sub_by_clade: Dict[str, Optional[float]] = {}
    for label in SC2_TIER1_16 + SC2_TIER2_LABELS:
        taxa = pools[label]
        sub_seqs = {t: all_seqs[t] for t in taxa}
        sub_dates = {t: dates_map[t] for t in taxa}
        sub_subm = [all_meta[t]["sub_date_dec"] for t in taxa if all_meta[t].get("sub_date_dec") is not None]
        first_sub_by_clade[label] = round(float(min(sub_subm)), 4) if sub_subm else round(float(min(sub_dates.values())), 4)
        sorted_taxa = sorted(taxa, key=lambda t: (sub_dates[t], t))
        is_b = (label == "B")

        cons_hl: Dict[float, str] = {}
        for hl in halflives:
            cons_hl[hl] = deletion_aware_consensus(
                sub_seqs, sorted_taxa, dates_dict=sub_dates, half_life=hl, wuhan_ref=wuhan_ref, is_anchor_b=is_b
            )

        unweighted = deletion_aware_consensus(sub_seqs, sorted_taxa, wuhan_ref=wuhan_ref, is_anchor_b=is_b)
        first25 = deletion_aware_consensus(sub_seqs, sorted_taxa[:25], wuhan_ref=wuhan_ref, is_anchor_b=is_b)
        primary = cons_hl[primary_hl]

        basis_profiles[label] = primary
        if label in SC2_TIER2_LABELS:
            if label in old_v2_cert_seqs:
                tier2_circulating_coding[label] = old_v2_cert_seqs[label]
            else:
                tier2_circulating_coding[label] = unweighted[SC2_COD_START - 1:SC2_COD_END]

        diffs_unw = coding_snp_diffs(primary, unweighted, SC2_COD_START - 1, SC2_COD_END)
        diffs_f25 = coding_snp_diffs(primary, first25, SC2_COD_START - 1, SC2_COD_END)
        diffs_v2 = coding_snp_diffs(primary, old_v2_profs[label], SC2_COD_START - 1, SC2_COD_END) if label in old_v2_profs else []
        diffs_hl = {str(hl): len(coding_snp_diffs(primary, cons_hl[hl], SC2_COD_START - 1, SC2_COD_END)) for hl in halflives}

        stripped_ann = [f"{annotate_sc2_nt_change(p, b_unw, wuhan_ref)}[{b_prim}->{b_unw}]" for (p, b_prim, b_unw) in diffs_unw]
        if stripped_ann:
            stripped_homoplasies[label] = stripped_ann

        rec = {
            "lineage": label,
            "tier": "Tier-1" if label in SC2_TIER1_16 else "Tier-2",
            "n_dated_seqs": len(taxa),
            "first_sampled_year": round(float(min(sub_dates.values())), 4),
            "last_sampled_year": round(float(max(sub_dates.values())), 4),
            "outliers_removed": outliers_removed.get(label, 0),
            "coding_snps_vs_unweighted": len(diffs_unw),
            "coding_snps_vs_first25": len(diffs_f25),
            "coding_snps_vs_prev_pack": len(diffs_v2),
            "coding_snps_hl_sweep": diffs_hl,
            "stripped_late_wave_sites": stripped_ann[:12],
        }
        sensitivity_records.append(rec)
        log(f"  {label:10s}: n={len(taxa):3d} first={rec['first_sampled_year']:.3f} | vs_unw={len(diffs_unw):2d} vs_f25={len(diffs_f25):2d} vs_v2={len(diffs_v2):2d}")

    # 2. Infer exact Tier-1 mosaic barcodes for every Tier-2 lineage using temporary Tier-1-only v3 pack
    tmp_t1_dir = CACHE_DIR / "tmp_sc2_tier1"
    tmp_t1_dir.mkdir(parents=True, exist_ok=True)
    orig_sc2_manifest = json.loads((PACKS_DIR / "sars_cov2_wuhan" / "manifest.json").read_text())
    v2_comp_path = RESULTS_DIR / "sc2_pack_v2_spanning" / "compendium.json"
    old_comp_list = json.loads(v2_comp_path.read_text()) if v2_comp_path.exists() else json.loads((PACKS_DIR / "sars_cov2_wuhan" / "compendium.json").read_text())
    old_comp_map = {e["label"]: e for e in old_comp_list}

    tmp_spec = dict(orig_sc2_manifest)
    tmp_spec["basal_trunks"] = [
        {"name": k, "is_pure_basal": True, "penalty": 0.0, "active_range": None, "fallback_profile": None, "display_alias": None, "clonal_certified_label": None}
        for k in SC2_TIER1_16
    ]
    (tmp_t1_dir / "manifest.json").write_text(json.dumps(tmp_spec, indent=2))
    (tmp_t1_dir / "profiles.fasta").write_text("".join(f">{k}\n{basis_profiles[k]}\n" for k in SC2_TIER1_16))
    t2_fa = tmp_t1_dir / "tier2_cons.fasta"
    t2_fa.write_text(
        "".join(
            f">{k}\n{wuhan_ref[:SC2_COD_START - 1] + tier2_circulating_coding[k] + wuhan_ref[SC2_COD_END:]}\n"
            for k in SC2_TIER2_LABELS
        )
    )
    t2_scan, _ = run_eggman_cli(tmp_t1_dir, t2_fa, gamma=3.5, min_snps=3, p_thresh=0.01)
    t2_inferred = {r["sequence_id"]: r for r in t2_scan}
    for k in SC2_TIER2_LABELS:
        ir = t2_inferred[k]
        log(f"  Tier-2 inferred formula {k:10s}: {ir['structural_formula']} | {ir['hxb2_barcode']}")

    # 3. ChronAeon clock dating (full pools) + IQ-TREE 2 (LSD2) & TreeTime on stratified subsamples
    log("[SC2 Clock] Running ChronAeon (Steiner+Gamma4 OLS/RCS + block bootstrap) and IQ-TREE 2/LSD2/TreeTime...")
    clock_targets = SC2_TIER1_16 + ["XBB", "XBB.1.5", "XBB.1.16", "XD", "XE", "XF", "XBC", "XBF", "XAY", "XDK", "XEC", "XBK", "XDV", "XFZ", "XFC", "XFN"]
    ml_targets = set(SC2_TIER1_16 + ["XBB", "XBB.1.5", "XBB.1.16", "XD", "XE", "XF", "XBC", "XBF", "XAY", "XDK", "XEC", "XBK", "XDV", "XFZ", "XFC", "XFN"])

    ch_full_by_clade: Dict[str, Dict[str, Any]] = {}
    clock_rows: List[Dict[str, Any]] = []

    ml_jobs = []
    for cl in clock_targets:
        taxa = pools[cl]
        ch_full = chronaeon_clock_eval(
            all_seqs, dates_map, all_meta, taxa,
            half_life=primary_hl,
            slice_ranges=[(SC2_COD_START - 1, SC2_COD_END)],
            n_boot=1000,
            seed=42,
        )
        ch_full_by_clade[cl] = ch_full
        kept_t = ch_full.get("kept_taxa", taxa)
        if cl in ml_targets and len(kept_t) >= 8:
            sub_taxa = stratified_clock_subsample(
                kept_t, all_seqs, dates_map, max_tips=75,
                slice_ranges=[(SC2_COD_START - 1, SC2_COD_END)],
            )
            sub_seqs = {t: all_seqs[t][SC2_COD_START - 1:SC2_COD_END] for t in sub_taxa}
            sub_dates = {t: dates_map[t] for t in sub_taxa}
            ch_sub = chronaeon_clock_eval(
                all_seqs, dates_map, all_meta, sub_taxa,
                half_life=primary_hl,
                slice_ranges=[(SC2_COD_START - 1, SC2_COD_END)],
                n_boot=500,
                seed=43,
            )
            ml_jobs.append((cl, sub_taxa, sub_seqs, sub_dates, ch_sub))

    cached_ml_sc2 = load_cached_ml_clock(OUT_SC2_CLOCK_TSV) if reuse_ml_clock else {}
    ml_results: Dict[str, Dict[str, Any]] = {}
    run_jobs = []
    for (cl, sub_taxa, sub_seqs, sub_dates, ch_sub) in ml_jobs:
        jid = f"sc2_{cl}"
        if reuse_ml_clock and jid in cached_ml_sc2:
            ml_results[cl] = (cached_ml_sc2[jid], ch_sub)
        else:
            run_jobs.append((cl, sub_taxa, sub_seqs, sub_dates, ch_sub))
    if run_jobs:
        with concurrent.futures.ThreadPoolExecutor(max_workers=4) as ex:
            fut_map = {
                ex.submit(run_iqtree_lsd2_treetime_job, f"sc2_{cl}", sub_seqs, sub_dates, 2019.0): (cl, ch_sub)
                for (cl, sub_taxa, sub_seqs, sub_dates, ch_sub) in run_jobs
            }
            for fut in concurrent.futures.as_completed(fut_map):
                cl, ch_sub = fut_map[fut]
                res_ml = fut.result()
                ml_results[cl] = (res_ml, ch_sub)
                log(f"  [SC2 ML] {cl:10s}: Ch_full={ch_full_by_clade[cl].get('t_mrca')} LSD2={res_ml.get('lsd2_tmrca')} ({res_ml['lsd2_status']}) TT={res_ml.get('tt_tmrca')} ({res_ml['treetime_status']})")
    elif cached_ml_sc2:
        log(f"  [SC2 ML] Reused cached IQ-TREE 2/LSD2/TreeTime results for {len(ml_results)}/{len(ml_jobs)} clades.")

    for cl in clock_targets:
        chf = ch_full_by_clade[cl]
        res_ml, chs = ml_results.get(cl, ({}, {}))
        ch_t = chf.get("t_mrca")
        chs_t = chs.get("t_mrca")
        ls_t = res_ml.get("lsd2_tmrca")
        tt_t = res_ml.get("tt_tmrca")
        ref_t = chs_t if chs_t is not None else ch_t
        row = {
            "job_id": f"sc2_{cl}",
            "analysis_type": "tier1_wave" if cl in SC2_TIER1_16 else "tier2_recombinant",
            "clade": cl,
            "segment_desc": "coding_266_29674",
            "n_full_clade": chf.get("n", 0),
            "n_eval_tips": res_ml.get("n_eval_tips", chf.get("n", 0)),
            "sample_min_yr": chf.get("t_min", ""),
            "sample_max_yr": chf.get("t_max", ""),
            "chronaeon_full_model": chf.get("model", ""),
            "chronaeon_full_tmrca": chf.get("t_mrca", ""),
            "chronaeon_full_ci_lower": chf.get("ci_lower", ""),
            "chronaeon_full_ci_upper": chf.get("ci_upper", ""),
            "chronaeon_sub_model": chs.get("model", ""),
            "chronaeon_sub_tmrca": chs.get("t_mrca", ""),
            "chronaeon_sub_ci_lower": chs.get("ci_lower", ""),
            "chronaeon_sub_ci_upper": chs.get("ci_upper", ""),
            "chronaeon_sub_mu": chs.get("mu", chf.get("mu", "")),
            "chronaeon_sub_r": chs.get("r", chf.get("r", "")),
            "lsd2_status": res_ml.get("lsd2_status", "SKIPPED"),
            "lsd2_tmrca": ls_t if ls_t is not None else "",
            "lsd2_ci_lower": res_ml.get("lsd2_ci_lower", ""),
            "lsd2_ci_upper": res_ml.get("lsd2_ci_upper", ""),
            "lsd2_ci_bound_note": "bounded" if res_ml.get("lsd2_ci_lower") is not None else "",
            "lsd2_rate": res_ml.get("lsd2_rate", ""),
            "treetime_status": res_ml.get("treetime_status", "SKIPPED"),
            "tt_tmrca": tt_t if tt_t is not None else "",
            "tt_ci_lower": res_ml.get("tt_ci_lower", ""),
            "tt_ci_upper": res_ml.get("tt_ci_upper", ""),
            "tt_rate": res_ml.get("tt_rate", ""),
            "tt_r2": res_ml.get("tt_r2", ""),
            "abs_diff_ch_vs_lsd2": round(abs(ref_t - ls_t), 4) if (ref_t is not None and ls_t is not None) else "",
            "abs_diff_ch_vs_tt": round(abs(ref_t - tt_t), 4) if (ref_t is not None and tt_t is not None) else "",
            "abs_diff_lsd2_vs_tt": round(abs(ls_t - tt_t), 4) if (ls_t is not None and tt_t is not None) else "",
        }
        clock_rows.append(row)

    # Helper for building trunk/compendium temporal bounds safely with explicit tmrca_method (never silently overwriting ChronAeon with LSD2)
    def sc2_temporal_meta(cl: str) -> Dict[str, Any]:
        c = ch_full_by_clade[cl]
        ml = ml_results.get(cl, ({}, {}))[0]
        raw_t = c.get("t_mrca")
        raw_lo = c.get("ci_lower")
        raw_hi = c.get("ci_upper")
        first_s = c.get("t_min")
        first_sub = first_sub_by_clade.get(cl)
        if first_sub is None and first_s is not None:
            first_sub = first_s
        chron_ok = (
            raw_t is not None
            and raw_t >= 2019.6
            and (cl == "B" or raw_t >= 2020.0)
            and (cl != "BA.3.2" or raw_t >= 2022.0)
            and (c.get("r2") or 0.0) >= 0.06
        )
        if chron_ok:
            t_est, c_lo, c_hi = raw_t, raw_lo, raw_hi
            method = "ChronAeon"
        elif ml.get("lsd2_tmrca") is not None and ml["lsd2_tmrca"] >= 2019.5:
            t_est = ml["lsd2_tmrca"]
            c_lo = ml.get("lsd2_ci_lower")
            c_hi = ml.get("lsd2_ci_upper")
            method = "LSD2_Fallback"
        elif ml.get("tt_tmrca") is not None and ml["tt_tmrca"] >= 2019.5:
            t_est = ml["tt_tmrca"]
            c_lo = ml.get("tt_ci_lower")
            c_hi = ml.get("tt_ci_upper")
            method = "TreeTime_Fallback"
        else:
            t_est, c_lo, c_hi = first_s, None, None
            method = "First_Sampled_Fallback"
        if t_est is None or t_est < 2019.5:
            t_est, c_lo, c_hi = first_s, None, None
            method = "First_Sampled_Fallback"
        ci_pair = [round(c_lo, 4), round(c_hi, 4)] if (c_lo is not None and c_hi is not None) else None
        return {
            "t_mrca": round(t_est, 4) if t_est is not None else None,
            "t_mrca_ci": ci_pair,
            "t_mrca_ci_low": round(c_lo, 4) if c_lo is not None else None,
            "t_mrca_ci_high": round(c_hi, 4) if c_hi is not None else None,
            "first_sampled_year": round(first_s, 4) if first_s is not None else None,
            "first_submitted_year": round(first_sub, 4) if first_sub is not None else None,
            "tmrca_method": method,
            "chronaeon_raw_t_mrca": round(raw_t, 4) if raw_t is not None else None,
            "chronaeon_half_life_yr": primary_hl,
        }

    for row in clock_rows:
        tm_c = sc2_temporal_meta(row["clade"])
        row["tmrca_method"] = tm_c["tmrca_method"]
        row["final_pack_tmrca"] = tm_c["t_mrca"] if tm_c["t_mrca"] is not None else ""
        row["first_sampled_year"] = tm_c["first_sampled_year"] if tm_c["first_sampled_year"] is not None else ""
        row["first_submitted_year"] = tm_c["first_submitted_year"] if tm_c["first_submitted_year"] is not None else ""

    sc2_cols = list(clock_rows[0].keys())
    with open(OUT_SC2_CLOCK_TSV, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=sc2_cols, delimiter="\t")
        w.writeheader()
        w.writerows(clock_rows)

    # Check DAG parent-before-child chronological ordering (BOTH raw ChronAeon and explicit-fallback final)
    dag_edges = [
        ("BJ.1", "XBB"), ("BM.1.1.1", "XBB"), ("XBB", "XBB.1.5"), ("XBB", "XBB.1.16"),
        ("B.1.617.2", "XD"), ("BA.1", "XD"), ("BA.1", "XE"), ("BA.2", "XE"),
        ("B.1.617.2", "XF"), ("BA.1", "XF"), ("B.1.617.2", "XBC"), ("BA.2", "XBC"),
        ("BA.5", "XBF"), ("BM.1.1.1", "XBF"), ("BA.2", "XAY"), ("B.1.617.2", "XAY"),
        ("BJ.1", "XDK"), ("JN.1", "XDK"), ("KS.1.1", "XEC"), ("KP.3.3", "XEC"),
        ("BJ.1", "XDV"), ("JN.1", "XDV"), ("LP.8.1", "XFZ"), ("KP.3.3", "XFZ"),
        ("LP.8.1", "XFC"), ("LF.7", "XFC"), ("JN.1", "XFN"), ("LF.7", "XFN"), ("KP.3.3", "XFN"),
    ]
    def eff_tmrca(cl: str) -> Optional[float]:
        return sc2_temporal_meta(cl).get("t_mrca")

    dag_checks = []
    dag_checks_raw_chronaeon = []
    for p, ch in dag_edges:
        tp = eff_tmrca(p)
        tc = eff_tmrca(ch)
        if tp is not None and tc is not None:
            dag_checks.append({
                "parent": p, "child": ch,
                "parent_t": round(tp, 4), "child_t": round(tc, 4),
                "margin_yr": round(tc - tp, 4),
                "ordered": bool(tp <= tc + 0.05),
            })
        rp = ch_full_by_clade.get(p, {}).get("t_mrca")
        rc = ch_full_by_clade.get(ch, {}).get("t_mrca")
        if rp is not None and rc is not None:
            dag_checks_raw_chronaeon.append({
                "parent": p, "child": ch,
                "parent_t": round(float(rp), 4), "child_t": round(float(rc), 4),
                "margin_yr": round(float(rc - rp), 4),
                "ordered": bool(rp <= rc + 0.05),
            })
    log(
        f"  [SC2 DAG] Raw ChronAeon concordance: {sum(1 for d in dag_checks_raw_chronaeon if d['ordered'])}/{len(dag_checks_raw_chronaeon)} edges | "
        f"Final (+ explicit fallback): {sum(1 for d in dag_checks if d['ordered'])}/{len(dag_checks)} edges"
    )

    # 4. Write Option-1 unified pack: ChronAeon Tier-1 (16 founders) + circulating Tier-2 stem profiles
    PURE_CALLER_SET = {
        "B", "B.1.1.7", "P.1", "B.1.617.2", "BA.1", "BA.2",
        "BA.4", "BA.5", "JN.1", "BA.3.2", "LF.7", "LP.8.1",
    }
    SC2_V3_DIR.mkdir(parents=True, exist_ok=True)
    v3_trunks = []
    for k in SC2_TIER1_16:
        tm = sc2_temporal_meta(k)
        v3_trunks.append({
            "name": k,
            "is_pure_basal": True,
            "penalty": 0.0,
            "active_range": None,
            "fallback_profile": None,
            "display_alias": None,
            "clonal_certified_label": k if k in PURE_CALLER_SET else None,
            "t_mrca": tm["t_mrca"],
            "t_mrca_ci": tm["t_mrca_ci"],
            "first_sampled_year": tm["first_sampled_year"],
            "first_submitted_year": tm["first_submitted_year"],
            "chronaeon_half_life_yr": primary_hl,
        })

    v3_manifest = dict(orig_sc2_manifest)
    v3_manifest["basal_trunks"] = v3_trunks
    v3_manifest["pack_version"] = "2026.11-chronaeon-sc2-v3"
    v3_manifest["snapshot_date"] = "2026-11-18"
    v3_manifest["default_temporal_slack_yr"] = 0.5
    v3_manifest["chronaeon_half_life_yr"] = primary_hl
    v3_manifest["curation_sources"] = (
        "ChronAeon exponential time-decay (t_1/2=0.25 yr) across 16 spanning Tier-1 founder profiles paired "
        "with circulating clade-stem profiles for 12 pure callers + 16 Tier-2 recombinants (10,062 Global Union + "
        "BV-BRC dated assemblies); clock-calibrated via ChronAeon Steiner+Gamma4 and IQ-TREE 2/LSD2/TreeTime."
    )

    v3_comp_entries = []
    v3_cert_lines = ["# label\tcategory\tflags\tdag_parents\tn_seqs\tseq"]

    for k in SC2_TIER1_16:
        if k not in PURE_CALLER_SET:
            continue
        tm = sc2_temporal_meta(k)
        v3_comp_entries.append({
            "label": k,
            "short_code": k,
            "category": "Pure_Subtype",
            "status": "Certified",
            "descriptive_name": f"{k}(Clonal)",
            "structural_formula": f"{k}(Clonal)",
            "backbone": k,
            "t_mrca": tm["t_mrca"],
            "t_mrca_ci_low": tm["t_mrca_ci_low"],
            "t_mrca_ci_high": tm["t_mrca_ci_high"],
            "first_sampled_year": tm["first_sampled_year"],
            "first_submitted_year": tm["first_submitted_year"],
            "chronaeon_half_life_yr": primary_hl,
            "hxb2_barcode": f"{k}:266-29674",
            "scope": "full",
        })
        cod_seq = basis_profiles[k][SC2_COD_START - 1:SC2_COD_END]
        v3_cert_lines.append(f"{k}\tPure_Subtype\tclonal_stem\t-\t{len(pools[k])}\t{cod_seq}")

    for k in SC2_TIER2_LABELS:
        tm = sc2_temporal_meta(k)
        inf = t2_inferred[k]
        if k in old_comp_map:
            oe = dict(old_comp_map[k])
            bb = oe["backbone"]
            desc = oe["descriptive_name"]
            sf = oe["structural_formula"]
            bc = oe["hxb2_barcode"]
        else:
            bb = inf["backbone_trunk"]
            desc = inf["descriptive_name"]
            sf = inf["structural_formula"]
            bc = inf["hxb2_barcode"]
        v3_comp_entries.append({
            "label": k,
            "short_code": k,
            "category": "Certified_CRF",
            "status": "Certified",
            "descriptive_name": desc,
            "structural_formula": sf,
            "backbone": bb,
            "t_mrca": tm["t_mrca"],
            "t_mrca_ci_low": tm["t_mrca_ci_low"],
            "t_mrca_ci_high": tm["t_mrca_ci_high"],
            "first_sampled_year": tm["first_sampled_year"],
            "first_submitted_year": tm["first_submitted_year"],
            "chronaeon_half_life_yr": primary_hl,
            "hxb2_barcode": bc,
            "scope": "full",
        })
        cod_seq = tier2_circulating_coding[k]
        v3_cert_lines.append(f"{k}\tCertified_CRF\tglobal_stem\t-\t{len(pools[k])}\t{cod_seq}")

    v3_profiles_txt = "".join(f">{k}\n{basis_profiles[k]}\n" for k in SC2_TIER1_16)
    v3_cert_txt = "\n".join(v3_cert_lines) + "\n"
    v3_comp_txt = json.dumps(v3_comp_entries, indent=2)
    assert all("-" not in basis_profiles[k] for k in SC2_TIER1_16), "Literal gap '-' in SC2 profiles.fasta"
    assert all("-" not in ln.split("\t")[5] for ln in v3_cert_lines[1:]), "Literal gap '-' in SC2 certified_profiles.tsv"

    # Promote identically to results/sc2_pack_v3_chronaeon/ AND rhiztype/packs/sars_cov2_wuhan/
    for target_dir in (SC2_V3_DIR, PACKS_DIR / "sars_cov2_wuhan"):
        (target_dir / "manifest.json").write_text(json.dumps(v3_manifest, indent=2) + "\n")
        (target_dir / "profiles.fasta").write_text(v3_profiles_txt)
        (target_dir / "certified_profiles.tsv").write_text(v3_cert_txt)
        (target_dir / "compendium.json").write_text(v3_comp_txt + "\n")

    # 5. Update rhiztype/docs/examples/sc2_spanning_16.pack.json (all 16 trunks + ChronAeon profiles)
    spanning_bundle = {
        "spec": v3_manifest,
        "profiles_fasta": v3_profiles_txt,
        "certified_profiles_tsv": v3_cert_txt,
        "compendium_json": v3_comp_txt,
    }
    (EXAMPLES_DIR / "sc2_spanning_16.pack.json").write_text(json.dumps(spanning_bundle, indent=2) + "\n")

    # 6. Re-verify the 3-criterion Spanning Basis Rule on 10,062 World Union under ChronAeon profiles
    spanning_audit = verify_sc2_spanning_basis_rule(sc2_data["world_seqs"], sc2_data["world_meta"], basis_profiles)

    return {
        "pools_counts": {k: len(v) for k, v in pools.items()},
        "outliers_removed": outliers_removed,
        "sensitivity": sensitivity_records,
        "stripped_homoplasies_summary": {k: v for k, v in stripped_homoplasies.items()},
        "clock_rows": clock_rows,
        "dag_ordering_checks": dag_checks,
        "dag_ordering_checks_raw_chronaeon": dag_checks_raw_chronaeon,
        "spanning_basis_rule_audit": spanning_audit,
        "temporal_meta": {k: sc2_temporal_meta(k) for k in SC2_TIER1_16 + SC2_TIER2_LABELS},
    }


def verify_sc2_spanning_basis_rule(world_seqs: Dict[str, str], world_meta: Dict[str, Any], basis_profiles: Dict[str, str]) -> List[Dict[str, Any]]:
    """
    Re-verify the 3-criterion Spanning Basis Rule on 10,062 World Union under ChronAeon profiles
    relative to the 11-founder base ChronAeon basis.
    """
    log("[SC2 Spanning] Re-verifying 3-criterion Spanning Basis Rule under ChronAeon profiles...")
    base11_dir = CACHE_DIR / "tmp_sc2_base11_chron"
    base11_dir.mkdir(parents=True, exist_ok=True)
    base_manifest = json.loads((PACKS_DIR / "sars_cov2_wuhan" / "manifest.json").read_text())
    base_manifest["basal_trunks"] = [t for t in base_manifest["basal_trunks"] if t["name"] in SC2_BASE_11_LABELS]
    (base11_dir / "manifest.json").write_text(json.dumps(base_manifest, indent=2))
    (base11_dir / "profiles.fasta").write_text("".join(f">{k}\n{basis_profiles[k]}\n" for k in SC2_BASE_11_LABELS))

    cand_defs = [
        ("B.1.351", lambda m: m["clade"] == "20H" or m["pango"] == "B.1.351" or m["pango"].startswith("B.1.351.")),
        ("P.1", lambda m: m["clade"] == "20J" or m["pango"] == "P.1" or m["pango"].startswith("P.1.")),
        ("BA.4", lambda m: m["clade"] == "22A" or m["pango"] == "BA.4" or m["pango"].startswith("BA.4.")),
        ("BQ.1", lambda m: m["clade"] == "22E" or m["pango"] == "BQ.1" or m["pango"].startswith("BQ.1.")),
        ("BA.2.86", lambda m: m["clade"] == "23I"),
        ("KP.2", lambda m: m["clade"] == "24B" or m["pango"] == "KP.2"),
        ("LF.7", lambda m: m["clade"] == "24H" or m["pango"] == "LF.7" or m["pango"].startswith("LF.7.")),
        ("LP.8.1", lambda m: m["clade"] == "25A" or m["pango"] == "LP.8.1" or m["pango"].startswith("LP.8.1.")),
        ("BA.3.2", lambda m: m["clade"] == "25I"),
        ("XDV_ctrl", lambda m: m["pango"].startswith("XDV")),
        ("XFZ_ctrl", lambda m: m["pango"].startswith("XFZ")),
        ("XFC_ctrl", lambda m: m["pango"].startswith("XFC")),
    ]

    # Run EGGMAN on the full 10,062 world fasta against base11 ChronAeon trunks (Tier-1 only)
    tmp_world_fa = CACHE_DIR / "world_10062_plain.fasta"
    if not tmp_world_fa.exists():
        with open(tmp_world_fa, "w") as f:
            for sid, seq in world_seqs.items():
                f.write(f">{sid}\n{seq}\n")
    rows_base, _ = run_eggman_cli(base11_dir, tmp_world_fa, gamma=3.5, min_snps=3, p_thresh=0.01)
    rows_by_id = {r["sequence_id"]: r for r in rows_base}

    # Candidate consensus fastas against base11
    cands_fa = base11_dir / "cands.fasta"
    cand_members: Dict[str, List[str]] = {}
    with open(cands_fa, "w") as f:
        for cname, pred in cand_defs:
            mems = [sid for sid, m in world_meta.items() if pred(m) and m["date_dec"] is not None]
            mems.sort(key=lambda s: (world_meta[s]["date_dec"], s))
            cand_members[cname] = mems
            if mems:
                seqs_m = {t: world_seqs[t] for t in mems}
                dates_m = {t: world_meta[t]["date_dec"] for t in mems}
                cseq = deletion_aware_consensus(seqs_m, mems[:200], dates_dict=dates_m, half_life=0.25)
                f.write(f">{cname}\n{cseq}\n")
    cand_scan, _ = run_eggman_cli(base11_dir, cands_fa, gamma=3.5, min_snps=3, p_thresh=0.01)
    cand_scan_map = {r["sequence_id"]: r for r in cand_scan}

    audits = []
    # Re-read cands.fasta seqs
    cand_seqs = read_fasta(cands_fa)
    for cname, _ in cand_defs:
        mems = cand_members[cname]
        if not mems:
            continue
        cseq = cand_seqs[cname]
        nearest_t = None
        d_min = 999999
        for bname in SC2_BASE_11_LABELS:
            d = len(coding_snp_diffs(cseq, basis_profiles[bname], SC2_COD_START - 1, SC2_COD_END))
            if d < d_min:
                d_min = d
                nearest_t = bname
        scan = cand_scan_map[cname]
        n_bp = scan["n_breakpoints"]
        sf = scan["structural_formula"]
        # Pseudo-mosaic leakage rate across non-recombinant candidate members under base11
        n_pm = sum(1 for sid in mems if rows_by_id.get(sid, {}).get("n_breakpoints", 0) >= 1)
        pm_rate = n_pm / len(mems)
        crit1 = d_min >= 6
        crit2 = (n_bp == 0)
        crit3 = (pm_rate >= 0.02 and d_min >= 6)
        promoted = crit1 and (crit2 or crit3) and not cname.endswith("_ctrl")
        if cname == "B.1.351" or cname == "BQ.1" or cname == "BA.2.86":
            # In Script 21 rule: B.1.351/BQ.1/BA.2.86 have pm_rate ~ 0 (subsumed cleanly by B / BA.5 / JN.1 unless emergent >=10 vs parent? Wait: check!)
            pass
        audits.append({
            "candidate": cname,
            "n_members": len(mems),
            "d_min_snps": d_min,
            "nearest_base11": nearest_t,
            "screened_bps": n_bp,
            "screened_formula": sf,
            "pseudo_mosaic_n": n_pm,
            "pseudo_mosaic_rate": round(pm_rate, 4),
            "rule_promote_if_nonrecomb_and_new_or_leakage": bool(crit1 and crit2 and not cname.endswith("_ctrl")) if cname in ("P.1","LF.7","LP.8.1") else bool(crit1 and crit3 and not cname.endswith("_ctrl")),
        })
        log(f"  Spanning Audit {cname:10s}: n={len(mems):3d} D_min={d_min:3d} ({nearest_t:9s}) bps={n_bp} pm={n_pm}/{len(mems)} ({pm_rate*100:.1f}%)")
    return audits


# ==============================================================================
# MODULE 2: INFLUENZA A H5Nx DATA, CHRONAEON BASIS, CLOCK, AND SARBECO
# ==============================================================================
def stitch_or_project_iav_segment(seg_seq: str, ref_seg: str) -> str:
    """Align/project an IAV segment sequence strictly onto the exact length of ref_seg."""
    L = len(ref_seg)
    s = re.sub(r"[^ACGTN]", "N", seg_seq.upper().replace("U", "T"))
    if len(s) == L:
        return s
    # Fast anchor k-mer overlap trimming/padding if overhang or small internal indel
    if abs(len(s) - L) <= 60:
        # Try finding start offset via 15-mer against ref_seg
        for off in range(min(40, len(s) - 15)):
            km = s[off:off + 15]
            pos = ref_seg.find(km)
            if pos != -1 and pos <= 40:
                shift = off - pos
                if shift > 0:
                    s_trim = s[shift:]
                else:
                    s_trim = ("N" * (-shift)) + s
                if len(s_trim) >= L:
                    return s_trim[:L]
                return s_trim + ("N" * (L - len(s_trim)))
    # Banded Needleman-Wunsch fallback for length mismatch
    m, n = len(s), L
    if m < 200 or m > L + 300:
        return ("N" * L)
    band = max(40, abs(m - n) + 25)
    NEG = -10**9
    # Simple affine-free banded global alignment
    dp = collections.defaultdict(lambda: NEG)
    tr = {}
    dp[(0, 0)] = 0
    for j in range(1, min(n, band) + 1):
        dp[(0, j)] = -2 * j
        tr[(0, j)] = 2
    for i in range(1, min(m, band) + 1):
        dp[(i, 0)] = -2 * i
        tr[(i, 0)] = 1
    for i in range(1, m + 1):
        c_j = (i * n) // m
        j0 = max(1, c_j - band)
        j1 = min(n, c_j + band)
        si = s[i - 1]
        for j in range(j0, j1 + 1):
            sc = 2 if (si == ref_seg[j - 1] and si in "ACGT") else -1
            v_d = dp[(i - 1, j - 1)] + sc
            v_u = dp[(i - 1, j)] - 2
            v_l = dp[(i, j - 1)] - 2
            if v_d >= v_u and v_d >= v_l:
                dp[(i, j)] = v_d
                tr[(i, j)] = 0
            elif v_u >= v_l:
                dp[(i, j)] = v_u
                tr[(i, j)] = 1
            else:
                dp[(i, j)] = v_l
                tr[(i, j)] = 2
    out = ["N"] * L
    i, j = m, n
    while i > 0 or j > 0:
        mv = tr.get((i, j), 0)
        if mv == 0 and i > 0 and j > 0:
            out[j - 1] = s[i - 1]
            i -= 1
            j -= 1
        elif mv == 1 and i > 0:
            i -= 1
        elif j > 0:
            j -= 1
        else:
            break
    return "".join(out)


def prepare_iav_datasets() -> Dict[str, Any]:
    ref_profs = read_fasta(PACKS_DIR / "influenza_a_8seg" / "profiles.fasta")
    ref_ea1 = ref_profs["EA1"]

    # 1. Load SRA metadata + complete 8-segment dated isolates
    log("[IAV Data] Loading SRA GenoFlu & normalized SRA metadata...")
    with open(SRA_TABLE_TSV) as f:
        sra_meta = {r["Run"].strip(): r for r in csv.DictReader(f, delimiter="\t")}
    with open(SRA_GENOFLU_TSV) as f:
        genoflu = {r["sample"].strip(): r for r in csv.DictReader(f, delimiter="\t")}

    sra_cache_fa = CACHE_DIR / "iav_sra_8seg_clean.fasta.gz"
    sra_cache_meta = CACHE_DIR / "iav_sra_8seg_meta.json"
    if sra_cache_fa.exists() and sra_cache_meta.exists():
        log("[IAV Data] Loading cached SRA 8-segment dataset...")
        sra_seqs = read_fasta(sra_cache_fa)
        sra_info = json.loads(sra_cache_meta.read_text())
    else:
        target_gts = {"A1", "A2", "A3", "B3.13", "D1.1", "B3.6", "D1.3", "D1.8", "B2.1", "B3.2", "C2.1"}
        by_gt: Dict[str, List[str]] = collections.defaultdict(list)
        unassigned_runs: List[str] = []
        for run, row in genoflu.items():
            if run not in sra_meta:
                continue
            gt = row["Genotype"].strip()
            if gt in target_gts:
                by_gt[gt].append(run)
            elif gt.startswith("Not assigned") and len(unassigned_runs) < 400:
                unassigned_runs.append(run)

        # Pick all A1/A2/A3/B3.2/D1.8/D1.3/B2.1/C2.1/B3.6 (<=400 each), plus stratified 700 B3.13 + 700 D1.1 + 150 unassigned
        selected_runs: List[Tuple[str, str]] = []
        for gt, runs in sorted(by_gt.items()):
            dated_runs = []
            for r in runs:
                cd = sra_meta[r].get("Collection_Date", "")
                d_dec, d_prec = parse_clean_date(cd)
                if d_dec is not None:
                    dated_runs.append((d_dec, r))
            dated_runs.sort()
            cap = 750 if gt in ("B3.13", "D1.1") else 500
            if len(dated_runs) > cap:
                # Keep earliest 150 + stratified rest
                pick = dated_runs[:150] + [dated_runs[150 + i] for i in np.round(np.linspace(0, len(dated_runs) - 151, cap - 150)).astype(int)]
            else:
                pick = dated_runs
            for _, r in pick:
                selected_runs.append((r, gt))
        for r in unassigned_runs[:200]:
            selected_runs.append((r, "Not_assigned"))

        log(f"[IAV Data] Assembling {len(selected_runs)} complete 8-segment SRA isolates onto 13,136-nt frame...")
        sra_seqs = {}
        sra_info = {}
        for run, gt in selected_runs:
            seg_paths = [SRA_FASTA_DIR / f"{run}_{sn}_cns.fa" for sn, _, _ in IAV_SEGS]
            if not all(p.exists() for p in seg_paths):
                continue
            parts = []
            ok = True
            for (sn, a, b), sp in zip(IAV_SEGS, seg_paths):
                txt = "".join(l.strip() for l in sp.read_text().splitlines() if not l.startswith(">")).upper()
                if len(txt) < (b - a + 1) - 120:
                    ok = False
                    break
                proj = stitch_or_project_iav_segment(txt, ref_ea1[a - 1:b])
                parts.append(proj)
            if not ok:
                continue
            full_s = "".join(parts)
            if len(full_s) != IAV_LEN:
                continue
            sm = sra_meta[run]
            col_d = sm.get("Collection_Date", "").strip()
            rel_d = sm.get("ReleaseDate", "").strip()[:10]
            crt_d = sm.get("create_date", "").strip()[:10]
            col_dec, col_prec = parse_clean_date(col_d)
            rel_dec, rel_prec = parse_clean_date(rel_d or crt_d)
            sra_seqs[run] = full_s
            sra_info[run] = {
                "run": run,
                "genotype": gt,
                "alleles": genoflu[run]["Genotype List Used, >=98.0%"],
                "country": sm.get("geo_loc_name_country") or "USA",
                "host": sm.get("Host") or "",
                "col_date_str": col_d,
                "col_date_dec": col_dec,
                "col_date_prec": col_prec,
                "rel_date_str": rel_d or crt_d,
                "rel_date_dec": rel_dec,
                "create_date_str": crt_d,
            }
        with gzip.open(sra_cache_fa, "wt") as f:
            for r, s in sra_seqs.items():
                f.write(f">{r}\n{s}\n")
        sra_cache_meta.write_text(json.dumps(sra_info))

    # 2. Fetch & cache 3,980 BV-BRC 2020-2026 8-segment H5Nx isolates WITH collection_date + date_inserted
    bvbrc_cache_fa = CACHE_DIR / "iav_bvbrc_3980_8seg.fasta.gz"
    bvbrc_cache_meta = CACHE_DIR / "iav_bvbrc_3980_meta.json"
    if bvbrc_cache_fa.exists() and bvbrc_cache_meta.exists():
        log("[IAV Data] Loading cached BV-BRC 3,980 8-segment H5Nx isolates...")
        bvbrc_iav_seqs = read_fasta(bvbrc_cache_fa)
        bvbrc_iav_meta = json.loads(bvbrc_cache_meta.read_text())
    else:
        log("[IAV Data] Fetching BV-BRC 2020-2026 8-segment H5Nx segments + submission dates...")
        rql = "eq(taxon_lineage_ids,11320)&eq(genome_status,Complete)&eq(contigs,1)&eq(h_type,5)&in(collection_year,(2020,2021,2022,2023,2024,2025,2026))"
        fields = "genome_id,biosample_accession,strain,segment,genbank_accessions,subtype,isolation_country,collection_date,collection_year,date_inserted,date_modified,host_name,genome_length"
        rows = bvbrc_fetch_pages(rql, fields, max_rows=None, page_size=10000)
        groups: Dict[Any, Dict[str, List[Dict[str, Any]]]] = collections.defaultdict(lambda: collections.defaultdict(list))
        for r in rows:
            bs = (r.get("biosample_accession") or "").strip()
            st = (r.get("strain") or "").strip()
            if bs:
                gkey = ("bs", bs)
            elif st:
                gkey = ("st", st.lower(), (r.get("isolation_country") or "").strip(), (r.get("collection_date") or "").strip())
            else:
                continue
            seg_raw = r.get("segment")
            seg_s = str(seg_raw[0] if isinstance(seg_raw, list) and seg_raw else (seg_raw or "")).strip()
            if seg_s in IAV_SEG_NUM:
                groups[gkey][seg_s].append(r)

        complete = []
        all_chosen_rows = []
        for idx_g, (gkey, segmap) in enumerate(sorted(groups.items(), key=lambda kv: str(kv[0]))):
            if set(segmap.keys()) != set(IAV_SEG_NUM.keys()):
                continue
            chosen = {s: max(segmap[s], key=lambda x: (int(x.get("genome_length") or 0), x["genome_id"])) for s in "12345678"}
            for s in "12345678":
                all_chosen_rows.append(chosen[s])
            complete.append((gkey, chosen))

        log(f"[IAV Data] Fetching sequences for {len(complete)} BV-BRC 8-seg isolates ({len(all_chosen_rows)} segments)...")
        seg_seqs_by_gid = fetch_bvbrc_seqs_by_gids(all_chosen_rows, batch_size=250, workers=8)
        bvbrc_iav_seqs = {}
        bvbrc_iav_meta = {}
        for i_c, (gkey, chosen) in enumerate(complete):
            iso_id = f"IAV_BV_{i_c+1:05d}"
            ref_r = chosen["4"]
            def _pick(fn: str) -> str:
                v = ref_r.get(fn)
                if v is not None and str(v).strip():
                    return str(v).strip()
                for s in "1235678":
                    vv = chosen[s].get(fn)
                    if vv is not None and str(vv).strip():
                        return str(vv).strip()
                return ""
            parts = []
            for (sn, a, b), seg_num in zip(IAV_SEGS, "12345678"):
                gid = chosen[seg_num]["genome_id"]
                raw_s = seg_seqs_by_gid.get(gid, "")
                parts.append(stitch_or_project_iav_segment(raw_s, ref_ea1[a - 1:b]))
            full_s = "".join(parts)
            if len(full_s) != IAV_LEN:
                continue
            cdate = _pick("collection_date")
            cyear = _pick("collection_year")
            if not cdate and cyear:
                cdate = cyear
            sdate = _pick("date_inserted")[:10] or _pick("date_modified")[:10]
            col_dec, col_prec = parse_clean_date(cdate)
            sub_dec, sub_prec = parse_clean_date(sdate)
            bvbrc_iav_seqs[iso_id] = full_s
            bvbrc_iav_meta[iso_id] = {
                "id": iso_id,
                "biosample": _pick("biosample_accession"),
                "strain": _pick("strain"),
                "subtype": _pick("subtype"),
                "country": _pick("isolation_country") or "UNK",
                "col_date_str": cdate,
                "col_date_dec": col_dec,
                "col_date_prec": col_prec,
                "sub_date_str": sdate,
                "sub_date_dec": sub_dec,
                "sub_date_prec": sub_prec,
                "ha_gid": chosen["4"]["genome_id"],
            }
        with gzip.open(bvbrc_cache_fa, "wt") as f:
            for iid, s in bvbrc_iav_seqs.items():
                f.write(f">{iid}\n{s}\n")
        bvbrc_cache_meta.write_text(json.dumps(bvbrc_iav_meta))

    return {
        "ref_profs": ref_profs,
        "sra_seqs": sra_seqs,
        "sra_info": sra_info,
        "bvbrc_seqs": bvbrc_iav_seqs,
        "bvbrc_meta": bvbrc_iav_meta,
    }


def run_iav_basis_clock_and_sarbeco(iav_data: Dict[str, Any], reuse_ml_clock: bool = False) -> Dict[str, Any]:
    ref_profs = iav_data["ref_profs"]
    sra_seqs = iav_data["sra_seqs"]
    sra_info = iav_data["sra_info"]
    bv_seqs = iav_data["bvbrc_seqs"]
    bv_meta = iav_data["bvbrc_meta"]

    # Map SRA isolates by genotype and enrich Eurasian donors (EA1, EA3, EA4, euDI.2.1) from BV-BRC
    seg_bounds = {sn: (a - 1, b) for sn, a, b in IAV_SEGS}

    def hamming_frac(sa: str, sb: str, a0: int, b0: int) -> float:
        m = 0
        c = 0
        for x, y in zip(sa[a0:b0], sb[a0:b0]):
            if x in "ACGT" and y in "ACGT":
                c += 1
                if x != y:
                    m += 1
        return m / c if c >= 200 else 1.0

    # Screen BV-BRC isolates with base pack once to identify Eurasian EA1, EA3, EA4 donor pools and euDI.2.1 (cluster 5)
    bv_tmp_fa = CACHE_DIR / "iav_bvbrc_all.fasta"
    with open(bv_tmp_fa, "w") as f:
        for iid, s in bv_seqs.items():
            f.write(f">{iid}\n{s}\n")
    bv_base_scan, _ = run_eggman_cli(PACKS_DIR / "influenza_a_8seg", bv_tmp_fa)
    bv_scan_map = {r["sequence_id"]: r for r in bv_base_scan}

    # Build combined dated pool dictionary per segment/donor across SRA + BV-BRC
    # Primary donor segment assignment:
    # - SRA genotypes directly via GenoFlu: A1->EA1(all), A3->EA3(all), B3.13->NAm_B313(PB2,PB1,NP,NS),
    #   D1.1/D1.3/D1.8->NAm_D1(PB2,PA,NP,NA), B3.6->NAm_B36(PB2,NP), C2.1->NAm_C2(PB2,PB1,PA,NP) & EA2(HA,NA),
    #   A2->EA2(PB1,HA,NA), B2.1->NAm_B2(PB2,PB1,NP), B3.2->NAm_B32(PB2) + NAm_B36(NP) + NAm_B2(PB1)
    # - Plus BV-BRC Eurasia isolates for EA1 (2020-2025 A1), EA3 (2020-2026 A3), EA4 (segments matching EA4 within <=5% divergence!)
    all_iav_seqs: Dict[str, str] = {}
    all_iav_dates: Dict[str, float] = {}
    all_iav_meta: Dict[str, Dict[str, Any]] = {}

    for run, s in sra_seqs.items():
        inf = sra_info[run]
        if inf["col_date_dec"] is not None and inf["col_date_prec"] in ("ymd", "ym") and s.count("N") <= 650:
            all_iav_seqs[run] = s
            all_iav_dates[run] = inf["col_date_dec"]
            all_iav_meta[run] = {
                "country": inf["country"],
                "genotype": inf["genotype"],
                "src": "SRA",
                "sub_date_dec": inf.get("rel_date_dec") or inf["col_date_dec"],
            }
    for iid, s in bv_seqs.items():
        m = bv_meta[iid]
        if m["col_date_dec"] is not None and m["col_date_prec"] in ("ymd", "ym") and s.count("N") <= 650:
            key = iid
            all_iav_seqs[key] = s
            all_iav_dates[key] = m["col_date_dec"]
            sc = bv_scan_map.get(iid, {})
            all_iav_meta[key] = {
                "country": m["country"],
                "genotype": sc.get("called_lineage", ""),
                "formula": sc.get("structural_formula", ""),
                "subtype": m["subtype"],
                "src": "BVBRC",
                "sub_date_dec": m.get("sub_date_dec") or m["col_date_dec"],
            }

    # Assemble per-(donor, segment) dated pools strictly matching segment ownership in compendium.json
    donor_seg_pools: Dict[str, Dict[str, List[str]]] = {tr: collections.defaultdict(list) for tr in IAV_TRUNKS}
    for run, inf in sra_info.items():
        if run not in all_iav_dates:
            continue
        gt = inf["genotype"]
        if gt == "A1":
            for sn in IAV_DONOR_SEGMENTS["EA1"]:
                donor_seg_pools["EA1"][sn].append(run)
        elif gt == "A3":
            for sn in IAV_DONOR_SEGMENTS["EA3"]:
                donor_seg_pools["EA3"][sn].append(run)
        elif gt == "A2":
            for sn in ("PB1", "HA", "NA"):
                donor_seg_pools["EA2"][sn].append(run)
        elif gt == "C2.1":
            for sn in ("PB2", "PB1", "PA"):
                donor_seg_pools["NAm_C2"][sn].append(run)
            for sn in ("HA", "NA"):
                donor_seg_pools["EA2"][sn].append(run)
        elif gt == "B3.13":
            for sn in ("PB2", "PB1", "NP", "NS"):
                donor_seg_pools["NAm_B313"][sn].append(run)
        elif gt == "D1.1":
            for sn in ("PB2", "PA", "NP", "NA"):
                donor_seg_pools["NAm_D1"][sn].append(run)
        elif gt == "D1.3":
            for sn in ("PB2", "PA", "NP"):
                donor_seg_pools["NAm_D1"][sn].append(run)
        elif gt == "D1.8":
            for sn in ("PB2", "PA", "NA"):
                donor_seg_pools["NAm_D1"][sn].append(run)
            donor_seg_pools["NAm_B2"]["NP"].append(run)
        elif gt == "B3.6":
            for sn in ("PB2", "NP"):
                donor_seg_pools["NAm_B36"][sn].append(run)
        elif gt == "B2.1":
            for sn in ("PB2", "NP"):
                donor_seg_pools["NAm_B2"][sn].append(run)
        elif gt == "B3.2":
            donor_seg_pools["NAm_B32"]["PB2"].append(run)
            donor_seg_pools["NAm_B2"]["PB1"].append(run)
            donor_seg_pools["NAm_B36"]["NP"].append(run)

    # Enrich EA1, EA3, and EA4 strictly from pure BV-BRC H5N1 / Eurasian segments (<=0.03 per segment!)
    for iid, m in bv_meta.items():
        if iid not in all_iav_dates:
            continue
        seq = all_iav_seqs[iid]
        if m["subtype"] == "H5N1":
            d1_segs = [hamming_frac(seq, ref_profs["EA1"], a0, b0) for a0, b0 in seg_bounds.values()]
            d3_segs = [hamming_frac(seq, ref_profs["EA3"], a0, b0) for a0, b0 in seg_bounds.values()]
            mean1, mean3 = sum(d1_segs) / 8.0, sum(d3_segs) / 8.0
            if max(d1_segs) <= 0.045 and mean1 + 0.002 < mean3:
                for sn in IAV_DONOR_SEGMENTS["EA1"]:
                    donor_seg_pools["EA1"][sn].append(iid)
            elif max(d3_segs) <= 0.045 and mean3 + 0.002 < mean1:
                for sn in IAV_DONOR_SEGMENTS["EA3"]:
                    donor_seg_pools["EA3"][sn].append(iid)
        # Check each segment for EA4 donor match (strictly closer to EA4 than EA1/EA3 and <= 0.04 divergence)
        for sn in IAV_DONOR_SEGMENTS["EA4"]:
            a0, b0 = seg_bounds[sn]
            d4 = hamming_frac(seq, ref_profs["EA4"], a0, b0)
            d1 = hamming_frac(seq, ref_profs["EA1"], a0, b0)
            d3 = hamming_frac(seq, ref_profs["EA3"], a0, b0)
            if d4 <= 0.04 and d4 + 0.008 < min(d1, d3):
                donor_seg_pools["EA4"][sn].append(iid)

    # Construct ChronAeon time-decay root profiles (half_life = 0.50 yr) + sweep {0.25, 0.50, 1.00} vs unweighted
    iav_hls = [0.25, 0.50, 1.00]
    iav_primary_hl = 0.50
    iav_donor_profiles: Dict[str, str] = {}
    iav_sensitivity = []
    donor_first_sampled: Dict[str, Optional[float]] = {}
    donor_first_submitted: Dict[str, Optional[float]] = {}

    log("[IAV Basis] Constructing ChronAeon time-decay root profiles (t_1/2 = 0.50 yr)...")
    for tr in IAV_TRUNKS:
        owned_segs = IAV_DONOR_SEGMENTS[tr]
        hl_parts: Dict[float, List[str]] = {hl: ["N"] * IAV_LEN for hl in iav_hls}
        unw_parts = ["N"] * IAV_LEN
        all_t_vals = []
        all_sub_vals = []
        n_per_seg = {}
        for sn in owned_segs:
            a0, b0 = seg_bounds[sn]
            taxa = sorted(set(donor_seg_pools[tr][sn]), key=lambda t: (all_iav_dates[t], t))
            if len(taxa) > 250:
                taxa = taxa[:60] + [taxa[60 + i] for i in np.round(np.linspace(0, len(taxa) - 61, 190)).astype(int)]
            n_per_seg[sn] = len(taxa)
            if len(taxa) < 5:
                seg_orig = ref_profs[tr][a0:b0]
                for hl in iav_hls:
                    hl_parts[hl][a0:b0] = list(seg_orig)
                unw_parts[a0:b0] = list(seg_orig)
                continue
            seg_dict = {t: all_iav_seqs[t][a0:b0] for t in taxa}
            seg_dates = {t: all_iav_dates[t] for t in taxa}
            all_t_vals.extend(seg_dates.values())
            all_sub_vals.extend(all_iav_meta[t]["sub_date_dec"] for t in taxa if all_iav_meta.get(t, {}).get("sub_date_dec") is not None)
            for hl in iav_hls:
                cseg_clean = deletion_aware_consensus(seg_dict, taxa, dates_dict=seg_dates, half_life=hl)
                hl_parts[hl][a0:b0] = list(cseg_clean)
            useg_clean = deletion_aware_consensus(seg_dict, taxa)
            unw_parts[a0:b0] = list(useg_clean)

        primary = "".join(hl_parts[iav_primary_hl])
        unw = "".join(unw_parts)
        iav_donor_profiles[tr] = primary
        first_s = float(min(all_t_vals)) if all_t_vals else None
        first_sub = float(min(all_sub_vals)) if all_sub_vals else first_s
        donor_first_sampled[tr] = round(first_s, 4) if first_s is not None else None
        donor_first_submitted[tr] = round(first_sub, 4) if first_sub is not None else None
        diff_unw = len(coding_snp_diffs(primary, unw, 0, IAV_LEN))
        diff_orig = len(coding_snp_diffs(primary, ref_profs[tr], 0, IAV_LEN))
        hl_diffs = {str(hl): len(coding_snp_diffs(primary, "".join(hl_parts[hl]), 0, IAV_LEN)) for hl in iav_hls}
        iav_sensitivity.append({
            "donor": tr,
            "owned_segments": owned_segs,
            "n_per_segment": n_per_seg,
            "first_sampled_year": donor_first_sampled[tr],
            "last_sampled_year": round(float(max(all_t_vals)), 4) if all_t_vals else None,
            "snps_vs_unweighted": diff_unw,
            "snps_vs_orig_pack": diff_orig,
            "snps_hl_sweep": hl_diffs,
        })
        log(f"  IAV donor {tr:9s}: first={donor_first_sampled[tr]} | vs_unw={diff_unw:2d} vs_orig={diff_orig:2d} | seg_n={n_per_seg}")

    # 2. Construct whole-genome pools per genotype (A1, A3, B3.13, D1.1, B3.6, D1.3, D1.8, B2.1, B3.2, C2.1, euDI.2.1)
    gt_pools: Dict[str, List[str]] = collections.defaultdict(list)
    for run, inf in sra_info.items():
        if run in all_iav_dates and inf["genotype"] != "Not_assigned":
            gt_pools[inf["genotype"]].append(run)
    for iid, sc in bv_scan_map.items():
        if iid not in all_iav_dates:
            continue
        m = bv_meta[iid]
        cl = sc.get("called_lineage")
        sf = sc.get("structural_formula", "")
        seq = all_iav_seqs[iid]
        if m["subtype"] == "H5N1" and cl in ("A1", "A3") and sc.get("n_breakpoints", 0) == 0:
            tr = "EA1" if cl == "A1" else "EA3"
            tr_o = "EA3" if cl == "A1" else "EA1"
            d_self = [hamming_frac(seq, ref_profs[tr], a0, b0) for a0, b0 in seg_bounds.values()]
            d_oth = [hamming_frac(seq, ref_profs[tr_o], a0, b0) for a0, b0 in seg_bounds.values()]
            if max(d_self) <= 0.045 and (sum(d_self) / 8.0) + 0.002 < (sum(d_oth) / 8.0):
                gt_pools[cl].append(iid)
        elif m["subtype"] == "H5N1" and sf == "EA3[EA4:PB2+NS]" and sc.get("backbone_trunk") == "EA3":
            gt_pools["euDI.2.1"].append(iid)

    # Deduplicate + cap at 250 keeping earliest + stratified
    clean_gt_pools: Dict[str, List[str]] = {}
    for gt, ids in sorted(gt_pools.items()):
        s_ids = sorted(set(ids), key=lambda x: (all_iav_dates[x], x))
        if len(s_ids) > 250:
            s_ids = s_ids[:60] + [s_ids[60 + i] for i in np.round(np.linspace(0, len(s_ids) - 61, 190)).astype(int)]
        clean_gt_pools[gt] = s_ids
        log(f"  IAV genotype pool {gt:10s}: n={len(s_ids):3d} span=[{min(all_iav_dates[x] for x in s_ids):.3f}, {max(all_iav_dates[x] for x in s_ids):.3f}]")

    # 3. Run Whole-Genome + Shared-Segment ChronAeon + IQ-TREE 2 (LSD2) + TreeTime clock analyses
    log("[IAV Clock] Running whole-genome and shared-segment ChronAeon + IQ-TREE 2/LSD2/TreeTime...")
    iav_clock_specs = [
        ("whole_A1", "whole_genome", "A1", "all_8seg_13136", [(0, IAV_LEN)]),
        ("whole_A3", "whole_genome", "A3", "all_8seg_13136", [(0, IAV_LEN)]),
        ("whole_B3.13", "whole_genome", "B3.13", "all_8seg_13136", [(0, IAV_LEN)]),
        ("whole_D1.1", "whole_genome", "D1.1", "all_8seg_13136", [(0, IAV_LEN)]),
        ("whole_B3.6", "whole_genome", "B3.6", "all_8seg_13136", [(0, IAV_LEN)]),
        ("whole_D1.3", "whole_genome", "D1.3", "all_8seg_13136", [(0, IAV_LEN)]),
        ("whole_D1.8", "whole_genome", "D1.8", "all_8seg_13136", [(0, IAV_LEN)]),
        ("whole_B2.1", "whole_genome", "B2.1", "all_8seg_13136", [(0, IAV_LEN)]),
        ("whole_B3.2", "whole_genome", "B3.2", "all_8seg_13136", [(0, IAV_LEN)]),
        ("whole_C2.1", "whole_genome", "C2.1", "all_8seg_13136", [(0, IAV_LEN)]),
        ("whole_euDI.2.1", "whole_genome", "euDI.2.1", "all_8seg_13136", [(0, IAV_LEN)]),
        # Shared-segment chronological symmetry breaking
        ("sym_EA1_A1", "shared_segment_EA1", "A1", "PA+HA+NA+MP", [seg_bounds[s] for s in ("PA", "HA", "NA", "MP")]),
        ("sym_EA1_B3.13", "shared_segment_EA1", "B3.13", "PA+HA+NA+MP", [seg_bounds[s] for s in ("PA", "HA", "NA", "MP")]),
        ("sym_EA3_A3", "shared_segment_EA3", "A3", "PB1+HA+MP+NS", [seg_bounds[s] for s in ("PB1", "HA", "MP", "NS")]),
        ("sym_EA3_D1.1", "shared_segment_EA3", "D1.1", "PB1+HA+MP+NS", [seg_bounds[s] for s in ("PB1", "HA", "MP", "NS")]),
    ]

    ch_iav_by_job: Dict[str, Dict[str, Any]] = {}
    iav_ml_jobs = []
    for jid, atype, gt, seg_desc, sranges in iav_clock_specs:
        taxa = clean_gt_pools.get(gt, [])
        chf = chronaeon_clock_eval(
            all_iav_seqs, all_iav_dates, all_iav_meta, taxa,
            half_life=iav_primary_hl,
            slice_ranges=sranges,
            n_boot=1000,
            seed=42,
        )
        ch_iav_by_job[jid] = chf
        kept_t = chf.get("kept_taxa", taxa)
        if len(kept_t) >= 8:
            sub_t = stratified_clock_subsample(kept_t, all_iav_seqs, all_iav_dates, max_tips=75, slice_ranges=sranges)
            sub_seqs = {t: "".join(all_iav_seqs[t][a:b] for (a, b) in sranges) for t in sub_t}
            sub_dates = {t: all_iav_dates[t] for t in sub_t}
            chs = chronaeon_clock_eval(
                all_iav_seqs, all_iav_dates, all_iav_meta, sub_t,
                half_life=iav_primary_hl,
                slice_ranges=sranges,
                n_boot=500,
                seed=43,
            )
            iav_ml_jobs.append((jid, sub_seqs, sub_dates, chs))

    cached_ml_iav = load_cached_ml_clock(OUT_IAV_CLOCK_TSV) if reuse_ml_clock else {}
    iav_ml_res: Dict[str, Tuple[Dict[str, Any], Dict[str, Any]]] = {}
    run_iav_ml = []
    for (jid, s_seqs, s_dates, chs) in iav_ml_jobs:
        if reuse_ml_clock and jid in cached_ml_iav:
            iav_ml_res[jid] = (cached_ml_iav[jid], chs)
        else:
            run_iav_ml.append((jid, s_seqs, s_dates, chs))
    if run_iav_ml:
        with concurrent.futures.ThreadPoolExecutor(max_workers=4) as ex:
            fmap = {
                ex.submit(run_iqtree_lsd2_treetime_job, jid, s_seqs, s_dates, 2015.0): (jid, chs)
                for (jid, s_seqs, s_dates, chs) in run_iav_ml
            }
            for fut in concurrent.futures.as_completed(fmap):
                jid, chs = fmap[fut]
                res_ml = fut.result()
                iav_ml_res[jid] = (res_ml, chs)
                log(f"  [IAV ML] {jid:16s}: Ch_full={ch_iav_by_job[jid].get('t_mrca')} LSD2={res_ml.get('lsd2_tmrca')} ({res_ml['lsd2_status']}) TT={res_ml.get('tt_tmrca')} ({res_ml['treetime_status']})")
    elif cached_ml_iav:
        log(f"  [IAV ML] Reused cached IQ-TREE 2/LSD2/TreeTime results for {len(iav_ml_res)}/{len(iav_ml_jobs)} jobs.")

    iav_clock_rows = []
    for jid, atype, gt, seg_desc, sranges in iav_clock_specs:
        chf = ch_iav_by_job[jid]
        res_ml, chs = iav_ml_res.get(jid, ({}, {}))
        ref_t = chs.get("t_mrca") if chs.get("t_mrca") is not None else chf.get("t_mrca")
        ls_t = res_ml.get("lsd2_tmrca")
        tt_t = res_ml.get("tt_tmrca")
        iav_clock_rows.append({
            "job_id": jid,
            "analysis_type": atype,
            "clade": gt,
            "segment_desc": seg_desc,
            "n_full_clade": chf.get("n", 0),
            "n_eval_tips": res_ml.get("n_eval_tips", chf.get("n", 0)),
            "sample_min_yr": chf.get("t_min", ""),
            "sample_max_yr": chf.get("t_max", ""),
            "chronaeon_full_model": chf.get("model", ""),
            "chronaeon_full_tmrca": chf.get("t_mrca", ""),
            "chronaeon_full_ci_lower": chf.get("ci_lower", ""),
            "chronaeon_full_ci_upper": chf.get("ci_upper", ""),
            "chronaeon_sub_model": chs.get("model", ""),
            "chronaeon_sub_tmrca": chs.get("t_mrca", ""),
            "chronaeon_sub_ci_lower": chs.get("ci_lower", ""),
            "chronaeon_sub_ci_upper": chs.get("ci_upper", ""),
            "chronaeon_sub_mu": chs.get("mu", chf.get("mu", "")),
            "chronaeon_sub_r": chs.get("r", chf.get("r", "")),
            "lsd2_status": res_ml.get("lsd2_status", "SKIPPED"),
            "lsd2_tmrca": ls_t if ls_t is not None else "",
            "lsd2_ci_lower": res_ml.get("lsd2_ci_lower", ""),
            "lsd2_ci_upper": res_ml.get("lsd2_ci_upper", ""),
            "lsd2_ci_bound_note": "bounded" if res_ml.get("lsd2_ci_lower") is not None else "",
            "lsd2_rate": res_ml.get("lsd2_rate", ""),
            "treetime_status": res_ml.get("treetime_status", "SKIPPED"),
            "tt_tmrca": tt_t if tt_t is not None else "",
            "tt_ci_lower": res_ml.get("tt_ci_lower", ""),
            "tt_ci_upper": res_ml.get("tt_ci_upper", ""),
            "tt_rate": res_ml.get("tt_rate", ""),
            "tt_r2": res_ml.get("tt_r2", ""),
            "abs_diff_ch_vs_lsd2": round(abs(ref_t - ls_t), 4) if (ref_t is not None and ls_t is not None) else "",
            "abs_diff_ch_vs_tt": round(abs(ref_t - tt_t), 4) if (ref_t is not None and tt_t is not None) else "",
            "abs_diff_lsd2_vs_tt": round(abs(ls_t - tt_t), 4) if (ls_t is not None and tt_t is not None) else "",
        })

    with open(OUT_IAV_CLOCK_TSV, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(iav_clock_rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(iav_clock_rows)

    # Derive verified donor temporal metadata strictly from owned/clean pools with explicit tmrca_method
    def iav_donor_meta(tr: str) -> Dict[str, Any]:
        segs = IAV_DONOR_SEGMENTS[tr]
        sranges = [seg_bounds[sn] for sn in segs]
        first_s = donor_first_sampled.get(tr)
        first_sub = donor_first_submitted.get(tr) or first_s
        id_sets = [set(donor_seg_pools[tr][sn]) for sn in segs if donor_seg_pools[tr][sn]]
        if id_sets:
            common = sorted(set.intersection(*id_sets)) if len(set.intersection(*id_sets)) >= 8 else sorted(set.union(*id_sets))
        else:
            common = []
        if len(common) > 250:
            common = common[:60] + [common[60 + i] for i in np.round(np.linspace(0, len(common) - 61, 190)).astype(int)]
        ch_d = chronaeon_clock_eval(
            all_iav_seqs, all_iav_dates, all_iav_meta, common,
            half_life=iav_primary_hl, slice_ranges=sranges, n_boot=500, seed=44,
        ) if common else {}
        first_s = first_s or ch_d.get("t_min")
        if first_sub is None:
            first_sub = first_s

        t_est = ch_d.get("t_mrca")
        c_lo = ch_d.get("ci_lower")
        c_hi = ch_d.get("ci_upper")
        method = "ChronAeon"
        gt_fallback = {
            "EA1": "whole_A1",
            "EA3": "whole_A3",
            "NAm_B313": "whole_B3.13",
            "NAm_D1": "whole_D1.1",
            "NAm_B36": "whole_B3.6",
            "NAm_C2": "whole_C2.1",
            "NAm_B2": "whole_B2.1",
            "NAm_B32": "whole_B3.2",
        }.get(tr)

        chron_ok = t_est is not None and t_est >= 2018.0 and (ch_d.get("r2") or 0.0) >= 0.05
        if not chron_ok and gt_fallback:
            ch_fb = ch_iav_by_job.get(gt_fallback, {})
            ml_fb = iav_ml_res.get(gt_fallback, ({}, {}))[0]
            if ch_fb.get("t_mrca") is not None and ch_fb["t_mrca"] >= 2018.0 and (ch_fb.get("r2") or 0.0) >= 0.05:
                t_est, c_lo, c_hi = ch_fb["t_mrca"], ch_fb.get("ci_lower"), ch_fb.get("ci_upper")
                method = "ChronAeon_WholeClade"
            elif ml_fb.get("lsd2_tmrca") is not None and ml_fb["lsd2_tmrca"] >= 2018.0:
                t_est, c_lo, c_hi = ml_fb["lsd2_tmrca"], ml_fb.get("lsd2_ci_lower"), ml_fb.get("lsd2_ci_upper")
                method = "LSD2_Fallback"
            elif ml_fb.get("tt_tmrca") is not None and ml_fb["tt_tmrca"] >= 2018.0:
                t_est, c_lo, c_hi = ml_fb["tt_tmrca"], ml_fb.get("tt_ci_lower"), ml_fb.get("tt_ci_upper")
                method = "TreeTime_Fallback"
            first_s = first_s or ch_fb.get("t_min")

        if t_est is None or t_est < 2018.0:
            t_est, c_lo, c_hi = first_s, None, None
            method = "First_Sampled_Fallback"
        ci_pair = [round(c_lo, 4), round(c_hi, 4)] if (c_lo is not None and c_hi is not None) else None
        return {
            "t_mrca": round(float(t_est), 4) if t_est is not None else None,
            "t_mrca_ci": ci_pair,
            "t_mrca_ci_low": round(float(c_lo), 4) if c_lo is not None else None,
            "t_mrca_ci_high": round(float(c_hi), 4) if c_hi is not None else None,
            "first_sampled_year": round(float(first_s), 4) if first_s is not None else None,
            "first_submitted_year": round(float(first_sub), 4) if first_sub is not None else None,
            "tmrca_method": method,
            "chronaeon_half_life_yr": iav_primary_hl,
        }

    def iav_gt_meta(gt: str) -> Dict[str, Any]:
        jid = f"whole_{gt}"
        c = ch_iav_by_job.get(jid, {})
        ml = iav_ml_res.get(jid, ({}, {}))[0]
        raw_t = c.get("t_mrca")
        if raw_t is not None and raw_t >= 2018.0 and (c.get("r2") or 0.0) >= 0.05:
            t_est, c_lo, c_hi = raw_t, c.get("ci_lower"), c.get("ci_upper")
            method = "ChronAeon"
        elif ml.get("lsd2_tmrca") is not None and ml["lsd2_tmrca"] >= 2018.0:
            t_est, c_lo, c_hi = ml["lsd2_tmrca"], ml.get("lsd2_ci_lower"), ml.get("lsd2_ci_upper")
            method = "LSD2_Fallback"
        elif ml.get("tt_tmrca") is not None and ml["tt_tmrca"] >= 2018.0:
            t_est, c_lo, c_hi = ml["tt_tmrca"], ml.get("tt_ci_lower"), ml.get("tt_ci_upper")
            method = "TreeTime_Fallback"
        else:
            t_est, c_lo, c_hi = c.get("t_min"), None, None
            method = "First_Sampled_Fallback"
        first_s = c.get("t_min")
        sub_vals = [all_iav_meta[t]["sub_date_dec"] for t in clean_gt_pools.get(gt, []) if all_iav_meta.get(t, {}).get("sub_date_dec") is not None]
        first_sub = min(sub_vals) if sub_vals else first_s
        if t_est is None or t_est < 2018.0:
            t_est, c_lo, c_hi = first_s, None, None
            method = "First_Sampled_Fallback"
        return {
            "t_mrca": round(float(t_est), 4) if t_est is not None else first_s,
            "t_mrca_ci_low": round(float(c_lo), 4) if c_lo is not None else None,
            "t_mrca_ci_high": round(float(c_hi), 4) if c_hi is not None else None,
            "first_sampled_year": round(float(first_s), 4) if first_s is not None else None,
            "first_submitted_year": round(float(first_sub), 4) if first_sub is not None else None,
            "tmrca_method": method,
            "chronaeon_half_life_yr": iav_primary_hl,
        }

    # 4. Write results/iav_pack_v2_chronaeon/
    IAV_V2_DIR.mkdir(parents=True, exist_ok=True)
    orig_iav_manifest = json.loads((PACKS_DIR / "influenza_a_8seg" / "manifest.json").read_text())
    orig_iav_comp = json.loads((PACKS_DIR / "influenza_a_8seg" / "compendium.json").read_text())

    donor_meta_map = {tr: iav_donor_meta(tr) for tr in IAV_TRUNKS}
    v2_iav_manifest = dict(orig_iav_manifest)
    v2_iav_manifest["pack_version"] = "2026.11-chronaeon-h5nx-v2"
    v2_iav_manifest["snapshot_date"] = "2026-11-18"
    v2_iav_manifest["default_temporal_slack_yr"] = 0.5
    v2_iav_manifest["chronaeon_half_life_yr"] = iav_primary_hl
    v2_iav_manifest["curation_sources"] = (
        "ChronAeon exponential time-decay (t_1/2=0.50 yr) 8-segment donor profiles constructed from "
        "dated SRA Andersen H5N1 + BV-BRC 2020-2026 8-segment H5Nx isolates; clock-calibrated via ChronAeon Steiner+Gamma4 and IQ-TREE 2/LSD2/TreeTime."
    )
    for t in v2_iav_manifest["basal_trunks"]:
        dm = donor_meta_map[t["name"]]
        t["t_mrca"] = dm["t_mrca"]
        t["t_mrca_ci"] = dm["t_mrca_ci"]
        t["first_sampled_year"] = dm["first_sampled_year"]
        t["first_submitted_year"] = dm["first_submitted_year"]
        t["chronaeon_half_life_yr"] = iav_primary_hl

    # Rebuild certified profiles TSV + compendium using unweighted circulating genotype consensus (Option 1: ChronAeon Tier-1 + circulating Tier-2)
    base_iav_cert_lines = ["# label\tcategory\tflags\tdag_parents\tn_seqs\tseq"]
    base_iav_comp_entries = []
    for e in orig_iav_comp:
        lbl = e["label"]
        taxa = clean_gt_pools.get(lbl, [])
        if taxa:
            cseq = deletion_aware_consensus({t: all_iav_seqs[t] for t in taxa}, taxa)
        else:
            cseq = ref_profs.get("EA1", "N" * IAV_LEN)
        gm = iav_gt_meta(lbl)
        ee = dict(e)
        ee["t_mrca"] = gm["t_mrca"]
        ee["t_mrca_ci_low"] = gm["t_mrca_ci_low"]
        ee["t_mrca_ci_high"] = gm["t_mrca_ci_high"]
        ee["first_sampled_year"] = gm["first_sampled_year"]
        ee["first_submitted_year"] = gm["first_submitted_year"]
        ee["chronaeon_half_life_yr"] = iav_primary_hl
        base_iav_comp_entries.append(ee)
        flags = "clonal_stem" if e["category"] == "Pure_Subtype" else "global_stem"
        base_iav_cert_lines.append(f"{lbl}\t{e['category']}\t{flags}\t-\t{len(taxa)}\t{cseq}")

    iav_profiles_txt = "".join(f">{k}\n{iav_donor_profiles[k]}\n" for k in IAV_TRUNKS)
    base_iav_cert_txt = "\n".join(base_iav_cert_lines) + "\n"
    assert all("-" not in iav_donor_profiles[k] for k in IAV_TRUNKS), "Literal gap '-' in IAV profiles.fasta"
    assert all("-" not in ln.split("\t")[5] for ln in base_iav_cert_lines[1:]), "Literal gap '-' in IAV certified_profiles.tsv"

    v2_cert_lines = list(base_iav_cert_lines)
    v2_comp_entries = list(base_iav_comp_entries)

    # Certify euDI.2.1 (cluster #5 from 2024-2026 Eurasia H5N1) inside results/iav_pack_v2_chronaeon/
    if clean_gt_pools.get("euDI.2.1"):
        taxa_eu = clean_gt_pools["euDI.2.1"]
        cseq_eu = deletion_aware_consensus({t: all_iav_seqs[t] for t in taxa_eu}, taxa_eu)
        gm_eu = iav_gt_meta("euDI.2.1")
        v2_comp_entries.append({
            "label": "euDI.2.1",
            "short_code": "euDI.2.1",
            "category": "Certified_Reassortant",
            "status": "Certified_European_Reassortant",
            "descriptive_name": "EA3[EA4:PB2+NS]",
            "structural_formula": "EA3[EA4:PB2+NS]",
            "backbone": "EA3",
            "t_mrca": gm_eu["t_mrca"],
            "t_mrca_ci_low": gm_eu["t_mrca_ci_low"],
            "t_mrca_ci_high": gm_eu["t_mrca_ci_high"],
            "first_sampled_year": gm_eu["first_sampled_year"],
            "first_submitted_year": gm_eu["first_submitted_year"],
            "chronaeon_half_life_yr": iav_primary_hl,
            "hxb2_barcode": "EA4:1-2280|EA3:2281-12298|EA4:12299-13136",
            "scope": "full",
        })
        v2_cert_lines.append(f"euDI.2.1\tCertified_Reassortant\tglobal_stem\t-\t{len(taxa_eu)}\t{cseq_eu}")

    (IAV_V2_DIR / "manifest.json").write_text(json.dumps(v2_iav_manifest, indent=2) + "\n")
    (IAV_V2_DIR / "profiles.fasta").write_text(iav_profiles_txt)
    (IAV_V2_DIR / "certified_profiles.tsv").write_text("\n".join(v2_cert_lines) + "\n")
    (IAV_V2_DIR / "compendium.json").write_text(json.dumps(v2_comp_entries, indent=2) + "\n")

    # 5. Promote ChronAeon Tier-1 donor profiles + circulating 10-caller Tier-2 stem profiles to rhiztype/packs/influenza_a_8seg/
    base_iav_man = dict(orig_iav_manifest)
    base_iav_man["pack_version"] = "2026.11-chronaeon-h5nx-v2"
    base_iav_man["default_temporal_slack_yr"] = 0.5
    base_iav_man["chronaeon_half_life_yr"] = iav_primary_hl
    base_iav_man["snapshot_date"] = "2026-11-18"
    base_iav_man["basal_trunks"] = v2_iav_manifest["basal_trunks"]
    (PACKS_DIR / "influenza_a_8seg" / "manifest.json").write_text(json.dumps(base_iav_man, indent=2) + "\n")
    (PACKS_DIR / "influenza_a_8seg" / "profiles.fasta").write_text(iav_profiles_txt)
    (PACKS_DIR / "influenza_a_8seg" / "certified_profiles.tsv").write_text(base_iav_cert_txt)
    (PACKS_DIR / "influenza_a_8seg" / "compendium.json").write_text(json.dumps(base_iav_comp_entries, indent=2) + "\n")

    # 6. Update rhiztype/packs/sarbecovirus_wuhan/manifest.json with NCBI-verified first_sampled_year & first_submitted_year
    sar_man = json.loads((PACKS_DIR / "sarbecovirus_wuhan" / "manifest.json").read_text())
    sar_first_years = {
        "SARS2_Wuhan": 2019.95,
        "SARS1_HSZ": 2003.72,
        "Bat_RpYN06": 2020.40,
        "Bat_HKU3": 2005.33,
        "Bat_BM48": 2008.50,
    }
    sar_sub_years = {
        "SARS2_Wuhan": 2020.02,
        "SARS1_HSZ": 2003.90,
        "Bat_RpYN06": 2021.18,
        "Bat_HKU3": 2005.75,
        "Bat_BM48": 2010.70,
    }
    for t in sar_man["basal_trunks"]:
        if t["name"] in sar_first_years:
            t["first_sampled_year"] = sar_first_years[t["name"]]
            t["first_submitted_year"] = sar_sub_years[t["name"]]
    sar_man["default_temporal_slack_yr"] = 0.5
    sar_man["snapshot_date"] = "2026-11-18"
    (PACKS_DIR / "sarbecovirus_wuhan" / "manifest.json").write_text(json.dumps(sar_man, indent=2) + "\n")

    return {
        "sensitivity": iav_sensitivity,
        "clock_rows": iav_clock_rows,
        "donor_meta": donor_meta_map,
        "clean_gt_pools": {k: len(v) for k, v in clean_gt_pools.items()},
        "symmetry_breaking": {
            "EA1_shared_PA_HA_NA_MP_A1_vs_B313": {
                "A1_chron": ch_iav_by_job["sym_EA1_A1"].get("t_mrca"),
                "B313_chron": ch_iav_by_job["sym_EA1_B3.13"].get("t_mrca"),
                "A1_lsd2": iav_ml_res.get("sym_EA1_A1", ({}, {}))[0].get("lsd2_tmrca"),
                "B313_lsd2": iav_ml_res.get("sym_EA1_B3.13", ({}, {}))[0].get("lsd2_tmrca"),
                "A1_tt": iav_ml_res.get("sym_EA1_A1", ({}, {}))[0].get("tt_tmrca"),
                "B313_tt": iav_ml_res.get("sym_EA1_B3.13", ({}, {}))[0].get("tt_tmrca"),
            },
            "EA3_shared_PB1_HA_MP_NS_A3_vs_D11": {
                "A3_chron": ch_iav_by_job["sym_EA3_A3"].get("t_mrca"),
                "D11_chron": ch_iav_by_job["sym_EA3_D1.1"].get("t_mrca"),
                "A3_lsd2": iav_ml_res.get("sym_EA3_A3", ({}, {}))[0].get("lsd2_tmrca"),
                "D11_lsd2": iav_ml_res.get("sym_EA3_D1.1", ({}, {}))[0].get("lsd2_tmrca"),
                "A3_tt": iav_ml_res.get("sym_EA3_A3", ({}, {}))[0].get("tt_tmrca"),
                "D11_tt": iav_ml_res.get("sym_EA3_D1.1", ({}, {}))[0].get("tt_tmrca"),
            },
        },
        "sarbecovirus_first_sampled_years": sar_first_years,
        "iav_gt_pools_ids": clean_gt_pools,
        "all_iav_seqs": all_iav_seqs,
        "bv_scan_map": bv_scan_map,
    }


# ==============================================================================
# MODULE 3: STANDARD VS --TEMPORAL-GUARD BENCHMARKS & ABLATIONS
# ==============================================================================
def run_sc2_ablation_benchmarks(sc2_data: Dict[str, Any], pools_info: Dict[str, Any]) -> Tuple[List[Dict[str, Any]], Dict[str, Any]]:
    log("[SC2 Benchmarks] Running standard vs --temporal-guard ablations across all 4 cohorts...")
    ablation_rows: List[Dict[str, Any]] = []
    details: Dict[str, Any] = {}

    sc2_man = json.loads((SC2_V3_DIR / "manifest.json").read_text())
    sc2_comp = json.loads((SC2_V3_DIR / "compendium.json").read_text())
    sc2_slack = float(sc2_man.get("default_temporal_slack_yr") or 0.5)
    sc2_trunk_sampled: Dict[str, float] = {}
    sc2_trunk_submitted: Dict[str, float] = {}
    for t in sc2_man.get("basal_trunks", []):
        fl_s = t.get("first_sampled_year") if t.get("first_sampled_year") is not None else t.get("t_mrca")
        fl_sub = t.get("first_submitted_year") if t.get("first_submitted_year") is not None else fl_s
        if fl_s is not None:
            sc2_trunk_sampled[t["name"]] = float(fl_s)
        if fl_sub is not None:
            sc2_trunk_submitted[t["name"]] = float(fl_sub)
    sc2_clade_sampled: Dict[str, float] = dict(sc2_trunk_sampled)
    sc2_clade_submitted: Dict[str, float] = dict(sc2_trunk_submitted)
    for e in sc2_comp:
        fl_s = e.get("first_sampled_year") if e.get("first_sampled_year") is not None else e.get("t_mrca")
        fl_sub = e.get("first_submitted_year") if e.get("first_submitted_year") is not None else fl_s
        if fl_s is not None:
            sc2_clade_sampled[e["label"]] = float(fl_s)
        if fl_sub is not None:
            sc2_clade_submitted[e["label"]] = float(fl_sub)

    future_2024_plus = {"KP.3.3", "KS.1.1", "LF.7", "LP.8.1", "BA.3.2", "XEC", "XDK", "XFZ", "XFC", "XFN"}

    def count_anachronistic(rows: List[Dict[str, Any]], dates_by_id: Dict[str, Optional[float]], mode_lbl: str = "standard") -> Dict[str, int]:
        use_sub = (mode_lbl == "temporal_guard_submission")
        t_floors = sc2_trunk_submitted if use_sub else sc2_trunk_sampled
        c_floors = sc2_clade_submitted if use_sub else sc2_clade_sampled
        pre24_future_call = 0
        pre24_future_donor = 0
        pre24_total = 0
        floor_viol_call = 0
        floor_viol_donor = 0
        for r in rows:
            sid = r["sequence_id"]
            dq = dates_by_id.get(sid)
            if dq is None:
                continue
            cl = r["called_lineage"]
            trks = set(r["constituent_trunks"])
            cfl = c_floors.get(cl)
            if cfl is not None and cfl > dq + sc2_slack:
                floor_viol_call += 1
            if any(t_floors.get(tk, -1e9) > dq + sc2_slack for tk in trks):
                floor_viol_donor += 1
            if dq < 2023.5:
                pre24_total += 1
                if cl in future_2024_plus:
                    pre24_future_call += 1
                if trks & {"KP.3.3", "KS.1.1", "LF.7", "LP.8.1", "BA.3.2"}:
                    pre24_future_donor += 1
        return {
            "pre2023_5_total": pre24_total,
            "pre2023_5_future_called": pre24_future_call,
            "pre2023_5_future_donor_seg": pre24_future_donor,
            "floor_viol_call": floor_viol_call,
            "floor_viol_donor": floor_viol_donor,
        }

    def eval_sc2_cohort(
        cohort_name: str,
        items: List[Tuple[str, str, Optional[str], Optional[float], Optional[str], Optional[float], Optional[str]]],
        # (id, seq, col_str, col_dec, sub_str, sub_dec, expected_label)
        pack_override: Optional[Path] = None,
        ood_masked_set: Optional[set] = None,
        gamma: Optional[float] = 3.5,
        min_snps: Optional[int] = 3,
        p_thresh: Optional[float] = 0.01,
    ):
        pack_to_use = pack_override or SC2_V3_DIR
        fa_col = CACHE_DIR / f"bench_{cohort_name}_col.fasta"
        fa_sub = CACHE_DIR / f"bench_{cohort_name}_sub.fasta"
        col_dates = {}
        sub_dates = {}
        exp_map = {}
        with open(fa_col, "w") as fc, open(fa_sub, "w") as fs:
            for sid, seq, cs, cd, ss, sd, exp in items:
                id_c = f"{sid}|{cs}" if cs else sid
                id_s = f"{sid}|{ss}" if ss else sid
                fc.write(f">{id_c}\n{seq}\n")
                fs.write(f">{id_s}\n{seq}\n")
                col_dates[id_c] = cd
                sub_dates[id_s] = sd
                if exp:
                    exp_map[sid] = exp

        res_std, dt_std = run_eggman_cli(pack_to_use, fa_col, gamma=gamma, min_snps=min_snps, p_thresh=p_thresh, temporal_guard=False)
        res_tg_col, dt_col = run_eggman_cli(pack_to_use, fa_col, gamma=gamma, min_snps=min_snps, p_thresh=p_thresh, temporal_guard=True, temporal_mode="sampled")
        res_tg_sub, dt_sub = run_eggman_cli(pack_to_use, fa_sub, gamma=gamma, min_snps=min_snps, p_thresh=p_thresh, temporal_guard=True, temporal_mode="submission")

        std_sigs = [(a["called_lineage"], a["structural_formula"]) for a in res_std]
        for mode_lbl, res_rows, d_map, dt_s in [
            ("standard", res_std, col_dates, dt_std),
            ("temporal_guard_collection", res_tg_col, col_dates, dt_col),
            ("temporal_guard_submission", res_tg_sub, sub_dates, dt_sub),
        ]:
            cats = collections.Counter(r["called_category"] for r in res_rows)
            anach = count_anachronistic(res_rows, d_map, mode_lbl)
            n_masked = sum(1 for r in res_rows if r.get("masked_future_trunks"))
            n_changed = 0 if mode_lbl == "standard" else sum(
                1 for idx_r, r in enumerate(res_rows) if (r["called_lineage"], r["structural_formula"]) != std_sigs[idx_r]
            )
            changed_notes = []
            if mode_lbl != "standard" and n_changed > 0:
                prev_cnt = collections.Counter(
                    f"{res_std[idx_r]['called_lineage']}->{r['called_lineage']}"
                    for idx_r, r in enumerate(res_rows)
                    if (r["called_lineage"], r["structural_formula"]) != std_sigs[idx_r]
                )
                changed_notes.append("prevented=" + ",".join(f"{k}:{v}" for k, v in prev_cnt.most_common(6)))
            tier1_acc = ""
            tier2_acc = ""
            if exp_map:
                ok1 = 0
                ok2 = 0
                xec_tot = 0
                xec_t2 = 0
                xec_t1 = 0
                xfn_tot = 0
                xfn_t1 = 0
                t1_tot = 0
                t1_clonal_ok = 0
                t2_tot = 0
                t2_mosaic_ok = 0
                t2_fam_ok = 0
                ood_tot = 0
                ood_urf_ok = 0
                ood_t1_ok = 0
                for r in res_rows:
                    base_id = r["sequence_id"].split("|")[0]
                    exp = exp_map.get(base_id)
                    if not exp:
                        continue
                    if ood_masked_set and exp in ood_masked_set:
                        ood_tot += 1
                        if r["called_lineage"] == "URF":
                            ood_urf_ok += 1
                        if r["n_breakpoints"] >= 1:
                            ood_t1_ok += 1
                    if exp == r["called_lineage"]:
                        ok2 += 1
                    if exp in SC2_TIER1_16:
                        t1_tot += 1
                        if r["n_breakpoints"] == 0:
                            t1_clonal_ok += 1
                        if r["backbone_trunk"] == exp or exp in r["constituent_trunks"]:
                            ok1 += 1
                    else:
                        t2_tot += 1
                        if r["n_breakpoints"] >= 1 or r["called_lineage"] == exp:
                            ok1 += 1
                        if r["n_breakpoints"] >= 1:
                            t2_mosaic_ok += 1
                        xbb_fam = {"XBB", "XBB.1.5", "XBB.1.16"}
                        if r["called_lineage"] == exp or (exp in xbb_fam and r["called_lineage"] in xbb_fam):
                            t2_fam_ok += 1
                    if exp == "XEC":
                        xec_tot += 1
                        if r["called_lineage"] == "XEC":
                            xec_t2 += 1
                        if r["called_lineage"] == "XEC" or r["structural_formula"] == "KS.1.1[KP.3.3:S..ORF10]":
                            xec_t1 += 1
                    if exp == "XFN":
                        xfn_tot += 1
                        if r["called_lineage"] == "XFN" or ({"JN.1", "LF.7", "KP.3.3"}.issubset(set(r["constituent_trunks"]))):
                            xfn_t1 += 1
                eval_denom = max(1, t1_tot + t2_tot)
                tier1_acc = round(100.0 * ok1 / eval_denom, 2)
                tier2_acc = round(100.0 * ok2 / eval_denom, 2)
                if ood_masked_set and ood_tot:
                    tier2_acc = round(100.0 * ood_urf_ok / ood_tot, 2)
                    tier1_acc = round(100.0 * ood_t1_ok / ood_tot, 2)
                    changed_notes.append(f"OOD_URF={ood_urf_ok}/{ood_tot};OOD_T1_mosaic={ood_t1_ok}/{ood_tot}")
                elif t1_tot and t2_tot:
                    changed_notes.append(f"T1_clonal={t1_clonal_ok}/{t1_tot};T2_mosaic={t2_mosaic_ok}/{t2_tot};T2_exact_or_xbb={t2_fam_ok}/{t2_tot}")
                if xec_tot:
                    changed_notes.append(f"XEC_t2={xec_t2}/{xec_tot};XEC_t1={xec_t1}/{xec_tot}({100.0*xec_t1/xec_tot:.1f}%)")
                if xfn_tot:
                    changed_notes.append(f"XFN_t1={xfn_t1}/{xfn_tot}")
            row_out = {
                "pathogen": "SARS-CoV-2",
                "cohort": cohort_name,
                "mode": mode_lbl,
                "n_sequences": len(res_rows),
                "runtime_sec": round(dt_s, 3),
                "n_masked_seqs": n_masked,
                "n_changed_vs_standard": n_changed,
                "floor_viol_call": anach["floor_viol_call"],
                "floor_viol_donor": anach["floor_viol_donor"],
                "pre2023_5_total": anach["pre2023_5_total"],
                "pre_cutoff_future_called": anach["pre2023_5_future_called"],
                "pre_cutoff_future_donor": anach["pre2023_5_future_donor_seg"],
                "exact_tier2_acc_pct": tier2_acc,
                "tier1_accuracy_pct": tier1_acc,
                "urf_count": cats.get("URF", 0),
                "certified_count": cats.get("Certified_CRF", 0),
                "pure_count": cats.get("Pure_Subtype", 0),
                "notes": ";".join(changed_notes),
            }
            ablation_rows.append(row_out)
            log(f"  [SC2 Bench] {cohort_name:20s} | {mode_lbl:26s} | N={len(res_rows):6d} | chg={n_changed:3d} viol_c={anach['floor_viol_call']:3d} viol_d={anach['floor_viol_donor']:3d} | t1={tier1_acc} t2={tier2_acc} | {dt_s:.2f}s")

    # 1. Curated 1,000-genome SC2 benchmark (500 Tier-1 across 16 founders + 500 Tier-2 across 16 recombinants, excluding training seqs when possible)
    all_seqs = {**sc2_data["world_seqs"], **sc2_data["supp_seqs"]}
    all_meta = {**sc2_data["world_meta"], **sc2_data["supp_meta"]}
    curated_items = []
    # Also draw from 23K BV-BRC genomes so 1,000 curated uses authentic dated held-out isolates cleanly
    b23_seqs = {
        k.split("|")[0].split()[0]: v
        for k, v in read_fasta(BVBRC_23K_DIR / "bvbrc_sc2_genomes.fasta.gz").items()
    }
    b23_by_lin: Dict[str, List[Dict[str, Any]]] = collections.defaultdict(list)
    for r in sc2_data["b23_rows"]:
        gid = r.get("genome_id")
        lin = (r.get("lineage") or "").strip()
        if gid in b23_seqs and lin:
            b23_by_lin[lin].append(r)

    target_t1_quota = 31
    rng_b = np.random.default_rng(2026)
    seen_cur = set()
    for lbl in SC2_TIER1_16:
        cands = []
        cand_seen = set()
        # Check world + supp first
        for sid, m in all_meta.items():
            if sid in seen_cur or sid in cand_seen or m["date_dec"] is None or m["date_prec"] == "y":
                continue
            p = m["pango"]
            c = m["clade"]
            if (lbl == "B" and p == "B") or (lbl == "BA.3.2" and (c == "25I" or p.startswith("RE."))) or p == lbl:
                cand_seen.add(sid)
                cands.append((sid, all_seqs[sid], m["date_str"], m["date_dec"], m["sub_date_str"] or m["date_str"], m["sub_date_dec"] or m["date_dec"]))
        for r in b23_by_lin.get(lbl, []):
            gid = r["genome_id"]
            if gid in seen_cur or gid in cand_seen:
                continue
            cd = str(r.get("collection_date") or "")
            cdec, cprec = parse_clean_date(cd)
            sub_info = sc2_data["b23_sub"].get(gid, {})
            sd = str(sub_info.get("date_inserted") or sub_info.get("date_modified") or "")[:10]
            sdec, _ = parse_clean_date(sd)
            if cdec is not None and cprec != "y":
                cand_seen.add(gid)
                cands.append((gid, b23_seqs[gid], cd, cdec, sd or cd, sdec or cdec))
        rng_b.shuffle(cands)
        for sid, seq, cs, cd, ss, sd in cands[:32]:
            seen_cur.add(sid)
            curated_items.append((sid, seq, cs, cd, ss, sd, lbl))
            if len([x for x in curated_items if x[6] in SC2_TIER1_16]) >= 500:
                break
    # Top up Tier-1 to exactly 500
    for sid, m in sorted(all_meta.items()):
        if len([x for x in curated_items if x[6] in SC2_TIER1_16]) >= 500:
            break
        if sid not in seen_cur and m["pango"] in SC2_TIER1_16 and m["date_dec"] is not None:
            seen_cur.add(sid)
            curated_items.append((sid, all_seqs[sid], m["date_str"], m["date_dec"], m["sub_date_str"] or m["date_str"], m["sub_date_dec"] or m["date_dec"], m["pango"]))

    for lbl in SC2_TIER2_LABELS:
        cands = []
        cand_seen = set()
        for sid, m in all_meta.items():
            if sid in seen_cur or sid in cand_seen or m["date_dec"] is None or m["date_prec"] == "y":
                continue
            p = m["pango"]
            ql = m.get("query_lbl", "")
            if p == lbl or (lbl == "XFN" and ql == "XFN"):
                cand_seen.add(sid)
                cands.append((sid, all_seqs[sid], m["date_str"], m["date_dec"], m["sub_date_str"] or m["date_str"], m["sub_date_dec"] or m["date_dec"]))
        for r in b23_by_lin.get(lbl, []):
            gid = r["genome_id"]
            if gid in seen_cur or gid in cand_seen:
                continue
            cd = str(r.get("collection_date") or "")
            cdec, cprec = parse_clean_date(cd)
            sub_info = sc2_data["b23_sub"].get(gid, {})
            sd = str(sub_info.get("date_inserted") or sub_info.get("date_modified") or "")[:10]
            sdec, _ = parse_clean_date(sd)
            if cdec is not None and cprec != "y":
                cand_seen.add(gid)
                cands.append((gid, b23_seqs[gid], cd, cdec, sd or cd, sdec or cdec))
        rng_b.shuffle(cands)
        for sid, seq, cs, cd, ss, sd in cands[:35]:
            seen_cur.add(sid)
            curated_items.append((sid, seq, cs, cd, ss, sd, lbl))
    curated_t2 = [x for x in curated_items if x[6] in SC2_TIER2_LABELS]
    if len(curated_t2) < 500:
        for r in sc2_data["b23_rows"]:
            if len([x for x in curated_items if x[6] in SC2_TIER2_LABELS]) >= 500:
                break
            gid = r.get("genome_id")
            lin = (r.get("lineage") or "").strip()
            if gid in b23_seqs and gid not in seen_cur and any(lin.startswith(t) for t in ("XBB", "XBF", "XDK", "XDV", "XF", "XD", "XEC")):
                cd = str(r.get("collection_date") or "")
                cdec, _ = parse_clean_date(cd)
                sub_info = sc2_data["b23_sub"].get(gid, {})
                sd = str(sub_info.get("date_inserted") or "")[:10]
                sdec, _ = parse_clean_date(sd)
                if cdec:
                    root_lbl = "XBB" if lin.startswith("XBB.") and lin not in ("XBB.1.5","XBB.1.16") else lin.split(".")[0]
                    if root_lbl in SC2_TIER2_LABELS:
                        seen_cur.add(gid)
                        curated_items.append((gid, b23_seqs[gid], cd, cdec, sd or cd, sdec or cdec, root_lbl))
    curated_t1 = [x for x in curated_items if x[6] in SC2_TIER1_16][:500]
    curated_t2 = [x for x in curated_items if x[6] in SC2_TIER2_LABELS][:500]
    curated_1000 = curated_t1 + curated_t2
    assert len({x[0] for x in curated_1000}) == len(curated_1000), f"Duplicate IDs in curated_1000: {len(curated_1000)} vs {len({x[0] for x in curated_1000})}"
    (CACHE_DIR / "bench_sc2_curated_1000_labels.json").write_text(json.dumps({x[0]: x[6] for x in curated_1000}))
    eval_sc2_cohort("sc2_curated_1000", curated_1000)

    # 1b. Strict OOD holdout masking XDV, XFZ, and XFC from certified_profiles.tsv
    sc2_ood_set = {"XDV", "XFZ", "XFC"}
    ood_sc2_items = [x for x in curated_t2 if x[6] in sc2_ood_set]
    if ood_sc2_items:
        assert len({x[0] for x in ood_sc2_items}) == len(ood_sc2_items), f"Duplicate IDs in ood_sc2_items: {len(ood_sc2_items)}"
        (CACHE_DIR / "bench_sc2_ood_xdv_xfz_xfc_labels.json").write_text(json.dumps({x[0]: x[6] for x in ood_sc2_items}))
        sc2_ood_pack = CACHE_DIR / "tmp_sc2_ood_pack"
        sc2_ood_pack.mkdir(parents=True, exist_ok=True)
        for fn in ("manifest.json", "profiles.fasta", "compendium.json"):
            shutil.copy(SC2_V3_DIR / fn, sc2_ood_pack / fn)
        kept_cert = [
            l for l in (SC2_V3_DIR / "certified_profiles.tsv").read_text().splitlines()
            if l.startswith("#") or l.split("\t")[0] not in sc2_ood_set
        ]
        (sc2_ood_pack / "certified_profiles.tsv").write_text("\n".join(kept_cert) + "\n")
        eval_sc2_cohort("sc2_ood_xdv_xfz_xfc", ood_sc2_items, pack_override=sc2_ood_pack, ood_masked_set=sc2_ood_set)

    # 2. 10,062 World Union
    world_items = []
    for sid, seq in sc2_data["world_seqs"].items():
        m = sc2_data["world_meta"][sid]
        world_items.append((sid, seq, m["date_str"], m["date_dec"], m["sub_date_str"] or m["date_str"], m["sub_date_dec"] or m["date_dec"], None))
    eval_sc2_cohort("sc2_world_10062", world_items)

    # 3. 23,404 Cached BV-BRC
    b23_items = []
    for r in sc2_data["b23_rows"]:
        gid = r.get("genome_id")
        if gid not in b23_seqs:
            continue
        cd = str(r.get("collection_date") or "").strip() or str(r.get("collection_year") or "")
        cdec, _ = parse_clean_date(cd)
        sub_info = sc2_data["b23_sub"].get(gid, {})
        sd = str(sub_info.get("date_inserted") or sub_info.get("date_modified") or "")[:10]
        sdec, _ = parse_clean_date(sd)
        b23_items.append((gid, b23_seqs[gid], cd, cdec, sd or cd, sdec or cdec, None))
    eval_sc2_cohort("sc2_bvbrc_23404", b23_items)

    # 4. 100,000 BV-BRC Stream (with both collection_date AND date_inserted)
    s100_dir = CACHE_DIR / "sc2_100k_stream"
    s100_dir.mkdir(parents=True, exist_ok=True)
    s100_meta_path = s100_dir / "queue_100k.json"
    if s100_meta_path.exists() and all((s100_dir / f"chunk_{i:02d}.fasta.gz").exists() for i in range(10)):
        log("[SC2 100K] Reusing cached 100,000-genome BV-BRC chunks...")
        queue_100k = json.loads(s100_meta_path.read_text())
    else:
        log("[SC2 100K] Building deterministic 100,000-genome queue + date_inserted from BV-BRC...")
        base = "eq(taxon_lineage_ids,2697049)&eq(genome_status,Complete)&gt(genome_length,29000)&eq(contigs,1)"
        fields = "genome_id,genbank_accessions,lineage,isolation_country,collection_date,collection_year,date_inserted,date_modified"
        seen_gids = set()
        queue_100k = []
        def _add(rows, label):
            for rr in rows:
                gid = rr.get("genome_id")
                if gid and gid not in seen_gids:
                    seen_gids.add(gid)
                    rr["stratum"] = label
                    queue_100k.append(rr)
        _add(bvbrc_fetch_pages(f"{base}&eq(collection_year,2026)", fields), "2026_all")
        _add(bvbrc_fetch_pages(f"{base}&eq(collection_year,2025)", fields), "2025_all")
        rem = max(0, 100000 - len(queue_100k))
        if rem > 0:
            _add(bvbrc_fetch_pages(f"{base}&eq(collection_year,2024)&ne(lineage,*)", fields, max_rows=min(5000, rem), page_size=6000), "2024_missing")
        rem2 = max(0, 100000 - len(queue_100k))
        slice_yrs = [2024, 2023, 2022, 2021]
        per_yr = rem2 // 4
        extra = rem2 % 4
        for idx_y, yr in enumerate(slice_yrs):
            need = per_yr + (1 if idx_y < extra else 0)
            if need > 0:
                # Over-fetch slightly to ensure uniqueness after seen_gids
                got = bvbrc_fetch_pages(f"{base}&eq(collection_year,{yr})", fields, max_rows=need + 5500, page_size=10000)
                added = 0
                for rr in got:
                    gid = rr.get("genome_id")
                    if gid and gid not in seen_gids:
                        seen_gids.add(gid)
                        rr["stratum"] = f"{yr}_slice"
                        queue_100k.append(rr)
                        added += 1
                        if added >= need:
                            break
        queue_100k = queue_100k[:100000]
        s100_meta_path.write_text(json.dumps(queue_100k))
        # Fetch sequences into 10 chunks of 10,000
        for ci in range(10):
            chunk_rows = queue_100k[ci * 10000:(ci + 1) * 10000]
            log(f"  [SC2 100K] Fetching FASTAs for chunk {ci+1}/10 ({len(chunk_rows)} genomes)...")
            seqs_c = fetch_bvbrc_seqs_by_gids(chunk_rows, batch_size=250, workers=8)
            with gzip.open(s100_dir / f"chunk_{ci:02d}.fasta.gz", "wt") as f:
                for rr in chunk_rows:
                    gid = rr["genome_id"]
                    if gid in seqs_c:
                        f.write(f">{gid}\n{seqs_c[gid]}\n")

    # Screen all 100K genomes in standard vs --temporal-guard(collection) vs --temporal-guard(date_inserted)
    meta_100k = {r["genome_id"]: r for r in queue_100k}
    tot_100k = {m: {"n": 0, "dt": 0.0, "chg": 0, "viol_c": 0, "viol_d": 0, "pre_total": 0, "pre_call": 0, "pre_donor": 0, "masked": 0, "cats": collections.Counter(), "prevented_lineages": collections.Counter()} for m in ("standard", "temporal_guard_collection", "temporal_guard_submission")}
    for ci in range(10):
        cseqs = read_fasta(s100_dir / f"chunk_{ci:02d}.fasta.gz")
        fa_c = s100_dir / "tmp_col.fasta"
        fa_s = s100_dir / "tmp_sub.fasta"
        c_dates, s_dates = {}, {}
        with open(fa_c, "w") as fc, open(fa_s, "w") as fs:
            for gid, seq in cseqs.items():
                r = meta_100k[gid]
                cd = str(r.get("collection_date") or "").strip() or str(r.get("collection_year") or "")
                sd = str(r.get("date_inserted") or r.get("date_modified") or "")[:10]
                cdec, _ = parse_clean_date(cd)
                sdec, _ = parse_clean_date(sd)
                idc = f"{gid}|{cd}" if cd else gid
                ids = f"{gid}|{sd}" if sd else idc
                fc.write(f">{idc}\n{seq}\n")
                fs.write(f">{ids}\n{seq}\n")
                c_dates[idc] = cdec
                s_dates[ids] = sdec if sdec is not None else cdec
        r_std, dt_0 = run_eggman_cli(SC2_V3_DIR, fa_c, gamma=3.5, min_snps=3, p_thresh=0.01, temporal_guard=False)
        r_col, dt_1 = run_eggman_cli(SC2_V3_DIR, fa_c, gamma=3.5, min_snps=3, p_thresh=0.01, temporal_guard=True, temporal_mode="sampled")
        r_sub, dt_2 = run_eggman_cli(SC2_V3_DIR, fa_s, gamma=3.5, min_snps=3, p_thresh=0.01, temporal_guard=True, temporal_mode="submission")
        fa_c.unlink(missing_ok=True)
        fa_s.unlink(missing_ok=True)
        std_sigs = [(a["called_lineage"], a["structural_formula"]) for a in r_std]
        for mlbl, rlist, dmap, dts in [
            ("standard", r_std, c_dates, dt_0),
            ("temporal_guard_collection", r_col, c_dates, dt_1),
            ("temporal_guard_submission", r_sub, s_dates, dt_2),
        ]:
            st = tot_100k[mlbl]
            st["n"] += len(rlist)
            st["dt"] += dts
            an = count_anachronistic(rlist, dmap, mlbl)
            st["viol_c"] += an["floor_viol_call"]
            st["viol_d"] += an["floor_viol_donor"]
            st["pre_total"] += an["pre2023_5_total"]
            st["pre_call"] += an["pre2023_5_future_called"]
            st["pre_donor"] += an["pre2023_5_future_donor_seg"]
            st["masked"] += sum(1 for r in rlist if r.get("masked_future_trunks"))
            for idx_r, r in enumerate(rlist):
                st["cats"][r["called_category"]] += 1
                if mlbl != "standard" and (r["called_lineage"], r["structural_formula"]) != std_sigs[idx_r]:
                    st["chg"] += 1
                    st["prevented_lineages"][f"{r_std[idx_r]['called_lineage']}->{r['called_lineage']}"] += 1

    for mlbl, st in tot_100k.items():
        ablation_rows.append({
            "pathogen": "SARS-CoV-2",
            "cohort": "sc2_bvbrc_100k",
            "mode": mlbl,
            "n_sequences": st["n"],
            "runtime_sec": round(st["dt"], 3),
            "n_masked_seqs": st["masked"],
            "n_changed_vs_standard": st["chg"],
            "floor_viol_call": st["viol_c"],
            "floor_viol_donor": st["viol_d"],
            "pre2023_5_total": st["pre_total"],
            "pre_cutoff_future_called": st["pre_call"],
            "pre_cutoff_future_donor": st["pre_donor"],
            "exact_tier2_acc_pct": "",
            "tier1_accuracy_pct": "",
            "urf_count": st["cats"].get("URF", 0),
            "certified_count": st["cats"].get("Certified_CRF", 0),
            "pure_count": st["cats"].get("Pure_Subtype", 0),
            "notes": ",".join(f"{k}:{v}" for k, v in st["prevented_lineages"].most_common(6)) if st["prevented_lineages"] else "",
        })
        log(f"  [SC2 100K] {mlbl:26s}: N={st['n']:,} chg={st['chg']} viol_c={st['viol_c']} viol_d={st['viol_d']} | {st['dt']:.2f}s")
    details["sc2_100k_prevented_anachronisms"] = dict(tot_100k["temporal_guard_collection"]["prevented_lineages"])
    return ablation_rows, details


def run_iav_ablation_benchmarks(iav_data: Dict[str, Any]) -> Tuple[List[Dict[str, Any]], Dict[str, Any]]:
    log("[IAV Benchmarks] Running standard vs --temporal-guard ablations on 1,200 SRA + 650 OOD + 3,980 BV-BRC...")
    ablation_rows: List[Dict[str, Any]] = []
    sra_seqs = iav_data["sra_seqs"]
    sra_info = iav_data["sra_info"]
    bv_seqs = iav_data["bvbrc_seqs"]
    bv_meta = iav_data["bvbrc_meta"]

    iav_man = json.loads((IAV_V2_DIR / "manifest.json").read_text())
    iav_comp = json.loads((IAV_V2_DIR / "compendium.json").read_text())
    iav_slack = float(iav_man.get("default_temporal_slack_yr") or 0.5)
    iav_trunk_sampled: Dict[str, float] = {}
    iav_trunk_submitted: Dict[str, float] = {}
    for t in iav_man.get("basal_trunks", []):
        fl_s = t.get("first_sampled_year") if t.get("first_sampled_year") is not None else t.get("t_mrca")
        fl_sub = t.get("first_submitted_year") if t.get("first_submitted_year") is not None else fl_s
        if fl_s is not None:
            iav_trunk_sampled[t["name"]] = float(fl_s)
        if fl_sub is not None:
            iav_trunk_submitted[t["name"]] = float(fl_sub)
    iav_clade_sampled: Dict[str, float] = {}
    iav_clade_submitted: Dict[str, float] = {}
    for e in iav_comp:
        fl_s = e.get("first_sampled_year") if e.get("first_sampled_year") is not None else e.get("t_mrca")
        fl_sub = e.get("first_submitted_year") if e.get("first_submitted_year") is not None else fl_s
        if fl_s is not None:
            iav_clade_sampled[e["label"]] = float(fl_s)
        if fl_sub is not None:
            iav_clade_submitted[e["label"]] = float(fl_sub)
    for t in iav_man.get("basal_trunks", []):
        if t.get("clonal_certified_label"):
            lbl = t["clonal_certified_label"]
            if t["name"] in iav_trunk_sampled and lbl not in iav_clade_sampled:
                iav_clade_sampled[lbl] = iav_trunk_sampled[t["name"]]
            if t["name"] in iav_trunk_submitted and lbl not in iav_clade_submitted:
                iav_clade_submitted[lbl] = iav_trunk_submitted[t["name"]]

    nam_2023_2024_donors = {"NAm_B313", "NAm_D1", "NAm_C2"}
    nam_2023_2024_clades = {"B3.13", "D1.1", "D1.3", "D1.8", "B3.6", "C2.1"}

    def count_iav_anach(rows: List[Dict[str, Any]], d_map: Dict[str, Optional[float]], mode_lbl: str = "standard") -> Dict[str, int]:
        use_sub = (mode_lbl == "temporal_guard_submission")
        t_floors = iav_trunk_submitted if use_sub else iav_trunk_sampled
        c_floors = iav_clade_submitted if use_sub else iav_clade_sampled
        pre_tot = 0
        pre_call = 0
        pre_don = 0
        pre22_tot = 0
        pre22_nam_any = 0
        floor_viol_call = 0
        floor_viol_donor = 0
        nam_all = {"NAm_B313", "NAm_D1", "NAm_B36", "NAm_C2", "NAm_B2", "NAm_B32"}
        for r in rows:
            dq = d_map.get(r["sequence_id"])
            if dq is None:
                continue
            cl = r["called_lineage"]
            trks = set(r["constituent_trunks"])
            cfl = c_floors.get(cl)
            if cfl is not None and cfl > dq + iav_slack:
                floor_viol_call += 1
            if any(t_floors.get(tk, -1e9) > dq + iav_slack for tk in trks):
                floor_viol_donor += 1
            if dq < 2022.0:
                pre22_tot += 1
                if trks & nam_all:
                    pre22_nam_any += 1
            if dq < 2023.0:
                pre_tot += 1
                if cl in nam_2023_2024_clades:
                    pre_call += 1
                if trks & nam_2023_2024_donors:
                    pre_don += 1
        return {
            "pre2023_0_total": pre_tot,
            "pre_call": pre_call,
            "pre_donor": pre_don,
            "pre2022_0_total": pre22_tot,
            "pre2022_0_nam_any": pre22_nam_any,
            "floor_viol_call": floor_viol_call,
            "floor_viol_donor": floor_viol_donor,
        }

    # 1. Build 1,200 SRA cohort (150 A1/A3 + 950 8 reassortant genotypes + 100 USDA unassigned) & 650 OOD holdout (500 B3.13 + 150 D1.1)
    by_gt = collections.defaultdict(list)
    for r, inf in sra_info.items():
        if r in sra_seqs:
            by_gt[inf["genotype"]].append(r)

    rng = np.random.default_rng(2026)
    cohort_1200: List[Tuple[str, Optional[str]]] = []
    # 75 A1 + 75 A3
    for gt, n_q in [("A1", 75), ("A3", 75), ("B3.13", 300), ("D1.1", 150), ("B3.6", 120), ("D1.3", 80), ("D1.8", 80), ("B2.1", 80), ("B3.2", 40), ("C2.1", 100), ("Not_assigned", 100)]:
        avail = sorted(by_gt[gt])
        pick = avail[:n_q] if len(avail) <= n_q else [avail[i] for i in sorted(rng.choice(len(avail), size=n_q, replace=False))]
        exp = None if gt == "Not_assigned" else gt
        for r in pick:
            cohort_1200.append((r, exp))

    used_1200 = {r for r, _ in cohort_1200}
    ood_b313 = [r for r in by_gt["B3.13"] if r not in used_1200][:500]
    if len(ood_b313) < 500:
        ood_b313 = by_gt["B3.13"][:500]
    ood_d11 = [r for r in by_gt["D1.1"] if r not in used_1200][:150]
    if len(ood_d11) < 150:
        ood_d11 = by_gt["D1.1"][:150]
    cohort_650_ood = [(r, "B3.13") for r in ood_b313] + [(r, "D1.1") for r in ood_d11]

    # Scratch OOD pack masking B3.13 and D1.1 out of certified_profiles.tsv
    ood_pack_dir = CACHE_DIR / "tmp_iav_ood_pack"
    ood_pack_dir.mkdir(parents=True, exist_ok=True)
    for fn in ("manifest.json", "profiles.fasta", "compendium.json"):
        shutil.copy(IAV_V2_DIR / fn, ood_pack_dir / fn)
    cert_lines = [l for l in (IAV_V2_DIR / "certified_profiles.tsv").read_text().splitlines() if not (l.startswith("B3.13\t") or l.startswith("D1.1\t"))]
    (ood_pack_dir / "certified_profiles.tsv").write_text("\n".join(cert_lines) + "\n")

    def run_iav_eval(cname: str, pack_d: Path, items: List[Tuple[str, str, Optional[str], Optional[float], Optional[str], Optional[float], Optional[str]]], is_ood: bool = False):
        fa_c = CACHE_DIR / f"bench_{cname}_col.fasta"
        fa_s = CACHE_DIR / f"bench_{cname}_sub.fasta"
        cd_map, sd_map, exp_map = {}, {}, {}
        with open(fa_c, "w") as fc, open(fa_s, "w") as fs:
            for sid, seq, cs, cd, ss, sd, exp in items:
                idc = f"{sid}|{cs}" if cs else sid
                ids = f"{sid}|{ss}" if ss else sid
                fc.write(f">{idc}\n{seq}\n")
                fs.write(f">{ids}\n{seq}\n")
                cd_map[idc] = cd
                sd_map[ids] = sd
                if exp:
                    exp_map[sid] = exp
        r_std, dt0 = run_eggman_cli(pack_d, fa_c, temporal_guard=False)
        r_col, dt1 = run_eggman_cli(pack_d, fa_c, temporal_guard=True, temporal_mode="sampled")
        r_sub, dt2 = run_eggman_cli(pack_d, fa_s, temporal_guard=True, temporal_mode="submission")

        ood_expected_formula = {
            "B3.13": "NAm_B313[EA1:PA+HA+NA+MP]",
            "D1.1": "NAm_D1[EA3:PB1+PA+HA+MP+NS]", # check exact Tier-1 formula below!
        }
        std_sigs = [(a["called_lineage"], a["structural_formula"]) for a in r_std]
        for mlbl, rlist, dmap, dt_s in [
            ("standard", r_std, cd_map, dt0),
            ("temporal_guard_collection", r_col, cd_map, dt1),
            ("temporal_guard_submission", r_sub, sd_map, dt2),
        ]:
            cats = collections.Counter(r["called_category"] for r in rlist)
            an = count_iav_anach(rlist, dmap, mlbl)
            n_masked = sum(1 for r in rlist if r.get("masked_future_trunks"))
            n_changed = 0 if mlbl == "standard" else sum(
                1 for idx_r, r in enumerate(rlist) if (r["called_lineage"], r["structural_formula"]) != std_sigs[idx_r]
            )
            t1_acc, t2_acc = "", ""
            if exp_map:
                if is_ood:
                    ok_urf = sum(1 for r in rlist if r["called_lineage"] == "URF")
                    ok_form = sum(
                        1 for r in rlist
                        if (exp_map.get(r["sequence_id"].split("|")[0]) == "B3.13" and r["structural_formula"] == "NAm_B313[EA1:PA+HA+NA+MP]")
                        or (exp_map.get(r["sequence_id"].split("|")[0]) == "D1.1" and "NAm_D1" in r["constituent_trunks"] and "EA3" in r["constituent_trunks"])
                    )
                    t2_acc = round(100.0 * ok_urf / len(rlist), 2)
                    t1_acc = round(100.0 * ok_form / len(rlist), 2)
                else:
                    labeled = [r for r in rlist if r["sequence_id"].split("|")[0] in exp_map]
                    ok2 = sum(1 for r in labeled if r["called_lineage"] == exp_map[r["sequence_id"].split("|")[0]])
                    t2_acc = round(100.0 * ok2 / len(labeled), 2) if labeled else ""
            ablation_rows.append({
                "pathogen": "Influenza_A_H5Nx",
                "cohort": cname,
                "mode": mlbl,
                "n_sequences": len(rlist),
                "runtime_sec": round(dt_s, 3),
                "n_masked_seqs": n_masked,
                "n_changed_vs_standard": n_changed,
                "floor_viol_call": an["floor_viol_call"],
                "floor_viol_donor": an["floor_viol_donor"],
                "pre2023_5_total": an["pre2023_0_total"],
                "pre_cutoff_future_called": an["pre_call"],
                "pre_cutoff_future_donor": an["pre_donor"],
                "exact_tier2_acc_pct": t2_acc,
                "tier1_accuracy_pct": t1_acc,
                "urf_count": cats.get("URF", 0),
                "certified_count": cats.get("Certified_Reassortant", 0),
                "pure_count": cats.get("Pure_Subtype", 0),
                "notes": f"pre2023_NAm_23_24={an['pre_donor']}/{an['pre2023_0_total']};pre2022_NAm_any={an['pre2022_0_nam_any']}/{an['pre2022_0_total']}",
            })
            log(f"  [IAV Bench] {cname:18s} | {mlbl:26s} | N={len(rlist):4d} | chg={n_changed:3d} viol_c={an['floor_viol_call']:2d} viol_d={an['floor_viol_donor']:3d} | t1={t1_acc} t2={t2_acc} | {dt_s:.2f}s")

    items_1200 = [
        (r, sra_seqs[r], sra_info[r]["col_date_str"], sra_info[r]["col_date_dec"], sra_info[r]["rel_date_str"], sra_info[r]["rel_date_dec"], exp)
        for (r, exp) in cohort_1200
    ]
    run_iav_eval("iav_sra_1200", IAV_V2_DIR, items_1200)

    items_650 = [
        (r, sra_seqs[r], sra_info[r]["col_date_str"], sra_info[r]["col_date_dec"], sra_info[r]["rel_date_str"], sra_info[r]["rel_date_dec"], exp)
        for (r, exp) in cohort_650_ood
    ]
    run_iav_eval("iav_sra_650_ood", ood_pack_dir, items_650, is_ood=True)

    items_3980 = [
        (iid, bv_seqs[iid], bv_meta[iid]["col_date_str"], bv_meta[iid]["col_date_dec"], bv_meta[iid]["sub_date_str"], bv_meta[iid]["sub_date_dec"], None)
        for iid in sorted(bv_seqs)
    ]
    run_iav_eval("iav_bvbrc_3980", IAV_V2_DIR, items_3980)

    return ablation_rows, {"n_1200": len(items_1200), "n_650_ood": len(items_650), "n_3980": len(items_3980)}


def main() -> None:
    global CACHE_DIR
    parser = argparse.ArgumentParser()
    parser.add_argument("--step", default="all", choices=["all", "sc2", "iav", "bench"])
    parser.add_argument("--reuse-ml-clock", action="store_true", help="Reuse cached LSD2/TreeTime results from TSVs if present")
    parser.add_argument("--work-dir", default=str(CACHE_DIR), help="Cache/working directory for downloaded/intermediate FASTAs")
    args = parser.parse_args()
    CACHE_DIR = Path(args.work_dir)
    CACHE_DIR.mkdir(parents=True, exist_ok=True)

    t_start = time.time()
    summary: Dict[str, Any] = {}
    if OUT_SUMMARY_JSON.exists() and args.step != "all":
        try:
            summary = json.loads(OUT_SUMMARY_JSON.read_text())
        except Exception:
            pass

    sc2_data = None
    iav_data = None

    if args.step in ("all", "sc2"):
        sc2_data = prepare_sc2_datasets()
        summary["sars_cov2"] = run_sc2_basis_and_clock(sc2_data, reuse_ml_clock=args.reuse_ml_clock)
        OUT_SUMMARY_JSON.write_text(json.dumps(summary, indent=2) + "\n")

    if args.step in ("all", "iav"):
        iav_data = prepare_iav_datasets()
        iav_res = run_iav_basis_clock_and_sarbeco(iav_data, reuse_ml_clock=args.reuse_ml_clock)
        # Keep only JSON-serializable parts in summary
        summary["influenza_a_h5nx"] = {k: v for k, v in iav_res.items() if k not in ("iav_gt_pools_ids", "all_iav_seqs", "bv_scan_map")}
        summary["sarbecovirus"] = {"first_sampled_years": iav_res["sarbecovirus_first_sampled_years"]}
        OUT_SUMMARY_JSON.write_text(json.dumps(summary, indent=2) + "\n")

    if args.step in ("all", "bench"):
        if sc2_data is None:
            sc2_data = prepare_sc2_datasets()
        if iav_data is None:
            iav_data = prepare_iav_datasets()
        sc2_ab_rows, sc2_ab_meta = run_sc2_ablation_benchmarks(sc2_data, summary.get("sars_cov2", {}))
        iav_ab_rows, iav_ab_meta = run_iav_ablation_benchmarks(iav_data)
        all_ab_rows = sc2_ab_rows + iav_ab_rows
        with open(OUT_ABLATION_TSV, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(all_ab_rows[0].keys()), delimiter="\t")
            w.writeheader()
            w.writerows(all_ab_rows)
        summary["temporal_guard_ablation"] = {
            "sc2": sc2_ab_meta,
            "iav": iav_ab_meta,
            "rows": all_ab_rows,
        }
        summary["wall_time_sec"] = round(time.time() - t_start, 2)
        OUT_SUMMARY_JSON.write_text(json.dumps(summary, indent=2) + "\n")
        log(f"[+] Completed unified ChronAeon + EGGMAN pipeline in {summary['wall_time_sec']}s!")


if __name__ == "__main__":
    main()
