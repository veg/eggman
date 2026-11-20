#!/usr/bin/env python3
"""
24_bvbrc_100k_streaming_census_and_screen.py
=============================================
Live BV-BRC archive census and large-scale in-memory EGGMAN streaming screen:
  Module A: Exact live BV-BRC archive census for SARS-CoV-2 & Influenza A
            -> results/bvbrc_archive_census.json
  Module B: In-memory 100,000-genome complete SARS-CoV-2 streaming screen
            -> results/bvbrc_100k_sc2_stream_summary.json
            -> results/bvbrc_100k_sc2_novel_and_mislabeled.tsv
            -> data/bvbrc_stream_hits/sc2_novel_hits.fasta.gz
  Module C: Exhaustive complete 8-segment 2020-2026 H5Nx (h_type=5) IAV stream
            -> results/bvbrc_iav_h5_8seg_stream_summary.json
            -> results/bvbrc_iav_h5_novel_constellations.tsv
            -> data/bvbrc_stream_hits/iav_h5_novel_hits.fasta.gz
"""

import argparse
import concurrent.futures
import csv
import gzip
import json
import os
import random
import re
import subprocess
import tempfile
import time
import urllib.error
import urllib.request
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path("/Users/sergei/Projects/HIV-recombination")
RESULTS_DIR = ROOT / "results"
HITS_DIR = ROOT / "data" / "bvbrc_stream_hits"

EGGMAN_BIN = "/Users/sergei/Projects/HIV-recombination/rhiztype/target/release/eggman"
SC2_PACK = "/Users/sergei/Projects/HIV-recombination/rhiztype/packs/sars_cov2_wuhan"
IAV_PACK = "/Users/sergei/Projects/HIV-recombination/rhiztype/packs/influenza_a_8seg"
DEFAULT_CACHE_DIR = "/Users/sergei/.gemini/antigravity-cli/brain/d7d3d768-93f4-475c-81ac-ccb936c7805d/scratch/chronaeon_cache"


def load_diverged_post_wuhan_trunks(pack_dir, min_coding_snps=45):
    """Identify diverged post-Wuhan trunks purely from pack coding SNP distance to Wuhan B."""
    prof_fa = Path(pack_dir) / "profiles.fasta"
    seqs = {}
    cur = None
    buf = []
    with open(prof_fa, "r", encoding="utf-8") as f:
        for ln in f:
            ln = ln.strip()
            if not ln:
                continue
            if ln.startswith(">"):
                if cur:
                    seqs[cur] = "".join(buf).upper()
                cur = ln[1:].split()[0]
                buf = []
            else:
                buf.append(ln)
        if cur:
            seqs[cur] = "".join(buf).upper()
    b_seq = seqs.get("B", "")
    if not b_seq:
        raise RuntimeError(f"Missing ancestral B profile in {prof_fa}")
    b_cod = b_seq[265:29674]
    diverged = set()
    for name, seq in seqs.items():
        if name == "B":
            continue
        s_cod = seq[265:29674]
        n_diff = sum(
            1 for x, y in zip(b_cod, s_cod)
            if x in "ACGT" and y in "ACGT" and x != y
        )
        if n_diff >= min_coding_snps:
            diverged.add(name)
    return diverged


DIVERGED_POST_WUHAN_TRUNKS = load_diverged_post_wuhan_trunks(SC2_PACK)

OUT_CENSUS = RESULTS_DIR / "bvbrc_archive_census.json"
OUT_SC2_JSON = RESULTS_DIR / "bvbrc_100k_sc2_stream_summary.json"
OUT_SC2_TSV = RESULTS_DIR / "bvbrc_100k_sc2_novel_and_mislabeled.tsv"
OUT_SC2_FASTA_GZ = HITS_DIR / "sc2_novel_hits.fasta.gz"

OUT_IAV_JSON = RESULTS_DIR / "bvbrc_iav_h5_8seg_stream_summary.json"
OUT_IAV_TSV = RESULTS_DIR / "bvbrc_iav_h5_novel_constellations.tsv"
OUT_IAV_FASTA_GZ = HITS_DIR / "iav_h5_novel_hits.fasta.gz"

BVBRC_GENOME_URL = "https://www.bv-brc.org/api/genome/"
BVBRC_SEQ_URL = "https://www.bv-brc.org/api/genome_sequence/"
HEADERS_JSON = {"Accept": "application/json", "User-Agent": "curl/8.7.1"}
HEADERS_FASTA_POST = {
    "Content-Type": "application/rqlquery+x-www-form-urlencoded",
    "Accept": "application/dna+fasta",
    "User-Agent": "curl/8.7.1",
}

GID_RE = re.compile(r"\|\s*([0-9]+\.[0-9]+)\s*\]\s*$")
ACCN_RE = re.compile(r"^>accn\|([^\s|]+)")

XDV_ALIAS_ROOTS = {"NB", "PQ", "RC", "RH", "RK", "SV", "SZ", "TM"}
XBB_ALIAS_ROOTS = {"EG", "FL", "FE", "FU", "FY", "GJ", "GN", "GW", "HA", "HK", "HV", "JD", "JG"}

IAV_SEG_NAMES = {
    "1": "PB2",
    "2": "PB1",
    "3": "PA",
    "4": "HA",
    "5": "NP",
    "6": "NA",
    "7": "MP",
    "8": "NS",
}


class ThrottledError(RuntimeError):
    """Raised when BV-BRC returns HTTP 429 or 503."""


def http_get_json_with_headers(url, timeout=45, retries=5):
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
            wait = 1.5 * (2 ** attempt) + random.uniform(0.1, 0.6)
            if e.code in (429, 503):
                wait *= 2.0
            time.sleep(wait)
        except Exception as e:
            last_err = e
            time.sleep(1.5 * (2 ** attempt) + random.uniform(0.1, 0.5))
    raise RuntimeError(f"GET failed after {retries} attempts for {url}: {last_err}")


def bvbrc_count(rql):
    url = f"{BVBRC_GENOME_URL}?{rql}&limit(1)"
    cr, _ = http_get_json_with_headers(url, timeout=30, retries=5)
    if "/" not in cr:
        raise RuntimeError(f"Missing Content-Range total for {url}: {cr!r}")
    return int(cr.split("/")[-1])


def bvbrc_fetch_pages(base_rql, select_fields, max_rows=None, page_size=10000, seen_global=None):
    """
    Fetch BV-BRC genome metadata cleanly using keyset pagination on +genome_id.
    Handles BV-BRC's inclusive gt(genome_id, ...) overlap by tracking seen_gids.
    """
    if seen_global is None:
        seen_global = set()
    collected = []
    seen_local = set()
    last_gid = None
    offset_fallback = 0
    use_offset = False

    while True:
        if max_rows is not None and len(collected) >= max_rows:
            break
        cur_limit = page_size
        if max_rows is not None:
            # Over-fetch modestly when skipping previously seen IDs, capped at page_size
            cur_limit = min(page_size, max(500, (max_rows - len(collected)) * 2))

        if use_offset:
            url = (
                f"{BVBRC_GENOME_URL}?{base_rql}&sort(+genome_id)"
                f"&limit({cur_limit},{offset_fallback})&select({select_fields})"
            )
        else:
            gt_clause = f"&gt(genome_id,{last_gid})" if last_gid else ""
            url = (
                f"{BVBRC_GENOME_URL}?{base_rql}{gt_clause}&sort(+genome_id)"
                f"&limit({cur_limit})&select({select_fields})"
            )

        _, rows = http_get_json_with_headers(url, timeout=60, retries=5)
        if not rows:
            break

        new_in_page = 0
        new_unique_global = 0
        for r in rows:
            gid = r.get("genome_id")
            if not gid:
                continue
            if gid not in seen_local:
                seen_local.add(gid)
                new_in_page += 1
            if gid not in seen_global:
                seen_global.add(gid)
                collected.append(r)
                new_unique_global += 1
                if max_rows is not None and len(collected) >= max_rows:
                    break

        if max_rows is not None and len(collected) >= max_rows:
            break

        if len(rows) < cur_limit:
            break

        if new_in_page == 0:
            # Keyset stalled cleanly; switch to offset or terminate
            if not use_offset:
                use_offset = True
                offset_fallback = len(seen_local)
            else:
                break
        else:
            if use_offset:
                offset_fallback += len(rows)
            else:
                last_gid = rows[-1].get("genome_id")
                if not last_gid:
                    break

    return collected


def post_fasta_batch(genome_ids, retries=5):
    body = f"in(genome_id,({','.join(genome_ids)}))&limit(1000)".encode("utf-8")
    last_err = None
    saw_throttle = False
    for attempt in range(retries):
        try:
            req = urllib.request.Request(BVBRC_SEQ_URL, data=body, headers=HEADERS_FASTA_POST)
            with urllib.request.urlopen(req, timeout=75) as resp:
                return resp.read().decode("utf-8", errors="replace"), saw_throttle
        except urllib.error.HTTPError as e:
            last_err = e
            if e.code in (429, 503):
                saw_throttle = True
                time.sleep(2.0 * (2 ** attempt) + random.uniform(0.2, 0.8))
            else:
                time.sleep(1.2 * (2 ** attempt) + random.uniform(0.1, 0.5))
        except Exception as e:
            last_err = e
            time.sleep(1.5 * (2 ** attempt) + random.uniform(0.1, 0.5))
    if saw_throttle:
        raise ThrottledError(f"Persistent HTTP 429/503 throttling after {retries} retries: {last_err}")
    raise RuntimeError(f"POST genome_sequence failed after {retries} retries for {len(genome_ids)} IDs: {last_err}")


def fetch_fasta_batches_adaptive(batches, workers=8):
    """
    Fetch batches of genome_ids concurrently. Starts at `workers` (default 8)
    and automatically downgrades to 4 workers if BV-BRC emits HTTP 429/503.
    Raises immediately if BV-BRC persistently throttles even at 4 workers.
    """
    results = [None] * len(batches)
    remaining = list(range(len(batches)))
    cur_workers = workers

    while remaining:
        throttled_indices = []
        batch_retries = 2 if cur_workers > 4 else 5
        with concurrent.futures.ThreadPoolExecutor(max_workers=cur_workers) as ex:
            fut_map = {ex.submit(post_fasta_batch, batches[idx], batch_retries): idx for idx in remaining}
            for fut in concurrent.futures.as_completed(fut_map):
                idx = fut_map[fut]
                try:
                    txt, saw_thr = fut.result()
                    results[idx] = txt
                    if saw_thr and cur_workers > 4:
                        cur_workers = 4
                except ThrottledError as te:
                    if cur_workers <= 4:
                        raise RuntimeError(f"BV-BRC persistently throttled even at max_workers=4: {te}") from te
                    throttled_indices.append(idx)

        if not throttled_indices:
            break
        print(f"  [Throttle] HTTP 429/503 detected; falling back to max_workers=4 for {len(throttled_indices)} batches...", flush=True)
        cur_workers = 4
        remaining = throttled_indices
        time.sleep(3.0)

    return results, cur_workers


def parse_fasta_and_join(fasta_texts, valid_gids, acc_to_gid):
    """
    Parse BV-BRC FASTA responses and join headers to metadata genome_ids.
    Primary: trailing `| <gid>]` inside `[...]`.
    Fallback: `accn|<acc>` matched against `genbank_accessions`.
    """
    seq_by_gid = {}
    total_headers = 0
    joined_headers = 0
    unmatched_examples = []

    for txt in fasta_texts:
        if not txt:
            continue
        cur_gid = None
        buf = []

        def flush_record():
            if cur_gid is not None and buf:
                seq = "".join(buf).upper()
                if cur_gid not in seq_by_gid or len(seq) > len(seq_by_gid[cur_gid]):
                    seq_by_gid[cur_gid] = seq

        for raw_line in txt.splitlines():
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                flush_record()
                total_headers += 1
                cur_gid = None
                buf = []

                m = GID_RE.search(line)
                if m and m.group(1) in valid_gids:
                    cur_gid = m.group(1)
                elif "[" in line and "|" in line:
                    bracket_part = line.rsplit("[", 1)[-1].rstrip("]").strip()
                    cand = bracket_part.rsplit("|", 1)[-1].strip()
                    if cand in valid_gids:
                        cur_gid = cand

                if cur_gid is None:
                    am = ACCN_RE.match(line)
                    if am:
                        acc = am.group(1).strip()
                        cand = acc_to_gid.get(acc) or acc_to_gid.get(acc.split(".")[0])
                        if cand in valid_gids:
                            cur_gid = cand

                if cur_gid is not None:
                    joined_headers += 1
                elif len(unmatched_examples) < 5:
                    unmatched_examples.append(line[:180])
            else:
                if cur_gid is not None:
                    buf.append(line)

        flush_record()

    return seq_by_gid, total_headers, joined_headers, unmatched_examples


def run_eggman_json(pack_dir, fasta_path, gamma=None, min_snps=None, p_thresh=None):
    fd, tmp_json = tempfile.mkstemp(suffix=".eggman.json")
    os.close(fd)
    cmd = [EGGMAN_BIN, "--pack", str(pack_dir), "--json", tmp_json, "-q"]
    if gamma is not None:
        cmd.extend(["--gamma", str(gamma)])
    if min_snps is not None:
        cmd.extend(["--min-snps", str(min_snps)])
    if p_thresh is not None:
        cmd.extend(["--p-thresh", str(p_thresh)])
    cmd.append(str(fasta_path))

    t0 = time.time()
    try:
        subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, text=True)
        dt = time.time() - t0
        with open(tmp_json, "r", encoding="utf-8") as jf:
            data = json.load(jf)
        rows = data if isinstance(data, list) else data.get("results", [])
        return rows, dt
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"EGGMAN failed ({' '.join(cmd)}): {e.stderr}") from e
    finally:
        if os.path.exists(tmp_json):
            os.remove(tmp_json)


def min_bp_fisher_p(res_row):
    bps = res_row.get("breakpoints") or []
    pvals = [float(bp["p_fisher"]) for bp in bps if bp.get("p_fisher") is not None]
    return min(pvals) if pvals else None


def is_wuhan_b_backfill(res_row):
    """
    Sequence-intrinsic Wuhan-B amplicon backfill artifact detector:
    triggers solely from decoded segment identities when a near-verbatim Wuhan B tract
    (identity_pct >= 99.8) coexists with a diverged post-Wuhan Omicron+ trunk
    (>= 45 coding SNPs from Wuhan B), with zero inspection of year or reported_pango.
    """
    form = res_row.get("structural_formula") or ""
    segs = res_row.get("segments") or []
    has_verbatim_b = any(
        s.get("trunk") == "B" and float(s.get("identity_pct") or 0.0) >= 99.8
        for s in segs
    )
    has_diverged_partner = any(
        s.get("trunk") in DIVERGED_POST_WUHAN_TRUNKS
        for s in segs
    )
    if not (has_verbatim_b and has_diverged_partner):
        return False, False
    seg_backfill = ("[B:" in form) or ("| B:" in form)
    maj_backfill = form.startswith("B[")
    return seg_backfill, maj_backfill


# =============================================================================
# MODULE A: EXACT LIVE BV-BRC ARCHIVE CENSUS
# =============================================================================
def run_module_a():
    print("\n" + "=" * 88, flush=True)
    print("MODULE A: EXACT LIVE BV-BRC ARCHIVE CENSUS (SARS-CoV-2 & Influenza A)", flush=True)
    print("=" * 88, flush=True)
    t0 = time.time()

    sc2_base = "eq(taxon_lineage_ids,2697049)"
    sc2_comp = f"{sc2_base}&eq(genome_status,Complete)&gt(genome_length,29000)&eq(contigs,1)"

    sc2_census = {
        "taxon_lineage_id": "2697049",
        "all_genomes": bvbrc_count(sc2_base),
        "complete_status": bvbrc_count(f"{sc2_base}&eq(genome_status,Complete)"),
        "complete_1contig_gt29000": bvbrc_count(sc2_comp),
        "lineage_x_star_all": bvbrc_count(f"{sc2_base}&eq(lineage,X*)"),
        "lineage_x_star_complete_1contig_gt29000": bvbrc_count(f"{sc2_comp}&eq(lineage,X*)"),
        "lineage_missing_all": bvbrc_count(f"{sc2_base}&ne(lineage,*)"),
        "lineage_missing_complete_1contig_gt29000": bvbrc_count(f"{sc2_comp}&ne(lineage,*)"),
        "by_year": {},
    }

    for yr in range(2019, 2027):
        yr_all = bvbrc_count(f"{sc2_base}&eq(collection_year,{yr})")
        yr_comp = bvbrc_count(f"{sc2_comp}&eq(collection_year,{yr})")
        yr_x = bvbrc_count(f"{sc2_comp}&eq(collection_year,{yr})&eq(lineage,X*)")
        yr_miss = bvbrc_count(f"{sc2_comp}&eq(collection_year,{yr})&ne(lineage,*)")
        sc2_census["by_year"][str(yr)] = {
            "all": yr_all,
            "complete_1contig_gt29000": yr_comp,
            "lineage_x_star_complete_1contig_gt29000": yr_x,
            "lineage_missing_complete_1contig_gt29000": yr_miss,
        }
        print(
            f"  [SC2 {yr}] all={yr_all:,} | complete_1c_gt29k={yr_comp:,} | "
            f"X*={yr_x:,} | missing_lineage={yr_miss:,}",
            flush=True,
        )

    iav_base = "eq(taxon_lineage_ids,11320)"
    iav_comp = f"{iav_base}&eq(genome_status,Complete)&eq(contigs,1)"
    iav_census = {
        "taxon_lineage_id": "11320",
        "all_records": bvbrc_count(iav_base),
        "complete_1contig": bvbrc_count(iav_comp),
        "h5_2020_2026_complete_1contig": bvbrc_count(
            f"{iav_comp}&eq(h_type,5)&in(collection_year,(2020,2021,2022,2023,2024,2025,2026))"
        ),
        "by_segment": {},
        "by_subtype": {},
    }

    for seg in range(1, 9):
        s_all = bvbrc_count(f"{iav_base}&eq(segment,{seg})")
        s_comp = bvbrc_count(f"{iav_comp}&eq(segment,{seg})")
        iav_census["by_segment"][str(seg)] = {
            "segment_name": IAV_SEG_NAMES[str(seg)],
            "all": s_all,
            "complete_1contig": s_comp,
        }
        print(f"  [IAV Seg {seg} ({IAV_SEG_NAMES[str(seg)]})] all={s_all:,} | complete_1c={s_comp:,}", flush=True)

    for sub in ["H1N1", "H3N2", "H5N1", "H5N2", "H5N6", "H5N8", "H7N9", "H9N2"]:
        sub_all = bvbrc_count(f"{iav_base}&eq(subtype,{sub})")
        sub_comp = bvbrc_count(f"{iav_comp}&eq(subtype,{sub})")
        sub_ha_all = bvbrc_count(f"{iav_base}&eq(subtype,{sub})&eq(segment,4)")
        sub_ha_comp = bvbrc_count(f"{iav_comp}&eq(subtype,{sub})&eq(segment,4)")
        iav_census["by_subtype"][sub] = {
            "all_segments": sub_all,
            "complete_1contig_all_segments": sub_comp,
            "segment_4_ha_all": sub_ha_all,
            "segment_4_ha_complete_1contig": sub_ha_comp,
        }
        print(
            f"  [IAV {sub}] all_segs={sub_all:,} (comp_1c={sub_comp:,}) | "
            f"HA_seg4={sub_ha_all:,} (comp_1c={sub_ha_comp:,})",
            flush=True,
        )

    dt = round(time.time() - t0, 2)
    out = {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "runtime_sec": dt,
        "source_api": BVBRC_GENOME_URL,
        "sars_cov_2": sc2_census,
        "influenza_a": iav_census,
    }
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    with open(OUT_CENSUS, "w", encoding="utf-8") as f:
        json.dump(out, f, indent=2)
    print(f"[+] Module A complete in {dt:.2f}s -> {OUT_CENSUS}", flush=True)
    return out


# =============================================================================
# MODULE B: 100,000-GENOME IN-MEMORY SARS-CoV-2 STREAMING SCREEN
# =============================================================================
def build_sc2_100k_queue():
    print("\n" + "=" * 88, flush=True)
    print("MODULE B (STEP 1): BUILDING DETERMINISTIC 100,000-GENOME SARS-CoV-2 QUEUE", flush=True)
    print("=" * 88, flush=True)
    t0 = time.time()
    base = "eq(taxon_lineage_ids,2697049)&eq(genome_status,Complete)&gt(genome_length,29000)&eq(contigs,1)"
    fields = "genome_id,genbank_accessions,lineage,isolation_country,collection_date,collection_year,genome_length,strain"

    seen_gids = set()
    queue = []
    stratum_counts = {}

    def add_rows(rows, label):
        for r in rows:
            r["stratum"] = label
            queue.append(r)
        stratum_counts[label] = len(rows)
        print(f"  [Queue] {label:22s}: +{len(rows):6d} (total={len(queue):6d})", flush=True)

    # 1. All 2026 and 2025 complete 1-contig >29kb genomes
    rows_2026 = bvbrc_fetch_pages(f"{base}&eq(collection_year,2026)", fields, max_rows=None, page_size=10000, seen_global=seen_gids)
    add_rows(rows_2026, "2026_all")

    rows_2025 = bvbrc_fetch_pages(f"{base}&eq(collection_year,2025)", fields, max_rows=None, page_size=10000, seen_global=seen_gids)
    add_rows(rows_2025, "2025_all")

    # 2. Top up to 100,000: up to 5,000 from 2024 ne(lineage,*) + deterministic slices of 2024, 2023, 2022, 2021
    target_total = 100000
    rem = max(0, target_total - len(queue))
    if rem > 0:
        n_miss_2024 = min(5000, rem)
        rows_24_miss = bvbrc_fetch_pages(
            f"{base}&eq(collection_year,2024)&ne(lineage,*)",
            fields,
            max_rows=n_miss_2024,
            page_size=6000,
            seen_global=seen_gids,
        )
        add_rows(rows_24_miss, "2024_missing_lineage")

    rem2 = max(0, target_total - len(queue))
    slice_years = [2024, 2023, 2022, 2021]
    base_per_yr = rem2 // len(slice_years)
    extra = rem2 % len(slice_years)
    for idx, yr in enumerate(slice_years):
        need = base_per_yr + (1 if idx < extra else 0)
        if need <= 0:
            continue
        rows_yr = bvbrc_fetch_pages(
            f"{base}&eq(collection_year,{yr})",
            fields,
            max_rows=need,
            page_size=10000,
            seen_global=seen_gids,
        )
        add_rows(rows_yr, f"{yr}_slice")

    if len(queue) > target_total:
        queue = queue[:target_total]
    dt = round(time.time() - t0, 2)
    print(f"[+] Queue built: {len(queue):,} unique genome_ids in {dt:.2f}s", flush=True)
    if len(queue) != target_total:
        raise RuntimeError(f"Expected {target_total} unique genomes, got {len(queue)}")
    return queue, stratum_counts, dt


def init_mode_stats():
    return {
        "n_screened": 0,
        "runtime_sec": 0.0,
        "categories": Counter(),
        "called_lineages": Counter(),
        "by_year": defaultdict(Counter),
        "by_stratum": defaultdict(Counter),
        "total_mosaics": 0,
        "wuhan_b_backfill_mosaics": 0,
        "wuhan_b_majority_backfill_mosaics": 0,
        "genuine_mosaics": 0,
        "pango_audit": defaultdict(Counter),
        "pango_25_26_audit": defaultdict(Counter),
        "xbb_literal_total": 0,
        "xbb_literal_recovered": 0,
        "xbb_alias_total": 0,
        "xbb_alias_recovered": 0,
        "xdv_literal_total": 0,
        "xdv_literal_recovered": 0,
        "xdv_alias_total": 0,
        "xdv_alias_recovered": 0,
        "urf_clusters": defaultdict(list),
        "mislabeled_non_x": [],
    }


def update_mode_stats(stats, rows, meta_by_gid, n_bases_by_gid, n_ambig_by_gid, dt):
    stats["runtime_sec"] += dt
    for r in rows:
        sid = r["sequence_id"]
        gid = sid.split("|", 1)[0]
        m = meta_by_gid[gid]
        lin = (m.get("lineage") or "").strip()
        lin_label = lin if lin else "MISSING"
        lin_root = lin.split(".")[0] if lin else "MISSING"
        yr = m.get("collection_year")
        yr_str = str(yr) if yr is not None else "UNK"
        stratum = m.get("stratum", "UNK")

        nbp = int(r.get("n_breakpoints") or 0)
        cat = r.get("called_category") or "UNK"
        cal = r.get("called_lineage") or "UNK"
        form = r.get("structural_formula") or ""
        epi_form = r.get("epidemiological_formula") or ""
        barc = r.get("hxb2_barcode") or ""
        mapped_nt = int(r.get("mapped_coding_nt") or 0)
        trunks = r.get("constituent_trunks") or []
        min_p = min_bp_fisher_p(r)
        n_cnt = n_bases_by_gid.get(gid, 0)
        ambig_cnt = n_ambig_by_gid.get(gid, 0)

        seg_backfill, maj_backfill = is_wuhan_b_backfill(r)
        has_backfill = seg_backfill

        stats["n_screened"] += 1
        stats["categories"][cat] += 1
        stats["called_lineages"][cal] += 1

        stats["by_year"][yr_str]["total"] += 1
        stats["by_year"][yr_str][f"cat:{cat}"] += 1
        stats["by_stratum"][stratum]["total"] += 1
        stats["by_stratum"][stratum][f"cat:{cat}"] += 1

        if nbp > 0:
            stats["total_mosaics"] += 1
            stats["by_year"][yr_str]["mosaic_total"] += 1
            stats["by_stratum"][stratum]["mosaic_total"] += 1
            if has_backfill:
                stats["wuhan_b_backfill_mosaics"] += 1
                stats["by_year"][yr_str]["mosaic_wuhan_b_backfill"] += 1
                stats["by_stratum"][stratum]["mosaic_wuhan_b_backfill"] += 1
            else:
                stats["genuine_mosaics"] += 1
                stats["by_year"][yr_str]["mosaic_genuine"] += 1
                stats["by_stratum"][stratum]["mosaic_genuine"] += 1
            if maj_backfill:
                stats["wuhan_b_majority_backfill_mosaics"] += 1

        # Alias-masked and literal XBB / XDV recovery
        if lin_root.startswith("XBB"):
            stats["xbb_literal_total"] += 1
            if cal.startswith("XBB") or ("BJ.1" in form and "BM.1.1.1" in form):
                stats["xbb_literal_recovered"] += 1
        if lin_root in XBB_ALIAS_ROOTS:
            stats["xbb_alias_total"] += 1
            if cal.startswith("XBB") or ("BJ.1" in form and "BM.1.1.1" in form):
                stats["xbb_alias_recovered"] += 1
        if lin_root.startswith("XDV"):
            stats["xdv_literal_total"] += 1
            if nbp > 0 and ("BJ.1" in form and "JN.1" in form):
                stats["xdv_literal_recovered"] += 1
        if lin_root in XDV_ALIAS_ROOTS:
            stats["xdv_alias_total"] += 1
            if nbp > 0 and ("BJ.1" in form and "JN.1" in form):
                stats["xdv_alias_recovered"] += 1

        # Per-Pango-lineage audit for non-X* lineages (overall and 2025-2026)
        if lin and not lin_root.startswith("X"):
            for bucket in (
                [stats["pango_audit"][lin], stats["pango_25_26_audit"][lin]]
                if yr in (2025, 2026)
                else [stats["pango_audit"][lin]]
            ):
                bucket["total"] += 1
                if nbp > 0:
                    bucket["mosaic"] += 1
                    if has_backfill:
                        bucket["wuhan_b_backfill"] += 1
                    else:
                        bucket["genuine_mosaic"] += 1
                        bucket["f:" + form] += 1
                    bucket["fall:" + form] += 1

        is_known_x_alias = (lin_root in XBB_ALIAS_ROOTS) or (lin_root in XDV_ALIAS_ROOTS)

        # Mislabeled non-X* clinical submissions
        if (
            lin
            and not lin_root.startswith("X")
            and nbp >= 1
            and (not has_backfill)
            and n_cnt < 300
            and (min_p is not None and min_p <= 1e-5)
        ):
            stats["mislabeled_non_x"].append({
                "genome_id": gid,
                "acc": (m.get("genbank_accessions") or "").split(",")[0].strip(),
                "pango": lin,
                "is_known_x_alias": is_known_x_alias,
                "year": yr,
                "country": m.get("isolation_country") or "",
                "collection_date": m.get("collection_date") or "",
                "called_category": cat,
                "called_lineage": cal,
                "structural_formula": form,
                "epidemiological_formula": epi_form,
                "hxb2_barcode": barc,
                "n_breakpoints": nbp,
                "min_fisher_p": min_p,
                "mapped_coding_nt": mapped_nt,
                "n_bases": n_cnt,
                "n_non_acgt": ambig_cnt,
                "constituent_trunks": trunks,
            })

        # Unannotated / novel URF clusters (>=3 genomes sharing non-backfill structural_formula, <300 N bases)
        if cat == "URF" and nbp >= 1 and (not has_backfill) and n_cnt < 300:
            stats["urf_clusters"][form].append({
                "genome_id": gid,
                "seq_id": sid,
                "acc": (m.get("genbank_accessions") or "").split(",")[0].strip(),
                "pango": lin_label,
                "is_known_x_alias": is_known_x_alias,
                "year": yr,
                "country": m.get("isolation_country") or "",
                "collection_date": m.get("collection_date") or "",
                "epidemiological_formula": epi_form,
                "hxb2_barcode": barc,
                "n_breakpoints": nbp,
                "min_fisher_p": min_p,
                "mapped_coding_nt": mapped_nt,
                "n_bases": n_cnt,
                "n_non_acgt": ambig_cnt,
                "constituent_trunks": trunks,
            })


def summarize_sc2_mode(stats):
    n = stats["n_screened"]
    dt = stats["runtime_sec"]

    def build_lineage_audit(audit_dict):
        out = []
        high_pseudo_2526 = []
        for lin, cnt in audit_dict.items():
            tot = cnt["total"]
            if tot < 30:
                continue
            mos = cnt["mosaic"]
            wbf = cnt["wuhan_b_backfill"]
            gen = cnt["genuine_mosaic"]
            pct_mos = round(100.0 * mos / tot, 2)
            pct_gen = round(100.0 * gen / tot, 2)
            top_gen_f = max(((k[2:], v) for k, v in cnt.items() if k.startswith("f:")), key=lambda x: x[1], default=("", 0))
            top_all_f = max(((k[5:], v) for k, v in cnt.items() if k.startswith("fall:")), key=lambda x: x[1], default=("", 0))
            lin_root = lin.split(".")[0]
            is_alias = (lin_root in XBB_ALIAS_ROOTS) or (lin_root in XDV_ALIAS_ROOTS)
            rec = {
                "lineage": lin,
                "total": tot,
                "mosaic_total": mos,
                "mosaic_pct": pct_mos,
                "wuhan_b_backfill": wbf,
                "genuine_mosaic": gen,
                "genuine_mosaic_pct": pct_gen,
                "is_known_x_alias": is_alias,
                "top_genuine_formula": top_gen_f,
                "top_overall_formula": top_all_f,
            }
            out.append(rec)
            if pct_mos > 50.0 or pct_gen > 50.0:
                high_pseudo_2526.append(rec)
        out.sort(key=lambda x: (-x["genuine_mosaic_pct"], -x["mosaic_pct"], -x["total"]))
        high_pseudo_2526.sort(key=lambda x: (-x["mosaic_pct"], -x["total"]))
        return out, high_pseudo_2526

    audit_all, high_pseudo_all = build_lineage_audit(stats["pango_audit"])
    audit_25_26, high_pseudo_25_26 = build_lineage_audit(stats["pango_25_26_audit"])

    urf_list = []
    n_3parent_clusters = 0
    n_3parent_genomes = 0
    for form, items in sorted(stats["urf_clusters"].items(), key=lambda x: -len(x[1])):
        if len(items) < 3:
            continue
        trunks = items[0]["constituent_trunks"]
        is_3parent = len(trunks) >= 3
        if is_3parent:
            n_3parent_clusters += 1
            n_3parent_genomes += len(items)
        pangos = Counter(x["pango"] for x in items).most_common(8)
        countries = Counter(x["country"] for x in items if x["country"]).most_common(6)
        year_counts = dict(sorted(Counter(str(x["year"]) for x in items if x["year"] is not None).items()))
        years = sorted(year_counts.keys())
        ns = sorted(x["n_bases"] for x in items)
        ps = sorted(x["min_fisher_p"] for x in items if x["min_fisher_p"] is not None)
        n_alias = sum(1 for x in items if x["is_known_x_alias"])
        n_x_reported = sum(1 for x in items if x["pango"].startswith("X"))
        urf_list.append({
            "structural_formula": form,
            "epidemiological_formula": items[0]["epidemiological_formula"],
            "constituent_trunks": trunks,
            "n_parents": len(trunks),
            "is_3parent_or_more": is_3parent,
            "n_genomes": len(items),
            "n_2025_2026": year_counts.get("2025", 0) + year_counts.get("2026", 0),
            "year_counts": year_counts,
            "n_known_x_alias_reported": n_alias,
            "n_literal_x_reported": n_x_reported,
            "n_non_x_non_alias_reported": len(items) - n_alias - n_x_reported,
            "median_n_bases": ns[len(ns) // 2],
            "min_n_bases": ns[0],
            "min_fisher_p": ps[0] if ps else None,
            "median_fisher_p": ps[len(ps) // 2] if ps else None,
            "top_reported_pango": pangos,
            "top_countries": countries,
            "years": years,
            "representative_barcode": items[0]["hxb2_barcode"],
            "example_accs": [x["acc"] for x in items[:5] if x["acc"]],
        })

    mislabeled = stats["mislabeled_non_x"]
    mislabeled_non_alias = [x for x in mislabeled if not x["is_known_x_alias"]]
    mislabeled_alias = [x for x in mislabeled if x["is_known_x_alias"]]

    return {
        "n_screened": n,
        "eggman_runtime_sec": round(dt, 3),
        "throughput_seq_per_sec": round(n / dt, 1) if dt > 0 else 0.0,
        "categories": dict(stats["categories"]),
        "top_called_lineages": stats["called_lineages"].most_common(25),
        "mosaic_breakdown": {
            "total_mosaics_n_bp_ge_1": stats["total_mosaics"],
            "wuhan_b_segment_backfill_post2021": stats["wuhan_b_backfill_mosaics"],
            "wuhan_b_majority_backfill_post2021": stats["wuhan_b_majority_backfill_mosaics"],
            "genuine_biological_mosaics": stats["genuine_mosaics"],
            "genuine_mosaic_pct": round(100.0 * stats["genuine_mosaics"] / n, 3) if n else 0.0,
        },
        "by_year": {k: dict(v) for k, v in sorted(stats["by_year"].items())},
        "by_stratum": {k: dict(v) for k, v in stats["by_stratum"].items()},
        "alias_recovery": {
            "xbb_literal_total": stats["xbb_literal_total"],
            "xbb_literal_recovered": stats["xbb_literal_recovered"],
            "xbb_literal_pct": round(100.0 * stats["xbb_literal_recovered"] / stats["xbb_literal_total"], 2) if stats["xbb_literal_total"] else None,
            "xbb_alias_masked_total": stats["xbb_alias_total"],
            "xbb_alias_masked_recovered": stats["xbb_alias_recovered"],
            "xbb_alias_masked_pct": round(100.0 * stats["xbb_alias_recovered"] / stats["xbb_alias_total"], 2) if stats["xbb_alias_total"] else None,
            "xdv_literal_total": stats["xdv_literal_total"],
            "xdv_literal_recovered": stats["xdv_literal_recovered"],
            "xdv_literal_pct": round(100.0 * stats["xdv_literal_recovered"] / stats["xdv_literal_total"], 2) if stats["xdv_literal_total"] else None,
            "xdv_alias_masked_total": stats["xdv_alias_total"],
            "xdv_alias_masked_recovered": stats["xdv_alias_recovered"],
            "xdv_alias_masked_pct": round(100.0 * stats["xdv_alias_recovered"] / stats["xdv_alias_total"], 2) if stats["xdv_alias_total"] else None,
        },
        "pango_pseudo_mosaic_audit_all_n_ge_30": audit_all,
        "pango_pseudo_mosaic_audit_2025_2026_n_ge_30": audit_25_26,
        "clonal_waves_all_gt_50pct_mosaic": high_pseudo_all,
        "clonal_waves_2025_2026_gt_50pct_mosaic": high_pseudo_25_26,
        "mislabeled_non_x_summary": {
            "total_qualifying": len(mislabeled),
            "known_x_alias_count": len(mislabeled_alias),
            "true_non_alias_count": len(mislabeled_non_alias),
            "top_non_alias_reported_pango": Counter(x["pango"] for x in mislabeled_non_alias).most_common(25),
            "top_non_alias_formulas": Counter(x["structural_formula"] for x in mislabeled_non_alias).most_common(25),
        },
        "novel_urf_clusters_ge_3": {
            "n_clusters": len(urf_list),
            "total_genomes_in_clusters": sum(c["n_genomes"] for c in urf_list),
            "n_3parent_or_more_clusters": n_3parent_clusters,
            "n_3parent_or_more_genomes": n_3parent_genomes,
            "clusters": urf_list,
        },
    }


def read_cached_chunk_fasta_gz(path):
    seqs = {}
    cur = None
    buf = []
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as f:
        for raw_ln in f:
            ln = raw_ln.strip()
            if not ln:
                continue
            if ln.startswith(">"):
                if cur is not None and buf:
                    seqs[cur] = "".join(buf).upper()
                cur = ln[1:].split()[0].split("|")[0]
                buf = []
            else:
                buf.append(ln)
        if cur is not None and buf:
            seqs[cur] = "".join(buf).upper()
    return seqs


def run_module_b(cache_dir=None):
    t_mod_start = time.time()
    cache_s100 = Path(cache_dir) / "sc2_100k_stream" if cache_dir else None
    use_cached_sc2 = (
        cache_s100 is not None
        and (cache_s100 / "queue_100k.json").exists()
        and all((cache_s100 / f"chunk_{i:02d}.fasta.gz").exists() for i in range(10))
    )
    if use_cached_sc2:
        t_q0 = time.time()
        with open(cache_s100 / "queue_100k.json", "r", encoding="utf-8") as qf:
            queue = json.load(qf)
        for r in queue:
            if r.get("stratum") == "2024_missing":
                r["stratum"] = "2024_missing_lineage"
        stratum_counts = dict(Counter(r.get("stratum", "UNK") for r in queue))
        dt_queue = round(time.time() - t_q0, 2)
        print(f"\n[+] Loaded cached 100k SC2 queue ({len(queue):,} genomes) from {cache_s100} in {dt_queue:.2f}s", flush=True)
    else:
        queue, stratum_counts, dt_queue = build_sc2_100k_queue()

    print("\n" + "=" * 88, flush=True)
    print("MODULE B (STEP 2): STREAMING 100,000 GENOMES IN 10,000-GENOME MACRO-CHUNKS", flush=True)
    print("=" * 88, flush=True)

    meta_by_gid = {r["genome_id"]: r for r in queue}
    n_bases_by_gid = {}
    n_ambig_by_gid = {}

    stats_def = init_mode_stats()
    stats_cov = init_mode_stats()

    # Keep small bounded sequences ONLY for novel non-backfill URF candidates (<300 N bases)
    # Keyed by (formula, gid) -> (n_bases, mapped_nt, seq, meta_row, eggman_row)
    kept_urf_seqs_cov = defaultdict(list)

    macro_size = 10000
    batch_size = 500
    cur_workers = 8
    chunk_reports = []
    total_stream_fetch_sec = 0.0
    total_joined = 0

    for c_idx in range(0, len(queue), macro_size):
        chunk_num = (c_idx // macro_size) + 1
        chunk_rows = queue[c_idx : c_idx + macro_size]
        chunk_gids = [r["genome_id"] for r in chunk_rows]
        valid_gids = set(chunk_gids)

        if use_cached_sc2:
            t_fetch0 = time.time()
            seq_by_gid = read_cached_chunk_fasta_gz(cache_s100 / f"chunk_{chunk_num - 1:02d}.fasta.gz")
            seq_by_gid = {gid: seq_by_gid[gid] for gid in chunk_gids if gid in seq_by_gid}
            dt_fetch = time.time() - t_fetch0
            hdr_cnt = len(seq_by_gid)
            joined_cnt = len(seq_by_gid)
            unmatched = []
        else:
            acc_to_gid = {}
            for r in chunk_rows:
                accs = (r.get("genbank_accessions") or "").split(",")
                for a in accs:
                    a_clean = a.strip()
                    if a_clean:
                        acc_to_gid[a_clean] = r["genome_id"]
                        acc_to_gid[a_clean.split(".")[0]] = r["genome_id"]

            batches = [chunk_gids[i : i + batch_size] for i in range(0, len(chunk_gids), batch_size)]
            t_fetch0 = time.time()
            fasta_texts, cur_workers = fetch_fasta_batches_adaptive(batches, workers=cur_workers)
            dt_fetch = time.time() - t_fetch0
            seq_by_gid, hdr_cnt, joined_cnt, unmatched = parse_fasta_and_join(fasta_texts, valid_gids, acc_to_gid)
            del fasta_texts
        total_stream_fetch_sec += dt_fetch

        join_rate_vs_expected = len(seq_by_gid) / len(chunk_rows)
        join_rate_vs_headers = (joined_cnt / hdr_cnt) if hdr_cnt else 0.0
        total_joined += len(seq_by_gid)

        print(
            f"  [SC2 Chunk {chunk_num:2d}/10] fetched={len(seq_by_gid)}/{len(chunk_rows)} "
            f"({100*join_rate_vs_expected:.2f}% expected, {100*join_rate_vs_headers:.2f}% headers) "
            f"in {dt_fetch:.2f}s (workers={cur_workers})",
            flush=True,
        )
        if unmatched:
            print(f"    Unmatched header samples: {unmatched[:2]}", flush=True)

        if chunk_num == 1 and (join_rate_vs_expected < 0.99 or join_rate_vs_headers < 0.99):
            raise RuntimeError(
                f"ABORT: Chunk 1 header-to-metadata join rate below 99% "
                f"(expected={100*join_rate_vs_expected:.2f}%, headers={100*join_rate_vs_headers:.2f}%)"
            )

        # Write temporary normalized FASTA for EGGMAN, run Default & CoV-tuned, then delete immediately
        fd, tmp_fa = tempfile.mkstemp(suffix=".sc2.fasta")
        os.close(fd)
        try:
            with open(tmp_fa, "w", encoding="utf-8") as tf:
                for r in chunk_rows:
                    gid = r["genome_id"]
                    seq = seq_by_gid.get(gid)
                    if not seq:
                        continue
                    n_cnt = seq.count("N")
                    ambig_cnt = sum(1 for c in seq if c not in "ACGT")
                    n_bases_by_gid[gid] = n_cnt
                    n_ambig_by_gid[gid] = ambig_cnt
                    acc = (r.get("genbank_accessions") or "NA").split(",")[0].strip().replace(" ", "_") or "NA"
                    tf.write(f">{gid}|{acc}\n{seq}\n")

            def_rows, dt_def = run_eggman_json(SC2_PACK, tmp_fa)
            cov_rows, dt_cov = run_eggman_json(SC2_PACK, tmp_fa, gamma=3.5, min_snps=3, p_thresh=0.01)
        finally:
            if os.path.exists(tmp_fa):
                os.remove(tmp_fa)

        update_mode_stats(stats_def, def_rows, meta_by_gid, n_bases_by_gid, n_ambig_by_gid, dt_def)
        update_mode_stats(stats_cov, cov_rows, meta_by_gid, n_bases_by_gid, n_ambig_by_gid, dt_cov)

        # Retain sequences in memory ONLY for non-backfill URF candidates (<300 N bases) up to 15 per formula
        for er in cov_rows:
            if er.get("called_category") == "URF" and int(er.get("n_breakpoints") or 0) >= 1:
                gid = er["sequence_id"].split("|", 1)[0]
                m = meta_by_gid[gid]
                form = er.get("structural_formula") or ""
                seg_bf, _ = is_wuhan_b_backfill(er)
                n_cnt = n_bases_by_gid.get(gid, 99999)
                if (not seg_bf) and n_cnt < 300 and gid in seq_by_gid:
                    bucket = kept_urf_seqs_cov[form]
                    item = (n_cnt, -int(er.get("mapped_coding_nt") or 0), gid, er["sequence_id"], seq_by_gid[gid], er)
                    bucket.append(item)
                    if len(bucket) > 15:
                        bucket.sort(key=lambda x: (x[0], x[1], x[2]))
                        del bucket[15:]

        del seq_by_gid
        print(
            f"    EGGMAN Default: {len(def_rows)} seqs in {dt_def:.2f}s ({len(def_rows)/dt_def:.0f} seq/s) | "
            f"CoV-tuned: {len(cov_rows)} seqs in {dt_cov:.2f}s ({len(cov_rows)/dt_cov:.0f} seq/s)",
            flush=True,
        )
        chunk_reports.append({
            "chunk": chunk_num,
            "requested": len(chunk_rows),
            "joined_unique": len(def_rows),
            "headers_received": hdr_cnt,
            "headers_joined": joined_cnt,
            "join_rate_expected": round(join_rate_vs_expected, 5),
            "join_rate_headers": round(join_rate_vs_headers, 5),
            "fetch_sec": round(dt_fetch, 2),
            "eggman_default_sec": round(dt_def, 3),
            "eggman_cov_tuned_sec": round(dt_cov, 3),
        })

    sum_def = summarize_sc2_mode(stats_def)
    sum_cov = summarize_sc2_mode(stats_cov)

    # Write representative sequences (up to 10 per novel URF cluster with >=3 genomes) to sc2_novel_hits.fasta.gz
    HITS_DIR.mkdir(parents=True, exist_ok=True)
    novel_clusters_cov = sum_cov["novel_urf_clusters_ge_3"]["clusters"]
    saved_rep_count = 0
    saved_cluster_count = 0
    with gzip.open(OUT_SC2_FASTA_GZ, "wt", encoding="utf-8") as gz:
        for cl in novel_clusters_cov:
            form = cl["structural_formula"]
            reps = sorted(kept_urf_seqs_cov.get(form, []), key=lambda x: (x[0], x[1], x[2]))[:10]
            if reps:
                saved_cluster_count += 1
            for n_cnt, _, gid, sid, seq, er in reps:
                m = meta_by_gid[gid]
                lin = (m.get("lineage") or "MISSING").strip()
                yr = m.get("collection_year") or "UNK"
                ctry = (m.get("isolation_country") or "UNK").replace(" ", "_")
                gz.write(f">{sid} pango={lin} year={yr} country={ctry} n_bases={n_cnt} formula={form}\n{seq}\n")
                saved_rep_count += 1

    # Write results/bvbrc_100k_sc2_novel_and_mislabeled.tsv
    # Includes: (1) NOVEL_URF_CLUSTER summary rows, (2) NOVEL_URF_MEMBER rows (up to 10 reps/cluster),
    # and (3) MISLABELED_NON_X rows (sorted with non-alias entries first, then min_fisher_p ascending)
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    cov_mislabeled_by_gid = {x["genome_id"]: x for x in stats_cov["mislabeled_non_x"]}
    def_mislabeled_gids = {x["genome_id"] for x in stats_def["mislabeled_non_x"]}

    tsv_cols = [
        "record_type",
        "cluster_size",
        "genome_id",
        "genbank_accession",
        "reported_pango",
        "is_known_x_alias",
        "collection_year",
        "collection_date",
        "isolation_country",
        "called_category",
        "called_lineage",
        "structural_formula",
        "epidemiological_formula",
        "n_parents",
        "n_breakpoints",
        "min_fisher_p",
        "also_called_default_mode",
        "mapped_coding_nt",
        "n_bases",
        "hxb2_barcode",
    ]
    tsv_rows_written = 0
    with open(OUT_SC2_TSV, "w", newline="", encoding="utf-8") as tf:
        w = csv.DictWriter(tf, fieldnames=tsv_cols, delimiter="\t")
        w.writeheader()
        for cl in novel_clusters_cov:
            form = cl["structural_formula"]
            w.writerow({
                "record_type": "NOVEL_URF_CLUSTER",
                "cluster_size": cl["n_genomes"],
                "genome_id": f"CLUSTER_N={cl['n_genomes']}",
                "genbank_accession": ",".join(cl["example_accs"]),
                "reported_pango": ";".join(f"{p}:{c}" for p, c in cl["top_reported_pango"][:4]),
                "is_known_x_alias": cl["n_known_x_alias_reported"] > 0,
                "collection_year": ",".join(cl["years"]),
                "collection_date": "",
                "isolation_country": ";".join(f"{ct}:{c}" for ct, c in cl["top_countries"][:3]),
                "called_category": "URF",
                "called_lineage": "URF",
                "structural_formula": form,
                "epidemiological_formula": cl["epidemiological_formula"],
                "n_parents": cl["n_parents"],
                "n_breakpoints": cl["representative_barcode"].count("|"),
                "min_fisher_p": cl["min_fisher_p"],
                "also_called_default_mode": "",
                "mapped_coding_nt": "",
                "n_bases": cl["median_n_bases"],
                "hxb2_barcode": cl["representative_barcode"],
            })
            tsv_rows_written += 1
            for n_cnt, _, gid, sid, _, er in sorted(kept_urf_seqs_cov.get(form, []), key=lambda x: (x[0], x[1], x[2]))[:10]:
                m = meta_by_gid[gid]
                lin = (m.get("lineage") or "MISSING").strip()
                lin_root = lin.split(".")[0]
                w.writerow({
                    "record_type": "NOVEL_URF_REP",
                    "cluster_size": cl["n_genomes"],
                    "genome_id": gid,
                    "genbank_accession": (m.get("genbank_accessions") or "").split(",")[0].strip(),
                    "reported_pango": lin,
                    "is_known_x_alias": (lin_root in XBB_ALIAS_ROOTS) or (lin_root in XDV_ALIAS_ROOTS),
                    "collection_year": m.get("collection_year") or "",
                    "collection_date": m.get("collection_date") or "",
                    "isolation_country": m.get("isolation_country") or "",
                    "called_category": er.get("called_category"),
                    "called_lineage": er.get("called_lineage"),
                    "structural_formula": form,
                    "epidemiological_formula": er.get("epidemiological_formula") or "",
                    "n_parents": len(er.get("constituent_trunks") or []),
                    "n_breakpoints": er.get("n_breakpoints"),
                    "min_fisher_p": min_bp_fisher_p(er),
                    "also_called_default_mode": gid in def_mislabeled_gids,
                    "mapped_coding_nt": er.get("mapped_coding_nt"),
                    "n_bases": n_cnt,
                    "hxb2_barcode": er.get("hxb2_barcode") or "",
                })
                tsv_rows_written += 1

        # Mislabeled non-X* clinical submissions (non-alias first, then alias; sorted by min_fisher_p)
        sorted_mislabeled = sorted(
            cov_mislabeled_by_gid.values(),
            key=lambda x: (x["is_known_x_alias"], x["min_fisher_p"], x["n_bases"], x["genome_id"]),
        )
        for item in sorted_mislabeled:
            w.writerow({
                "record_type": "MISLABELED_NON_X",
                "cluster_size": "",
                "genome_id": item["genome_id"],
                "genbank_accession": item["acc"],
                "reported_pango": item["pango"],
                "is_known_x_alias": item["is_known_x_alias"],
                "collection_year": item["collection_year"] if "collection_year" in item else (item["year"] or ""),
                "collection_date": item["collection_date"],
                "isolation_country": item["country"],
                "called_category": item["called_category"],
                "called_lineage": item["called_lineage"],
                "structural_formula": item["structural_formula"],
                "epidemiological_formula": item["epidemiological_formula"],
                "n_parents": len(item["constituent_trunks"]),
                "n_breakpoints": item["n_breakpoints"],
                "min_fisher_p": item["min_fisher_p"],
                "also_called_default_mode": item["genome_id"] in def_mislabeled_gids,
                "mapped_coding_nt": item["mapped_coding_nt"],
                "n_bases": item["n_bases"],
                "hxb2_barcode": item["hxb2_barcode"],
            })
            tsv_rows_written += 1

    total_dt = round(time.time() - t_mod_start, 2)
    out = {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "queue_size_requested": len(queue),
        "total_genomes_screened": total_joined,
        "stratum_counts": stratum_counts,
        "timings_sec": {
            "metadata_queue_sec": dt_queue,
            "fasta_stream_download_sec": round(total_stream_fetch_sec, 2),
            "eggman_default_sec": sum_def["eggman_runtime_sec"],
            "eggman_cov_tuned_sec": sum_cov["eggman_runtime_sec"],
            "total_module_b_sec": total_dt,
        },
        "chunk_reports": chunk_reports,
        "outputs": {
            "novel_urf_fasta_gz": str(OUT_SC2_FASTA_GZ),
            "novel_urf_clusters_saved_to_fasta": saved_cluster_count,
            "novel_urf_representatives_saved_to_fasta": saved_rep_count,
            "novel_and_mislabeled_tsv": str(OUT_SC2_TSV),
            "tsv_rows_written": tsv_rows_written,
        },
        "eggman_default": sum_def,
        "eggman_cov_tuned": sum_cov,
    }
    with open(OUT_SC2_JSON, "w", encoding="utf-8") as f:
        json.dump(out, f, indent=2)

    print(
        f"[+] Module B complete in {total_dt:.2f}s | screened={total_joined:,} | "
        f"Default={sum_def['throughput_seq_per_sec']} seq/s | CoV-tuned={sum_cov['throughput_seq_per_sec']} seq/s",
        flush=True,
    )
    return out


# =============================================================================
# MODULE C: EXHAUSTIVE 2020-2026 H5Nx (h_type=5) 8-SEGMENT IAV STREAM
# =============================================================================
def run_module_c():
    print("\n" + "=" * 88, flush=True)
    print("MODULE C: EXHAUSTIVE 2020-2026 H5Nx (h_type=5) COMPLETE 8-SEGMENT IAV STREAM", flush=True)
    print("=" * 88, flush=True)
    t_mod_start = time.time()

    rql = (
        "eq(taxon_lineage_ids,11320)&eq(genome_status,Complete)&eq(contigs,1)"
        "&eq(h_type,5)&in(collection_year,(2020,2021,2022,2023,2024,2025,2026))"
    )
    fields = (
        "genome_id,biosample_accession,strain,segment,genbank_accessions,"
        "subtype,isolation_country,collection_date,collection_year,host_name,genome_length"
    )
    expected_total = bvbrc_count(rql)
    t_meta0 = time.time()
    rows = bvbrc_fetch_pages(rql, fields, max_rows=None, page_size=10000)
    dt_meta = round(time.time() - t_meta0, 2)
    print(f"  [IAV Meta] Fetched {len(rows):,} / {expected_total:,} segment rows in {dt_meta:.2f}s", flush=True)

    groups = defaultdict(lambda: defaultdict(list))
    skipped_empty_id = 0
    invalid_seg = 0
    grouped_by_biosample = 0
    grouped_by_strain = 0
    seen_group_mode = {}

    for r in rows:
        bs = (r.get("biosample_accession") or "").strip()
        st = (r.get("strain") or "").strip()
        if bs:
            gkey = ("bs", bs)
        elif st:
            ctry = (r.get("isolation_country") or "").strip()
            cdate = (r.get("collection_date") or "").strip()
            gkey = ("st", st.lower(), ctry, cdate)
        else:
            skipped_empty_id += 1
            continue

        seg = r.get("segment")
        if isinstance(seg, list):
            seg_str = str(seg[0]).strip() if seg else ""
        else:
            seg_str = str(seg or "").strip()

        if seg_str not in IAV_SEG_NAMES:
            invalid_seg += 1
            continue

        if gkey not in seen_group_mode:
            seen_group_mode[gkey] = gkey[0]
            if gkey[0] == "bs":
                grouped_by_biosample += 1
            else:
                grouped_by_strain += 1

        groups[gkey][seg_str].append(r)

    complete_isolates = []
    duplicate_slot_count = 0
    for idx_g, (gkey, segmap) in enumerate(sorted(groups.items(), key=lambda x: str(x[0]))):
        if set(segmap.keys()) != set(IAV_SEG_NAMES.keys()):
            continue
        chosen = {}
        for s in "12345678":
            cands = segmap[s]
            if len(cands) > 1:
                duplicate_slot_count += 1
            chosen[s] = max(cands, key=lambda x: (int(x.get("genome_length") or 0), x["genome_id"]))

        # Derive isolate-level metadata primarily from HA ("4"), falling back across segments
        ref_row = chosen["4"]
        def pick_field(fn):
            val = ref_row.get(fn)
            if val is not None and str(val).strip():
                return val
            for s in "1235678":
                v = chosen[s].get(fn)
                if v is not None and str(v).strip():
                    return v
            return ""

        bs = str(pick_field("biosample_accession")).strip()
        strain = str(pick_field("strain")).strip()
        subtype = str(pick_field("subtype")).strip() or "H5Nx"
        country = str(pick_field("isolation_country")).strip() or "UNK"
        cdate = str(pick_field("collection_date")).strip()
        cyear = pick_field("collection_year")
        host = str(pick_field("host_name")).strip() or "UNK"
        ha_acc = str(chosen["4"].get("genbank_accessions") or "").split(",")[0].strip()

        id_tok = re.sub(r"[^A-Za-z0-9_.-]+", "_", bs or strain)[:60]
        uid = f"IAV_{len(complete_isolates):05d}|{id_tok}"
        complete_isolates.append({
            "uid": uid,
            "biosample_accession": bs,
            "strain": strain,
            "subtype": subtype,
            "country": country,
            "collection_date": cdate,
            "collection_year": cyear,
            "host_name": host,
            "ha_accession": ha_acc,
            "segments": chosen,
        })

    print(
        f"  [IAV Assembly] groups={len(groups):,} (biosample={grouped_by_biosample:,}, strain={grouped_by_strain:,}) | "
        f"complete_8seg={len(complete_isolates):,} | duplicate_slots_resolved={duplicate_slot_count:,} | "
        f"skipped_empty_id={skipped_empty_id}",
        flush=True,
    )

    # Stream FASTA sequences for all 8 * N_complete segments
    needed_gids = []
    valid_gids = set()
    acc_to_gid = {}
    for iso in complete_isolates:
        for s in "12345678":
            sr = iso["segments"][s]
            gid = sr["genome_id"]
            needed_gids.append(gid)
            valid_gids.add(gid)
            for a in (sr.get("genbank_accessions") or "").split(","):
                ac = a.strip()
                if ac:
                    acc_to_gid[ac] = gid
                    acc_to_gid[ac.split(".")[0]] = gid

    batches = [needed_gids[i : i + 500] for i in range(0, len(needed_gids), 500)]
    print(f"  [IAV FASTA] Streaming {len(needed_gids):,} segment sequences in {len(batches)} batches...", flush=True)
    t_fetch0 = time.time()
    fasta_texts, cur_workers = fetch_fasta_batches_adaptive(batches, workers=8)
    dt_fetch = round(time.time() - t_fetch0, 2)
    seq_by_gid, hdr_cnt, joined_cnt, unmatched = parse_fasta_and_join(fasta_texts, valid_gids, acc_to_gid)
    del fasta_texts

    join_rate_expected = len(seq_by_gid) / len(needed_gids) if needed_gids else 0.0
    join_rate_headers = (joined_cnt / hdr_cnt) if hdr_cnt else 0.0
    print(
        f"  [IAV FASTA] Joined {len(seq_by_gid):,}/{len(needed_gids):,} ({100*join_rate_expected:.2f}% expected, "
        f"{100*join_rate_headers:.2f}% headers) in {dt_fetch:.2f}s (workers={cur_workers})",
        flush=True,
    )
    if join_rate_expected < 0.99 or join_rate_headers < 0.99:
        raise RuntimeError(
            f"ABORT: IAV segment FASTA join rate below 99% "
            f"(expected={100*join_rate_expected:.2f}%, headers={100*join_rate_headers:.2f}%)"
        )

    # Assemble concatenated 8-segment sequences and run EGGMAN in batches of 1000
    assembled_isolates = []
    missing_seg_isolates = 0
    for iso in complete_isolates:
        parts = []
        ok = True
        for s in "12345678":
            gid = iso["segments"][s]["genome_id"]
            sq = seq_by_gid.get(gid)
            if not sq:
                ok = False
                break
            parts.append(sq)
        if not ok:
            missing_seg_isolates += 1
            continue
        concat = "".join(parts)
        iso["concat_len"] = len(concat)
        iso["n_bases"] = concat.count("N")
        iso["n_non_acgt"] = sum(1 for c in concat if c not in "ACGT")
        iso["concat_seq"] = concat
        assembled_isolates.append(iso)
    del seq_by_gid

    iso_by_uid = {iso["uid"]: iso for iso in assembled_isolates}
    egg_results = []
    dt_eggman = 0.0
    for b_start in range(0, len(assembled_isolates), 1000):
        batch_isos = assembled_isolates[b_start : b_start + 1000]
        fd, tmp_fa = tempfile.mkstemp(suffix=".iav.fasta")
        os.close(fd)
        try:
            with open(tmp_fa, "w", encoding="utf-8") as tf:
                for iso in batch_isos:
                    tf.write(f">{iso['uid']}\n{iso['concat_seq']}\n")
            erows, dt_b = run_eggman_json(IAV_PACK, tmp_fa)
            dt_eggman += dt_b
            egg_results.extend(erows)
        finally:
            if os.path.exists(tmp_fa):
                os.remove(tmp_fa)

    print(
        f"  [IAV EGGMAN] Screened {len(egg_results):,} 8-segment isolates in {dt_eggman:.3f}s "
        f"({len(egg_results)/dt_eggman:.1f} iso/s)",
        flush=True,
    )

    # Summarize all assembled and high-coverage (mapped_coding_nt >= 11500)
    all_cats = Counter()
    all_calls = Counter()
    hc_cats = Counter()
    hc_calls = Counter()
    hc_subtypes = Counter()
    hc_countries = Counter()
    hc_hosts = Counter()
    hc_years = defaultdict(Counter)
    hc_urf_clusters = defaultdict(list)

    n_high_cov = 0
    for er in egg_results:
        uid = er["sequence_id"]
        iso = iso_by_uid[uid]
        cat = er.get("called_category") or "UNK"
        cal = er.get("called_lineage") or "UNK"
        mapped_nt = int(er.get("mapped_coding_nt") or 0)
        nbp = int(er.get("n_breakpoints") or 0)
        form = er.get("structural_formula") or ""
        epi_form = er.get("epidemiological_formula") or ""
        barc = er.get("hxb2_barcode") or ""
        trunks = er.get("constituent_trunks") or []
        min_p = min_bp_fisher_p(er)

        all_cats[cat] += 1
        all_calls[cal] += 1

        if mapped_nt >= 11500:
            n_high_cov += 1
            hc_cats[cat] += 1
            hc_calls[cal] += 1
            hc_subtypes[iso["subtype"]] += 1
            hc_countries[iso["country"]] += 1
            hc_hosts[iso["host_name"]] += 1
            yr_s = str(iso["collection_year"]) if iso["collection_year"] is not None else "UNK"
            hc_years[yr_s]["total"] += 1
            hc_years[yr_s][f"cat:{cat}"] += 1
            hc_years[yr_s][f"call:{cal}"] += 1

            if cat == "URF" and nbp >= 1:
                hc_urf_clusters[form].append({
                    "uid": uid,
                    "biosample": iso["biosample_accession"],
                    "strain": iso["strain"],
                    "ha_acc": iso["ha_accession"],
                    "subtype": iso["subtype"],
                    "country": iso["country"],
                    "year": iso["collection_year"],
                    "collection_date": iso["collection_date"],
                    "host": iso["host_name"],
                    "mapped_coding_nt": mapped_nt,
                    "n_bases": iso["n_bases"],
                    "n_non_acgt": iso["n_non_acgt"],
                    "structural_formula": form,
                    "epidemiological_formula": epi_form,
                    "backbone_trunk": er.get("backbone_trunk") or "",
                    "constituent_trunks": trunks,
                    "n_segments": er.get("n_segments"),
                    "n_breakpoints": nbp,
                    "min_fisher_p": min_p,
                    "hxb2_barcode": barc,
                })

    # Cluster novel high-coverage URF constellations (>= 3 isolates)
    novel_constellations = []
    HITS_DIR.mkdir(parents=True, exist_ok=True)
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    saved_iav_reps = 0

    with gzip.open(OUT_IAV_FASTA_GZ, "wt", encoding="utf-8") as gz:
        for form, items in sorted(hc_urf_clusters.items(), key=lambda x: -len(x[1])):
            if len(items) < 3:
                continue
            ps = sorted(x["min_fisher_p"] for x in items if x["min_fisher_p"] is not None)
            m_nts = sorted(x["mapped_coding_nt"] for x in items)
            subtypes = Counter(x["subtype"] for x in items).most_common(6)
            countries = Counter(x["country"] for x in items).most_common(6)
            hosts = Counter(x["host"] for x in items).most_common(8)
            years = sorted(set(str(x["year"]) for x in items if x["year"] is not None))
            epi_forms = Counter(x["epidemiological_formula"] for x in items).most_common(3)
            cl_rec = {
                "structural_formula": form,
                "epidemiological_formula": epi_forms[0][0] if epi_forms else "",
                "backbone_trunk": items[0]["backbone_trunk"],
                "constituent_trunks": items[0]["constituent_trunks"],
                "n_isolates": len(items),
                "n_segments_in_formula": items[0]["n_segments"],
                "n_breakpoints": items[0]["n_breakpoints"],
                "hxb2_barcode": items[0]["hxb2_barcode"],
                "median_mapped_coding_nt": m_nts[len(m_nts) // 2],
                "min_fisher_p": ps[0] if ps else None,
                "median_fisher_p": ps[len(ps) // 2] if ps else None,
                "subtypes": subtypes,
                "years": years,
                "top_countries": countries,
                "top_hosts": hosts,
                "example_strains": [x["strain"] for x in items[:5] if x["strain"]],
                "example_biosamples": [x["biosample"] for x in items[:5] if x["biosample"]],
                "example_ha_accs": [x["ha_acc"] for x in items[:5] if x["ha_acc"]],
            }
            novel_constellations.append(cl_rec)

            # Save up to 10 best-covered representatives to iav_h5_novel_hits.fasta.gz
            reps = sorted(items, key=lambda x: (x["n_bases"], -x["mapped_coding_nt"], x["uid"]))[:10]
            for rep in reps:
                iso = iso_by_uid[rep["uid"]]
                st_clean = rep["strain"].replace(" ", "_")
                ct_clean = rep["country"].replace(" ", "_")
                hst_clean = rep["host"].replace(" ", "_")
                gz.write(
                    f">{rep['uid']} biosample={rep['biosample'] or 'NA'} strain={st_clean} "
                    f"subtype={rep['subtype']} year={rep['year']} country={ct_clean} host={hst_clean} "
                    f"mapped_nt={rep['mapped_coding_nt']} formula={form}\n{iso['concat_seq']}\n"
                )
                saved_iav_reps += 1

    # Write results/bvbrc_iav_h5_novel_constellations.tsv
    iav_tsv_cols = [
        "rank",
        "n_isolates",
        "structural_formula",
        "epidemiological_formula",
        "backbone_trunk",
        "constituent_trunks",
        "n_breakpoints",
        "median_mapped_coding_nt",
        "min_fisher_p",
        "median_fisher_p",
        "subtypes",
        "years",
        "top_countries",
        "top_hosts",
        "example_strains",
        "example_biosamples",
        "example_ha_accs",
        "hxb2_barcode",
    ]
    with open(OUT_IAV_TSV, "w", newline="", encoding="utf-8") as tf:
        w = csv.DictWriter(tf, fieldnames=iav_tsv_cols, delimiter="\t")
        w.writeheader()
        for rank, cl in enumerate(novel_constellations, start=1):
            w.writerow({
                "rank": rank,
                "n_isolates": cl["n_isolates"],
                "structural_formula": cl["structural_formula"],
                "epidemiological_formula": cl["epidemiological_formula"],
                "backbone_trunk": cl["backbone_trunk"],
                "constituent_trunks": ",".join(cl["constituent_trunks"]),
                "n_breakpoints": cl["n_breakpoints"],
                "median_mapped_coding_nt": cl["median_mapped_coding_nt"],
                "min_fisher_p": cl["min_fisher_p"],
                "median_fisher_p": cl["median_fisher_p"],
                "subtypes": ";".join(f"{k}:{v}" for k, v in cl["subtypes"]),
                "years": ",".join(cl["years"]),
                "top_countries": ";".join(f"{k}:{v}" for k, v in cl["top_countries"]),
                "top_hosts": ";".join(f"{k}:{v}" for k, v in cl["top_hosts"]),
                "example_strains": ";".join(cl["example_strains"]),
                "example_biosamples": ";".join(cl["example_biosamples"]),
                "example_ha_accs": ";".join(cl["example_ha_accs"]),
                "hxb2_barcode": cl["hxb2_barcode"],
            })

    # Clean up concatenated sequences from memory
    for iso in assembled_isolates:
        iso.pop("concat_seq", None)

    years_summary = {}
    for yr_s, cnt in sorted(hc_years.items()):
        top_yr_calls = sorted(
            [(k[5:], v) for k, v in cnt.items() if k.startswith("call:")],
            key=lambda x: -x[1],
        )[:12]
        years_summary[yr_s] = {
            "total": cnt["total"],
            "categories": {k[4:]: v for k, v in cnt.items() if k.startswith("cat:")},
            "top_calls": top_yr_calls,
        }

    total_dt = round(time.time() - t_mod_start, 2)
    out = {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "segment_rows_expected": expected_total,
        "segment_rows_fetched": len(rows),
        "assembly_stats": {
            "skipped_empty_biosample_and_strain": skipped_empty_id,
            "invalid_segment_rows": invalid_seg,
            "unique_isolate_groups": len(groups),
            "grouped_by_biosample": grouped_by_biosample,
            "grouped_by_strain": grouped_by_strain,
            "complete_8seg_isolates": len(complete_isolates),
            "duplicate_segment_slots_resolved": duplicate_slot_count,
            "fasta_segments_requested": len(needed_gids),
            "fasta_segments_joined": len(needed_gids) - missing_seg_isolates * 8,
            "fasta_join_rate_expected": round(join_rate_expected, 5),
            "fasta_join_rate_headers": round(join_rate_headers, 5),
            "assembled_8seg_screened": len(egg_results),
            "high_coverage_mapped_ge_11500": n_high_cov,
        },
        "timings_sec": {
            "metadata_fetch_sec": dt_meta,
            "fasta_stream_download_sec": dt_fetch,
            "eggman_sec": round(dt_eggman, 3),
            "throughput_iso_per_sec": round(len(egg_results) / dt_eggman, 1) if dt_eggman > 0 else 0.0,
            "total_module_c_sec": total_dt,
        },
        "all_assembled_summary": {
            "categories": dict(all_cats),
            "calls": all_calls.most_common(30),
        },
        "high_coverage_summary": {
            "n_isolates": n_high_cov,
            "categories": dict(hc_cats),
            "calls": hc_calls.most_common(30),
            "subtypes": hc_subtypes.most_common(20),
            "top_countries": hc_countries.most_common(25),
            "top_hosts": hc_hosts.most_common(30),
            "by_year": years_summary,
        },
        "novel_urf_constellations_ge_3": {
            "n_constellations": len(novel_constellations),
            "total_isolates_in_constellations": sum(c["n_isolates"] for c in novel_constellations),
            "n_total_high_cov_urf_isolates": hc_cats.get("URF", 0),
            "n_distinct_high_cov_urf_formulas": len(hc_urf_clusters),
            "representatives_saved_to_fasta": saved_iav_reps,
            "constellations": novel_constellations,
        },
        "outputs": {
            "summary_json": str(OUT_IAV_JSON),
            "novel_constellations_tsv": str(OUT_IAV_TSV),
            "novel_hits_fasta_gz": str(OUT_IAV_FASTA_GZ),
        },
    }
    with open(OUT_IAV_JSON, "w", encoding="utf-8") as f:
        json.dump(out, f, indent=2)

    print(
        f"[+] Module C complete in {total_dt:.2f}s | complete_8seg={len(egg_results):,} | "
        f"high_cov(>=11500)={n_high_cov:,} | novel_URF_constellations(>=3)={len(novel_constellations)}",
        flush=True,
    )
    return out


def main():
    parser = argparse.ArgumentParser(description="BV-BRC live census & 100k SC2 + 8-seg H5Nx EGGMAN streaming screen")
    parser.add_argument("--module", choices=["all", "A", "B", "C"], default="all")
    parser.add_argument("--cache-dir", default=DEFAULT_CACHE_DIR, help="Local cache dir for sc2_100k_stream & iav_bvbrc_3980")
    parser.add_argument("--reuse-census", action="store_true", help="Reuse cached bvbrc_archive_census.json in Module A if present")
    args = parser.parse_args()

    t0 = time.time()
    if args.module in ("all", "A"):
        if args.reuse_census and OUT_CENSUS.exists():
            print(f"[+] Reusing cached BV-BRC archive census from {OUT_CENSUS}", flush=True)
        else:
            run_module_a()
    if args.module in ("all", "B"):
        run_module_b(cache_dir=args.cache_dir)
    if args.module in ("all", "C"):
        run_module_c()
    print(f"\n[DONE] Total elapsed: {time.time() - t0:.2f}s", flush=True)


if __name__ == "__main__":
    main()
