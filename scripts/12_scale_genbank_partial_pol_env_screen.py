#!/usr/bin/env python3
"""
Script 12: Large-Scale GenBank + LANL Partial `pol` (2000..4500 nt) and
Partial `env` (2000..3500 nt) Harvesting, Patient-Clone Deduplication,
Profile-Locked Alignment, and RhizAeon Circulating Recombinant Discovery.

Scans:
  - All 21,257 intact GenBank near-full `pol` sequences (2,000..4,500 nt)
  - All 83,794 intact GenBank near-full `env` sequences (2,000..3,500 nt)
  - Combined with the 2,385 strictly partial-only `pol` and 4,723 strictly partial-only
    `env` sequences from the LANL curated alignments (`Script 11`).
"""

import json
import os
import re
import subprocess
import sys
import tempfile
import time
import urllib.parse
import urllib.request
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
s08b = import_module("08b_fetch_genbank_cluster_metadata")
s09 = import_module("09_generate_segment_phylogenies")
s10 = import_module("10_recombination_viability_and_establishment_sieve")
s11 = import_module("11_screen_partial_pol_and_env")

PROJECT_ROOT = Path("/Users/sergei/Projects/HIV-recombination")
DATA_DIR = PROJECT_ROOT / "data"
RESULTS_DIR = PROJECT_ROOT / "results"
TREES_DIR = RESULTS_DIR / "trees"
EGGMAN_BIN = PROJECT_ROOT / "rhiztype" / "target" / "release" / "eggman"
MAFFT_BIN = "/usr/local/bin/mafft"

COUNTRY_MAP = {
    "united states": "US", "usa": "US", "china": "CN", "thailand": "TH", "south africa": "ZA",
    "uganda": "UG", "kenya": "KE", "tanzania": "TZ", "cameroon": "CM", "brazil": "BR",
    "russia": "RU", "russian federation": "RU", "spain": "ES", "united kingdom": "GB",
    "france": "FR", "germany": "DE", "italy": "IT", "canada": "CA", "india": "IN",
    "malaysia": "MY", "viet nam": "VN", "vietnam": "VN", "myanmar": "MM", "japan": "JP",
    "south korea": "KR", "korea": "KR", "australia": "AU", "sweden": "SE", "belgium": "BE",
    "switzerland": "CH", "portugal": "PT", "argentina": "AR", "peru": "PE", "chile": "CL",
    "mexico": "MX", "cuba": "CU", "nigeria": "NG", "senegal": "SN", "ghana": "GH",
    "democratic republic of the congo": "CD", "drc": "CD", "zambia": "ZM", "malawi": "MW",
    "botswana": "BW", "zimbabwe": "ZW", "rwanda": "RW", "ethiopia": "ET", "cyprus": "CY",
    "greece": "GR", "poland": "PL", "ukraine": "UA", "belarus": "BY", "kazakhstan": "KZ",
    "uzbekistan": "UZ", "pakistan": "PK", "indonesia": "ID", "philippines": "PH", "taiwan": "TW",
    "hong kong": "HK", "singapore": "SG", "denmark": "DK", "netherlands": "NL", "austria": "AT",
    "finland": "FI", "luxembourg": "LU", "romania": "RO", "angola": "AO", "equatorial guinea": "GQ",
}


def parse_genbank_fasta_header(hdr_line: str):
    """
    Extract (accession, patient_key, country, year, title) from a raw GenBank FASTA header:
      >PZ982380.1 HIV-1 isolate IP1_Oct16_PR-IN_1_3s from USA pol protein (pol) gene, partial cds
    """
    raw = hdr_line.lstrip(">").strip()
    parts = raw.split(None, 1)
    acc_ver = parts[0]
    acc = acc_ver.split(".")[0].strip()
    title = parts[1] if len(parts) > 1 else ""
    tl = title.lower()

    # Extract country if "from <Country>" is in title
    country = "Unknown"
    m_cty = re.search(r"\bfrom\s+([A-Za-z\s]+?)(?:\s+(?:gag|pol|vif|vpr|tat|rev|vpu|env|nef|complete|partial|clone|isolate|strain|in\s+\d{4}|,|$))", title, re.IGNORECASE)
    if m_cty:
        c_raw = m_cty.group(1).strip().lower()
        country = COUNTRY_MAP.get(c_raw, c_raw[:12].title())

    # Extract year if present
    year = np.nan
    m_yr = re.search(r"\b(19[789]\d|20[012]\d)\b", title)
    if m_yr:
        year = int(m_yr.group(1))

    # Extract isolate / clone / strain identifier
    iso = ""
    m_iso = re.search(r"\b(?:isolate|strain|patient|subject|participant)\s+([^\s,;]+)", title, re.IGNORECASE)
    if m_iso:
        iso = m_iso.group(1).strip()

    # Normalize isolate to patient root so multi-clone SGA dumps from 1 patient collapse cleanly:
    # e.g. V704_1716_090_GP_64_1 -> V704_1716; IP1_Oct16_PR-IN_1_3s -> IP1; CAP256_W30_c12 -> CAP256
    acc_prefix = acc[:4]  # Study submission block prefix
    if iso:
        clean_iso = re.sub(r"^(?:HIV-?1[-_]?|HIV[-_]?)", "", iso, flags=re.IGNORECASE)
        tokens = [t for t in re.split(r"[_.\-/]+", clean_iso) if t]
        if len(tokens) >= 2 and len(tokens[0]) <= 5 and any(ch.isdigit() for ch in tokens[1]):
            pat_root = f"{tokens[0]}_{tokens[1]}"
        elif len(tokens) >= 1 and len(tokens[0]) >= 3:
            pat_root = tokens[0]
        else:
            pat_root = clean_iso[:8]
        # Further strip trailing clone/timepoint markers like w12, d14, m6, c1, cl2, p1, sga3
        pat_root = re.sub(r"[-_](?:w|wk|week|d|day|m|mo|yr|v|vis|c|cl|clone|sga|seq|env|pol|gp\d+|p\d+)\d*$", "", pat_root, flags=re.IGNORECASE)
        pat_key = f"{acc_prefix}:{pat_root.upper()}"
    else:
        # No isolate token in title -> group by 5-char accession block + country + first 35 chars of title
        pat_key = f"{acc[:5]}:{country}:{re.sub(r'[^a-z0-9]', '', tl)[:24]}"

    return acc, pat_key, country, year, title


def harvest_and_dedup_genbank_region(
    region_name: str,
    query_str: str,
    min_valid_len: int,
    exclude_accs: set,
    out_fasta: Path,
    out_meta_csv: Path,
):
    if out_fasta.exists() and out_meta_csv.exists() and out_fasta.stat().st_size > 100_000:
        meta_df = pd.read_csv(out_meta_csv)
        print(f"[+] Using cached patient-deduplicated GenBank {region_name} dataset: {len(meta_df):,} unique patient isolates ({out_fasta.name})")
        return meta_df

    print("\n" + "=" * 88)
    print(f"HARVESTING & PATIENT-DEDUPLICATING GENBANK {region_name.upper()} AT SCALE")
    print("=" * 88)

    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?" + urllib.parse.urlencode({
        "db": "nuccore", "term": query_str, "retmode": "json", "usehistory": "y", "retmax": 0
    })
    with urllib.request.urlopen(url, timeout=30) as r:
        res = json.loads(r.read().decode())["esearchresult"]
        total_count = int(res["count"])
        webenv = res["webenv"]
        qk = res["querykey"]

    print(f"NCBI Entrez matched {total_count:,} intact {region_name} sequences. Streaming in batches of 5,000...")

    best_by_pat_key = {}  # pat_key -> (valid_len, acc, country, year, title, seq)
    n_streamed = 0
    n_already_in_full_or_lanl = 0
    n_filtered_bad = 0

    batch_size = 5000
    for retstart in range(0, total_count, batch_size):
        f_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?" + urllib.parse.urlencode({
            "db": "nuccore",
            "query_key": qk,
            "WebEnv": webenv,
            "rettype": "fasta",
            "retmode": "text",
            "retstart": retstart,
            "retmax": batch_size,
        })
        fasta_chunk = ""
        for attempt in range(4):
            try:
                with urllib.request.urlopen(f_url, timeout=90) as resp:
                    fasta_chunk = resp.read().decode("utf-8", errors="ignore")
                break
            except Exception as e:
                print(f"  [!] Batch {retstart} attempt {attempt+1} error: {e}; retrying in {2*(attempt+1)}s...")
                time.sleep(2 * (attempt + 1))

        if not fasta_chunk:
            continue

        cur_hdr = None
        cur_seq = []
        for line in fasta_chunk.splitlines():
            if line.startswith(">"):
                if cur_hdr is not None:
                    n_streamed += 1
                    seq_s = "".join(cur_seq).upper()
                    acc, pk, cty, yr, ttl = parse_genbank_fasta_header(cur_hdr)
                    tl = ttl.lower()
                    if acc in exclude_accs:
                        n_already_in_full_or_lanl += 1
                    elif any(w in tl for w in ("nonfunctional", "defective", "hypermutated", "unverified", "protein-like", "pseudogene", "synthetic", "vector", "construct", "patent", "group o", "group n", "group p", "siv")):
                        n_filtered_bad += 1
                    else:
                        acgt_len = sum(seq_s.count(b) for b in "ACGT")
                        if acgt_len >= min_valid_len:
                            prev = best_by_pat_key.get(pk)
                            if prev is None or acgt_len > prev[0]:
                                best_by_pat_key[pk] = (acgt_len, acc, cty, yr, ttl, seq_s)
                cur_hdr = line.strip()
                cur_seq = []
            else:
                cur_seq.append(line.strip())

        if cur_hdr is not None:
            n_streamed += 1
            seq_s = "".join(cur_seq).upper()
            acc, pk, cty, yr, ttl = parse_genbank_fasta_header(cur_hdr)
            tl = ttl.lower()
            if acc in exclude_accs:
                n_already_in_full_or_lanl += 1
            elif any(w in tl for w in ("nonfunctional", "defective", "hypermutated", "unverified", "protein-like", "pseudogene", "synthetic", "vector", "construct", "patent", "group o", "group n", "group p", "siv")):
                n_filtered_bad += 1
            else:
                acgt_len = sum(seq_s.count(b) for b in "ACGT")
                if acgt_len >= min_valid_len:
                    prev = best_by_pat_key.get(pk)
                    if prev is None or acgt_len > prev[0]:
                        best_by_pat_key[pk] = (acgt_len, acc, cty, yr, ttl, seq_s)

        print(f"  Streamed {min(retstart + batch_size, total_count):6,d}/{total_count:,d} | "
              f"Excl (Full/LANL): {n_already_in_full_or_lanl:5,d} | "
              f"Unique New Patient Keys: {len(best_by_pat_key):5,d}")
        time.sleep(0.35)

    meta_rows = []
    with open(out_fasta, "w") as f_out:
        for pk, (acgt_len, acc, cty, yr, ttl, seq_s) in best_by_pat_key.items():
            f_out.write(f">{acc}\n{seq_s}\n")
            meta_rows.append({
                "accession": acc,
                "pat_key": pk,
                "country": cty,
                "year": yr,
                "acgt_len": acgt_len,
                "title": ttl,
            })

    meta_df = pd.DataFrame(meta_rows)
    meta_df.to_csv(out_meta_csv, index=False)
    print(f"[+] Saved {len(meta_df):,} new patient-deduplicated GenBank {region_name} sequences to {out_fasta.name}")
    return meta_df


# Worker globals for alignment + RhizAeon window screening
G_WIN_PROFILE_PATH = None
G_TRUNK_MAT_WIN = None
G_TRUNK_NAMES = None
G_HXB2_COORDS_WIN = None
G_CM240_WIN = None
G_MASK_BOOL = None
G_WIN_START = None
G_WIN_END = None
G_MIN_VALID_NT = None


def build_window_profile_fasta(full_profile_path: Path, win_start: int, win_end: int, out_profile_path: Path):
    """Slice the 15-trunk 9,719-nt HXB2 profile to [win_start, win_end] (1-based HXB2)."""
    recs = s3.read_fasta(full_profile_path)
    with open(out_profile_path, "w") as f:
        for hdr, seq in recs:
            sub_seq = seq[win_start - 1 : win_end]
            f.write(f">{hdr}\n{sub_seq}\n")
    return out_profile_path


def screen_partial_pairs_with_eggman(raw_pairs, w_start: int, w_end: int, min_aln_nt: int, chunk_size: int = 4000):
    """Screen partial pol/env sequences via rebuilt eggman binary (--gamma 10.5 --min-snps 7 --min-span 200 --p-thresh 1e-5)."""
    assert EGGMAN_BIN.exists(), f"Missing eggman binary: {EGGMAN_BIN}"
    out = {}
    for idx in range(0, len(raw_pairs), chunk_size):
        chunk = raw_pairs[idx : idx + chunk_size]
        with tempfile.TemporaryDirectory() as td:
            in_fa = Path(td) / "chunk.fasta"
            out_js = Path(td) / "chunk.json"
            with open(in_fa, "w") as f:
                for acc, seq in chunk:
                    clean_s = str(seq).upper().replace("-", "").replace(".", "")
                    if len(clean_s) >= 200:
                        f.write(f">{acc}\n{clean_s}\n")
            subprocess.run(
                [
                    str(EGGMAN_BIN),
                    str(in_fa),
                    "--gamma", "10.5",
                    "--min-snps", "7",
                    "--min-span", "200",
                    "--p-thresh", "1e-5",
                    "--json", str(out_js),
                    "-q",
                ],
                check=True,
            )
            recs = json.loads(out_js.read_text())
        for r in recs:
            cat = str(r.get("called_category") or "")
            if cat in ("Insufficient_Data", "Non_GroupM_Outlier"):
                continue
            s0 = int(r.get("hxb2_span_start") or 0)
            s1 = int(r.get("hxb2_span_end") or 0)
            ov = min(s1, w_end) - max(s0, w_start) + 1
            if ov < min_aln_nt or int(r.get("mapped_coding_nt") or 0) < min_aln_nt:
                continue
            segs = r.get("segments") or []
            if not segs:
                continue
            min_id = min(float(s.get("identity_pct") or 0.0) for s in segs)
            if min_id < 83.0:
                continue
            max_div = round(max(0.0, 1.0 - min_id / 100.0), 4)
            bc = str(r.get("basal_barcode") or r.get("hxb2_barcode") or "")
            formula, norm_bc = s11.project_barcode_to_window(bc, w_start, w_end, min_seg_nt=180)
            if not norm_bc:
                continue
            out[str(r["sequence_id"])] = {
                "valid_nt": int(r["mapped_coding_nt"]),
                "max_div": max_div,
                "egg_category": cat,
                "egg_lineage": str(r.get("called_lineage") or ""),
                "window_formula": formula,
                "window_barcode": norm_bc,
                "basal_formula": str(r.get("basal_formula") or ""),
            }
    return out


def align_partial_candidates_to_window(acc_seq_pairs, win_profile_path: Path, expected_len: int, mask_bool: np.ndarray):
    """Align selected candidate partial sequences to the [w_start, w_end] 15-trunk profile via mafft --keeplength."""
    if not acc_seq_pairs:
        return {}
    out = {}
    for start in range(0, len(acc_seq_pairs), 400):
        batch = acc_seq_pairs[start : start + 400]
        with tempfile.NamedTemporaryFile("w", suffix=".fasta", delete=False) as tf:
            frag_path = tf.name
            for acc, seq in batch:
                clean_s = str(seq).upper().replace("-", "").replace(".", "")
                tf.write(f">{acc}\n{clean_s}\n")
        try:
            proc = subprocess.run(
                [MAFFT_BIN, "--quiet", "--6merpair", "--keeplength", "--thread", "4", "--addfragments", frag_path, str(win_profile_path)],
                capture_output=True, text=True, check=True,
            )
            aln_text = proc.stdout
        finally:
            if os.path.exists(frag_path):
                os.remove(frag_path)
        cur_id, cur_s = None, []
        for line in aln_text.splitlines():
            if line.startswith(">"):
                if cur_id is not None and not cur_id.startswith("Profile_"):
                    arr = s3.seq_to_int_array("".join(cur_s).upper())
                    if len(arr) == expected_len:
                        arr[mask_bool] = -1
                        out[cur_id] = arr
                cur_id = line[1:].strip().split()[0]
                cur_s = []
            else:
                cur_s.append(line.strip())
        if cur_id is not None and not cur_id.startswith("Profile_"):
            arr = s3.seq_to_int_array("".join(cur_s).upper())
            if len(arr) == expected_len:
                arr[mask_bool] = -1
                out[cur_id] = arr
    return out


def main():
    t0_all = time.time()
    TREES_DIR.mkdir(parents=True, exist_ok=True)
    full_profile_path, trunk_mat_full, trunk_names, hxb2_coords_full, cm240_full = s8.build_15trunk_hxb2_profile()
    selected_refs, crf_g0_tail, a1_ae_refs, a1_ag_refs, ref_hxb2_coords = s09.load_reference_panel_coding()

    gb_full = pd.read_csv(RESULTS_DIR / "genbank_31k_rhizaeon_classifications.csv")
    lanl_full = pd.read_csv(RESULTS_DIR / "lanl_all_full_genomes_rhizaeon_validation.csv")
    lanl_pol = pd.read_csv(RESULTS_DIR / "partial_pol_rhizaeon_screen.csv")
    lanl_env = pd.read_csv(RESULTS_DIR / "partial_env_rhizaeon_screen.csv")

    full_genome_accs = set(gb_full["accession"].astype(str)) | set(
        lanl_full["sequence_id"].apply(lambda x: str(x).split(".")[-1])
    )
    exclude_pol_accs = full_genome_accs | set(lanl_pol["accession"].astype(str))
    exclude_env_accs = full_genome_accs | set(lanl_env["accession"].astype(str))

    comp_df = pd.read_csv(RESULTS_DIR / "hiv1_descriptive_nomenclature_compendium.csv")
    novel_df = pd.read_csv(RESULTS_DIR / "genbank_novel_circulating_recombinants.csv")
    reg_mosaic_df = novel_df[novel_df["phylogenetic_classification"] == "Candidate_Regional_Mosaic"].copy() if "phylogenetic_classification" in novel_df.columns else pd.DataFrame()

    # Align validated Candidate_Regional_Mosaic (RU_G0B_gagMA) genomes onto coding 790..9417 for partial-carrier checks
    all_reg_accs = set()
    for m_str in reg_mosaic_df.get("member_accessions", pd.Series([], dtype=str)):
        for a in str(m_str).split(","):
            if a.strip():
                all_reg_accs.add(a.strip())
    reg_coding_map = {}
    if all_reg_accs:
        raw_pairs = []
        with open(s8.GENBANK_FASTA, "r") as f:
            cur_acc, cur_seq = None, []
            for line in f:
                if line.startswith(">"):
                    if cur_acc in all_reg_accs:
                        raw_pairs.append((cur_acc, "".join(cur_seq).upper()))
                    cur_acc = line[1:].strip().split()[0].split(".")[0]
                    cur_seq = []
                elif cur_acc in all_reg_accs:
                    cur_seq.append(line.strip())
            if cur_acc in all_reg_accs:
                raw_pairs.append((cur_acc, "".join(cur_seq).upper()))
        loop_mask_full = np.zeros(len(hxb2_coords_full), dtype=bool)
        for h_lo, h_hi in s8.HYPERVARIABLE_MASK_RANGES:
            loop_mask_full |= (hxb2_coords_full >= h_lo) & (hxb2_coords_full <= h_hi)
        reg_coding_map = s8.align_sequences_to_coding_arrays(raw_pairs, full_profile_path, loop_mask_full)

    not_clause = (
        "NOT (nonfunctional[Title] OR defective[Title] OR hypermutated[Title] OR "
        "unverified[Title] OR synthetic[Title] OR patent[Title] OR vector[Title] OR "
        "construct[Title] OR \"complete genome\"[Title] OR \"protein-like\"[Title] OR pseudogene[Title])"
    )
    q_pol = f"txid11676[Organism:exp] AND 2000:4500[SLEN] AND (pol[Gene] OR pol[Title] OR protease[Title] OR \"reverse transcriptase\"[Title] OR integrase[Title]) {not_clause}"
    q_env = f"txid11676[Organism:exp] AND 2000:3500[SLEN] AND (env[Gene] OR env[Title] OR envelope[Title] OR gp160[Title]) {not_clause}"

    configs = [
        ("pol", 2085, 5096, q_pol, 1850, exclude_pol_accs, lanl_pol, 1500),
        ("env", 6225, 8795, q_env, 1850, exclude_env_accs, lanl_env, 1300),
    ]

    for region_name, w_start, w_end, q_str, min_raw_len, excl_accs, lanl_win_df, min_aln_nt in configs:
        t_reg = time.time()
        gb_fasta = DATA_DIR / f"genbank_partial_{region_name}_dedup.fasta"
        gb_meta_csv = DATA_DIR / f"genbank_partial_{region_name}_dedup_metadata.csv"

        meta_df = harvest_and_dedup_genbank_region(
            region_name, q_str, min_raw_len, excl_accs, gb_fasta, gb_meta_csv
        )
        meta_by_acc = {str(r["accession"]): r for _, r in meta_df.iterrows()}

        win_profile_path = DATA_DIR / f"hxb2_15trunk_{region_name}_{w_start}_{w_end}.fasta"
        build_window_profile_fasta(full_profile_path, w_start, w_end, win_profile_path)

        idx0 = w_start - 790
        idx1 = w_end - 790 + 1
        hxb2_coords_win = np.arange(w_start, w_end + 1, dtype=np.int32)
        expected_len = w_end - w_start + 1

        mask_bool = np.zeros(expected_len, dtype=bool)
        if region_name == "env":
            for s_l, e_l in s11.ENV_LOOP_RANGES:
                mask_bool |= ((hxb2_coords_win >= s_l) & (hxb2_coords_win <= e_l))

        # 1. Screen cached GenBank partial FASTA via rebuilt eggman binary
        raw_gb_recs = s3.read_fasta(gb_fasta)
        raw_gb_by_acc = {acc: seq.upper().replace("-", "").replace(".", "") for acc, seq in raw_gb_recs}
        print(f"\n[*] Screening {len(raw_gb_recs):,} GenBank partial {region_name} sequences via {EGGMAN_BIN.name} --gamma 10.5 ...")
        gb_egg = screen_partial_pairs_with_eggman(raw_gb_recs, w_start, w_end, min_aln_nt)
        print(f"    -> {len(gb_egg):,} validated GenBank partial {region_name} isolates ({time.time() - t_reg:.1f}s)")

        # 2. Screen strictly partial-only LANL sequences via rebuilt eggman binary
        lanl_fasta = DATA_DIR / f"lanl_all_2022_{region_name}.fasta"
        lanl_recs = s3.read_fasta(lanl_fasta)
        hxb2_rec = [seq for hdr, seq in lanl_recs if "K03455" in hdr][0]
        hxb2_cols = s8.get_hxb2_all_cols(hxb2_rec)
        lanl_meta_by_acc = {str(r["accession"]): r for _, r in lanl_win_df[~lanl_win_df["in_full_genome_set"]].iterrows()}

        seq_arrays = {}
        lanl_raw_pairs = []
        seen_lanl = set()
        for hdr, seq in lanl_recs:
            acc = hdr.split(".")[-1].strip()
            if acc in full_genome_accs or acc not in lanl_meta_by_acc or acc in seen_lanl:
                continue
            seen_lanl.add(acc)
            arr = s3.seq_to_int_array(seq.upper())[hxb2_cols]
            if len(arr) == expected_len:
                arr[mask_bool] = -1
                seq_arrays[acc] = arr
            lanl_raw_pairs.append((acc, seq.upper().replace("-", "").replace(".", "")))

        print(f"[*] Screening {len(lanl_raw_pairs):,} strictly partial LANL {region_name} sequences via {EGGMAN_BIN.name} --gamma 10.5 ...")
        lanl_egg = screen_partial_pairs_with_eggman(lanl_raw_pairs, w_start, w_end, min_aln_nt)
        print(f"    -> {len(lanl_egg):,} validated LANL partial {region_name} isolates")

        # Build unified table of all validated partial isolates
        gb_screen_rows = []
        for acc, res in gb_egg.items():
            m_row = meta_by_acc.get(acc, {})
            gb_screen_rows.append({
                "source": "GenBank_Scale_Partial",
                "region": region_name,
                "accession": acc,
                "country": m_row.get("country", "Unknown"),
                "year": m_row.get("year", np.nan),
                "pat_key": m_row.get("pat_key", acc),
                "title": m_row.get("title", ""),
                "in_full_genome_set": False,
                "valid_nt": res["valid_nt"],
                "egg_category": res["egg_category"],
                "egg_lineage": res["egg_lineage"],
                "window_formula": res["window_formula"],
                "window_barcode": res["window_barcode"],
            })
        for acc, res in lanl_egg.items():
            lr = lanl_meta_by_acc[acc]
            gb_screen_rows.append({
                "source": "LANL_Partial_Only",
                "region": region_name,
                "accession": acc,
                "country": str(lr["country"]),
                "year": lr["year"],
                "pat_key": f"LANL:{acc}",
                "title": str(lr["header"]),
                "in_full_genome_set": False,
                "valid_nt": res["valid_nt"],
                "egg_category": res["egg_category"],
                "egg_lineage": res["egg_lineage"],
                "window_formula": res["window_formula"],
                "window_barcode": res["window_barcode"],
            })

        # Build Known CRF window lookup map (2023 LANL + 8 post-2023 promoted CRFs from compendium_bundle.json)
        bundle_path = PROJECT_ROOT / "rhiztype" / "crates" / "eggman-core" / "data" / "compendium_bundle.json"
        bundle_entries = json.loads(bundle_path.read_text())["compendium"] if bundle_path.exists() else []
        comp_lbls_seen = set(comp_df["legacy_label"].astype(str))
        known_ref_pairs = [(str(r["legacy_label"]), str(r["hxb2_barcode"])) for _, r in comp_df.iterrows()]
        for be in bundle_entries:
            lbl = str(be.get("label", ""))
            if lbl not in comp_lbls_seen and be.get("status") == "Verified_Recombinant_Lineage" and be.get("hxb2_barcode"):
                known_ref_pairs.append((lbl, str(be["hxb2_barcode"])))

        known_win_map = {"_SETS_": []}
        for k_lbl, k_bc in known_ref_pairs:
            wf, wbc = s11.project_barcode_to_window(k_bc, w_start, w_end, min_seg_nt=180)
            if "[" in wf:
                if wf not in known_win_map:
                    known_win_map[wf] = k_lbl
                else:
                    known_win_map[wf] += f"/{k_lbl}"
                k_segs = s10.parse_barcode_map(wbc)
                k_dom = wf.split("[")[0]
                k_donors = {(trk, s3.hxb2_to_short_gene(s_c, e_c)) for s_c, e_c, trk in k_segs if trk != k_dom}
                known_win_map["_SETS_"].append((k_lbl, k_dom, k_donors))

        reg_win_info = []
        for _, r in reg_mosaic_df.iterrows():
            wf, wbc = s11.project_barcode_to_window(str(r["representative_hxb2_barcode"]), w_start, w_end, min_seg_nt=180)
            w_segs = s10.parse_barcode_map(wbc)
            w_bps = [bp for bp, _, _ in s10.get_internal_breakpoints(w_segs)]
            mem_accs = [x.strip() for x in str(r["member_accessions"]).split(",") if x.strip()]
            rep_arrs = [reg_coding_map[a][idx0:idx1] for a in mem_accs if a in reg_coding_map]
            reg_win_info.append({
                "label": "RU_G0B_gagMA",
                "full_formula": str(r["rhizaeon_formula"]),
                "win_formula": wf,
                "win_bps": w_bps,
                "rep_arrays": rep_arrs,
            })

        # First annotation pass (eggman + known window formulas + RU_G0B_gagMA direct calls)
        for r in gb_screen_rows:
            segs = s10.parse_barcode_map(r["window_barcode"])
            formula = r["window_formula"]
            dom_trunk = formula.split("[")[0] if "[" in formula else formula
            matched_known = ""
            if r["egg_category"] in ("Certified_CRF", "Ancestral_Recombinant_Subtype"):
                matched_known = f"Eggman:{r['egg_lineage']}"
            elif r["egg_category"] == "Pure_Subtype":
                matched_known = f"Clonal:{r['egg_lineage']}"
            elif formula in known_win_map:
                matched_known = known_win_map[formula]
            elif "[" in formula:
                q_donors = {(trk, s3.hxb2_to_short_gene(s_c, e_c)) for s_c, e_c, trk in segs if trk != dom_trunk}
                for k_lbl, k_dom, k_donors in known_win_map["_SETS_"]:
                    if k_dom == dom_trunk and q_donors.issubset(k_donors) and len(q_donors) > 0:
                        matched_known = f"Attenuated:{k_lbl}"
                        break

            matched_reg = r["egg_lineage"] if r["egg_lineage"] == "RU_G0B_gagMA" else ""
            r["dominant_trunk"] = dom_trunk
            r["n_segments"] = len(segs)
            r["n_breakpoints"] = len(segs) - 1
            r["matched_known_crf_window"] = matched_known
            r["matched_novel15_lineage"] = matched_reg
            r["novel15_match_mode"] = "Eggman_Tier2_Certified" if matched_reg else ""
            r["novel15_genetic_div"] = np.nan

        comb_df = pd.DataFrame(gb_screen_rows)

        # Surface all candidate multi-sequence URF formulas (egg_category == "URF", not matched_known, not matched_reg)
        unid_pre = comb_df[
            (comb_df["egg_category"] == "URF")
            & (comb_df["n_segments"] > 1)
            & (comb_df["matched_known_crf_window"] == "")
            & (comb_df["matched_novel15_lineage"] == "")
        ]
        vc_unid = unid_pre["window_formula"].value_counts()
        cand_formulas = set(vc_unid[vc_unid >= 2].index)
        cand_accs = sorted(set(unid_pre[unid_pre["window_formula"].isin(cand_formulas)]["accession"].astype(str)))

        # Align multi-patient URF candidates + RU_G0B_gagMA formula carriers + fetch 08b flatfiles!
        gb_cand_pairs = [(a, raw_gb_by_acc[a]) for a in cand_accs if a in raw_gb_by_acc and a not in seq_arrays]
        if gb_cand_pairs:
            print(f"[+] Aligning {len(gb_cand_pairs)} candidate partial {region_name} URF sequences onto HXB2 window profile...")
            seq_arrays.update(align_partial_candidates_to_window(gb_cand_pairs, win_profile_path, expected_len, mask_bool))

        if cand_accs:
            print(f"[+] Fetching/checking 08b GenBank flatfiles for {len(cand_accs)} candidate partial {region_name} URF accessions...")
            s08b.ensure_flatfiles_cached(cand_accs)

        ff_df = pd.read_csv(DATA_DIR / "genbank_cluster_flatfiles.tsv", sep="\t", dtype=str).fillna("")
        ff_map = {str(r["accession"]): r.to_dict() for _, r in ff_df.iterrows()}

        # Check flatfile titles/official_crf_note + strict patient root + RU_G0B_gagMA carrier alignment
        def strict_patient_root(row):
            acc = str(row["accession"])
            ttl = str(row["title"])
            src = str(row["source"])
            cty = str(row["country"])
            if src == "LANL_Partial_Only":
                parts = ttl.split(".")
                strain = parts[3] if len(parts) >= 5 else acc
                root = re.split(r"[_.\-#/]+", strain)[0].upper()
                root = re.sub(r"(?:V\d+|W\d+|D\d+|M\d+)$", "", root)
                if len(root) < 2:
                    root = strain[:6].upper()
                return f"{cty}:{root}"
            m = re.search(r"\b(?:isolate|strain|clone|patient|subject|participant)\s+([^\s,;]+)", ttl, re.IGNORECASE)
            if m:
                iso = m.group(1).strip()
                iso = re.sub(r"^(?:HIV-?1[-_.]?|HIV[-_.]?)", "", iso, flags=re.IGNORECASE)
                toks = [t for t in re.split(r"[_.\-#/]+", iso) if t]
                if len(toks) >= 2 and toks[0].upper() in ("TDR", "CN", "CY", "FRESH", "NV", "V704", "V703", "CAP", "RV217", "RV144", "P", "PAT", "ID", "SO"):
                    root = f"{toks[0]}_{toks[1]}".upper()
                elif len(toks) >= 1:
                    root = toks[0].upper()
                else:
                    root = iso[:8].upper()
                root = re.sub(r"(?:V\d+|W\d+|D\d+|M\d+)$", "", root)
                return f"{acc[:4]}:{root}"
            return f"{acc[:5]}:{cty}"

        # Enrich comb_df with flatfile annotations and check known CRFs revealed by flatfiles
        crf_notes_col, patients_col, clones_col, comments_col, authors_col, ff_iso_col, ff_title_col = [], [], [], [], [], [], []
        for idx_r, r in comb_df.iterrows():
            acc = str(r["accession"])
            ff_r = ff_map.get(acc, {})
            note_crf = s8._norm_str(ff_r.get("official_crf_note", ""))
            ff_t = s8._norm_str(ff_r.get("title", ""))
            crf_notes_col.append(note_crf)
            patients_col.append(s8._norm_str(ff_r.get("patient", "")))
            clones_col.append(s8._norm_str(ff_r.get("clone", "")))
            comments_col.append(s8._norm_str(ff_r.get("comment", "")))
            authors_col.append(s8._norm_str(ff_r.get("authors", "")))
            ff_iso_col.append(s8._norm_str(ff_r.get("isolate", "")))
            ff_title_col.append(ff_t)

            if not r["matched_known_crf_window"] and r["egg_category"] == "URF":
                t_crf = re.search(r"\b(CRF\s*\d{2,3}[_\-A-Za-z0-9]*)", ff_t + " " + note_crf, re.I)
                if t_crf:
                    comb_df.at[idx_r, "matched_known_crf_window"] = f"Flatfile:{t_crf.group(1).replace(' ', '')}"

        comb_df["official_crf_note"] = crf_notes_col
        comb_df["patient"] = patients_col
        comb_df["clone"] = clones_col
        comb_df["comment"] = comments_col
        comb_df["authors"] = authors_col
        comb_df["ff_isolate"] = ff_iso_col
        comb_df["ff_title"] = ff_title_col
        comb_df["strict_pat"] = [strict_patient_root(r) for _, r in comb_df.iterrows()]
        comb_df.to_csv(RESULTS_DIR / f"scale_partial_{region_name}_all_classifications.csv", index=False)

        print("\n" + "-" * 88)
        print(f"UNIFIED GENBANK + LANL PARTIAL-ONLY {region_name.upper()} RESULTS (N={len(comb_df):,} validated partial isolates)")
        print("-" * 88)
        cat_vc = comb_df["egg_category"].value_counts()
        for k_c, v_c in cat_vc.items():
            print(f"    * {k_c:<32s}: {v_c:6,d}")
        reg_hits = comb_df[comb_df["matched_novel15_lineage"] != ""]
        print(f"[A] Partial-{region_name.upper()} Carriers of Candidate Regional Mosaic (RU_G0B_gagMA): {len(reg_hits)}")
        for lin, grp in reg_hits.groupby("matched_novel15_lineage"):
            print(f"    * {lin}: +{len(grp)} partial isolates | Countries: {Counter(grp['country']).most_common(4)}")

        unid = comb_df[
            (comb_df["egg_category"] == "URF")
            & (comb_df["n_segments"] > 1)
            & (comb_df["matched_known_crf_window"] == "")
            & (comb_df["matched_novel15_lineage"] == "")
        ].copy()

        novel_clusters = []
        iqtree_tasks = []
        pending_clusters = []
        for old_f in TREES_DIR.glob(f"partial_{region_name}_*"):
            old_f.unlink()

        for w_form, grp in unid.groupby("window_formula"):
            if len(grp) < 3:
                continue
            dedup = grp.drop_duplicates("strict_pat").copy()
            if len(dedup) < 3:
                continue

            kept_rows = []
            kept_arrs = []
            kept_tokens = []
            for _, r in dedup.iterrows():
                acc = str(r["accession"])
                arr_cand = seq_arrays.get(acc)
                if arr_cand is None:
                    continue
                d_row = r.to_dict()
                if d_row.get("ff_title"):
                    d_row["title"] = d_row["ff_title"]
                tok_c = s8._extract_dedup_tokens(d_row)
                is_dup = False
                for ka, tok_r, r_prev in zip(kept_arrs, kept_tokens, kept_rows):
                    same_study = bool(tok_c["study_sig"] and tok_r["study_sig"] and tok_c["study_sig"] == tok_r["study_sig"])
                    if tok_c["patient_id"] and tok_r["patient_id"] and tok_c["patient_id"] == tok_r["patient_id"] and (same_study or d_row.get("country") == r_prev.get("country")):
                        is_dup = True
                        break
                    if tok_c["clone_id"] and tok_r["clone_id"] and tok_c["clone_id"] == tok_r["clone_id"] and (same_study or (tok_c["isolate_norm"] and tok_c["isolate_norm"] == tok_r["isolate_norm"])):
                        is_dup = True
                        break
                    if tok_c["isolate_norm"] and tok_r["isolate_norm"] and tok_c["isolate_norm"] == tok_r["isolate_norm"] and same_study:
                        is_dup = True
                        break
                    if tok_c["pair_key"] and tok_r["pair_key"] and tok_c["pair_key"] == tok_r["pair_key"]:
                        is_dup = True
                        break
                    both = (arr_cand >= 0) & (ka >= 0)
                    if np.sum(both) >= 800:
                        div = float(np.sum(arr_cand[both] != ka[both])) / int(np.sum(both))
                        if div < 0.008:
                            is_dup = True
                            break
                if not is_dup:
                    kept_rows.append(d_row)
                    kept_arrs.append(arr_cand)
                    kept_tokens.append(tok_c)

            if len(kept_rows) < 3:
                continue

            # Homologous breakpoint clustering (+-250 nt)
            bp_clusters = []
            for r, arr in zip(kept_rows, kept_arrs):
                segs = s10.parse_barcode_map(r["window_barcode"])
                bps = [bp for bp, _, _ in s10.get_internal_breakpoints(segs)]
                if not bps:
                    continue
                placed = False
                for cl in bp_clusters:
                    ref_bps = cl[0][1]
                    if len(bps) == len(ref_bps) and all(abs(b1 - b2) <= 250 for b1, b2 in zip(bps, ref_bps)):
                        cl.append((r, bps, arr, segs))
                        placed = True
                        break
                if not placed:
                    bp_clusters.append([(r, bps, arr, segs)])

            if not bp_clusters:
                continue
            best_cl = max(bp_clusters, key=len)
            if len(best_cl) < 3:
                continue

            h_rows = [x[0] for x in best_cl]
            h_arrs = [x[2] for x in best_cl]
            h_segs0 = best_cl[0][3]
            bp_medians = [np.median([x[1][j] for x in best_cl]) for j in range(len(best_cl[0][1]))]

            # Study independence check via flatfile author signatures + accession prefixes
            study_keys = set()
            for r in h_rows:
                auth_sig = str(r.get("authors") or "")[:40].strip().lower()
                acc_pre = str(r["accession"])[:4]
                study_keys.add(auth_sig if auth_sig else acc_pre)
            n_acc_prefixes = len({str(r["accession"])[:4] for r in h_rows})
            if len(study_keys) < 2 and n_acc_prefixes < 2:
                continue

            # Mean pairwise locus coding p-distance across homologous cluster
            pair_divs = []
            for i in range(len(h_arrs)):
                for j in range(i + 1, len(h_arrs)):
                    both = (h_arrs[i] >= 0) & (h_arrs[j] >= 0)
                    if np.sum(both) >= 800:
                        pair_divs.append(float(np.mean(h_arrs[i][both] != h_arrs[j][both])))
            mean_div_pct = round(100.0 * float(np.mean(pair_divs)), 2) if pair_divs else 999.0

            cl_idx = len(pending_clusters) + 1
            safe_form = re.sub(r"[^A-Za-z0-9]+", "_", w_form).strip("_")
            n_segs_queued = 0
            for seg_i, (h0, h1, tr) in enumerate(h_segs0, 1):
                if (h1 - h0 + 1) < 250:
                    continue
                c0_win = max(0, int(np.searchsorted(hxb2_coords_win, h0)))
                c1_win = min(expected_len, int(np.searchsorted(hxb2_coords_win, h1, side="right")))
                c0_ref = max(0, int(np.searchsorted(ref_hxb2_coords, h0)))
                c1_ref = min(len(ref_hxb2_coords), int(np.searchsorted(ref_hxb2_coords, h1, side="right")))
                if (c1_win - c0_win) < 250 or (c1_ref - c0_ref) != (c1_win - c0_win):
                    continue

                fasta_path = TREES_DIR / f"partial_{region_name}_c{cl_idx:02d}_{safe_form}_seg{seg_i}_{tr}_{h0}_{h1}.fasta"
                c_trunks = {seg_t for _, _, seg_t in h_segs0}
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
                        s_slice = s09.arr_to_fasta_slice(ref_arr, c0_ref, c1_ref)
                        if (len(s_slice) - s_slice.count("-")) >= 150:
                            f.write(f">{ref_name}\n{s_slice}\n")
                    if use_ae_a1:
                        for ref_name, ref_tr, ref_arr in a1_ae_refs:
                            s_slice = s09.arr_to_fasta_slice(ref_arr, c0_ref, c1_ref)
                            if (len(s_slice) - s_slice.count("-")) >= 150:
                                f.write(f">{ref_name}\n{s_slice}\n")
                    elif use_ag_a1:
                        for ref_name, ref_tr, ref_arr in a1_ag_refs:
                            s_slice = s09.arr_to_fasta_slice(ref_arr, c0_ref, c1_ref)
                            if (len(s_slice) - s_slice.count("-")) >= 150:
                                f.write(f">{ref_name}\n{s_slice}\n")
                    if h1 > 8200 and h0 < 8900:
                        for ref_name, ref_tr, ref_arr in crf_g0_tail:
                            s_slice = s09.arr_to_fasta_slice(ref_arr, c0_ref, c1_ref)
                            if (len(s_slice) - s_slice.count("-")) >= 150:
                                f.write(f">{ref_name}\n{s_slice}\n")
                    for r, arr in zip(h_rows, h_arrs):
                        s_slice = s09.arr_to_fasta_slice(arr, c0_win, c1_win)
                        if (len(s_slice) - s_slice.count("-")) >= 150:
                            qn = f"Query_{r['accession']}"
                            q_names.append(qn)
                            f.write(f">{qn}\n{s_slice}\n")
                if len(q_names) >= 3:
                    iqtree_tasks.append((cl_idx, w_form, seg_i, tr, h0, h1, fasta_path, q_names))
                    n_segs_queued += 1

            pending_clusters.append({
                "cl_idx": cl_idx,
                "w_form": w_form,
                "grp": grp,
                "h_rows": h_rows,
                "n_studies": max(len(study_keys), n_acc_prefixes),
                "bp_medians": bp_medians,
                "mean_div_pct": mean_div_pct,
                "n_segs_queued": n_segs_queued,
            })

        # Run IQ-TREE 2 GTR+G on all partial candidate segments >= 250 nt
        partial_tree_results = []
        if iqtree_tasks:
            n_workers = min(8, cpu_count())
            print(f"[+] Running IQ-TREE 2 (GTR+G, 1000 UFBoot) on {len(iqtree_tasks)} partial {region_name} segments across {len(pending_clusters)} candidate clusters...")
            with Pool(processes=n_workers) as pool:
                for res in pool.imap_unordered(s09.run_iqtree_task, iqtree_tasks):
                    partial_tree_results.append(res)

        pt_df = pd.DataFrame(partial_tree_results)
        if not pt_df.empty:
            pt_df = pt_df.sort_values(["lineage_rank", "segment_idx"]).reset_index(drop=True)
            pt_df.to_csv(RESULTS_DIR / f"scale_partial_{region_name}_phylogenetic_validation.csv", index=False)

        for pc in pending_clusters:
            cl_idx = pc["cl_idx"]
            sub_t = pt_df[pt_df["lineage_rank"] == cl_idx] if not pt_df.empty else pd.DataFrame()
            n_seg = len(sub_t)
            n_mono = int(sub_t["is_strictly_monophyletic"].sum()) if n_seg else 0
            min_ufb = float(sub_t["lineage_cluster_ufboot_pct"].min()) if n_seg and sub_t["lineage_cluster_ufboot_pct"].notna().all() else np.nan
            max_pat = float(sub_t["mean_intracluster_patristic_dist"].max()) if n_seg and sub_t["mean_intracluster_patristic_dist"].notna().all() else np.nan
            mean_div_pct = pc["mean_div_pct"]
            strict_mono = bool(
                n_seg >= 2
                and n_mono == n_seg
                and pd.notna(min_ufb) and min_ufb >= 80.0
                and pd.notna(mean_div_pct) and mean_div_pct <= 6.0
                and pd.notna(max_pat) and max_pat <= 0.10
            )
            h_df = pd.DataFrame(pc["h_rows"])
            ctys = ", ".join(f"{k}({v})" for k, v in Counter(h_df["country"]).most_common(5))
            srcs = ", ".join(f"{k}({v})" for k, v in Counter(h_df["source"]).most_common(2))
            accs_str = ", ".join(h_df["accession"].astype(str))
            bp_str = ", ".join(str(int(round(m))) for m in pc["bp_medians"])
            novel_clusters.append({
                "region": region_name,
                "window_formula": pc["w_form"],
                "phylogenetic_classification": "Candidate_Partial_Mosaic" if strict_mono else "Convergent_Hotspot_Polyphyletic",
                "all_segments_strictly_monophyletic": strict_mono,
                "n_raw_partial_isolates": len(pc["grp"]),
                "n_independent_patients": len(h_df),
                "n_distinct_studies": pc["n_studies"],
                "mean_pairwise_div_pct": mean_div_pct,
                "n_segments_tested": n_seg,
                "n_segments_monophyletic": n_mono,
                "min_lineage_ufboot_pct": round(min_ufb, 1) if pd.notna(min_ufb) else np.nan,
                "max_segment_patristic_dist": round(max_pat, 5) if pd.notna(max_pat) else np.nan,
                "median_breakpoints_hxb2": bp_str,
                "countries": ctys,
                "sources": srcs,
                "member_accessions": accs_str,
                "representative_barcode": h_df.iloc[0]["window_barcode"],
            })

        clust_cols = [
            "region", "window_formula", "phylogenetic_classification", "all_segments_strictly_monophyletic",
            "n_raw_partial_isolates", "n_independent_patients", "n_distinct_studies", "mean_pairwise_div_pct",
            "n_segments_tested", "n_segments_monophyletic", "min_lineage_ufboot_pct", "max_segment_patristic_dist",
            "median_breakpoints_hxb2", "countries", "sources", "member_accessions", "representative_barcode",
        ]
        clust_df = pd.DataFrame(novel_clusters, columns=clust_cols)
        if len(clust_df) > 0:
            clust_df = clust_df.sort_values(
                ["all_segments_strictly_monophyletic", "n_independent_patients", "mean_pairwise_div_pct"],
                ascending=[False, False, True],
            )
        clust_df.to_csv(RESULTS_DIR / f"scale_partial_{region_name}_novel_circulating_lineages.csv", index=False)

        n_pass = int(clust_df["all_segments_strictly_monophyletic"].sum()) if len(clust_df) else 0
        print(f"\n[B] Partial-{region_name.upper()} Homologous Multi-Study Clusters Tested: {len(clust_df)} | Strict Monophyletic (Candidate_Partial_Mosaic): {n_pass}")
        for _, r in clust_df.iterrows():
            print(f"    * {r['window_formula']:40s} | {r['phylogenetic_classification']:31s} | N={r['n_independent_patients']:2d} | div={r['mean_pairwise_div_pct']:.2f}% | mono={r['n_segments_monophyletic']}/{r['n_segments_tested']} | minUFB={r['min_lineage_ufboot_pct']} | maxPat={r['max_segment_patristic_dist']} | {r['countries']}")

    print(f"\n[+] Step 12 Partial Scale Screen complete in {time.time() - t0_all:.1f}s")


if __name__ == "__main__":
    main()
