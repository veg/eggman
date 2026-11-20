#!/usr/bin/env python3
"""
08b_fetch_genbank_cluster_metadata.py
=====================================
Fetches authentic NCBI GenBank flatfiles (`rettype=gb&retmode=text`) via E-utilities
for:
  1. All multi-sequence URF / candidate novel CRF clusters in GenBank 31K & Panel 2
  2. All multi-sequence partial pol / env candidate clusters & Panel 3
  3. All LANL 2023 reference genomes (for exact NCBI first_submitted_year)
  4. Any extra accessions passed via CLI `--accs` or helper `ensure_flatfiles_cached(accs)`

Parses full flatfile metadata into `data/genbank_cluster_flatfiles.tsv`:
  accession, locus_date, first_submitted_date, first_submitted_year,
  collection_date, country, host, isolate, patient, clone, title,
  authors, comment, note, official_crf_note, source_qualifiers_json
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import time
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional

import pandas as pd

PROJECT_ROOT = Path("/Users/sergei/Projects/HIV-recombination")
DATA_DIR = PROJECT_ROOT / "data"
RESULTS_DIR = PROJECT_ROOT / "results"
OUT_TSV = DATA_DIR / "genbank_cluster_flatfiles.tsv"

EFETCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
MONTH_MAP = {
    "JAN": 1, "FEB": 2, "MAR": 3, "APR": 4, "MAY": 5, "JUN": 6,
    "JUL": 7, "AUG": 8, "SEP": 9, "OCT": 10, "NOV": 11, "DEC": 12,
}

TSV_COLUMNS = [
    "accession",
    "locus_date",
    "first_submitted_date",
    "first_submitted_year",
    "collection_date",
    "collection_year",
    "country",
    "host",
    "isolate",
    "patient",
    "clone",
    "title",
    "authors",
    "comment",
    "note",
    "official_crf_note",
]


def parse_gb_date_to_decimal(d_str: str) -> Optional[float]:
    s = (d_str or "").strip()
    if not s:
        return None
    m1 = re.match(r"^(\d{1,2})-([A-Za-z]{3})-(\d{4})$", s)
    if m1:
        day = int(m1.group(1))
        mon = MONTH_MAP.get(m1.group(2).upper())
        yr = int(m1.group(3))
        if mon:
            return round(yr + ((mon - 1) * 30.4375 + (day - 0.5)) / 365.25, 4)
    m2 = re.match(r"^([A-Za-z]{3})-(\d{4})$", s)
    if m2:
        mon = MONTH_MAP.get(m2.group(1).upper())
        yr = int(m2.group(2))
        if mon:
            return round(yr + (mon - 0.5) / 12.0, 4)
    m3 = re.match(r"^(\d{4})-(\d{2})-(\d{2})$", s)
    if m3:
        yr, mon, day = int(m3.group(1)), int(m3.group(2)), int(m3.group(3))
        if 1 <= mon <= 12 and 1 <= day <= 31:
            return round(yr + ((mon - 1) * 30.4375 + (day - 0.5)) / 365.25, 4)
    m4 = re.match(r"^(\d{4})-(\d{2})$", s)
    if m4:
        yr, mon = int(m4.group(1)), int(m4.group(2))
        if 1 <= mon <= 12:
            return round(yr + (mon - 0.5) / 12.0, 4)
    m5 = re.match(r"^(\d{4})$", s)
    if m5:
        return float(int(m5.group(1))) + 0.5
    return None


def extract_official_crf_mentions(text: str) -> str:
    if not text:
        return ""
    found = []
    for m in re.findall(r"\bCRF[_\s\-]*0*(\d{1,3})[_\-]?([A-Za-z0-9]+)?\b", text, re.IGNORECASE):
        num = int(m[0])
        suf = (m[1] or "").strip()
        # Ignore vague words right after CRF numbers
        if suf.lower() in {"in", "of", "and", "with", "from", "to", "for", "is", "are", "was", "or", "genome", "strains", "strain", "sequence", "sequences", "like", "related"}:
            suf = ""
        tok = f"CRF{num:02d}_{suf}" if suf else f"CRF{num:02d}"
        if tok not in found:
            found.append(tok)
    return ";".join(found)


def parse_single_gb_record(raw_rec: str) -> Optional[Dict[str, Any]]:
    lines = raw_rec.splitlines()
    if not lines:
        return None

    acc = ""
    locus_date = ""
    def_lines: List[str] = []
    comment_lines: List[str] = []
    ref_blocks: List[Dict[str, str]] = []
    source_quals: Dict[str, List[str]] = {}

    section = ""
    cur_ref: Optional[Dict[str, str]] = None
    cur_ref_field = ""
    in_source_feat = False
    cur_qual_key = ""
    cur_qual_val: List[str] = []

    def _flush_qual():
        nonlocal cur_qual_key, cur_qual_val
        if cur_qual_key:
            val = " ".join(cur_qual_val).strip()
            if val.startswith('"') and val.endswith('"'):
                val = val[1:-1]
            source_quals.setdefault(cur_qual_key, []).append(val)
        cur_qual_key = ""
        cur_qual_val = []

    for line in lines:
        if not line.strip():
            continue
        top_key = line[:12].strip() if not line.startswith(" ") else ""

        if top_key == "LOCUS":
            section = "LOCUS"
            parts = line.split()
            if len(parts) >= 2:
                acc = parts[1].strip().split(".")[0]
            m_ld = re.search(r"(\d{2}-[A-Za-z]{3}-\d{4})\s*$", line)
            if m_ld:
                locus_date = m_ld.group(1).upper()
            continue

        if top_key == "VERSION":
            section = "VERSION"
            parts = line.split()
            if len(parts) >= 2:
                acc = parts[1].strip().split(".")[0]
            continue

        if top_key == "ACCESSION":
            if not acc:
                parts = line.split()
                if len(parts) >= 2:
                    acc = parts[1].strip().split(".")[0]
            section = "ACCESSION"
            continue

        if top_key == "DEFINITION":
            section = "DEFINITION"
            def_lines.append(line[12:].strip())
            continue

        if top_key == "COMMENT":
            section = "COMMENT"
            comment_lines.append(line[12:].strip())
            continue

        if top_key == "REFERENCE":
            section = "REFERENCE"
            cur_ref = {"authors": "", "title": "", "journal": ""}
            ref_blocks.append(cur_ref)
            cur_ref_field = ""
            continue

        if top_key == "FEATURES":
            section = "FEATURES"
            in_source_feat = False
            continue

        if top_key in ("ORIGIN", "CONTIG", "BASE", "//"):
            _flush_qual()
            in_source_feat = False
            section = top_key
            continue

        if top_key:
            _flush_qual()
            in_source_feat = False
            section = top_key
            continue

        # Continuation lines
        if section == "DEFINITION":
            def_lines.append(line[12:].strip())
        elif section == "COMMENT":
            comment_lines.append(line[12:].strip())
        elif section == "REFERENCE" and cur_ref is not None:
            sub_k = line[:12].strip()
            if sub_k in ("AUTHORS", "CONSRTM"):
                cur_ref_field = "authors"
                cur_ref["authors"] += (" " if cur_ref["authors"] else "") + line[12:].strip()
            elif sub_k == "TITLE":
                cur_ref_field = "title"
                cur_ref["title"] += (" " if cur_ref["title"] else "") + line[12:].strip()
            elif sub_k == "JOURNAL":
                cur_ref_field = "journal"
                cur_ref["journal"] += (" " if cur_ref["journal"] else "") + line[12:].strip()
            elif sub_k in ("PUBMED", "REMARK"):
                cur_ref_field = ""
            elif cur_ref_field:
                cur_ref[cur_ref_field] += " " + line[12:].strip()
        elif section == "FEATURES":
            feat_hdr = line[5:21].strip() if len(line) > 5 else ""
            if feat_hdr:
                _flush_qual()
                in_source_feat = (feat_hdr == "source")
                continue
            if in_source_feat:
                sline = line[21:] if len(line) > 21 else line.strip()
                if sline.startswith("/"):
                    _flush_qual()
                    if "=" in sline:
                        qk, qv = sline[1:].split("=", 1)
                        cur_qual_key = qk.strip()
                        cur_qual_val = [qv.strip()]
                    else:
                        cur_qual_key = sline[1:].strip()
                        cur_qual_val = ["True"]
                elif cur_qual_key:
                    cur_qual_val.append(sline.strip())

    _flush_qual()
    if not acc:
        return None

    # Find earliest Submitted journal date, or earliest publication/locus date
    sub_dates = []
    authors_all = []
    titles_all = []
    for rb in ref_blocks:
        if rb.get("authors"):
            authors_all.append(rb["authors"])
        if rb.get("title") and rb["title"] != "Direct Submission":
            titles_all.append(rb["title"])
        j = rb.get("journal", "")
        m_sub = re.search(r"Submitted\s*\((\d{1,2}-[A-Za-z]{3}-\d{4})\)", j)
        if m_sub:
            ds = m_sub.group(1).upper()
            dec = parse_gb_date_to_decimal(ds)
            if dec is not None:
                sub_dates.append((dec, ds))

    if sub_dates:
        sub_dates.sort()
        first_sub_dec, first_sub_str = sub_dates[0]
    else:
        first_sub_str = locus_date
        first_sub_dec = parse_gb_date_to_decimal(locus_date)

    col_dates = source_quals.get("collection_date", [])
    col_date_str = col_dates[0] if col_dates else ""
    col_dec = parse_gb_date_to_decimal(col_date_str)
    country = (source_quals.get("geo_loc_name") or source_quals.get("country") or [""])[0]
    host = (source_quals.get("host") or [""])[0]
    isolate = (source_quals.get("isolate") or source_quals.get("strain") or [""])[0]
    patient = (source_quals.get("patient") or source_quals.get("isolation_source") or [""])[0]
    clone = (source_quals.get("clone") or [""])[0]
    note_str = " | ".join(source_quals.get("note", []))
    def_str = " ".join(def_lines).strip()
    comment_str = " ".join(comment_lines).strip()
    title_joined = " | ".join(titles_all) if titles_all else def_str
    authors_str = " | ".join(authors_all[:2])

    crf_text = " ".join([def_str, title_joined, comment_str, note_str])
    crf_note = extract_official_crf_mentions(crf_text)

    # Also infer patient token from note/comment if not explicit /patient qualifier
    if not source_quals.get("patient"):
        m_pat = re.search(r"\bpatient\s*[:#]?\s*([A-Za-z0-9_\-]+)", note_str + " " + comment_str + " " + title_joined, re.IGNORECASE)
        if m_pat:
            patient = f"patient:{m_pat.group(1)}"

    return {
        "accession": acc,
        "locus_date": locus_date,
        "first_submitted_date": first_sub_str,
        "first_submitted_year": round(first_sub_dec, 4) if first_sub_dec is not None else "",
        "collection_date": col_date_str,
        "collection_year": round(col_dec, 4) if col_dec is not None else "",
        "country": country,
        "host": host,
        "isolate": isolate,
        "patient": patient,
        "clone": clone,
        "title": title_joined,
        "authors": authors_str,
        "comment": comment_str,
        "note": note_str,
        "official_crf_note": crf_note,
    }


def fetch_gb_batch(accs: List[str], retries: int = 4) -> List[Dict[str, Any]]:
    ids = ",".join(accs)
    params = urllib.parse.urlencode({"db": "nuccore", "id": ids, "rettype": "gb", "retmode": "text"})
    url = f"{EFETCH_URL}?{params}"
    last_err = None
    txt = ""
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "HIV-recomb-audit/1.0"})
            with urllib.request.urlopen(req, timeout=90) as resp:
                txt = resp.read().decode("utf-8", errors="replace")
            break
        except Exception as e:
            last_err = e
            time.sleep(1.5 * (2 ** attempt))
    if not txt:
        raise RuntimeError(f"NCBI efetch failed for {len(accs)} accessions: {last_err}")

    rows = []
    for raw_rec in txt.split("\n//\n"):
        if not raw_rec.strip():
            continue
        parsed = parse_single_gb_record(raw_rec)
        if parsed:
            rows.append(parsed)
    return rows


def ensure_flatfiles_cached(accessions: Iterable[str], out_tsv: Path = OUT_TSV) -> pd.DataFrame:
    acc_set = sorted({str(a).strip().split(".")[0] for a in accessions if str(a).strip() and str(a).strip() != "nan"})
    if out_tsv.exists() and out_tsv.stat().st_size > 100:
        existing_df = pd.read_csv(out_tsv, sep="\t", dtype=str).fillna("")
    else:
        existing_df = pd.DataFrame(columns=TSV_COLUMNS)

    done = set(existing_df["accession"].astype(str))
    missing = [a for a in acc_set if a not in done]
    print(f"[*] GenBank flatfile cache: {len(done)} already cached, {len(missing)} missing to fetch.", flush=True)

    new_rows = []
    batch_size = 150
    for i in range(0, len(missing), batch_size):
        chunk = missing[i : i + batch_size]
        rows = fetch_gb_batch(chunk)
        new_rows.extend(rows)
        print(f"    Fetched batch {i // batch_size + 1}/{(len(missing) + batch_size - 1) // batch_size}: +{len(rows)} records", flush=True)
        if i + batch_size < len(missing):
            time.sleep(0.45)

    if new_rows:
        comb = pd.concat([existing_df, pd.DataFrame(new_rows)], ignore_index=True)
    else:
        comb = existing_df
    comb = comb.drop_duplicates(subset=["accession"], keep="last").sort_values("accession").reset_index(drop=True)
    comb.to_csv(out_tsv, sep="\t", index=False)
    print(f"[+] Saved {len(comb)} GenBank flatfile records to {out_tsv}", flush=True)
    return comb


def collect_target_accessions() -> List[str]:
    accs = set()

    cls_p = RESULTS_DIR / "genbank_31k_rhizaeon_classifications.csv"
    if cls_p.exists():
        gb = pd.read_csv(cls_p)
        u = gb[gb["recomb_category"] == "Unidentified_Novel_Recombinant"]
        vc = u["rhizaeon_formula"].value_counts()
        multi_forms = set(vc[vc >= 2].index)
        accs.update(u[u["rhizaeon_formula"].isin(multi_forms)]["accession"].astype(str))

    nov_p = RESULTS_DIR / "genbank_novel_circulating_recombinants.csv"
    if nov_p.exists():
        nov = pd.read_csv(nov_p)
        for s in nov["member_accessions"].dropna():
            for a in str(s).split(","):
                if a.strip():
                    accs.add(a.strip())

    for reg in ("pol", "env"):
        cp = RESULTS_DIR / f"scale_partial_{reg}_novel_circulating_lineages.csv"
        ap = RESULTS_DIR / f"scale_partial_{reg}_all_classifications.csv"
        if cp.exists() and ap.exists():
            cand = pd.read_csv(cp)
            c_forms = set(cand["window_formula"].astype(str))
            allc = pd.read_csv(ap)
            sub = allc[
                allc["window_formula"].isin(c_forms)
                & allc["matched_known_crf_window"].isna()
                & allc["matched_novel15_lineage"].isna()
            ]
            accs.update(sub["accession"].astype(str))
            for s in cand["example_accessions"].dropna():
                for a in str(s).split(","):
                    if a.strip():
                        accs.add(a.strip())

    bench_p = DATA_DIR / "benchmark_rhiztype_vs_comet.fasta"
    if bench_p.exists():
        for line in open(bench_p):
            if line.startswith(">P2|") or line.startswith(">P3|"):
                accs.add(line[1:].strip().split("|")[-1].split(".")[0])

    ref_meta_p = DATA_DIR / "lanl_ref_2023_metadata.csv"
    if ref_meta_p.exists():
        rm = pd.read_csv(ref_meta_p)
        for a in rm["accession"].dropna():
            a_s = str(a).strip().split(".")[0]
            if a_s and re.match(r"^[A-Z0-9_]+$", a_s):
                accs.add(a_s)

    return sorted(accs)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--accs", nargs="*", default=None, help="Optional explicit accessions to fetch")
    args = parser.parse_args()

    targets = collect_target_accessions()
    if args.accs:
        targets = sorted(set(targets) | set(args.accs))
    print(f"[*] Target accessions collected for GenBank flatfile fetch: {len(targets)}")
    df = ensure_flatfiles_cached(targets)
    n_crf = int((df["official_crf_note"].astype(str) != "").sum())
    n_pat = int((df["patient"].astype(str) != "").sum())
    n_clone = int((df["clone"].astype(str) != "").sum())
    print(f"[+] Flatfile summary: {len(df)} records | official_crf_note={n_crf} | patient_info={n_pat} | clone_info={n_clone}")


if __name__ == "__main__":
    main()
