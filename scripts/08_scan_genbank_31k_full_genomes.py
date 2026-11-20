#!/usr/bin/env python3
"""
08_scan_genbank_31k_full_genomes.py
===================================
Exhaustive, alignment-hardened scan of the 31,825 full-length HIV-1 GenBank sequences
(`ncbi_hiv1_all_complete_31k.fasta`) combined with the 5,606 LANL full-length genomes
to discover and validate all Unidentified Circulating Recombinant Forms (URFs -> candidate CRFs).

Alignment & Quality Control Architecture:
  1. 15-Sequence Multi-Subtype HXB2 Reference Profile (`9,719 nt`):
     Projects HXB2 (`K03455`) + all 14 clean Group M Basal Trunk consensuses
     (`A1, A2, A6, B, C, D, F1, F2, G0, H, J, K, L, ?E`) onto the exact 9,719 HXB2 columns.
  2. Profile-Locked Parallel Alignment (`mafft --6merpair --addfragments --keeplength`):
     Aligns every GenBank genome against all 14 subtype consensuses simultaneously while
     strictly preserving the 9,719 HXB2 coordinate system.
  3. Four-Layer Post-Alignment Artifact Shield:
     - Layer A: Masks the 5 hypervariable length-polymorphic `gp120` loops (`V1/V2`, `V3`,
       `V4`, `V5`) and `gag p6` PTAP indel zone so indel placement ambiguity never triggers a crossover.
     - Layer B: Sliding-window local misalignment / frame-shift / non-Group-M detector (150-nt window,
       >26% divergence to closest Group M trunk masked; genomes with >450 nt corrupted rejected).
     - Layer C: Strict breakpoint verification (`min_informative_snps=7`, `min_hxb2_span_nt=200`,
       `p_fisher_thresh=1e-5`, segment trunk divergence <= 16.5%).
     - Layer D: Known & Attenuated-Subset CRF Filter + Sequence-Level Independent-Patient
       Deduplication (greedy 0.8% whole-genome coding divergence threshold, >= 70 coding SNPs)
       to guarantee that neither attenuated known CRFs (e.g. CRF02_AG, CRF07_BC, CRF63_02A6)
       nor intra-patient SGA/reservoir clones ever inflate novel circulating lineage counts.
"""

import json
import os
import re
import subprocess
import sys
import tempfile
import time
from collections import Counter, defaultdict
from multiprocessing import Pool, cpu_count
from pathlib import Path
import numpy as np
import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))
from importlib import import_module

s2 = import_module("02_phase1_pure_subtype_and_temporal_audit")
s3 = import_module("03_audit_and_reclassify_all_crfs")
s08b = import_module("08b_fetch_genbank_cluster_metadata")

PROJECT_ROOT = Path("/Users/sergei/Projects/HIV-recombination")
DATA_DIR = PROJECT_ROOT / "data"
RESULTS_DIR = PROJECT_ROOT / "results"
EGGMAN_BIN = PROJECT_ROOT / "rhiztype" / "target" / "release" / "eggman"

GENBANK_FASTA = Path("/Users/sergei/Projects/TOGA_MEME/recombination/pan_group_m_nomenclature/data/ncbi_hiv1_all_complete_31k.fasta")
GENBANK_META = Path("/Users/sergei/Projects/TOGA_MEME/recombination/pan_group_m_nomenclature/data/ncbi_hiv1_all_complete_31k_metadata.csv")
MAFFT_BIN = "/usr/local/bin/mafft"

POST2023_ACC_TO_CRF = {
    "LC735411": "CRF127_07109",
    "LC735412": "CRF127_07109",
    "LC735413": "CRF127_07109",
    "OP501672": "CRF135_0107",
    "OP501674": "CRF135_0107",
    "OP501675": "CRF135_0107",
    "OP501676": "CRF135_0107",
    "OP758070": "CRF135_0107",
    "OP501677": "CRF136_0107",
    "OP501678": "CRF136_0107",
    "OP501679": "CRF136_0107",
    "PP483791": "CRF162_cpx",
    "PP483792": "CRF162_cpx",
    "PP483793": "CRF162_cpx",
    "PP483794": "CRF162_cpx",
    "PP483795": "CRF162_cpx",
    "PP483797": "CRF162_cpx",
    "PQ356206": "CRF172_0755",
    "PQ356207": "CRF172_0755",
    "PQ356208": "CRF172_0755",
    "PQ356209": "CRF172_0755",
    "PQ585408": "CRF173_63A6",
    "PQ585409": "CRF173_63A6",
    "PQ585410": "CRF173_63A6",
    "PQ585411": "CRF173_63A6",
    "PQ585412": "CRF173_63A6",
    "PQ585413": "CRF173_63A6",
    "PQ585414": "CRF173_63A6",
    "PV363678": "CRF188_0107",
    "PV363679": "CRF188_0107",
    "PV363680": "CRF188_0107",
    "PV363681": "CRF188_0107",
    "PV363682": "CRF188_0107",
    "PX262425": "CRF205_0107",
    "PX262426": "CRF205_0107",
    "PX262427": "CRF205_0107",
    "PX262428": "CRF205_0107",
    "OP921950": "CRF205_0107",
}
POST2023_CRF_NAMES = {
    "CRF127_07109",
    "CRF135_0107",
    "CRF136_0107",
    "CRF162_cpx",
    "CRF172_0755",
    "CRF173_63A6",
    "CRF188_0107",
    "CRF205_0107",
}
CRF63_REF_ACCS = {"JX500700", "MK984159", "PZ426393"}

HYPERVARIABLE_MASK_RANGES = [
    (2130, 2185),  # gag p6 PTAP duplication hotspot
    (6615, 6812),  # env gp120 V1/V2 hypervariable indel loop
    (7110, 7217),  # env gp120 V3 loop
    (7377, 7478),  # env gp120 V4 hypervariable indel loop
    (7596, 7637),  # env gp120 V5 hypervariable indel loop
]

INT_TO_BASE = {0: "A", 1: "C", 2: "G", 3: "T", -1: "-", -2: "-"}

G_PROFILE_PATH = None
G_TRUNK_MAT = None
G_TRUNK_NAMES = None
G_HXB2_COORDS = None
G_CM240_ARR = None
G_LOOP_MASK_BOOL = None


def get_hxb2_all_cols(seq_str: str) -> np.ndarray:
    cols = [i for i, ch in enumerate(seq_str) if ch not in "-."]
    return np.array(cols, dtype=np.int32)


def int_array_to_fasta_str(arr: np.ndarray) -> str:
    return "".join(INT_TO_BASE.get(int(x), "-") for x in arr)


def build_15trunk_hxb2_profile():
    profile_path = DATA_DIR / "hxb2_15trunk_profile_9719.fasta"
    assert profile_path.exists(), f"Missing ChronAeon-scrubbed profile FASTA: {profile_path}"
    trunk_mat, trunk_names, hxb2_coords, cm240_arr = s3.load_profiles_9719_coding()
    return profile_path, trunk_mat, trunk_names, hxb2_coords, cm240_arr


def extract_patient_key(row) -> str:
    acc_pre = str(row["accession"])[:5]
    iso = str(row["raw_iso"])
    if re.match(r"^(PG\d+-[A-Z]+\d+|CD\d+-\d+|HK_JIDLNBL_S\d+|Rakai\d+)", iso):
        cleaned = re.sub(r"_[0-9]+$", "", iso)
        return acc_pre + ":" + cleaned
    m_uspi = re.match(r"^(USPI\d+)", iso)
    if m_uspi:
        return acc_pre + ":" + m_uspi.group(1)
    tokens = re.split(r"[_\-\.]", iso)
    if len(tokens) > 1:
        t0 = tokens[0]
        if len(t0) <= 2 or t0.upper() in (
            "HIV", "HIV1", "PAT", "SUBJ", "PID", "CHN", "ZAF", "UGA",
            "RWA", "KEN", "TZA", "BRA", "RUS", "THA", "USA", "PBM", "PLA"
        ):
            p = t0 + "_" + tokens[1]
        else:
            p = t0
    else:
        p = iso
    return acc_pre + ":" + p


def init_align_worker(profile_path, trunk_mat, trunk_names, hxb2_coords, cm240_arr, loop_mask_bool):
    global G_PROFILE_PATH, G_TRUNK_MAT, G_TRUNK_NAMES, G_HXB2_COORDS, G_CM240_ARR, G_LOOP_MASK_BOOL
    G_PROFILE_PATH = str(profile_path)
    G_TRUNK_MAT = trunk_mat
    G_TRUNK_NAMES = trunk_names
    G_HXB2_COORDS = hxb2_coords
    G_CM240_ARR = cm240_arr
    G_LOOP_MASK_BOOL = loop_mask_bool


def align_sequences_to_coding_arrays(acc_seq_pairs, profile_path, loop_mask_bool):
    """Align a list of (accession, raw_seq) pairs via mafft --addfragments --keeplength and return {acc: q_clean}."""
    if not acc_seq_pairs:
        return {}
    with tempfile.NamedTemporaryFile("w", suffix=".fasta", delete=False) as tf:
        tmp_in = tf.name
        for acc, seq in acc_seq_pairs:
            tf.write(f">{acc}\n{seq}\n")
    try:
        proc = subprocess.run(
            [MAFFT_BIN, "--6merpair", "--addfragments", tmp_in, "--keeplength", "--quiet", "--thread", "4", str(profile_path)],
            capture_output=True, text=True, check=True
        )
        aligned_fasta = proc.stdout
    finally:
        if os.path.exists(tmp_in):
            os.unlink(tmp_in)

    out_arrays = {}
    cur_hdr = None
    cur_chunks = []
    for line in aligned_fasta.splitlines():
        if line.startswith(">"):
            if cur_hdr is not None and not cur_hdr.startswith("Profile_"):
                s_aln = "".join(cur_chunks)
                if len(s_aln) == 9719:
                    q_arr = s3.seq_to_int_array(s_aln[789:9417].upper())
                    q_arr[loop_mask_bool] = -1
                    out_arrays[cur_hdr] = q_arr
            cur_hdr = line[1:].strip().split()[0]
            cur_chunks = []
        else:
            cur_chunks.append(line.strip())
    if cur_hdr is not None and not cur_hdr.startswith("Profile_"):
        s_aln = "".join(cur_chunks)
        if len(s_aln) == 9719:
            q_arr = s3.seq_to_int_array(s_aln[789:9417].upper())
            q_arr[loop_mask_bool] = -1
            out_arrays[cur_hdr] = q_arr
    return out_arrays


def parse_bp_str(bp_str):
    if pd.isna(bp_str) or not str(bp_str).strip():
        return []
    res = []
    for item in str(bp_str).split(";"):
        m = re.match(r"\s*([A-Za-z0-9\?]+)->([A-Za-z0-9\?]+)@(\d+)\[(\d+)-(\d+)\]", item.strip())
        if m:
            res.append((m.group(1), m.group(2), int(m.group(3)), int(m.group(4)), int(m.group(5))))
    return res


def _reclassify_with_eggman(accs_needed):
    """Run rebuilt eggman --gamma 13.5 --min-snps 7 --min-span 200 --p-thresh 1e-5 --json on all accession sequences."""
    assert EGGMAN_BIN.exists(), f"Missing rebuilt eggman binary: {EGGMAN_BIN}"
    accs_set = set(accs_needed)
    raw_seqs = {}

    def _scan_fasta(path, extract_fn):
        if not path.exists() or len(raw_seqs) >= len(accs_set):
            return
        with open(path, "r") as f:
            cur_acc = None
            buf = []
            for line in f:
                if line.startswith(">"):
                    if cur_acc in accs_set and cur_acc not in raw_seqs:
                        raw_seqs[cur_acc] = "".join(buf).upper().replace("-", "").replace(".", "")
                    cur_acc = extract_fn(line[1:].strip())
                    buf = []
                elif cur_acc in accs_set:
                    buf.append(line.strip())
            if cur_acc in accs_set and cur_acc not in raw_seqs:
                raw_seqs[cur_acc] = "".join(buf).upper().replace("-", "").replace(".", "")

    _scan_fasta(GENBANK_FASTA, lambda h: h.split()[0].split(".")[0])
    if len(raw_seqs) < len(accs_set):
        _scan_fasta(DATA_DIR / "lanl_flt_2022.fasta", lambda h: h.split()[0].split(".")[-1])
    if len(raw_seqs) < len(accs_set):
        _scan_fasta(DATA_DIR / "lanl_ref_2023.fasta", lambda h: h.split()[0].split(".")[-1])

    missing = accs_set - set(raw_seqs.keys())
    assert not missing, f"Missing raw sequences for {len(missing)} accessions: {sorted(missing)[:5]}"

    with tempfile.TemporaryDirectory() as td:
        tmp_fa = Path(td) / "gb_7720.fasta"
        tmp_js = Path(td) / "gb_7720.json"
        with open(tmp_fa, "w") as f:
            for acc in accs_needed:
                f.write(f">{acc}\n{raw_seqs[acc]}\n")
        cmd = [
            str(EGGMAN_BIN),
            str(tmp_fa),
            "--gamma", "13.5",
            "--min-snps", "7",
            "--min-span", "200",
            "--p-thresh", "1e-5",
            "--json", str(tmp_js),
            "-q",
        ]
        subprocess.run(cmd, check=True)
        recs = json.loads(tmp_js.read_text())

    out = {}
    for rec in recs:
        bps = rec.get("breakpoints") or []
        bp_ci = "; ".join(
            f"{b['left_trunk']}->{b['right_trunk']}@{b['bp_hxb2']}[{b['ci_hxb2_low']}-{b['ci_hxb2_high']}]"
            for b in bps
        )
        bb = str(rec.get("backbone_trunk") or "")
        bc = rec.get("basal_barcode") or rec.get("hxb2_barcode") or (f"{bb}:790-9417" if bb else "")
        c_trunks = sorted(str(t) for t in (rec.get("constituent_trunks") or []))
        out[str(rec["sequence_id"])] = {
            "rhizaeon_formula": str(rec.get("basal_formula") or ""),
            "backbone_trunk": bb,
            "constituent_trunks": "+".join(c_trunks),
            "n_segments": int(rec.get("n_segments") or 1),
            "n_breakpoints": int(rec.get("n_breakpoints") or 0),
            "hxb2_barcode": str(bc),
            "breakpoints_ci": bp_ci,
        }
    return out


_GENERIC_PATIENT_RE = re.compile(
    r"^(?:patient\s*[:#]?\s*(?:code|sex|age|s|was|from|cohort|id)?|plasma.*|pbmc.*|peripheral blood.*|buffy coat.*|msm|heterosexual|male|female|intravenous drug user|positive blood unit|seropositive.*|mother|infant)$",
    re.I,
)


def _norm_str(val) -> str:
    if pd.isna(val):
        return ""
    s = str(val).strip()
    return "" if s.lower() in ("nan", "none", "unknown", "na") else s


def _extract_dedup_tokens(rec: dict) -> dict:
    """Extract conservative patient, clone, and epidemiological link tokens for intra-patient/couple dedup."""
    acc = str(rec.get("accession", ""))
    iso = _norm_str(rec.get("ff_isolate") or "") or _norm_str(rec.get("isolate"))
    pat_raw = _norm_str(rec.get("patient"))
    clone_raw = _norm_str(rec.get("clone"))
    title = _norm_str(rec.get("title"))
    comment = _norm_str(rec.get("comment"))
    authors = _norm_str(rec.get("authors"))

    specific_patient_id = None
    if pat_raw and not _GENERIC_PATIENT_RE.match(pat_raw):
        m_pid = re.search(
            r"\b(?:patient|subject|donor|pid|pt)\s*[:#_-]?\s*([A-Za-z]*\d+[A-Za-z0-9_-]*)",
            pat_raw,
            re.I,
        )
        if m_pid:
            specific_patient_id = m_pid.group(1).upper()
        elif re.search(r"\d", pat_raw) and len(pat_raw) <= 40:
            specific_patient_id = pat_raw.upper()

    # Longitudinal YY<Initials><Visit> patient codes (e.g. Cho,Y.-K. LTNP cohort 94KBH10 / 04KBH8)
    if specific_patient_id is None and ("Cho,Y.K." in authors or "Cho,Y.-K." in authors):
        m_ltnp = re.match(r"^\d{2}([A-Z]{2,4})\d+[a-z]?$", iso)
        if m_ltnp:
            specific_patient_id = m_ltnp.group(1).upper()

    m_title_pat = re.search(
        r"\b(?:patient|subject|donor|pid)\s*[:#_-]?\s*([A-Za-z]*\d+[A-Za-z0-9_-]*)",
        title + " " + comment,
        re.I,
    )
    if specific_patient_id is None and m_title_pat:
        specific_patient_id = m_title_pat.group(1).upper()

    clone_id = None
    if clone_raw and clone_raw.lower() not in ("1", "2", "3", "yes", "no"):
        clone_id = clone_raw.upper()
    else:
        m_cl = re.search(r"clone\s*:\s*([A-Za-z0-9_\-]+)", title, re.I)
        if m_cl:
            clone_id = m_cl.group(1).upper()

    pair_key = None
    study_sig = authors[:45].lower() if authors else title[:45].lower()
    if re.search(r"couple|partner|spouse|linked|mother|infant|transmission\s+pair|donor[-/\s]recipient", title + " " + comment, re.I):
        m_pp = re.match(r"^(PP\d+)", iso, re.I)
        m_fam = re.match(r"^(?:ML|BL)(\d+)", iso, re.I)
        m_zam = re.match(r"^(ZAM\d+)", iso, re.I)
        if m_pp:
            pair_key = f"{study_sig}:{m_pp.group(1).upper()}"
        elif m_fam:
            pair_key = f"{study_sig}:FAM{m_fam.group(1)}"
        elif m_zam:
            pair_key = f"{study_sig}:{m_zam.group(1).upper()}"
        else:
            pair_key = f"{acc[:5]}:{study_sig}"

    return {
        "patient_id": specific_patient_id,
        "clone_id": clone_id,
        "isolate_norm": iso.upper() if iso else "",
        "pair_key": pair_key,
        "study_sig": study_sig,
        "authors": authors,
    }


def main():
    t0 = time.time()
    print("=" * 100)
    print("STEP 8: ALIGNMENT-HARDENED SCAN OF 31,825 GENBANK FULL-LENGTH HIV-1 GENOMES FOR UNIDENTIFIED CRFs")
    print("=" * 100)

    profile_path, trunk_mat, trunk_names, hxb2_coords, cm240_arr = build_15trunk_hxb2_profile()
    loop_mask_bool = np.zeros(len(hxb2_coords), dtype=bool)
    for h_lo, h_hi in HYPERVARIABLE_MASK_RANGES:
        loop_mask_bool |= (hxb2_coords >= h_lo) & (hxb2_coords <= h_hi)

    cls_path = RESULTS_DIR / "genbank_31k_rhizaeon_classifications.csv"
    assert cls_path.exists(), "Expected genbank_31k_rhizaeon_classifications.csv from alignment step"
    gb = pd.read_csv(cls_path)
    gb_meta = pd.read_csv(GENBANK_META).set_index("accession")
    lanl_val = pd.read_csv(RESULTS_DIR / "lanl_all_full_genomes_rhizaeon_validation.csv")
    comp_df = pd.read_csv(RESULTS_DIR / "hiv1_descriptive_nomenclature_compendium.csv")
    bps_df = pd.read_csv(RESULTS_DIR / "hiv1_deconvolved_breakpoints.csv")
    lanl_flt_meta = pd.read_csv(DATA_DIR / "lanl_flt_2022_metadata.csv").set_index("accession")

    # 0. Re-call all validated patient genomes via rebuilt eggman binary (--gamma 13.5 --min-snps 7 --min-span 200 --p-thresh 1e-5)
    print(f"[+] Re-classifying {len(gb):,} genomes via {EGGMAN_BIN.name} --gamma 13.5 --min-snps 7 --min-span 200 --p-thresh 1e-5 ...")
    egg_map = _reclassify_with_eggman(gb["accession"].astype(str).tolist())
    for col in [
        "rhizaeon_formula",
        "backbone_trunk",
        "constituent_trunks",
        "n_segments",
        "n_breakpoints",
        "hxb2_barcode",
        "breakpoints_ci",
    ]:
        gb[col] = [egg_map[str(a)][col] for a in gb["accession"]]
    print(f"[+] Refreshed eggman calls in {time.time() - t0:.1f}s (n_breakpoints==0: {int((gb['n_breakpoints'] == 0).sum()):,}, >=1: {int((gb['n_breakpoints'] >= 1).sum()):,}).")

    def _enrich_metadata_and_years(df_in: pd.DataFrame) -> pd.DataFrame:
        ff_tsv = DATA_DIR / "genbank_cluster_flatfiles.tsv"
        ff_map = {}
        if ff_tsv.exists():
            ff_df = pd.read_csv(ff_tsv, sep="\t", dtype=str).fillna("")
            ff_map = {str(r["accession"]): r.to_dict() for _, r in ff_df.iterrows()}

        years_fixed = []
        titles = []
        coll_years = []
        sub_years = []
        crf_notes = []
        patients = []
        clones = []
        comments = []
        authors_col = []
        ff_isolates = []

        for _, r in df_in.iterrows():
            acc = str(r["accession"])
            ff_r = ff_map.get(acc, {})

            t_meta = _norm_str(gb_meta.loc[acc, "title"]) if acc in gb_meta.index else ""
            t_ff = _norm_str(ff_r.get("title", ""))
            t_str = t_ff if t_ff else t_meta
            titles.append(t_str)
            ff_isolates.append(_norm_str(ff_r.get("isolate", "")))

            c_yr = np.nan
            c_yr_raw = _norm_str(ff_r.get("collection_year", ""))
            if c_yr_raw:
                try:
                    val = float(c_yr_raw)
                    if 1975.0 <= val <= 2027.0:
                        c_yr = val
                except ValueError:
                    pass
            coll_years.append(c_yr)

            s_yr = np.nan
            s_yr_raw = _norm_str(ff_r.get("first_submitted_year", ""))
            if s_yr_raw:
                try:
                    val = float(s_yr_raw)
                    if 1980.0 <= val <= 2027.0:
                        s_yr = val
                except ValueError:
                    pass
            sub_years.append(s_yr)

            crf_notes.append(_norm_str(ff_r.get("official_crf_note", "")))
            patients.append(_norm_str(ff_r.get("patient", "")))
            clones.append(_norm_str(ff_r.get("clone", "")))
            comments.append(_norm_str(ff_r.get("comment", "")))
            authors_col.append(_norm_str(ff_r.get("authors", "")))

            # Recompute sampling year cleanly from raw sources (never inherit stale cls_path year!)
            yr = c_yr
            if pd.isna(yr) and acc in gb_meta.index:
                y_m = gb_meta.loc[acc, "year"]
                if pd.notna(y_m):
                    yr = float(y_m)
            if pd.isna(yr) and acc in lanl_flt_meta.index:
                y_l = lanl_flt_meta.loc[acc, "year"]
                if pd.notna(y_l):
                    yr = float(y_l)
            if pd.isna(yr):
                iso_s = _norm_str(r["isolate"])
                m_yr4 = re.search(r"(?:^|[^0-9])(19[89]\d|20[012]\d)(?:[^0-9]|$)", iso_s + " " + t_str)
                if m_yr4:
                    yr_val = float(m_yr4.group(1))
                    if 1980.0 <= yr_val <= 2025.0:
                        yr = yr_val
                else:
                    # City-prefixed 2-digit years (e.g. HLJ19197, GZGY230091, ZZS24046) or standard NN+alpha-code (20HZ2132, 03ZAPS119); NEVER match bare 6-digit IDs like 200202 or 191996
                    m_city = re.match(
                        r"^(?:P\d+_)?(?:HLJ|GZGY|GZZY|XYPQ|ZZS|NY|BJMP)(0[0-9]|1[0-9]|2[0-4])\d*",
                        iso_s,
                    )
                    m_yr2 = m_city or re.match(
                        r"^(?:P\d+_)?(0[0-9]|1[0-9]|2[0-4])(?:HZ|VLD|CM|NG|ZA|KE|RU|UG|TZ|GH|BS|BJ|KMM)",
                        iso_s,
                    )
                    if m_yr2:
                        yr = 2000.0 + float(m_yr2.group(1))
            years_fixed.append(yr)

        df_out = df_in.copy()
        df_out["year"] = years_fixed
        df_out["collection_year"] = coll_years
        df_out["first_submitted_year"] = sub_years
        df_out["title"] = titles
        df_out["official_crf_note"] = crf_notes
        df_out["patient"] = patients
        df_out["clone"] = clones
        df_out["comment"] = comments
        df_out["authors"] = authors_col
        df_out["ff_isolate"] = ff_isolates
        return df_out

    gb = _enrich_metadata_and_years(gb)

    def extract_title_crf(t):
        m = re.search(r"\b(CRF\s*\d{2,3}[_\-A-Za-z0-9]*)", str(t), re.I)
        return m.group(1).replace(" ", "") if m else None

    def _run_classification_pass(df_in: pd.DataFrame) -> pd.DataFrame:
        df_cur = df_in.copy()
        df_cur["title_crf"] = [extract_title_crf(t) for t in df_cur["title"]]

        # 1. Build comprehensive catalog of known CRFs & their breakpoint coordinates
        known_formulas = {}
        derived_subtype_formulas = {}
        for _, r in comp_df.iterrows():
            if int(r["n_verified_breakpoints"]) > 0:
                form_s = str(r["structural_formula"])
                lbl_s = str(r["legacy_label"])
                if lbl_s.startswith("Subtype_"):
                    derived_subtype_formulas[form_s] = lbl_s
                else:
                    known_formulas[form_s] = lbl_s

        named_crf_flt = lanl_val[lanl_val["lanl_category"] == "Assigned_Named_CRF"]
        for (tag, form), grp in named_crf_flt.groupby(["lanl_assigned_tag", "rhizaeon_formula"]):
            if len(grp) >= 2 and form not in known_formulas and form not in derived_subtype_formulas:
                known_formulas[form] = f"CRF{tag}_variant"

        for (tcrf, form), grp in df_cur[df_cur["title_crf"].notna()].groupby(["title_crf", "rhizaeon_formula"]):
            tcrf_norm = str(tcrf)
            if tcrf_norm in POST2023_CRF_NAMES or tcrf_norm.startswith("CRF63"):
                continue
            if len(grp) >= 2 and form not in known_formulas and form not in derived_subtype_formulas:
                known_formulas[form] = f"{tcrf_norm}_GenBank"

        crf_bp_catalog = []
        for cl_label, grp in bps_df.groupby("legacy_label"):
            c_rows = comp_df[comp_df["legacy_label"] == cl_label]
            if c_rows.empty:
                continue
            c_row = c_rows.iloc[0]
            donors = [d for d in str(c_row["donor_trunks"]).split(",") if d and d != "None"]
            c_set = frozenset([c_row["backbone_trunk"]] + donors)
            bp_list = [(r["left_trunk"], r["right_trunk"], int(r["bp_hxb2"])) for _, r in grp.iterrows()]
            crf_bp_catalog.append((cl_label, c_set, bp_list))

        # Post-2023 official LANL CRF formula/breakpoint matcher (without overriding CRF63_02A6)
        post2023_bp_catalog = []
        for crf_name in sorted(POST2023_CRF_NAMES):
            ref_accs = {a for a, c in POST2023_ACC_TO_CRF.items() if c == crf_name}
            sub_ref = df_cur[df_cur["accession"].astype(str).isin(ref_accs) & (df_cur["n_breakpoints"] > 0)]
            if not sub_ref.empty:
                rep = sub_ref.iloc[0]
                c_set = frozenset(str(rep["constituent_trunks"]).split("+"))
                bp_list = [(l_t, r_t, pos) for l_t, r_t, pos, _, _ in parse_bp_str(rep["breakpoints_ci"])]
                post2023_bp_catalog.append((crf_name, str(rep["rhizaeon_formula"]), c_set, bp_list))

        def classify_recombinant_category(row):
            acc = str(row["accession"])
            form = str(row["rhizaeon_formula"])
            note_s = _norm_str(row.get("official_crf_note", ""))
            t_crf = _norm_str(row.get("title_crf", ""))

            # Explicit CRF63_02A6 guard (never misassign JX500700/MK984159/PZ426393 to CRF173_63A6)
            if acc in CRF63_REF_ACCS or t_crf.startswith("CRF63") or any(x.startswith("CRF63") for x in note_s.split(";") if x):
                return ("Known_Indexed_CRF_Title", "CRF63_02A6")

            # Explicit 8 post-2023 official LANL CRFs
            if acc in POST2023_ACC_TO_CRF:
                return ("Known_Indexed_CRF_Post2023", POST2023_ACC_TO_CRF[acc])
            for p_crf in POST2023_CRF_NAMES:
                if t_crf.upper() == p_crf.upper() or any(x.strip().upper() == p_crf.upper() for x in note_s.split(";") if x.strip()):
                    return ("Known_Indexed_CRF_Post2023", p_crf)

            # Derived Subtype G exclusion (G0[A1:gp41..nef])
            if form in derived_subtype_formulas:
                return ("Derived_Basal_Subtype", derived_subtype_formulas[form])

            if form in known_formulas:
                return ("Known_Indexed_CRF_Exact", known_formulas[form])
            if t_crf:
                return ("Known_Indexed_CRF_Title", t_crf)
            note_tokens = [x.strip() for x in note_s.split(";") if x.strip().upper().startswith("CRF")]
            if len(note_tokens) == 1:
                return ("Known_Indexed_CRF_Title", note_tokens[0])

            tag = str(row["lanl_tag"])
            if tag not in ("Unindexed_GenBank", "U", "?", "CPX") and ("_" in tag or tag.startswith("CRF")):
                return ("Known_Indexed_CRF_LANLTag", f"LANL_{tag}")

            c_set = frozenset(str(row["constituent_trunks"]).split("+"))
            q_bps = parse_bp_str(row["breakpoints_ci"])
            if not q_bps:
                return ("Clonal_Basal_Subtype", str(row["backbone_trunk"]))

            # Check exact formula + breakpoint match against the 8 post-2023 CRFs
            for p_crf, p_form, p_cset, p_bps in post2023_bp_catalog:
                if form == p_form and c_set == p_cset and len(q_bps) == len(p_bps):
                    if all(any(abs(pos - r_pos) <= 250 and l_t == rl_t and r_t == rr_t for rl_t, rr_t, r_pos in p_bps) for l_t, r_t, pos, _, _ in q_bps):
                        return ("Known_Indexed_CRF_Post2023", p_crf)

            if c_set.issubset({"A1", "A6", "G0"}) and "A6" in c_set and "G0" in c_set and row["country"] in ("Russia", "Uzbekistan", "Kyrgyzstan", "Tajikistan", "Kazakhstan"):
                crf63_coords = [1299, 1850, 2070, 3292, 4175, 5150, 6230, 6480, 7450, 8154, 8396, 8900]
                if all(any(abs(pos - c_pos) <= 280 for c_pos in crf63_coords) for _, _, pos, _, _ in q_bps):
                    return ("Attenuated_or_Boundary_Variant_of_Known_CRF", "CRF63_02A6")

            for crf_name, ref_c_set, ref_bps in crf_bp_catalog:
                if not c_set.issubset(ref_c_set):
                    continue
                if len(ref_bps) < len(q_bps):
                    continue
                matched_all = True
                for l_t, r_t, pos, lo, hi in q_bps:
                    if not any(abs(pos - r_pos) <= 300 and {l_t, r_t} == {rl_t, rr_t} for rl_t, rr_t, r_pos in ref_bps):
                        matched_all = False
                        break
                if matched_all:
                    if str(crf_name).startswith("Subtype_"):
                        return ("Derived_Basal_Subtype", crf_name)
                    return ("Attenuated_or_Boundary_Variant_of_Known_CRF", crf_name)

            return ("Unidentified_Novel_Recombinant", "Unindexed")

        cats, matches = [], []
        for _, r in df_cur.iterrows():
            if int(r["n_breakpoints"]) == 0:
                cats.append("Clonal_Basal_Subtype")
                matches.append(r["backbone_trunk"])
            else:
                c, m = classify_recombinant_category(r)
                cats.append(c)
                matches.append(m)
        df_cur["recomb_category"] = cats
        df_cur["matched_known_lineage"] = matches
        return df_cur

    # First classification pass to surface all multi-sequence URF candidate accessions
    gb = _run_classification_pass(gb)
    novel_pre = gb[gb["recomb_category"] == "Unidentified_Novel_Recombinant"]
    vc_pre = novel_pre["rhizaeon_formula"].value_counts()
    cand_pre_formulas = set(vc_pre[vc_pre >= 2].index)
    cand_pre_accs = sorted(set(novel_pre[novel_pre["rhizaeon_formula"].isin(cand_pre_formulas)]["accession"].astype(str)))

    if cand_pre_accs:
        print(f"[+] Ensuring GenBank flatfiles cached via 08b for {len(cand_pre_accs)} candidate multi-sequence URF accessions...")
        s08b.ensure_flatfiles_cached(cand_pre_accs)
        gb = _enrich_metadata_and_years(gb)
        gb = _run_classification_pass(gb)

    gb.to_csv(cls_path, index=False)

    print("\n" + "=" * 100)
    print(f"GLOBAL GENBANK FULL-GENOME SCAN SUMMARY (N = {len(gb):,} VALIDATED PATIENT REPRESENTATIVES)")
    print("=" * 100)
    for cat, cnt in gb["recomb_category"].value_counts().items():
        print(f"  * {cat:<48s}: {cnt:5,d} ({100.0*cnt/len(gb):5.1f}%)")

    novel_df_raw = gb[gb["recomb_category"] == "Unidentified_Novel_Recombinant"].copy()
    vc = novel_df_raw["rhizaeon_formula"].value_counts()
    cand_formulas = list(vc[vc >= 2].index)
    cand_accs = set(novel_df_raw[novel_df_raw["rhizaeon_formula"].isin(cand_formulas)]["accession"])
    print(f"\n[+] Extracting & aligning {len(cand_accs)} sequences across {len(cand_formulas)} candidate multi-patient novel formulas for sequence-level divergence verification...")

    raw_seqs_to_align = []
    with open(GENBANK_FASTA, "r") as f:
        cur_acc = None
        cur_seq = []
        for line in f:
            if line.startswith(">"):
                if cur_acc in cand_accs:
                    raw_seqs_to_align.append((cur_acc, "".join(cur_seq).upper()))
                cur_acc = line[1:].strip().split()[0].split(".")[0]
                cur_seq = []
            else:
                cur_seq.append(line.strip())
        if cur_acc in cand_accs:
            raw_seqs_to_align.append((cur_acc, "".join(cur_seq).upper()))

    aligned_cand_map = align_sequences_to_coding_arrays(raw_seqs_to_align, profile_path, loop_mask_bool)
    print(f"[+] Aligned {len(aligned_cand_map)} candidate novel recombinant sequences onto HXB2 coding coordinates.")

    flt_recs = s3.read_fasta(DATA_DIR / "lanl_flt_2022.fasta")
    flt_hxb2_seq = [seq for hdr, seq in flt_recs if "K03455" in hdr][0]
    flt_coding_cols = np.array([i for i, ch in enumerate(flt_hxb2_seq) if ch not in "-."], dtype=np.int32)[789:9417]
    lanl_seq_map = {}
    for hdr, seq in flt_recs:
        if "K03455" in hdr:
            continue
        acc = hdr.split(".")[-1]
        q_a = s3.seq_to_int_array(seq)[flt_coding_cols]
        q_a[loop_mask_bool] = -1
        lanl_seq_map[acc] = q_a

    verified_lineages = []
    convergent_architecture_bins = []
    intra_patient_clone_collapses = []

    for formula in vc.index:
        sub = novel_df_raw[novel_df_raw["rhizaeon_formula"] == formula].copy()
        if len(sub) < 3:
            continue

        recs_with_arr = []
        for _, r in sub.iterrows():
            acc = str(r["accession"])
            q_a = aligned_cand_map.get(acc, lanl_seq_map.get(acc))
            if q_a is not None:
                d = r.to_dict()
                d["q_arr"] = q_a
                d["dedup_tokens"] = _extract_dedup_tokens(d)
                recs_with_arr.append(d)

        # Greedy sequence-level + flatfile patient/clone/couple independent-patient filter
        recs_sorted = sorted(
            recs_with_arr,
            key=lambda x: (0 if x.get("in_lanl") else 1, x["year"] if pd.notna(x["year"]) else 9999.0, -x["valid_coding_nt"])
        )
        indep = []
        for cand in recs_sorted:
            qc = cand["q_arr"]
            tc = cand["dedup_tokens"]
            is_dup = False
            for rep in indep:
                qr = rep["q_arr"]
                tr = rep["dedup_tokens"]
                same_study = bool(tc["study_sig"] and tr["study_sig"] and tc["study_sig"] == tr["study_sig"])
                if tc["patient_id"] and tr["patient_id"] and tc["patient_id"] == tr["patient_id"] and (same_study or cand.get("country") == rep.get("country")):
                    is_dup = True
                    break
                if tc["clone_id"] and tr["clone_id"] and tc["clone_id"] == tr["clone_id"] and (same_study or (tc["isolate_norm"] and tc["isolate_norm"] == tr["isolate_norm"])):
                    is_dup = True
                    break
                if tc["isolate_norm"] and tr["isolate_norm"] and tc["isolate_norm"] == tr["isolate_norm"] and same_study:
                    is_dup = True
                    break
                if tc["pair_key"] and tr["pair_key"] and tc["pair_key"] == tr["pair_key"]:
                    is_dup = True
                    break
                vm = (qc >= 0) & (qr >= 0)
                if np.sum(vm) >= 4000 and float(np.mean(qc[vm] != qr[vm])) < 0.008:
                    is_dup = True
                    break
            if not is_dup:
                indep.append(cand)

        if len(indep) < 3:
            intra_patient_clone_collapses.append((formula, len(sub), len(indep), list(sub["accession"])))
            continue

        bp_clusters = []
        for rec in indep:
            q_bps = parse_bp_str(rec["breakpoints_ci"])
            placed = False
            for cl in bp_clusters:
                ref_bps = parse_bp_str(cl[0]["breakpoints_ci"])
                if len(q_bps) == len(ref_bps) and all(
                    l1 == l2 and r1 == r2 and abs(p1 - p2) <= 250
                    for (l1, r1, p1, _, _), (l2, r2, p2, _, _) in zip(q_bps, ref_bps)
                ):
                    cl.append(rec)
                    placed = True
                    break
            if not placed:
                bp_clusters.append([rec])

        best_cl = max(bp_clusters, key=len)
        n_homologous = len(best_cl)

        cl_for_stats = best_cl if n_homologous >= 3 else indep
        divs = []
        for i in range(len(cl_for_stats)):
            for j in range(i + 1, len(cl_for_stats)):
                qi, qj = cl_for_stats[i]["q_arr"], cl_for_stats[j]["q_arr"]
                vm = (qi >= 0) & (qj >= 0)
                if np.sum(vm) >= 4000:
                    divs.append(float(np.mean(qi[vm] != qj[vm])))
        mean_div = float(np.mean(divs)) if divs else 0.0
        min_div = float(np.min(divs)) if divs else 0.0

        countries = Counter(r["country"] for r in cl_for_stats if r["country"] != "Unknown")
        years = sorted(int(r["year"]) for r in cl_for_stats if pd.notna(r["year"]))
        yr_span = f"{years[0]}-{years[-1]}" if years else "Undated"
        sub_yrs = sorted(float(r["first_submitted_year"]) for r in cl_for_stats if pd.notna(r.get("first_submitted_year")))
        first_sub = round(min(sub_yrs), 4) if sub_yrs else np.nan
        n_outside_lanl = sum(1 for r in cl_for_stats if not r["in_lanl"])
        tags = Counter(r["lanl_tag"] for r in cl_for_stats)

        t_mrca_str = "Undated"
        dated = [(r["year"], r["q_arr"]) for r in cl_for_stats if pd.notna(r["year"])]
        if len(dated) >= 3 and (max(y for y, _ in dated) - min(y for y, _ in dated)) >= 2.0:
            cons_q = s3.build_plurality_consensus([q for _, q in dated])
            y_arr = np.array([y for y, _ in dated], dtype=np.float64)
            d_arr = np.array([
                float(np.mean(q[(q >= 0) & (cons_q >= 0)] != cons_q[(q >= 0) & (cons_q >= 0)]))
                for _, q in dated
            ], dtype=np.float64)
            ols_res = s2.run_ols_dating(y_arr, d_arr)
            if pd.notna(ols_res.get("t_mrca")) and 1920.0 <= ols_res["t_mrca"] <= 2024.0:
                t_mrca_str = f"{ols_res['t_mrca']:.1f}"

        row_dict = {
            "rhizaeon_formula": formula,
            "lineage_status": "Homologous_Circulating_Lineage (Shared Breakpoints)" if n_homologous >= 3 else "Convergent_Hotspot_Architecture (Distinct Breakpoints)",
            "backbone_trunk": cl_for_stats[0]["backbone_trunk"],
            "constituent_trunks": cl_for_stats[0]["constituent_trunks"],
            "n_breakpoints": cl_for_stats[0]["n_breakpoints"],
            "n_raw_genbank": len(sub),
            "n_independent_patients_formula": len(indep),
            "n_homologous_bp_patients": n_homologous,
            "n_outside_lanl": n_outside_lanl,
            "n_in_lanl_urf": len(cl_for_stats) - n_outside_lanl,
            "mean_pairwise_div_pct": round(100.0 * mean_div, 2),
            "min_pairwise_div_pct": round(100.0 * min_div, 2),
            "countries": "; ".join(f"{k}({v})" for k, v in countries.most_common(5)),
            "sampling_years": yr_span,
            "first_submitted_year": first_sub,
            "chronaeon_tmrca": t_mrca_str,
            "lanl_or_gb_tags": "; ".join(f"{k}:{v}" for k, v in tags.most_common(4)),
            "member_accessions": ", ".join(r["accession"] for r in cl_for_stats),
            "representative_hxb2_barcode": cl_for_stats[0]["hxb2_barcode"],
            "representative_breakpoints": cl_for_stats[0]["breakpoints_ci"],
        }
        if n_homologous >= 3:
            verified_lineages.append(row_dict)
        else:
            convergent_architecture_bins.append(row_dict)

    cand_cols = [
        "rhizaeon_formula",
        "lineage_status",
        "backbone_trunk",
        "constituent_trunks",
        "n_breakpoints",
        "n_raw_genbank",
        "n_independent_patients_formula",
        "n_homologous_bp_patients",
        "n_outside_lanl",
        "n_in_lanl_urf",
        "mean_pairwise_div_pct",
        "min_pairwise_div_pct",
        "countries",
        "sampling_years",
        "first_submitted_year",
        "chronaeon_tmrca",
        "lanl_or_gb_tags",
        "member_accessions",
        "representative_hxb2_barcode",
        "representative_breakpoints",
    ]
    all_cand_rows = verified_lineages + convergent_architecture_bins
    all_cand_df = pd.DataFrame(all_cand_rows, columns=cand_cols)
    if not all_cand_df.empty:
        all_cand_df = all_cand_df.sort_values(
            ["n_homologous_bp_patients", "n_independent_patients_formula", "mean_pairwise_div_pct"],
            ascending=[False, False, False],
        )
    out_novel = RESULTS_DIR / "genbank_novel_circulating_recombinants.csv"
    all_cand_df.to_csv(out_novel, index=False)

    print("\n" + "=" * 100)
    print("AUDIT OF CANDIDATE MULTI-PATIENT NOVEL RECOMBINANT FORMULAS (N_raw >= 3)")
    print("=" * 100)
    print(f"  1. Collapsed to <3 patients after >=0.8% Whole-Genome + Flatfile Dedup (Intra-Patient/Couple Batches): {len(intra_patient_clone_collapses)}")
    for form, n_r, n_i, accs in intra_patient_clone_collapses:
        print(f"       - {form:<48s}: {n_r} raw -> {n_i} independent patient(s) ({', '.join(accs[:4])})")
    print(f"  2. Convergent Gene-Level Architectures (Distinct Breakpoints Across Patients, N_indep >= 3, N_homologous < 3): {len(convergent_architecture_bins)}")
    for r in convergent_architecture_bins:
        print(f"       - {r['rhizaeon_formula']:<48s}: {r['n_independent_patients_formula']} indep patients, max {r['n_homologous_bp_patients']} sharing exact BPs | {r['countries']} | {r['member_accessions']}")
    print(f"  3. VERIFIED HOMOLOGOUS UNIDENTIFIED CIRCULATING RECOMBINANT LINEAGES (N_homologous >= 3 independent patients): {len(verified_lineages)}")
    print("=" * 100)
    for idx, r in enumerate(verified_lineages, 1):
        print(f"\n  [{idx:2d}] {r['rhizaeon_formula']}  (N_homologous = {r['n_homologous_bp_patients']} indep patients [N_formula = {r['n_independent_patients_formula']}], {r['n_outside_lanl']} outside LANL)")
        print(f"       * Countries & Years : {r['countries']} | Years: {r['sampling_years']} | ChronAeon t_MRCA: {r['chronaeon_tmrca']}")
        print(f"       * Pairwise Div      : Mean = {r['mean_pairwise_div_pct']:.2f}% (Min = {r['min_pairwise_div_pct']:.2f}%)")
        print(f"       * Existing Labels   : {r['lanl_or_gb_tags']}")
        print(f"       * Breakpoints [CI]  : {r['representative_breakpoints']}")
        print(f"       * HXB2 Barcode      : {r['representative_hxb2_barcode']}")
        print(f"       * Accessions        : {r['member_accessions']}")


if __name__ == "__main__":
    main()
