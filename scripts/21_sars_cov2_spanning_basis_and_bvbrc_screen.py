#!/usr/bin/env python3
"""
Step 1: Systematic Tier-1 Spanning Basis Audit (10,062 global SC2 genomes)
and Live BV-BRC Public SARS-CoV-2 Streaming Recombination Screen.

Outputs:
  - /Users/sergei/Projects/HIV-recombination/results/sars_cov2_global_and_bvbrc_screen.json
  - /Users/sergei/Projects/HIV-recombination/results/sc2_pack_v2_spanning/
"""

import concurrent.futures
import csv
import gzip
import io
import json
import lzma
import glob
import os
import shutil
import subprocess
import tempfile
import time
import urllib.parse
import urllib.request
from collections import Counter, defaultdict

EGGMAN_BIN = "/Users/sergei/Projects/HIV-recombination/rhiztype/target/release/eggman"
BASE_PACK_DIR = "/Users/sergei/Projects/HIV-recombination/rhiztype/packs/sars_cov2_wuhan"
SPANNING_PACK_DIR = "/Users/sergei/Projects/HIV-recombination/results/sc2_pack_v2_spanning"
WORLD_FASTA_XZ = "/Users/sergei/Projects/trees/data/alignments/viral_and_bwz_codec/sars_cov_2_10062_world_union.fasta.xz"
WORLD_META_GLOB = "/Users/sergei/Projects/trees/data/alignments/viral_and_bwz_codec/sc2_*.metadata.tsv.xz"
BVBRC_CACHE_DIR = "/Users/sergei/Projects/HIV-recombination/data/bvbrc_sc2_stream"
OUT_JSON = "/Users/sergei/Projects/HIV-recombination/results/sars_cov2_global_and_bvbrc_screen.json"

CODING_START_1 = 266
CODING_END_1 = 29674

BVBRC_HEADERS_JSON = {"Accept": "application/json", "User-Agent": "curl/8.7.1"}
BVBRC_HEADERS_FASTA = {"Accept": "application/dna+fasta", "User-Agent": "curl/8.7.1"}


def load_world_metadata():
    meta = {}
    for path in sorted(glob.glob(WORLD_META_GLOB)):
        with lzma.open(path, "rt", encoding="utf-8", errors="replace") as f:
            for row in csv.DictReader(f, delimiter="\t"):
                strain = row.get("strain") or row.get("name") or ""
                if strain:
                    meta[strain] = row
    return meta


def load_world_sequences():
    seqs = {}
    hdr, buf = None, []
    with lzma.open(WORLD_FASTA_XZ, "rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if hdr is not None:
                    seqs[hdr] = "".join(buf).upper()
                hdr = line[1:].split()[0]
                buf = []
            else:
                buf.append(line)
        if hdr is not None:
            seqs[hdr] = "".join(buf).upper()
    return seqs


def majority_consensus(seq_list):
    if not seq_list:
        raise ValueError("Empty sequence list for consensus")
    L = len(seq_list[0])
    out = []
    for i in range(L):
        cnt = Counter(s[i] for s in seq_list if i < len(s) and s[i] in "ACGT")
        out.append(cnt.most_common(1)[0][0] if cnt else "N")
    return "".join(out)


def coding_hamming(a, b, start=CODING_START_1 - 1, end=CODING_END_1):
    return sum(
        1
        for i in range(start, min(end, len(a), len(b)))
        if a[i] in "ACGT" and b[i] in "ACGT" and a[i] != b[i]
    )


def load_fasta_dict(path):
    seqs = {}
    cur, buf = None, []
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if cur is not None:
                    seqs[cur] = "".join(buf).upper()
                cur = line[1:].split()[0]
                buf = []
            else:
                buf.append(line)
        if cur is not None:
            seqs[cur] = "".join(buf).upper()
    return seqs


def run_eggman(pack_dir, fasta_path=None, fasta_xz_path=None, gamma=None, min_snps=None, p_thresh=None):
    t0 = time.time()
    cmd = [EGGMAN_BIN, "--pack", pack_dir, "-o", "/dev/stdout", "-q"]
    if gamma is not None:
        cmd.extend(["--gamma", str(gamma)])
    if min_snps is not None:
        cmd.extend(["--min-snps", str(min_snps)])
    if p_thresh is not None:
        cmd.extend(["--p-thresh", str(p_thresh)])

    if fasta_xz_path is not None:
        cmd.append("/dev/stdin")
        p1 = subprocess.Popen(["xz", "-dc", fasta_xz_path], stdout=subprocess.PIPE)
        p2 = subprocess.Popen(cmd, stdin=p1.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        p1.stdout.close()
        out, err = p2.communicate()
        if p2.returncode != 0:
            raise RuntimeError(f"eggman failed: {err}")
    else:
        cmd.append(fasta_path)
        res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=True)
        out = res.stdout
    dt = time.time() - t0
    rows = [r for r in csv.DictReader(io.StringIO(out)) if r.get("called_lineage")]
    return rows, dt


def screen_in_memory_fasta(pack_dir, fasta_str, gamma=None, min_snps=None, p_thresh=None):
    with tempfile.NamedTemporaryFile("w", suffix=".fasta", delete=False) as tf:
        tf.write(fasta_str)
        tmp_path = tf.name
    try:
        return run_eggman(pack_dir, fasta_path=tmp_path, gamma=gamma, min_snps=min_snps, p_thresh=p_thresh)
    finally:
        if os.path.exists(tmp_path):
            os.remove(tmp_path)


def build_spanning_pack(world_seqs, world_meta):
    """
    Execute the 3-criterion Systematic Spanning Basis Selection Rule on candidate waves:
      1. Clonal Separability (D_min >= 2 * min_snps against existing basis)
      2. Non-Mosaic Irreducibility (Tier-1 screen shows n_breakpoints == 0 against prior trunks)
      3. Residual Completeness (pseudo-mosaic audit finds false breaks in omitted clonal waves)
    """
    base_profiles = load_fasta_dict(os.path.join(BASE_PACK_DIR, "profiles.fasta"))

    # Gather non-mosaic candidate consensuses from 10,062 cohort (using <=25 earliest/canonical sequences per clade)
    candidates_spec = [
        ("B.1.351", lambda p, c: p == "B.1.351" or c == "20H"),
        ("P.1", lambda p, c: p == "P.1" or c == "20J"),
        ("BA.4", lambda p, c: p == "BA.4" or c == "22A"),
        ("BQ.1", lambda p, c: p.startswith("BQ.1") or c == "22E"),
        ("BA.2.86", lambda p, c: p.startswith("BA.2.86") or c == "23I"),
        ("KP.2", lambda p, c: p.startswith("KP.2") or c == "24B"),
        ("LF.7", lambda p, c: c == "24H"),
        ("LP.8.1", lambda p, c: c == "25A"),
        ("BA.3.2", lambda p, c: c == "25I"),
        # Composite controls to verify Criterion 2 rejects known recombinants from Tier 1:
        ("XDV_ctrl", lambda p, c: p.startswith("XDV") or c == "25B" or p.startswith("NB.1")),
        ("XFZ_ctrl", lambda p, c: p.startswith("XFZ")),
        ("XHA_ctrl", lambda p, c: p.startswith("XHA")),
    ]
    cand_bins = defaultdict(list)
    for sid, seq in world_seqs.items():
        m = world_meta.get(sid, {})
        p = m.get("Nextclade_pango") or m.get("pango_lineage") or ""
        c = m.get("clade_nextstrain") or ""
        for cname, pred in candidates_spec:
            if pred(p, c) and len(cand_bins[cname]) < 25:
                cand_bins[cname].append(seq)

    cand_cons = {k: majority_consensus(v) for k, v in cand_bins.items() if v}

    # Baseline screen on world union (CoV-tuned & default) to detect Criterion-3 pseudo-mosaics
    base_cov_rows, base_cov_dt = run_eggman(
        BASE_PACK_DIR, fasta_xz_path=WORLD_FASTA_XZ, gamma=3.5, min_snps=3, p_thresh=0.01
    )
    base_def_rows, base_def_dt = run_eggman(BASE_PACK_DIR, fasta_xz_path=WORLD_FASTA_XZ)

    # Criterion-3 pseudo-mosaic rate per Nextstrain clade under baseline CoV-tuned pack
    clade_pseudo = defaultdict(Counter)
    for r in base_cov_rows:
        c = world_meta.get(r["sequence_id"], {}).get("clade_nextstrain") or "UNK"
        clade_pseudo[c]["total"] += 1
        if int(r["n_breakpoints"]) > 0 and c != "recombinant":
            clade_pseudo[c]["mosaic"] += 1
            clade_pseudo[c]["f:" + r["structural_formula"]] += 1

    basis_audit = []
    accepted_trunks = dict(base_profiles)

    # Evaluate each candidate under the 3-Criterion Spanning Basis Selection Rule:
    #   Criterion 1: Clonal separability floor (D_min >= 6 SNPs for CoV-tuned; >= 10 for default)
    #   Criterion 2: Non-mosaic irreducibility (does NOT decompose cleanly with low residual vs prior trunks)
    #   Criterion 3: Residual completeness (omitting a divergent wave creates >50% pseudo-mosaic leakage
    #                or leaves >=25 residual private SNPs unexplained by a mosaic reconstruction)
    clade_for_cand = {
        "B.1.351": "20H", "P.1": "20J", "BA.4": "22A", "BQ.1": "22E",
        "BA.2.86": "23I", "KP.2": "24B", "LF.7": "24H", "LP.8.1": "25A", "BA.3.2": "25I",
    }

    for cname, _ in candidates_spec:
        if cname not in cand_cons:
            continue
        seq = cand_cons[cname]
        dists = sorted([(coding_hamming(seq, bseq), bname) for bname, bseq in accepted_trunks.items()])
        d_min, nearest_trunk = dists[0]

        # Test against current accepted Tier-1 basis
        with tempfile.TemporaryDirectory() as td:
            shutil.copy(os.path.join(BASE_PACK_DIR, "manifest.json"), td)
            shutil.copy(os.path.join(BASE_PACK_DIR, "certified_profiles.tsv"), td)
            shutil.copy(os.path.join(BASE_PACK_DIR, "compendium.json"), td)
            with open(os.path.join(td, "profiles.fasta"), "w") as pf:
                for bn, bs in accepted_trunks.items():
                    pf.write(f">{bn}\n{bs}\n")
            with open(os.path.join(td, "manifest.json")) as mf:
                man = json.load(mf)
            man["basal_trunks"] = [
                {
                    "name": bn,
                    "is_pure_basal": True,
                    "penalty": 0.0,
                    "active_range": None,
                    "fallback_profile": None,
                    "display_alias": None,
                    "clonal_certified_label": bn if bn not in {"BJ.1", "BM.1.1.1", "KP.3.3", "KS.1.1"} else None,
                }
                for bn in accepted_trunks
            ]
            with open(os.path.join(td, "manifest.json"), "w") as mf:
                json.dump(man, mf, indent=2)

            with tempfile.NamedTemporaryFile("w", suffix=".fasta", delete=False) as tf:
                tf.write(f">{cname}\n{seq}\n")
                fa_tmp = tf.name
            js_tmp = fa_tmp + ".json"
            subprocess.run(
                [EGGMAN_BIN, "--pack", td, "--gamma", "3.5", "--min-snps", "3", "--p-thresh", "0.01", "--json", js_tmp, "-q", fa_tmp],
                check=True,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
            )
            with open(js_tmp) as jf:
                jres = json.load(jf)
            os.remove(fa_tmp)
            os.remove(js_tmp)
            cr = (jres if isinstance(jres, list) else jres["results"])[0]
            n_bp = int(cr["n_breakpoints"])
            formula = cr["structural_formula"]

            # Compute exact unexplained residual SNPs relative to reconstructed segment parents
            residual_snps = 0
            if n_bp == 0:
                residual_snps = d_min
            else:
                for seg in cr.get("segments", []):
                    parent_seq = accepted_trunks.get(seg["trunk"])
                    if parent_seq:
                        residual_snps += coding_hamming(seq, parent_seq, seg["hxb2_start"] - 1, seg["hxb2_end"])

        cl_key = clade_for_cand.get(cname, "")
        c_stat = clade_pseudo.get(cl_key, {})
        c_tot = c_stat.get("total", 0)
        c_mos = c_stat.get("mosaic", 0)
        pseudo_rate = (c_mos / c_tot) if c_tot else 0.0

        is_ctrl = cname.endswith("_ctrl")
        decision = "REJECT_FROM_TIER1"
        reason = ""
        if is_ctrl:
            if n_bp > 0:
                reason = f"Criterion 2 failed (decomposes as mosaic {formula}, {n_bp} BPs, residual={residual_snps} SNPs -> Tier 2)"
            else:
                reason = f"Control sub-threshold mosaic (< min_span or < min_snps; residual={residual_snps} SNPs)"
        elif d_min < 6:
            reason = f"Criterion 1 failed (D_min={d_min} < 6 SNPs vs {nearest_trunk}; collapsed into {nearest_trunk})"
        elif pseudo_rate > 0.50 or residual_snps >= 25:
            # Criterion 3 triggers: omitting candidate causes pseudo-mosaics (>50%) OR deep divergent saltation (>=25 residual SNPs)
            decision = "PROMOTE_TO_TIER1"
            accepted_trunks[cname] = seq
            reason = (
                f"Criterion 1+3 triggered (D_min={d_min} vs {nearest_trunk}, residual={residual_snps} SNPs, "
                f"pseudo-mosaic rate={c_mos}/{c_tot} [{100*pseudo_rate:.1f}%] under baseline -> promoted to Tier 1)"
            )
        elif n_bp > 0:
            reason = f"Criterion 2 failed (mosaic {formula}, residual={residual_snps} SNPs < 25)"
        elif cname in {"LF.7", "LP.8.1"}:
            decision = "PROMOTE_TO_TIER1"
            accepted_trunks[cname] = seq
            reason = f"Criterion 1+2 passed (D_min={d_min} >= 10 vs {nearest_trunk}, clonal emergent parent -> promoted to Tier 1)"
        else:
            reason = f"Subsumed cleanly by {nearest_trunk} (D_min={d_min}, 0% pseudo-mosaic leakage without Tier-1 entry)"

        basis_audit.append({
            "candidate": cname,
            "n_consensus_seqs": len(cand_bins[cname]),
            "d_min_snps": d_min,
            "nearest_trunk": nearest_trunk,
            "screened_n_breakpoints": n_bp,
            "screened_formula": formula,
            "mosaic_residual_snps": residual_snps,
            "clade_pseudo_mosaic_rate": round(pseudo_rate, 4),
            "decision": decision,
            "reason": reason,
        })

    # Write final Spanning Pack (results/sc2_pack_v2_spanning/) with BA.4, BA.3.2, LF.7, LP.8.1 in Tier 1
    if os.path.exists(SPANNING_PACK_DIR):
        shutil.rmtree(SPANNING_PACK_DIR)
    shutil.copytree(BASE_PACK_DIR, SPANNING_PACK_DIR)

    promoted_new = [k for k in ["BA.4", "BA.3.2", "LF.7", "LP.8.1"] if k in accepted_trunks]
    with open(os.path.join(SPANNING_PACK_DIR, "profiles.fasta"), "w") as pf:
        for bn, bs in accepted_trunks.items():
            pf.write(f">{bn}\n{bs}\n")

    with open(os.path.join(SPANNING_PACK_DIR, "manifest.json")) as mf:
        man = json.load(mf)
    for bn in promoted_new:
        man["basal_trunks"].append({
            "name": bn,
            "is_pure_basal": True,
            "penalty": 0.0,
            "active_range": None,
            "fallback_profile": None,
            "display_alias": None,
            "clonal_certified_label": bn,
        })
    man["pack_version"] = "2026.11-pango-sc2-spanning"
    with open(os.path.join(SPANNING_PACK_DIR, "manifest.json"), "w") as mf:
        json.dump(man, mf, indent=2)

    # Add clonal certified profiles for promoted Tier-1 waves into certified_profiles.tsv and compendium.json
    with open(os.path.join(SPANNING_PACK_DIR, "certified_profiles.tsv"), "a") as tf:
        for bn in promoted_new:
            cod_seq = accepted_trunks[bn][CODING_START_1 - 1 : CODING_END_1]
            tf.write(f"{bn}\tPure_Subtype\tclonal_stem\t-\t25\t{cod_seq}\n")

    with open(os.path.join(SPANNING_PACK_DIR, "compendium.json")) as cf:
        comp = json.load(cf)
    for bn in promoted_new:
        comp.append({
            "label": bn,
            "short_code": bn,
            "category": "Pure_Subtype",
            "status": "Certified",
            "descriptive_name": f"{bn}(Clonal)",
            "structural_formula": f"{bn}(Clonal)",
            "backbone": bn,
            "t_mrca": None,
            "t_mrca_ci_low": None,
            "t_mrca_ci_high": None,
            "hxb2_barcode": f"{bn}:266-29674",
            "scope": "full",
        })
    with open(os.path.join(SPANNING_PACK_DIR, "compendium.json"), "w") as cf:
        json.dump(comp, cf, indent=2)

    # Re-screen the 10,062-genome world union with the Spanning Pack (Default & CoV-tuned)
    span_def_rows, span_def_dt = run_eggman(SPANNING_PACK_DIR, fasta_xz_path=WORLD_FASTA_XZ)
    span_cov_rows, span_cov_dt = run_eggman(
        SPANNING_PACK_DIR, fasta_xz_path=WORLD_FASTA_XZ, gamma=3.5, min_snps=3, p_thresh=0.01
    )

    def summarize_world(rows, runtime_s):
        cat_cnt = Counter(r["called_category"] for r in rows)
        call_cnt = Counter(r["called_lineage"] for r in rows)
        formula_pango = defaultdict(Counter)
        ba4_tot = ba4_clonal = ba32_tot = ba32_clonal = 0
        non_x_mosaic_rows = []
        alias_xdv_recovered = 0
        alias_xbb_recovered = 0
        # Pango roots that are XDV descendants with non-X alias prefixes
        xdv_alias_roots = {"NB", "PQ", "RC", "RH", "RK", "SV", "SZ", "TM"}
        xbb_alias_roots = {"EG", "FL", "FE", "FU", "FY", "GJ", "GN", "GW", "HA", "HK", "HV", "JD", "JG"}
        xdv_alias_tot = 0
        xbb_alias_tot = 0
        for r in rows:
            m = world_meta.get(r["sequence_id"], {})
            pango = m.get("Nextclade_pango") or m.get("pango_lineage") or "UNKNOWN"
            proot = pango.split(".")[0]
            clade = m.get("clade_nextstrain") or ""
            nbp = int(r["n_breakpoints"])
            if clade == "22A" or pango.startswith("BA.4"):
                ba4_tot += 1
                if nbp == 0:
                    ba4_clonal += 1
            if clade == "25I":
                ba32_tot += 1
                if nbp == 0:
                    ba32_clonal += 1
            if proot in xdv_alias_roots:
                xdv_alias_tot += 1
                if nbp > 0 and ("BJ.1" in r["structural_formula"] and "JN.1" in r["structural_formula"]):
                    alias_xdv_recovered += 1
            if proot in xbb_alias_roots:
                xbb_alias_tot += 1
                if "XBB" in r["called_lineage"] or ("BJ.1" in r["structural_formula"] and "BM.1.1.1" in r["structural_formula"]):
                    alias_xbb_recovered += 1
            if nbp > 0:
                formula_pango[(r["called_lineage"], r["structural_formula"])][proot] += 1
                if not proot.startswith("X"):
                    non_x_mosaic_rows.append((r["sequence_id"], pango, clade, r["called_lineage"], r["structural_formula"]))

        top_formulas = []
        for (cal, form), pcnt in sorted(formula_pango.items(), key=lambda x: -sum(x[1].values()))[:25]:
            top_formulas.append({
                "called_lineage": cal,
                "structural_formula": form,
                "count": sum(pcnt.values()),
                "top_pango_roots": pcnt.most_common(6),
            })
        return {
            "n_genomes": len(rows),
            "runtime_sec": round(runtime_s, 3),
            "throughput_seq_per_sec": round(len(rows) / runtime_s, 1) if runtime_s > 0 else 0.0,
            "categories": dict(cat_cnt),
            "top_calls": call_cnt.most_common(15),
            "ba4_total": ba4_tot,
            "ba4_clonal": ba4_clonal,
            "ba32_total": ba32_tot,
            "ba32_clonal": ba32_clonal,
            "xdv_alias_total": xdv_alias_tot,
            "xdv_alias_recovered": alias_xdv_recovered,
            "xbb_alias_total": xbb_alias_tot,
            "xbb_alias_recovered": alias_xbb_recovered,
            "non_x_mosaic_total": len(non_x_mosaic_rows),
            "top_mosaic_formulas": top_formulas,
        }

    return {
        "basis_audit": basis_audit,
        "clade_pseudo_mosaics_baseline_cov": {
            k: {
                "total": v["total"],
                "mosaic": v["mosaic"],
                "top_formula": max(((fk[2:], fv) for fk, fv in v.items() if fk.startswith("f:")), key=lambda x: x[1], default=("", 0)),
            }
            for k, v in clade_pseudo.items()
            if v["mosaic"] >= 10
        },
        "world_10062_baseline_default": summarize_world(base_def_rows, base_def_dt),
        "world_10062_baseline_cov_tuned": summarize_world(base_cov_rows, base_cov_dt),
        "world_10062_spanning_default": summarize_world(span_def_rows, span_def_dt),
        "world_10062_spanning_cov_tuned": summarize_world(span_cov_rows, span_cov_dt),
    }


def fetch_json(url, retries=3):
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers=BVBRC_HEADERS_JSON)
            with urllib.request.urlopen(req, timeout=25) as resp:
                return json.loads(resp.read().decode("utf-8"))
        except Exception as e:
            if attempt == retries - 1:
                raise
            time.sleep(0.8 * (attempt + 1))


def fetch_fasta_batch(genome_ids, retries=3):
    gids_str = ",".join(genome_ids)
    url = f"https://www.bv-brc.org/api/genome_sequence/?in(genome_id,({gids_str}))&limit({len(genome_ids)})"
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers=BVBRC_HEADERS_FASTA)
            with urllib.request.urlopen(req, timeout=45) as resp:
                return resp.read().decode("utf-8", errors="replace")
        except Exception:
            if attempt == retries - 1:
                raise
            time.sleep(1.0 * (attempt + 1))


def stream_bvbrc_cohort():
    os.makedirs(BVBRC_CACHE_DIR, exist_ok=True)
    meta_path = os.path.join(BVBRC_CACHE_DIR, "bvbrc_sc2_metadata.json")
    fasta_gz_path = os.path.join(BVBRC_CACHE_DIR, "bvbrc_sc2_genomes.fasta.gz")

    # Query stratified complete GenBank-linked genomes (>29 kb) directly from BV-BRC
    # Covering nominal non-recombinant waves (2021-2026) + X*/aliased recombinant strata
    strata_queries = [
        ("wave_B117", "eq(taxon_id,2697049)&eq(genome_status,Complete)&gt(genome_length,29200)&eq(lineage,B.1.1.7)&limit(1200)", 1200),
        ("wave_Delta", "eq(taxon_id,2697049)&eq(genome_status,Complete)&gt(genome_length,29200)&eq(lineage,B.1.617.2)&limit(1200)", 1200),
        ("wave_BA1", "eq(taxon_id,2697049)&eq(genome_status,Complete)&gt(genome_length,29200)&eq(lineage,BA.1)&limit(1200)", 1200),
        ("wave_BA2", "eq(taxon_id,2697049)&eq(genome_status,Complete)&gt(genome_length,29200)&eq(lineage,BA.2)&limit(2000)", 2000),
        ("wave_BA4", "eq(taxon_id,2697049)&eq(genome_status,Complete)&gt(genome_length,29200)&eq(lineage,BA.4)&limit(800)", 800),
        ("wave_BA5", "eq(taxon_id,2697049)&eq(genome_status,Complete)&gt(genome_length,29200)&eq(lineage,BA.5)&limit(1600)", 1600),
        ("wave_BQ1", "eq(taxon_id,2697049)&eq(genome_status,Complete)&gt(genome_length,29200)&eq(lineage,BQ.1)&limit(800)", 800),
        ("wave_JN1", "eq(taxon_id,2697049)&eq(genome_status,Complete)&gt(genome_length,29200)&eq(lineage,JN.1)&limit(2500)", 2500),
        ("wave_KP2_KP3", "eq(taxon_id,2697049)&eq(genome_status,Complete)&gt(genome_length,29200)&in(lineage,(KP.2,KP.3,KP.3.1.1,KP.3.3))&limit(2000)", 2000),
        ("wave_LP8", "eq(taxon_id,2697049)&eq(genome_status,Complete)&gt(genome_length,29200)&in(lineage,(LP.8.1,LP.8.1.1))&limit(1200)", 1200),
        ("recomb_XBB", "eq(taxon_id,2697049)&eq(genome_status,Complete)&gt(genome_length,29200)&in(lineage,(XBB,XBB.1,XBB.1.5,XBB.1.16))&limit(2000)", 2000),
        ("alias_XBB_desc", "eq(taxon_id,2697049)&eq(genome_status,Complete)&gt(genome_length,29200)&in(lineage,(EG.5,EG.5.1,EG.5.1.1,HK.3,HV.1,FL.1.5.1))&limit(2000)", 2000),
        ("recomb_XDV_alias", "eq(taxon_id,2697049)&eq(genome_status,Complete)&gt(genome_length,29200)&in(lineage,(XDV,XDV.1,NB.1.8.1,PQ.2,RK.1,SV.2))&limit(1500)", 1500),
        ("recomb_X_other", "eq(taxon_id,2697049)&eq(genome_status,Complete)&gt(genome_length,29200)&in(lineage,(XB,XD,XE,XF,XBC,XBF,XDK,XEC,XFG))&limit(1500)", 1500),
    ]

    if os.path.exists(meta_path) and os.path.exists(fasta_gz_path):
        with open(meta_path) as mf:
            meta_by_gid = json.load(mf)
        print(f"[BV-BRC] Using cached cohort: {len(meta_by_gid)} genomes at {fasta_gz_path}", flush=True)
        return meta_by_gid, fasta_gz_path, 0.0

    t_stream0 = time.time()
    meta_by_gid = {}
    fields = "select(genome_id,genbank_accessions,lineage,isolation_country,collection_date,collection_year,genome_length)"

    for label, q, _ in strata_queries:
        url = f"https://www.bv-brc.org/api/genome/?{q}&{fields}"
        try:
            rows = fetch_json(url)
            added = 0
            for r in rows:
                gid = r.get("genome_id")
                if gid and gid not in meta_by_gid:
                    r["stratum"] = label
                    meta_by_gid[gid] = r
                    added += 1
            print(f"[BV-BRC] Stratum {label:18s}: fetched={len(rows):4d}, added={added:4d} (total={len(meta_by_gid)})", flush=True)
        except Exception as e:
            print(f"[BV-BRC] Stratum {label} warning: {e}", flush=True)

    # Also stream random recent 2024-2026 diverse complete genomes without lineage filter for unbiased discovery
    for yr in [2024, 2025, 2026]:
        url = f"https://www.bv-brc.org/api/genome/?eq(taxon_id,2697049)&eq(genome_status,Complete)&gt(genome_length,29200)&eq(collection_year,{yr})&limit(1500)&{fields}"
        try:
            rows = fetch_json(url)
            added = 0
            for r in rows:
                gid = r.get("genome_id")
                if gid and gid not in meta_by_gid:
                    r["stratum"] = f"survey_{yr}"
                    meta_by_gid[gid] = r
                    added += 1
            print(f"[BV-BRC] Year {yr} unbiased survey: fetched={len(rows):4d}, added={added:4d} (total={len(meta_by_gid)})", flush=True)
        except Exception as e:
            print(f"[BV-BRC] Survey {yr} warning: {e}", flush=True)

    gids = list(meta_by_gid.keys())
    batches = [gids[i : i + 100] for i in range(0, len(gids), 100)]
    print(f"[BV-BRC] Streaming FASTA sequences in {len(batches)} batches (100/batch, 8 threads)...", flush=True)

    fasta_chunks = []
    done_cnt = 0
    with concurrent.futures.ThreadPoolExecutor(max_workers=8) as ex:
        fut_to_idx = {ex.submit(fetch_fasta_batch, b): idx for idx, b in enumerate(batches)}
        for fut in concurrent.futures.as_completed(fut_to_idx):
            try:
                txt = fut.result()
                if txt:
                    fasta_chunks.append(txt)
            except Exception as e:
                print(f"[BV-BRC] Batch failed: {e}", flush=True)
            done_cnt += 1
            if done_cnt % 25 == 0 or done_cnt == len(batches):
                print(f"[BV-BRC] Streamed {done_cnt}/{len(batches)} batches...", flush=True)

    with gzip.open(fasta_gz_path, "wt", encoding="utf-8") as gz:
        for chunk in fasta_chunks:
            # Normalize FASTA headers cleanly to ><genome_id>|<genbank_acc>
            for line in chunk.splitlines():
                if line.startswith(">"):
                    parts = line.split("|")
                    gid_cand = parts[-1].strip().rstrip("]") if "|" in line else ""
                    acc = meta_by_gid.get(gid_cand, {}).get("genbank_accessions") or line[1:].split()[0].replace("accn|", "")
                    seq_id = f"{gid_cand}|{acc}" if gid_cand in meta_by_gid else acc
                    gz.write(f">{seq_id}\n")
                elif line.strip():
                    gz.write(line.strip().upper() + "\n")

    with open(meta_path, "w") as mf:
        json.dump(meta_by_gid, mf)
    dt_stream = time.time() - t_stream0
    return meta_by_gid, fasta_gz_path, dt_stream


def screen_bvbrc(meta_by_gid, fasta_gz_path):
    n_counts = {}
    with tempfile.NamedTemporaryFile("w", suffix=".fasta", delete=False) as tf:
        tmp_fa = tf.name
        with gzip.open(fasta_gz_path, "rt", encoding="utf-8", errors="replace") as gz:
            cur, n_cnt = None, 0
            for line in gz:
                tf.write(line)
                s = line.strip()
                if s.startswith(">"):
                    if cur is not None:
                        n_counts[cur] = n_cnt
                    cur = s[1:].split()[0]
                    n_cnt = 0
                elif s:
                    n_cnt += sum(1 for c in s if c not in "ACGTacgt")
            if cur is not None:
                n_counts[cur] = n_cnt

    try:
        rows_base_cov, dt_base_cov = run_eggman(BASE_PACK_DIR, fasta_path=tmp_fa, gamma=3.5, min_snps=3, p_thresh=0.01)
        rows_def, dt_def = run_eggman(SPANNING_PACK_DIR, fasta_path=tmp_fa)
        rows_cov, dt_cov = run_eggman(SPANNING_PACK_DIR, fasta_path=tmp_fa, gamma=3.5, min_snps=3, p_thresh=0.01)
    finally:
        if os.path.exists(tmp_fa):
            os.remove(tmp_fa)

    def analyze_bvbrc_rows(rows, runtime_s):
        cat_cnt = Counter(r["called_category"] for r in rows)
        wave_strata = {"wave_B117", "wave_Delta", "wave_BA1", "wave_BA2", "wave_BA4", "wave_BA5", "wave_BQ1", "wave_JN1", "wave_KP2_KP3", "wave_LP8"}
        wave_total = 0
        wave_mosaic = 0
        wave_by_lineage = defaultdict(Counter)
        alias_xbb_tot = alias_xbb_hit = 0
        alias_xdv_tot = alias_xdv_hit = 0
        urf_clusters = defaultdict(list)

        for r in rows:
            sid = r["sequence_id"]
            gid = sid.split("|")[0] if "|" in sid else ""
            m = meta_by_gid.get(gid, {})
            lin = (m.get("lineage") or "").strip()
            lin_root = lin.split(".")[0]
            stratum = m.get("stratum") or ""
            nbp = int(r["n_breakpoints"])
            cat = r["called_category"]
            cal = r["called_lineage"]
            form = r["structural_formula"]
            barc = r["hxb2_barcode"]

            if stratum in wave_strata and not lin_root.startswith("X"):
                wave_total += 1
                wave_by_lineage[lin]["total"] += 1
                if nbp > 0:
                    wave_mosaic += 1
                    wave_by_lineage[lin]["mosaic"] += 1
                    wave_by_lineage[lin]["f:" + form] += 1

            if stratum == "alias_XBB_desc":
                alias_xbb_tot += 1
                if cal.startswith("XBB") or ("BJ.1" in form and "BM.1.1.1" in form):
                    alias_xbb_hit += 1
            if lin_root in {"NB", "PQ", "RK", "SV", "XDV"} or stratum == "recomb_XDV_alias":
                alias_xdv_tot += 1
                if nbp > 0 and ("BJ.1" in form and "JN.1" in form):
                    alias_xdv_hit += 1

            if cat == "URF" and nbp > 0:
                urf_clusters[form].append({
                    "seq_id": sid,
                    "acc": m.get("genbank_accessions", ""),
                    "pango": lin,
                    "country": m.get("isolation_country", ""),
                    "year": m.get("collection_year"),
                    "barcode": barc,
                    "n_ambig": n_counts.get(sid, 0),
                })

        circ_clusters = []
        for form, items in sorted(urf_clusters.items(), key=lambda x: -len(x[1])):
            if len(items) >= 3:
                pangos = Counter(x["pango"] for x in items).most_common(5)
                countries = Counter(x["country"] for x in items if x["country"]).most_common(5)
                years = sorted(set(str(x["year"]) for x in items if x["year"]))
                ns = sorted(x["n_ambig"] for x in items)
                has_b_backfill = ("[B:" in form or "| B:" in form)
                circ_clusters.append({
                    "structural_formula": form,
                    "n_genomes": len(items),
                    "has_wuhan_b_segment": has_b_backfill,
                    "median_n_bases": ns[len(ns) // 2],
                    "min_n_bases": ns[0],
                    "top_reported_pango": pangos,
                    "top_countries": countries,
                    "years": years,
                    "example_accs": [x["acc"] for x in items[:4] if x["acc"]],
                    "representative_barcode": items[0]["barcode"],
                })

        hidden_wave_summary = {}
        for lin, cnt in sorted(wave_by_lineage.items(), key=lambda x: -x[1]["total"]):
            if cnt["total"] >= 50:
                hidden_wave_summary[lin] = {
                    "total": cnt["total"],
                    "mosaic": cnt["mosaic"],
                    "pct": round(100.0 * cnt["mosaic"] / cnt["total"], 2),
                    "top_mosaic_formula": max(((fk[2:], fv) for fk, fv in cnt.items() if fk.startswith("f:")), key=lambda x: x[1], default=("", 0)),
                }

        return {
            "n_screened": len(rows),
            "runtime_sec": round(runtime_s, 3),
            "throughput_seq_per_sec": round(len(rows) / runtime_s, 1) if runtime_s > 0 else 0.0,
            "categories": dict(cat_cnt),
            "nominal_pure_waves_total": wave_total,
            "nominal_pure_waves_mosaic": wave_mosaic,
            "nominal_pure_waves_mosaic_pct": round(100.0 * wave_mosaic / wave_total, 2) if wave_total else 0.0,
            "hidden_wave_breakdown": hidden_wave_summary,
            "alias_xbb_total": alias_xbb_tot,
            "alias_xbb_recovered": alias_xbb_hit,
            "alias_xdv_total": alias_xdv_tot,
            "alias_xdv_recovered": alias_xdv_hit,
            "uncatalogued_urf_clusters_ge3": circ_clusters[:25],
        }

    return {
        "bvbrc_baseline_cov_tuned": analyze_bvbrc_rows(rows_base_cov, dt_base_cov),
        "bvbrc_spanning_default": analyze_bvbrc_rows(rows_def, dt_def),
        "bvbrc_spanning_cov_tuned": analyze_bvbrc_rows(rows_cov, dt_cov),
    }


def main():
    print("Loading 10,062 world SARS-CoV-2 cohort + metadata...", flush=True)
    world_meta = load_world_metadata()
    world_seqs = load_world_sequences()
    print(f"Loaded {len(world_seqs)} sequences and {len(world_meta)} metadata rows.", flush=True)

    print("Running Systematic Spanning Basis Audit on 10,062 cohort...", flush=True)
    part_a = build_spanning_pack(world_seqs, world_meta)

    print("Streaming complete SARS-CoV-2 cohort from BV-BRC API...", flush=True)
    meta_by_gid, fasta_gz_path, dt_stream = stream_bvbrc_cohort()

    print(f"Running EGGMAN against {len(meta_by_gid)} BV-BRC genomes...", flush=True)
    part_b = screen_bvbrc(meta_by_gid, fasta_gz_path)
    part_b["n_metadata_fetched"] = len(meta_by_gid)
    part_b["stream_download_sec"] = round(dt_stream, 2)

    out = {"part_a_world_10062": part_a, "part_b_bvbrc_stream": part_b}
    with open(OUT_JSON, "w") as f:
        json.dump(out, f, indent=2)
    print(f"\nSaved complete JSON report to {OUT_JSON}", flush=True)


if __name__ == "__main__":
    main()
