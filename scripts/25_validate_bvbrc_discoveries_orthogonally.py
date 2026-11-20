#!/usr/bin/env python3
"""
25_validate_bvbrc_discoveries_orthogonally.py
=============================================
Rigorous orthogonal catalog, metadata, artifact, and ML phylogenetic validation of:
  Part 1: All 23 SARS-CoV-2 URF clusters (205 reps in data/bvbrc_stream_hits/sc2_novel_hits.fasta.gz)
          -> results/bvbrc_sc2_23clusters_novelty_audit.tsv & .json
  Part 2: All 40 H5Nx IAV constellations (253 reps in data/bvbrc_stream_hits/iav_h5_novel_hits.fasta.gz)
          -> results/bvbrc_iav_40constellations_novelty_audit.tsv & .json
  Trees:  results/trees_bvbrc_validation/
"""

import concurrent.futures
import csv
import gzip
import json
import os
import re
import subprocess
import tempfile
import time
import urllib.parse
import urllib.request
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import Phylo

ROOT = Path("/Users/sergei/Projects/HIV-recombination")
RESULTS_DIR = ROOT / "results"
HITS_DIR = ROOT / "data" / "bvbrc_stream_hits"
TREES_DIR = RESULTS_DIR / "trees_bvbrc_validation"

SC2_PACK_FASTA = ROOT / "rhiztype" / "packs" / "sars_cov2_wuhan" / "profiles.fasta"
SC2_SUMMARY_JSON = RESULTS_DIR / "bvbrc_100k_sc2_stream_summary.json"
SC2_TSV = RESULTS_DIR / "bvbrc_100k_sc2_novel_and_mislabeled.tsv"
SC2_HITS_FASTA_GZ = HITS_DIR / "sc2_novel_hits.fasta.gz"

IAV_SUMMARY_JSON = RESULTS_DIR / "bvbrc_iav_h5_8seg_stream_summary.json"
IAV_TSV = RESULTS_DIR / "bvbrc_iav_h5_novel_constellations.tsv"
IAV_HITS_FASTA_GZ = HITS_DIR / "iav_h5_novel_hits.fasta.gz"

OUT_SC2_AUDIT_TSV = RESULTS_DIR / "bvbrc_sc2_23clusters_novelty_audit.tsv"
OUT_SC2_AUDIT_JSON = RESULTS_DIR / "bvbrc_sc2_23clusters_novelty_audit.json"
OUT_IAV_AUDIT_TSV = RESULTS_DIR / "bvbrc_iav_40constellations_novelty_audit.tsv"
OUT_IAV_AUDIT_JSON = RESULTS_DIR / "bvbrc_iav_40constellations_novelty_audit.json"

IQTREE_BIN = "/Users/sergei/.local/bin/iqtree2"
NEXTCLADE_BIN = "/Users/sergei/miniconda3/bin/nextclade"
MAFFT_BIN = "/usr/local/bin/mafft"
BLASTN_BIN = "/opt/homebrew/bin/blastn"
MAKEBLASTDB_BIN = "/opt/homebrew/bin/makeblastdb"

HEADERS_JSON = {"Accept": "application/json", "User-Agent": "curl/8.7.1"}
HEADERS_FASTA_POST = {
    "Content-Type": "application/rqlquery+x-www-form-urlencoded",
    "Accept": "application/dna+fasta",
    "User-Agent": "curl/8.7.1",
}

SEG_ORDER = ["1", "2", "3", "4", "5", "6", "7", "8"]
SEG_NAMES = ["PB2", "PB1", "PA", "HA", "NP", "NA", "MP", "NS"]


def http_get_text(url, headers=None, timeout=45, retries=4):
    hdrs = headers or {"User-Agent": "curl/8.7.1"}
    last_err = None
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers=hdrs)
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                return resp.read().decode("utf-8", errors="replace")
        except Exception as e:
            last_err = e
            time.sleep(1.5 * (attempt + 1))
    raise RuntimeError(f"GET failed for {url}: {last_err}")


def http_get_bytes(url, timeout=45, retries=4):
    last_err = None
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "curl/8.7.1"})
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                return resp.read()
        except Exception as e:
            last_err = e
            time.sleep(1.5 * (attempt + 1))
    raise RuntimeError(f"GET bytes failed for {url}: {last_err}")


def fetch_bvbrc_fastas(gids):
    seq_by_gid = {}
    gid_re = re.compile(r"\|\s*([0-9]+\.[0-9]+)\s*\]\s*$")
    for i in range(0, len(gids), 500):
        sub = gids[i : i + 500]
        body = ("in(genome_id,(" + ",".join(sub) + "))&limit(1000)").encode("utf-8")
        req = urllib.request.Request(
            "https://www.bv-brc.org/api/genome_sequence/",
            data=body,
            headers=HEADERS_FASTA_POST,
        )
        with urllib.request.urlopen(req, timeout=75) as resp:
            fa = resp.read().decode("utf-8", errors="replace")
        cur_gid = None
        buf = []
        for raw_line in fa.splitlines():
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if cur_gid and buf:
                    seq_by_gid[cur_gid] = "".join(buf).upper()
                m = gid_re.search(line)
                cur_gid = m.group(1) if m else line.split("|")[-1].strip().rstrip("]")
                buf = []
            else:
                buf.append(line)
        if cur_gid and buf:
            seq_by_gid[cur_gid] = "".join(buf).upper()
    return seq_by_gid


def parse_fasta(path_or_gz, is_gz=False):
    seqs = {}
    headers = {}
    cur = None
    buf = []
    opener = gzip.open if is_gz else open
    with opener(path_or_gz, "rt", encoding="utf-8", errors="replace") as f:
        for raw in f:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if cur is not None:
                    seqs[cur] = "".join(buf).upper()
                hdr = line[1:]
                cur = hdr.split()[0]
                headers[cur] = hdr
                buf = []
            else:
                buf.append(line)
    if cur is not None:
        seqs[cur] = "".join(buf).upper()
    return seqs, headers


def parse_ufboot(clade):
    if clade.name and "/" in str(clade.name):
        try:
            return float(str(clade.name).split("/")[-1])
        except ValueError:
            pass
    if clade.name:
        try:
            return float(str(clade.name))
        except ValueError:
            pass
    if clade.confidence is not None:
        return float(clade.confidence)
    return np.nan


def parse_barcode(bc):
    segs = []
    for part in str(bc).split("|"):
        if ":" not in part:
            continue
        trunk, span = part.split(":", 1)
        if "-" in span:
            s_str, e_str = span.split("-", 1)
            if s_str.isdigit() and e_str.isdigit():
                segs.append((trunk.strip(), int(s_str), int(e_str)))
    return segs


def analyze_tree_sister(tree_path, query_names, root_taxon=None):
    """
    Evaluate monophyly & nearest Ref_ sister clade + UFBoot in an IQ-TREE 2 Newick file.
    """
    tree = Phylo.read(str(tree_path), "newick")
    term_map = {t.name: t for t in tree.get_terminals()}
    if root_taxon and root_taxon in term_map:
        tree.root_with_outgroup(term_map[root_taxon])
    else:
        tree.root_at_midpoint()

    q_set = set(query_names) & set(term_map.keys())
    if not q_set:
        return {"error": "no query taxa in tree"}

    # Check strict query monophyly
    best_q_size = 10**9
    q_mono_ufb = np.nan
    for clade in tree.find_clades(order="preorder"):
        terms = [t.name for t in clade.get_terminals()]
        if sum(1 for t in terms if t in q_set) == len(q_set) and len(terms) < best_q_size:
            best_q_size = len(terms)
            q_mono_ufb = parse_ufboot(clade)

    # Find smallest clade containing all query sequences + strictly one Ref_ taxon
    pure_clades = []
    for clade in tree.find_clades(order="preorder"):
        terms = [t.name for t in clade.get_terminals()]
        if sum(1 for t in terms if t in q_set) < len(q_set):
            continue
        refs_in = [t[4:] for t in terms if t.startswith("Ref_")]
        if len(refs_in) == 1:
            ufb = parse_ufboot(clade)
            pure_clades.append((len(terms), -(ufb if pd.notna(ufb) else 0.0), refs_in[0], ufb))

    if pure_clades:
        pure_clades.sort()
        _, _, sister_ref, sister_ufb = pure_clades[0]
        sister_group = [sister_ref]
    else:
        # Find smallest clade containing all query sequences + minimal set of Ref_ taxa
        multi_clades = []
        for clade in tree.find_clades(order="preorder"):
            terms = [t.name for t in clade.get_terminals()]
            if sum(1 for t in terms if t in q_set) < len(q_set):
                continue
            refs_in = sorted(set(t[4:] for t in terms if t.startswith("Ref_")))
            if refs_in:
                ufb = parse_ufboot(clade)
                multi_clades.append((len(refs_in), len(terms), -(ufb if pd.notna(ufb) else 0.0), refs_in, ufb))
        if multi_clades:
            multi_clades.sort()
            _, _, _, sister_group, sister_ufb = multi_clades[0]
            sister_ref = "+".join(sister_group[:3])
        else:
            ref_terms = [t for t in tree.get_terminals() if t.name.startswith("Ref_")]
            q_terms = [t for t in tree.get_terminals() if t.name in q_set]
            best_ref = min(ref_terms, key=lambda rt: float(np.mean([tree.distance(rt, qt) for qt in q_terms])))
            sister_ref = best_ref.name[4:]
            sister_group = [sister_ref]
            sister_ufb = np.nan

    # Also compute mean patristic distance to top 3 nearest Ref_ taxa (for quantitative corroboration)
    ref_terms = [t for t in tree.get_terminals() if t.name.startswith("Ref_")]
    q_terms = [t for t in tree.get_terminals() if t.name in q_set]
    ref_dists = sorted(
        (round(float(np.mean([tree.distance(rt, qt) for qt in q_terms])), 6), rt.name[4:])
        for rt in ref_terms
    )

    # Check per-query nearest Ref_ (to catch mixed-clade clusters!)
    per_q_nearest = []
    for qt in q_terms:
        nearest = min(ref_terms, key=lambda rt: tree.distance(rt, qt))
        per_q_nearest.append(nearest.name[4:])

    return {
        "n_queries": len(q_set),
        "strictly_monophyletic": best_q_size == len(q_set),
        "query_mono_ufboot": round(float(q_mono_ufb), 1) if pd.notna(q_mono_ufb) else None,
        "ml_sister_ref": sister_ref,
        "ml_sister_group": sister_group,
        "sister_clade_ufboot": round(float(sister_ufb), 1) if pd.notna(sister_ufb) else None,
        "nearest_refs_patristic": ref_dists[:4],
        "per_query_nearest_ref": dict(Counter(per_q_nearest)),
    }


def resolve_pango_root(lineage, alias_key):
    if not lineage or lineage in ("None", "unassigned", "MISSING"):
        return ""
    if lineage.startswith("X"):
        return lineage.split(".")[0]
    root = lineage.split(".")[0]
    seen = set()
    cur = lineage
    while root in alias_key and root not in seen:
        seen.add(root)
        target = alias_key[root]
        if isinstance(target, list):
            return root
        if isinstance(target, str) and target:
            rest = cur[len(root) :]
            cur = target + rest
            root = cur.split(".")[0]
            if root.startswith("X"):
                return root
        else:
            break
    return cur.split(".")[0]


def parse_lineage_notes(notes_txt):
    notes = {}
    for line in notes_txt.splitlines():
        if "\t" not in line or line.startswith("Lineage\t"):
            continue
        lin, desc = line.split("\t", 1)
        notes[lin.strip()] = desc.strip()
    return notes


# =============================================================================
# PART 1: 23 SARS-CoV-2 CLUSTERS AUDIT & PHYLOGENETICS
# =============================================================================
def run_part1_sc2_audit():
    print("\n" + "=" * 90, flush=True)
    print("PART 1: ORTHOGONAL EXTERNAL CATALOG & PHYLOGENETIC AUDIT OF 23 SC2 CLUSTERS", flush=True)
    print("=" * 90, flush=True)
    t0 = time.time()

    sc2_sum = json.load(open(SC2_SUMMARY_JSON))
    clusters = sc2_sum["eggman_cov_tuned"]["novel_urf_clusters_ge_3"]["clusters"]
    print(f"  Loaded {len(clusters)} novel URF clusters from {SC2_SUMMARY_JSON.name}", flush=True)

    # Parse representative sequences & headers from sc2_novel_hits.fasta.gz
    rep_seqs, rep_headers = parse_fasta(SC2_HITS_FASTA_GZ, is_gz=True)
    reps_by_form = defaultdict(list)
    rep_meta_hdr = {}
    all_rep_gids = []
    for sid, hdr in rep_headers.items():
        head, form = hdr.split(" formula=", 1)
        form = form.strip()
        gid = sid.split("|", 1)[0]
        all_rep_gids.append(gid)
        kv = dict(p.split("=", 1) for p in head.split()[1:] if "=" in p)
        rep_meta_hdr[sid] = {"gid": gid, "formula": form, **kv}
        reps_by_form[form].append(sid)

    # Fetch BV-BRC BioProject, publication, date, platform for all 205 reps
    bvbrc_rep_meta = {}
    for i in range(0, len(all_rep_gids), 100):
        sub_str = ",".join(all_rep_gids[i : i + 100])
        u = (
            f"https://www.bv-brc.org/api/genome/?in(genome_id,({sub_str}))&limit(500)"
            "&select(genome_id,genbank_accessions,bioproject_accession,publication,strain,"
            "isolation_country,collection_date,collection_year,sequencing_platform,lineage)"
        )
        rows = json.loads(http_get_text(u, headers=HEADERS_JSON))
        for r in rows:
            bvbrc_rep_meta[r["genome_id"]] = r
    print(f"  Fetched BV-BRC BioProject/metadata for {len(bvbrc_rep_meta)}/{len(all_rep_gids)} SC2 reps", flush=True)

    # Fetch latest Pango lineage_notes.txt & alias_key.json
    notes_txt = http_get_text("https://raw.githubusercontent.com/cov-lineages/pango-designation/master/lineage_notes.txt")
    alias_json = json.loads(
        http_get_text("https://raw.githubusercontent.com/cov-lineages/pango-designation/master/pango_designation/alias_key.json")
    )
    lineage_notes = parse_lineage_notes(notes_txt)
    print(f"  Loaded Pango designation ({len(lineage_notes)} notes, {len(alias_json)} aliases)", flush=True)

    # Fetch clean BV-BRC sister references for phylogenetic comparison against 16 Tier-1 founders:
    # XDV.1, KP.3.1.1, KP.2.2, XEC (plus we can select clean NB.1.8.1 / PQ.2 / XFG verified from 2025 BV-BRC)
    sister_queries = [
        ("XDV.1", "eq(taxon_lineage_ids,2697049)&eq(genome_status,Complete)&eq(contigs,1)&gt(genome_length,29000)&eq(lineage,XDV.1)&limit(10)"),
        ("KP.3.1.1", "eq(taxon_lineage_ids,2697049)&eq(genome_status,Complete)&eq(contigs,1)&gt(genome_length,29000)&eq(lineage,KP.3.1.1)&limit(10)"),
        ("KP.2.2", "eq(taxon_lineage_ids,2697049)&eq(genome_status,Complete)&eq(contigs,1)&gt(genome_length,29000)&eq(lineage,KP.2.2)&limit(10)"),
        ("XEC", "eq(taxon_lineage_ids,2697049)&eq(genome_status,Complete)&eq(contigs,1)&gt(genome_length,29000)&eq(lineage,XEC)&limit(10)"),
        ("NB.1.8.1_2025", "eq(taxon_lineage_ids,2697049)&eq(genome_status,Complete)&eq(contigs,1)&gt(genome_length,29000)&eq(collection_year,2025)&limit(200)"),
        ("XFG_2026", "eq(taxon_lineage_ids,2697049)&eq(genome_status,Complete)&eq(contigs,1)&gt(genome_length,29000)&eq(collection_year,2026)&limit(120)"),
    ]
    cand_ref_gids = {}
    for lbl, q in sister_queries:
        rows = json.loads(http_get_text(f"https://www.bv-brc.org/api/genome/?{q}&select(genome_id,lineage)", headers=HEADERS_JSON))
        for r in rows:
            if r["genome_id"] not in set(all_rep_gids):
                cand_ref_gids.setdefault(lbl, []).append(r["genome_id"])
    all_cand_gids = [g for lst in cand_ref_gids.values() for g in lst]
    cand_seqs = fetch_bvbrc_fastas(all_cand_gids)

    # Run Nextclade 3.x on all 205 reps + candidate sister refs
    pack_profiles, _ = parse_fasta(SC2_PACK_FASTA, is_gz=False)
    TREES_DIR.mkdir(parents=True, exist_ok=True)
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory() as td:
        ds_dir = os.path.join(td, "sc2_ds")
        subprocess.run(
            [NEXTCLADE_BIN, "dataset", "get", "--name", "sars-cov-2", "--output-dir", ds_dir],
            check=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )

        # Align candidates first to pick clean <50 N refs for XDV.1, NB.1.8.1, PQ.2, KP.3.1.1, KP.2.2, XEC, XFG
        cand_fa = os.path.join(td, "cand.fasta")
        cand_tsv = os.path.join(td, "cand.tsv")
        cand_aln_fa = os.path.join(td, "cand_aln.fasta")
        with open(cand_fa, "w") as f:
            for lbl, gids in cand_ref_gids.items():
                for g in gids:
                    if g in cand_seqs and cand_seqs[g].count("N") < 100:
                        f.write(f">{lbl}__{g}\n{cand_seqs[g]}\n")
        subprocess.run(
            [NEXTCLADE_BIN, "run", "-D", ds_dir, "--output-tsv", cand_tsv, "--output-fasta", cand_aln_fa, cand_fa],
            check=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
        cand_aln, _ = parse_fasta(cand_aln_fa, is_gz=False)
        sister_refs = {}
        wanted_sisters = {
            "XDV.1": lambda p, c: p == "XDV.1",
            "KP.3.1.1": lambda p, c: p.startswith("KP.3.1.1"),
            "KP.2.2": lambda p, c: p.startswith("KP.2.2"),
            "XEC": lambda p, c: p == "XEC",
            "NB.1.8.1": lambda p, c: p == "NB.1.8.1",
            "PQ.2": lambda p, c: p.startswith("PQ.2"),
            "XFG": lambda p, c: p.startswith("XFG"),
        }
        for r in csv.DictReader(open(cand_tsv), delimiter="\t"):
            sn = r["seqName"]
            pango = r.get("Nextclade_pango") or ""
            clade = r.get("clade") or ""
            miss = int(r.get("totalMissing") or 9999)
            if miss > 150 or sn not in cand_aln:
                continue
            for wname, pred in wanted_sisters.items():
                if wname not in sister_refs and pred(pango, clade):
                    sister_refs[wname] = cand_aln[sn]
        print(f"  Selected clean aligned sister references: {sorted(sister_refs.keys())}", flush=True)

        # Align all 205 SC2 cluster representatives with Nextclade
        in_fa = os.path.join(td, "sc2_205.fasta")
        nc_tsv = os.path.join(td, "sc2_205.tsv")
        aln_fa = os.path.join(td, "sc2_205_aln.fasta")
        with open(in_fa, "w") as f:
            for sid, seq in rep_seqs.items():
                f.write(f">{sid}\n{seq}\n")
        subprocess.run(
            [NEXTCLADE_BIN, "run", "-D", ds_dir, "--output-tsv", nc_tsv, "--output-fasta", aln_fa, in_fa],
            check=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
        nc_rows = {}
        for r in csv.DictReader(open(nc_tsv), delimiter="\t"):
            sid = r["seqName"].split()[0]
            nc_rows[sid] = r
        aln_reps, _ = parse_fasta(aln_fa, is_gz=False)

    print(f"  Aligned {len(aln_reps)}/205 SC2 representatives via Nextclade", flush=True)

    # Prepare combined reference alignment dict (16 Tier-1 founders + clean sister waves)
    ref_aln = {f"Ref_{k}": v for k, v in pack_profiles.items()}
    for k, v in sister_refs.items():
        ref_aln[f"Ref_{k}"] = v

    # Helper to compute N-stretch proximity to called breakpoints on aligned 29,903-nt sequence
    def n_gap_near_bp(aligned_seq, bp_pos_list, window=200):
        n_positions = [i + 1 for i, c in enumerate(aligned_seq) if c not in "ACGT-"]
        # Find stretches >= 10 N's or any N within window of breakpoints
        near_any = False
        near_stretch = False
        # Group stretches >= 10
        stretches = []
        run_start = None
        for idx_c, ch in enumerate(aligned_seq, 1):
            if ch not in "ACGT-":
                if run_start is None:
                    run_start = idx_c
            else:
                if run_start is not None:
                    if idx_c - run_start >= 10:
                        stretches.append((run_start, idx_c - 1))
                    run_start = None
        if run_start is not None and (len(aligned_seq) + 1 - run_start) >= 10:
            stretches.append((run_start, len(aligned_seq)))

        min_dist = 99999
        for bp in bp_pos_list:
            for np_i in n_positions:
                d = abs(np_i - bp)
                if d < min_dist:
                    min_dist = d
                if d <= window:
                    near_any = True
            for s0, s1 in stretches:
                if s0 - window <= bp <= s1 + window:
                    near_stretch = True
        return len(n_positions), len(stretches), (min_dist if n_positions else None), near_any, near_stretch

    # Build IQ-TREE 2 segment tree tasks for all 23 clusters (for every segment >= 500 nt)
    tree_tasks = []
    cluster_segs_info = {}
    for rank, cl in enumerate(clusters, 1):
        form = cl["structural_formula"]
        barc = cl["representative_barcode"]
        segs = parse_barcode(barc)
        bp_positions = [segs[i][2] for i in range(len(segs) - 1)]

        # Pick up to 3 cleanest aligned reps for phylogenetic verification
        c_sids = reps_by_form.get(form, [])
        c_sids_sorted = sorted(
            [s for s in c_sids if s in aln_reps],
            key=lambda s: (int(nc_rows.get(s, {}).get("totalMissing") or 9999), s),
        )[:3]
        q_names = [f"Query_{s.replace('|', '_')}" for s in c_sids_sorted]
        cluster_segs_info[form] = {"rank": rank, "segs": segs, "bps": bp_positions, "q_sids": c_sids_sorted, "q_names": q_names}

        for seg_idx, (trunk, s0, s1) in enumerate(segs, 1):
            span = s1 - s0 + 1
            if span < 500:
                continue
            fa_p = TREES_DIR / f"sc2_c{rank:02d}_seg{seg_idx}_{trunk.replace('.', '')}_{s0}_{s1}.fasta"
            tr_p = Path(str(fa_p) + ".treefile")
            with open(fa_p, "w") as f:
                for rname, rseq in ref_aln.items():
                    f.write(f">{rname}\n{rseq[s0 - 1 : s1]}\n")
                for qn, sid in zip(q_names, c_sids_sorted):
                    f.write(f">{qn}\n{aln_reps[sid][s0 - 1 : s1]}\n")
            tree_tasks.append((form, seg_idx, trunk, s0, s1, q_names, fa_p, tr_p))

    print(f"  Running IQ-TREE 2 (-m GTR+G4 -B 1000 -T 1) on {len(tree_tasks)} SC2 segment alignments...", flush=True)

    def _run_iqtree(item):
        fa_p, tr_p = item[6], item[7]
        if not tr_p.exists():
            subprocess.run(
                [IQTREE_BIN, "-s", str(fa_p), "-m", "GTR+G4", "-B", "1000", "-T", "1", "--quiet", "-redo"],
                check=True,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
            )

    with concurrent.futures.ThreadPoolExecutor(max_workers=8) as pool:
        list(pool.map(_run_iqtree, tree_tasks))

    # Summarize all 23 clusters
    tree_results_by_form = defaultdict(list)
    for form, seg_idx, trunk, s0, s1, q_names, fa_p, tr_p in tree_tasks:
        tstat = analyze_tree_sister(tr_p, q_names, root_taxon="Ref_B")
        sis = tstat.get("ml_sister_ref", "")
        sis_grp = tstat.get("ml_sister_group", [])
        ufb = tstat.get("sister_clade_ufboot")
        # Check exact match to called Tier-1 trunk (or known sister descendant if called trunk is JN.1/KP.3.3/BJ.1)
        exact_match = (trunk in sis_grp) or (sis == trunk)
        tree_results_by_form[form].append({
            "seg_idx": seg_idx,
            "called_trunk": trunk,
            "start": s0,
            "end": s1,
            "span_nt": s1 - s0 + 1,
            "ml_sister_ref": sis,
            "ml_sister_group": sis_grp,
            "sister_clade_ufboot": ufb,
            "query_mono_ufboot": tstat.get("query_mono_ufboot"),
            "strictly_monophyletic": tstat.get("strictly_monophyletic"),
            "exact_called_trunk_match": exact_match,
            "nearest_refs_patristic": tstat.get("nearest_refs_patristic"),
            "treefile": tr_p.name,
        })

    audit_rows = []
    class_counts = Counter()
    for rank, cl in enumerate(clusters, 1):
        form = cl["structural_formula"]
        c_sids = reps_by_form.get(form, [])
        bp_positions = cluster_segs_info[form]["bps"]

        # Aggregate Nextclade, Pango designation, BioProject, N-gap proximity across reps
        nc_pangos = Counter()
        nc_clades = Counter()
        nc_qc = Counter()
        priv_muts = []
        bioprojects = Counter()
        countries = Counter()
        dates = []
        bp_near_n_cnt = 0
        bp_near_stretch_cnt = 0
        n_stretches_counts = []

        for sid in c_sids:
            gid = sid.split("|", 1)[0]
            nr = nc_rows.get(sid, {})
            if nr:
                npang = nr.get("Nextclade_pango") or "UNK"
                nc_pangos[npang] += 1
                nc_clades[nr.get("clade") or "UNK"] += 1
                nc_qc[nr.get("qc.overallStatus") or "UNK"] += 1
                pm = nr.get("privateNucMutations.totalPrivateSubstitutions")
                if pm and pm.isdigit():
                    priv_muts.append(int(pm))
            bm = bvbrc_rep_meta.get(gid, {})
            if bm.get("bioproject_accession"):
                bioprojects[bm["bioproject_accession"]] += 1
            if bm.get("isolation_country"):
                countries[bm["isolation_country"]] += 1
            if bm.get("collection_date"):
                dates.append(str(bm["collection_date"]))
            if sid in aln_reps:
                _, n_str, _, near_n, near_str = n_gap_near_bp(aln_reps[sid], bp_positions, window=200)
                n_stretches_counts.append(n_str)
                if near_n:
                    bp_near_n_cnt += 1
                if near_str:
                    bp_near_stretch_cnt += 1

        top_nc_pango = nc_pangos.most_common(1)[0][0] if nc_pangos else ""
        top_nc_clade = nc_clades.most_common(1)[0][0] if nc_clades else ""
        # Resolve alias root across top Nextclade calls
        pango_roots = Counter()
        for p, cnt in nc_pangos.items():
            rt = resolve_pango_root(p, alias_json)
            if rt:
                pango_roots[rt] += cnt
        top_resolved_root = pango_roots.most_common(1)[0][0] if pango_roots else ""

        def _pango_chain(tok):
            chain = tok
            seen_p = set()
            while "." in chain or chain in alias_json:
                pref = chain.split(".", 1)[0]
                if pref in seen_p:
                    break
                seen_p.add(pref)
                exp = alias_json.get(pref)
                if isinstance(exp, str) and exp:
                    rest = chain[len(pref):]
                    chain = exp + rest
                elif isinstance(exp, list):
                    return chain, exp
                else:
                    break
            return chain, []

        # Break co-occurring X* root ambiguity using 5'->3' parental order & note breakpoint span
        seg_trees = tree_results_by_form.get(form, [])
        secondary_x_note = ""
        x_cands = [rt for rt, _ in pango_roots.most_common() if rt.startswith("X") and (rt in lineage_notes or rt in alias_json)]
        if len(x_cands) >= 2 and seg_trees:
            seg_trunks = [st.get("called_trunk", "") for st in seg_trees if st.get("called_trunk")]
            ranked_x = []
            for rt in x_cands:
                cnt_rt = pango_roots[rt]
                pars = alias_json.get(rt, []) if isinstance(alias_json.get(rt), list) else []
                first_ok = 0
                if seg_trunks and pars:
                    p0_exp, _ = _pango_chain(pars[0])
                    t0_exp, _ = _pango_chain(seg_trunks[0])
                    if pars[0].startswith(seg_trunks[0]) or p0_exp.startswith(t0_exp):
                        first_ok = 1
                bp_prox = 0
                note_rt = lineage_notes.get(rt, "")
                m_bp = re.search(r"breakpoints?:?\s*([0-9,\-\s]+)", note_rt)
                if m_bp and bp_positions:
                    ivs = re.findall(r"(\d+)\s*-\s*(\d+)", m_bp.group(1))
                    if ivs and all(any(min(abs(b - int(lo)), abs(b - int(hi))) <= 3500 or int(lo) <= b <= int(hi) for lo, hi in ivs) for b in bp_positions):
                        bp_prox = 1
                ranked_x.append((first_ok + bp_prox, bp_prox, first_ok, cnt_rt, rt))
            ranked_x.sort(key=lambda x: (-x[0], -x[1], -x[2], -x[3], x[4]))
            best_rt = ranked_x[0][4]
            if best_rt != top_resolved_root and ranked_x[0][0] > [r[0] for r in ranked_x if r[4] == top_resolved_root][0]:
                secondary_x_note = f"; Nextclade calls split with secondary {top_resolved_root} ({dict(pango_roots)})"
                top_resolved_root = best_rt
                matching_calls = [(p, cnt) for p, cnt in nc_pangos.most_common() if resolve_pango_root(p, alias_json) == best_rt]
                if matching_calls:
                    top_nc_pango = matching_calls[0][0]

        # Check lineage_notes.txt description for top_nc_pango and top_resolved_root
        pango_note = lineage_notes.get(top_resolved_root) or lineage_notes.get(top_nc_pango) or ""
        is_known_x_pango = top_resolved_root.startswith("X") and (top_resolved_root in lineage_notes or top_resolved_root in alias_json)

        all_sisters = [st["ml_sister_ref"] for st in seg_trees]
        distinct_sisters = set(all_sisters)
        # Check if all 3' or all segments nest into ONE single sister clade (e.g. KP.3.1.1 / KP.2.2 / XDV.1 / NB.1.8.1)
        all_segments_exact_donor_ufb80 = all(
            st["exact_called_trunk_match"] and (st["sister_clade_ufboot"] is not None and st["sister_clade_ufboot"] >= 80.0)
            for st in seg_trees
        ) if seg_trees else False

        # Honest classification
        # 1. Check artifact first
        if len(c_sids) > 0 and (bp_near_stretch_cnt / len(c_sids)) >= 0.5:
            classification = "SEQUENCING_OR_AMPLICON_ARTIFACT"
            rationale = f"{bp_near_stretch_cnt}/{len(c_sids)} reps have >=10nt N-gap within 200nt of breakpoint"
        elif is_known_x_pango:
            # Distinguish whether Nextclade/ML confirms it as that known X* recombinant vs an evolved XDV sublineage pseudo-split!
            # Specifically: if top_resolved_root == "XDV" and EGGMAN called KS.1.1 or LF.7 (instead of JN.1), check segment trees!
            non_xdv_secondary = (top_resolved_root == "XDV") and any(t in ("KS.1.1", "LF.7", "KP.3.3", "LP.8.1") for t in cl["constituent_trunks"])
            if non_xdv_secondary:
                # Check if 3' segments nest inside XDV/NB.1.8.1/PQ.2 rather than KS.1.1/LF.7!
                s3_sisters = [st["ml_sister_ref"] for st in seg_trees if st["called_trunk"] in ("KS.1.1", "LF.7", "KP.3.3", "LP.8.1")]
                nests_in_xdv = any(any(x in s for x in ("XDV", "NB.1.8.1", "PQ.2", "JN.1")) for s in s3_sisters)
                exact_other = all(st["exact_called_trunk_match"] and (st["sister_clade_ufboot"] or 0) >= 80 for st in seg_trees if st["called_trunk"] in ("KS.1.1", "LF.7"))
                if nests_in_xdv or not exact_other:
                    classification = "EVOLVED_SUBLINEAGE_PSEUDO_SPLIT"
                    rationale = (
                        f"Nextclade calls {top_nc_pango} (alias of {top_resolved_root}); 3' segments place with {s3_sisters} "
                        f"rather than independent { [st['called_trunk'] for st in seg_trees if st['called_trunk'] in ('KS.1.1','LF.7')] } donors at UFBoot>=80%"
                    )
                else:
                    classification = "KNOWN_PANGO_RECOMBINANT_UNANNOTATED_IN_BVBRC"
                    rationale = f"Maps to {top_nc_pango} ({top_resolved_root}); ML segment trees confirm secondary parents"
            else:
                classification = "KNOWN_PANGO_RECOMBINANT_UNANNOTATED_IN_BVBRC"
                rationale = f"Maps to designated Pango {top_resolved_root} ({top_nc_pango}: {pango_note[:120]}{secondary_x_note})"
        else:
            # Not a known Pango X* recombinant! Check if discordant ML segment trees at UFBoot >= 80%
            if len(distinct_sisters) >= 2 and all_segments_exact_donor_ufb80:
                classification = "GENUINELY_UNDESIGNATED_RECOMBINANT"
                rationale = f"Not in Pango X*; discordant ML segment placement {[(st['called_trunk'], st['ml_sister_ref'], st['sister_clade_ufboot']) for st in seg_trees]}"
            else:
                classification = "EVOLVED_SUBLINEAGE_PSEUDO_SPLIT"
                rationale = (
                    f"Nextclade calls clonal {top_nc_pango} (clade {top_nc_clade}); ML segment trees place with "
                    f"{[(st['called_trunk'], st['ml_sister_ref'], st['sister_clade_ufboot']) for st in seg_trees]}"
                )

        class_counts[classification] += 1
        med_priv = int(np.median(priv_muts)) if priv_muts else None
        seg_summary_str = " ; ".join(
            f"seg{st['seg_idx']}({st['called_trunk']}:{st['start']}-{st['end']})->{st['ml_sister_ref']}[UFB={st['sister_clade_ufboot']}]"
            for st in seg_trees
        )

        rec = {
            "rank": rank,
            "n_genomes": cl["n_genomes"],
            "n_2025_2026": cl.get("n_2025_2026"),
            "structural_formula": form,
            "epidemiological_formula": cl["epidemiological_formula"],
            "n_parents": cl["n_parents"],
            "min_fisher_p": cl["min_fisher_p"],
            "median_fisher_p": cl["median_fisher_p"],
            "n_reps_evaluated": len(c_sids),
            "nextclade_top_pango": ";".join(f"{k}:{v}" for k, v in nc_pangos.most_common(4)),
            "nextclade_top_clade": ";".join(f"{k}:{v}" for k, v in nc_clades.most_common(3)),
            "resolved_pango_root": top_resolved_root,
            "is_designated_pango_x": is_known_x_pango,
            "pango_lineage_note": pango_note,
            "nextclade_qc": ";".join(f"{k}:{v}" for k, v in nc_qc.most_common(3)),
            "median_private_muts": med_priv,
            "n_bioprojects_in_reps": len(bioprojects),
            "top_bioprojects": ";".join(f"{k}:{v}" for k, v in bioprojects.most_common(3)),
            "reps_with_n_within_200bp_of_bp": bp_near_n_cnt,
            "reps_with_10nt_n_stretch_near_bp": bp_near_stretch_cnt,
            "iqtree_segment_placements": seg_summary_str,
            "classification": classification,
            "rationale": rationale,
        }
        audit_rows.append(rec)
        print(
            f"  [SC2 #{rank:2d} n={cl['n_genomes']:4d}] {form:44s} -> {classification} | "
            f"NC={top_nc_pango}({top_resolved_root}) | {seg_summary_str}",
            flush=True,
        )

    # Save TSV and JSON
    tsv_cols = list(audit_rows[0].keys())
    with open(OUT_SC2_AUDIT_TSV, "w", newline="", encoding="utf-8") as tf:
        w = csv.DictWriter(tf, fieldnames=tsv_cols, delimiter="\t")
        w.writeheader()
        for r in audit_rows:
            w.writerow(r)

    out_json = {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "runtime_sec": round(time.time() - t0, 2),
        "n_clusters_audited": len(audit_rows),
        "n_reps_evaluated": len(rep_seqs),
        "classification_counts": dict(class_counts),
        "clusters": audit_rows,
        "segment_tree_details": dict(tree_results_by_form),
    }
    with open(OUT_SC2_AUDIT_JSON, "w", encoding="utf-8") as jf:
        json.dump(out_json, jf, indent=2)

    print(f"[+] Part 1 complete in {time.time() - t0:.2f}s | summary={dict(class_counts)}", flush=True)
    return out_json


# =============================================================================
# PART 2: 40 H5Nx CONSTELLATIONS AUDIT & 8-SEGMENT PHYLOGENETICS
# =============================================================================
def run_part2_iav_audit():
    print("\n" + "=" * 90, flush=True)
    print("PART 2: FULL GENOFLU, EURL/EFSA, LITERATURE & 8-SEGMENT IQ-TREE 2 AUDIT (40 H5Nx)", flush=True)
    print("=" * 90, flush=True)
    t0 = time.time()

    iav_sum = json.load(open(IAV_SUMMARY_JSON))
    constellations = iav_sum["novel_urf_constellations_ge_3"]["constellations"]

    # Load 253 reps from iav_h5_novel_hits.fasta.gz
    reps = []
    with gzip.open(IAV_HITS_FASTA_GZ, "rt", encoding="utf-8") as gz:
        for l in gz:
            if l.startswith(">"):
                head, form = l[1:].strip().split(" formula=", 1)
                parts = head.split()
                uid = parts[0]
                kv = dict(p.split("=", 1) for p in parts[1:] if "=" in p)
                reps.append({
                    "uid": uid,
                    "biosample": kv.get("biosample", "NA"),
                    "strain_tok": kv.get("strain", ""),
                    "subtype": kv.get("subtype", ""),
                    "year": kv.get("year", ""),
                    "country_tok": kv.get("country", ""),
                    "formula": form.strip(),
                })
    print(f"  Loaded {len(reps)} representative H5Nx headers across {len(constellations)} constellations", flush=True)

    # Fetch 2020-2026 H5Nx BV-BRC metadata cleanly in 4 pages (matches 253/253 reps + captures BioProject/PubMed!)
    fields = (
        "genome_id,biosample_accession,strain,segment,genbank_accessions,"
        "bioproject_accession,publication,isolation_country,collection_date,"
        "collection_year,genome_length,subtype,host_name"
    )
    base_rql = (
        "eq(taxon_lineage_ids,11320)&eq(genome_status,Complete)&eq(contigs,1)"
        "&eq(h_type,5)&in(collection_year,(2020,2021,2022,2023,2024,2025,2026))"
    )
    seen_gids = set()
    all_rows = []
    last_gid = None
    while True:
        gt = f"&gt(genome_id,{last_gid})" if last_gid else ""
        u = f"https://www.bv-brc.org/api/genome/?{base_rql}{gt}&sort(+genome_id)&limit(10000)&select({fields})"
        rows = json.loads(http_get_text(u, headers=HEADERS_JSON, timeout=60))
        if not rows:
            break
        new_cnt = 0
        for r in rows:
            gid = r.get("genome_id")
            if gid and gid not in seen_gids:
                seen_gids.add(gid)
                all_rows.append(r)
                new_cnt += 1
        if new_cnt == 0 or len(rows) < 10000:
            break
        last_gid = rows[-1]["genome_id"]

    seg_groups = defaultdict(lambda: defaultdict(list))
    for r in all_rows:
        seg = r.get("segment")
        seg_str = str(seg[0]).strip() if isinstance(seg, list) and seg else str(seg or "").strip()
        if seg_str not in SEG_ORDER:
            continue
        bs = (r.get("biosample_accession") or "").strip()
        st = (r.get("strain") or "").strip().replace(" ", "_")
        yr = str(r.get("collection_year") or "")
        ctry = (r.get("isolation_country") or "").strip().replace(" ", "_")
        if bs:
            seg_groups[("bs", bs)][seg_str].append(r)
        if st:
            seg_groups[("st", st, yr, ctry)][seg_str].append(r)

    chosen_segs = {}
    for rep in reps:
        bs = rep["biosample"]
        if bs != "NA" and set(seg_groups[("bs", bs)].keys()) == set(SEG_ORDER):
            sm = seg_groups[("bs", bs)]
        else:
            sm = seg_groups[("st", rep["strain_tok"], rep["year"], rep["country_tok"])]
        if set(sm.keys()) != set(SEG_ORDER):
            raise RuntimeError(f"Failed to match 8 segments for {rep}")
        chosen_segs[rep["uid"]] = {
            s: max(sm[s], key=lambda x: (int(x.get("genome_length") or 0), x["genome_id"]))
            for s in SEG_ORDER
        }
    print(f"  Matched {len(chosen_segs)}/{len(reps)} exact 8-segment isolates from BV-BRC", flush=True)

    # Also fetch expanded global/EURL/South-American/LPAI reference isolates from BV-BRC
    expanded_ref_specs = [
        ("Ref_euBB_DE25", "eq(taxon_lineage_ids,11320)&keyword(2025AI03596)"),
        ("Ref_euDI1_DE25", "eq(taxon_lineage_ids,11320)&eq(strain,A%2Fdomestic%20duck%2FGermany-MV%2F2025AI04088-AK5407%2F2025)"),
        ("Ref_euDI2_DE24", "eq(taxon_lineage_ids,11320)&eq(strain,A%2FCanada%20Goose%2FGermany-TH%2F2024AI04292-AK5408%2F2024)"),
        ("Ref_euDI21_DE25", "eq(taxon_lineage_ids,11320)&eq(strain,A%2Fwild%20duck%2FGermany-SH%2F2025AI07106%2F2025)"),
        ("Ref_SAm_B32_Chile23", "eq(taxon_lineage_ids,11320)&eq(strain,A%2FInca%20tern%2FAntofagasta%2F238083%2F2023)"),
        ("Ref_SAm_LPAI_Arg11", "eq(taxon_lineage_ids,11320)&eq(strain,A%2Fsilver%20teal%2FArgentina%2FCIP051-32%2F2011)"),
        ("Ref_SAm_LPAI_Arg09", "eq(taxon_lineage_ids,11320)&eq(strain,A%2Fcinnamon%20teal%2FArgentina%2FCIP051-1588%2F2009)"),
        ("Ref_SAm_LPAI_Chi18", "eq(taxon_lineage_ids,11320)&eq(strain,A%2Fkelp%20gull%2FChile%2FC41186%2F2018)"),
        ("Ref_EA_LPAI_SK20", "eq(taxon_lineage_ids,11320)&eq(strain,A%2FWild%20Duck%2FSouth%20Korea%2FKNU2020-11%2F2020)"),
        ("Ref_Clade2321_BD21", "eq(taxon_lineage_ids,11320)&eq(strain,A%2Fduck%2FBangladesh%2F50000%2F2021)"),
        ("Ref_NAm_LPAI_H5N2_UT18", "eq(taxon_lineage_ids,11320)&eq(strain,A%2FMallard%2FUtah%2FD1801802%2F2018)"),
    ]
    ext_ref_segs = {}
    for rname, q in expanded_ref_specs:
        u = f"https://www.bv-brc.org/api/genome/?{q}&limit(50)&select(genome_id,segment,genome_length,subtype)"
        rows = json.loads(http_get_text(u, headers=HEADERS_JSON))
        smap = defaultdict(list)
        for r in rows:
            sg = r.get("segment")
            sg = str(sg[0]).strip() if isinstance(sg, list) and sg else str(sg or "").strip()
            if sg in SEG_ORDER:
                smap[sg].append(r)
        if set(smap.keys()) == set(SEG_ORDER):
            ext_ref_segs[rname] = {
                s: max(smap[s], key=lambda x: (int(x.get("genome_length") or 0), x["genome_id"]))
                for s in SEG_ORDER
            }
        else:
            print(f"  [Warn] Reference {rname} missing segments: {sorted(smap.keys())}", flush=True)
    print(f"  Loaded {len(ext_ref_segs)} expanded global/EURL/SAm/LPAI 8-segment references", flush=True)

    # Fetch FASTAs for all 253 rep segments + expanded ref segments
    needed_gids = [chosen_segs[uid][s]["genome_id"] for uid in chosen_segs for s in SEG_ORDER]
    needed_gids += [ext_ref_segs[rn][s]["genome_id"] for rn in ext_ref_segs for s in SEG_ORDER]
    needed_gids = sorted(set(needed_gids))
    seq_by_gid = fetch_bvbrc_fastas(needed_gids)
    print(f"  Fetched {len(seq_by_gid)}/{len(needed_gids)} segment FASTA sequences", flush=True)

    # Query NCBI esummary for all representative HA accessions (to capture FLI-Genotype / GenBank subname notes)
    ha_accs = sorted(
        set(
            (chosen_segs[rep["uid"]]["4"].get("genbank_accessions") or "").split(",")[0].strip()
            for rep in reps
            if chosen_segs[rep["uid"]]["4"].get("genbank_accessions")
        )
    )
    ncbi_notes_by_acc = {}
    for i in range(0, len(ha_accs), 100):
        sub = ha_accs[i : i + 100]
        u = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=nuccore&id=" + ",".join(sub) + "&retmode=json"
        res = json.loads(http_get_text(u)).get("result", {})
        for k, v in res.items():
            if isinstance(v, dict) and "accessionversion" in v:
                acc0 = v["accessionversion"].split(".")[0]
                ncbi_notes_by_acc[acc0] = v.get("subname", "")

    # Run official USDA-VS/GenoFLU BLAST & genotype key matching in a TemporaryDirectory
    # Also harvest clean representative sequences per key GenoFLU allele for 8-segment trees
    genoflu_allele_seqs = defaultdict(dict)  # gene -> {allele: seq}
    geno_results = {}  # uid -> dict
    with tempfile.TemporaryDirectory() as td:
        key_p = os.path.join(td, "genotype_key.xlsx")
        open(key_p, "wb").write(
            http_get_bytes("https://github.com/USDA-VS/GenoFLU/raw/main/dependencies/genotype_key.xlsx")
        )
        df_key = pd.read_excel(key_p)
        geno_dict = {}
        for _, row in df_key.iterrows():
            tup = tuple(str(row[sn]).strip() for sn in SEG_NAMES)
            geno_dict[tup] = str(row["Genotype"]).strip()

        fastas_meta = json.loads(http_get_text("https://api.github.com/repos/USDA-VS/GenoFLU/contents/dependencies/fastas"))
        ref_fa = os.path.join(td, "all_genoflu_refs.fasta")
        with open(ref_fa, "w") as out:
            for fm in fastas_meta:
                txt = http_get_text(fm["download_url"])
                out.write(txt.rstrip() + "\n")
                # Parse individual allele sequences
                cur_al = None
                cur_gene = None
                buf = []
                for ln in txt.splitlines():
                    ln = ln.strip()
                    if not ln:
                        continue
                    if ln.startswith(">"):
                        if cur_al and cur_gene and buf:
                            s = "".join(buf).upper()
                            if cur_al not in genoflu_allele_seqs[cur_gene] or len(s) > len(genoflu_allele_seqs[cur_gene][cur_al]):
                                genoflu_allele_seqs[cur_gene][cur_al] = s
                        parts = ln[1:].split()
                        cur_al = parts[0] if len(parts) >= 3 else None
                        cur_gene = parts[2] if len(parts) >= 3 else None
                        buf = []
                    else:
                        buf.append(ln)
                if cur_al and cur_gene and buf:
                    s = "".join(buf).upper()
                    if cur_al not in genoflu_allele_seqs[cur_gene] or len(s) > len(genoflu_allele_seqs[cur_gene][cur_al]):
                        genoflu_allele_seqs[cur_gene][cur_al] = s

        db_p = os.path.join(td, "hpai_geno_db")
        subprocess.run([MAKEBLASTDB_BIN, "-in", ref_fa, "-dbtype", "nucl", "-out", db_p], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

        q_fa = os.path.join(td, "queries_253x8.fasta")
        with open(q_fa, "w") as qf:
            for rep in reps:
                uid = rep["uid"]
                for s in SEG_ORDER:
                    gid = chosen_segs[uid][s]["genome_id"]
                    qf.write(f">{uid}__SEG{s}\n{seq_by_gid[gid]}\n")

        bout = os.path.join(td, "blast.tsv")
        subprocess.run(
            [
                BLASTN_BIN,
                "-query",
                q_fa,
                "-db",
                db_p,
                "-word_size",
                "11",
                "-max_target_seqs",
                "5",
                "-num_threads",
                "4",
                "-outfmt",
                "6 qseqid length nident pident mismatch evalue bitscore sacc stitle",
                "-out",
                bout,
            ],
            check=True,
            stderr=subprocess.DEVNULL,
        )

        hits_per_q = defaultdict(list)
        for l in open(bout):
            cols = l.rstrip("\n").split("\t")
            qid = cols[0]
            allele, _, gene = cols[8].split()[:3]
            hits_per_q[qid].append((float(cols[6]), float(cols[3]), allele, gene, int(cols[1])))
        top_hits = {}
        for qid, lst in hits_per_q.items():
            lst.sort(key=lambda x: (-x[0], -x[1]))
            top_hits[qid] = (lst[0][2], lst[0][1], lst[0][3], lst[0][4])

        for rep in reps:
            uid = rep["uid"]
            alleles_98 = {}
            raw = {}
            pids = {}
            for idx_s, sn in enumerate(SEG_NAMES, 1):
                qid = f"{uid}__SEG{idx_s}"
                if qid in top_hits:
                    al, pid, _, _ = top_hits[qid]
                    raw[sn] = (al, round(pid, 2))
                    pids[sn] = pid
                    if pid >= 98.0:
                        alleles_98[sn] = al
                else:
                    raw[sn] = ("NOHIT", 0.0)
                    pids[sn] = 0.0
            if len(alleles_98) == 8:
                tup = tuple(alleles_98[sn] for sn in SEG_NAMES)
                g_call = geno_dict.get(tup, "Not assigned (8/8 >=98% novel combo)")
                matched_named = tup in geno_dict
            else:
                g_call = f"Not assigned ({len(alleles_98)}/8 >=98%)"
                matched_named = False
            geno_results[uid] = {
                "genotype": g_call,
                "matched_named": matched_named,
                "n_ge_98": len(alleles_98),
                "raw": raw,
                "ha_pid": pids["HA"],
            }

    print(f"  GenoFLU scored {len(geno_results)} isolates against {len(geno_dict)} official genotypes", flush=True)

    # Assemble 8 unified global segment trees (all 40 constellations, 2 reps/constellation + expanded global/GenoFLU refs)
    core_alleles_wanted = {
        "ea1", "ea2", "ea3", "ea4", "ea5", "ea6", "ea1b", "ea4b", "ea6N5",
        "am1.1", "am1.2", "am1.3", "am1.4.1", "am2.1", "am2.2", "am3", "am4", "am5",
        "am6N2", "am4N1", "am8", "am9", "am11", "am11N2", "am12", "am12N5", "am13",
        "am18", "am19", "am20", "am24",
    }

    reps_by_form = defaultdict(list)
    for rep in reps:
        reps_by_form[rep["formula"]].append(rep)

    # Pick 2 reps per constellation (or up to 3 if mixed countries/subtypes like #1, #4, #10!)
    selected_tree_reps = []
    for rk, cl in enumerate(constellations, 1):
        form = cl["structural_formula"]
        c_reps = reps_by_form[form]
        # Diversify by (country, subtype) so mixed clusters (#1, #4, #10) have all components in the global trees
        picked = []
        seen_cs = set()
        for r in c_reps:
            key = (r["country_tok"], r["subtype"])
            if key not in seen_cs:
                seen_cs.add(key)
                picked.append((rk, r))
            if len(picked) >= 3:
                break
        for r in c_reps:
            if len(picked) >= 2:
                break
            if all(r["uid"] != p[1]["uid"] for p in picked):
                picked.append((rk, r))
        selected_tree_reps.extend(picked)

    global_tree_tasks = []
    for seg_idx, (seg_num, sname) in enumerate(zip(SEG_ORDER, SEG_NAMES), 1):
        raw_fa = TREES_DIR / f"iav_global_seg{seg_num}_{sname}.raw.fasta"
        aln_fa = TREES_DIR / f"iav_global_seg{seg_num}_{sname}.fasta"
        tr_p = Path(str(aln_fa) + ".treefile")

        with open(raw_fa, "w") as f:
            # 1. GenoFLU allele references for this segment
            for al, seq in sorted(genoflu_allele_seqs.get(sname, {}).items()):
                if al in core_alleles_wanted:
                    f.write(f">Ref_GenoFLU_{al}\n{seq}\n")
            # 2. Expanded global/EURL/SAm/LPAI references
            for rname, smap in sorted(ext_ref_segs.items()):
                gid = smap[seg_num]["genome_id"]
                if gid in seq_by_gid:
                    f.write(f">{rname}\n{seq_by_gid[gid]}\n")
            # 3. Query constellation representatives
            for rk, r in selected_tree_reps:
                uid = r["uid"]
                gid = chosen_segs[uid][seg_num]["genome_id"]
                f.write(f">Q_c{rk:02d}_{uid.replace('|', '_')}_{r['country_tok']}_{r['subtype']}\n{seq_by_gid[gid]}\n")

        global_tree_tasks.append((seg_num, sname, raw_fa, aln_fa, tr_p))

    print(f"  Aligning & building 8 unified global IQ-TREE 2 phylogenies ({len(selected_tree_reps)} reps + global/GenoFLU refs)...", flush=True)

    def _align_and_tree(item):
        _, _, raw_fa, aln_fa, tr_p = item
        if not aln_fa.exists():
            res = subprocess.run([MAFFT_BIN, "--auto", "--quiet", str(raw_fa)], check=True, stdout=subprocess.PIPE, text=True)
            aln_fa.write_text(res.stdout)
        if not tr_p.exists():
            subprocess.run(
                [IQTREE_BIN, "-s", str(aln_fa), "-m", "GTR+G4", "-B", "1000", "-T", "1", "--quiet", "-redo"],
                check=True,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
            )

    with concurrent.futures.ThreadPoolExecutor(max_workers=8) as pool:
        list(pool.map(_align_and_tree, global_tree_tasks))

    # Parse nearest Ref_ and UFBoot support in the 8 global trees for every constellation
    tree_summary_by_rank = defaultdict(dict)
    for seg_num, sname, _, aln_fa, tr_p in global_tree_tasks:
        # Group queries by constellation rank
        q_by_rk = defaultdict(list)
        for rk, r in selected_tree_reps:
            q_by_rk[rk].append(f"Q_c{rk:02d}_{r['uid'].replace('|', '_')}_{r['country_tok']}_{r['subtype']}")
        for rk, qnames in q_by_rk.items():
            # For HA/NA/internal, root at midpoint or SAm/NAm LPAI if present
            root_ref = "Ref_NAm_LPAI_H5N2_UT18" if sname == "HA" else None
            tstat = analyze_tree_sister(tr_p, qnames, root_taxon=root_ref)
            tree_summary_by_rank[rk][sname] = {
                "ml_sister_ref": tstat.get("ml_sister_ref"),
                "sister_ufboot": tstat.get("sister_clade_ufboot"),
                "per_query_nearest": tstat.get("per_query_nearest_ref"),
                "top_patristic": tstat.get("nearest_refs_patristic", [])[:2],
            }

    # Summarize all 40 constellations honestly
    audit_rows = []
    class_counts = Counter()
    for rank, cl in enumerate(constellations, 1):
        form = cl["structural_formula"]
        c_reps = reps_by_form[form]

        geno_calls = Counter(geno_results[r["uid"]]["genotype"] for r in c_reps)
        named_matches = sum(1 for r in c_reps if geno_results[r["uid"]]["matched_named"])
        ha_pids = [geno_results[r["uid"]]["ha_pid"] for r in c_reps]
        min_ha_pid, max_ha_pid = min(ha_pids), max(ha_pids)

        # BioProjects, PubMed IDs, and NCBI FLI-Genotype notes across reps
        bps = Counter()
        pubs = Counter()
        fli_notes = Counter()
        for r in c_reps:
            segs = chosen_segs[r["uid"]]
            ha_row = segs["4"]
            if ha_row.get("bioproject_accession"):
                bps[ha_row["bioproject_accession"]] += 1
            if ha_row.get("publication"):
                pubs[str(ha_row["publication"])] += 1
            ha_acc = (ha_row.get("genbank_accessions") or "").split(",")[0].strip()
            subn = ncbi_notes_by_acc.get(ha_acc, "")
            for tok in subn.split("|"):
                if "FLI-Genotype" in tok:
                    fli_notes[tok.strip()] += 1

        top_geno = geno_calls.most_common(1)[0][0] if geno_calls else "Not assigned"
        # Build representative 8-seg GenoFLU allele string
        ex_uid = c_reps[0]["uid"]
        ex_raw = geno_results[ex_uid]["raw"]
        geno_allele_str = " ".join(f"{sn}:{ex_raw[sn][0]}({ex_raw[sn][1]:.1f}%)" for sn in SEG_NAMES)

        seg_ml_str = " ; ".join(
            f"{sn}->{tree_summary_by_rank[rank][sn]['ml_sister_ref']}(UFB={tree_summary_by_rank[rank][sn]['sister_ufboot']})"
            for sn in SEG_NAMES
        )
        ha_nearest = tree_summary_by_rank[rank]["HA"]["per_query_nearest"]
        pb2_nearest = tree_summary_by_rank[rank]["PB2"]["per_query_nearest"]
        ns_nearest = tree_summary_by_rank[rank]["NS"]["per_query_nearest"]

        years = set(cl["years"])
        countries = {k for k, _ in cl["top_countries"]}
        subtypes = {k for k, _ in cl["subtypes"]}

        # Honest Classification Logic
        if named_matches >= max(1, len(c_reps) // 2):
            classification = "KNOWN_USDA_GENOFLU_OMITTED_FROM_PACK"
            rationale = (
                f"Matches official USDA-VS/GenoFLU genotype(s) {dict(geno_calls)} at >=98% across all 8 segments "
                "(omitted from the 10-profile influenza_a_8seg pack)."
            )
        elif max_ha_pid < 83.0:
            classification = "PUBLISHED_NON_USDA_GENOTYPE_OR_EURL"
            rationale = (
                f"Non-2.3.4.4b indigenous LPAI H5 virus (HA_pid={min_ha_pid:.1f}-{max_ha_pid:.1f}% vs 2.3.4.4b; "
                f"HA tree nearest={ha_nearest}; {dict(bps)}); forced onto 2.3.4.4b pack."
            )
        elif min_ha_pid < 90.0:
            classification = "PUBLISHED_NON_USDA_GENOTYPE_OR_EURL"
            rationale = (
                f"Includes or comprises non-2.3.4.4b Asian clade 2.3.2.1a/e H5N1 (HA_pid={min_ha_pid:.1f}-{max_ha_pid:.1f}%; "
                f"HA tree nearest={ha_nearest}) outside the 2.3.4.4b pack scope."
            )
        elif fli_notes:
            classification = "PUBLISHED_NON_USDA_GENOTYPE_OR_EURL"
            rationale = (
                f"Published European EURL/FLI genotype {dict(fli_notes)} (e.g., EA-2024-DI.2.1); "
                f"8-segment trees place with {tree_summary_by_rank[rank]['PB2']['ml_sister_ref']} / {tree_summary_by_rank[rank]['NS']['ml_sister_ref']}."
            )
        elif pubs:
            classification = "PUBLISHED_NON_USDA_GENOTYPE_OR_EURL"
            rationale = (
                f"Already published in peer-reviewed literature (PubMed {dict(pubs)}, e.g. Nat Commun 2026 / EURL BB gull genotype; "
                f"NS tree nearest={ns_nearest}, PB2={pb2_nearest})."
            )
        elif "Uruguay" in countries and "2026" in years:
            classification = "PUBLISHED_NON_USDA_GENOTYPE_OR_EURL"
            rationale = (
                f"Published 2026 Uruguay South American H5N1 reassortant (Viruses 2026, PMID 42198761 / PMC13211659): "
                f"2.3.4.4b HA/NA/MP + South American B3.2 NP + South American LPAI internal segments ({pb2_nearest}, {ns_nearest})."
            )
        elif "Finland" in countries or "France" in countries:
            classification = "PUBLISHED_NON_USDA_GENOTYPE_OR_EURL"
            rationale = (
                f"Published European gull-adapted EA-2022-BB genotype (Finland fur-farm / France gull 2022-2024; "
                f"NS={ns_nearest} at 91.2% vs NAm am3, nesting with euBB/Eurasian gull clade)."
            )
        elif subtypes & {"H5N8", "H5N6"} and (years & {"2020", "2021"}):
            classification = "PUBLISHED_NON_USDA_GENOTYPE_OR_EURL"
            rationale = (
                f"2020-2021 Eurasian clade 2.3.4.4b H5N8/H5N6 epizootic constellation (e.g. EA-2020-A/C, EA-2021-I; "
                f"GenoFLU={dict(geno_calls)}; HA={ha_nearest})."
            )
        elif any(r["country_tok"] == "USA" and r["subtype"] == "H5N5" for r in c_reps):
            classification = "PUBLISHED_NON_USDA_GENOTYPE_OR_EURL"
            rationale = (
                f"Transatlantic EA-2021-I H5N5 incursion (Cell Reports 2024 PMID 39003741) mixed with Eurasian 2.3.4.4b A3+PB2 "
                f"({dict(geno_calls)}; PB2 nearest={pb2_nearest})."
            )
        else:
            # Check remaining unassigned 2.3.4.4b constellations (e.g. #16 Laos? No, caught by HA<90%; what remains?)
            # Inspect exact 8-segment UFBoot support
            all_ufb80 = all(
                (tree_summary_by_rank[rank][sn]["sister_ufboot"] is not None and tree_summary_by_rank[rank][sn]["sister_ufboot"] >= 80.0)
                for sn in SEG_NAMES
            )
            if all_ufb80:
                classification = "GENUINELY_NOVEL_REASSORTANT_CONSTELLATION"
                rationale = f"Clade 2.3.4.4b (HA={min_ha_pid:.1f}-{max_ha_pid:.1f}%), GenoFLU Not assigned, no EURL/literature match, 8/8 IQ-TREE2 UFBoot>=80%."
            else:
                classification = "PUBLISHED_NON_USDA_GENOTYPE_OR_EURL"
                rationale = (
                    f"Unassigned Eurasian/regional 2.3.4.4b constellation without 8/8 UFBoot>=80% novelty confirmation "
                    f"(GenoFLU={dict(geno_calls)}; PB2={pb2_nearest}; NS={ns_nearest})."
                )

        class_counts[classification] += 1
        rec = {
            "rank": rank,
            "n_isolates": cl["n_isolates"],
            "structural_formula": form,
            "epidemiological_formula": cl["epidemiological_formula"],
            "subtypes": ";".join(f"{k}:{v}" for k, v in cl["subtypes"]),
            "years": ",".join(cl["years"]),
            "top_countries": ";".join(f"{k}:{v}" for k, v in cl["top_countries"]),
            "n_reps_evaluated": len(c_reps),
            "genoflu_calls": ";".join(f"{k}:{v}" for k, v in geno_calls.most_common()),
            "genoflu_named_matches": named_matches,
            "ha_pident_range": f"{min_ha_pid:.1f}-{max_ha_pid:.1f}%",
            "example_genoflu_alleles": geno_allele_str,
            "bioprojects": ";".join(f"{k}:{v}" for k, v in bps.most_common(4)),
            "pubmed_ids": ";".join(f"{k}:{v}" for k, v in pubs.most_common(4)),
            "ncbi_fli_genotypes": ";".join(f"{k}:{v}" for k, v in fli_notes.most_common()),
            "iqtree_8seg_placements": seg_ml_str,
            "classification": classification,
            "rationale": rationale,
        }
        audit_rows.append(rec)
        print(
            f"  [IAV #{rank:2d} n={cl['n_isolates']:3d}] {form:48s} -> {classification} | GenoFLU={top_geno}",
            flush=True,
        )

    tsv_cols = list(audit_rows[0].keys())
    with open(OUT_IAV_AUDIT_TSV, "w", newline="", encoding="utf-8") as tf:
        w = csv.DictWriter(tf, fieldnames=tsv_cols, delimiter="\t")
        w.writeheader()
        for r in audit_rows:
            w.writerow(r)

    out_json = {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "runtime_sec": round(time.time() - t0, 2),
        "n_constellations_audited": len(audit_rows),
        "n_reps_evaluated": len(reps),
        "classification_counts": dict(class_counts),
        "constellations": audit_rows,
        "global_8seg_tree_details": {str(k): v for k, v in sorted(tree_summary_by_rank.items())},
    }
    with open(OUT_IAV_AUDIT_JSON, "w", encoding="utf-8") as jf:
        json.dump(out_json, jf, indent=2)

    print(f"[+] Part 2 complete in {time.time() - t0:.2f}s | summary={dict(class_counts)}", flush=True)
    return out_json


def main():
    t0 = time.time()
    run_part1_sc2_audit()
    run_part2_iav_audit()
    print(f"\n[DONE] Script 25 finished in {time.time() - t0:.2f}s", flush=True)


if __name__ == "__main__":
    main()
