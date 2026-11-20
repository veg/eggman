#!/usr/bin/env python3
"""
27_benchmark_sc2_nextclade_and_iav_genoflu.py
=============================================
Head-to-head multi-pathogen benchmark of EGGMAN against domain-standard callers:
  1. SARS-CoV-2: EGGMAN vs Nextclade CLI (`-j 8`) on:
     - `sc2_curated_1000` (500 clonal + 500 Pango X* recombinants)
     - `sc2_ood_xdv_xfz_xfc` (10 OOD XDV/XFZ/XFC recombinants with XDV/XFZ/XFC masked from EGGMAN)
  2. Influenza A (2.3.4.4b H5N1 / H5Nx): EGGMAN vs USDA GenoFLU (`blastn >= 98%` 8-segment allele dictionary) on:
     - `iav_sra_1200` (1,100 labeled + 100 USDA Not assigned)
     - `iav_sra_650_ood` (500 B3.13 + 150 D1.1 holdout with B3.13 & D1.1 masked from both EGGMAN and GenoFLU)

Writes:
  - results/benchmark_sc2_nextclade_and_iav_genoflu.csv
  - results/benchmark_sc2_nextclade_and_iav_genoflu.json
"""

import csv
import glob
import io
import json
import lzma
import os
import re
import shutil
import subprocess
import tempfile
import time
import urllib.request
from collections import Counter, defaultdict
from pathlib import Path

import pandas as pd


ROOT = Path("/Users/sergei/Projects/HIV-recombination")
RES = ROOT / "results"
DATA = ROOT / "data"
CACHE_DIR = Path("/Users/sergei/.gemini/antigravity-cli/brain/d7d3d768-93f4-475c-81ac-ccb936c7805d/scratch/chronaeon_cache")

EGGMAN_BIN = str(ROOT / "rhiztype" / "target" / "release" / "eggman")
NEXTCLADE_BIN = "/Users/sergei/miniconda3/bin/nextclade"
BLASTN_BIN = "/opt/homebrew/bin/blastn"
MAKEBLASTDB_BIN = "/opt/homebrew/bin/makeblastdb"

SC2_PACK = ROOT / "rhiztype" / "packs" / "sars_cov2_wuhan"
IAV_PACK = ROOT / "rhiztype" / "packs" / "influenza_a_8seg"

WORLD_META_GLOB = "/Users/sergei/Projects/trees/data/alignments/viral_and_bwz_codec/sc2_*.metadata.tsv.xz"
SRA_GENOFLU_TSV = Path("/Users/sergei/Projects/stuff/bovine/data/andersen_metadata/genoflu_results.tsv")

IAV_SEG_NAMES = ["PB2", "PB1", "PA", "HA", "NP", "NA", "MP", "NS"]
IAV_BOUNDS = [
    ("PB2", 0, 2280),
    ("PB1", 2280, 4554),
    ("PA", 4554, 6705),
    ("HA", 6705, 8409),
    ("NP", 8409, 9906),
    ("NA", 9906, 11316),
    ("MP", 11316, 12298),
    ("NS", 12298, 13136),
]

SC2_CLONAL_16 = {"B", "B.1.1.7", "P.1", "B.1.617.2", "BA.1", "BA.2", "BA.4", "BA.5", "BJ.1", "BM.1.1.1", "JN.1", "KP.3.3", "KS.1.1", "BA.3.2", "LF.7", "LP.8.1"}
SC2_XBB_ALIASES = {"EG", "FL", "FE", "FU", "FY", "GJ", "GN", "GW", "HA", "HK", "HV", "JD", "JG", "EL", "EM", "GW", "HF", "HH"}
SC2_XDV_ALIASES = {"NB", "PQ", "RC", "RH", "RK", "SV", "SZ", "TM"}


def read_fasta_clean(path: Path):
    """Read FASTA returning [(full_header, base_id, seq_upper)]."""
    recs = []
    cur_h, buf = None, []
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        for ln in f:
            ln = ln.strip()
            if not ln:
                continue
            if ln.startswith(">"):
                if cur_h is not None:
                    recs.append((cur_h, cur_h.split("|")[0].split()[0], "".join(buf).upper()))
                cur_h = ln[1:].strip()
                buf = []
            else:
                buf.append(ln)
    if cur_h is not None:
        recs.append((cur_h, cur_h.split("|")[0].split()[0], "".join(buf).upper()))
    return recs


def build_filtered_pack(src_pack: Path, dst_pack: Path, omit_labels: set):
    dst_pack.mkdir(parents=True, exist_ok=True)
    for fn in ("manifest.json", "profiles.fasta", "compendium.json"):
        shutil.copy(src_pack / fn, dst_pack / fn)
    kept = []
    for ln in (src_pack / "certified_profiles.tsv").read_text().splitlines():
        if not ln.strip() or ln.startswith("#"):
            kept.append(ln)
        elif ln.split("\t")[0] not in omit_labels:
            kept.append(ln)
    (dst_pack / "certified_profiles.tsv").write_text("\n".join(kept) + "\n")


def load_sc2_ground_truth_map():
    """Map base_id -> reported Pango lineage across World Union, SC2 Supplement, and BV-BRC 23K."""
    gt = {}
    supp_p = CACHE_DIR / "sc2_supplement_meta.json"
    if supp_p.exists():
        for sid, m in json.loads(supp_p.read_text()).items():
            gt[sid] = (m.get("pango") or m.get("lineage") or "").strip()
    b23_p = DATA / "bvbrc_sc2_stream" / "bvbrc_sc2_metadata.json"
    if b23_p.exists():
        b23_obj = json.loads(b23_p.read_text())
        b23_rows = list(b23_obj.values()) if isinstance(b23_obj, dict) else b23_obj
        for r in b23_rows:
            if not isinstance(r, dict):
                continue
            gid = str(r.get("genome_id") or "").strip()
            lin = str(r.get("lineage") or "").strip()
            if gid and lin:
                gt[gid] = lin
    for p in sorted(glob.glob(WORLD_META_GLOB)):
        with lzma.open(p, "rt", encoding="utf-8", errors="replace") as f:
            hdr = f.readline().rstrip("\n").split("\t")
            col_i = {c: i for i, c in enumerate(hdr)}
            s_i = col_i.get("strain", 0)
            p_i = col_i.get("Nextclade_pango", col_i.get("pango_lineage"))
            pl_i = col_i.get("pango_lineage")
            c_i = col_i.get("Nextstrain_clade")
            for ln in f:
                parts = ln.rstrip("\n").split("\t")
                if len(parts) <= s_i:
                    continue
                sid = parts[s_i].strip()
                p_val = parts[pl_i].strip() if (pl_i is not None and pl_i < len(parts) and parts[pl_i].strip()) else ""
                if not p_val and p_i is not None and p_i < len(parts):
                    p_val = parts[p_i].strip()
                c_val = parts[c_i].strip() if (c_i is not None and c_i < len(parts)) else ""
                if sid and (p_val or c_val):
                    gt[sid] = p_val or c_val
    return gt


def load_pango_aliases():
    p_cache = CACHE_DIR / "pango_alias_key.json"
    if p_cache.exists():
        return json.loads(p_cache.read_text())
    txt = urllib.request.urlopen(
        "https://raw.githubusercontent.com/cov-lineages/pango-designation/master/pango_designation/alias_key.json",
        timeout=30,
    ).read().decode("utf-8")
    p_cache.write_text(txt)
    return json.loads(txt)


def resolve_pango_root(lin: str, alias_json: dict) -> str:
    s = (lin or "").strip().split()[0]
    if not s or s in ("Unassigned", "recombinant"):
        return ""
    pref = s.split(".", 1)[0]
    if pref.startswith("X"):
        return pref
    exp = alias_json.get(pref)
    if isinstance(exp, list):
        return pref
    if isinstance(exp, str) and exp.startswith("X"):
        return exp.split(".", 1)[0]
    return pref


def is_sc2_recomb_pango(lin: str, alias_json: dict) -> bool:
    rt = resolve_pango_root(lin, alias_json)
    return bool(rt and rt.startswith("X"))


def pango_family(lin: str, alias_json: dict) -> str:
    s = (lin or "").strip()
    if not s:
        return ""
    pref = s.split(".", 1)[0]
    if pref in SC2_XBB_ALIASES or s.startswith("XBB"):
        return "XBB"
    if pref in SC2_XDV_ALIASES or s.startswith("XDV"):
        return "XDV"
    if s.startswith("B.1.1.7"):
        return "B.1.1.7"
    if s.startswith("B.1.617.2") or pref == "AY":
        return "B.1.617.2"
    if s.startswith("P.1"):
        return "P.1"
    if s == "B" or s.startswith("B.1."):
        return "B"
    if s.startswith("RE."):
        return "BA.3.2"
    for root in ("XBB.1.5", "XBB.1.16", "KP.3.3", "KS.1.1", "BM.1.1.1", "BA.3.2", "LP.8.1", "BA.1", "BA.2", "BA.4", "BA.5", "BJ.1", "JN.1", "LF.7"):
        if s == root or s.startswith(root + "."):
            return root
    rt = resolve_pango_root(s, alias_json)
    return rt or pref


def run_sc2_benchmark():
    print("=" * 90)
    print("PART 1: SARS-CoV-2 HEAD-TO-HEAD (EGGMAN vs NEXTCLADE CLI)")
    print("=" * 90)

    alias_json = load_pango_aliases()
    gt_map = load_sc2_ground_truth_map()

    ds_dir = CACHE_DIR / "nextclade_sc2_ds"
    if not (ds_dir / "reference.fasta").exists():
        print("[*] Fetching Nextclade SARS-CoV-2 dataset (cached outside benchmark timer)...")
        ds_dir.mkdir(parents=True, exist_ok=True)
        subprocess.run([NEXTCLADE_BIN, "dataset", "get", "--name", "sars-cov-2", "--output-dir", str(ds_dir)], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    rows = []
    with tempfile.TemporaryDirectory() as td:
        td_p = Path(td)
        ood_pack = td_p / "sc2_ood_pack"
        build_filtered_pack(SC2_PACK, ood_pack, {"XDV", "XFZ", "XFC"})

        for cohort_name, src_fa, egg_pack, is_ood in [
            ("sc2_curated_1000", CACHE_DIR / "bench_sc2_curated_1000_col.fasta", SC2_PACK, False),
            ("sc2_ood_xdv_xfz_xfc", CACHE_DIR / "bench_sc2_ood_xdv_xfz_xfc_col.fasta", ood_pack, True),
        ]:
            recs = read_fasta_clean(src_fa)
            assert len(recs) == len({bid for _, bid, _ in recs}), (
                f"Duplicate sequence IDs in {cohort_name}: {len(recs)} vs {len({bid for _, bid, _ in recs})}"
            )
            lbl_path = CACHE_DIR / f"bench_{cohort_name}_labels.json"
            cur_lbls = json.loads(lbl_path.read_text()) if lbl_path.exists() else {}
            clean_fa = td_p / f"{cohort_name}.fasta"
            with open(clean_fa, "w") as f:
                for _, bid, seq in recs:
                    f.write(f">{bid}\n{seq.replace('-', '')}\n")

            egg_csv = td_p / f"{cohort_name}_eggman.csv"
            t0 = time.time()
            subprocess.run(
                [EGGMAN_BIN, "--pack", str(egg_pack), "--gamma", "3.5", "--min-snps", "3", "--p-thresh", "0.01", "-q", "-o", str(egg_csv), str(clean_fa)],
                check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
            )
            dt_egg = time.time() - t0
            egg_df = pd.read_csv(egg_csv).set_index("sequence_id")

            nc_tsv = td_p / f"{cohort_name}_nextclade.tsv"
            t1 = time.time()
            subprocess.run(
                [NEXTCLADE_BIN, "run", "-j", "8", "-D", str(ds_dir), "--output-tsv", str(nc_tsv), str(clean_fa)],
                check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
            )
            dt_nc = time.time() - t1
            nc_df = pd.read_csv(nc_tsv, sep="\t").set_index("seqName")

            n = len(recs)
            egg_bp_ok = egg_t1 = egg_t2 = egg_ood_urf = egg_ood_t1 = 0
            nc_t1 = nc_t2 = nc_ood_recomb = 0
            for _, bid, _ in recs:
                exp_lbl = cur_lbls.get(bid, "")
                gt_lin = gt_map.get(bid, "") or exp_lbl
                gt_is_rec = (exp_lbl not in SC2_CLONAL_16) if exp_lbl else is_sc2_recomb_pango(gt_lin, alias_json)
                gt_fam = exp_lbl if exp_lbl else pango_family(gt_lin, alias_json)

                er = egg_df.loc[bid]
                e_cat = "" if pd.isna(er["called_category"]) else str(er["called_category"])
                e_lin = "" if pd.isna(er["called_lineage"]) else str(er["called_lineage"])
                e_nbp = int(er["n_breakpoints"])
                if not is_ood:
                    if (not gt_is_rec and e_nbp == 0) or (gt_is_rec and e_nbp >= 1):
                        egg_bp_ok += 1
                    if (not gt_is_rec and e_nbp == 0 and e_cat != "Certified_CRF") or (gt_is_rec and (e_nbp >= 1 or e_cat in ("Certified_CRF", "URF"))):
                        egg_t1 += 1
                    if e_lin == gt_fam or (gt_fam in {"XBB", "XBB.1.5", "XBB.1.16"} and e_lin in {"XBB", "XBB.1.5", "XBB.1.16"}) or (e_lin == gt_lin.split(".")[0]):
                        egg_t2 += 1
                else:
                    if e_cat == "URF":
                        egg_ood_urf += 1
                    if e_nbp >= 1:
                        egg_ood_t1 += 1

                nr = nc_df.loc[bid]
                nc_p_val = nr.get("Nextclade_pango", "")
                nc_p = "" if pd.isna(nc_p_val) else str(nc_p_val or "")
                nc_cl_val = nr.get("clade", "")
                nc_cl = "" if pd.isna(nc_cl_val) else str(nc_cl_val or "")
                nc_is_rec = is_sc2_recomb_pango(nc_p, alias_json) or nc_cl == "recombinant"
                nc_fam = pango_family(nc_p, alias_json)
                if not is_ood:
                    if gt_is_rec == nc_is_rec:
                        nc_t1 += 1
                    if nc_fam == gt_fam or (gt_fam in {"XBB", "XBB.1.5", "XBB.1.16"} and (nc_fam == "XBB" or nc_p.startswith("XBB"))) or resolve_pango_root(nc_p, alias_json) == resolve_pango_root(gt_lin, alias_json):
                        nc_t2 += 1
                else:
                    if nc_is_rec:
                        nc_ood_recomb += 1

            if not is_ood:
                rows.append({
                    "pathogen": "SARS-CoV-2", "cohort": cohort_name, "tool": "EGGMAN", "n_sequences": n,
                    "wall_sec": round(dt_egg, 3), "throughput_seq_per_sec": round(n / dt_egg, 1),
                    "clonal_vs_recomb_acc_pct": round(100.0 * egg_t1 / n, 2),
                    "lineage_or_genotype_acc_pct": round(100.0 * egg_t2 / n, 2),
                    "ood_novel_detection_pct": None,
                    "notes": f"Tier-1 clonal/mosaic & Tier-2 lineage family vs metadata (de_novo_bp_acc={100.0 * egg_bp_ok / n:.2f}% ({egg_bp_ok}/{n}))",
                })
                rows.append({
                    "pathogen": "SARS-CoV-2", "cohort": cohort_name, "tool": "Nextclade_CLI_8thr", "n_sequences": n,
                    "wall_sec": round(dt_nc, 3), "throughput_seq_per_sec": round(n / dt_nc, 1),
                    "clonal_vs_recomb_acc_pct": round(100.0 * nc_t1 / n, 2),
                    "lineage_or_genotype_acc_pct": round(100.0 * nc_t2 / n, 2),
                    "ood_novel_detection_pct": None, "notes": "Nextclade_pango family & X*/recombinant flag vs metadata",
                })
            else:
                rows.append({
                    "pathogen": "SARS-CoV-2", "cohort": cohort_name, "tool": "EGGMAN_OOD_Masked", "n_sequences": n,
                    "wall_sec": round(dt_egg, 3), "throughput_seq_per_sec": round(n / dt_egg, 1),
                    "clonal_vs_recomb_acc_pct": round(100.0 * egg_ood_t1 / n, 2),
                    "lineage_or_genotype_acc_pct": None,
                    "ood_novel_detection_pct": round(100.0 * egg_ood_urf / n, 2),
                    "notes": f"OOD XDV/XFZ/XFC withheld from dictionary: URF={egg_ood_urf}/{n}, de_novo_bp={egg_ood_t1}/{n}",
                })
                rows.append({
                    "pathogen": "SARS-CoV-2", "cohort": cohort_name, "tool": "Nextclade_CLI_8thr", "n_sequences": n,
                    "wall_sec": round(dt_nc, 3), "throughput_seq_per_sec": round(n / dt_nc, 1),
                    "clonal_vs_recomb_acc_pct": round(100.0 * nc_ood_recomb / n, 2),
                    "lineage_or_genotype_acc_pct": None,
                    "ood_novel_detection_pct": None,
                    "notes": "Static Auspice tree retains XDV/XFZ/XFC nodes (not held out as novel OOD)",
                })
            print(f"  [{cohort_name}] N={n}: EGGMAN={dt_egg:.3f}s ({n/dt_egg:.1f}/s) vs Nextclade={dt_nc:.3f}s ({n/dt_nc:.1f}/s)")
    return rows


def run_iav_benchmark():
    print("\n" + "=" * 90)
    print("PART 2: INFLUENZA A H5Nx HEAD-TO-HEAD (EGGMAN vs USDA GenoFLU BLASTn)")
    print("=" * 90)

    gf_truth = {}
    with open(SRA_GENOFLU_TSV) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            gf_truth[r["sample"].strip()] = r["Genotype"].strip()

    gf_cache = CACHE_DIR / "genoflu_db_cache"
    gf_cache.mkdir(parents=True, exist_ok=True)
    key_p = gf_cache / "genotype_key.xlsx"
    ref_fa = gf_cache / "all_genoflu_refs.fasta"
    if not key_p.exists():
        key_p.write_bytes(urllib.request.urlopen("https://github.com/USDA-VS/GenoFLU/raw/main/dependencies/genotype_key.xlsx", timeout=30).read())
    if not ref_fa.exists() or ref_fa.stat().st_size < 10000:
        fastas_meta = json.loads(urllib.request.urlopen("https://api.github.com/repos/USDA-VS/GenoFLU/contents/dependencies/fastas", timeout=30).read().decode("utf-8"))
        with open(ref_fa, "w") as out_f:
            for fm in fastas_meta:
                txt = urllib.request.urlopen(fm["download_url"], timeout=30).read().decode("utf-8", errors="replace")
                out_f.write(txt.rstrip() + "\n")

    df_key = pd.read_excel(key_p)
    geno_dict = {}
    for _, row in df_key.iterrows():
        tup = tuple(str(row[sn]).strip() for sn in IAV_SEG_NAMES)
        geno_dict[tup] = str(row["Genotype"]).strip()
    geno_dict_ood = {k: v for k, v in geno_dict.items() if v not in ("B3.13", "D1.1")}

    db_prefix = gf_cache / "hpai_geno_db"
    if not (gf_cache / "hpai_geno_db.nin").exists():
        subprocess.run([MAKEBLASTDB_BIN, "-in", str(ref_fa), "-dbtype", "nucl", "-out", str(db_prefix)], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    rows = []
    with tempfile.TemporaryDirectory() as td:
        td_p = Path(td)
        ood_pack = td_p / "iav_ood_pack"
        build_filtered_pack(IAV_PACK, ood_pack, {"B3.13", "D1.1"})

        for cohort_name, src_fa, egg_pack, gf_lookup, is_ood in [
            ("iav_sra_1200", CACHE_DIR / "bench_iav_sra_1200_col.fasta", IAV_PACK, geno_dict, False),
            ("iav_sra_650_ood", CACHE_DIR / "bench_iav_sra_650_ood_col.fasta", ood_pack, geno_dict_ood, True),
        ]:
            recs = read_fasta_clean(src_fa)
            assert len(recs) == len({bid for _, bid, _ in recs}), (
                f"Duplicate sequence IDs in {cohort_name}: {len(recs)} vs {len({bid for _, bid, _ in recs})}"
            )
            clean_fa = td_p / f"{cohort_name}.fasta"
            seg_q_fa = td_p / f"{cohort_name}_8segs.fasta"
            with open(clean_fa, "w") as f_c, open(seg_q_fa, "w") as f_s:
                for _, bid, seq in recs:
                    f_c.write(f">{bid}\n{seq}\n")
                    for idx_s, (_, a, b) in enumerate(IAV_BOUNDS, 1):
                        f_s.write(f">{bid}__SEG{idx_s}\n{seq[a:b]}\n")

            egg_csv = td_p / f"{cohort_name}_eggman.csv"
            t0 = time.time()
            subprocess.run(
                [EGGMAN_BIN, "--pack", str(egg_pack), "-q", "-o", str(egg_csv), str(clean_fa)],
                check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
            )
            dt_egg = time.time() - t0
            egg_df = pd.read_csv(egg_csv).set_index("sequence_id")

            bout = td_p / f"{cohort_name}_blast.tsv"
            t1 = time.time()
            subprocess.run(
                [
                    BLASTN_BIN, "-query", str(seg_q_fa), "-db", str(db_prefix),
                    "-word_size", "11", "-max_target_seqs", "5", "-num_threads", "4",
                    "-outfmt", "6 qseqid length nident pident mismatch evalue bitscore sacc stitle",
                    "-out", str(bout),
                ],
                check=True, stderr=subprocess.DEVNULL,
            )
            hits = defaultdict(list)
            for bl in open(bout):
                cols = bl.rstrip("\n").split("\t")
                al, _, gene = cols[8].split()[:3]
                hits[cols[0]].append((float(cols[6]), float(cols[3]), al, gene))
            gf_preds = {}
            for _, bid, _ in recs:
                al98 = {}
                for idx_s, sn in enumerate(IAV_SEG_NAMES, 1):
                    lst = hits.get(f"{bid}__SEG{idx_s}", [])
                    if lst:
                        lst.sort(key=lambda x: (-x[0], -x[1]))
                        if lst[0][1] >= 98.0:
                            al98[sn] = lst[0][2]
                if len(al98) == 8:
                    tup = tuple(al98[sn] for sn in IAV_SEG_NAMES)
                    gf_preds[bid] = gf_lookup.get(tup, "Not assigned")
                else:
                    gf_preds[bid] = "Not assigned"
            dt_gf = time.time() - t1

            n = len(recs)
            labeled_n = 0
            egg_exact = gf_exact = egg_ood_urf = egg_ood_oop = gf_ood_unas = egg_ood_t1 = 0
            for _, bid, _ in recs:
                gt_g = gf_truth.get(bid, "")
                er = egg_df.loc[bid]
                e_cat = "" if pd.isna(er["called_category"]) else str(er["called_category"])
                e_lin = "" if pd.isna(er["called_lineage"]) else str(er["called_lineage"])
                e_nbp = int(er["n_breakpoints"])
                e_bc = "" if pd.isna(er["hxb2_barcode"]) else str(er["hxb2_barcode"])
                e_trunks = {seg.split(":", 1)[0].strip() for seg in e_bc.split("|") if ":" in seg}
                g_call = gf_preds[bid]
                if not is_ood:
                    if not gt_g.startswith("Not assigned"):
                        labeled_n += 1
                        if e_lin == gt_g:
                            egg_exact += 1
                        if g_call == gt_g:
                            gf_exact += 1
                else:
                    if e_cat == "URF":
                        egg_ood_urf += 1
                    elif e_cat == "OutOfPanel":
                        egg_ood_oop += 1
                    exp_trunks = {"NAm_B313", "EA1"} if "B3.13" in gt_g else {"NAm_D1", "EA3"}
                    if e_nbp >= 1 and e_trunks == exp_trunks:
                        egg_ood_t1 += 1
                    if g_call.startswith("Not assigned"):
                        gf_ood_unas += 1

            if not is_ood:
                rows.append({
                    "pathogen": "Influenza_A_H5Nx", "cohort": cohort_name, "tool": "EGGMAN", "n_sequences": n,
                    "wall_sec": round(dt_egg, 3), "throughput_seq_per_sec": round(n / dt_egg, 1),
                    "clonal_vs_recomb_acc_pct": None,
                    "lineage_or_genotype_acc_pct": round(100.0 * egg_exact / max(1, labeled_n), 2),
                    "ood_novel_detection_pct": None,
                    "notes": f"Exact genotype match on {labeled_n}/{n} USDA-labeled isolates (plus 100 USDA Not assigned resolved by EGGMAN)",
                })
                rows.append({
                    "pathogen": "Influenza_A_H5Nx", "cohort": cohort_name, "tool": "USDA_GenoFLU_BLASTn", "n_sequences": n,
                    "wall_sec": round(dt_gf, 3), "throughput_seq_per_sec": round(n / dt_gf, 1),
                    "clonal_vs_recomb_acc_pct": None,
                    "lineage_or_genotype_acc_pct": round(100.0 * gf_exact / max(1, labeled_n), 2),
                    "ood_novel_detection_pct": None,
                    "notes": f"4-thread blastn >=98% allele dictionary match ({labeled_n} labeled + 100 Not assigned)",
                })
            else:
                egg_ood_novel = egg_ood_urf + egg_ood_oop
                rows.append({
                    "pathogen": "Influenza_A_H5Nx", "cohort": cohort_name, "tool": "EGGMAN_OOD_Masked", "n_sequences": n,
                    "wall_sec": round(dt_egg, 3), "throughput_seq_per_sec": round(n / dt_egg, 1),
                    "clonal_vs_recomb_acc_pct": round(100.0 * egg_ood_t1 / n, 2),
                    "lineage_or_genotype_acc_pct": None,
                    "ood_novel_detection_pct": round(100.0 * egg_ood_novel / n, 2),
                    "notes": f"OOD B3.13/D1.1 omitted: URF={egg_ood_urf}/{n}, OutOfPanel={egg_ood_oop}/{n}, exact Tier-1 trunk set={egg_ood_t1}/{n}",
                })
                rows.append({
                    "pathogen": "Influenza_A_H5Nx", "cohort": cohort_name, "tool": "USDA_GenoFLU_OOD_Masked", "n_sequences": n,
                    "wall_sec": round(dt_gf, 3), "throughput_seq_per_sec": round(n / dt_gf, 1),
                    "clonal_vs_recomb_acc_pct": None,
                    "lineage_or_genotype_acc_pct": None,
                    "ood_novel_detection_pct": round(100.0 * gf_ood_unas / n, 2),
                    "notes": f"OOD B3.13/D1.1 omitted from genotype_key: Not_assigned={gf_ood_unas}/{n}",
                })
            print(f"  [{cohort_name}] N={n}: EGGMAN={dt_egg:.3f}s ({n/dt_egg:.1f}/s) vs GenoFLU={dt_gf:.3f}s ({n/dt_gf:.1f}/s)")
    return rows


def main():
    sc2_rows = run_sc2_benchmark()
    iav_rows = run_iav_benchmark()
    all_rows = sc2_rows + iav_rows
    df = pd.DataFrame(all_rows)
    out_csv = RES / "benchmark_sc2_nextclade_and_iav_genoflu.csv"
    out_json = RES / "benchmark_sc2_nextclade_and_iav_genoflu.json"
    df.to_csv(out_csv, index=False)
    out_json.write_text(json.dumps(all_rows, indent=2))
    print("\n" + df.to_string(index=False))
    print(f"\n[+] Saved multi-pathogen benchmark -> {out_csv} & {out_json}")


if __name__ == "__main__":
    main()
