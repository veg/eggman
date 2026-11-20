#!/usr/bin/env python3
"""
23_validate_usda_h5n1_and_sieve_controls.py
===========================================
Computes three empirical validation tables requested in BARS v2 review:
  1. Stage 1 parent-pair-weighted diagnostic-SNP density null and
     simulation-calibrated EGGMAN Tier-1 detectability null across HIV-1 coding zones
     (`results/sieve_stage1_snp_and_sim_null.csv`).
  2. Stage 2 joint multivariable logistic regression (`tat1 <-> tat2`, `RRE <-> rev2`,
     and `8150-8550` gp41 return-crossover hotspot) and hotspot-stratified tests
     (`results/sieve_stage2_tat_rre_multivariable.csv`).
  3. Orthogonal 8-segment IQ-TREE 2 (`GTR+G`, 1,000 UFBoot) validation of
     USDA GenoFlu "Not assigned" H5N1 field assemblies (`D1.8` and unindexed novel URF
     constellations) against donor reference alleles
     (`results/usda_h5n1_unassigned_iqtree_validation.csv`).
"""

import concurrent.futures
import csv
import io
import json
import os
import re
import shutil
import subprocess
import tempfile
import urllib.request
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import Phylo
from scipy import optimize, stats


ROOT = Path("/Users/sergei/Projects/HIV-recombination")
RES = ROOT / "results"
DATA = ROOT / "data"
CACHE_GENOFLU = Path("/Users/sergei/.gemini/antigravity-cli/brain/d7d3d768-93f4-475c-81ac-ccb936c7805d/scratch/chronaeon_cache/genoflu_db_cache")
EGGMAN_BIN = str(ROOT / "rhiztype" / "target" / "release" / "eggman")
IQTREE_BIN = "/Users/sergei/.local/bin/iqtree2"
MAFFT_BIN = shutil.which("mafft") or "/usr/local/bin/mafft"
BLASTN_BIN = "/opt/homebrew/bin/blastn"
MAKEBLASTDB_BIN = "/opt/homebrew/bin/makeblastdb"

FASTA_DIR_H5N1 = Path("/Users/sergei/Projects/stuff/bovine/data/andersen_avian_influenza/fasta")
GENOFLU_TSV = Path("/Users/sergei/Projects/stuff/bovine/data/andersen_metadata/genoflu_results.tsv")

H0, H1 = 790, 9417
CODING_LEN = H1 - H0 + 1  # 8,628 nt
HYPERVARIABLE_LOOPS = [(6615, 6692), (6693, 6812), (7110, 7217), (7377, 7478), (7602, 7634)]

ZONES = [
    ("gp120 Core (6450-7650)", [(6450, 7650)]),
    ("Gag-CA Core (1280-1850)", [(1280, 1850)]),
    ("Pol-RT Core (2650-3450)", [(2650, 3450)]),
    ("All 3 Rigid Cores Combined", [(1280, 1850), (2650, 3450), (6450, 7650)]),
    ("All 4 Domain Hinges Combined", [(1100, 1260), (1980, 2620), (4950, 6350), (8150, 8950)]),
    ("gag/pol & PR/RT Hinge (1980-2620)", [(1980, 2620)]),
]


def parse_bc(bc: str):
    segs = []
    for p in str(bc).split("|"):
        if ":" not in p:
            continue
        t, r = p.split(":", 1)
        if "-" in r and r.replace("-", "").isdigit():
            s, e = r.split("-", 1)
        elif "-" in t and t.replace("-", "").isdigit():
            s, e = t.split("-", 1)
            t = r
        else:
            continue
        segs.append((int(s), int(e), t.strip()))
    return segs


def trunk_at(segs, pos: int) -> str:
    for s, e, t in segs:
        if s <= pos <= e:
            return t
    return min(segs, key=lambda x: min(abs(pos - x[0]), abs(pos - x[1])))[2]


def has_bp_in(segs, lo: int, hi: int) -> int:
    for i in range(len(segs) - 1):
        bp = int(round(0.5 * (segs[i][1] + segs[i + 1][0])))
        if lo <= bp <= hi:
            return 1
    return 0


def bp_positions(segs_list):
    out = []
    for segs in segs_list:
        for i in range(len(segs) - 1):
            out.append(int(round(0.5 * (segs[i][1] + segs[i + 1][0]))))
    return np.array(out, dtype=int)


def load_cohorts():
    comp_df = pd.read_csv(RES / "hiv1_descriptive_nomenclature_compendium.csv")
    novel_df = pd.read_csv(RES / "genbank_novel_circulating_recombinants.csv")
    man_df = pd.read_csv(RES / "supp_manifest_sieve_749_lineages.csv")

    crf_comp = comp_df[comp_df["legacy_category"].str.startswith("CRF") & (comp_df["n_verified_segments"] > 1)]
    seen_labels = set(crf_comp["legacy_label"].astype(str).str.strip())
    est_segs = [parse_bc(b) for b in crf_comp["hxb2_barcode"]]

    bundle_path = ROOT / "rhiztype" / "crates" / "eggman-core" / "data" / "compendium_bundle.json"
    for ent in json.loads(bundle_path.read_text())["compendium"]:
        lbl = str(ent.get("label", "")).strip()
        if not str(ent.get("category", "")).startswith("CRF") or lbl in seen_labels:
            continue
        segs = parse_bc(ent.get("hxb2_barcode", ""))
        if len(segs) > 1:
            seen_labels.add(lbl)
            est_segs.append(segs)

    est_segs += [
        parse_bc(b)
        for b in novel_df[novel_df["lineage_status"].astype(str).str.startswith("Candidate_Regional_Mosaic")][
            "representative_hxb2_barcode"
        ]
    ]
    urf_segs = [parse_bc(b) for b in man_df[man_df["in_singleton_urf_461"] == True]["hxb2_barcode"]]
    return est_segs, urf_segs, man_df


def load_profile_seqs():
    seqs = {}
    cur, buf = None, []
    with open(DATA / "hxb2_15trunk_profile_9719.fasta") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if cur:
                    seqs[cur] = "".join(buf).upper()
                cur = line[1:].split()[0]
                buf = []
            else:
                buf.append(line)
        if cur:
            seqs[cur] = "".join(buf).upper()
    return seqs


def run_stage1_controls(est_segs, urf_segs, man_df):
    print("\n" + "=" * 90)
    print("PART 1: STAGE 1 DIAGNOSTIC-SNP DENSITY NULL & SIMULATION-CALIBRATED DETECTABILITY NULL")
    print("=" * 90)
    seqs = load_profile_seqs()

    # Build empirical adjacent parental pair weights across the 749-lineage manifest
    pair_weights = Counter()
    for bc in man_df["hxb2_barcode"]:
        segs = parse_bc(bc)
        for i in range(len(segs) - 1):
            a, b = segs[i][2], segs[i + 1][2]
            t1 = "Profile_Trunk_" + ("E" if a == "?E" else a)
            t2 = "Profile_Trunk_" + ("E" if b == "?E" else b)
            if t1 in seqs and t2 in seqs and t1 != t2:
                pair_weights[(min(t1, t2), max(t1, t2))] += 1

    # Empirical parent-pair-weighted diagnostic SNP profile across HXB2 790..9417
    loop_mask = np.zeros(CODING_LEN, dtype=bool)
    for lo, hi in HYPERVARIABLE_LOOPS:
        loop_mask[lo - H0 : hi - H0 + 1] = True

    snp_weight = np.zeros(CODING_LEN, dtype=np.float64)
    for (t1, t2), w in pair_weights.items():
        s1, s2 = seqs[t1], seqs[t2]
        has_e = (t1 == "Profile_Trunk_E") or (t2 == "Profile_Trunk_E")
        for idx in range(CODING_LEN):
            if loop_mask[idx]:
                continue
            h = H0 + idx
            if has_e and h < 5096:
                continue
            c1, c2 = s1[h - 1], s2[h - 1]
            if c1 in "ACGT" and c2 in "ACGT" and c1 != c2:
                snp_weight[idx] += w
    tot_snp_w = float(snp_weight.sum())

    # Simulate 1,000 single-crossover recombinants uniformly across HXB2 790..9417
    # with 2% terminal mutation noise on top of real Tier-1 parent consensus templates
    rng = np.random.default_rng(20261117)
    pairs = list(pair_weights.keys())
    p_weights = np.array([pair_weights[p] for p in pairs], dtype=np.float64)
    p_weights /= p_weights.sum()
    bases = ["A", "C", "G", "T"]
    hxb2_tpl = list(seqs["Profile_HXB2_K03455"])

    sim_records = []
    mut_rate = 0.02
    for sim_i in range(1000):
        t1, t2 = pairs[rng.choice(len(pairs), p=p_weights)]
        lo_bp = 5200 if ("Profile_Trunk_E" in (t1, t2)) else H0 + 140
        bp = int(rng.integers(lo_bp, H1 - 140))
        s1, s2 = seqs[t1], seqs[t2]
        arr = list(hxb2_tpl)
        for h in range(H0, H1 + 1):
            src = s1 if h <= bp else s2
            ch = src[h - 1]
            if ch not in "ACGT":
                ch = hxb2_tpl[h - 1]
            if ch in "ACGT" and rng.random() < mut_rate:
                ch = bases[(bases.index(ch) + int(rng.integers(1, 4))) % 4]
            arr[h - 1] = ch
        sim_records.append((f"SIM_{sim_i:04d}|bp={bp}|{t1.split('_')[-1]}_{t2.split('_')[-1]}", "".join(arr[H0 - 1 : H1])))

    with tempfile.NamedTemporaryFile("w", suffix=".fasta", delete=False) as tf:
        for hdr, s in sim_records:
            tf.write(f">{hdr}\n{s}\n")
        sim_fa = tf.name

    bps_urf = bp_positions(urf_segs)
    bps_est = bp_positions(est_segs)

    sim_regimes = [
        ("FullGenome_Screen_g13.5_p1e-5_s7_span200", ["--gamma", "13.5", "--min-snps", "7", "--min-span", "200", "--p-thresh", "1e-5"]),
        ("Default_g7.0", ["--gamma", "7.0"]),
    ]

    rows = []
    try:
        for reg_name, reg_flags in sim_regimes:
            out_csv = subprocess.check_output([EGGMAN_BIN] + reg_flags + ["-o", "/dev/stdout", "-q", sim_fa], text=True)
            sim_called_bps = []
            for r in csv.DictReader(io.StringIO(out_csv)):
                if int(r.get("n_breakpoints", 0) or 0) >= 1:
                    segs = parse_bc(r["hxb2_barcode"])
                    for i in range(len(segs) - 1):
                        sim_called_bps.append(int(round(0.5 * (segs[i][1] + segs[i + 1][0]))))
            sim_called_bps = np.array(sim_called_bps, dtype=int)
            print(f"[+] [{reg_name}] Simulated 1,000 uniform recombinants -> {len(sim_called_bps)} Tier-1 EGGMAN breakpoints called.")

            for z_name, ivs in ZONES:
                z_len = sum(e - s + 1 for s, e in ivs)
                unif_f = z_len / CODING_LEN
                snp_f = float(sum(snp_weight[s - H0 : e - H0 + 1].sum() for s, e in ivs)) / tot_snp_w

                sim_k = int(sum(np.sum((sim_called_bps >= s) & (sim_called_bps <= e)) for s, e in ivs))
                sim_f = sim_k / len(sim_called_bps)

                urf_k = int(sum(np.sum((bps_urf >= s) & (bps_urf <= e)) for s, e in ivs))
                urf_f = urf_k / len(bps_urf)

                est_k = int(sum(np.sum((bps_est >= s) & (bps_est <= e)) for s, e in ivs))
                est_f = est_k / len(bps_est)

                p_est_unif = float(stats.binomtest(est_k, len(bps_est), unif_f).pvalue)
                p_est_snp = float(stats.binomtest(est_k, len(bps_est), snp_f).pvalue)
                p_est_sim = float(stats.binomtest(est_k, len(bps_est), sim_f).pvalue)
                p_urf_snp = float(stats.binomtest(urf_k, len(bps_urf), snp_f).pvalue)
                p_urf_sim = float(stats.binomtest(urf_k, len(bps_urf), sim_f).pvalue)

                rows.append(
                    {
                        "sim_regime": reg_name,
                        "zone": z_name,
                        "nt_span": z_len,
                        "uniform_share": round(unif_f, 4),
                        "snp_density_null_share": round(snp_f, 4),
                        "sim_detectability_null_share": round(sim_f, 4),
                        "sim_called_k": sim_k,
                        "sim_called_n": len(sim_called_bps),
                        "urf_k": urf_k,
                        "urf_n": len(bps_urf),
                        "urf_share": round(urf_f, 4),
                        "urf_enrichment_vs_unif": round(urf_f / unif_f, 3),
                        "urf_enrichment_vs_snp": round(urf_f / snp_f, 3),
                        "urf_enrichment_vs_sim": round(urf_f / sim_f, 3),
                        "urf_p_vs_snp": p_urf_snp,
                        "urf_p_vs_sim": p_urf_sim,
                        "est_k": est_k,
                        "est_n": len(bps_est),
                        "est_share": round(est_f, 4),
                        "est_enrichment_vs_unif": round(est_f / unif_f, 3),
                        "est_enrichment_vs_snp": round(est_f / snp_f, 3),
                        "est_enrichment_vs_sim": round(est_f / sim_f, 3),
                        "est_p_vs_unif": p_est_unif,
                        "est_p_vs_snp": p_est_snp,
                        "est_p_vs_sim": p_est_sim,
                    }
                )
                print(
                    f"  * [{reg_name[:10]}] {z_name:34s}: Unif={100*unif_f:4.1f}%, SNP={100*snp_f:4.1f}%, Sim={100*sim_f:4.1f}% | "
                    f"Est={100*est_f:4.1f}% (vsSNP {est_f/snp_f:.2f}x p={p_est_snp:.2e}; vsSim {est_f/sim_f:.2f}x p={p_est_sim:.2e})"
                )
    finally:
        os.remove(sim_fa)

    out_path = RES / "sieve_stage1_snp_and_sim_null.csv"
    pd.DataFrame(rows).to_csv(out_path, index=False)
    print(f"[+] Saved Stage 1 SNP & Simulation Nulls -> {out_path}")


def run_stage2_controls(est_segs, urf_segs):
    print("\n" + "=" * 90)
    print("PART 2: STAGE 2 TAT vs RRE MULTIVARIABLE REGRESSION & 8150-8550 HOTSPOT CONDITIONING")
    print("=" * 90)

    records = []
    for y_val, cohort in [(1, est_segs), (0, urf_segs)]:
        for segs in cohort:
            tat = int(trunk_at(segs, 5950) == trunk_at(segs, 8400))
            rre = int(trunk_at(segs, 7880) == trunk_at(segs, 8500))
            bp83 = has_bp_in(segs, 8150, 8550)
            records.append((y_val, len(segs) - 1, tat, rre, bp83))
    arr = np.array(records, dtype=float)
    y = arr[:, 0]
    X = np.column_stack([np.ones(len(y)), arr[:, 1], arr[:, 2], arr[:, 3], arr[:, 4]])

    def nll(b):
        z = X @ b
        p = 1.0 / (1.0 + np.exp(-z))
        return -np.sum(y * np.log(p + 1e-15) + (1.0 - y) * np.log(1.0 - p + 1e-15))

    res = optimize.minimize(nll, np.zeros(X.shape[1]), method="BFGS")
    b = res.x
    p_hat = 1.0 / (1.0 + np.exp(-(X @ b)))
    W = p_hat * (1.0 - p_hat)
    cov = np.linalg.inv(X.T @ (W[:, None] * X))
    se = np.sqrt(np.diag(cov))
    z_stat = b / se
    p_val = 2.0 * (1.0 - stats.norm.cdf(np.abs(z_stat)))

    out_rows = []
    names = ["intercept", "n_breakpoints", "tat1_tat2_coancestry", "rre_rev2_coancestry", "has_bp_8150_8550"]
    for nm, bi, sei, zi, pi in zip(names, b, se, z_stat, p_val):
        out_rows.append(
            {
                "analysis": "Joint_Multivariable_Logistic",
                "stratum_or_term": nm,
                "n_established": int(y.sum()),
                "n_singleton_urf": int(len(y) - y.sum()),
                "est_rate_or_beta": round(float(bi), 5),
                "urf_rate_or_se": round(float(sei), 5),
                "odds_ratio": round(float(np.exp(bi)), 4),
                "z_stat": round(float(zi), 3),
                "p_value": float(pi),
            }
        )
        print(f"  [Logit] {nm:22s}: OR={np.exp(bi):.3f}, beta={bi:.4f}, SE={sei:.4f}, Z={zi:.2f}, p={pi:.4e}")

    for cond_val, cname in [(0, "No_BP_8150_8550"), (1, "Has_BP_8150_8550")]:
        e_sub = [s for s in est_segs if has_bp_in(s, 8150, 8550) == cond_val]
        u_sub = [s for s in urf_segs if has_bp_in(s, 8150, 8550) == cond_val]
        for fname, pos1, pos2 in [("tat1_tat2_5950_8400", 5950, 8400), ("rre_rev2_7880_8500", 7880, 8500)]:
            ek = sum(trunk_at(s, pos1) == trunk_at(s, pos2) for s in e_sub)
            uk = sum(trunk_at(s, pos1) == trunk_at(s, pos2) for s in u_sub)
            or_v, pv = stats.fisher_exact([[ek, len(e_sub) - ek], [uk, len(u_sub) - uk]])
            out_rows.append(
                {
                    "analysis": f"Stratum_{cname}",
                    "stratum_or_term": fname,
                    "n_established": len(e_sub),
                    "n_singleton_urf": len(u_sub),
                    "est_rate_or_beta": round(ek / len(e_sub), 4),
                    "urf_rate_or_se": round(uk / len(u_sub), 4),
                    "odds_ratio": round(float(or_v), 4),
                    "z_stat": np.nan,
                    "p_value": float(pv),
                }
            )
            print(
                f"  [{cname:16s} | {fname}] Est={ek}/{len(e_sub)} ({100*ek/len(e_sub):.1f}%) vs "
                f"URF={uk}/{len(u_sub)} ({100*uk/len(u_sub):.1f}%) | OR={or_v:.2f}, p={pv:.4e}"
            )

    out_path = RES / "sieve_stage2_tat_rre_multivariable.csv"
    pd.DataFrame(out_rows).to_csv(out_path, index=False)
    print(f"[+] Saved Stage 2 Multivariable & Hotspot Controls -> {out_path}")


def run_usda_h5n1_iqtree_validation():
    print("\n" + "=" * 90)
    print("PART 3: 8-SEGMENT IQ-TREE 2 VALIDATION OF USDA-UNASSIGNED H5N1 CONSTELLATIONS")
    print("=" * 90)

    # Load 10 USDA reference donor profiles & segment coordinates
    bounds = [
        ("PB2", 0, 2280),
        ("PB1", 2280, 4554),
        ("PA", 4554, 6705),
        ("HA", 6705, 8409),
        ("NP", 8409, 9906),
        ("NA", 9906, 11316),
        ("MP", 11316, 12298),
        ("NS", 12298, 13136),
    ]
    profiles = {}
    cur, buf = None, []
    with open(ROOT / "rhiztype" / "packs" / "influenza_a_8seg" / "profiles.fasta") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if cur:
                    profiles[cur] = "".join(buf).upper()
                cur = line[1:].split()[0]
                buf = []
            else:
                buf.append(line)
        if cur:
            profiles[cur] = "".join(buf).upper()

    # Collect complete 8-segment "Not assigned" SRA samples
    not_ass = []
    with open(GENOFLU_TSV) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if "Not assigned" in (r.get("Genotype") or ""):
                not_ass.append(r["sample"])

    segs_only = [s for s, _, _ in bounds]
    d18_candidates = []
    novel_candidates = []

    # Screen complete 8-seg isolates matching exact reference segment lengths (for ungapped clean MSA)
    with tempfile.NamedTemporaryFile("w", suffix=".fasta", delete=False) as tf:
        kept = []
        for srr in not_ass:
            if not all((FASTA_DIR_H5N1 / f"{srr}_{sg}_cns.fa").exists() for sg in segs_only):
                continue
            parts = []
            ok = True
            for sg, a, b in bounds:
                with open(FASTA_DIR_H5N1 / f"{srr}_{sg}_cns.fa") as sf:
                    seq = "".join(l.strip() for l in sf if not l.startswith(">") and l.strip()).upper()
                if len(seq) != (b - a) or seq.count("N") > 0.05 * len(seq):
                    ok = False
                    break
                parts.append(seq)
            if ok:
                tf.write(f">{srr}\n" + "".join(parts) + "\n")
                kept.append(srr)
        tmp_fa = tf.name

    try:
        out_csv = subprocess.check_output(
            [EGGMAN_BIN, "--pack", str(ROOT / "rhiztype" / "packs" / "influenza_a_8seg"), "-o", "/dev/stdout", "-q", tmp_fa],
            text=True,
        )
    finally:
        os.remove(tmp_fa)

    for r in csv.DictReader(io.StringIO(out_csv)):
        if r["called_lineage"] == "D1.8" and len(d18_candidates) < 10:
            d18_candidates.append(r["sequence_id"])
        elif (
            r["called_lineage"] == "URF"
            and r["structural_formula"] == "EA3[NAm_B313:PB2 | NAm_D1:PA+NP+NA]"
            and len(novel_candidates) < 5
        ):
            novel_candidates.append(r["sequence_id"])

    print(f"[+] Selected {len(d18_candidates)} USDA-unassigned D1.8 isolates and {len(novel_candidates)} novel EA3[NAm_B313:PB2 | NAm_D1:PA+NP+NA] isolates.")

    cohorts = [
        (
            "D1.8_Unassigned",
            "NAm_D1[EA3:PB1+HA+MP+NS | NAm_B2:NP]",
            d18_candidates,
            {
                "PB2": "NAm_D1",
                "PB1": "EA3",
                "PA": "NAm_D1",
                "HA": "EA3",
                "NP": "NAm_B2",
                "NA": "NAm_D1",
                "MP": "EA3",
                "NS": "EA3",
            },
        ),
        (
            "Novel_URF_EA3_B313_D1",
            "EA3[NAm_B313:PB2 | NAm_D1:PA+NP+NA]",
            novel_candidates,
            {
                "PB2": "NAm_B313",
                "PB1": "EA3",
                "PA": "NAm_D1",
                "HA": "EA3",
                "NP": "NAm_D1",
                "NA": "NAm_D1",
                "MP": "EA3",
                "NS": "EA3",
            },
        ),
    ]

    h5_tree_dir = RES / "trees_h5n1"
    h5_tree_dir.mkdir(parents=True, exist_ok=True)
    for stale in h5_tree_dir.glob("*.ckp.gz"):
        tf = Path(str(stale).replace(".ckp.gz", ".treefile"))
        if not tf.exists():
            stale.unlink(missing_ok=True)

    results_rows = []

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

    # Load official USDA GenoFLU reference allele sequences & match queries via blastn
    CACHE_GENOFLU.mkdir(parents=True, exist_ok=True)
    geno_ref_fa = CACHE_GENOFLU / "all_genoflu_refs.fasta"
    if not geno_ref_fa.exists() or geno_ref_fa.stat().st_size < 10000:
        fastas_meta = json.loads(urllib.request.urlopen("https://api.github.com/repos/USDA-VS/GenoFLU/contents/dependencies/fastas", timeout=30).read().decode("utf-8"))
        with open(geno_ref_fa, "w") as out_gf:
            for fm in fastas_meta:
                txt = urllib.request.urlopen(fm["download_url"], timeout=30).read().decode("utf-8", errors="replace")
                out_gf.write(txt.rstrip() + "\n")
    genoflu_allele_seqs = defaultdict(dict)  # gene -> {allele: seq}
    cur_al, cur_gene, buf_gf = None, None, []
    for ln in geno_ref_fa.read_text().splitlines():
        ln = ln.strip()
        if not ln:
            continue
        if ln.startswith(">"):
            if cur_al and cur_gene and buf_gf:
                s_gf = "".join(buf_gf).upper()
                if cur_al not in genoflu_allele_seqs[cur_gene] or len(s_gf) > len(genoflu_allele_seqs[cur_gene][cur_al]):
                    genoflu_allele_seqs[cur_gene][cur_al] = s_gf
            parts = ln[1:].split()
            cur_al = parts[0] if len(parts) >= 3 else None
            cur_gene = parts[2] if len(parts) >= 3 else None
            buf_gf = []
        else:
            buf_gf.append(ln)
    if cur_al and cur_gene and buf_gf:
        s_gf = "".join(buf_gf).upper()
        if cur_al not in genoflu_allele_seqs[cur_gene] or len(s_gf) > len(genoflu_allele_seqs[cur_gene][cur_al]):
            genoflu_allele_seqs[cur_gene][cur_al] = s_gf

    # Run blastn to identify official GenoFLU allele nearest to each query segment
    top_geno_allele_by_qseg = {}
    with tempfile.TemporaryDirectory() as td_gf:
        db_p = os.path.join(td_gf, "gf_db")
        subprocess.run([MAKEBLASTDB_BIN, "-in", str(geno_ref_fa), "-dbtype", "nucl", "-out", db_p], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        q_fa_gf = os.path.join(td_gf, "q16.fasta")
        with open(q_fa_gf, "w") as qf:
            for cname, _, srrs, _ in cohorts:
                for sg, _, _ in bounds:
                    for srr in srrs:
                        with open(FASTA_DIR_H5N1 / f"{srr}_{sg}_cns.fa") as sf:
                            qseq = "".join(l.strip() for l in sf if not l.startswith(">") and l.strip()).upper()
                        qf.write(f">{cname}__{srr}__{sg}\n{qseq}\n")
        bout_gf = os.path.join(td_gf, "blast.tsv")
        subprocess.run([BLASTN_BIN, "-query", q_fa_gf, "-db", db_p, "-word_size", "11", "-max_target_seqs", "5", "-num_threads", "4", "-outfmt", "6 qseqid bitscore pident stitle", "-out", bout_gf], check=True, stderr=subprocess.DEVNULL)
        hits_gf = defaultdict(list)
        for bl in open(bout_gf):
            c_b = bl.rstrip("\n").split("\t")
            al_b, _, g_b = c_b[3].split()[:3]
            hits_gf[c_b[0]].append((float(c_b[1]), float(c_b[2]), al_b, g_b))
        for qk, lst_h in hits_gf.items():
            lst_h.sort(key=lambda x: (-x[0], -x[1]))
            top_geno_allele_by_qseg[qk] = (lst_h[0][2], lst_h[0][1])

    core_gf_alleles = {"ea1", "ea2", "ea3", "ea4", "am1.1", "am1.2", "am2.1", "am2.2", "am3", "am4", "am5", "am11", "am12", "am19", "am20"}

    tasks = []
    for cname, formula, srrs, exp_map in cohorts:
        for sg, a, b in bounds:
            raw_fa = h5_tree_dir / f"{cname}_{sg}.raw.fasta"
            fa_path = h5_tree_dir / f"{cname}_{sg}.aln.fasta"
            tree_path = Path(str(fa_path) + ".treefile")

            valid_refs = [rk for rk, rseq in profiles.items() if rseq[a:b].count("N") < 0.2 * (b - a)]
            matched_alleles = {top_geno_allele_by_qseg[f"{cname}__{srr}__{sg}"][0] for srr in srrs if f"{cname}__{srr}__{sg}" in top_geno_allele_by_qseg}
            gf_alleles_use = sorted((core_gf_alleles | matched_alleles) & set(genoflu_allele_seqs.get(sg, {}).keys()))
            q_names = []
            with open(raw_fa, "w") as f:
                for rk in valid_refs:
                    f.write(f">Ref_{rk}\n{profiles[rk][a:b]}\n")
                for al in gf_alleles_use:
                    f.write(f">Ref_GenoFLU_{al}\n{genoflu_allele_seqs[sg][al]}\n")
                for srr in srrs:
                    with open(FASTA_DIR_H5N1 / f"{srr}_{sg}_cns.fa") as sf:
                        qseq = "".join(l.strip() for l in sf if not l.startswith(">") and l.strip()).upper()
                    qn = f"Query_{srr}"
                    q_names.append(qn)
                    f.write(f">{qn}\n{qseq}\n")
            tasks.append((cname, formula, sg, exp_map[sg], valid_refs, q_names, raw_fa, fa_path, tree_path, sorted(matched_alleles)))

    def _run_tree(item):
        raw_p, aln_p, tr_p = item[6], item[7], item[8]
        if not aln_p.exists():
            res = subprocess.run([MAFFT_BIN, "--auto", "--quiet", str(raw_p)], check=True, stdout=subprocess.PIPE, text=True)
            aln_p.write_text(res.stdout)
        if not tr_p.exists():
            subprocess.run(
                [IQTREE_BIN, "-s", str(aln_p), "-m", "GTR+G4", "-B", "1000", "-T", "1", "--quiet", "-redo"],
                check=True,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
            )

    with concurrent.futures.ThreadPoolExecutor(max_workers=8) as pool:
        list(pool.map(_run_tree, tasks))

    for cname, formula, sg, exp_donor, valid_refs, q_names, _, fa_path, tree_path, matched_gf in tasks:
        tree = Phylo.read(str(tree_path), "newick")
        tree.root_at_midpoint()
        q_set = set(q_names)

        # Find smallest clade containing all query sequences + strictly one Ref_<pack_profile> donor (may also include Ref_GenoFLU_*)
        pure_pack_clades = []
        pure_gf_clades = []
        for clade in tree.find_clades(order="preorder"):
            terms = [t.name for t in clade.get_terminals()]
            q_in = [t for t in terms if t in q_set]
            if len(q_in) < len(q_set):
                continue
            pack_refs = [t[4:] for t in terms if t.startswith("Ref_") and not t.startswith("Ref_GenoFLU_")]
            gf_refs = [t[12:] for t in terms if t.startswith("Ref_GenoFLU_")]
            ufb = parse_ufboot(clade)
            if len(pack_refs) == 1:
                pure_pack_clades.append((len(terms), -(ufb if pd.notna(ufb) else 0.0), pack_refs[0], ufb))
            if len(gf_refs) == 1 and len(pack_refs) <= 1:
                pure_gf_clades.append((len(terms), -(ufb if pd.notna(ufb) else 0.0), gf_refs[0], ufb))

        if pure_pack_clades:
            pure_pack_clades.sort()
            _, _, sister_donor, ufb_val = pure_pack_clades[0]
        else:
            ref_terms = [t for t in tree.get_terminals() if t.name.startswith("Ref_") and not t.name.startswith("Ref_GenoFLU_")]
            q_terms = [t for t in tree.get_terminals() if t.name in q_set]
            best_ref = min(ref_terms, key=lambda rt: np.mean([tree.distance(rt, qt) for qt in q_terms]))
            sister_donor = best_ref.name[4:]
            ufb_val = np.nan

        if pure_gf_clades:
            pure_gf_clades.sort()
            _, _, sister_gf, gf_ufb = pure_gf_clades[0]
        else:
            gf_terms = [t for t in tree.get_terminals() if t.name.startswith("Ref_GenoFLU_")]
            q_terms = [t for t in tree.get_terminals() if t.name in q_set]
            best_gf = min(gf_terms, key=lambda rt: np.mean([tree.distance(rt, qt) for qt in q_terms]))
            sister_gf = best_gf.name[12:]
            gf_ufb = np.nan

        best_q_size = 999999
        q_ufb = np.nan
        for clade in tree.find_clades(order="preorder"):
            terms = [t.name for t in clade.get_terminals()]
            if sum(1 for t in terms if t in q_set) == len(q_set) and len(terms) < best_q_size:
                best_q_size = len(terms)
                q_ufb = parse_ufboot(clade)

        concordant_exact = sister_donor == exp_donor
        results_rows.append(
            {
                "cohort": cname,
                "eggman_formula": formula,
                "segment": sg,
                "n_queries": len(q_set),
                "n_pack_refs": len(valid_refs),
                "expected_donor": exp_donor,
                "ml_sister_donor": sister_donor,
                "exact_donor_concordant": concordant_exact,
                "sister_clade_ufboot_pct": round(float(ufb_val), 1) if pd.notna(ufb_val) else np.nan,
                "blast_top_genoflu_alleles": ";".join(matched_gf),
                "ml_sister_genoflu_allele": sister_gf,
                "genoflu_sister_ufboot_pct": round(float(gf_ufb), 1) if pd.notna(gf_ufb) else np.nan,
                "query_monophyly_ufboot_pct": round(float(q_ufb), 1) if pd.notna(q_ufb) else np.nan,
                "strictly_monophyletic_queries": best_q_size == len(q_set),
                "treefile": str(tree_path.name),
            }
        )
        print(
            f"  [{cname:21s} | {sg:3s}] Expected={exp_donor:8s}, ML_Pack={sister_donor:8s} "
            f"(Match={str(concordant_exact):5s}, UFB={ufb_val}%) | ML_GenoFLU={sister_gf:6s} (UFB={gf_ufb}%)",
            flush=True,
        )

    out_path = RES / "usda_h5n1_unassigned_iqtree_validation.csv"
    pd.DataFrame(results_rows).to_csv(out_path, index=False)
    print(f"[+] Saved 16-segment USDA H5N1 ML tree validation -> {out_path}")


def main():
    est_segs, urf_segs, man_df = load_cohorts()
    run_stage1_controls(est_segs, urf_segs, man_df)
    run_stage2_controls(est_segs, urf_segs)
    run_usda_h5n1_iqtree_validation()


if __name__ == "__main__":
    main()
