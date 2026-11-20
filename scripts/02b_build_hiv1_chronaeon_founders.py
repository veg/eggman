#!/usr/bin/env python3
"""
02b_build_hiv1_chronaeon_founders.py
====================================
Rebuild HIV-1 Group M's 13 scrubbed basal subtype founders (+ `?E` on HXB2 5097..9417)
in `rhiztype/crates/eggman-core/data/profiles_9719.tsv` and
`data/hxb2_15trunk_profile_9719.fasta` from dated non-mosaic LANL Group M sequences
(`lanl_ref_2023.fasta` + `lanl_flt_2022.fasta`) via deletion-aware ChronAeon
exponential temporal decay (t_{1/2} = 6.93 yr, w_s = exp(-ln(2) * (t_s - t_min) / 6.93)).

Scrubbing rules (`C1d-A`):
  1. Exclude A3, A4, A7, A8 (only A1 -> A1, A2 -> A2, A6 -> A6).
  2. Exclude any LANL nominal pure sequence flagged as non-clonal / cryptic recombinant
     in `results/lanl_all_full_genomes_rhizaeon_validation.csv` (except Subtype G, which
     feeds G0 after masking the ancestral A1 gp41..nef tract).
  3. Mask all non-ghost contaminated tracts identified in Phase 1
     (`results/phase1_contaminated_pure_references.csv`: PT2695 6627..8880,
     93BR020_1 7074..7651, 95IN21068 7233..7438, AF286236 8976..9417).
  4. Mask all Subtype G sequences across HXB2 8262..8880 (ancestral A1 gp41..nef tract)
     and emit '?' in G0 across 8262..8880.
  5. Deletion-aware ChronAeon aggregation at each HXB2 column j (1..9719):
     - Leading/trailing missing sequence ends outside [first_valid, last_valid] have weight 0.
     - Masked/ambiguous bases (N/?) have weight 0.
     - Internal gaps ('-'/'.') count toward W_gap(j); valid A/C/G/T count toward W_base(j).
     - If W_gap(j) / (W_gap(j) + W_base(j)) > 0.50 -> emit '-'.
     - Otherwise emit argmax_{b in {A,C,G,T}} W_b(j).
  6. Keep `HXB2` (K03455, 9,719 nt) and `CM240` (U54771, 9,719 nt), and emit `?E`
     from dated CRF01_AE sequences via deletion-aware ChronAeon on 5097..9417 ('-' outside).
"""

import json
import math
from collections import defaultdict
from pathlib import Path
import numpy as np
import pandas as pd

PROJECT_ROOT = Path("/Users/sergei/Projects/HIV-recombination")
DATA_DIR = PROJECT_ROOT / "data"
RESULTS_DIR = PROJECT_ROOT / "results"
CORE_DATA_DIR = PROJECT_ROOT / "rhiztype" / "crates" / "eggman-core" / "data"

HALF_LIFE_YR = 6.93
GAMMA_DECAY = math.log(2.0) / HALF_LIFE_YR  # ~0.1000212 / yr

BASE_TO_INT = {ord("A"): 0, ord("C"): 1, ord("G"): 2, ord("T"): 3, ord("U"): 3}
INT_TO_BASE = {0: "A", 1: "C", 2: "G", 3: "T"}

PURE_TRUNKS = ["A1", "A2", "A6", "B", "C", "D", "F1", "F2", "G0", "H", "J", "K", "L"]


def read_fasta(path: Path):
    records = []
    hdr = None
    chunks = []
    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if hdr is not None:
                    records.append((hdr, "".join(chunks).upper()))
                hdr = line[1:].strip()
                chunks = []
            else:
                chunks.append(line)
        if hdr is not None:
            records.append((hdr, "".join(chunks).upper()))
    return records


def project_alignment_to_hxb2_9719(fasta_path: Path):
    """Map every MSA record onto exact 9,719 HXB2 non-gap columns."""
    recs = read_fasta(fasta_path)
    hxb2_seqs = [seq for hdr, seq in recs if "K03455" in hdr]
    if len(hxb2_seqs) != 1:
        raise RuntimeError(f"Expected 1 HXB2 K03455 record in {fasta_path}, got {len(hxb2_seqs)}")
    cols_9719 = np.array([i for i, ch in enumerate(hxb2_seqs[0]) if ch not in "-."], dtype=np.int32)
    if len(cols_9719) != 9719:
        raise RuntimeError(f"Expected 9,719 HXB2 columns in {fasta_path}, got {len(cols_9719)}")
    projected = []
    for hdr, seq in recs:
        b = np.frombuffer(seq.encode("ascii", errors="replace"), dtype=np.uint8)[cols_9719]
        projected.append((hdr, b))
    return projected


def deletion_aware_chronaeon_consensus(byte_rows, years, mask_ranges_per_seq, default_char="-"):
    """
    Aggregate aligned 9,719-nt byte rows via deletion-aware ChronAeon exponential decay.
    - Leading/trailing '-' outside each sequence's sequenced [first_valid..last_valid] span
      and masked/ambiguous positions get weight 0.
    - Internal '-' within [first_valid..last_valid] contributes to gap weight W_gap.
    """
    n_seqs = len(byte_rows)
    if n_seqs == 0:
        raise RuntimeError("No sequences supplied to deletion_aware_chronaeon_consensus")
    t_min = min(years)
    weights = np.array([math.exp(-GAMMA_DECAY * (t - t_min)) for t in years], dtype=np.float64)

    # Encode states: 0..3 = A,C,G,T ; -1 = internal deletion gap ; -2 = missing/masked (ignored)
    mat = np.full((n_seqs, 9719), -2, dtype=np.int8)
    for row_i, (arr, masks) in enumerate(zip(byte_rows, mask_ranges_per_seq)):
        is_base = (
            (arr == ord("A"))
            | (arr == ord("C"))
            | (arr == ord("G"))
            | (arr == ord("T"))
            | (arr == ord("U"))
        )
        valid_pos = np.where(is_base)[0]
        if len(valid_pos) == 0:
            continue
        p0, p1 = int(valid_pos[0]), int(valid_pos[-1])
        sub = arr[p0 : p1 + 1]
        row_enc = np.full(p1 - p0 + 1, -2, dtype=np.int8)
        for ch_ord, b_idx in BASE_TO_INT.items():
            row_enc[sub == ch_ord] = b_idx
        row_enc[(sub == ord("-")) | (sub == ord("."))] = -1
        mat[row_i, p0 : p1 + 1] = row_enc
        for h_start, h_end in masks:
            mat[row_i, (h_start - 1) : h_end] = -2

    w_gap = weights @ (mat == -1).astype(np.float64)
    w_bases = np.vstack([weights @ (mat == b).astype(np.float64) for b in range(4)])
    w_base_tot = np.sum(w_bases, axis=0)
    w_tot = w_gap + w_base_tot
    best_b = np.argmax(w_bases, axis=0)

    out = []
    for j in range(9719):
        if w_tot[j] <= 0.0:
            out.append(default_char)
        elif (w_gap[j] / w_tot[j]) > 0.50:
            out.append("-")
        else:
            out.append(INT_TO_BASE[int(best_b[j])])
    return "".join(out), t_min, max(years), len(years)


def main():
    print("=" * 92)
    print("STEP 2B: REBUILDING SCRUBBED HIV-1 CHRONAEON FOUNDERS (t_1/2 = 6.93 yr, DELETION-AWARE)")
    print("=" * 92)

    ref_meta = pd.read_csv(DATA_DIR / "lanl_ref_2023_metadata.csv")
    flt_meta = pd.read_csv(DATA_DIR / "lanl_flt_2022_metadata.csv")
    ref_year_map = dict(zip(ref_meta["seq_id"], ref_meta["year"]))
    ref_acc_map = dict(zip(ref_meta["seq_id"], ref_meta["accession"]))
    flt_year_map = dict(zip(flt_meta["seq_id"], flt_meta["year"]))
    flt_acc_map = dict(zip(flt_meta["seq_id"], flt_meta["accession"]))

    # Identify accessions flagged as non-clonal / cryptic mosaic in LANL validation (except Subtype G's ancestral A1 tract)
    bad_accs = set()
    val_path = RESULTS_DIR / "lanl_all_full_genomes_rhizaeon_validation.csv"
    if val_path.exists():
        val_df = pd.read_csv(val_path)
        pure_nom = {"A1", "A2", "A6", "B", "C", "D", "F1", "F2", "H", "J", "K", "L"}
        bad = val_df[
            val_df["lanl_assigned_tag"].isin(pure_nom)
            & (val_df["validation_verdict"] != "Concordant_Clonal_Subtype")
        ]
        bad_accs = set(bad["sequence_id"].apply(lambda s: str(s).split(".")[-1]))
        print(f"[+] Excluding {len(bad_accs)} cryptic non-clonal LANL nominal-pure accessions from founder pools")

    # Contaminated intervals from Phase 1 audit
    contam_by_acc = defaultdict(list)
    contam_df = pd.read_csv(RESULTS_DIR / "phase1_contaminated_pure_references.csv")
    for _, row in contam_df.iterrows():
        if bool(row["is_ghost"]):
            continue
        acc = str(row["recombinant_id"]).split(".")[-1]
        contam_by_acc[acc].append((int(row["hxb2_start"]), int(row["hxb2_end"])))

    ref_proj = project_alignment_to_hxb2_9719(DATA_DIR / "lanl_ref_2023.fasta")
    flt_proj = project_alignment_to_hxb2_9719(DATA_DIR / "lanl_flt_2022.fasta")

    hxb2_str = None
    cm240_str = None
    pools = defaultdict(list)  # tag -> list of (acc, seq_id, yr, byte_arr, masks)
    seen_acc_by_tag = defaultdict(set)

    def ingest(records, year_map, acc_map, source_name):
        nonlocal hxb2_str, cm240_str
        for hdr, arr in records:
            if "K03455" in hdr and hxb2_str is None:
                hxb2_str = arr.tobytes().decode("ascii")
                continue
            if "U54771" in hdr and cm240_str is None:
                cm240_str = arr.tobytes().decode("ascii")
            acc = str(acc_map.get(hdr, hdr.split(".")[-1]))
            yr = year_map.get(hdr)
            if pd.isna(yr) or yr is None:
                continue
            yr = float(yr)
            parts = hdr.split(".")
            tag = parts[1] if parts[0] == "Ref" else parts[0]
            if tag not in {"A1", "A2", "A6", "B", "C", "D", "F1", "F2", "G", "H", "J", "K", "L", "01_AE"}:
                continue
            if acc in bad_accs:
                continue
            if acc in seen_acc_by_tag[tag]:
                continue
            seen_acc_by_tag[tag].add(acc)
            masks = list(contam_by_acc.get(acc, []))
            if tag == "G":
                masks.append((8262, 8880))
            pools[tag].append((acc, hdr, yr, arr, masks, source_name))

    ingest(ref_proj, ref_year_map, ref_acc_map, "LANL_Ref2023")
    ingest(flt_proj, flt_year_map, flt_acc_map, "LANL_Flt2022")
    assert hxb2_str is not None and len(hxb2_str) == 9719
    assert cm240_str is not None and len(cm240_str) == 9719

    # Read previous profiles_9719.tsv for exact SNP diff auditing
    old_profiles = {}
    old_tsv = CORE_DATA_DIR / "profiles_9719.tsv"
    if old_tsv.exists():
        with open(old_tsv) as f:
            for line in f:
                cols = line.rstrip("\n").split("\t")
                if len(cols) >= 2:
                    old_profiles[cols[0]] = cols[1].strip()

    new_profiles = {"HXB2": hxb2_str}
    audit_rows = []

    for trunk in PURE_TRUNKS:
        src_tag = "G" if trunk == "G0" else trunk
        entries = pools[src_tag]
        arrs = [e[3] for e in entries]
        yrs = [e[2] for e in entries]
        masks = [e[4] for e in entries]
        cons, t_min, t_max, n_seqs = deletion_aware_chronaeon_consensus(
            arrs, yrs, masks, default_char=("?" if trunk == "G0" else "-")
        )
        if trunk == "G0":
            cons = cons[:8261] + ("?" * (8880 - 8261)) + cons[8880:]
        new_profiles[trunk] = cons

        old_s = old_profiles.get(trunk, "")
        snp_diffs_cod = 0
        indel_diffs_cod = 0
        if len(old_s) == 9719:
            for i in range(789, 9417):
                a, b = old_s[i], cons[i]
                if a in "ACGT" and b in "ACGT" and a != b:
                    snp_diffs_cod += 1
                elif (a in "ACGT") != (b in "ACGT"):
                    indel_diffs_cod += 1
        acgt_cod = sum(cons[789:9417].count(x) for x in "ACGT")
        gap_cod = cons[789:9417].count("-")
        q_cod = cons[789:9417].count("?")
        audit_rows.append({
            "trunk": trunk,
            "source_tag": src_tag,
            "n_dated_seqs": n_seqs,
            "t_min": t_min,
            "t_max": t_max,
            "coding_acgt_nt": acgt_cod,
            "coding_gap_nt": gap_cod,
            "coding_masked_q_nt": q_cod,
            "coding_snp_diffs_vs_old": snp_diffs_cod,
            "coding_indel_diffs_vs_old": indel_diffs_cod,
        })
        print(
            f"  {trunk:4s} (N={n_seqs:4d}, {t_min:.0f}-{t_max:.0f}) | coding ACGT={acgt_cod:4d} gap={gap_cod:3d} ?={q_cod:3d} | "
            f"diffs vs old={snp_diffs_cod:3d} SNPs, {indel_diffs_cod:3d} indels"
        )

    # Build ?E on 5097..9417 via deletion-aware ChronAeon from dated 01_AE sequences ('-' outside 5097..9417)
    ae_entries = pools["01_AE"]
    ae_cons, ae_tmin, ae_tmax, ae_n = deletion_aware_chronaeon_consensus(
        [e[3] for e in ae_entries],
        [e[2] for e in ae_entries],
        [e[4] for e in ae_entries],
        default_char="-",
    )
    ghost_e = ("-" * 5096) + ae_cons[5096:9417] + ("-" * (9719 - 9417))
    assert len(ghost_e) == 9719
    new_profiles["?E"] = ghost_e
    new_profiles["CM240"] = cm240_str

    old_e = old_profiles.get("?E", "")
    e_snp = sum(1 for i in range(5096, 9417) if i < len(old_e) and old_e[i] in "ACGT" and ghost_e[i] in "ACGT" and old_e[i] != ghost_e[i])
    e_ind = sum(1 for i in range(5096, 9417) if i < len(old_e) and (old_e[i] in "ACGT") != (ghost_e[i] in "ACGT"))
    audit_rows.append({
        "trunk": "?E",
        "source_tag": "01_AE(5097..9417)",
        "n_dated_seqs": ae_n,
        "t_min": ae_tmin,
        "t_max": ae_tmax,
        "coding_acgt_nt": sum(ghost_e[5096:9417].count(x) for x in "ACGT"),
        "coding_gap_nt": ghost_e[5096:9417].count("-"),
        "coding_masked_q_nt": 0,
        "coding_snp_diffs_vs_old": e_snp,
        "coding_indel_diffs_vs_old": e_ind,
    })
    print(f"  ?E   (N={ae_n:4d}, {ae_tmin:.0f}-{ae_tmax:.0f}) | 5097..9417 ACGT={sum(ghost_e[5096:9417].count(x) for x in 'ACGT'):4d} gap={ghost_e[5096:9417].count('-'):3d} | diffs vs old={e_snp:3d} SNPs, {e_ind:3d} indels")

    # Write profiles_9719.tsv
    row_order = ["HXB2"] + PURE_TRUNKS + ["?E", "CM240"]
    with open(old_tsv, "w") as f:
        for k in row_order:
            f.write(f"{k}\t{new_profiles[k]}\n")
    print(f"[+] Wrote {len(row_order)} scrubbed ChronAeon profiles (9,719 nt) to {old_tsv}")

    # Write data/hxb2_15trunk_profile_9719.fasta (matching expected MAFFT profile headers)
    prof_fa = DATA_DIR / "hxb2_15trunk_profile_9719.fasta"
    with open(prof_fa, "w") as f:
        f.write(f">Profile_HXB2_K03455\n{new_profiles['HXB2']}\n")
        for k in PURE_TRUNKS:
            f.write(f">Profile_Trunk_{k}\n{new_profiles[k].replace('?', '-')}\n")
        f.write(f">Profile_Trunk_E\n{new_profiles['?E']}\n")
    print(f"[+] Wrote 15-trunk MAFFT anchor profile FASTA to {prof_fa}")

    audit_df = pd.DataFrame(audit_rows)
    audit_path = RESULTS_DIR / "hiv1_chronaeon_founder_audit.tsv"
    audit_df.to_csv(audit_path, sep="\t", index=False)
    print(f"[+] Wrote founder audit table to {audit_path}")


if __name__ == "__main__":
    main()
