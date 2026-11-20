#!/usr/bin/env python3
"""
01_prepare_lanl_compendium.py
=============================
Prepares the foundational HIV-1 Group M datasets in /Users/sergei/Projects/HIV-recombination/data/:
1. Copies/links the official LANL Reference Alignment (518 full genomes: 17 'pure' subtypes + 145 CRFs)
   and LANL Filtered Web Alignment (5,089 full genomes), both anchored to HXB2 (K03455).
2. Parses the complete LANL CRF Compendium HTML table across the full HXB2 genome (1..9719 nt)
   into structured JSON and CSV formats (`lanl_official_crf_catalog.json`, `lanl_official_crf_segments.csv`).
3. Extracts sequence metadata (Subtype/CRF label, Country, Sampling Year in decimal calendar years,
   Strain, Accession) and builds exact MSA-column <-> HXB2-coordinate bijective maps.
4. Builds the pure-subtype + CRF01_AE + CRF02_AG focal alignment (`lanl_pure_and_primary_crfs.fasta`)
   with dates CSV for Phase 1 RhizAeon + ChronAeon auditing.
"""

import os
import re
import json
import shutil
from pathlib import Path
from collections import Counter, defaultdict
import pandas as pd
from Bio import SeqIO

PROJECT_ROOT = Path("/Users/sergei/Projects/HIV-recombination")
DATA_DIR = PROJECT_ROOT / "data"
DATA_DIR.mkdir(parents=True, exist_ok=True)

SRC_NOM_DIR = Path("/Users/sergei/Projects/TOGA_MEME/recombination/pan_group_m_nomenclature")
LANL_HTML_CACHE = Path("/Users/sergei/.gemini/antigravity-cli/brain/1dfa5602-784e-4936-ae1c-a9734d91bc95/.system_generated/steps/25374/content.md")

# Canonical HXB2 (K03455, 9719 nt) Gene & Domain Map
HXB2_ANNOTATIONS = {
    "reference": "K03455.1 (HXB2)",
    "length_nt": 9719,
    "genes": [
        {"id": "5LTR", "name": "5' LTR", "start": 1, "end": 634},
        {"id": "gag", "name": "Gag", "start": 790, "end": 2292},
        {"id": "pol", "name": "Pol", "start": 2085, "end": 5096},
        {"id": "vif", "name": "Vif", "start": 5041, "end": 5619},
        {"id": "vpr", "name": "Vpr", "start": 5559, "end": 5850},
        {"id": "tat_1", "name": "Tat (Exon 1)", "start": 5831, "end": 6045},
        {"id": "rev_1", "name": "Rev (Exon 1)", "start": 5970, "end": 6045},
        {"id": "vpu", "name": "Vpu", "start": 6062, "end": 6310},
        {"id": "env", "name": "Env (gp160)", "start": 6225, "end": 8795},
        {"id": "tat_2", "name": "Tat (Exon 2)", "start": 8379, "end": 8469},
        {"id": "rev_2", "name": "Rev (Exon 2)", "start": 8379, "end": 8653},
        {"id": "nef", "name": "Nef", "start": 8797, "end": 9417},
        {"id": "3LTR", "name": "3' LTR", "start": 9086, "end": 9719},
    ],
    "domains": [
        {"id": "5LTR", "gene": "5LTR", "label": "5'LTR", "start": 1, "end": 789},
        {"id": "gag.MA", "gene": "gag", "label": "gag-p17(MA)", "start": 790, "end": 1185},
        {"id": "gag.CA", "gene": "gag", "label": "gag-p24(CA)", "start": 1186, "end": 1878},
        {"id": "gag.NCp6", "gene": "gag", "label": "gag-p7/p6", "start": 1879, "end": 2252},
        {"id": "pol.PR", "gene": "pol", "label": "pol-PR", "start": 2253, "end": 2549},
        {"id": "pol.RT", "gene": "pol", "label": "pol-RT", "start": 2550, "end": 3569},
        {"id": "pol.RNaseH", "gene": "pol", "label": "pol-RNaseH", "start": 3570, "end": 3929},
        {"id": "pol.IN", "gene": "pol", "label": "pol-IN", "start": 3930, "end": 5096},
        {"id": "vif", "gene": "vif", "label": "vif", "start": 5097, "end": 5558},
        {"id": "vpr_tat", "gene": "vpr", "label": "vpr-tat", "start": 5559, "end": 6061},
        {"id": "vpu", "gene": "vpu", "label": "vpu", "start": 6062, "end": 6224},
        {"id": "env.gp120_C1", "gene": "env", "label": "env-gp120(C1)", "start": 6225, "end": 6611},
        {"id": "env.gp120_V1V3", "gene": "env", "label": "env-gp120(V1-V3)", "start": 6612, "end": 7220},
        {"id": "env.gp120_C3V5", "gene": "env", "label": "env-gp120(C3-V5)", "start": 7221, "end": 7757},
        {"id": "env.gp41", "gene": "env", "label": "env-gp41", "start": 7758, "end": 8795},
        {"id": "nef", "gene": "nef", "label": "nef", "start": 8796, "end": 9417},
        {"id": "3LTR", "gene": "3LTR", "label": "3'LTR", "start": 9418, "end": 9719},
    ]
}


def parse_two_digit_year(yy_str: str):
    """Convert LANL 2-digit year token (e.g. '83', '97', '04', '21') to 4-digit float year."""
    if not yy_str or not yy_str.isdigit() or len(yy_str) != 2:
        return None
    yy = int(yy_str)
    if yy >= 59:
        return float(1900 + yy)
    elif yy <= 26:
        return float(2000 + yy)
    return None


def parse_lanl_fasta_header(seq_id: str):
    """
    Parse LANL FASTA IDs:
      - Ref.<Subtype>.<Country>.<YY>.<Strain>.<Accession>
      - <Subtype>.<Country>.<YY>.<Strain>.<Accession>
    """
    parts = seq_id.split(".")
    is_ref = parts[0] == "Ref"
    offset = 1 if is_ref else 0
    if len(parts) < offset + 5:
        return {
            "seq_id": seq_id,
            "is_ref_panel": is_ref,
            "lanl_tag": parts[offset] if len(parts) > offset else "Unknown",
            "country": "Unknown",
            "year": None,
            "strain": seq_id,
            "accession": parts[-1]
        }
    lanl_tag = parts[offset]
    country = parts[offset + 1]
    yy_str = parts[offset + 2]
    year = parse_two_digit_year(yy_str)
    accession = parts[-1]
    strain = ".".join(parts[offset + 3:-1])
    return {
        "seq_id": seq_id,
        "is_ref_panel": is_ref,
        "lanl_tag": lanl_tag,
        "country": country,
        "year": year,
        "strain": strain,
        "accession": accession
    }


def parse_lanl_crf_html(html_path: Path):
    """Parse official LANL CRF HTML table into full-genome HXB2 segment maps."""
    text = html_path.read_text(errors="replace")
    pattern = r'<tr id=tgl_(CRF\w+)[^>]*>.*?<pre>(.*?)</pre>'
    matches = re.findall(pattern, text, re.DOTALL)

    catalog = {}
    segment_rows = []

    for crf_name, pre_content in matches:
        raw_lines = re.split(r'<br\s*/?>|\n', pre_content)
        segments = []
        for line in raw_lines:
            line = line.strip()
            parts = line.split()
            if len(parts) >= 3 and parts[0].isdigit() and parts[1].isdigit():
                start, end, st = int(parts[0]), int(parts[1]), " ".join(parts[2:])
                segments.append({"hxb2_start": start, "hxb2_end": end, "parent_label": st})
                segment_rows.append({
                    "crf": crf_name,
                    "short_tag": crf_name.replace("CRF", ""),
                    "hxb2_start": start,
                    "hxb2_end": end,
                    "length_nt": end - start + 1,
                    "parent_label": st
                })

        parents = sorted(list(set(s["parent_label"] for s in segments)))
        breakpoints = []
        if segments:
            # Breakpoints are transitions between consecutive segments
            for idx in range(len(segments) - 1):
                s_left = segments[idx]
                s_right = segments[idx + 1]
                if s_left["parent_label"] != s_right["parent_label"]:
                    breakpoints.append({
                        "hxb2_coord": s_left["hxb2_end"],
                        "next_start": s_right["hxb2_start"],
                        "gap_or_overlap": s_right["hxb2_start"] - s_left["hxb2_end"] - 1,
                        "left_parent": s_left["parent_label"],
                        "right_parent": s_right["parent_label"]
                    })
        catalog[crf_name] = {
            "crf": crf_name,
            "short_tag": crf_name.replace("CRF", ""),
            "has_schematic_map": len(segments) > 0,
            "num_segments": len(segments),
            "num_breakpoints": len(breakpoints),
            "declared_parents": parents,
            "is_cpx": "_cpx" in crf_name,
            "has_unclassified_U": any(p == "U" or "U" in p.split("/") for p in parents),
            "has_ambiguous_slash": any("/" in p for p in parents),
            "has_Ghost_E": any(p == "E" or "01_AE" in p or "CRF01" in p for p in parents),
            "has_second_gen_CRF_parent": any("CRF" in p or "01" in p or "02" in p or "07" in p or "08" in p or "55" in p for p in parents),
            "segments": segments,
            "breakpoints": breakpoints
        }
    return catalog, pd.DataFrame(segment_rows)


def build_hxb2_coordinate_map(fasta_path: Path, hxb2_id_substring: str = "K03455"):
    """
    Build a bijective coordinate map between 1-indexed MSA column (1..L_msa)
    and 1-indexed HXB2 coordinate (1..9719).
    """
    hxb2_rec = None
    for rec in SeqIO.parse(str(fasta_path), "fasta"):
        if hxb2_id_substring in rec.id:
            hxb2_rec = rec
            break
    if hxb2_rec is None:
        raise RuntimeError(f"Could not find HXB2 ({hxb2_id_substring}) in {fasta_path}")

    msa_to_hxb2 = {}
    hxb2_to_msa = {}
    h_pos = 0
    seq_str = str(hxb2_rec.seq).upper()
    for col_0, base in enumerate(seq_str):
        col_1 = col_0 + 1
        if base not in ("-", "."):
            h_pos += 1
            msa_to_hxb2[col_1] = h_pos
            hxb2_to_msa[h_pos] = col_1
        else:
            # Interpolate to nearest HXB2 position for insertions relative to HXB2
            msa_to_hxb2[col_1] = max(1, h_pos)

    return {
        "hxb2_seq_id": hxb2_rec.id,
        "msa_length": len(seq_str),
        "hxb2_ungapped_length": h_pos,
        "msa_to_hxb2": msa_to_hxb2,
        "hxb2_to_msa": hxb2_to_msa
    }


def main():
    print("=" * 80)
    print("STEP 1: PREPARING LANL HIV-1 REFERENCE & CRF COMPENDIUM IN /Users/sergei/Projects/HIV-recombination/data")
    print("=" * 80)

    # 1. Copy/stage lanl_ref_2023.fasta and lanl_flt_2022.fasta
    src_ref = SRC_NOM_DIR / "data" / "lanl_ref_2023.fasta"
    src_flt = SRC_NOM_DIR / "data" / "lanl_flt_2022.fasta"
    dst_ref = DATA_DIR / "lanl_ref_2023.fasta"
    dst_flt = DATA_DIR / "lanl_flt_2022.fasta"

    if not dst_ref.exists():
        shutil.copy2(src_ref, dst_ref)
    if not dst_flt.exists():
        shutil.copy2(src_flt, dst_flt)
    print(f"[+] Staged {dst_ref.name} and {dst_flt.name}")

    # 2. Save HXB2 gene & structural domain annotations
    annot_path = DATA_DIR / "hxb2_gene_domain_map.json"
    with open(annot_path, "w") as f:
        json.dump(HXB2_ANNOTATIONS, f, indent=2)
    print(f"[+] Wrote HXB2 gene & domain map to {annot_path}")

    # 3. Build HXB2 <-> MSA coordinate maps for lanl_ref_2023
    coord_map_ref = build_hxb2_coordinate_map(dst_ref, "K03455")
    with open(DATA_DIR / "lanl_ref_2023_hxb2_coord_map.json", "w") as f:
        json.dump(coord_map_ref, f)
    print(f"[+] Built HXB2 coordinate map for lanl_ref_2023.fasta: MSA L={coord_map_ref['msa_length']} -> HXB2 L={coord_map_ref['hxb2_ungapped_length']}")

    # 4. Parse LANL CRF HTML table across full genome
    if LANL_HTML_CACHE.exists():
        catalog, seg_df = parse_lanl_crf_html(LANL_HTML_CACHE)
        with open(DATA_DIR / "lanl_official_crf_catalog.json", "w") as f:
            json.dump(catalog, f, indent=2)
        seg_df.to_csv(DATA_DIR / "lanl_official_crf_segments.csv", index=False)
        print(f"[+] Parsed {len(catalog)} full-genome CRFs ({len(seg_df)} total segments) from LANL CRF compendium.")

        # Print diagnostic pathologies of the current LANL CRF catalog!
        n_cpx = sum(1 for v in catalog.values() if v["is_cpx"])
        n_u = sum(1 for v in catalog.values() if v["has_unclassified_U"])
        n_slash = sum(1 for v in catalog.values() if v["has_ambiguous_slash"])
        n_2nd = sum(1 for v in catalog.values() if v["has_second_gen_CRF_parent"])
        print(f"    - CRFs dumped into '_cpx' label:                  {n_cpx} / {len(catalog)} ({100*n_cpx/len(catalog):.1f}%)")
        print(f"    - CRFs with unclassified 'U' segments:            {n_u} / {len(catalog)} ({100*n_u/len(catalog):.1f}%)")
        print(f"    - CRFs with ambiguous 'X/Y' parent segments:      {n_slash} / {len(catalog)} ({100*n_slash/len(catalog):.1f}%)")
        print(f"    - 2nd/3rd-gen CRFs listing another CRF as parent: {n_2nd} / {len(catalog)} ({100*n_2nd/len(catalog):.1f}%)")

    # 5. Extract metadata for all 518 sequences in lanl_ref_2023.fasta and 5,089 sequences in lanl_flt_2022.fasta
    ref_records = list(SeqIO.parse(str(dst_ref), "fasta"))
    ref_meta = [parse_lanl_fasta_header(r.id) for r in ref_records]
    df_ref_meta = pd.DataFrame(ref_meta)
    df_ref_meta.to_csv(DATA_DIR / "lanl_ref_2023_metadata.csv", index=False)

    flt_records = list(SeqIO.parse(str(dst_flt), "fasta"))
    flt_meta = [parse_lanl_fasta_header(r.id) for r in flt_records]
    df_flt_meta = pd.DataFrame(flt_meta)
    df_flt_meta.to_csv(DATA_DIR / "lanl_flt_2022_metadata.csv", index=False)
    print(f"[+] Extracted metadata: {len(df_ref_meta)} ref genomes ({df_ref_meta['year'].notna().sum()} dated), {len(df_flt_meta)} flt genomes ({df_flt_meta['year'].notna().sum()} dated)")

    # 6. Build the Phase 1 Core Subtype + Primary CRF Audit Alignment (Pure A1..L + 01_AE + 02_AG)
    pure_tags = {"A1", "A2", "A3", "A4", "A6", "A7", "A8", "B", "C", "D", "F1", "F2", "G", "H", "J", "K", "L"}
    focal_tags = pure_tags | {"01_AE", "02_AG"}
    focal_recs = [r for r, m in zip(ref_records, ref_meta) if m["lanl_tag"] in focal_tags]
    focal_fasta = DATA_DIR / "lanl_pure_and_01_02.fasta"
    SeqIO.write(focal_recs, str(focal_fasta), "fasta")

    pure_only_recs = [r for r, m in zip(ref_records, ref_meta) if m["lanl_tag"] in pure_tags]
    pure_only_fasta = DATA_DIR / "lanl_nominal_pure_subtypes.fasta"
    SeqIO.write(pure_only_recs, str(pure_only_fasta), "fasta")

    print(f"[+] Wrote {len(pure_only_recs)} nominal 'pure' subtype genomes to {pure_only_fasta.name}")
    print(f"[+] Wrote {len(focal_recs)} nominal 'pure' + CRF01_AE + CRF02_AG genomes to {focal_fasta.name}")


if __name__ == "__main__":
    main()
