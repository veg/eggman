# Reproducible Analysis & Benchmarking Scripts

This directory archives the Python analysis, validation, surveillance streaming, and publication-figure scripts supporting the **EGGMAN** manuscript.

## Pipeline Overview

- **HIV-1 Reference Scrubbing, ChronAeon Weighting, & LANL Audit (`01`–`07`, `16`, `19`)**
  - `01_prepare_lanl_compendium.py`, `02_phase1_pure_subtype_and_temporal_audit.py`, `02b_build_hiv1_chronaeon_founders.py`: Basal trunk scrubbing (`G0`, `L`, `?E`), `ChronAeon` emergence profiling, and pure-subtype audit.
  - `03_audit_and_reclassify_all_crfs.py`, `05_reconstruct_recombination_dag.py`, `06_validate_all_lanl_full_genomes.py`: Bridged audit of 986 LANL CRF junctions, demotions, and 3-generation recombination DAG.
  - `16_build_rhiztype_certified_profiles.py`, `16_validate_west_african_crf02_rt.py`, `19_chronaeon_lsd2_treetime_validation.py`: Tier-2 certified catalog construction and `IQ-TREE 2` / `LSD2` / `TreeTime` clock validation.
- **GenBank & Stanford HIVDB Archive Audits (`08`–`12`, `15`)**
  - `08_scan_genbank_31k_full_genomes.py`, `08b_fetch_genbank_cluster_metadata.py`, `09_generate_segment_phylogenies.py`: Full-length GenBank deduplication (`31,825` seqs / `7,720` patients), 3-class lineage/hotspot separation (`RU_G0B_gagMA` vs. convergent hotspots), and segment `IQ-TREE 2` validation.
  - `11_screen_partial_pol_and_env.py`, `12_scale_genbank_partial_pol_env_screen.py`, `15_audit_stanford_hivdb_subtypes.py`: Partial-only `pol`/`env` census (`38,250` seqs / `20,797` hosts) and Stanford HIVDB (`120,468` clinical records + `720` references).
- **HIV-1 Benchmarks, Null Simulations, & Two-Stage Sieve (`10`, `14`, `17`, `18`, `23`)**
  - `14_benchmark_rhiztype_vs_comet.py`, `17_benchmark_strict_and_jphmm.py`, `18_simulate_tier1_null_fpr.py`: 12-panel HIV-1 head-to-head vs. `COMET v2.4` and 28-worker `jpHMM`, 500-mosaic synthetic breakpoint localization, LOSO/OOD holdouts, and 1,300 non-recombinant null simulations.
  - `10_recombination_viability_and_establishment_sieve.py`, `23_validate_usda_h5n1_and_sieve_controls.py`: National `H_circ` supply models (WLS, HC3, NB), discrete Hurwitz-zeta MLE (`N = 761`), Stage 1 folding-core sieve (3 nulls including 1,000-sim detectability null), and Stage 2 transmission sieve (CMH, permutations, logistic regression).
- **Multi-Pathogen Packs, Streaming Surveillance, & Comparators (`21`, `24`–`27`)**
  - `21_sars_cov2_spanning_basis_and_bvbrc_screen.py`, `24_bvbrc_100k_streaming_census_and_screen.py`, `25_validate_bvbrc_discoveries_orthogonally.py`: SARS-CoV-2 16-founder Spanning Basis Rule, 100,000-genome live BV-BRC stream (including post-2021 Wuhan-`B` backfill QA filter), 3,980-isolate BV-BRC H5Nx stream, and orthogonal ML verification.
  - `26_unified_chronaeon_temporal_basis_sc2_iav.py`, `27_benchmark_sc2_nextclade_and_iav_genoflu.py`: `ChronAeon` SC2/IAV dating, `--temporal-guard` (`sampled`/`submission`) ablations, and head-to-head benchmarks against `USDA GenoFLU` and `Nextclade CLI`.
- **Publication Figures & Audit (`20`, `22`, `audit_figures_standard.py`)**
  - `20_export_supplementary_manifests.py`, `22_generate_science_figures.py`, `audit_figures_standard.py`: Supplementary manifest generation and automated 1:1 publication figure rendering & overlap auditing.
