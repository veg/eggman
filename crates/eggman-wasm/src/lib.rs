use eggman_core::{
    CompendiumEntry, EggmanEngine, GeneLandmark, PathogenProfilePack, ScreeningConfig,
    SegmentedAssemblyResult, SequenceResult, TemporalFloorMode, TemporalGuardConfig,
};
use serde::Serialize;
use wasm_bindgen::prelude::*;

#[derive(Debug, Clone, Serialize)]
pub struct PackHitScores {
    pub hiv1_rate: f64,
    pub iav_rate: f64,
    pub sc2_rate: f64,
    pub sarbeco_rate: f64,
}

#[derive(Debug, Clone, Serialize)]
pub struct PathogenDetectionResult {
    pub detected_pack: String,
    pub pathogen_name: String,
    pub reference_name: String,
    pub is_segmented: bool,
    pub confidence: String,
    pub reason: String,
    pub raw_record_count: usize,
    pub hit_scores: PackHitScores,
}

#[derive(Debug, Clone, Serialize)]
pub struct PackMetaPayload {
    pub pack_id: String,
    pub engine_version: String,
    pub pack_version: String,
    pub snapshot_date: String,
    pub curation_sources: String,
    pub certification_criteria: String,
    pub pathogen_name: String,
    pub reference_name: String,
    pub ref_full_len: usize,
    pub coding_start: usize,
    pub coding_end: usize,
    pub is_segmented: bool,
    pub landmarks: Vec<GeneLandmark>,
    pub trunk_names: Vec<String>,
    pub certified_count: usize,
    pub compendium: Vec<CompendiumEntry>,
}

#[derive(Debug, Clone, Serialize)]
pub struct MultiPackScreenResponse {
    pub active_pack: String,
    pub detection: PathogenDetectionResult,
    pub pack_meta: PackMetaPayload,
    pub assembly: SegmentedAssemblyResult,
    pub assembled_fasta: String,
    pub results: Vec<SequenceResult>,
}

#[wasm_bindgen]
pub struct WasmEggman {
    hiv1: EggmanEngine,
    iav: EggmanEngine,
    sc2: EggmanEngine,
    sarbeco: EggmanEngine,
    custom: Option<EggmanEngine>,
    custom_pack: Option<PathogenProfilePack>,
}

impl WasmEggman {
    fn engine_for(&self, pack_id: &str) -> (&'static str, &EggmanEngine) {
        match pack_id.trim().to_ascii_lowercase().as_str() {
            "custom" => {
                if let Some(ref eng) = self.custom {
                    ("custom", eng)
                } else {
                    ("hiv1", &self.hiv1)
                }
            }
            "iav" | "flu" | "h5n1" | "influenza" | "influenza_a_8seg" => ("iav", &self.iav),
            "sc2" | "sars2" | "sars-cov-2" | "sars_cov2" | "sars_cov2_wuhan" => ("sc2", &self.sc2),
            "sarbeco" | "sarbecovirus" | "sarbecovirus_wuhan" => ("sarbeco", &self.sarbeco),
            _ => ("hiv1", &self.hiv1),
        }
    }

    fn build_pack_meta(&self, pack_id: &str) -> PackMetaPayload {
        let (canonical_id, eng) = self.engine_for(pack_id);
        let spec = eng.spec();
        PackMetaPayload {
            pack_id: canonical_id.to_string(),
            engine_version: env!("CARGO_PKG_VERSION").to_string(),
            pack_version: spec.pack_version.clone(),
            snapshot_date: spec.snapshot_date.clone(),
            curation_sources: spec.curation_sources.clone(),
            certification_criteria: spec.certification_criteria.clone(),
            pathogen_name: spec.pathogen_name.clone(),
            reference_name: spec.reference_name.clone(),
            ref_full_len: spec.ref_full_len,
            coding_start: spec.coding_start,
            coding_end: spec.coding_end,
            is_segmented: spec.is_segmented,
            landmarks: spec.landmarks.clone(),
            trunk_names: eng.trunk_names.clone(),
            certified_count: eng.certified_lineage_count(),
            compendium: eng.compendium.clone(),
        }
    }

    fn detect_pathogen_internal(&self, fasta_text: &str) -> PathogenDetectionResult {
        let raw_records: Vec<(String, String)> = EggmanEngine::parse_fasta(fasta_text)
            .into_iter()
            .filter_map(|(hdr, seq)| {
                let clean: String = seq
                    .chars()
                    .filter(|&c| c != '|' && c != '-' && c != '.')
                    .collect();
                if clean.is_empty() {
                    None
                } else {
                    Some((hdr, clean))
                }
            })
            .collect();

        let raw_record_count = raw_records.len();
        if raw_records.is_empty() {
            return PathogenDetectionResult {
                detected_pack: "hiv1".to_string(),
                pathogen_name: self.hiv1.pathogen_name().to_string(),
                reference_name: self.hiv1.reference_name().to_string(),
                is_segmented: false,
                confidence: "default".to_string(),
                reason: "Empty input; defaulting to HIV-1 Group M".to_string(),
                raw_record_count: 0,
                hit_scores: PackHitScores {
                    hiv1_rate: 0.0,
                    iav_rate: 0.0,
                    sc2_rate: 0.0,
                    sarbeco_rate: 0.0,
                },
            };
        }

        let sample_n = raw_records.len().min(16);
        let mut hiv_m = 0usize;
        let mut iav_m = 0usize;
        let mut sc2_m = 0usize;
        let mut sar_m = 0usize;
        let mut sum_valid = 0usize;
        let mut max_sc2_prof_rate = 0.0f64;
        let mut wildlife_cue = false;
        let mut flu_cue = false;
        let mut max_len = 0usize;

        for (hdr, seq) in raw_records.iter().take(sample_n) {
            max_len = max_len.max(seq.len());
            let low_hdr = hdr.to_ascii_lowercase();
            if low_hdr.contains("ratg13")
                || low_hdr.contains("banal")
                || low_hdr.contains("covzc")
                || low_hdr.contains("prc31")
                || low_hdr.contains("hku3")
                || low_hdr.contains("wiv16")
                || low_hdr.contains("rpshaanxi")
                || low_hdr.contains("bm48")
                || low_hdr.contains("sarbeco")
                || low_hdr.contains("pangolin")
                || low_hdr.contains("rhinolophus")
            {
                wildlife_cue = true;
            }
            if low_hdr.contains("h5n1")
                || low_hdr.contains("influenza")
                || low_hdr.contains("|pb2")
                || low_hdr.contains("|ha")
                || low_hdr.contains("|na")
                || low_hdr.contains("iav-")
            {
                flu_cue = true;
            }

            let (hm, hv, _) = self.hiv1.kmer_hit_stats(seq);
            let (im, _, _) = self.iav.kmer_hit_stats(seq);
            let (sm, sv, s_prof) = self.sc2.kmer_hit_stats(seq);
            let (rm, _, _) = self.sarbeco.kmer_hit_stats(seq);

            let valid = hv.max(sv).max(1);
            hiv_m += hm;
            iav_m += im;
            sc2_m += sm;
            sar_m += rm;
            sum_valid += valid;
            let p_rate = (s_prof as f64) / (valid as f64);
            if p_rate > max_sc2_prof_rate {
                max_sc2_prof_rate = p_rate;
            }
        }

        let denom = sum_valid.max(1) as f64;
        let hiv1_rate = (hiv_m as f64) / denom;
        let iav_rate = (iav_m as f64) / denom;
        let sc2_rate = (sc2_m as f64) / denom;
        let sarbeco_rate = (sar_m as f64) / denom;
        let cov_rate = sc2_rate.max(sarbeco_rate);

        let (detected_pack, confidence, reason) = if iav_rate >= 0.015
            && iav_rate > hiv1_rate * 1.25
            && iav_rate > cov_rate * 1.25
        {
            (
                "iav".to_string(),
                if iav_rate >= 0.15 { "high" } else { "medium" }.to_string(),
                format!(
                    "8-segment Influenza A 15-mer index matched {:.1}% of sampled kmers across {} FASTA record(s)",
                    iav_rate * 100.0,
                    raw_record_count
                ),
            )
        } else if cov_rate >= 0.015 && cov_rate > hiv1_rate * 1.25 && cov_rate > iav_rate * 1.25 {
            // Distinguish clinical SARS-CoV-2 from divergent wildlife Sarbecovirus
            let is_wildlife = wildlife_cue
                || sar_m > sc2_m.saturating_add(sum_valid / 200)
                || (sc2_rate < 0.34 && max_sc2_prof_rate < 0.36);
            if is_wildlife {
                (
                    "sarbeco".to_string(),
                    "high".to_string(),
                    format!(
                        "Wildlife Sarbecovirus signatures detected ({:.1}% Sarbeco vs {:.1}% SC2 15-mer match; max clinical SC2 profile {:.1}%)",
                        sarbeco_rate * 100.0,
                        sc2_rate * 100.0,
                        max_sc2_prof_rate * 100.0
                    ),
                )
            } else {
                (
                    "sc2".to_string(),
                    if sc2_rate >= 0.25 { "high" } else { "medium" }.to_string(),
                    format!(
                        "Clinical SARS-CoV-2 Wuhan-Hu-1 index matched {:.1}% of sampled kmers across {} genome(s)",
                        sc2_rate * 100.0,
                        raw_record_count
                    ),
                )
            }
        } else if hiv1_rate >= 0.008 && hiv1_rate >= iav_rate && hiv1_rate >= cov_rate {
            (
                "hiv1".to_string(),
                if hiv1_rate >= 0.06 { "high" } else { "medium" }.to_string(),
                format!(
                    "HIV-1 Group M HXB2 15-mer index matched {:.1}% of sampled kmers across {} sequence(s)",
                    hiv1_rate * 100.0,
                    raw_record_count
                ),
            )
        } else if flu_cue {
            (
                "iav".to_string(),
                "low".to_string(),
                "Header keywords matched Influenza A; very low 15-mer identity".to_string(),
            )
        } else if max_len >= 22_000 {
            if wildlife_cue {
                (
                    "sarbeco".to_string(),
                    "low".to_string(),
                    "~29 kb wildlife coronavirus genome inferred from sequence length/header".to_string(),
                )
            } else {
                (
                    "sc2".to_string(),
                    "low".to_string(),
                    "~29.9 kb coronavirus genome inferred from sequence length".to_string(),
                )
            }
        } else {
            (
                "hiv1".to_string(),
                "low".to_string(),
                "Low 15-mer overlap across packs; defaulting to HIV-1 Group M".to_string(),
            )
        };

        let (_, eng) = self.engine_for(&detected_pack);
        PathogenDetectionResult {
            detected_pack,
            pathogen_name: eng.pathogen_name().to_string(),
            reference_name: eng.reference_name().to_string(),
            is_segmented: eng.spec().is_segmented,
            confidence,
            reason,
            raw_record_count,
            hit_scores: PackHitScores {
                hiv1_rate,
                iav_rate,
                sc2_rate,
                sarbeco_rate,
            },
        }
    }

    fn make_config(
        gamma: f64,
        min_span_nt: usize,
        min_informative_snps: usize,
        p_fisher_thresh: f64,
        mask_hypervariable_loops: bool,
    ) -> ScreeningConfig {
        ScreeningConfig {
            gamma: if gamma.is_finite() && gamma > 0.0 {
                Some(gamma)
            } else {
                None
            },
            min_span_nt: if min_span_nt > 0 {
                Some(min_span_nt)
            } else {
                None
            },
            min_informative_snps: if min_informative_snps > 0 {
                min_informative_snps
            } else {
                5
            },
            p_fisher_thresh: if p_fisher_thresh.is_finite()
                && p_fisher_thresh > 0.0
                && p_fisher_thresh <= 1.0
            {
                p_fisher_thresh
            } else {
                1e-4
            },
            mask_hypervariable_loops,
        }
    }
}

#[wasm_bindgen]
impl WasmEggman {
    #[wasm_bindgen(constructor)]
    pub fn new() -> Self {
        Self {
            hiv1: EggmanEngine::new(),
            iav: EggmanEngine::from_builtin_id("iav").expect("valid embedded IAV pack"),
            sc2: EggmanEngine::from_builtin_id("sc2").expect("valid embedded SC2 pack"),
            sarbeco: EggmanEngine::from_builtin_id("sarbeco").expect("valid embedded Sarbeco pack"),
            custom: None,
            custom_pack: None,
        }
    }

    /// Load a custom `PathogenProfilePack` JSON (or bare `PathogenSpec` + optional FASTA profiles string) into slot `"custom"`.
    pub fn load_custom_pack(&mut self, json_str: &str, fasta_override: &str) -> String {
        match PathogenProfilePack::from_json_and_fasta_str(json_str, fasta_override) {
            Ok(pack) => match EggmanEngine::from_pathogen_pack(pack.clone()) {
                Ok(eng) => {
                    self.custom = Some(eng);
                    self.custom_pack = Some(pack);
                    let meta = self.build_pack_meta("custom");
                    serde_json::json!({ "ok": true, "meta": meta }).to_string()
                }
                Err(e) => serde_json::json!({ "ok": false, "error": e }).to_string(),
            },
            Err(e) => serde_json::json!({ "ok": false, "error": e }).to_string(),
        }
    }

    /// Export any built-in pack (`hiv1`, `iav`, `sc2`, `sarbeco`) or loaded `"custom"` pack as a self-contained JSON bundle.
    pub fn export_pack_bundle_json(&self, pack_id: &str) -> String {
        let low = pack_id.trim().to_ascii_lowercase();
        if low == "custom" {
            if let Some(ref pack) = self.custom_pack {
                return serde_json::to_string_pretty(pack).unwrap_or_else(|_| "{}".to_string());
            }
        }
        match PathogenProfilePack::from_builtin_id(&low) {
            Ok(pack) => serde_json::to_string_pretty(&pack).unwrap_or_else(|_| "{}".to_string()),
            Err(e) => serde_json::json!({ "error": e }).to_string(),
        }
    }

    /// Inspect FASTA records and automatically infer whether the dataset is HIV-1, IAV, SARS-CoV-2, or Sarbecovirus.
    pub fn detect_pathogen(&self, fasta_text: &str) -> String {
        let det = self.detect_pathogen_internal(fasta_text);
        serde_json::to_string(&det).unwrap_or_else(|_| "{}".to_string())
    }

    /// Return metadata and reference compendium JSON for a built-in pack (`hiv1`, `iav`, `sc2`, `sarbeco`).
    pub fn get_pack_meta_json(&self, pack_id: &str) -> String {
        let meta = self.build_pack_meta(pack_id);
        serde_json::to_string(&meta).unwrap_or_else(|_| "{}".to_string())
    }

    /// Screen a FASTA string using `"auto"` pathogen detection or an explicit pack id,
    /// automatically assembling multi-record Influenza A segments (`seg_mode`: `"auto"`, `"single_isolate"`, `"per_record"`).
    pub fn screen_with_pack(
        &self,
        pack_id: &str,
        fasta_text: &str,
        seg_mode: &str,
        gamma: f64,
        min_span_nt: usize,
        min_informative_snps: usize,
        p_fisher_thresh: f64,
        mask_hypervariable_loops: bool,
    ) -> String {
        self.screen_with_pack_and_temporal(
            pack_id,
            fasta_text,
            seg_mode,
            gamma,
            min_span_nt,
            min_informative_snps,
            p_fisher_thresh,
            mask_hypervariable_loops,
            false,
            "sampled",
            f64::NAN,
            f64::NAN,
        )
    }

    /// Screen a FASTA string with optional ChronAeon non-lookahead temporal guarding
    /// (`temporal_mode`: `"sampled"` or `"submission"`; pass `<= 0` or `NaN` for `query_year` to parse from FASTA headers,
    /// and `< 0` or `NaN` for `temporal_slack_yr` to use pack/mode defaults).
    #[allow(clippy::too_many_arguments)]
    pub fn screen_with_pack_and_temporal(
        &self,
        pack_id: &str,
        fasta_text: &str,
        seg_mode: &str,
        gamma: f64,
        min_span_nt: usize,
        min_informative_snps: usize,
        p_fisher_thresh: f64,
        mask_hypervariable_loops: bool,
        temporal_guard: bool,
        temporal_mode: &str,
        query_year: f64,
        temporal_slack_yr: f64,
    ) -> String {
        let detection = self.detect_pathogen_internal(fasta_text);
        let chosen_id = if pack_id.trim().eq_ignore_ascii_case("auto") || pack_id.trim().is_empty()
        {
            detection.detected_pack.as_str()
        } else {
            pack_id.trim()
        };
        let (canonical_id, engine) = self.engine_for(chosen_id);
        let cfg = Self::make_config(
            gamma,
            min_span_nt,
            min_informative_snps,
            p_fisher_thresh,
            mask_hypervariable_loops,
        );
        let mode = match temporal_mode.trim().to_ascii_lowercase().as_str() {
            "submission" | "submitted" => TemporalFloorMode::Submission,
            _ => TemporalFloorMode::Sampled,
        };
        let temporal_cfg = TemporalGuardConfig {
            enabled: temporal_guard,
            query_year: if query_year.is_finite() && query_year > 0.0 {
                Some(query_year)
            } else {
                None
            },
            slack_yr: if temporal_slack_yr.is_finite() && temporal_slack_yr >= 0.0 {
                Some(temporal_slack_yr)
            } else {
                None
            },
            mode,
        };
        let (assembly, results, assembled_fasta) = engine.classify_fasta_with_segmentation_temporal(
            fasta_text,
            seg_mode,
            &cfg,
            &temporal_cfg,
        );
        let payload = MultiPackScreenResponse {
            active_pack: canonical_id.to_string(),
            detection,
            pack_meta: self.build_pack_meta(canonical_id),
            assembly,
            assembled_fasta,
            results,
        };
        serde_json::to_string(&payload).unwrap_or_else(|_| "{}".to_string())
    }

    /// Legacy HIV-1 FASTA screening entry point preserved for backward compatibility.
    pub fn screen_fasta(
        &self,
        fasta_text: &str,
        gamma: f64,
        min_span_nt: usize,
        min_informative_snps: usize,
        p_fisher_thresh: f64,
        mask_hypervariable_loops: bool,
    ) -> String {
        let records = EggmanEngine::parse_fasta(fasta_text);
        let cfg = Self::make_config(
            gamma,
            min_span_nt,
            min_informative_snps,
            p_fisher_thresh,
            mask_hypervariable_loops,
        );
        let results: Vec<_> = records
            .iter()
            .map(|(hdr, seq)| self.hiv1.classify_sequence(hdr, seq, &cfg))
            .collect();
        serde_json::to_string(&results).unwrap_or_else(|_| "[]".to_string())
    }

    /// Return the embedded HIV-1 171-lineage reference compendium as a JSON array.
    pub fn get_compendium_json(&self) -> String {
        serde_json::to_string(&self.hiv1.compendium).unwrap_or_else(|_| "[]".to_string())
    }
}

#[cfg(test)]
mod wasm_multipathogen_tests {
    use super::*;

    #[test]
    fn autodetects_all_four_empirical_demo_fastas_and_stiches_segments() {
        let wasm = WasmEggman::new();

        let demos = [
            (
                "hiv1",
                include_str!("../../../docs/examples/demos.fasta"),
                5usize,
            ),
            (
                "hiv1",
                include_str!("../../../docs/examples/partial_env.fasta"),
                5usize,
            ),
            (
                "iav",
                include_str!("../../../docs/examples/iav_reassortment.fasta"),
                7usize,
            ),
            (
                "iav",
                include_str!("../../../docs/examples/iav_raw_8segments.fasta"),
                2usize,
            ),
            (
                "sc2",
                include_str!("../../../docs/examples/sars_cov2_recomb.fasta"),
                7usize,
            ),
            (
                "sarbeco",
                include_str!("../../../docs/examples/sarbecovirus_wildlife.fasta"),
                6usize,
            ),
        ];

        for (expected_pack, fasta, expected_isolates) in demos {
            let det = wasm.detect_pathogen_internal(fasta);
            assert_eq!(
                det.detected_pack, expected_pack,
                "failed detection on expected pack {} (got {:?})",
                expected_pack, det
            );

            let out_json = wasm.screen_with_pack("auto", fasta, "auto", 0.0, 0, 5, 1e-4, true);
            let parsed: serde_json::Value = serde_json::from_str(&out_json).unwrap();
            assert_eq!(parsed["active_pack"], expected_pack);
            let res_arr = parsed["results"].as_array().unwrap();
            assert_eq!(res_arr.len(), expected_isolates);

            if expected_pack == "iav" && expected_isolates == 2 {
                assert_eq!(parsed["assembly"]["was_multi_record_assembled"], true);
                assert_eq!(res_arr[0]["called_lineage"], "B3.13");
                assert_eq!(res_arr[1]["called_lineage"], "D1.1");
            }

            // Verify pure sequence-only detection when FASTA headers are completely anonymized
            let mut anon_fasta = String::new();
            for (i, (_, seq)) in EggmanEngine::parse_fasta(fasta).into_iter().enumerate() {
                anon_fasta.push_str(&format!(">ANON_{}\n{}\n", i + 1, seq));
            }
            let anon_det = wasm.detect_pathogen_internal(&anon_fasta);
            assert_eq!(
                anon_det.detected_pack, expected_pack,
                "failed headerless sequence-only detection for {} (got {:?})",
                expected_pack, anon_det
            );
        }
    }

    #[test]
    fn custom_pack_bundle_export_and_reload_works() {
        let mut wasm = WasmEggman::new();
        let bundle = wasm.export_pack_bundle_json("sc2");
        let load_res: serde_json::Value = serde_json::from_str(&wasm.load_custom_pack(&bundle, "")).unwrap();
        assert_eq!(load_res["ok"], true, "custom pack load failed: {:?}", load_res);
        assert_eq!(load_res["meta"]["pack_id"], "custom");

        let fasta = include_str!("../../../docs/examples/sars_cov2_recomb.fasta");
        let out_json = wasm.screen_with_pack("custom", fasta, "auto", 3.5, 0, 3, 0.01, false);
        let parsed: serde_json::Value = serde_json::from_str(&out_json).unwrap();
        assert_eq!(parsed["active_pack"], "custom");
        assert_eq!(parsed["results"].as_array().unwrap().len(), 7);
    }

    #[test]
    fn wasm_temporal_guard_masks_future_clades() {
        let wasm = WasmEggman::new();
        let fasta = include_str!("../../../docs/examples/sars_cov2_recomb.fasta");
        let out_json = wasm.screen_with_pack_and_temporal(
            "sc2",
            fasta,
            "auto",
            0.0,
            0,
            5,
            1e-4,
            true,
            true,
            "submission",
            2020.5,
            0.0,
        );
        let parsed: serde_json::Value = serde_json::from_str(&out_json).unwrap();
        let res_arr = parsed["results"].as_array().unwrap();
        assert_eq!(res_arr.len(), 7);
        let masked = res_arr[0]["masked_future_trunks"].as_array().unwrap();
        assert!(
            !masked.is_empty(),
            "expected post-2020.5 SARS-CoV-2 founders to be masked under submission temporal guard"
        );
    }
}
