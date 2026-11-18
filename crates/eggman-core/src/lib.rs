//! Core reference-anchored viral subtyping and mosaic deconvolution engine.
//!
//! Defaults to HIV-1 Group M anchored on HXB2 (`K03455`, 9,719 nt), and supports
//! arbitrary user-supplied pathogen profile packs (contiguous genomes such as SARS-CoV-2
//! on Wuhan-Hu-1, or segmented genomes such as Influenza A 8-segment constellations)
//! via seeded 15-mer diagonal anchor hashing + banded Needleman-Wunsch coordinate projection.

use serde::{Deserialize, Serialize};
use std::cmp::Reverse;
use std::collections::{BTreeMap, HashMap, HashSet};
use std::fs;
use std::path::Path;

const PROFILES_TSV: &str = include_str!("../data/profiles_9719.tsv");
const FOUNDER_PROFILES_TSV: &str = include_str!("../data/founder_profiles_9719.tsv");
const CERTIFIED_PROFILES_TSV: &str = include_str!("../data/certified_profiles_8628.tsv");
const COMPENDIUM_JSON: &str = include_str!("../data/compendium_bundle.json");

const IAV_MANIFEST_JSON: &str = include_str!("../../../packs/influenza_a_8seg/manifest.json");
const IAV_PROFILES_FASTA: &str = include_str!("../../../packs/influenza_a_8seg/profiles.fasta");
const IAV_CERTIFIED_TSV: &str = include_str!("../../../packs/influenza_a_8seg/certified_profiles.tsv");
const IAV_COMPENDIUM_JSON: &str = include_str!("../../../packs/influenza_a_8seg/compendium.json");

const SC2_MANIFEST_JSON: &str = include_str!("../../../packs/sars_cov2_wuhan/manifest.json");
const SC2_PROFILES_FASTA: &str = include_str!("../../../packs/sars_cov2_wuhan/profiles.fasta");
const SC2_CERTIFIED_TSV: &str = include_str!("../../../packs/sars_cov2_wuhan/certified_profiles.tsv");
const SC2_COMPENDIUM_JSON: &str = include_str!("../../../packs/sars_cov2_wuhan/compendium.json");

const SARBECO_MANIFEST_JSON: &str = include_str!("../../../packs/sarbecovirus_wuhan/manifest.json");
const SARBECO_PROFILES_FASTA: &str = include_str!("../../../packs/sarbecovirus_wuhan/profiles.fasta");

pub const HXB2_FULL_LEN: usize = 9719;
pub const CODING_START_HXB2: usize = 790;
pub const CODING_END_HXB2: usize = 9417;
pub const CODING_LEN: usize = CODING_END_HXB2 - CODING_START_HXB2 + 1; // 8628
pub const PURE_BASAL_COUNT: usize = 13; // A1, A2, A6, B, C, D, F1, F2, G0, H, J, K, L (excludes ?E)

/// Hypervariable indel-rich loops in `gag` (p6 PTAP) and `env` (`V1/V2`, `V3`, `V4`, `V5`)
const HYPERVARIABLE_MASK_RANGES: &[(usize, usize)] = &[
    (2130, 2185),
    (6615, 6812),
    (7110, 7217),
    (7377, 7478),
    (7596, 7637),
];

fn default_true() -> bool {
    true
}

fn default_min_mapped_nt() -> usize {
    100
}

fn default_epi_founder_penalty() -> f64 {
    0.005
}

/// Gene, domain, or RNA segment landmark along a pathogen's 1-based reference coordinates.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GeneLandmark {
    pub short_name: String,
    #[serde(default)]
    pub parent_gene: Option<String>,
    #[serde(default)]
    pub parent_min_span: usize,
    pub start: usize,
    pub end: usize,
}

/// Coordinate rule classifying partial amplicon locus scopes (e.g. `@pol`, `@env`, `@S`).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LocusScopeRule {
    pub name: String,
    pub min_start: usize,
    #[serde(default)]
    pub min_end: usize,
    pub max_end: usize,
    #[serde(default)]
    pub min_span_nt: usize,
}

/// Definition of a basal evolutionary trunk (or domain-restricted ancestral ghost trunk).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TrunkSpec {
    pub name: String,
    #[serde(default = "default_true")]
    pub is_pure_basal: bool,
    #[serde(default)]
    pub penalty: f64,
    #[serde(default)]
    pub active_range: Option<(usize, usize)>,
    #[serde(default)]
    pub fallback_profile: Option<String>,
    #[serde(default)]
    pub display_alias: Option<String>,
    #[serde(default)]
    pub clonal_certified_label: Option<String>,
    #[serde(default)]
    pub t_mrca: Option<f64>,
    #[serde(default)]
    pub t_mrca_ci: Option<(f64, f64)>,
    #[serde(default)]
    pub first_sampled_year: Option<f64>,
    #[serde(default)]
    pub chronaeon_half_life_yr: Option<f64>,
}

impl TrunkSpec {
    /// Admissible emergence floor `min(t_mrca_ci_low, t_mrca, first_sampled_year)`.
    /// Returns `None` if no finite temporal bounds exist for this trunk.
    pub fn emergence_floor(&self) -> Option<f64> {
        min_temporal_floor(
            self.t_mrca_ci.map(|(lo, _)| lo),
            self.t_mrca,
            self.first_sampled_year,
        )
    }
}

fn min_temporal_floor(
    t_mrca_ci_low: Option<f64>,
    t_mrca: Option<f64>,
    first_sampled_year: Option<f64>,
) -> Option<f64> {
    let mut floor: Option<f64> = None;
    for candidate in [t_mrca_ci_low, t_mrca, first_sampled_year]
        .into_iter()
        .flatten()
    {
        if candidate.is_finite() {
            floor = Some(match floor {
                Some(prev) => prev.min(candidate),
                None => candidate,
            });
        }
    }
    floor
}

/// Rule for rescuing unattributed ghost segments that match a known extinct/unsampled donor lineage.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GhostRelabelRule {
    pub ghost_trunk: String,
    pub donor_profile: String,
    pub partner_trunk: String,
    pub min_ref_coord: usize,
    pub min_margin_per_site: f64,
    pub min_span_nt: usize,
}

/// Pathogen architecture, reference coordinate frame, gene map, and lineage classification rules.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PathogenSpec {
    pub pathogen_name: String,
    pub reference_name: String,
    pub ref_full_len: usize,
    pub coding_start: usize,
    pub coding_end: usize,
    #[serde(default)]
    pub max_query_len: usize,
    #[serde(default = "default_min_mapped_nt")]
    pub min_mapped_coding_nt: usize,
    #[serde(default)]
    pub full_span_min_nt: usize,
    #[serde(default)]
    pub full_mapped_min_nt: usize,
    #[serde(default)]
    pub is_segmented: bool,
    #[serde(default)]
    pub forbid_intra_segment_recomb: bool,
    #[serde(default)]
    pub use_hiv1_legacy_gene_mapper: bool,
    #[serde(default)]
    pub hypervariable_masks: Vec<(usize, usize)>,
    #[serde(default)]
    pub landmarks: Vec<GeneLandmark>,
    #[serde(default)]
    pub locus_scopes: Vec<LocusScopeRule>,
    pub basal_trunks: Vec<TrunkSpec>,
    #[serde(default)]
    pub ghost_relabel_rule: Option<GhostRelabelRule>,
    #[serde(default)]
    pub align_founder_keys: Vec<String>,
    #[serde(default)]
    pub epi_founder_order: Vec<(String, String)>,
    #[serde(default = "default_epi_founder_penalty")]
    pub epi_founder_penalty: f64,
    #[serde(default)]
    pub epi_certified_replacements: Vec<(String, String, String)>,
    #[serde(default)]
    pub pack_version: String,
    #[serde(default)]
    pub snapshot_date: String,
    #[serde(default)]
    pub curation_sources: String,
    #[serde(default)]
    pub certification_criteria: String,
    #[serde(default)]
    pub default_temporal_slack_yr: Option<f64>,
    #[serde(default)]
    pub chronaeon_half_life_yr: Option<f64>,
}

impl PathogenSpec {
    /// Effective admissible lookahead tolerance (years) for `--temporal-guard`.
    /// Defaults to `3.0` yr for HIV-1 packs and `0.5` yr for non-HIV external packs.
    pub fn effective_temporal_slack_yr(&self) -> f64 {
        if let Some(s) = self.default_temporal_slack_yr {
            if s.is_finite() && s >= 0.0 {
                return s;
            }
        }
        if self.pathogen_name.to_ascii_lowercase().contains("hiv") {
            3.0
        } else {
            0.5
        }
    }

    /// Embedded HIV-1 Group M reference specification (HXB2 9,719 nt / coding 790..9417).
    pub fn hiv1_group_m_default() -> Self {
        let pure_meta: [(&str, f64, Option<(f64, f64)>, f64); 13] = [
            ("A1", 1959.7, Some((1954.1, 1964.1)), 1985.0),
            ("A2", 1987.0, None, 1987.0),
            ("A6", 1982.8, Some((1971.4, 1989.3)), 1999.0),
            ("B", 1967.8, Some((1965.1, 1970.2)), 1979.0),
            ("C", 1951.4, Some((1945.0, 1956.5)), 1986.0),
            ("D", 1966.1, Some((1958.8, 1971.2)), 1983.0),
            ("F1", 1968.5, Some((1952.7, 1977.3)), 1989.0),
            ("F2", 1983.8, Some((1974.3, 1989.0)), 1995.0),
            ("G0", 1961.5, Some((1953.4, 1967.5)), 1987.0),
            ("H", 1982.0, Some((1967.0, 1987.7)), 1990.0),
            ("J", 1986.0, Some((1978.1, 1989.6)), 1993.0),
            ("K", 1996.0, None, 1996.0),
            ("L", 1983.0, None, 1983.0),
        ];
        let mut basal_trunks: Vec<TrunkSpec> = pure_meta
            .iter()
            .map(|&(name, t_mrca, t_mrca_ci, first_sampled)| TrunkSpec {
                name: name.to_string(),
                is_pure_basal: true,
                penalty: 0.0,
                active_range: None,
                fallback_profile: None,
                display_alias: if name == "G0" {
                    Some("G".to_string())
                } else {
                    None
                },
                clonal_certified_label: if name == "G0" {
                    Some("Subtype_G".to_string())
                } else {
                    None
                },
                t_mrca: Some(t_mrca),
                t_mrca_ci,
                first_sampled_year: Some(first_sampled),
                chronaeon_half_life_yr: Some(6.93),
            })
            .collect();
        basal_trunks.push(TrunkSpec {
            name: "?E".to_string(),
            is_pure_basal: false,
            penalty: 0.035,
            active_range: Some((5097, CODING_END_HXB2)),
            fallback_profile: Some("CM240".to_string()),
            display_alias: None,
            clonal_certified_label: None,
            t_mrca: Some(1974.7),
            t_mrca_ci: Some((1970.6, 1978.0)),
            first_sampled_year: Some(1990.0),
            chronaeon_half_life_yr: Some(6.93),
        });

        let landmarks = vec![
            GeneLandmark { short_name: "5'LTR".into(), parent_gene: None, parent_min_span: 0, start: 1, end: 789 },
            GeneLandmark { short_name: "gag.MA".into(), parent_gene: Some("gag".into()), parent_min_span: 900, start: 790, end: 1185 },
            GeneLandmark { short_name: "gag.CA".into(), parent_gene: Some("gag".into()), parent_min_span: 900, start: 1186, end: 1878 },
            GeneLandmark { short_name: "gag.NC".into(), parent_gene: Some("gag".into()), parent_min_span: 900, start: 1879, end: 2292 },
            GeneLandmark { short_name: "PR".into(), parent_gene: Some("pol".into()), parent_min_span: 1800, start: 2293, end: 2549 },
            GeneLandmark { short_name: "RT".into(), parent_gene: Some("pol".into()), parent_min_span: 1800, start: 2550, end: 4229 },
            GeneLandmark { short_name: "IN".into(), parent_gene: Some("pol".into()), parent_min_span: 1800, start: 4230, end: 5096 },
            GeneLandmark { short_name: "vif".into(), parent_gene: None, parent_min_span: 0, start: 5097, end: 5559 },
            GeneLandmark { short_name: "vpr".into(), parent_gene: None, parent_min_span: 0, start: 5560, end: 5850 },
            GeneLandmark { short_name: "tat1".into(), parent_gene: None, parent_min_span: 0, start: 5851, end: 6062 },
            GeneLandmark { short_name: "vpu".into(), parent_gene: None, parent_min_span: 0, start: 6063, end: 6225 },
            GeneLandmark { short_name: "gp120".into(), parent_gene: Some("env".into()), parent_min_span: 1800, start: 6226, end: 7757 },
            GeneLandmark { short_name: "gp41".into(), parent_gene: Some("env".into()), parent_min_span: 1800, start: 7758, end: 8795 },
            GeneLandmark { short_name: "nef".into(), parent_gene: None, parent_min_span: 0, start: 8796, end: 9417 },
            GeneLandmark { short_name: "3'LTR".into(), parent_gene: None, parent_min_span: 0, start: 9418, end: HXB2_FULL_LEN },
        ];

        Self {
            pathogen_name: "HIV-1 Group M".to_string(),
            reference_name: "HXB2".to_string(),
            ref_full_len: HXB2_FULL_LEN,
            coding_start: CODING_START_HXB2,
            coding_end: CODING_END_HXB2,
            max_query_len: 20_000,
            min_mapped_coding_nt: 100,
            full_span_min_nt: 6800,
            full_mapped_min_nt: 6000,
            is_segmented: false,
            forbid_intra_segment_recomb: false,
            use_hiv1_legacy_gene_mapper: true,
            hypervariable_masks: HYPERVARIABLE_MASK_RANGES.to_vec(),
            landmarks,
            locus_scopes: vec![
                LocusScopeRule {
                    name: "pol".to_string(),
                    min_start: 1800,
                    min_end: 2293,
                    max_end: 5400,
                    min_span_nt: 0,
                },
                LocusScopeRule {
                    name: "env".to_string(),
                    min_start: 5900,
                    min_end: 0,
                    max_end: 9100,
                    min_span_nt: 0,
                },
            ],
            basal_trunks,
            ghost_relabel_rule: Some(GhostRelabelRule {
                ghost_trunk: "?E".to_string(),
                donor_profile: "CM240".to_string(),
                partner_trunk: "A1".to_string(),
                min_ref_coord: 5000,
                min_margin_per_site: 0.85,
                min_span_nt: 50,
            }),
            align_founder_keys: vec![
                "CRF01_AE".into(),
                "CRF02_AG".into(),
                "Subtype_G".into(),
                "CRF06_cpx".into(),
                "CRF07_BC".into(),
            ],
            epi_founder_order: [
                "CRF01_AE",
                "CRF02_AG",
                "Subtype_G",
                "CRF06_cpx",
                "CRF07_BC",
                "CRF08_BC",
                "CRF11_cpx",
                "CRF12_BF",
                "CRF13_cpx",
                "CRF15_01B",
                "CRF33_01B",
                "CRF55_01B",
                "CRF63_02A6",
            ]
            .into_iter()
            .map(|k| (k.to_string(), k.to_string()))
            .collect(),
            epi_founder_penalty: 0.005,
            epi_certified_replacements: vec![(
                "G0".to_string(),
                "Subtype_G".to_string(),
                "G".to_string(),
            )],
            pack_version: "2026.11-lanl2023".to_string(),
            snapshot_date: "2026-11-16".to_string(),
            curation_sources: "LANL HIV Sequence Compendium 2023 (55 pure-subtype refs; 139 verified CRFs, 6 clonal-collapsed CRFs demoted) + NCBI GenBank clinical census (278,829 sequences; 38 novel circulating clusters)".to_string(),
            certification_criteria: "4-gate guarded caller: >=95% span coverage, binomial parent dominance (n_c - n_u >= 5), stem synapomorphies (net_stem >= 8, p <= 1e-4), and Kadane residual scan (< 80 nt / max S < 12)".to_string(),
            default_temporal_slack_yr: Some(3.0),
            chronaeon_half_life_yr: Some(6.93),
        }
    }

    pub fn coding_len(&self) -> usize {
        self.coding_end.saturating_sub(self.coding_start) + 1
    }

    pub fn effective_max_query_len(&self) -> usize {
        if self.max_query_len > 0 {
            self.max_query_len
        } else {
            (self.ref_full_len * 3 / 2).max(20_000)
        }
    }

    pub fn effective_full_span_min_nt(&self) -> usize {
        if self.full_span_min_nt > 0 {
            self.full_span_min_nt
        } else {
            (self.coding_len() * 3) / 4
        }
    }

    pub fn effective_full_mapped_min_nt(&self) -> usize {
        if self.full_mapped_min_nt > 0 {
            self.full_mapped_min_nt
        } else {
            (self.coding_len() * 2) / 3
        }
    }

    pub fn classify_locus_scope(
        &self,
        h_start: usize,
        h_end: usize,
        span_width: usize,
        mapped_coding_nt: usize,
    ) -> String {
        if span_width >= self.effective_full_span_min_nt()
            && mapped_coding_nt >= self.effective_full_mapped_min_nt()
        {
            return "full".to_string();
        }
        for rule in &self.locus_scopes {
            if h_start >= rule.min_start
                && h_end >= rule.min_end
                && h_end <= rule.max_end
                && span_width >= rule.min_span_nt
            {
                return rule.name.clone();
            }
        }
        "partial".to_string()
    }

    fn nearest_landmark(&self, pos: usize) -> Option<&GeneLandmark> {
        if let Some(exact) = self
            .landmarks
            .iter()
            .find(|l| pos >= l.start && pos <= l.end)
        {
            return Some(exact);
        }
        self.landmarks.iter().min_by_key(|l| {
            if pos < l.start {
                l.start - pos
            } else if pos > l.end {
                pos - l.end
            } else {
                0
            }
        })
    }

    pub fn interval_to_short_gene(&self, h0: usize, h1: usize) -> String {
        if self.use_hiv1_legacy_gene_mapper {
            return hxb2_to_short_gene(h0, h1);
        }
        if self.landmarks.is_empty() {
            return format!("{}-{}", h0, h1);
        }
        if self.is_segmented {
            let mut hit_names: Vec<String> = Vec::new();
            for lm in &self.landmarks {
                let ovl_s = h0.max(lm.start);
                let ovl_e = h1.min(lm.end);
                if ovl_e >= ovl_s {
                    let ovl_len = ovl_e - ovl_s + 1;
                    let lm_len = lm.end.saturating_sub(lm.start) + 1;
                    let min_ovl = (lm_len / 4).min(30).max(1);
                    if ovl_len >= min_ovl && hit_names.last() != Some(&lm.short_name) {
                        hit_names.push(lm.short_name.clone());
                    }
                }
            }
            if !hit_names.is_empty() {
                return hit_names.join("+");
            }
            let mid = (h0 + h1) / 2;
            if let Some(lm) = self.nearest_landmark(mid) {
                return lm.short_name.clone();
            }
            return format!("{}-{}", h0, h1);
        }

        let span = h1.saturating_sub(h0);
        let p0 = (h0 + 15).min(h1);
        let p1 = h1.saturating_sub(15).max(h0);
        let (Some(l0), Some(l1)) = (self.nearest_landmark(p0), self.nearest_landmark(p1)) else {
            return format!("{}-{}", h0, h1);
        };
        if l0.short_name == l1.short_name {
            if let Some(ref parent) = l0.parent_gene {
                if l0.parent_min_span > 0 && span > l0.parent_min_span {
                    return parent.clone();
                }
            }
            return l0.short_name.clone();
        }
        if l0.parent_gene.is_some() && l0.parent_gene == l1.parent_gene {
            if span > l0.parent_min_span.max(1) {
                return l0.parent_gene.clone().unwrap();
            }
        }
        let g0_norm = if l0.parent_gene.is_some() && l0.parent_gene != l1.parent_gene {
            l0.parent_gene.as_deref().unwrap()
        } else {
            l0.short_name.as_str()
        };
        let g1_norm = if l1.parent_gene.is_some() && l1.parent_gene != l0.parent_gene {
            l1.parent_gene.as_deref().unwrap()
        } else {
            l1.short_name.as_str()
        };
        format!("{}..{}", g0_norm, g1_norm)
    }

    pub fn display_lineage_name(&self, raw: &str) -> String {
        if let Some(alias) = self
            .basal_trunks
            .iter()
            .find(|t| t.name == raw)
            .and_then(|t| t.display_alias.as_ref())
        {
            return alias.clone();
        }
        EggmanEngine::display_lineage_name(raw)
    }
}

/// Complete pathogen reference pack containing genome metadata and reference profile alignments.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PathogenProfilePack {
    pub spec: PathogenSpec,
    #[serde(default)]
    pub profiles_tsv: Option<String>,
    #[serde(default)]
    pub profiles_fasta: Option<String>,
    #[serde(default)]
    pub founder_profiles_tsv: Option<String>,
    #[serde(default)]
    pub founder_profiles_fasta: Option<String>,
    #[serde(default)]
    pub certified_profiles_tsv: Option<String>,
    #[serde(default)]
    pub compendium_json: Option<String>,
}

fn resolve_pack_field(
    base_dir: Option<&Path>,
    val: Option<String>,
    default_filenames: &[&str],
) -> Result<Option<String>, String> {
    if let Some(s) = val {
        let trimmed = s.trim();
        let is_inline = trimmed.contains('\n')
            || trimmed.contains('\r')
            || trimmed.contains('\t')
            || trimmed.starts_with('>')
            || trimmed.starts_with('{')
            || trimmed.starts_with('[')
            || trimmed.starts_with('#');
        if is_inline {
            return Ok(Some(s));
        }
        if let Some(dir) = base_dir {
            let joined = dir.join(trimmed);
            if joined.is_file() {
                return fs::read_to_string(&joined)
                    .map(Some)
                    .map_err(|e| format!("failed to read pack asset {:?}: {}", joined, e));
            }
        }
        let direct = Path::new(trimmed);
        if direct.is_file() {
            return fs::read_to_string(direct)
                .map(Some)
                .map_err(|e| format!("failed to read pack asset {:?}: {}", direct, e));
        }
        return Err(format!("referenced pack asset not found: '{}'", trimmed));
    }
    if let Some(dir) = base_dir {
        for &fname in default_filenames {
            let cand = dir.join(fname);
            if cand.is_file() {
                return fs::read_to_string(&cand)
                    .map(Some)
                    .map_err(|e| format!("failed to read pack asset {:?}: {}", cand, e));
            }
        }
    }
    Ok(None)
}

impl PathogenProfilePack {
    /// Embedded HIV-1 Group M pack using built-in HXB2 profiles and certified compendium.
    pub fn hiv1_group_m_default() -> Self {
        Self {
            spec: PathogenSpec::hiv1_group_m_default(),
            profiles_tsv: Some(PROFILES_TSV.to_string()),
            profiles_fasta: None,
            founder_profiles_tsv: Some(FOUNDER_PROFILES_TSV.to_string()),
            founder_profiles_fasta: None,
            certified_profiles_tsv: Some(CERTIFIED_PROFILES_TSV.to_string()),
            compendium_json: Some(COMPENDIUM_JSON.to_string()),
        }
    }

    /// Embedded 8-segment Influenza A (2.3.4.4b H5N1 / LPAI) pack.
    pub fn influenza_a_8seg_default() -> Result<Self, String> {
        let spec: PathogenSpec = serde_json::from_str(IAV_MANIFEST_JSON)
            .map_err(|e| format!("invalid embedded IAV manifest: {}", e))?;
        Ok(Self {
            spec,
            profiles_tsv: None,
            profiles_fasta: Some(IAV_PROFILES_FASTA.to_string()),
            founder_profiles_tsv: None,
            founder_profiles_fasta: None,
            certified_profiles_tsv: Some(IAV_CERTIFIED_TSV.to_string()),
            compendium_json: Some(IAV_COMPENDIUM_JSON.to_string()),
        })
    }

    /// Embedded clinical SARS-CoV-2 (Wuhan-Hu-1) pack.
    pub fn sars_cov2_wuhan_default() -> Result<Self, String> {
        let spec: PathogenSpec = serde_json::from_str(SC2_MANIFEST_JSON)
            .map_err(|e| format!("invalid embedded SARS-CoV-2 manifest: {}", e))?;
        Ok(Self {
            spec,
            profiles_tsv: None,
            profiles_fasta: Some(SC2_PROFILES_FASTA.to_string()),
            founder_profiles_tsv: None,
            founder_profiles_fasta: None,
            certified_profiles_tsv: Some(SC2_CERTIFIED_TSV.to_string()),
            compendium_json: Some(SC2_COMPENDIUM_JSON.to_string()),
        })
    }

    /// Embedded wildlife Sarbecovirus (Wuhan-Hu-1 coordinates) pack.
    pub fn sarbecovirus_wuhan_default() -> Result<Self, String> {
        let spec: PathogenSpec = serde_json::from_str(SARBECO_MANIFEST_JSON)
            .map_err(|e| format!("invalid embedded Sarbecovirus manifest: {}", e))?;
        Ok(Self {
            spec,
            profiles_tsv: None,
            profiles_fasta: Some(SARBECO_PROFILES_FASTA.to_string()),
            founder_profiles_tsv: None,
            founder_profiles_fasta: None,
            certified_profiles_tsv: None,
            compendium_json: None,
        })
    }

    /// Resolve a built-in pathogen pack by short identifier (`hiv1`, `iav`, `sc2`, `sarbeco`).
    pub fn from_builtin_id(pack_id: &str) -> Result<Self, String> {
        match pack_id.trim().to_ascii_lowercase().as_str() {
            "hiv1" | "hiv-1" | "hiv" | "default" => Ok(Self::hiv1_group_m_default()),
            "iav" | "flu" | "h5n1" | "influenza" | "influenza_a_8seg" => {
                Self::influenza_a_8seg_default()
            }
            "sc2" | "sars2" | "sars-cov-2" | "sars_cov2" | "sars_cov2_wuhan" => {
                Self::sars_cov2_wuhan_default()
            }
            "sarbeco" | "sarbecovirus" | "sarbecovirus_wuhan" => {
                Self::sarbecovirus_wuhan_default()
            }
            other => Err(format!("unknown built-in pathogen pack '{}'", other)),
        }
    }

    /// Construct a `PathogenProfilePack` directly from an in-memory JSON bundle (`PathogenProfilePack` or `PathogenSpec`)
    /// plus an optional aligned reference FASTA override (useful for browser/WASM custom pack loading).
    pub fn from_json_and_fasta_str(
        json_str: &str,
        profiles_fasta_override: &str,
    ) -> Result<Self, String> {
        let mut pack: PathogenProfilePack = match serde_json::from_str::<PathogenProfilePack>(json_str) {
            Ok(p) => p,
            Err(pack_err) => match serde_json::from_str::<PathogenSpec>(json_str) {
                Ok(spec) => PathogenProfilePack {
                    spec,
                    profiles_tsv: None,
                    profiles_fasta: None,
                    founder_profiles_tsv: None,
                    founder_profiles_fasta: None,
                    certified_profiles_tsv: None,
                    compendium_json: None,
                },
                Err(spec_err) => {
                    return Err(format!(
                        "invalid JSON bundle: not a PathogenProfilePack ({}) or PathogenSpec ({})",
                        pack_err, spec_err
                    ));
                }
            },
        };
        if !profiles_fasta_override.trim().is_empty() {
            pack.profiles_fasta = Some(profiles_fasta_override.to_string());
        }
        Ok(pack)
    }

    /// Load a pathogen profile pack from either a JSON pack file or a directory containing `manifest.json` / `pack.json`.
    pub fn from_path(path: &Path) -> Result<Self, String> {
        let (manifest_path, base_dir): (std::path::PathBuf, Option<&Path>) = if path.is_dir() {
            let m1 = path.join("manifest.json");
            let m2 = path.join("pack.json");
            if m1.is_file() {
                (m1, Some(path))
            } else if m2.is_file() {
                (m2, Some(path))
            } else {
                return Err(format!(
                    "pack directory {:?} must contain manifest.json or pack.json",
                    path
                ));
            }
        } else if path.is_file() {
            (path.to_path_buf(), path.parent())
        } else {
            return Err(format!("pathogen pack path does not exist: {:?}", path));
        };

        let raw_json = fs::read_to_string(&manifest_path)
            .map_err(|e| format!("failed to read manifest {:?}: {}", manifest_path, e))?;

        let mut pack: PathogenProfilePack = match serde_json::from_str::<PathogenProfilePack>(&raw_json) {
            Ok(p) => p,
            Err(pack_err) => match serde_json::from_str::<PathogenSpec>(&raw_json) {
                Ok(spec) => PathogenProfilePack {
                    spec,
                    profiles_tsv: None,
                    profiles_fasta: None,
                    founder_profiles_tsv: None,
                    founder_profiles_fasta: None,
                    certified_profiles_tsv: None,
                    compendium_json: None,
                },
                Err(spec_err) => {
                    return Err(format!(
                        "invalid JSON in {:?}: not a PathogenProfilePack ({}) or PathogenSpec ({})",
                        manifest_path, pack_err, spec_err
                    ));
                }
            },
        };

        pack.profiles_tsv = resolve_pack_field(base_dir, pack.profiles_tsv.take(), &["profiles.tsv"])?;
        pack.profiles_fasta =
            resolve_pack_field(base_dir, pack.profiles_fasta.take(), &["profiles.fasta", "profiles.fa"])?;
        pack.founder_profiles_tsv = resolve_pack_field(
            base_dir,
            pack.founder_profiles_tsv.take(),
            &["founder_profiles.tsv"],
        )?;
        pack.founder_profiles_fasta = resolve_pack_field(
            base_dir,
            pack.founder_profiles_fasta.take(),
            &["founder_profiles.fasta", "founder_profiles.fa"],
        )?;
        pack.certified_profiles_tsv = resolve_pack_field(
            base_dir,
            pack.certified_profiles_tsv.take(),
            &["certified_profiles.tsv"],
        )?;
        pack.compendium_json =
            resolve_pack_field(base_dir, pack.compendium_json.take(), &["compendium.json", "compendium_bundle.json"])?;

        Ok(pack)
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CompendiumEntry {
    pub label: String,
    pub short_code: String,
    pub category: String,
    pub status: String,
    pub descriptive_name: String,
    pub structural_formula: String,
    pub backbone: String,
    #[serde(default)]
    pub t_mrca: Option<f64>,
    #[serde(default)]
    pub t_mrca_ci_low: Option<f64>,
    #[serde(default)]
    pub t_mrca_ci_high: Option<f64>,
    #[serde(default)]
    pub first_sampled_year: Option<f64>,
    #[serde(default)]
    pub chronaeon_half_life_yr: Option<f64>,
    pub hxb2_barcode: String,
    pub scope: String,
}

#[derive(Debug, Clone, Deserialize)]
struct CompendiumBundle {
    #[serde(default)]
    lanl_msa_len: usize,
    #[serde(default)]
    lanl_msa_to_hxb2: Vec<usize>,
    #[serde(default)]
    compendium: Vec<CompendiumEntry>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ScreeningConfig {
    /// Template-switch penalty gamma (default 7.0)
    pub gamma: Option<f64>,
    /// Minimum segment span in HXB2 nucleotides (default 120 for full, 150 for partial)
    pub min_span_nt: Option<usize>,
    /// Minimum subtype-informative SNPs per segment (default 5)
    pub min_informative_snps: usize,
    /// Fisher exact test p-value threshold (default 1e-4)
    pub p_fisher_thresh: f64,
    /// Mask hypervariable gp120 V1-V5 loops and gag p6 PTAP duplication zone
    pub mask_hypervariable_loops: bool,
}

impl ScreeningConfig {
    /// Reject parameter values that would make the screen meaningless or non-deterministic.
    pub fn validate(&self) -> Result<(), String> {
        if let Some(g) = self.gamma {
            if !g.is_finite() || g <= 0.0 {
                return Err(format!("gamma must be a finite number > 0 (got {})", g));
            }
        }
        if self.min_span_nt == Some(0) {
            return Err("min-span must be >= 1".to_string());
        }
        if !self.p_fisher_thresh.is_finite() || self.p_fisher_thresh <= 0.0 || self.p_fisher_thresh > 1.0 {
            return Err(format!(
                "p-thresh must be in (0, 1] (got {})",
                self.p_fisher_thresh
            ));
        }
        Ok(())
    }
}

impl Default for ScreeningConfig {
    fn default() -> Self {
        Self {
            gamma: None,
            min_span_nt: None,
            min_informative_snps: 5,
            p_fisher_thresh: 1e-4,
            mask_hypervariable_loops: true,
        }
    }
}

/// Optional non-lookahead temporal donor admissibility configuration.
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct TemporalGuardConfig {
    /// Enable temporal masking of Tier-1 trunks and Tier-2 clades emerging after the query timestamp.
    pub enabled: bool,
    /// Explicit decimal query year override; when `None`, extracted from the FASTA header.
    pub query_year: Option<f64>,
    /// Lookahead tolerance in years (`None` uses `PathogenSpec::effective_temporal_slack_yr()`).
    pub slack_yr: Option<f64>,
}

impl TemporalGuardConfig {
    pub fn validate(&self) -> Result<(), String> {
        if let Some(y) = self.query_year {
            if !y.is_finite() || !(1900.0..=2100.0).contains(&y) {
                return Err(format!(
                    "query-year must be a finite year in [1900, 2100] (got {})",
                    y
                ));
            }
        }
        if let Some(s) = self.slack_yr {
            if !s.is_finite() || !(0.0..=50.0).contains(&s) {
                return Err(format!(
                    "temporal-slack-yr must be a finite value in [0, 50] (got {})",
                    s
                ));
            }
        }
        Ok(())
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SegmentCall {
    pub seg_idx: usize,
    pub trunk: String,
    pub hxb2_start: usize,
    pub hxb2_end: usize,
    pub span_nt: usize,
    pub ci_hxb2_start: usize,
    pub ci_hxb2_end: usize,
    pub short_gene: String,
    pub identity_pct: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BreakpointCall {
    pub bp_hxb2: usize,
    pub ci_hxb2_low: usize,
    pub ci_hxb2_high: usize,
    pub ci_width_nt: usize,
    pub left_trunk: String,
    pub right_trunk: String,
    pub p_fisher: f64,
    pub n11: usize,
    pub n12: usize,
    pub n21: usize,
    pub n22: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LineageMatch {
    pub label: String,
    pub descriptive_name: String,
    pub status: String,
    pub match_type: String,
    pub hxb2_barcode: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LineageCallSupport {
    pub identity_pct: f64,
    pub compared_nt: usize,
    pub net_stem_snps: i32,
    pub p_value: f64,
    pub note: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SequenceResult {
    pub sequence_id: String,
    pub raw_length_nt: usize,
    pub mapped_coding_nt: usize,
    pub hxb2_span_start: usize,
    pub hxb2_span_end: usize,
    pub locus_scope: String,
    pub alignment_mode: String,
    pub called_lineage: String,
    pub called_category: String,
    pub call_support: Option<LineageCallSupport>,
    pub epidemiological_formula: String,
    pub epidemiological_barcode: String,
    pub epidemiological_segments: Vec<SegmentCall>,
    pub basal_formula: String,
    pub basal_barcode: String,
    pub structural_formula: String,
    pub descriptive_name: String,
    pub backbone_trunk: String,
    pub constituent_trunks: Vec<String>,
    pub is_clonal: bool,
    pub n_segments: usize,
    pub n_breakpoints: usize,
    pub hxb2_barcode: String,
    pub ci_bracketed_barcode: String,
    pub segments: Vec<SegmentCall>,
    pub breakpoints: Vec<BreakpointCall>,
    pub matches: Vec<LineageMatch>,
    #[serde(default)]
    pub query_year: Option<f64>,
    #[serde(default)]
    pub masked_future_trunks: Vec<String>,
}

/// Slot coverage metadata for a single RNA segment or genomic partition of a segmented isolate.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SegmentSlotInfo {
    pub seg_name: String,
    pub seg_index: usize,
    pub ref_start: usize,
    pub ref_end: usize,
    pub present: bool,
    pub source_header: Option<String>,
    pub raw_len_nt: usize,
    pub mapped_nt: usize,
}

/// Summary of how an isolate's FASTA record(s) map across segmented genome slots.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SegmentedIsolateSummary {
    pub isolate_id: String,
    pub n_segments_present: usize,
    pub total_segments: usize,
    pub is_assembled_from_parts: bool,
    pub total_raw_nt: usize,
    pub total_mapped_nt: usize,
    pub segments: Vec<SegmentSlotInfo>,
}

/// Ingestion summary returned alongside classification results when handling segmented pathogens.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SegmentedAssemblyResult {
    pub mode: String,
    pub raw_record_count: usize,
    pub isolate_count: usize,
    pub was_multi_record_assembled: bool,
    pub isolates: Vec<SegmentedIsolateSummary>,
}

struct CertifiedLineage {
    label: String,
    category: String,
    backbone_idx: usize,
    global_stem: bool,
    clonal_stem: bool,
    dag_parent_indices: Vec<usize>,
    segments: Vec<(String, usize, usize)>,
    segment_trunk_indices: Vec<usize>,
    emergence_floor: Option<f64>,
    profile: Vec<i8>,
}

struct RivalProfile {
    label: String,
    profile: Vec<i8>,
    emergence_floor: Option<f64>,
    trunk_deps: Vec<usize>,
}

struct LineageCallDecision {
    called_lineage: String,
    called_category: String,
    certified_compendium_label: Option<String>,
    call_support: LineageCallSupport,
    residual_tract: Option<(String, String, usize, usize)>,
}

#[deprecated(note = "renamed to EggmanEngine")]
pub type RhizTypeEngine = EggmanEngine;

pub struct EggmanEngine {
    spec: PathogenSpec,
    pure_basal_count: usize,
    segment_id_by_ref_pos: Vec<u16>,
    pub trunk_names: Vec<String>,
    /// Optional emergence floor `min(t_mrca_ci_low, t_mrca, first_sampled_year)` per trunk.
    trunk_floors: Vec<Option<f64>>,
    /// Lookup from trunk name or display alias -> index in `trunk_names`.
    trunk_alias_map: HashMap<String, usize>,
    /// Trunks x coding sites (0=A, 1=C, 2=G, 3=T, -1=indel, -2=forbidden domain)
    trunk_mat: Vec<Vec<i8>>,
    /// Full reference-length profiles for each trunk + anchor ref + key founders for k-mer anchoring and banded alignment
    full_profiles_ref: Vec<Vec<i8>>,
    /// Bitmask over 1..ref_full_len of bases present in any full profile (accelerates banded NW inner loop)
    base_mask_ref: Vec<u8>,
    ghost_donor_coding: Vec<i8>,
    /// 15-mer hash index -> list of (profile_idx, ref_0based_pos)
    kmer_index: HashMap<u32, Vec<(u16, u32)>>,
    lanl_msa_len: usize,
    lanl_msa_to_hxb2: Vec<usize>,
    pub compendium: Vec<CompendiumEntry>,
    loop_mask_coding: Vec<bool>,
    /// Guarded Certified Clade Caller panel
    certified_lineages: Vec<CertifiedLineage>,
    /// Rival profiles for residual recombination checking
    founder_rivals: Vec<RivalProfile>,
    /// Tier-1 Epidemiological segmentation panel
    epi_names: Vec<String>,
    epi_mat: Vec<Vec<i8>>,
    epi_penalties: Vec<f64>,
    epi_floors: Vec<Option<f64>>,
    epi_trunk_deps: Vec<Vec<usize>>,
}

fn base_to_int(b: u8) -> i8 {
    match b.to_ascii_uppercase() {
        b'A' => 0,
        b'C' => 1,
        b'G' => 2,
        b'T' | b'U' => 3,
        _ => -1,
    }
}

fn segment_trunk_at(segments: &[(String, usize, usize)], hxb2: usize) -> Option<&str> {
    segments
        .iter()
        .find(|(_, s, e)| hxb2 >= *s && hxb2 <= *e)
        .map(|(t, _, _)| t.as_str())
}

fn seq_to_int_vec(seq: &str) -> Vec<i8> {
    seq.bytes().map(base_to_int).collect()
}

fn normalize_profile_vec(seq_str: &str, spec: &PathogenSpec, label: &str) -> Result<Vec<i8>, String> {
    let v = seq_to_int_vec(seq_str);
    let cod_len = spec.coding_len();
    if v.len() == spec.ref_full_len {
        Ok(v)
    } else if v.len() == cod_len {
        let mut full = vec![-1i8; spec.ref_full_len];
        full[(spec.coding_start - 1)..spec.coding_end].copy_from_slice(&v);
        Ok(full)
    } else {
        Err(format!(
            "Profile '{}' has length {} (expected ref_full_len={} or coding_len={})",
            label, v.len(), spec.ref_full_len, cod_len
        ))
    }
}

fn parse_tsv_profiles(tsv: &str, spec: &PathogenSpec) -> Result<HashMap<String, Vec<i8>>, String> {
    let mut map = HashMap::new();
    for line in tsv.lines() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let mut parts = line.split('\t');
        if let (Some(name), Some(seq)) = (parts.next(), parts.next()) {
            let v = normalize_profile_vec(seq.trim(), spec, name.trim())?;
            map.insert(name.trim().to_string(), v);
        }
    }
    Ok(map)
}

fn parse_aligned_fasta_profiles(
    fasta_str: &str,
    expected_len: usize,
) -> Result<HashMap<String, Vec<i8>>, String> {
    let mut map = HashMap::new();
    for (hdr, seq) in EggmanEngine::parse_fasta(fasta_str) {
        let name = hdr
            .split_whitespace()
            .next()
            .unwrap_or(hdr.as_str())
            .to_string();
        let v = seq_to_int_vec(&seq);
        if v.len() != expected_len {
            return Err(format!(
                "FASTA profile '{}' has length {} != expected {}",
                name,
                v.len(),
                expected_len
            ));
        }
        map.insert(name, v);
    }
    Ok(map)
}

fn encode_kmer_15(slice: &[i8]) -> Option<u32> {
    if slice.len() != 15 {
        return None;
    }
    let mut val: u32 = 0;
    for &b in slice {
        if b < 0 || b > 3 {
            return None;
        }
        val = (val << 2) | (b as u32);
    }
    Some(val)
}

impl Default for EggmanEngine {
    fn default() -> Self {
        Self::new()
    }
}

impl EggmanEngine {
    pub fn new() -> Self {
        Self::from_pathogen_pack(PathogenProfilePack::hiv1_group_m_default())
            .expect("valid embedded HIV-1 Group M pack")
    }

    pub fn with_certified_profiles_tsv(certified_profiles_tsv: &str) -> Self {
        let mut pack = PathogenProfilePack::hiv1_group_m_default();
        pack.certified_profiles_tsv = Some(certified_profiles_tsv.to_string());
        Self::from_pathogen_pack(pack).expect("valid HIV-1 Group M pack with custom certified profiles")
    }

    /// Construct a `EggmanEngine` from an embedded built-in pack id (`hiv1`, `iav`, `sc2`, `sarbeco`).
    pub fn from_builtin_id(pack_id: &str) -> Result<Self, String> {
        let pack = PathogenProfilePack::from_builtin_id(pack_id)?;
        Self::from_pathogen_pack(pack)
    }

    /// Construct a `EggmanEngine` from a pack directory or JSON manifest path.
    pub fn from_pack_path(path: &Path) -> Result<Self, String> {
        let pack = PathogenProfilePack::from_path(path)?;
        Self::from_pathogen_pack(pack)
    }

    /// Construct a `EggmanEngine` from an in-memory `PathogenProfilePack`.
    pub fn from_pathogen_pack(pack: PathogenProfilePack) -> Result<Self, String> {
        let spec = pack.spec;
        if spec.ref_full_len == 0 {
            return Err("PathogenSpec.ref_full_len must be > 0".to_string());
        }
        if spec.coding_start == 0
            || spec.coding_end < spec.coding_start
            || spec.coding_end > spec.ref_full_len
        {
            return Err(format!(
                "Invalid coding span {}..={} for ref_full_len={}",
                spec.coding_start, spec.coding_end, spec.ref_full_len
            ));
        }
        if spec.basal_trunks.is_empty() {
            return Err("PathogenSpec.basal_trunks cannot be empty".to_string());
        }

        let coding_len = spec.coding_len();
        let coding_start = spec.coding_start;
        let coding_end = spec.coding_end;

        let mut raw_profiles: HashMap<String, Vec<i8>> = HashMap::new();
        if let Some(ref tsv) = pack.profiles_tsv {
            for (k, v) in parse_tsv_profiles(tsv, &spec)? {
                raw_profiles.insert(k, v);
            }
        }
        if let Some(ref fa) = pack.profiles_fasta {
            let fa_map = match parse_aligned_fasta_profiles(fa, spec.ref_full_len) {
                Ok(m) => m,
                Err(_) if coding_len < spec.ref_full_len => {
                    let cod_map = parse_aligned_fasta_profiles(fa, coding_len)?;
                    cod_map
                        .into_iter()
                        .map(|(k, v)| {
                            let mut full = vec![-1i8; spec.ref_full_len];
                            full[(coding_start - 1)..coding_end].copy_from_slice(&v);
                            (k, full)
                        })
                        .collect()
                }
                Err(e) => return Err(e),
            };
            for (k, v) in fa_map {
                raw_profiles.insert(k, v);
            }
        }
        if raw_profiles.is_empty() {
            return Err("No basal trunk profiles provided (profiles_tsv or profiles_fasta required)".to_string());
        }

        let mut founder_profiles: HashMap<String, Vec<i8>> = HashMap::new();
        if let Some(ref f_tsv) = pack.founder_profiles_tsv {
            for (k, v) in parse_tsv_profiles(f_tsv, &spec)? {
                founder_profiles.insert(k, v);
            }
        }
        if let Some(ref f_fa) = pack.founder_profiles_fasta {
            let fa_map = match parse_aligned_fasta_profiles(f_fa, spec.ref_full_len) {
                Ok(m) => m,
                Err(_) if coding_len < spec.ref_full_len => {
                    let cod_map = parse_aligned_fasta_profiles(f_fa, coding_len)?;
                    cod_map
                        .into_iter()
                        .map(|(k, v)| {
                            let mut full = vec![-1i8; spec.ref_full_len];
                            full[(coding_start - 1)..coding_end].copy_from_slice(&v);
                            (k, full)
                        })
                        .collect()
                }
                Err(e) => return Err(e),
            };
            for (k, v) in fa_map {
                founder_profiles.insert(k, v);
            }
        }

        // Order pure basal trunks before domain-restricted non-pure trunks
        let mut ordered_trunks: Vec<&TrunkSpec> = spec
            .basal_trunks
            .iter()
            .filter(|t| t.is_pure_basal)
            .collect();
        let pure_basal_count = ordered_trunks.len();
        if pure_basal_count == 0 {
            return Err("At least one pure basal trunk (is_pure_basal=true) is required".to_string());
        }
        ordered_trunks.extend(spec.basal_trunks.iter().filter(|t| !t.is_pure_basal));

        let mut trunk_names = Vec::with_capacity(ordered_trunks.len());
        let mut trunk_floors = Vec::with_capacity(ordered_trunks.len());
        let mut trunk_mat = Vec::with_capacity(ordered_trunks.len());
        let mut full_profiles_ref = Vec::new();

        if let Some(ref_full) = raw_profiles.get(&spec.reference_name) {
            full_profiles_ref.push(ref_full.clone());
        }

        for t_spec in &ordered_trunks {
            let full_arr = raw_profiles.get(&t_spec.name).ok_or_else(|| {
                format!("Missing trunk profile '{}' in profiles table", t_spec.name)
            })?;
            full_profiles_ref.push(full_arr.clone());
            let mut cod_arr = full_arr[(coding_start - 1)..coding_end].to_vec();
            if let Some((act_lo, act_hi)) = t_spec.active_range {
                for (idx, cell) in cod_arr.iter_mut().enumerate() {
                    let ref_pos = coding_start + idx;
                    if ref_pos < act_lo || ref_pos > act_hi {
                        *cell = -2;
                    }
                }
            }
            trunk_names.push(t_spec.name.clone());
            trunk_floors.push(t_spec.emergence_floor());
            trunk_mat.push(cod_arr);
        }

        // Index optional divergent founder lineages alongside basal trunks for higher alignment sensitivity
        for fk in &spec.align_founder_keys {
            if let Some(fprof) = founder_profiles.get(fk) {
                full_profiles_ref.push(fprof.clone());
            }
        }

        let mut base_mask_ref = vec![0u8; spec.ref_full_len];
        for prof in &full_profiles_ref {
            for (pos, &b) in prof.iter().enumerate() {
                if (0..=3).contains(&b) {
                    base_mask_ref[pos] |= 1u8 << (b as u8);
                }
            }
        }

        let ghost_donor_coding = if let Some(ref rule) = spec.ghost_relabel_rule {
            if let Some(full_p) = raw_profiles
                .get(&rule.donor_profile)
                .or_else(|| founder_profiles.get(&rule.donor_profile))
            {
                full_p[(coding_start - 1)..coding_end].to_vec()
            } else {
                return Err(format!(
                    "Missing ghost donor profile '{}' in profiles table",
                    rule.donor_profile
                ));
            }
        } else {
            Vec::new()
        };

        // Build 15-mer index across all full profiles (step 2 nt for compact fast lookup)
        let mut kmer_index: HashMap<u32, Vec<(u16, u32)>> = HashMap::new();
        for (p_idx, prof) in full_profiles_ref.iter().enumerate() {
            for pos in (0..(prof.len().saturating_sub(15))).step_by(2) {
                if let Some(code) = encode_kmer_15(&prof[pos..(pos + 15)]) {
                    let entry = kmer_index.entry(code).or_default();
                    let pos32 = pos as u32;
                    if entry.len() < 8 && !entry.iter().any(|&(_, rp)| rp == pos32) {
                        entry.push((p_idx as u16, pos32));
                    }
                }
            }
        }

        let bundle: CompendiumBundle = if let Some(ref comp_json) = pack.compendium_json {
            match serde_json::from_str::<CompendiumBundle>(comp_json) {
                Ok(b) => b,
                Err(_) => {
                    let entries: Vec<CompendiumEntry> = serde_json::from_str(comp_json)
                        .map_err(|e| format!("Invalid compendium_json: {}", e))?;
                    CompendiumBundle {
                        lanl_msa_len: 0,
                        lanl_msa_to_hxb2: Vec::new(),
                        compendium: entries,
                    }
                }
            }
        } else {
            CompendiumBundle {
                lanl_msa_len: 0,
                lanl_msa_to_hxb2: Vec::new(),
                compendium: Vec::new(),
            }
        };

        let comp_map: HashMap<&str, &CompendiumEntry> = bundle
            .compendium
            .iter()
            .map(|e| (e.label.as_str(), e))
            .collect();
        let mut trunk_alias_map: HashMap<String, usize> = HashMap::new();
        for (i, t_spec) in ordered_trunks.iter().enumerate() {
            trunk_alias_map.insert(t_spec.name.clone(), i);
            if let Some(ref alias) = t_spec.display_alias {
                trunk_alias_map.entry(alias.clone()).or_insert(i);
            }
        }

        // Load Tier-2 Certified Lineages
        let mut certified_lineages = Vec::new();
        let mut raw_parent_names: Vec<Vec<String>> = Vec::new();
        if let Some(ref cert_tsv) = pack.certified_profiles_tsv {
            for line in cert_tsv.lines() {
                let line = line.trim();
                if line.is_empty() || line.starts_with('#') {
                    continue;
                }
                let cols: Vec<&str> = line.split('\t').collect();
                if cols.len() < 6 {
                    continue;
                }
                let label = cols[0].to_string();
                let category = cols[1].to_string();
                let flags = cols[2];
                let parents_col = cols[3];
                let seq_cod = seq_to_int_vec(cols[5].trim());
                if seq_cod.len() != coding_len {
                    return Err(format!(
                        "Certified profile '{}' has length {} != expected coding_len {}",
                        label,
                        seq_cod.len(),
                        coding_len
                    ));
                }

                let entry_opt = comp_map.get(label.as_str()).copied();
                let (backbone_idx, segments, emergence_floor) = if let Some(entry) = entry_opt {
                    let b_idx = *trunk_alias_map.get(&entry.backbone).unwrap_or(&0);
                    let floor = min_temporal_floor(
                        entry.t_mrca_ci_low,
                        entry.t_mrca,
                        entry.first_sampled_year,
                    );
                    (b_idx, parse_barcode_segments(&entry.hxb2_barcode), floor)
                } else if !bundle.compendium.is_empty() {
                    return Err(format!(
                        "Certified lineage '{}' missing from compendium_json",
                        label
                    ));
                } else {
                    (0usize, Vec::new(), None)
                };

                let mut segment_trunk_indices = vec![backbone_idx];
                for (seg_trk, _, _) in &segments {
                    if let Some(&ti) = trunk_alias_map.get(seg_trk.as_str()) {
                        if !segment_trunk_indices.contains(&ti) {
                            segment_trunk_indices.push(ti);
                        }
                    }
                }

                let global_stem = flags.split(',').any(|f| f == "global_stem");
                let clonal_stem = flags.split(',').any(|f| f == "clonal_stem");
                let p_names = if parents_col == "-" || parents_col.is_empty() {
                    Vec::new()
                } else {
                    parents_col.split('+').map(|s| s.to_string()).collect()
                };
                raw_parent_names.push(p_names);

                certified_lineages.push(CertifiedLineage {
                    label,
                    category,
                    backbone_idx,
                    global_stem,
                    clonal_stem,
                    dag_parent_indices: Vec::new(),
                    segments,
                    segment_trunk_indices,
                    emergence_floor,
                    profile: seq_cod,
                });
            }
        }

        let cert_idx_map: HashMap<String, usize> = certified_lineages
            .iter()
            .enumerate()
            .map(|(i, c)| (c.label.clone(), i))
            .collect();
        for (i, p_names) in raw_parent_names.into_iter().enumerate() {
            for pn in p_names {
                if let Some(&pi) = cert_idx_map.get(&pn) {
                    certified_lineages[i].dag_parent_indices.push(pi);
                }
            }
        }

        let founder_meta = |fk: &str| -> (Option<f64>, Vec<usize>) {
            if let Some(&ci) = cert_idx_map.get(fk) {
                let c = &certified_lineages[ci];
                return (c.emergence_floor, c.segment_trunk_indices.clone());
            }
            if let Some(entry) = comp_map.get(fk) {
                let mut deps = Vec::new();
                if let Some(&bi) = trunk_alias_map.get(&entry.backbone) {
                    deps.push(bi);
                }
                for (seg_trk, _, _) in parse_barcode_segments(&entry.hxb2_barcode) {
                    if let Some(&ti) = trunk_alias_map.get(seg_trk.as_str()) {
                        if !deps.contains(&ti) {
                            deps.push(ti);
                        }
                    }
                }
                let floor = min_temporal_floor(
                    entry.t_mrca_ci_low,
                    entry.t_mrca,
                    entry.first_sampled_year,
                );
                return (floor, deps);
            }
            if let Some(&ti) = trunk_alias_map.get(fk) {
                return (trunk_floors[ti], vec![ti]);
            }
            (None, Vec::new())
        };

        // Build rival list for Guard 4 (basal trunks + founder lineages)
        let mut founder_rivals = Vec::new();
        for (i, (t_name, t_prof)) in trunk_names.iter().zip(trunk_mat.iter()).enumerate() {
            founder_rivals.push(RivalProfile {
                label: t_name.clone(),
                profile: t_prof.clone(),
                emergence_floor: trunk_floors[i],
                trunk_deps: vec![i],
            });
        }
        for (fk, _) in &spec.epi_founder_order {
            if let Some(fprof) = founder_profiles.get(fk) {
                let (emergence_floor, trunk_deps) = founder_meta(fk);
                founder_rivals.push(RivalProfile {
                    label: fk.clone(),
                    profile: fprof[(coding_start - 1)..coding_end].to_vec(),
                    emergence_floor,
                    trunk_deps,
                });
            }
        }

        // Build Tier-1 Epidemiological Segmentation Panel
        let mut epi_names = Vec::new();
        let mut epi_mat = Vec::new();
        let mut epi_penalties = Vec::new();
        let mut epi_floors = Vec::new();
        let mut epi_trunk_deps = Vec::new();
        let mut replaced_founders: HashSet<String> = HashSet::new();
        for (i, t_name) in trunk_names.iter().enumerate().take(pure_basal_count) {
            if let Some((_, cert_key, disp_name)) = spec
                .epi_certified_replacements
                .iter()
                .find(|(tr, _, _)| tr == t_name)
            {
                epi_names.push(disp_name.clone());
                let rep_prof = founder_profiles
                    .get(cert_key)
                    .map(|p| p[(coding_start - 1)..coding_end].to_vec())
                    .unwrap_or_else(|| trunk_mat[i].clone());
                epi_mat.push(rep_prof);
                let (c_floor, mut c_deps) = founder_meta(cert_key);
                if !c_deps.contains(&i) {
                    c_deps.push(i);
                }
                epi_floors.push(c_floor.or(trunk_floors[i]));
                epi_trunk_deps.push(c_deps);
                replaced_founders.insert(cert_key.clone());
            } else {
                epi_names.push(t_name.clone());
                epi_mat.push(trunk_mat[i].clone());
                epi_floors.push(trunk_floors[i]);
                epi_trunk_deps.push(vec![i]);
            }
            epi_penalties.push(0.0);
        }
        for (fk, disp_name) in &spec.epi_founder_order {
            if replaced_founders.contains(fk) {
                continue;
            }
            if let Some(fprof) = founder_profiles.get(fk) {
                epi_names.push(disp_name.clone());
                epi_mat.push(fprof[(coding_start - 1)..coding_end].to_vec());
                epi_penalties.push(spec.epi_founder_penalty);
                let (emergence_floor, trunk_deps) = founder_meta(fk);
                epi_floors.push(emergence_floor);
                epi_trunk_deps.push(trunk_deps);
            }
        }

        let mut loop_mask_coding = vec![false; coding_len];
        for &(lo, hi) in &spec.hypervariable_masks {
            for h in lo..=hi {
                if h >= coding_start && h <= coding_end {
                    loop_mask_coding[h - coding_start] = true;
                }
            }
        }

        let mut segment_id_by_ref_pos = vec![0u16; spec.ref_full_len];
        if spec.is_segmented {
            for (seg_idx, lm) in spec.landmarks.iter().enumerate() {
                let s = lm.start.max(1);
                let e = lm.end.min(spec.ref_full_len);
                for p in s..=e {
                    segment_id_by_ref_pos[p - 1] = (seg_idx + 1) as u16;
                }
            }
        }

        Ok(Self {
            spec,
            pure_basal_count,
            segment_id_by_ref_pos,
            trunk_names,
            trunk_floors,
            trunk_alias_map,
            trunk_mat,
            full_profiles_ref,
            base_mask_ref,
            ghost_donor_coding,
            kmer_index,
            lanl_msa_len: bundle.lanl_msa_len,
            lanl_msa_to_hxb2: bundle.lanl_msa_to_hxb2,
            compendium: bundle.compendium,
            loop_mask_coding,
            certified_lineages,
            founder_rivals,
            epi_names,
            epi_mat,
            epi_penalties,
            epi_floors,
            epi_trunk_deps,
        })
    }

    pub fn spec(&self) -> &PathogenSpec {
        &self.spec
    }

    pub fn pathogen_name(&self) -> &str {
        &self.spec.pathogen_name
    }

    pub fn reference_name(&self) -> &str {
        &self.spec.reference_name
    }

    pub fn basal_trunk_count(&self) -> usize {
        self.trunk_names.len()
    }

    pub fn certified_lineage_count(&self) -> usize {
        self.certified_lineages.len()
    }

    #[inline]
    fn segment_id_at_ref(&self, ref_pos_1based: usize) -> u16 {
        if !self.spec.is_segmented || ref_pos_1based == 0 {
            return 0;
        }
        self.segment_id_by_ref_pos
            .get(ref_pos_1based - 1)
            .copied()
            .unwrap_or(0)
    }

    fn effective_segment_min_span(&self, h0: usize, h1: usize, default_min: usize) -> usize {
        if !self.spec.is_segmented {
            return default_min;
        }
        let min_lm_len = self
            .spec
            .landmarks
            .iter()
            .filter(|lm| h0 <= lm.end && h1 >= lm.start)
            .map(|lm| lm.end.saturating_sub(lm.start) + 1)
            .min();
        if let Some(lm_len) = min_lm_len {
            default_min.min((lm_len / 2).max(1))
        } else {
            default_min
        }
    }

    /// Parse a FASTA string into `(header, sequence)` tuples.
    pub fn parse_fasta(fasta_text: &str) -> Vec<(String, String)> {
        let mut records = Vec::new();
        let mut cur_hdr: Option<String> = None;
        let mut cur_seq = String::new();

        for line in fasta_text.trim_start_matches('\u{FEFF}').split(['\n', '\r']) {
            let line = line.trim();
            if line.is_empty() {
                continue;
            }
            if let Some(rest) = line.strip_prefix('>') {
                if let Some(h) = cur_hdr.take() {
                    records.push((h, cur_seq.clone()));
                }
                cur_hdr = Some(rest.trim().to_string());
                cur_seq.clear();
            } else {
                if cur_hdr.is_none() {
                    cur_hdr = Some("unnamed_1".to_string());
                }
                for ch in line.chars() {
                    if !ch.is_whitespace() {
                        cur_seq.push(ch.to_ascii_uppercase());
                    }
                }
            }
        }
        if let Some(h) = cur_hdr {
            records.push((h, cur_seq));
        }
        records
    }

    /// Return `(matched_15mers, total_valid_15mers, best_profile_exact_15mers)` on the stronger strand.
    pub fn kmer_hit_stats(&self, raw_seq: &str) -> (usize, usize, usize) {
        let q: Vec<i8> = raw_seq
            .bytes()
            .filter(|&b| b != b'-' && b != b'.' && b != b'|')
            .map(base_to_int)
            .collect();
        if q.len() < 15 {
            return (0, 0, 0);
        }
        let n_prof = self.full_profiles_ref.len();
        let eval_strand = |seq: &[i8]| -> (usize, usize, usize) {
            let mut matched = 0usize;
            let mut total_valid = 0usize;
            let mut prof_hits = vec![0usize; n_prof];
            for pos in 0..=(seq.len() - 15) {
                if let Some(code) = encode_kmer_15(&seq[pos..(pos + 15)]) {
                    total_valid += 1;
                    if let Some(locs) = self.kmer_index.get(&code) {
                        matched += 1;
                        let mut seen_mask_p: u64 = 0;
                        for &(p_idx, _) in locs {
                            let pi = p_idx as usize;
                            if pi < n_prof {
                                if pi < 64 {
                                    let bit = 1u64 << pi;
                                    if seen_mask_p & bit == 0 {
                                        seen_mask_p |= bit;
                                        prof_hits[pi] += 1;
                                    }
                                } else {
                                    prof_hits[pi] += 1;
                                }
                            }
                        }
                    }
                }
            }
            let max_prof = prof_hits.into_iter().max().unwrap_or(0);
            (matched, total_valid, max_prof)
        };

        let fwd = eval_strand(&q);
        // Fast check on reverse complement if forward match rate < 15%
        if fwd.1 > 0 && (fwd.0 as f64) / (fwd.1 as f64) >= 0.15 {
            return fwd;
        }
        let rc: Vec<i8> = q
            .iter()
            .rev()
            .map(|&b| if (0..=3).contains(&b) { 3 - b } else { -1 })
            .collect();
        let rev = eval_strand(&rc);
        if rev.0 > fwd.0 { rev } else { fwd }
    }

    /// Count 15-mer hits per segmented landmark index (`0..landmarks.len()`) on the stronger strand.
    fn identify_segmented_landmark(&self, header: &str, clean_seq: &str) -> Option<usize> {
        let n_lm = self.spec.landmarks.len();
        if !self.spec.is_segmented || n_lm == 0 {
            return None;
        }
        let q: Vec<i8> = clean_seq.bytes().map(base_to_int).collect();
        if q.len() >= 15 {
            let count_strand = |seq: &[i8]| -> Vec<usize> {
                let mut counts = vec![0usize; n_lm];
                for pos in 0..=(seq.len() - 15) {
                    if let Some(code) = encode_kmer_15(&seq[pos..(pos + 15)]) {
                        if let Some(locs) = self.kmer_index.get(&code) {
                            let mut seen_seg: u32 = 0;
                            for &(_, r_pos) in locs {
                                let sid = self.segment_id_at_ref((r_pos as usize) + 1);
                                if sid > 0 {
                                    let idx = (sid - 1) as usize;
                                    if idx < n_lm && (seen_seg & (1u32 << idx.min(31))) == 0 {
                                        seen_seg |= 1u32 << idx.min(31);
                                        counts[idx] += 1;
                                    }
                                }
                            }
                        }
                    }
                }
                counts
            };
            let fwd_counts = count_strand(&q);
            let fwd_total: usize = fwd_counts.iter().sum();
            let counts = if fwd_total >= 4 {
                fwd_counts
            } else {
                let rc: Vec<i8> = q
                    .iter()
                    .rev()
                    .map(|&b| if (0..=3).contains(&b) { 3 - b } else { -1 })
                    .collect();
                let rev_counts = count_strand(&rc);
                let rev_total: usize = rev_counts.iter().sum();
                if rev_total > fwd_total { rev_counts } else { fwd_counts }
            };
            if let Some((best_idx, &best_hits)) = counts
                .iter()
                .enumerate()
                .max_by_key(|&(_, &c)| c)
            {
                if best_hits >= 2 {
                    return Some(best_idx);
                }
            }
        }
        // Fallback: unambiguous header segment keywords
        infer_flu_segment_from_header(header, &self.spec.landmarks)
    }

    /// Screen FASTA records with automatic multi-record assembly for segmented genomes (`Influenza A`).
    /// `seg_mode`: `"auto"`, `"single_isolate"`, or `"per_record"`.
    pub fn classify_fasta_with_segmentation(
        &self,
        fasta_text: &str,
        seg_mode: &str,
        cfg: &ScreeningConfig,
    ) -> (SegmentedAssemblyResult, Vec<SequenceResult>, String) {
        self.classify_fasta_with_segmentation_temporal(
            fasta_text,
            seg_mode,
            cfg,
            &TemporalGuardConfig::default(),
        )
    }

    /// Segmented FASTA ingestion and classification with optional non-lookahead temporal guarding.
    pub fn classify_fasta_with_segmentation_temporal(
        &self,
        fasta_text: &str,
        seg_mode: &str,
        cfg: &ScreeningConfig,
        temporal: &TemporalGuardConfig,
    ) -> (SegmentedAssemblyResult, Vec<SequenceResult>, String) {
        let raw_records: Vec<(String, String)> = Self::parse_fasta(fasta_text)
            .into_iter()
            .filter_map(|(hdr, seq)| {
                let clean: String = seq.chars().filter(|&c| c != '|').collect();
                if clean.is_empty() {
                    None
                } else {
                    Some((hdr, clean))
                }
            })
            .collect();

        let raw_record_count = raw_records.len();
        let n_lm = self.spec.landmarks.len();
        let effective_mode = if !self.spec.is_segmented {
            "per_record"
        } else {
            match seg_mode.trim().to_ascii_lowercase().as_str() {
                "per_record" | "single_isolate" => seg_mode.trim(),
                _ => "auto",
            }
        };

        // Determine if there are unconcatenated sub-genome segment records to assemble
        let has_short_segments = self.spec.is_segmented
            && effective_mode != "per_record"
            && raw_records.iter().any(|(_, s)| s.len() < 8_000);

        if !has_short_segments {
            let mut results = Vec::with_capacity(raw_records.len());
            let mut isolates = Vec::with_capacity(raw_records.len());
            let mut out_fasta = String::new();
            for (hdr, seq) in &raw_records {
                let res = self.classify_sequence_with_temporal(hdr, seq, cfg, temporal);
                let cov_vec = self.compute_landmark_coverage_vec(seq);
                let mut slots = Vec::with_capacity(n_lm);
                let mut present_cnt = 0usize;
                for (si, lm) in self.spec.landmarks.iter().enumerate() {
                    let mapped_nt = cov_vec.get(si).copied().unwrap_or(0);
                    let present = mapped_nt >= 50;
                    if present {
                        present_cnt += 1;
                    }
                    slots.push(SegmentSlotInfo {
                        seg_name: lm.short_name.clone(),
                        seg_index: si + 1,
                        ref_start: lm.start,
                        ref_end: lm.end,
                        present,
                        source_header: if present { Some(hdr.clone()) } else { None },
                        raw_len_nt: if present { mapped_nt } else { 0 },
                        mapped_nt,
                    });
                }
                isolates.push(SegmentedIsolateSummary {
                    isolate_id: hdr.clone(),
                    n_segments_present: present_cnt,
                    total_segments: n_lm,
                    is_assembled_from_parts: false,
                    total_raw_nt: res.raw_length_nt,
                    total_mapped_nt: res.mapped_coding_nt,
                    segments: slots,
                });
                results.push(res);
                out_fasta.push('>');
                out_fasta.push_str(hdr);
                out_fasta.push('\n');
                out_fasta.push_str(seq);
                out_fasta.push('\n');
            }
            let assembly = SegmentedAssemblyResult {
                mode: effective_mode.to_string(),
                raw_record_count,
                isolate_count: isolates.len(),
                was_multi_record_assembled: false,
                isolates,
            };
            return (assembly, results, out_fasta);
        }

        // Partition long pre-concatenated records from short segment records
        struct PartRec {
            hdr: String,
            seq: String,
            seg_idx: Option<usize>,
            strain_key: Option<String>,
        }

        let mut full_items: Vec<(String, String)> = Vec::new();
        let mut short_parts: Vec<PartRec> = Vec::new();
        for (hdr, seq) in raw_records {
            if seq.len() >= 8_000 {
                full_items.push((hdr, seq));
            } else {
                let seg_idx = self.identify_segmented_landmark(&hdr, &seq);
                let strain_key = extract_flu_isolate_key(&hdr);
                short_parts.push(PartRec {
                    hdr,
                    seq,
                    seg_idx,
                    strain_key,
                });
            }
        }

        // Group short segment records into candidate isolates
        struct AssembledGroup {
            isolate_id: String,
            parts: Vec<(usize, String, String)>, // (seg_idx, src_hdr, seq)
        }
        let mut groups: Vec<AssembledGroup> = Vec::new();
        let mut passthrough_short: Vec<(String, String)> = Vec::new();

        if effective_mode == "single_isolate" {
            let first_id = short_parts
                .iter()
                .find_map(|p| p.strain_key.clone())
                .or_else(|| {
                    short_parts.first().and_then(|p| {
                        p.hdr
                            .split_whitespace()
                            .next()
                            .map(|tok| format!("{}_assembled", tok))
                    })
                })
                .unwrap_or_else(|| "IAV_Assembled_Isolate".to_string());
            let mut grp = AssembledGroup {
                isolate_id: first_id,
                parts: Vec::new(),
            };
            for p in short_parts {
                if let Some(si) = p.seg_idx {
                    if let Some(existing) = grp.parts.iter_mut().find(|x| x.0 == si) {
                        if p.seq.len() > existing.2.len() {
                            *existing = (si, p.hdr, p.seq);
                        }
                    } else {
                        grp.parts.push((si, p.hdr, p.seq));
                    }
                } else {
                    passthrough_short.push((p.hdr, p.seq));
                }
            }
            if !grp.parts.is_empty() {
                groups.push(grp);
            }
        } else {
            // "auto" mode: group by normalized strain key first
            let mut keyed_idx: HashMap<String, usize> = HashMap::new();
            let mut unkeyed: Vec<PartRec> = Vec::new();

            for p in short_parts {
                match (&p.strain_key, p.seg_idx) {
                    (Some(k), Some(si)) => {
                        if let Some(&g_idx) = keyed_idx.get(k) {
                            if let Some(existing) = groups[g_idx].parts.iter_mut().find(|x| x.0 == si) {
                                // Collision on same strain key + same segment -> keep distinct suffix group or longer
                                if p.seq.len() > existing.2.len() {
                                    *existing = (si, p.hdr, p.seq);
                                }
                            } else {
                                groups[g_idx].parts.push((si, p.hdr, p.seq));
                            }
                        } else {
                            let g_idx = groups.len();
                            keyed_idx.insert(k.clone(), g_idx);
                            groups.push(AssembledGroup {
                                isolate_id: k.clone(),
                                parts: vec![(si, p.hdr, p.seq)],
                            });
                        }
                    }
                    _ => unkeyed.push(p),
                }
            }

            // If there are <= n_lm unkeyed records (e.g., 8 separate bare GenBank accessions)
            // whose identified segments are all mutually distinct, auto-merge them into one isolate.
            let assigned_unkeyed: Vec<usize> =
                unkeyed.iter().filter_map(|p| p.seg_idx).collect();
            let mut seen_set: HashSet<usize> = HashSet::new();
            let all_distinct = !assigned_unkeyed.is_empty()
                && assigned_unkeyed.len() == unkeyed.len()
                && unkeyed.len() <= n_lm
                && assigned_unkeyed.iter().all(|&s| seen_set.insert(s));

            if all_distinct && unkeyed.len() >= 2 {
                let first_tok = unkeyed[0]
                    .hdr
                    .split_whitespace()
                    .next()
                    .unwrap_or("IAV_Isolate");
                let mut grp = AssembledGroup {
                    isolate_id: format!("{}_assembled_{}seg", first_tok, unkeyed.len()),
                    parts: Vec::new(),
                };
                for p in unkeyed {
                    if let Some(si) = p.seg_idx {
                        grp.parts.push((si, p.hdr, p.seq));
                    }
                }
                groups.push(grp);
            } else {
                for p in unkeyed {
                    passthrough_short.push((p.hdr, p.seq));
                }
            }
        }

        let mut results = Vec::new();
        let mut isolates = Vec::new();
        let mut out_fasta = String::new();
        let mut was_multi_record_assembled = false;

        // 1. Pre-concatenated full records
        for (hdr, seq) in full_items {
            let res = self.classify_sequence_with_temporal(&hdr, &seq, cfg, temporal);
            let cov_vec = self.compute_landmark_coverage_vec(&seq);
            let mut slots = Vec::with_capacity(n_lm);
            let mut present_cnt = 0usize;
            for (si, lm) in self.spec.landmarks.iter().enumerate() {
                let mapped_nt = cov_vec.get(si).copied().unwrap_or(0);
                let present = mapped_nt >= 50;
                if present {
                    present_cnt += 1;
                }
                slots.push(SegmentSlotInfo {
                    seg_name: lm.short_name.clone(),
                    seg_index: si + 1,
                    ref_start: lm.start,
                    ref_end: lm.end,
                    present,
                    source_header: if present { Some(hdr.clone()) } else { None },
                    raw_len_nt: if present { mapped_nt } else { 0 },
                    mapped_nt,
                });
            }
            isolates.push(SegmentedIsolateSummary {
                isolate_id: hdr.clone(),
                n_segments_present: present_cnt,
                total_segments: n_lm,
                is_assembled_from_parts: false,
                total_raw_nt: res.raw_length_nt,
                total_mapped_nt: res.mapped_coding_nt,
                segments: slots,
            });
            results.push(res);
            out_fasta.push('>');
            out_fasta.push_str(&hdr);
            out_fasta.push('\n');
            out_fasta.push_str(&seq);
            out_fasta.push('\n');
        }

        // 2. Assembled segmented groups
        for grp in groups {
            if grp.parts.len() > 1 {
                was_multi_record_assembled = true;
            }
            let mut scaffold = vec![-1i8; self.spec.ref_full_len];
            let mut total_raw_nt = 0usize;
            let mut slot_meta: Vec<Option<(String, usize, usize)>> = vec![None; n_lm];

            for (si, src_hdr, seg_seq) in &grp.parts {
                let lm = &self.spec.landmarks[*si];
                let raw_nt = seg_seq.bytes().filter(|&b| b != b'-' && b != b'.').count();
                total_raw_nt += raw_nt;
                let (cod, _, _) = self.project_to_hxb2_coding(seg_seq);
                let lo = lm.start.saturating_sub(self.spec.coding_start);
                let hi = lm
                    .end
                    .min(self.spec.coding_end)
                    .saturating_sub(self.spec.coding_start);
                let mut mapped_here = 0usize;
                for c_idx in lo..=hi {
                    if let Some(&b) = cod.get(c_idx) {
                        if b >= 0 {
                            let ref_idx = (self.spec.coding_start - 1) + c_idx;
                            if ref_idx < scaffold.len() {
                                scaffold[ref_idx] = b;
                                mapped_here += 1;
                            }
                        }
                    }
                }
                slot_meta[*si] = Some((src_hdr.clone(), raw_nt, mapped_here));
            }

            let stitched_seq: String = scaffold
                .iter()
                .map(|&b| match b {
                    0 => 'A',
                    1 => 'C',
                    2 => 'G',
                    3 => 'T',
                    _ => 'N',
                })
                .collect();

            let eff_temporal = if temporal.enabled && temporal.query_year.is_none() {
                let part_year = parse_header_decimal_year(&grp.isolate_id).or_else(|| {
                    grp.parts
                        .iter()
                        .filter_map(|(_, h, _)| parse_header_decimal_year(h))
                        .next()
                });
                TemporalGuardConfig {
                    enabled: true,
                    query_year: part_year,
                    slack_yr: temporal.slack_yr,
                }
            } else {
                temporal.clone()
            };
            let mut res =
                self.classify_sequence_with_temporal(&grp.isolate_id, &stitched_seq, cfg, &eff_temporal);
            res.raw_length_nt = total_raw_nt;
            if grp.parts.len() > 1 {
                res.alignment_mode = format!(
                    "Coordinate-stitched ({}/{} RNA segments to {})",
                    grp.parts.len(),
                    n_lm,
                    self.spec.reference_name
                );
            }

            let mut slots = Vec::with_capacity(n_lm);
            let mut present_cnt = 0usize;
            for (si, lm) in self.spec.landmarks.iter().enumerate() {
                if let Some((src_h, r_len, m_len)) = &slot_meta[si] {
                    present_cnt += 1;
                    slots.push(SegmentSlotInfo {
                        seg_name: lm.short_name.clone(),
                        seg_index: si + 1,
                        ref_start: lm.start,
                        ref_end: lm.end,
                        present: true,
                        source_header: Some(src_h.clone()),
                        raw_len_nt: *r_len,
                        mapped_nt: *m_len,
                    });
                } else {
                    slots.push(SegmentSlotInfo {
                        seg_name: lm.short_name.clone(),
                        seg_index: si + 1,
                        ref_start: lm.start,
                        ref_end: lm.end,
                        present: false,
                        source_header: None,
                        raw_len_nt: 0,
                        mapped_nt: 0,
                    });
                }
            }

            isolates.push(SegmentedIsolateSummary {
                isolate_id: grp.isolate_id.clone(),
                n_segments_present: present_cnt,
                total_segments: n_lm,
                is_assembled_from_parts: grp.parts.len() > 1,
                total_raw_nt,
                total_mapped_nt: res.mapped_coding_nt,
                segments: slots,
            });
            results.push(res);
            out_fasta.push('>');
            out_fasta.push_str(&grp.isolate_id);
            out_fasta.push('\n');
            out_fasta.push_str(&stitched_seq);
            out_fasta.push('\n');
        }

        // 3. Any unmerged/single passthrough short records
        for (hdr, seq) in passthrough_short {
            let res = self.classify_sequence_with_temporal(&hdr, &seq, cfg, temporal);
            let cov_vec = self.compute_landmark_coverage_vec(&seq);
            let mut slots = Vec::with_capacity(n_lm);
            let mut present_cnt = 0usize;
            for (si, lm) in self.spec.landmarks.iter().enumerate() {
                let mapped_nt = cov_vec.get(si).copied().unwrap_or(0);
                let present = mapped_nt >= 50;
                if present {
                    present_cnt += 1;
                }
                slots.push(SegmentSlotInfo {
                    seg_name: lm.short_name.clone(),
                    seg_index: si + 1,
                    ref_start: lm.start,
                    ref_end: lm.end,
                    present,
                    source_header: if present { Some(hdr.clone()) } else { None },
                    raw_len_nt: if present { mapped_nt } else { 0 },
                    mapped_nt,
                });
            }
            isolates.push(SegmentedIsolateSummary {
                isolate_id: hdr.clone(),
                n_segments_present: present_cnt,
                total_segments: n_lm,
                is_assembled_from_parts: false,
                total_raw_nt: res.raw_length_nt,
                total_mapped_nt: res.mapped_coding_nt,
                segments: slots,
            });
            results.push(res);
            out_fasta.push('>');
            out_fasta.push_str(&hdr);
            out_fasta.push('\n');
            out_fasta.push_str(&seq);
            out_fasta.push('\n');
        }

        let assembly = SegmentedAssemblyResult {
            mode: effective_mode.to_string(),
            raw_record_count,
            isolate_count: isolates.len(),
            was_multi_record_assembled,
            isolates,
        };
        (assembly, results, out_fasta)
    }

    /// Compute valid mapped nucleotides within each landmark span for a sequence.
    fn compute_landmark_coverage_vec(&self, raw_seq: &str) -> Vec<usize> {
        let (cod, _, _) = self.project_to_hxb2_coding(raw_seq);
        self.spec
            .landmarks
            .iter()
            .map(|lm| {
                if lm.end < self.spec.coding_start || lm.start > self.spec.coding_end {
                    return 0;
                }
                let lo = lm.start.saturating_sub(self.spec.coding_start);
                let hi = lm
                    .end
                    .min(self.spec.coding_end)
                    .saturating_sub(self.spec.coding_start);
                (lo..=hi)
                    .filter(|&ci| cod.get(ci).copied().unwrap_or(-1) >= 0)
                    .count()
            })
            .collect()
    }

    /// Project any aligned or unaligned sequence onto the reference coding coordinate grid (`coding_start..=coding_end`).
    /// Returns `(coding_array, alignment_mode_string, raw_ungapped_len)`.
    pub fn project_to_hxb2_coding(&self, raw_seq: &str) -> (Vec<i8>, String, usize) {
        let seq_bytes = raw_seq.as_bytes();
        let total_len = seq_bytes.len();
        let ungapped: Vec<i8> = seq_bytes
            .iter()
            .filter(|&&b| b != b'-' && b != b'.')
            .map(|&b| base_to_int(b))
            .collect();
        let raw_ungapped_len = ungapped.len();
        let mut notes: Vec<String> = Vec::new();
        let ref_len = self.spec.ref_full_len;
        let cod_start = self.spec.coding_start;
        let cod_end = self.spec.coding_end;
        let cod_len = self.spec.coding_len();

        let layout_case = if total_len == ref_len {
            Some(1)
        } else if total_len == cod_len {
            Some(2)
        } else if self.lanl_msa_len > 0 && total_len == self.lanl_msa_len {
            Some(3)
        } else {
            None
        };

        if let Some(case) = layout_case {
            if self.prealigned_plausible(seq_bytes, case) {
                match case {
                    // Case 1: Exact full reference-anchored alignment
                    1 => {
                        let full_arr = seq_to_int_vec(raw_seq);
                        let cod = full_arr[(cod_start - 1)..cod_end].to_vec();
                        return (
                            cod,
                            format!(
                                "Pre-aligned ({} {} nt)",
                                self.spec.reference_name,
                                fmt_nt(ref_len)
                            ),
                            raw_ungapped_len,
                        );
                    }
                    // Case 2: Exact coding-span alignment
                    2 => {
                        let cod = seq_to_int_vec(raw_seq);
                        return (
                            cod,
                            format!(
                                "Pre-aligned ({} coding {} nt)",
                                self.spec.reference_name,
                                fmt_nt(cod_len)
                            ),
                            raw_ungapped_len,
                        );
                    }
                    // Case 3: Exact MSA coordinate table alignment (e.g. LANL 2023 11,531-nt MSA)
                    _ => {
                        let full_arr = seq_to_int_vec(raw_seq);
                        let mut cod = vec![-1i8; cod_len];
                        for (msa_col, &h_pos) in self.lanl_msa_to_hxb2.iter().enumerate() {
                            if h_pos >= cod_start && h_pos <= cod_end {
                                let c_idx = h_pos - cod_start;
                                if cod[c_idx] < 0 && full_arr[msa_col] >= 0 {
                                    cod[c_idx] = full_arr[msa_col];
                                }
                            }
                        }
                        return (
                            cod,
                            format!(
                                "Pre-aligned (LANL 2023 MSA {} nt)",
                                fmt_nt(self.lanl_msa_len)
                            ),
                            raw_ungapped_len,
                        );
                    }
                }
            }
            notes.push("length matched pre-aligned layout but anchors disagreed".to_string());
        }

        // Case 4: Unaligned or differently-aligned sequence -> Seed-and-extend coordinate projection onto reference
        let (full_ref, mut flags) = self.align_unaligned_to_ref(&ungapped);
        notes.append(&mut flags);
        let cod = full_ref[(cod_start - 1)..cod_end].to_vec();
        let mut mode = format!(
            "Auto-aligned (15-mer anchor + banded NW to {}",
            self.spec.reference_name
        );
        for n in &notes {
            mode.push_str("; ");
            mode.push_str(n);
        }
        mode.push(')');
        (cod, mode, raw_ungapped_len)
    }

    /// Sanity check for records whose length equals a pre-aligned layout: gapped records are trusted;
    /// gap-free ones are probed with 15-mers against the reference/trunk index at the positions the layout predicts.
    fn prealigned_plausible(&self, seq: &[u8], case: u8) -> bool {
        if seq.iter().any(|&b| b == b'-' || b == b'.') {
            return true;
        }
        let cod_start = self.spec.coding_start;
        let cod_end = self.spec.coding_end;
        let cod_len = self.spec.coding_len();
        // (first column, last column exclusive) to probe, and the expected 0-based reference position of a column
        let (lo, hi) = match case {
            1 => (cod_start - 1, cod_end),
            2 => (0, cod_len),
            _ => (0, seq.len()),
        };
        let expected = |c: usize| -> Option<usize> {
            match case {
                1 => Some(c),
                2 => Some(c + cod_start - 1),
                _ => {
                    let h = *self.lanl_msa_to_hxb2.get(c)?;
                    if h >= cod_start && h <= cod_end {
                        Some(h - 1)
                    } else {
                        None
                    }
                }
            }
        };
        const N_PROBES: usize = 16;
        let span = hi.saturating_sub(lo).saturating_sub(40);
        let (mut agree, mut disagree) = (0usize, 0usize);
        for i in 0..N_PROBES {
            let start = lo + i * span / (N_PROBES - 1);
            // the index holds even reference positions only
            let mut c = start;
            let target = loop {
                if c >= start + 8 || c + 15 > seq.len() {
                    break None;
                }
                match expected(c) {
                    Some(e) if e % 2 == 0 => break Some(e),
                    _ => c += 1,
                }
            };
            let Some(target) = target else { continue };
            let bases: Vec<i8> = seq[c..c + 15].iter().map(|&b| base_to_int(b)).collect();
            let Some(code) = encode_kmer_15(&bases) else { continue };
            if let Some(locs) = self.kmer_index.get(&code) {
                if locs.iter().any(|&(_, r)| (r as usize).abs_diff(target) <= 3) {
                    agree += 1;
                } else if !locs.is_empty() {
                    disagree += 1;
                }
            }
        }
        disagree <= agree
    }

    /// Collect seed hits for one strand and chain them into non-overlapping anchors `(r_pos, q_pos, k)`.
    /// Returns the anchors and the chain score (`-1` when no chain exists).
    fn seed_and_chain(&self, q: &[i8]) -> (Vec<(usize, usize, usize)>, i32) {
        let qlen = q.len();
        // 1. Collect seed hits: (r_pos, q_pos, prof_idx, k)
        let mut hits: Vec<(usize, usize, usize, usize)> = Vec::new();
        for q_pos in 0..=(qlen.saturating_sub(15)) {
            if let Some(code) = encode_kmer_15(&q[q_pos..(q_pos + 15)]) {
                if let Some(locs) = self.kmer_index.get(&code) {
                    for &(p_idx, r_pos) in locs {
                        hits.push((r_pos as usize, q_pos, p_idx as usize, 15));
                    }
                }
            }
        }

        // If very few exact 15-mer hits, fall back to 11-mer scanning against indexed profiles
        if hits.len() < 4 && qlen >= 80 {
            for (p_idx, prof) in self.full_profiles_ref.iter().enumerate().take(16) {
                let mut k11_map: HashMap<u32, usize> = HashMap::new();
                for r_pos in (0..(prof.len().saturating_sub(11))).step_by(3) {
                    let slice = &prof[r_pos..(r_pos + 11)];
                    if slice.iter().all(|&b| b >= 0 && b <= 3) {
                        let mut v = 0u32;
                        for &b in slice {
                            v = (v << 2) | (b as u32);
                        }
                        k11_map.entry(v).or_insert(r_pos);
                    }
                }
                for q_pos in (0..=(qlen.saturating_sub(11))).step_by(2) {
                    let slice = &q[q_pos..(q_pos + 11)];
                    if slice.iter().all(|&b| b >= 0 && b <= 3) {
                        let mut v = 0u32;
                        for &b in slice {
                            v = (v << 2) | (b as u32);
                        }
                        if let Some(&r_pos) = k11_map.get(&v) {
                            hits.push((r_pos, q_pos, p_idx, 11));
                        }
                    }
                }
            }
        }

        if hits.is_empty() {
            return (Vec::new(), -1);
        }

        // Deduplicate (r_pos, q_pos) keeping the lowest profile index (stable, so 15-mers win over 11-mers)
        hits.sort_by_key(|&(r, q_p, p_idx, _)| (q_p, r, p_idx));
        hits.dedup_by_key(|h| (h.0, h.1));

        // 2. Find dominant diagonal bin (r_pos - q_pos) to filter spurious off-diagonal repeats
        let mut diag_counts: BTreeMap<i32, usize> = BTreeMap::new();
        for &(r, qp, _, _) in &hits {
            let bin = ((r as i32) - (qp as i32)).div_euclid(250);
            *diag_counts.entry(bin).or_insert(0) += 1;
            *diag_counts.entry(bin - 1).or_insert(0) += 1;
            *diag_counts.entry(bin + 1).or_insert(0) += 1;
        }
        let best_bin = diag_counts
            .into_iter()
            .max_by_key(|&(bin, count)| (count, Reverse(bin.abs()), Reverse(bin)))
            .map(|(b, _)| b)
            .unwrap_or(0);
        let center_diag = best_bin * 250;

        let filtered: Vec<(usize, usize, usize)> = hits
            .into_iter()
            .filter(|&(r, qp, _, _)| ((r as i32 - qp as i32) - center_diag).abs() <= 600)
            .map(|(r, qp, _, k)| (r, qp, k))
            .collect();

        if filtered.is_empty() {
            return (Vec::new(), -1);
        }

        // 3. Collinear sparse DP chaining over anchors
        let n_hits = filtered.len();
        let mut dp: Vec<i32> = filtered.iter().map(|&(_, _, k)| k as i32).collect();
        let mut prev = vec![usize::MAX; n_hits];
        let max_jump = if self.spec.ref_full_len <= 10_000 {
            2500
        } else {
            (self.spec.ref_full_len as i32) / 3
        };

        for i in 0..n_hits {
            let (ri, qi, ki) = filtered[i];
            let ki = ki as i32;
            let start_j = i.saturating_sub(120);
            for j in start_j..i {
                let (rj, qj, _) = filtered[j];
                if rj < ri && qj < qi {
                    let dr = (ri - rj) as i32;
                    let dq = (qi - qj) as i32;
                    let indel = (dr - dq).abs();
                    if indel <= 180 && dr <= max_jump && dq <= max_jump {
                        let overlap = (ki - dr.min(dq)).max(0);
                        let gain = ki - overlap - (indel / 3);
                        let cand = dp[j] + gain;
                        if cand > dp[i] {
                            dp[i] = cand;
                            prev[i] = j;
                        }
                    }
                }
            }
        }

        let mut best_end = 0;
        for i in 1..n_hits {
            if dp[i] > dp[best_end] {
                best_end = i;
            }
        }
        let score = dp[best_end];

        let mut chain: Vec<(usize, usize, usize)> = Vec::new();
        let mut cur = best_end;
        while cur != usize::MAX {
            chain.push(filtered[cur]);
            cur = prev[cur];
        }
        chain.reverse();

        // Merge overlapping anchors into non-overlapping anchor points
        let mut clean_anchors: Vec<(usize, usize, usize)> = Vec::new();
        for (r, qp, k) in chain {
            if let Some(&(lr, lq, lk)) = clean_anchors.last() {
                if r >= lr + lk && qp >= lq + lk {
                    clean_anchors.push((r, qp, k));
                }
            } else {
                clean_anchors.push((r, qp, k));
            }
        }
        (clean_anchors, score)
    }

    /// Aligns an unaligned sequence `q` onto the pathogen's reference coordinate grid using 15-mer seed chaining
    /// and banded dynamic programming between collinear anchors. Both strands are chained and the
    /// stronger one is used. The second return value lists diagnostic notes for the mode string.
    fn align_unaligned_to_ref(&self, q: &[i8]) -> (Vec<i8>, Vec<String>) {
        let ref_len = self.spec.ref_full_len;
        let max_q = self.spec.effective_max_query_len();
        let qlen = q.len();
        let mut notes: Vec<String> = Vec::new();
        if qlen < 40 {
            return (vec![-1i8; ref_len], notes);
        }
        if qlen > max_q {
            notes.push(format!("rejected: >{} nt", max_q));
            return (vec![-1i8; ref_len], notes);
        }

        let (fwd_anchors, fwd_score) = self.seed_and_chain(q);
        let rc: Vec<i8> = q
            .iter()
            .rev()
            .map(|&b| if (0..=3).contains(&b) { 3 - b } else { -1 })
            .collect();
        let (rc_anchors, rc_score) = self.seed_and_chain(&rc);

        let (clean_anchors, q, reverse) = if rc_score > fwd_score {
            (rc_anchors, &rc[..], true)
        } else {
            (fwd_anchors, q, false)
        };
        if clean_anchors.len() < 3 {
            return (vec![-1i8; ref_len], notes);
        }
        if reverse {
            notes.push("reverse-complemented".to_string());
        }

        let mut aligned = vec![-1i8; ref_len];
        let mut proportional_len = 0usize;

        // Fill left unanchored tail if short/moderate
        let (first_r, first_q, _) = clean_anchors[0];
        let mut left_untouched = 0usize;
        if first_q > 0 && first_r > 0 {
            let take_q = first_q.min(400);
            let take_r = first_r.min(take_q + 45);
            left_untouched = first_q - take_q;
            proportional_len = proportional_len.max(self.banded_align_segment(
                &q[(first_q - take_q)..first_q],
                first_r - take_r,
                first_r,
                true,
                false,
                &mut aligned,
            ));
        }

        // Fill exact anchors and banded gaps between consecutive anchors
        for idx in 0..clean_anchors.len() {
            let (r, qp, k) = clean_anchors[idx];
            for off in 0..k {
                if r + off < ref_len && qp + off < qlen {
                    aligned[r + off] = q[qp + off];
                }
            }
            if idx + 1 < clean_anchors.len() {
                let (next_r, next_q, _) = clean_anchors[idx + 1];
                let gap_r_start = r + k;
                let gap_q_start = qp + k;
                if next_r > gap_r_start && next_q > gap_q_start {
                    proportional_len = proportional_len.max(self.banded_align_segment(
                        &q[gap_q_start..next_q],
                        gap_r_start,
                        next_r,
                        false,
                        false,
                        &mut aligned,
                    ));
                }
            }
        }

        // Fill right unanchored tail
        let &(last_r, last_q, last_k) = clean_anchors.last().unwrap();
        let rem_q_start = last_q + last_k;
        let rem_r_start = last_r + last_k;
        let mut right_untouched = 0usize;
        if rem_q_start < qlen && rem_r_start < ref_len {
            let take_q = (qlen - rem_q_start).min(400);
            let take_r = (ref_len - rem_r_start).min(take_q + 45);
            right_untouched = qlen - rem_q_start - take_q;
            proportional_len = proportional_len.max(self.banded_align_segment(
                &q[rem_q_start..(rem_q_start + take_q)],
                rem_r_start,
                rem_r_start + take_r,
                false,
                true,
                &mut aligned,
            ));
        }

        if proportional_len > 0 {
            notes.push(format!("proportional projection over {} nt", proportional_len));
        }
        if left_untouched > 0 {
            notes.push(format!(
                "5' tail truncated after 400 nt ({} nt not aligned)",
                left_untouched
            ));
        }
        if right_untouched > 0 {
            notes.push(format!(
                "3' tail truncated after 400 nt ({} nt not aligned)",
                right_untouched
            ));
        }

        (aligned, notes)
    }

    /// Banded Needleman-Wunsch projection of `q_sub` onto reference interval `[r_start..r_end]`,
    /// scoring matches against the consensus profile ensemble so any basal trunk aligns cleanly.
    /// `free_ref_start` / `free_ref_end` make the reference slack at that end free (semi-global), leaving
    /// unused reference columns at -1. Returns the interval length when the interval was too large for
    /// DP and was projected proportionally, else 0.
    fn banded_align_segment(
        &self,
        q_sub: &[i8],
        r_start: usize,
        r_end: usize,
        free_ref_start: bool,
        free_ref_end: bool,
        aligned_out: &mut [i8],
    ) -> usize {
        let m = q_sub.len();
        let n = r_end.saturating_sub(r_start);
        if m == 0 || n == 0 {
            return 0;
        }
        // If interval is huge, project proportionally to avoid quadratic memory
        if m > 1200 || n > 1200 {
            for j in 0..n {
                let qi = (j * m) / n;
                aligned_out[r_start + j] = q_sub[qi];
            }
            return n.max(m);
        }

        let diff = (m as isize - n as isize).unsigned_abs();
        let band = diff + 35;
        let neg_inf = -1_000_000i32;

        let cols = n + 1;
        let mut score = vec![neg_inf; (m + 1) * cols];
        let mut trace = vec![0u8; (m + 1) * cols]; // 1=diag, 2=up (insertion in q), 3=left (deletion in q)

        score[0] = 0;
        for i in 1..=m.min(band) {
            score[i * cols] = -(i as i32) * 3;
            trace[i * cols] = 2;
        }
        for j in 1..=n.min(band) {
            score[j] = if free_ref_start { 0 } else { -(j as i32) * 3 };
            trace[j] = 3;
        }
        if free_ref_start {
            for j in 1..=n {
                score[j] = 0;
                trace[j] = 3;
            }
        }

        for i in 1..=m {
            let center_j = (i * n) / m;
            let j_lo = center_j.saturating_sub(band).max(1);
            let j_hi = (center_j + band).min(n);
            let qb = q_sub[i - 1];

            for j in j_lo..=j_hi {
                let r_pos = r_start + j - 1;
                // Match if query base appears in reference, any basal trunk, or indexed founders at r_pos
                let matched = (0..=3).contains(&qb)
                    && (self.base_mask_ref[r_pos] & (1u8 << (qb as u8))) != 0;
                let sub_score = if matched { 2 } else { -1 };

                let diag = score[(i - 1) * cols + (j - 1)] + sub_score;
                let up = score[(i - 1) * cols + j] - 3;
                let left = score[i * cols + (j - 1)] - 3;

                let idx = i * cols + j;
                if diag >= up && diag >= left {
                    score[idx] = diag;
                    trace[idx] = 1;
                } else if up >= left {
                    score[idx] = up;
                    trace[idx] = 2;
                } else {
                    score[idx] = left;
                    trace[idx] = 3;
                }
            }
        }

        let mut i = m;
        let mut j = n;
        if free_ref_end {
            // Traceback from the best-scoring cell of the last row; later reference columns stay unaligned
            let j_lo = n.saturating_sub(band).max(1);
            let mut best = neg_inf;
            for jj in j_lo..=n {
                let s = score[m * cols + jj];
                if s >= best {
                    best = s;
                    j = jj;
                }
            }
        }
        while i > 0 && j > 0 {
            let t = trace[i * cols + j];
            if t == 1 {
                aligned_out[r_start + j - 1] = q_sub[i - 1];
                i -= 1;
                j -= 1;
            } else if t == 2 {
                i -= 1;
            } else if t == 3 {
                aligned_out[r_start + j - 1] = -1;
                j -= 1;
            } else {
                // Fallback if outside band
                aligned_out[r_start + j - 1] = q_sub[i - 1];
                i -= 1;
                j -= 1;
            }
        }
        0
    }

    /// Resolve admissible Tier-1 trunks for a query under `TemporalGuardConfig`.
    /// Returns `(query_year, temporal_limit, allowed_trunks, masked_future_trunks)`.
    fn resolve_temporal_mask(
        &self,
        header: &str,
        temporal: &TemporalGuardConfig,
    ) -> (Option<f64>, Option<f64>, Vec<bool>, Vec<String>) {
        let k = self.trunk_names.len();
        if !temporal.enabled {
            return (None, None, vec![true; k], Vec::new());
        }
        let query_year = temporal
            .query_year
            .or_else(|| parse_header_decimal_year(header));
        let Some(t_q) = query_year else {
            return (None, None, vec![true; k], Vec::new());
        };
        let slack = temporal
            .slack_yr
            .unwrap_or_else(|| self.spec.effective_temporal_slack_yr());
        let limit = t_q + slack;

        let mut allowed: Vec<bool> = self
            .trunk_floors
            .iter()
            .map(|fl| fl.is_none_or(|v| v <= limit))
            .collect();

        // Safety fallback: ensure at least min(2, k) trunks and >= 1 pure basal trunk remain active.
        let min_keep = k.min(2);
        let need_basal = self.pure_basal_count > 0
            && !(0..self.pure_basal_count).any(|b| allowed[b]);

        if need_basal {
            if let Some((best_b, _)) = (0..self.pure_basal_count)
                .map(|b| (b, self.trunk_floors[b].unwrap_or(f64::NEG_INFINITY)))
                .min_by(|a, b| a.1.total_cmp(&b.1).then_with(|| a.0.cmp(&b.0)))
            {
                allowed[best_b] = true;
            }
        }

        if allowed.iter().filter(|&&a| a).count() < min_keep {
            let mut masked: Vec<usize> = (0..k).filter(|&i| !allowed[i]).collect();
            masked.sort_by(|&a, &b| {
                let pa = a < self.pure_basal_count;
                let pb = b < self.pure_basal_count;
                pb.cmp(&pa)
                    .then_with(|| {
                        self.trunk_floors[a]
                            .unwrap_or(f64::INFINITY)
                            .total_cmp(&self.trunk_floors[b].unwrap_or(f64::INFINITY))
                    })
                    .then_with(|| a.cmp(&b))
            });
            for idx in masked {
                if allowed.iter().filter(|&&a| a).count() >= min_keep {
                    break;
                }
                allowed[idx] = true;
            }
        }

        let masked_future_trunks: Vec<String> = (0..k)
            .filter(|&i| !allowed[i])
            .map(|i| self.trunk_names[i].clone())
            .collect();

        (Some(t_q), Some(limit), allowed, masked_future_trunks)
    }

    fn is_compendium_entry_admissible(
        &self,
        entry: &CompendiumEntry,
        allowed_trunks: &[bool],
        temporal_limit: f64,
    ) -> bool {
        let floor = min_temporal_floor(
            entry.t_mrca_ci_low,
            entry.t_mrca,
            entry.first_sampled_year,
        );
        if floor.is_some_and(|fl| fl > temporal_limit) {
            return false;
        }
        if let Some(&bi) = self.trunk_alias_map.get(&entry.backbone) {
            if !allowed_trunks.get(bi).copied().unwrap_or(true) {
                return false;
            }
        }
        for (seg_trk, _, _) in parse_barcode_segments(&entry.hxb2_barcode) {
            if let Some(&ti) = self.trunk_alias_map.get(seg_trk.as_str()) {
                if !allowed_trunks.get(ti).copied().unwrap_or(true) {
                    return false;
                }
            }
        }
        true
    }

    /// Screen a single sequence `(header, raw_seq)` and return its full `SequenceResult`.
    pub fn classify_sequence(
        &self,
        header: &str,
        raw_seq: &str,
        cfg: &ScreeningConfig,
    ) -> SequenceResult {
        self.classify_sequence_with_temporal(
            header,
            raw_seq,
            cfg,
            &TemporalGuardConfig::default(),
        )
    }

    /// Screen a single sequence `(header, raw_seq)` with optional non-lookahead temporal guarding.
    pub fn classify_sequence_with_temporal(
        &self,
        header: &str,
        raw_seq: &str,
        cfg: &ScreeningConfig,
        temporal: &TemporalGuardConfig,
    ) -> SequenceResult {
        let (query_year, temporal_limit, allowed_trunks, masked_future_trunks) =
            self.resolve_temporal_mask(header, temporal);
        let allowed_slice = if temporal.enabled && temporal_limit.is_some() {
            Some(allowed_trunks.as_slice())
        } else {
            None
        };

        let (mut q_coding, align_mode, raw_len) = self.project_to_hxb2_coding(raw_seq);
        let coding_len = self.spec.coding_len();
        let coding_start = self.spec.coding_start;

        let mapped_coding_nt = q_coding.iter().filter(|&&b| b >= 0).count();

        // Apply artifact shield for hypervariable loops if enabled (before the span is determined)
        if cfg.mask_hypervariable_loops {
            for i in 0..coding_len {
                if self.loop_mask_coding[i] {
                    q_coding[i] = -1;
                }
            }
        }

        // Determine active coding span [active_start..active_end] where usable nucleotides exist
        let valid_positions: Vec<usize> = q_coding
            .iter()
            .enumerate()
            .filter_map(|(i, &b)| if b >= 0 { Some(i) } else { None })
            .collect();

        if valid_positions.len() < self.spec.min_mapped_coding_nt {
            let unclass_formula = if cfg.mask_hypervariable_loops {
                format!(
                    "Unclassifiable (<{} usable coding nt after loop mask)",
                    self.spec.min_mapped_coding_nt
                )
            } else {
                format!(
                    "Unclassifiable (<{} mapped coding nt)",
                    self.spec.min_mapped_coding_nt
                )
            };
            return SequenceResult {
                sequence_id: header.to_string(),
                raw_length_nt: raw_len,
                mapped_coding_nt,
                hxb2_span_start: 0,
                hxb2_span_end: 0,
                locus_scope: "insufficient_coverage".to_string(),
                alignment_mode: align_mode,
                called_lineage: "Unclassifiable".to_string(),
                called_category: "Unclassifiable".to_string(),
                call_support: None,
                epidemiological_formula: unclass_formula.clone(),
                epidemiological_barcode: String::new(),
                epidemiological_segments: vec![],
                basal_formula: unclass_formula.clone(),
                basal_barcode: String::new(),
                structural_formula: unclass_formula,
                descriptive_name: "Unclassifiable".to_string(),
                backbone_trunk: "None".to_string(),
                constituent_trunks: vec![],
                is_clonal: false,
                n_segments: 0,
                n_breakpoints: 0,
                hxb2_barcode: String::new(),
                ci_bracketed_barcode: String::new(),
                segments: vec![],
                breakpoints: vec![],
                matches: vec![],
                query_year,
                masked_future_trunks,
            };
        }

        let first_idx = *valid_positions.first().unwrap();
        let last_idx = *valid_positions.last().unwrap();
        let h_start = coding_start + first_idx;
        let h_end = coding_start + last_idx;

        // Determine locus scope via PathogenSpec rules
        let span_width = h_end - h_start + 1;
        let locus_scope = self
            .spec
            .classify_locus_scope(h_start, h_end, span_width, mapped_coding_nt);

        let gamma = cfg.gamma.unwrap_or(7.0);
        let min_span_nt = cfg.min_span_nt.unwrap_or(if locus_scope == "full" {
            120
        } else {
            150
        });

        // Slice to [first_idx..=last_idx] so leading/trailing unsequenced regions in partial genes don't dilute edges
        let sub_len = last_idx - first_idx + 1;
        let q_sub = &q_coding[first_idx..=last_idx];
        let col_to_hxb2: Vec<usize> = (h_start..=h_end).collect();

        let (segs, bps) = self.run_paff_viterbi(
            allowed_slice,
            q_sub,
            first_idx,
            &col_to_hxb2,
            sub_len,
            gamma,
            min_span_nt,
            cfg.min_informative_snps,
            cfg.p_fisher_thresh,
        );

        let (desc_name, backbone, struct_formula) =
            format_descriptive_name(&segs, &locus_scope);
        let hxb2_barcode = segs
            .iter()
            .map(|s| format!("{}:{}-{}", s.trunk, s.hxb2_start, s.hxb2_end))
            .collect::<Vec<_>>()
            .join("|");
        let ci_bracketed_barcode = segs
            .iter()
            .map(|s| {
                format!(
                    "{}:{}-{}[CI:{}-{}]",
                    s.trunk, s.hxb2_start, s.hxb2_end, s.ci_hxb2_start, s.ci_hxb2_end
                )
            })
            .collect::<Vec<_>>()
            .join("|");

        let mut constituent_set: HashSet<String> = HashSet::new();
        for s in &segs {
            constituent_set.insert(s.trunk.clone());
        }
        let mut constituent_trunks: Vec<String> = constituent_set.into_iter().collect();
        constituent_trunks.sort();

        let mut matches = self.match_compendium(
            &struct_formula,
            &segs,
            &bps,
            &locus_scope,
            allowed_slice,
            temporal_limit,
        );
        let decision = self.call_certified_lineage(
            &q_coding,
            h_start,
            h_end,
            &locus_scope,
            &segs,
            &bps,
            cfg,
            allowed_slice,
            temporal_limit,
        );
        if let Some(comp_label) = &decision.certified_compendium_label {
            self.prepend_certified_match(&mut matches, comp_label, &decision.call_support);
        }

        let (epidemiological_formula, epidemiological_barcode, epidemiological_segments) = self
            .build_epidemiological_map(
                &decision,
                &locus_scope,
                h_start,
                h_end,
                q_sub,
                first_idx,
                &col_to_hxb2,
                &segs,
                cfg,
                allowed_slice,
                temporal_limit,
            );

        let basal_formula = desc_name.clone();
        let basal_barcode = hxb2_barcode.clone();

        SequenceResult {
            sequence_id: header.to_string(),
            raw_length_nt: raw_len,
            mapped_coding_nt,
            hxb2_span_start: h_start,
            hxb2_span_end: h_end,
            locus_scope,
            alignment_mode: align_mode,
            called_lineage: decision.called_lineage,
            called_category: decision.called_category,
            call_support: Some(decision.call_support),
            epidemiological_formula,
            epidemiological_barcode,
            epidemiological_segments,
            basal_formula,
            basal_barcode,
            structural_formula: struct_formula,
            descriptive_name: desc_name,
            backbone_trunk: backbone,
            constituent_trunks,
            is_clonal: segs.len() == 1,
            n_segments: segs.len(),
            n_breakpoints: bps.len(),
            hxb2_barcode,
            ci_bracketed_barcode,
            segments: segs,
            breakpoints: bps,
            matches,
            query_year,
            masked_future_trunks,
        }
    }

    #[allow(clippy::too_many_arguments)]
    fn run_paff_viterbi(
        &self,
        allowed_trunks: Option<&[bool]>,
        q_sub: &[i8],
        offset_in_coding: usize,
        col_to_hxb2: &[usize],
        l_sub: usize,
        gamma: f64,
        min_hxb2_span_nt: usize,
        min_informative_snps: usize,
        p_fisher_thresh: f64,
    ) -> (Vec<SegmentCall>, Vec<BreakpointCall>) {
        let mut basal_penalties: Vec<f64> = self
            .spec
            .basal_trunks
            .iter()
            .filter(|t| t.is_pure_basal)
            .chain(self.spec.basal_trunks.iter().filter(|t| !t.is_pure_basal))
            .map(|t| t.penalty)
            .collect();
        basal_penalties.resize(self.trunk_names.len(), 0.0);
        let allow_ghost = self.spec.ghost_relabel_rule.as_ref().is_none_or(|rule| {
            self.trunk_alias_map.get(&rule.ghost_trunk).is_none_or(|&idx| {
                allowed_trunks.is_none_or(|m| m.get(idx).copied().unwrap_or(true))
            })
        });
        self.run_paff_viterbi_panel(
            &self.trunk_names,
            &self.trunk_mat,
            &basal_penalties,
            allowed_trunks,
            allow_ghost,
            q_sub,
            offset_in_coding,
            col_to_hxb2,
            l_sub,
            gamma,
            min_hxb2_span_nt,
            min_informative_snps,
            p_fisher_thresh,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn run_paff_viterbi_panel(
        &self,
        raw_names: &[String],
        raw_mat: &[Vec<i8>],
        raw_penalties: &[f64],
        allowed: Option<&[bool]>,
        allow_ghost_e_relabel: bool,
        q_sub: &[i8],
        offset_in_coding: usize,
        col_to_hxb2: &[usize],
        l_sub: usize,
        gamma: f64,
        min_hxb2_span_nt: usize,
        min_informative_snps: usize,
        p_fisher_thresh: f64,
    ) -> (Vec<SegmentCall>, Vec<BreakpointCall>) {
        let (names, mat, penalties): (Vec<&str>, Vec<&[i8]>, Vec<f64>) = match allowed {
            Some(mask) if mask.iter().any(|&a| !a) => {
                let mut ns = Vec::with_capacity(raw_names.len());
                let mut ms = Vec::with_capacity(raw_names.len());
                let mut ps = Vec::with_capacity(raw_names.len());
                for i in 0..raw_names.len() {
                    if mask.get(i).copied().unwrap_or(true) {
                        ns.push(raw_names[i].as_str());
                        ms.push(raw_mat[i].as_slice());
                        ps.push(raw_penalties[i]);
                    }
                }
                if ms.is_empty() {
                    (
                        raw_names.iter().map(|s| s.as_str()).collect(),
                        raw_mat.iter().map(|v| v.as_slice()).collect(),
                        raw_penalties.to_vec(),
                    )
                } else {
                    (ns, ms, ps)
                }
            }
            _ => (
                raw_names.iter().map(|s| s.as_str()).collect(),
                raw_mat.iter().map(|v| v.as_slice()).collect(),
                raw_penalties.to_vec(),
            ),
        };

        let k_trunks = names.len();
        let k_total = k_trunks + 1;
        let ghost_idx = k_trunks;
        let mismatch_cost = 1.0f64;
        let ghost_cost = 0.58f64;

        let mut dp = vec![0.0f64; k_total * l_sub];
        let mut back = vec![0u16; k_total * l_sub];

        // Site 0 initialization
        for k in 0..k_trunks {
            let base_c =
                self.site_emission_cost(q_sub[0], mat[k][offset_in_coding], mismatch_cost);
            dp[k * l_sub] = if q_sub[0] >= 0 {
                base_c + penalties[k]
            } else {
                base_c
            };
        }
        dp[ghost_idx * l_sub] = if q_sub[0] >= 0 { ghost_cost } else { 0.0 };

        for u in 1..l_sub {
            let mut best_prev_idx = 0usize;
            let mut best_prev_val = dp[u - 1];
            for k in 1..k_total {
                let v = dp[k * l_sub + (u - 1)];
                if v < best_prev_val {
                    best_prev_val = v;
                    best_prev_idx = k;
                }
            }
            let step_gamma = if self.spec.is_segmented {
                let prev_seg = self.segment_id_at_ref(col_to_hxb2[u - 1]);
                let curr_seg = self.segment_id_at_ref(col_to_hxb2[u]);
                if prev_seg != 0 && curr_seg != 0 && prev_seg != curr_seg {
                    0.0
                } else if self.spec.forbid_intra_segment_recomb {
                    1.0e6
                } else {
                    gamma
                }
            } else {
                gamma
            };
            let switch_cost = best_prev_val + step_gamma;
            let cod_pos = offset_in_coding + u;

            for k in 0..k_total {
                let stay_cost = dp[k * l_sub + (u - 1)];
                let emit = if k < k_trunks {
                    let base_c =
                        self.site_emission_cost(q_sub[u], mat[k][cod_pos], mismatch_cost);
                    if q_sub[u] >= 0 {
                        base_c + penalties[k]
                    } else {
                        base_c
                    }
                } else if q_sub[u] >= 0 {
                    ghost_cost
                } else {
                    0.0
                };
                if switch_cost < stay_cost {
                    dp[k * l_sub + u] = switch_cost + emit;
                    back[k * l_sub + u] = best_prev_idx as u16;
                } else {
                    dp[k * l_sub + u] = stay_cost + emit;
                    back[k * l_sub + u] = k as u16;
                }
            }
        }

        let mut path = vec![0usize; l_sub];
        let mut best_last = 0usize;
        let mut best_last_val = dp[l_sub - 1];
        for k in 1..k_total {
            let v = dp[k * l_sub + (l_sub - 1)];
            if v < best_last_val {
                best_last_val = v;
                best_last = k;
            }
        }
        path[l_sub - 1] = best_last;
        for u in (1..l_sub).rev() {
            path[u - 1] = back[path[u] * l_sub + u] as usize;
        }

        let mut segs: Vec<(usize, usize, usize)> = Vec::new();
        let mut cur_t = path[0];
        let mut start_u = 0usize;
        for (u, &t) in path.iter().enumerate().skip(1) {
            if t != cur_t {
                segs.push((start_u, u, cur_t));
                cur_t = t;
                start_u = u;
            }
        }
        segs.push((start_u, l_sub, cur_t));

        // Iterative statistical boundary pruning
        let mut changed = true;
        while changed && segs.len() > 1 {
            changed = false;
            let mut best_merge_idx: Option<usize> = None;
            let mut worst_score = -1.0f64;

            for i in 0..(segs.len() - 1) {
                let (s1, e1, t1) = segs[i];
                let (s2, e2, t2) = segs[i + 1];
                if t1 == t2 {
                    best_merge_idx = Some(i);
                    break;
                }
                let span1 = (col_to_hxb2[e1 - 1] + 1).saturating_sub(col_to_hxb2[s1]);
                let span2 = (col_to_hxb2[e2 - 1] + 1).saturating_sub(col_to_hxb2[s2]);
                let seg1_end_id = self.segment_id_at_ref(col_to_hxb2[e1 - 1]);
                let seg2_start_id = self.segment_id_at_ref(col_to_hxb2[s2]);
                let is_seg_boundary = self.spec.is_segmented
                    && seg1_end_id != 0
                    && seg2_start_id != 0
                    && seg1_end_id != seg2_start_id;
                let eff_min1 = if is_seg_boundary {
                    self.effective_segment_min_span(
                        col_to_hxb2[s1],
                        col_to_hxb2[e1 - 1],
                        min_hxb2_span_nt,
                    )
                } else {
                    min_hxb2_span_nt
                };
                let eff_min2 = if is_seg_boundary {
                    self.effective_segment_min_span(
                        col_to_hxb2[s2],
                        col_to_hxb2[e2 - 1],
                        min_hxb2_span_nt,
                    )
                } else {
                    min_hxb2_span_nt
                };

                let (pval, n11, _n12, _n21, n22, _, _) = self.test_boundary_significance_for(
                    &mat,
                    q_sub,
                    offset_in_coding,
                    ghost_idx,
                    s1,
                    e1,
                    e2,
                    t1,
                    t2,
                    min_informative_snps,
                );
                if span1 < eff_min1
                    || span2 < eff_min2
                    || pval > p_fisher_thresh
                    || n11 < min_informative_snps
                    || n22 < min_informative_snps
                {
                    let score = pval
                        + if span1 < eff_min1 || span2 < eff_min2 {
                            1.5
                        } else {
                            0.0
                        };
                    if score > worst_score {
                        worst_score = score;
                        best_merge_idx = Some(i);
                    }
                }
            }

            if let Some(i) = best_merge_idx {
                let (s1, e1, t1) = segs[i];
                let (s2, e2, t2) = segs[i + 1];
                let span1 = (col_to_hxb2[e1 - 1] + 1).saturating_sub(col_to_hxb2[s1]);
                let span2 = (col_to_hxb2[e2 - 1] + 1).saturating_sub(col_to_hxb2[s2]);

                let valid1 = Self::trunk_valid_fraction_for(
                    &mat,
                    t1,
                    offset_in_coding + s1,
                    offset_in_coding + e2,
                );
                let valid2 = Self::trunk_valid_fraction_for(
                    &mat,
                    t2,
                    offset_in_coding + s1,
                    offset_in_coding + e2,
                );

                let keep_t = if valid1 < 0.7 && valid2 >= 0.7 {
                    t2
                } else if valid2 < 0.7 && valid1 >= 0.7 {
                    t1
                } else if span1 >= span2 {
                    t1
                } else {
                    t2
                };
                segs[i] = (s1, e2, keep_t);
                segs.remove(i + 1);
                changed = true;
            }
        }

        // Merge identical adjacent states
        let mut state_merged: Vec<(usize, usize, usize)> = Vec::new();
        for seg in segs {
            if let Some(last) = state_merged.last_mut() {
                if last.2 == seg.2 {
                    last.1 = seg.1;
                    continue;
                }
            }
            state_merged.push(seg);
        }

        // Label each piece, then merge adjacent pieces whose final label is equal
        // (e.g. trunk ?E next to a ghost segment relabelled ?E).
        let label_of = |s: usize, e: usize, t: usize| -> String {
            if t < k_trunks {
                return names[t].to_string();
            }
            if allow_ghost_e_relabel {
                if let Some(ref rule) = self.spec.ghost_relabel_rule {
                    if !self.ghost_donor_coding.is_empty() {
                        let h0 = col_to_hxb2[s];
                        let mut valid_cnt = 0usize;
                        let mut match_cnt = 0usize;
                        for u in s..e {
                            let qb = q_sub[u];
                            let cb = self.ghost_donor_coding[offset_in_coding + u];
                            if qb >= 0 && cb >= 0 {
                                valid_cnt += 1;
                                if qb == cb {
                                    match_cnt += 1;
                                }
                            }
                        }
                        if valid_cnt > rule.min_span_nt
                            && h0 >= rule.min_ref_coord
                            && (match_cnt as f64 / valid_cnt as f64) >= rule.min_margin_per_site
                        {
                            return rule.ghost_trunk.clone();
                        }
                    }
                }
            }
            "Unattributed".to_string()
        };
        // (start, end, state, label, mixed states, span of the piece keeping the state)
        let mut merged: Vec<(usize, usize, usize, String, bool, usize)> = Vec::new();
        for &(s, e, t) in &state_merged {
            let label = label_of(s, e, t);
            let span = (col_to_hxb2[e - 1] + 1).saturating_sub(col_to_hxb2[s]);
            if let Some(last) = merged.last_mut() {
                if last.3 == label {
                    last.4 |= last.2 != t;
                    if span > last.5 {
                        last.2 = t;
                        last.5 = span;
                    }
                    last.1 = e;
                    continue;
                }
            }
            merged.push((s, e, t, label, false, span));
        }

        let mut breakpoints = Vec::new();
        for i in 0..(merged.len().saturating_sub(1)) {
            let (s1, e1, t1, _, _, _) = merged[i];
            let (_s2, e2, t2, _, _, _) = merged[i + 1];
            let (pval, n11, n12, n21, n22, ci_low_col, ci_high_col) = self
                .test_boundary_significance_for(
                    &mat,
                    q_sub,
                    offset_in_coding,
                    ghost_idx,
                    s1,
                    e1,
                    e2,
                    t1,
                    t2,
                    min_informative_snps,
                );
            let ci_low_h = col_to_hxb2[ci_low_col];
            let ci_high_h = col_to_hxb2[ci_high_col];
            breakpoints.push(BreakpointCall {
                bp_hxb2: col_to_hxb2[e1],
                ci_hxb2_low: ci_low_h,
                ci_hxb2_high: ci_high_h,
                ci_width_nt: ci_high_h.saturating_sub(ci_low_h) + 1,
                left_trunk: String::new(),
                right_trunk: String::new(),
                p_fisher: pval,
                n11,
                n12,
                n21,
                n22,
            });
        }

        let mut out_segs = Vec::new();
        for (i, (s, e, t, t_label, mixed, _)) in merged.iter().cloned().enumerate() {
            let h0 = col_to_hxb2[s];
            let h1 = col_to_hxb2[e - 1];
            let ci_start = if i > 0 {
                breakpoints[i - 1].ci_hxb2_low
            } else {
                h0
            };
            let ci_end = if i < breakpoints.len() {
                breakpoints[i].ci_hxb2_high
            } else {
                h1
            };

            // Compute identity % against assigned trunk
            let mut id_valid = 0usize;
            let mut id_match = 0usize;
            let t_check: Option<&[i8]> = if t < k_trunks && !mixed {
                Some(mat[t])
            } else if allow_ghost_e_relabel
                && self
                    .spec
                    .ghost_relabel_rule
                    .as_ref()
                    .map(|r| r.ghost_trunk == t_label)
                    .unwrap_or(false)
                && !self.ghost_donor_coding.is_empty()
            {
                Some(self.ghost_donor_coding.as_slice())
            } else {
                None
            };
            if let Some(ref_arr) = t_check {
                for u in s..e {
                    let qb = q_sub[u];
                    let rb = ref_arr[offset_in_coding + u];
                    if qb >= 0 && rb >= 0 {
                        id_valid += 1;
                        if qb == rb {
                            id_match += 1;
                        }
                    }
                }
            }
            let identity_pct = if id_valid > 0 {
                (100.0 * id_match as f64 / id_valid as f64 * 10.0).round() / 10.0
            } else {
                0.0
            };

            out_segs.push(SegmentCall {
                seg_idx: i + 1,
                trunk: t_label,
                hxb2_start: h0,
                hxb2_end: h1,
                span_nt: h1 - h0 + 1,
                ci_hxb2_start: ci_start,
                ci_hxb2_end: ci_end,
                short_gene: self.spec.interval_to_short_gene(h0, h1),
                identity_pct,
            });
        }

        for (i, bp) in breakpoints.iter_mut().enumerate() {
            bp.left_trunk = out_segs[i].trunk.clone();
            bp.right_trunk = out_segs[i + 1].trunk.clone();
        }

        (out_segs, breakpoints)
    }

    #[inline]
    fn site_emission_cost(&self, q_base: i8, ref_base: i8, mismatch_cost: f64) -> f64 {
        if q_base < 0 {
            return 0.0;
        }
        if ref_base == -2 {
            return 2.0;
        }
        if ref_base == -1 {
            return 0.35;
        }
        if q_base != ref_base {
            mismatch_cost
        } else {
            0.0
        }
    }

    fn trunk_valid_fraction_for(
        mat: &[&[i8]],
        t_idx: usize,
        cod_start: usize,
        cod_end: usize,
    ) -> f64 {
        if t_idx >= mat.len() || cod_end <= cod_start {
            return 1.0;
        }
        let slice = &mat[t_idx][cod_start..cod_end];
        let valid = slice.iter().filter(|&&b| b >= 0).count();
        valid as f64 / slice.len() as f64
    }

    #[allow(clippy::too_many_arguments)]
    fn test_boundary_significance_for(
        &self,
        mat: &[&[i8]],
        q_sub: &[i8],
        offset: usize,
        ghost_idx: usize,
        col_start: usize,
        col_bp: usize,
        col_end: usize,
        t_l: usize,
        t_r: usize,
        min_snps: usize,
    ) -> (f64, usize, usize, usize, usize, usize, usize) {
        let mut n11 = 0usize;
        let mut n12 = 0usize;
        let mut n21 = 0usize;
        let mut n22 = 0usize;

        if t_l != ghost_idx && t_r != ghost_idx {
            let ref_l = &mat[t_l];
            let ref_r = &mat[t_r];

            let valid_r_on_left = (col_start..col_bp)
                .filter(|&u| ref_r[offset + u] >= 0)
                .count() as f64
                / ((col_bp - col_start).max(1) as f64);

            if valid_r_on_left < 0.5 {
                for u in col_start..col_bp {
                    let qb = q_sub[u];
                    let rl = ref_l[offset + u];
                    if qb >= 0 && rl >= 0 {
                        let other_diff = mat
                            .iter()
                            .any(|tr| tr[offset + u] >= 0 && tr[offset + u] != rl);
                        if other_diff {
                            if qb == rl {
                                n11 += 1;
                            } else {
                                n12 += 1;
                            }
                        }
                    }
                }
            } else {
                for u in col_start..col_bp {
                    let qb = q_sub[u];
                    let rl = ref_l[offset + u];
                    let rr = ref_r[offset + u];
                    if qb >= 0 && rl >= 0 && rr >= 0 && rl != rr {
                        if qb == rl {
                            n11 += 1;
                        } else if qb == rr {
                            n12 += 1;
                        }
                    }
                }
            }

            let valid_l_on_right = (col_bp..col_end)
                .filter(|&u| ref_l[offset + u] >= 0)
                .count() as f64
                / ((col_end - col_bp).max(1) as f64);

            if valid_l_on_right < 0.5 {
                for u in col_bp..col_end {
                    let qb = q_sub[u];
                    let rr = ref_r[offset + u];
                    if qb >= 0 && rr >= 0 {
                        let other_diff = mat
                            .iter()
                            .any(|tr| tr[offset + u] >= 0 && tr[offset + u] != rr);
                        if other_diff {
                            if qb == rr {
                                n22 += 1;
                            } else {
                                n21 += 1;
                            }
                        }
                    }
                }
            } else {
                for u in col_bp..col_end {
                    let qb = q_sub[u];
                    let rl = ref_l[offset + u];
                    let rr = ref_r[offset + u];
                    if qb >= 0 && rl >= 0 && rr >= 0 && rl != rr {
                        if qb == rl {
                            n21 += 1;
                        } else if qb == rr {
                            n22 += 1;
                        }
                    }
                }
            }
        } else {
            let known_t = if t_l == ghost_idx { t_r } else { t_l };
            let ref_k = &mat[known_t];

            let mut match_k_l = 0usize;
            let mut orphan_l = 0usize;
            for u in col_start..col_bp {
                let qb = q_sub[u];
                let rk = ref_k[offset + u];
                if qb >= 0 && rk >= 0 {
                    if qb == rk {
                        match_k_l += 1;
                    } else if !mat.iter().any(|tr| tr[offset + u] == qb) {
                        orphan_l += 1;
                    }
                }
            }

            let mut match_k_r = 0usize;
            let mut orphan_r = 0usize;
            for u in col_bp..col_end {
                let qb = q_sub[u];
                let rk = ref_k[offset + u];
                if qb >= 0 && rk >= 0 {
                    if qb == rk {
                        match_k_r += 1;
                    } else if !mat.iter().any(|tr| tr[offset + u] == qb) {
                        orphan_r += 1;
                    }
                }
            }

            if t_l == ghost_idx {
                n11 = orphan_l;
                n12 = match_k_l;
                n21 = orphan_r;
                n22 = match_k_r;
            } else {
                n11 = match_k_l;
                n12 = orphan_l;
                n21 = match_k_r;
                n22 = orphan_r;
            }
        }

        let pval = if n11 >= min_snps && n22 >= min_snps && (n11 * n22 > n12 * n21) {
            fisher_exact_greater(n11, n12, n21, n22)
        } else {
            1.0
        };

        // Informative site CI brackets
        let mut left_info_last = col_bp;
        let mut right_info_first = col_bp;

        if t_l != ghost_idx && t_r != ghost_idx {
            let ref_l = &mat[t_l];
            let ref_r = &mat[t_r];
            for u in col_start..col_bp {
                let qb = q_sub[u];
                let rl = ref_l[offset + u];
                let rr = ref_r[offset + u];
                if qb >= 0 && rl >= 0 && rr >= 0 && rl != rr && (qb == rl || qb == rr) {
                    left_info_last = u;
                }
            }
            for u in col_bp..col_end {
                let qb = q_sub[u];
                let rl = ref_l[offset + u];
                let rr = ref_r[offset + u];
                if qb >= 0 && rl >= 0 && rr >= 0 && rl != rr && (qb == rl || qb == rr) {
                    right_info_first = u;
                    break;
                }
            }
        }

        (
            pval,
            n11,
            n12,
            n21,
            n22,
            left_info_last,
            right_info_first,
        )
    }

    fn display_lineage_name(lbl: &str) -> String {
        if lbl == "G0" {
            "G".to_string()
        } else if let Some(stripped) = lbl.strip_prefix("Subtype_") {
            stripped.to_string()
        } else {
            lbl.to_string()
        }
    }

    #[allow(clippy::too_many_arguments)]
    fn kadane_rival_scan(
        &self,
        q_coding: &[i8],
        cand_prof: &[i8],
        rivals: &[RivalProfile],
        skip_labels: &[&str],
        min_stem: i32,
        p_thresh: f64,
        allowed_trunks: Option<&[bool]>,
        temporal_limit: Option<f64>,
    ) -> Option<(String, usize, usize, i32, f64)> {
        let coding_len = self.spec.coding_len();
        let coding_start = self.spec.coding_start;
        let mut worst_tract: Option<(String, usize, usize, i32, f64)> = None;
        let mut best_p = 1.0f64;

        for rival in rivals {
            if skip_labels.iter().any(|&s| s == rival.label) {
                continue;
            }
            if let Some(lim) = temporal_limit {
                if rival.emergence_floor.is_some_and(|fl| fl > lim) {
                    continue;
                }
            }
            if let Some(allowed) = allowed_trunks {
                if rival
                    .trunk_deps
                    .iter()
                    .any(|&ti| !allowed.get(ti).copied().unwrap_or(true))
                {
                    continue;
                }
            }
            let mut cur = 0i32;
            let mut cur_start = 0usize;
            let mut best = 0i32;
            let mut best_s = 0usize;
            let mut best_e = 0usize;

            for pos in 0..coding_len {
                let qb = q_coding[pos];
                if qb < 0 {
                    continue;
                }
                let cb = cand_prof[pos];
                let rb = rival.profile[pos];
                if cb < 0 || rb < 0 || cb == rb {
                    continue;
                }
                let w = if qb == rb {
                    1
                } else if qb == cb {
                    -1
                } else {
                    0
                };
                if w == 0 {
                    continue;
                }
                if cur + w <= 0 {
                    cur = 0;
                    cur_start = pos + 1;
                } else {
                    if cur == 0 {
                        cur_start = pos;
                    }
                    cur += w;
                    if cur > best {
                        best = cur;
                        best_s = cur_start;
                        best_e = pos;
                    }
                }
            }

            if best >= min_stem && best_e >= best_s {
                let mut wins_rival = 0usize;
                let mut wins_cand = 0usize;
                for pos in best_s..=best_e {
                    let qb = q_coding[pos];
                    let cb = cand_prof[pos];
                    let rb = rival.profile[pos];
                    if qb >= 0 && cb >= 0 && rb >= 0 && cb != rb {
                        if qb == rb {
                            wins_rival += 1;
                        } else if qb == cb {
                            wins_cand += 1;
                        }
                    }
                }
                let tot = wins_rival + wins_cand;
                let net = wins_rival as i32 - wins_cand as i32;
                let span_nt = best_e - best_s + 1;
                if net >= min_stem && span_nt >= 80 && tot > 0 {
                    let p_val = binom_upper_tail_half(wins_rival, tot);
                    if p_val < p_thresh && p_val < best_p {
                        best_p = p_val;
                        worst_tract = Some((
                            rival.label.clone(),
                            coding_start + best_s,
                            coding_start + best_e,
                            net,
                            p_val,
                        ));
                    }
                }
            }
        }
        worst_tract
    }

    #[allow(clippy::too_many_arguments)]
    fn call_certified_lineage(
        &self,
        q_coding: &[i8],
        _h_start: usize,
        _h_end: usize,
        locus_scope: &str,
        basal_segs: &[SegmentCall],
        basal_bps: &[BreakpointCall],
        cfg: &ScreeningConfig,
        allowed_trunks: Option<&[bool]>,
        temporal_limit: Option<f64>,
    ) -> LineageCallDecision {
        let coding_len = self.spec.coding_len();
        let coding_start = self.spec.coding_start;

        // 1. Best pure basal trunk strictly across admissible 0..pure_basal_count (never domain-restricted ghosts)
        let mut best_basal_idx = (0..self.pure_basal_count)
            .find(|&b| allowed_trunks.is_none_or(|m| m.get(b).copied().unwrap_or(true)))
            .unwrap_or(0);
        let mut best_basal_id = -1.0f64;
        let mut best_basal_matched = 0usize;
        let mut best_basal_compared = 0usize;
        let mut basal_ids = vec![(-1.0f64, 0usize); self.pure_basal_count];
        for b in 0..self.pure_basal_count {
            if allowed_trunks.is_some_and(|m| !m.get(b).copied().unwrap_or(true)) {
                continue;
            }
            let mut matched = 0usize;
            let mut compared = 0usize;
            for (&qb, &rb) in q_coding.iter().zip(self.trunk_mat[b].iter()) {
                if qb >= 0 && rb >= 0 {
                    compared += 1;
                    if qb == rb {
                        matched += 1;
                    }
                }
            }
            if compared >= 50 {
                let id = matched as f64 / compared as f64;
                basal_ids[b] = (id, b);
                if id > best_basal_id {
                    best_basal_id = id;
                    best_basal_idx = b;
                    best_basal_matched = matched;
                    best_basal_compared = compared;
                }
            }
        }
        let best_basal_name = self.trunk_names[best_basal_idx].as_str();
        let best_basal_prof = &self.trunk_mat[best_basal_idx];

        struct SurvivingCand<'a> {
            cand: &'a CertifiedLineage,
            cand_id: f64,
            compared: usize,
            net_stem: i32,
            p_stem: f64,
        }

        let is_cand_admissible = |cand: &CertifiedLineage| -> bool {
            if let Some(lim) = temporal_limit {
                if cand.emergence_floor.is_some_and(|fl| fl > lim) {
                    return false;
                }
            }
            if let Some(allowed) = allowed_trunks {
                if !allowed.get(cand.backbone_idx).copied().unwrap_or(true)
                    || cand
                        .segment_trunk_indices
                        .iter()
                        .any(|&ti| !allowed.get(ti).copied().unwrap_or(true))
                {
                    return false;
                }
            }
            true
        };

        let mut survivors: Vec<SurvivingCand<'_>> = Vec::new();
        for cand in &self.certified_lineages {
            if !is_cand_admissible(cand) {
                continue;
            }
            // If a certified clonal descendant corresponds to a basal trunk (e.g., Subtype_G <-> G0),
            // test its stem against the closest non-self pure basal trunk.
            let is_clonal_alias = self
                .spec
                .basal_trunks
                .iter()
                .any(|t| t.name == best_basal_name && t.clonal_certified_label.as_deref() == Some(cand.label.as_str()));
            let stem_ref_idx = if is_clonal_alias {
                basal_ids
                    .iter()
                    .filter(|&&(id, b)| {
                        id >= 0.0
                            && self.trunk_names[b] != best_basal_name
                            && allowed_trunks.is_none_or(|m| m.get(b).copied().unwrap_or(true))
                    })
                    .max_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal))
                    .map(|&(_, b)| b)
                    .unwrap_or(best_basal_idx)
            } else {
                best_basal_idx
            };
            let stem_ref_prof = &self.trunk_mat[stem_ref_idx];

            let mut matched = 0usize;
            let mut compared = 0usize;
            let mut cand_wins = 0usize;
            let mut basal_wins = 0usize;

            for pos in 0..coding_len {
                let qb = q_coding[pos];
                if qb < 0 {
                    continue;
                }
                let cb = cand.profile[pos];
                if cb >= 0 {
                    compared += 1;
                    if qb == cb {
                        matched += 1;
                    }
                    let bb = stem_ref_prof[pos];
                    if bb >= 0 && cb != bb {
                        if qb == cb {
                            cand_wins += 1;
                        } else if qb == bb {
                            basal_wins += 1;
                        }
                    }
                }
            }
            if compared < self.spec.min_mapped_coding_nt {
                continue;
            }
            let cand_id = matched as f64 / compared as f64;
            if cand_id < 0.84 {
                continue;
            }

            // Guard 3: Binomial clade-stem synapomorphy test vs best of pure basal trunks
            let net_stem = cand_wins as i32 - basal_wins as i32;
            let min_net = if compared >= 4000 { 8 } else { 5 };
            let max_p_stem = if compared >= 4000 {
                cfg.p_fisher_thresh.min(1e-4)
            } else {
                cfg.p_fisher_thresh.max(0.01).min(0.01)
            };
            if net_stem < min_net {
                continue;
            }
            let p_stem = binom_upper_tail_half(cand_wins, cand_wins + basal_wins);
            if p_stem > max_p_stem {
                continue;
            }

            // Guard 1: Donor-Tract Completeness & Window Distinguishability
            if cand.category != "Pure_Subtype"
                && cand.category != "Ancestral_Recombinant_Subtype"
                && cand.segments.len() >= 2
            {
                let mut covered_segs = 0usize;
                let mut donor_ok = true;
                let mut all_cand_trunks: HashSet<&str> = HashSet::new();
                let mut covered_trunks: HashSet<&str> = HashSet::new();
                let bb_name = self.trunk_names[cand.backbone_idx].as_str();
                let bb_prof = &self.trunk_mat[cand.backbone_idx];

                for (seg_trk, s_h0, s_h1) in &cand.segments {
                    all_cand_trunks.insert(seg_trk.as_str());
                    let c0 = s_h0.saturating_sub(coding_start).min(coding_len);
                    let c1 = (s_h1 + 1).saturating_sub(coding_start).min(coding_len);
                    if c1 <= c0 {
                        continue;
                    }
                    let mut cov = 0usize;
                    let mut cw_bb = 0usize;
                    let mut bw_bb = 0usize;
                    let mut cw_best = 0usize;
                    let mut bw_best = 0usize;
                    for pos in c0..c1 {
                        let qb = q_coding[pos];
                        let cb = cand.profile[pos];
                        if qb < 0 || cb < 0 {
                            continue;
                        }
                        cov += 1;
                        let bp_b = bb_prof[pos];
                        if bp_b >= 0 && cb != bp_b {
                            if qb == cb {
                                cw_bb += 1;
                            } else if qb == bp_b {
                                bw_bb += 1;
                            }
                        }
                        let bp_best = best_basal_prof[pos];
                        if bp_best >= 0 && cb != bp_best {
                            if qb == cb {
                                cw_best += 1;
                            } else if qb == bp_best {
                                bw_best += 1;
                            }
                        }
                    }
                    if cov >= 120 {
                        covered_segs += 1;
                        covered_trunks.insert(seg_trk.as_str());
                        if seg_trk != bb_name && cw_bb + bw_bb >= 6 && cw_bb <= bw_bb {
                            donor_ok = false;
                            break;
                        }
                        if seg_trk != best_basal_name && cw_best + bw_best >= 6 && cw_best <= bw_best {
                            donor_ok = false;
                            break;
                        }
                    }
                }
                if !donor_ok {
                    continue;
                }
                if locus_scope == "full" {
                    if covered_segs < 2 {
                        continue;
                    }
                } else if !cand.global_stem && !cand.clonal_stem {
                    if covered_segs < 2 {
                        continue;
                    }
                    // Partial-amplicon Novel Donor-Trunk Completeness: every non-ancestral donor trunk must be covered
                    let mut parent_all_trunks: HashSet<&str> = HashSet::new();
                    for &pi in &cand.dag_parent_indices {
                        for (p_trk, _, _) in &self.certified_lineages[pi].segments {
                            parent_all_trunks.insert(p_trk.as_str());
                        }
                    }
                    if let Some(ref rule) = self.spec.ghost_relabel_rule {
                        if all_cand_trunks.contains(rule.partner_trunk.as_str())
                            && all_cand_trunks.contains(rule.ghost_trunk.as_str())
                        {
                            parent_all_trunks.insert(rule.ghost_trunk.as_str());
                        }
                    }
                    if all_cand_trunks
                        .iter()
                        .any(|&t| !covered_trunks.contains(t) && !parent_all_trunks.contains(t))
                    {
                        continue;
                    }
                    // Full Segment Completeness vs Global Founders: avoid promoting an incompletely covered
                    // derived/novel CRF over an established global_stem founder with the same covered trunk composition
                    if covered_segs < cand.segments.len() {
                        let shadowed_by_founder = self.certified_lineages.iter().any(|gs| {
                            if !gs.global_stem || gs.label == cand.label || !is_cand_admissible(gs) {
                                return false;
                            }
                            let mut gs_trunks: HashSet<&str> = HashSet::new();
                            for (g_trk, g_h0, g_h1) in &gs.segments {
                                let gc0 = g_h0.saturating_sub(coding_start).min(coding_len);
                                let gc1 = (g_h1 + 1).saturating_sub(coding_start).min(coding_len);
                                if gc1 <= gc0 {
                                    continue;
                                }
                                let mut gcov = 0usize;
                                for pos in gc0..gc1 {
                                    if q_coding[pos] >= 0 && gs.profile[pos] >= 0 {
                                        gcov += 1;
                                        if gcov >= 120 {
                                            gs_trunks.insert(g_trk.as_str());
                                            break;
                                        }
                                    }
                                }
                            }
                            gs_trunks == covered_trunks
                        });
                        if shadowed_by_founder {
                            continue;
                        }
                    }
                }
            }

            // Guard 2: DAG Parental Distinguishability
            let max_p_dag = if compared >= 4000 { 1e-3 } else { 0.01 };
            let mut dag_ok = true;
            for &p_idx in &cand.dag_parent_indices {
                let parent = &self.certified_lineages[p_idx];
                let mut cw_p = 0usize;
                let mut pw = 0usize;
                let mut disc_cov = 0usize;
                let mut cw_p_disc = 0usize;
                let mut pw_disc = 0usize;
                for pos in 0..coding_len {
                    let qb = q_coding[pos];
                    let cb = cand.profile[pos];
                    let pb = parent.profile[pos];
                    if qb >= 0 && cb >= 0 && pb >= 0 {
                        let is_disc_seg = locus_scope != "full"
                            && (segment_trunk_at(&cand.segments, coding_start + pos)
                                != segment_trunk_at(&parent.segments, coding_start + pos));
                        if is_disc_seg {
                            disc_cov += 1;
                        }
                        if cb != pb {
                            if qb == cb {
                                cw_p += 1;
                                if is_disc_seg {
                                    cw_p_disc += 1;
                                }
                            } else if qb == pb {
                                pw += 1;
                                if is_disc_seg {
                                    pw_disc += 1;
                                }
                            }
                        }
                    }
                }
                let net_p = cw_p as i32 - pw as i32;
                if net_p < 5 || binom_upper_tail_half(cw_p, cw_p + pw) > max_p_dag {
                    dag_ok = false;
                    break;
                }
                if locus_scope != "full" && (disc_cov < 120 || (cw_p_disc as i32 - pw_disc as i32) < 5) {
                    dag_ok = false;
                    break;
                }
            }
            if !dag_ok {
                continue;
            }

            survivors.push(SurvivingCand {
                cand,
                cand_id,
                compared,
                net_stem,
                p_stem,
            });
        }

        survivors.sort_by(|a, b| {
            b.cand_id
                .partial_cmp(&a.cand_id)
                .unwrap_or(std::cmp::Ordering::Equal)
                .then_with(|| b.net_stem.cmp(&a.net_stem))
                .then_with(|| b.cand.global_stem.cmp(&a.cand.global_stem))
                .then_with(|| a.cand.label.cmp(&b.cand.label))
        });

        if let Some(best) = survivors.first() {
            let min_stem_g4 = if best.compared >= 4000 { 12 } else { 10 };
            let mut skip_vec: Vec<&str> = vec![best.cand.label.as_str()];
            for t in &self.spec.basal_trunks {
                if t.clonal_certified_label.as_deref() == Some(best.cand.label.as_str()) {
                    skip_vec.push(t.name.as_str());
                }
            }
            if let Some((rival_lbl, r_h0, r_h1, net_r, p_r)) = self.kadane_rival_scan(
                q_coding,
                &best.cand.profile,
                &self.founder_rivals,
                &skip_vec,
                min_stem_g4,
                1e-4,
                allowed_trunks,
                temporal_limit,
            ) {
                let disp_base = self.spec.display_lineage_name(&best.cand.label);
                let disp_rival = self.spec.display_lineage_name(&rival_lbl);
                return LineageCallDecision {
                    called_lineage: "URF".to_string(),
                    called_category: "URF".to_string(),
                    certified_compendium_label: None,
                    call_support: LineageCallSupport {
                        identity_pct: (best.cand_id * 1000.0).round() / 10.0,
                        compared_nt: best.compared,
                        net_stem_snps: -net_r,
                        p_value: p_r,
                        note: format!(
                            "Rejected {} due to residual {} tract ({}-{} {}, net {:+} SNPs, p={:.1e})",
                            disp_base, disp_rival, r_h0, r_h1, self.spec.reference_name, net_r, p_r
                        ),
                    },
                    residual_tract: Some((disp_base, disp_rival, r_h0, r_h1)),
                };
            }

            let disp = self.spec.display_lineage_name(&best.cand.label);
            let stem_note = if best.cand.global_stem {
                format!("Global founder stem vs {} (+{} SNPs)", best_basal_name, best.net_stem)
            } else if best.cand.clonal_stem {
                format!("Clonal clade stem vs {} (+{} SNPs)", best_basal_name, best.net_stem)
            } else {
                format!("Certified mosaic clade vs {} (+{} SNPs)", best_basal_name, best.net_stem)
            };
            return LineageCallDecision {
                called_lineage: disp,
                called_category: best.cand.category.clone(),
                certified_compendium_label: Some(best.cand.label.clone()),
                call_support: LineageCallSupport {
                    identity_pct: (best.cand_id * 1000.0).round() / 10.0,
                    compared_nt: best.compared,
                    net_stem_snps: best.net_stem,
                    p_value: best.p_stem,
                    note: stem_note,
                },
                residual_tract: None,
            };
        }

        // 3. Fallback: Pure Basal or Unindexed URF
        if basal_segs.len() >= 2 {
            let min_bp_p = basal_bps
                .iter()
                .map(|b| b.p_fisher)
                .fold(1.0f64, |acc, p| acc.min(p));
            return LineageCallDecision {
                called_lineage: "URF".to_string(),
                called_category: "URF".to_string(),
                certified_compendium_label: None,
                call_support: LineageCallSupport {
                    identity_pct: (best_basal_id.max(0.0) * 1000.0).round() / 10.0,
                    compared_nt: best_basal_compared,
                    net_stem_snps: 0,
                    p_value: min_bp_p,
                    note: format!(
                        "Multi-trunk mosaic outside all {} certified profiles",
                        self.certified_lineages.len()
                    ),
                },
                residual_tract: None,
            };
        }

        let only_trunk = basal_segs
            .first()
            .map(|s| s.trunk.as_str())
            .unwrap_or("Unattributed");
        if only_trunk == "Unattributed" || only_trunk.starts_with('?') {
            return LineageCallDecision {
                called_lineage: "Unclassifiable".to_string(),
                called_category: "Unclassifiable".to_string(),
                certified_compendium_label: None,
                call_support: LineageCallSupport {
                    identity_pct: (best_basal_id.max(0.0) * 1000.0).round() / 10.0,
                    compared_nt: best_basal_compared,
                    net_stem_snps: 0,
                    p_value: 1.0,
                    note: format!("Unclassifiable basal segment ({})", only_trunk),
                },
                residual_tract: None,
            };
        }

        let single_idx = self
            .trunk_names
            .iter()
            .position(|t| t == only_trunk)
            .unwrap_or(best_basal_idx);
        let mut s_matched = 0usize;
        let mut s_compared = 0usize;
        for (&qb, &rb) in q_coding.iter().zip(self.trunk_mat[single_idx].iter()) {
            if qb >= 0 && rb >= 0 {
                s_compared += 1;
                if qb == rb {
                    s_matched += 1;
                }
            }
        }
        let s_id = if s_compared > 0 {
            s_matched as f64 / s_compared as f64
        } else {
            0.0
        };
        let min_stem_g4 = if s_compared >= 4000 { 12 } else { 10 };
        let basal_rivals = &self.founder_rivals[..self.trunk_names.len()];
        let skip_vec = [only_trunk];
        if let Some((rival_lbl, r_h0, r_h1, net_r, p_r)) = self.kadane_rival_scan(
            q_coding,
            &self.trunk_mat[single_idx],
            basal_rivals,
            &skip_vec,
            min_stem_g4,
            1e-4,
            allowed_trunks,
            temporal_limit,
        ) {
            let disp_base = self.spec.display_lineage_name(only_trunk);
            let disp_rival = self.spec.display_lineage_name(&rival_lbl);
            return LineageCallDecision {
                called_lineage: "URF".to_string(),
                called_category: "URF".to_string(),
                certified_compendium_label: None,
                call_support: LineageCallSupport {
                    identity_pct: (s_id * 1000.0).round() / 10.0,
                    compared_nt: s_compared,
                    net_stem_snps: -net_r,
                    p_value: p_r,
                    note: format!(
                        "Residual {} tract in {} background ({}-{} {}, net {:+} SNPs, p={:.1e})",
                        disp_rival, disp_base, r_h0, r_h1, self.spec.reference_name, net_r, p_r
                    ),
                },
                residual_tract: Some((disp_base, disp_rival, r_h0, r_h1)),
            };
        }

        let disp_pure = self.spec.display_lineage_name(only_trunk);
        let comp_lbl = self
            .spec
            .basal_trunks
            .iter()
            .find(|t| t.name == only_trunk)
            .and_then(|t| t.clonal_certified_label.clone());
        let final_cat = "Pure_Subtype".to_string();
        LineageCallDecision {
            called_lineage: disp_pure.clone(),
            called_category: final_cat,
            certified_compendium_label: comp_lbl,
            call_support: LineageCallSupport {
                identity_pct: (s_id * 1000.0).round() / 10.0,
                compared_nt: s_compared,
                net_stem_snps: best_basal_matched as i32,
                p_value: 0.0,
                note: format!("Basal pure subtype {} with no residual breakpoint", disp_pure),
            },
            residual_tract: None,
        }
    }

    fn prepend_certified_match(
        &self,
        matches: &mut Vec<LineageMatch>,
        comp_label: &str,
        support: &LineageCallSupport,
    ) {
        let Some(entry) = self.compendium.iter().find(|e| e.label == comp_label) else {
            return;
        };
        let cert_match = LineageMatch {
            label: entry.label.clone(),
            descriptive_name: entry.descriptive_name.clone(),
            status: entry.status.clone(),
            match_type: format!(
                "Certified Sequence Match ({:.1}% id, {:+} stem SNPs, p={:.1e})",
                support.identity_pct, support.net_stem_snps, support.p_value
            ),
            hxb2_barcode: entry.hxb2_barcode.clone(),
        };
        matches.retain(|m| {
            !(m.label == entry.label && !m.match_type.starts_with("Homologous"))
        });
        matches.insert(0, cert_match);
        matches.truncate(6);
    }

    #[allow(clippy::too_many_arguments)]
    fn build_epidemiological_map(
        &self,
        decision: &LineageCallDecision,
        locus_scope: &str,
        h_start: usize,
        h_end: usize,
        q_sub: &[i8],
        first_idx: usize,
        col_to_hxb2: &[usize],
        basal_segs: &[SegmentCall],
        cfg: &ScreeningConfig,
        allowed_trunks: Option<&[bool]>,
        temporal_limit: Option<f64>,
    ) -> (String, String, Vec<SegmentCall>) {
        let scope_suffix = if locus_scope == "full" {
            String::new()
        } else {
            format!("@{}", locus_scope)
        };
        if decision.called_category != "URF" {
            let lin = decision.called_lineage.clone();
            let formula = format!("{}(Clonal){}", lin, scope_suffix);
            let barcode = format!("{}:{}-{}", lin, h_start, h_end);
            let segs = vec![SegmentCall {
                seg_idx: 1,
                trunk: lin,
                hxb2_start: h_start,
                hxb2_end: h_end,
                span_nt: h_end - h_start + 1,
                ci_hxb2_start: h_start,
                ci_hxb2_end: h_end,
                short_gene: self.spec.interval_to_short_gene(h_start, h_end),
                identity_pct: decision.call_support.identity_pct,
            }];
            return (formula, barcode, segs);
        }

        // URF: run epidemiological Viterbi segmentation (respecting temporal admissibility)
        let epi_allowed: Option<Vec<bool>> = if allowed_trunks.is_some() || temporal_limit.is_some() {
            let mut mask: Vec<bool> = (0..self.epi_names.len())
                .map(|i| {
                    let time_ok = temporal_limit
                        .is_none_or(|lim| self.epi_floors[i].is_none_or(|fl| fl <= lim));
                    let trunk_ok = allowed_trunks.is_none_or(|m| {
                        self.epi_trunk_deps[i]
                            .iter()
                            .all(|&ti| m.get(ti).copied().unwrap_or(true))
                    });
                    time_ok && trunk_ok
                })
                .collect();
            if mask.iter().filter(|&&a| a).count() < self.epi_names.len().min(2) {
                for (i, slot) in mask.iter_mut().enumerate() {
                    if self.epi_trunk_deps[i]
                        .iter()
                        .all(|&ti| allowed_trunks.is_none_or(|m| m.get(ti).copied().unwrap_or(true)))
                    {
                        *slot = true;
                    }
                }
            }
            Some(mask)
        } else {
            None
        };

        let gamma = cfg.gamma.unwrap_or(7.0);
        let min_span_nt = cfg.min_span_nt.unwrap_or(if locus_scope == "full" {
            120
        } else {
            150
        });
        let sub_len = q_sub.len();
        let (mut epi_segs, _) = self.run_paff_viterbi_panel(
            &self.epi_names,
            &self.epi_mat,
            &self.epi_penalties,
            epi_allowed.as_deref(),
            false,
            q_sub,
            first_idx,
            col_to_hxb2,
            sub_len,
            gamma,
            min_span_nt,
            cfg.min_informative_snps,
            cfg.p_fisher_thresh,
        );

        // If default gamma kept only 1 epi segment, try relaxed gamma/span or splice residual tract
        if epi_segs.len() < 2 {
            let (relaxed_segs, _) = self.run_paff_viterbi_panel(
                &self.epi_names,
                &self.epi_mat,
                &self.epi_penalties,
                epi_allowed.as_deref(),
                false,
                q_sub,
                first_idx,
                col_to_hxb2,
                sub_len,
                3.5,
                80,
                4,
                0.01,
            );
            if relaxed_segs.len() >= 2 {
                epi_segs = relaxed_segs;
            }
        }

        if epi_segs.len() < 2 {
            if let Some((base_lbl, rival_lbl, r_h0, r_h1)) = &decision.residual_tract {
                let s0 = (*r_h0).max(h_start);
                let s1 = (*r_h1).min(h_end);
                if s1 > s0 {
                    let mut spliced: Vec<(String, usize, usize)> = Vec::new();
                    if s0.saturating_sub(h_start) >= 60 {
                        spliced.push((base_lbl.clone(), h_start, s0 - 1));
                        spliced.push((rival_lbl.clone(), s0, s1));
                    } else {
                        spliced.push((rival_lbl.clone(), h_start, s1));
                    }
                    if h_end.saturating_sub(s1) >= 60 {
                        spliced.push((base_lbl.clone(), s1 + 1, h_end));
                    } else if let Some(last) = spliced.last_mut() {
                        last.2 = h_end;
                    }
                    if spliced.len() >= 2 {
                        epi_segs = spliced
                            .into_iter()
                            .enumerate()
                            .map(|(i, (trk, a, b))| SegmentCall {
                                seg_idx: i + 1,
                                trunk: trk,
                                hxb2_start: a,
                                hxb2_end: b,
                                span_nt: b - a + 1,
                                ci_hxb2_start: a,
                                ci_hxb2_end: b,
                                short_gene: self.spec.interval_to_short_gene(a, b),
                                identity_pct: decision.call_support.identity_pct,
                            })
                            .collect();
                    }
                }
            } else if basal_segs.len() >= 2 {
                epi_segs = basal_segs
                    .iter()
                    .enumerate()
                    .map(|(i, s)| {
                        let mut c = s.clone();
                        c.seg_idx = i + 1;
                        if let Some((_, _, disp)) = self
                            .spec
                            .epi_certified_replacements
                            .iter()
                            .find(|(tr, _, _)| tr == &c.trunk)
                        {
                            c.trunk = disp.clone();
                        }
                        c
                    })
                    .collect();
            }
        }

        let (epi_desc, _bb, _sf) = format_descriptive_name(&epi_segs, locus_scope);
        let epi_bc = epi_segs
            .iter()
            .map(|s| format!("{}:{}-{}", s.trunk, s.hxb2_start, s.hxb2_end))
            .collect::<Vec<_>>()
            .join("|");
        (epi_desc, epi_bc, epi_segs)
    }

    /// Match a sequence's structural formula, segment coordinates, and breakpoints against the compendium.
    fn match_compendium(
        &self,
        struct_formula: &str,
        segs: &[SegmentCall],
        bps: &[BreakpointCall],
        locus_scope: &str,
        allowed_trunks: Option<&[bool]>,
        temporal_limit: Option<f64>,
    ) -> Vec<LineageMatch> {
        let mut matches = Vec::new();
        let base_formula = struct_formula.split('@').next().unwrap_or(struct_formula);
        let entry_ok = |entry: &CompendiumEntry| -> bool {
            if let Some(lim) = temporal_limit {
                if !self.is_compendium_entry_admissible(entry, allowed_trunks.unwrap_or(&[]), lim) {
                    return false;
                }
            } else if let Some(allowed) = allowed_trunks {
                if !self.is_compendium_entry_admissible(entry, allowed, f64::INFINITY) {
                    return false;
                }
            }
            true
        };

        // 1. Exact structural formula matches
        for entry in &self.compendium {
            if entry.category == "Demoted_LANL_Entry"
                || entry.status.starts_with("Collapsed_to_Clonal_Trunk")
                || !entry_ok(entry)
            {
                continue;
            }
            let clonal_alias_match = self.spec.basal_trunks.iter().any(|t| {
                format!("{}(Clonal)", t.name) == base_formula
                    && t.clonal_certified_label.as_deref() == Some(entry.label.as_str())
            });
            if entry.structural_formula == base_formula || clonal_alias_match {
                if (locus_scope == "full" && entry.scope == "full")
                    || (locus_scope == entry.scope)
                    || segs.len() == 1
                {
                    matches.push(LineageMatch {
                        label: entry.label.clone(),
                        descriptive_name: entry.descriptive_name.clone(),
                        status: entry.status.clone(),
                        match_type: if entry.structural_formula == base_formula {
                            "Exact Structural Formula".to_string()
                        } else {
                            "Basal Trunk Match (G0)".to_string()
                        },
                        hxb2_barcode: entry.hxb2_barcode.clone(),
                    });
                }
            }
        }

        // 2. Check coordinate-level mosaic concordance and breakpoint proximity (+/- 200 nt) against compendium barcodes
        if !bps.is_empty() && matches.len() < 4 {
            let mut scored_candidates: Vec<(f64, usize, usize, usize, &CompendiumEntry)> = Vec::new();

            for entry in &self.compendium {
                if entry.category == "Demoted_LANL_Entry"
                    || entry.status.starts_with("Collapsed_to_Clonal_Trunk")
                    || !entry_ok(entry)
                {
                    continue;
                }
                if matches.iter().any(|m| m.label == entry.label) {
                    continue;
                }
                let entry_segs = parse_barcode_segments(&entry.hxb2_barcode);
                let entry_bps = parse_barcode_breakpoints(&entry.hxb2_barcode);
                if entry_segs.is_empty() || entry_bps.is_empty() {
                    continue;
                }

                // Compute coordinate overlap concordance across query segments
                let mut total_ovl = 0usize;
                let mut same_trunk_ovl = 0usize;
                for q_seg in segs {
                    for (e_trk, e_s, e_e) in &entry_segs {
                        let ovl_s = q_seg.hxb2_start.max(*e_s);
                        let ovl_e = q_seg.hxb2_end.min(*e_e);
                        if ovl_e >= ovl_s {
                            let len = ovl_e - ovl_s + 1;
                            total_ovl += len;
                            if &q_seg.trunk == e_trk {
                                same_trunk_ovl += len;
                            }
                        }
                    }
                }
                if total_ovl < 400 {
                    continue;
                }
                let concordance = same_trunk_ovl as f64 / total_ovl as f64;

                let q_min_h = segs.first().map(|s| s.hxb2_start).unwrap_or(self.spec.coding_start);
                let q_max_h = segs.last().map(|s| s.hxb2_end).unwrap_or(self.spec.coding_end);
                let span_ebps = entry_bps
                    .iter()
                    .filter(|(_, _, pos)| *pos >= q_min_h.saturating_sub(150) && *pos <= q_max_h + 150)
                    .count()
                    .max(1);

                let mut shared = 0usize;
                let mut pos_shared = 0usize;
                for bp in bps {
                    if entry_bps.iter().any(|(l_tr, r_tr, pos)| {
                        l_tr == &bp.left_trunk
                            && r_tr == &bp.right_trunk
                            && (*pos as isize - bp.bp_hxb2 as isize).unsigned_abs() <= 200
                    }) {
                        shared += 1;
                    }
                    if entry_bps.iter().any(|(l_tr, r_tr, pos)| {
                        (l_tr == &bp.left_trunk || r_tr == &bp.right_trunk)
                            && (*pos as isize - bp.bp_hxb2 as isize).unsigned_abs() <= 200
                    }) {
                        pos_shared += 1;
                    }
                }

                if (shared == bps.len() && shared == span_ebps)
                    || (concordance >= 0.88 && shared >= 1)
                    || (concordance >= 0.91 && pos_shared >= 2)
                    || (shared >= 3 && shared * 4 >= span_ebps * 3)
                {
                    scored_candidates.push((concordance, shared, pos_shared, span_ebps, entry));
                }
            }

            scored_candidates.sort_by(|a, b| {
                b.0.partial_cmp(&a.0)
                    .unwrap_or(std::cmp::Ordering::Equal)
                    .then_with(|| b.1.cmp(&a.1))
                    .then_with(|| b.2.cmp(&a.2))
                    .then_with(|| a.4.label.cmp(&b.4.label))
            });

            for (concordance, shared, _pos_shared, span_ebps, entry) in scored_candidates.into_iter().take(4) {
                let match_type = if shared == bps.len() && shared == span_ebps {
                    format!(
                        "Homologous Breakpoints ({:.1}% coordinate concordance)",
                        concordance * 100.0
                    )
                } else {
                    format!(
                        "Homologous Mosaic ({:.1}% concordance, {}/{} breakpoints)",
                        concordance * 100.0,
                        shared,
                        span_ebps.max(shared)
                    )
                };
                matches.push(LineageMatch {
                    label: entry.label.clone(),
                    descriptive_name: entry.descriptive_name.clone(),
                    status: entry.status.clone(),
                    match_type,
                    hxb2_barcode: entry.hxb2_barcode.clone(),
                });
            }
        }

        if matches.is_empty() && bps.is_empty() {
            matches.push(LineageMatch {
                label: "Unindexed_Clonal_Lineage".to_string(),
                descriptive_name: struct_formula.to_string(),
                status: "Clonal call with no compendium entry for this trunk/scope".to_string(),
                match_type: "No Compendium Match".to_string(),
                hxb2_barcode: String::new(),
            });
        } else if matches.is_empty() {
            matches.push(LineageMatch {
                label: "Novel_URF_or_Unindexed_Mosaic".to_string(),
                descriptive_name: struct_formula.to_string(),
                status: "Unique or Unindexed Recombinant Architecture".to_string(),
                match_type: "De Novo".to_string(),
                hxb2_barcode: String::new(),
            });
        }

        matches.truncate(6);
        matches
    }
}

fn fmt_nt(n: usize) -> String {
    let s = n.to_string();
    let mut out = String::with_capacity(s.len() + s.len() / 3);
    for (i, ch) in s.chars().enumerate() {
        if i > 0 && (s.len() - i) % 3 == 0 {
            out.push(',');
        }
        out.push(ch);
    }
    out
}

fn parse_barcode_segments(barcode: &str) -> Vec<(String, usize, usize)> {
    let mut parsed_segs: Vec<(String, usize, usize)> = Vec::new();
    for part in barcode.split('|') {
        if let Some((trk, coords)) = part.split_once(':') {
            let coords_clean = coords.split('[').next().unwrap_or(coords);
            if let Some((s_str, e_str)) = coords_clean.split_once('-') {
                if let (Ok(s), Ok(e)) = (s_str.parse::<usize>(), e_str.parse::<usize>()) {
                    parsed_segs.push((trk.to_string(), s, e));
                }
            }
        }
    }
    parsed_segs
}

fn parse_barcode_breakpoints(barcode: &str) -> Vec<(String, String, usize)> {
    let parsed_segs = parse_barcode_segments(barcode);
    let mut bps = Vec::new();
    for i in 0..(parsed_segs.len().saturating_sub(1)) {
        bps.push((
            parsed_segs[i].0.clone(),
            parsed_segs[i + 1].0.clone(),
            parsed_segs[i + 1].1,
        ));
    }
    bps
}

/// Map an HXB2 nucleotide interval `[h0, h1]` to a concise HIV-1 gene/domain span label.
pub fn hxb2_to_short_gene(h0: usize, h1: usize) -> String {
    let span = h1.saturating_sub(h0);
    let pt_gene = |pos: usize| -> &'static str {
        if pos < 790 {
            "5'LTR"
        } else if pos <= 2292 {
            if pos <= 1185 {
                "gag.MA"
            } else if pos <= 1878 {
                "gag.CA"
            } else {
                "gag.NC"
            }
        } else if pos <= 5096 {
            if pos <= 2549 {
                "PR"
            } else if pos <= 4229 {
                "RT"
            } else {
                "IN"
            }
        } else if pos <= 5559 {
            "vif"
        } else if pos <= 5850 {
            "vpr"
        } else if pos <= 6062 {
            "tat1"
        } else if pos <= 6225 {
            "vpu"
        } else if pos <= 8795 {
            if pos <= 7757 {
                "gp120"
            } else {
                "gp41"
            }
        } else if pos <= 9417 {
            "nef"
        } else {
            "3'LTR"
        }
    };

    let g0 = pt_gene((h0 + 15).min(h1));
    let g1 = pt_gene(h1.saturating_sub(15).max(h0));

    if g0 == g1 {
        if (g0 == "PR" || g0 == "RT" || g0 == "IN") && span > 1800 {
            return "pol".to_string();
        }
        if (g0 == "gp120" || g0 == "gp41") && span > 1800 {
            return "env".to_string();
        }
        return g0.to_string();
    }
    if g0.starts_with("gag.") && g1.starts_with("gag.") {
        if span > 900 {
            return "gag".to_string();
        }
        return format!("{}..{}", g0, g1);
    }
    if g0 == "PR" && g1 == "IN" {
        if span > 1800 {
            return "pol".to_string();
        }
        return "PR..IN".to_string();
    }
    if g0 == "gp120" && g1 == "gp41" {
        if span > 1800 {
            return "env".to_string();
        }
        return "gp120..gp41".to_string();
    }
    let g0_norm = if g0.starts_with("gag.") && !g1.starts_with("gag.") {
        "gag"
    } else {
        g0
    };
    let g1_norm = if g1.starts_with("gag.") && !g0.starts_with("gag.") {
        "gag"
    } else {
        g1
    };
    format!("{}..{}", g0_norm, g1_norm)
}

fn format_descriptive_name(
    segs: &[SegmentCall],
    locus_scope: &str,
) -> (String, String, String) {
    let mut trunk_spans: HashMap<String, usize> = HashMap::new();
    for s in segs {
        *trunk_spans.entry(s.trunk.clone()).or_insert(0) += s.span_nt;
    }

    let mut known_spans: Vec<(String, usize)> = trunk_spans
        .iter()
        .filter(|(k, _)| !k.starts_with('?') && *k != "Unattributed")
        .map(|(k, &v)| (k.clone(), v))
        .collect();
    known_spans.sort_by(|a, b| b.1.cmp(&a.1).then_with(|| a.0.cmp(&b.0)));

    let backbone = if let Some((tr, _)) = known_spans.first() {
        tr.clone()
    } else {
        let mut all_spans: Vec<(String, usize)> =
            trunk_spans.iter().map(|(k, &v)| (k.clone(), v)).collect();
        all_spans.sort_by(|a, b| b.1.cmp(&a.1).then_with(|| a.0.cmp(&b.0)));
        all_spans
            .first()
            .map(|x| x.0.clone())
            .unwrap_or_else(|| "Unknown".to_string())
    };

    let scope_suffix = if locus_scope == "full" {
        String::new()
    } else {
        format!("@{}", locus_scope)
    };

    if trunk_spans.len() == 1 {
        let base = if locus_scope == "full" {
            format!("{}(Clonal)", backbone)
        } else {
            format!("{}(Clonal){}", backbone, scope_suffix)
        };
        return (base.clone(), backbone, base);
    }

    let mut donor_order: Vec<String> = Vec::new();
    let mut donor_genes: HashMap<String, Vec<String>> = HashMap::new();
    for s in segs {
        if s.trunk != backbone {
            if !donor_order.contains(&s.trunk) {
                donor_order.push(s.trunk.clone());
            }
            let entry = donor_genes.entry(s.trunk.clone()).or_default();
            if entry.last() != Some(&s.short_gene) {
                entry.push(s.short_gene.clone());
            }
        }
    }

    let donor_parts: Vec<String> = donor_order
        .into_iter()
        .map(|d| {
            let gs = &donor_genes[&d];
            format!("{}:{}", d, gs.join("+"))
        })
        .collect();

    let structural_formula = format!("{}[{}]", backbone, donor_parts.join(" | "));
    let descriptive_name = format!("{}{}", structural_formula, scope_suffix);
    (descriptive_name, backbone, structural_formula)
}

#[cfg(test)]
mod naming_tests {
    use super::*;

    fn seg(trunk: &str, start: usize, end: usize) -> SegmentCall {
        SegmentCall {
            seg_idx: 0,
            trunk: trunk.to_string(),
            hxb2_start: start,
            hxb2_end: end,
            span_nt: end - start + 1,
            ci_hxb2_start: start,
            ci_hxb2_end: end,
            short_gene: "pol".to_string(),
            identity_pct: 100.0,
        }
    }

    #[test]
    fn backbone_tiebreak_is_deterministic() {
        let a = seg("Unattributed", 800, 999);
        let b = seg("?E", 1000, 1199);
        for _ in 0..32 {
            let (_, bb1, _) = format_descriptive_name(&[a.clone(), b.clone()], "full");
            let (_, bb2, _) = format_descriptive_name(&[b.clone(), a.clone()], "full");
            assert_eq!(bb1, "?E");
            assert_eq!(bb2, "?E");
        }
    }
}

#[cfg(test)]
mod pathogen_pack_tests {
    use super::*;
    use std::sync::atomic::{AtomicUsize, Ordering};

    static TEMP_DIR_COUNTER: AtomicUsize = AtomicUsize::new(1);

    fn temp_pack_dir(prefix: &str) -> std::path::PathBuf {
        let id = TEMP_DIR_COUNTER.fetch_add(1, Ordering::Relaxed);
        let dir = std::env::temp_dir().join(format!(
            "eggman_{}_{}_{}",
            prefix,
            std::process::id(),
            id
        ));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        dir
    }

    fn tsv_profile_col(tsv: &str, label: &str, col_idx: usize) -> String {
        for line in tsv.lines() {
            let trimmed = line.trim();
            if trimmed.is_empty() || trimmed.starts_with('#') {
                continue;
            }
            let cols: Vec<&str> = trimmed.split('\t').collect();
            if cols.first().map(|s| s.trim()) == Some(label) {
                return cols[col_idx].trim().to_string();
            }
        }
        panic!("Missing profile '{}' in TSV table", label);
    }

    fn deterministic_dna(len: usize, seed: u64) -> Vec<u8> {
        const BASES: [u8; 4] = [b'A', b'C', b'G', b'T'];
        let mut x = seed;
        let mut seq = Vec::with_capacity(len);
        for _ in 0..len {
            x ^= x << 13;
            x ^= x >> 7;
            x ^= x << 17;
            seq.push(BASES[(x as usize) & 3]);
        }
        seq
    }

    fn mutate_periodic(anc: &[u8], period: usize, phase: usize, shift: usize) -> Vec<u8> {
        const BASES: [u8; 4] = [b'A', b'C', b'G', b'T'];
        anc.iter()
            .enumerate()
            .map(|(i, &b)| {
                if i % period == phase {
                    let cur = match b {
                        b'A' => 0,
                        b'C' => 1,
                        b'G' => 2,
                        _ => 3,
                    };
                    BASES[(cur + shift) & 3]
                } else {
                    b
                }
            })
            .collect()
    }

    fn simple_trunk(name: &str) -> TrunkSpec {
        TrunkSpec {
            name: name.to_string(),
            is_pure_basal: true,
            penalty: 0.0,
            active_range: None,
            fallback_profile: None,
            display_alias: None,
            clonal_certified_label: None,
            t_mrca: None,
            t_mrca_ci: None,
            first_sampled_year: None,
            chronaeon_half_life_yr: None,
        }
    }

    #[test]
    fn hiv1_default_pack_regression() {
        let engine = EggmanEngine::new();
        assert_eq!(engine.pathogen_name(), "HIV-1 Group M");
        assert_eq!(engine.reference_name(), "HXB2");
        assert_eq!(engine.basal_trunk_count(), 14);
        assert_eq!(engine.certified_lineage_count(), 155);

        // JSON round-trip of default PathogenSpec produces an equivalent engine
        let spec_json = serde_json::to_string_pretty(&PathogenSpec::hiv1_group_m_default()).unwrap();
        let parsed_spec: PathogenSpec = serde_json::from_str(&spec_json).unwrap();
        let mut rt_pack = PathogenProfilePack::hiv1_group_m_default();
        rt_pack.spec = parsed_spec;
        let rt_engine = EggmanEngine::from_pathogen_pack(rt_pack).unwrap();
        assert_eq!(rt_engine.pathogen_name(), "HIV-1 Group M");
        assert_eq!(rt_engine.basal_trunk_count(), 14);
        assert_eq!(rt_engine.certified_lineage_count(), 155);

        let cfg = ScreeningConfig::default();

        let seq_b = tsv_profile_col(PROFILES_TSV, "B", 1);
        let res_b = rt_engine.classify_sequence("Q_B", &seq_b, &cfg);
        assert_eq!(res_b.called_lineage, "B");
        assert_eq!(res_b.called_category, "Pure_Subtype");
        assert_eq!(res_b.structural_formula, "B(Clonal)");
        assert!(res_b.is_clonal);
        assert!(res_b.matches.iter().any(|m| m.label == "Subtype_B"));

        let seq_c = tsv_profile_col(PROFILES_TSV, "C", 1);
        let res_c = rt_engine.classify_sequence("Q_C", &seq_c, &cfg);
        assert_eq!(res_c.called_lineage, "C");
        assert_eq!(res_c.called_category, "Pure_Subtype");
        assert_eq!(res_c.structural_formula, "C(Clonal)");
        assert!(res_c.is_clonal);
        assert!(res_c.matches.iter().any(|m| m.label == "Subtype_C"));
        assert!(!res_c
            .matches
            .iter()
            .any(|m| matches!(m.label.as_str(), "CRF31_BC" | "CRF92_C2U" | "CRF146_BC")));

        let seq_g0 = tsv_profile_col(PROFILES_TSV, "G0", 1);
        let res_g0 = rt_engine.classify_sequence("Q_G0", &seq_g0, &cfg);
        assert_eq!(res_g0.called_lineage, "G");
        assert_eq!(res_g0.called_category, "Ancestral_Recombinant_Subtype");
        assert_eq!(res_g0.structural_formula, "G0(Clonal)");

        let seq_g = tsv_profile_col(CERTIFIED_PROFILES_TSV, "Subtype_G", 5);
        let res_g = rt_engine.classify_sequence("Q_G", &seq_g, &cfg);
        assert_eq!(res_g.called_lineage, "G");
        assert_eq!(res_g.called_category, "Ancestral_Recombinant_Subtype");
        assert_eq!(res_g.structural_formula, "G0[A1:gp41..nef]");

        let seq_01 = tsv_profile_col(CERTIFIED_PROFILES_TSV, "CRF01_AE", 5);
        let res_01 = rt_engine.classify_sequence("Q_CRF01", &seq_01, &cfg);
        assert_eq!(res_01.called_lineage, "CRF01_AE");
        assert_eq!(res_01.called_category, "Certified_CRF");

        let seq_02 = tsv_profile_col(CERTIFIED_PROFILES_TSV, "CRF02_AG", 5);
        let res_02 = rt_engine.classify_sequence("Q_CRF02", &seq_02, &cfg);
        assert_eq!(res_02.called_lineage, "CRF02_AG");
        assert_eq!(res_02.called_category, "Certified_CRF");
    }

    #[test]
    fn influenza_a_segmented_reassortment_pack() {
        // Concatenated 8-segment Influenza A coordinate frame (13,588 nt total)
        let segs = [
            ("PB2", 1usize, 2341usize),
            ("PB1", 2342, 4682),
            ("PA", 4683, 6915),
            ("HA", 6916, 8693),
            ("NP", 8694, 10258),
            ("NA", 10259, 11671),
            ("MP", 11672, 12698),
            ("NS", 12699, 13588),
        ];
        let total_len = 13_588usize;
        let landmarks: Vec<GeneLandmark> = segs
            .iter()
            .map(|&(name, start, end)| GeneLandmark {
                short_name: name.to_string(),
                parent_gene: None,
                parent_min_span: 0,
                start,
                end,
            })
            .collect();

        let spec = PathogenSpec {
            pathogen_name: "Influenza A".to_string(),
            reference_name: "FluA-8seg".to_string(),
            ref_full_len: total_len,
            coding_start: 1,
            coding_end: total_len,
            max_query_len: 20_000,
            min_mapped_coding_nt: 400,
            full_span_min_nt: 10_000,
            full_mapped_min_nt: 9_000,
            is_segmented: true,
            forbid_intra_segment_recomb: true,
            use_hiv1_legacy_gene_mapper: false,
            hypervariable_masks: Vec::new(),
            landmarks,
            locus_scopes: vec![],
            basal_trunks: vec![
                simple_trunk("H1N1pdm09"),
                simple_trunk("H3N2"),
                simple_trunk("H5N1"),
                simple_trunk("TRIG"),
            ],
            ghost_relabel_rule: None,
            align_founder_keys: Vec::new(),
            epi_founder_order: Vec::new(),
            epi_founder_penalty: 0.005,
            epi_certified_replacements: Vec::new(),
            pack_version: "test".to_string(),
            snapshot_date: "2026-11-16".to_string(),
            curation_sources: String::new(),
            certification_criteria: String::new(),
            default_temporal_slack_yr: None,
            chronaeon_half_life_yr: None,
        };

        let anc = deterministic_dna(total_len, 0xCAFE_BABE_1234_5678);
        let h1 = mutate_periodic(&anc, 20, 0, 1);
        let h3 = mutate_periodic(&anc, 20, 5, 2);
        let h5 = mutate_periodic(&anc, 20, 10, 3);
        let trig = mutate_periodic(&anc, 20, 15, 1);

        let profiles_fasta = format!(
            ">H1N1pdm09\n{}\n>H3N2\n{}\n>H5N1\n{}\n>TRIG\n{}\n",
            String::from_utf8_lossy(&h1),
            String::from_utf8_lossy(&h3),
            String::from_utf8_lossy(&h5),
            String::from_utf8_lossy(&trig),
        );

        let pack = PathogenProfilePack {
            spec,
            profiles_tsv: None,
            profiles_fasta: Some(profiles_fasta),
            founder_profiles_tsv: None,
            founder_profiles_fasta: None,
            certified_profiles_tsv: None,
            compendium_json: None,
        };
        let engine = EggmanEngine::from_pathogen_pack(pack).unwrap();
        assert_eq!(engine.pathogen_name(), "Influenza A");
        assert_eq!(engine.reference_name(), "FluA-8seg");
        assert_eq!(engine.basal_trunk_count(), 4);

        let cfg = ScreeningConfig::default();

        // 1. Pure H5N1 constellation
        let res_pure = engine.classify_sequence(
            "Flu_Pure_H5N1",
            std::str::from_utf8(&h5).unwrap(),
            &cfg,
        );
        assert!(res_pure.is_clonal);
        assert_eq!(res_pure.called_lineage, "H5N1");
        assert_eq!(res_pure.called_category, "Pure_Subtype");
        assert_eq!(res_pure.structural_formula, "H5N1(Clonal)");

        // 2. 6+2 reassortant: H1N1pdm09 internal genes + H5N1 HA & NA
        let mut q62 = h1.clone();
        q62[(6916 - 1)..8693].copy_from_slice(&h5[(6916 - 1)..8693]);
        q62[(10259 - 1)..11671].copy_from_slice(&h5[(10259 - 1)..11671]);
        let res_62 = engine.classify_sequence(
            "Flu_6plus2",
            std::str::from_utf8(&q62).unwrap(),
            &cfg,
        );
        assert!(!res_62.is_clonal);
        assert_eq!(res_62.called_lineage, "URF");
        assert_eq!(res_62.structural_formula, "H1N1pdm09[H5N1:HA+NA]");

        // 3. Triple reassortant: H5N1 backbone [H1N1pdm09:PB2+NP | TRIG:MP]
        let mut q3 = h5.clone();
        q3[0..2341].copy_from_slice(&h1[0..2341]); // PB2
        q3[(8694 - 1)..10258].copy_from_slice(&h1[(8694 - 1)..10258]); // NP
        q3[(11672 - 1)..12698].copy_from_slice(&trig[(11672 - 1)..12698]); // MP
        let res_3 = engine.classify_sequence(
            "Flu_Triple",
            std::str::from_utf8(&q3).unwrap(),
            &cfg,
        );
        assert!(!res_3.is_clonal);
        assert_eq!(res_3.called_lineage, "URF");
        assert_eq!(res_3.structural_formula, "H5N1[H1N1pdm09:PB2+NP | TRIG:MP]");

        // 4. Unaligned reassortant (5'/3' trimmed + 3-nt deletion inside HA)
        let mut unaligned = Vec::with_capacity(q62.len());
        unaligned.extend_from_slice(&q62[200..7500]);
        unaligned.extend_from_slice(&q62[7503..(total_len - 100)]);
        let res_un = engine.classify_sequence(
            "Flu_6plus2_Unaligned",
            std::str::from_utf8(&unaligned).unwrap(),
            &cfg,
        );
        assert_eq!(
            res_un.alignment_mode,
            "Auto-aligned (15-mer anchor + banded NW to FluA-8seg)"
        );
        assert!(!res_un.is_clonal);
        assert_eq!(res_un.structural_formula, "H1N1pdm09[H5N1:HA+NA]");
    }

    #[test]
    fn sars_cov2_wuhan_29903_directory_pack() {
        const REF_LEN: usize = 29_903;
        const COD_START: usize = 266;
        const COD_END: usize = 29_674;
        const COD_LEN: usize = COD_END - COD_START + 1;

        let spec = PathogenSpec {
            pathogen_name: "SARS-CoV-2".to_string(),
            reference_name: "Wuhan-Hu-1".to_string(),
            ref_full_len: REF_LEN,
            coding_start: COD_START,
            coding_end: COD_END,
            max_query_len: 0, // tests dynamic effective_max_query_len() > 20 kb
            min_mapped_coding_nt: 500,
            full_span_min_nt: 22_000,
            full_mapped_min_nt: 20_000,
            is_segmented: false,
            forbid_intra_segment_recomb: false,
            use_hiv1_legacy_gene_mapper: false,
            hypervariable_masks: Vec::new(),
            landmarks: vec![
                GeneLandmark { short_name: "ORF1a".into(), parent_gene: None, parent_min_span: 0, start: 266, end: 13468 },
                GeneLandmark { short_name: "ORF1b".into(), parent_gene: None, parent_min_span: 0, start: 13468, end: 21555 },
                GeneLandmark { short_name: "S".into(), parent_gene: None, parent_min_span: 0, start: 21563, end: 25384 },
                GeneLandmark { short_name: "ORF3a".into(), parent_gene: None, parent_min_span: 0, start: 25393, end: 26220 },
                GeneLandmark { short_name: "E".into(), parent_gene: None, parent_min_span: 0, start: 26245, end: 26472 },
                GeneLandmark { short_name: "M".into(), parent_gene: None, parent_min_span: 0, start: 26523, end: 27191 },
                GeneLandmark { short_name: "N".into(), parent_gene: None, parent_min_span: 0, start: 28274, end: 29533 },
            ],
            locus_scopes: vec![LocusScopeRule {
                name: "S".to_string(),
                min_start: 21500,
                min_end: 25300,
                max_end: 25500,
                min_span_nt: 3000,
            }],
            basal_trunks: vec![
                simple_trunk("BA.2"),
                simple_trunk("BJ.1"),
                simple_trunk("BA.5"),
            ],
            ghost_relabel_rule: None,
            align_founder_keys: Vec::new(),
            epi_founder_order: Vec::new(),
            epi_founder_penalty: 0.005,
            epi_certified_replacements: Vec::new(),
            pack_version: "test".to_string(),
            snapshot_date: "2026-11-16".to_string(),
            curation_sources: String::new(),
            certification_criteria: String::new(),
            default_temporal_slack_yr: None,
            chronaeon_half_life_yr: None,
        };

        let anc = deterministic_dna(REF_LEN, 0x5A25_C0F2_2022_0001);
        let ba2 = mutate_periodic(&anc, 16, 0, 1);
        let bj1 = mutate_periodic(&anc, 16, 5, 2);
        let ba5 = mutate_periodic(&anc, 16, 10, 3);

        // Certified XBB: BA.2 backbone (266..21562 and 25385..29674) + BJ.1 Spike (21563..25384)
        // plus unique clade-stem synapomorphies across coding genome.
        let mut xbb_full = ba2.clone();
        xbb_full[(21563 - 1)..25384].copy_from_slice(&bj1[(21563 - 1)..25384]);
        for k in 0..30usize {
            let pos = COD_START + 200 + k * 900;
            if pos < COD_END && pos % 16 != 0 && pos % 16 != 5 && pos % 16 != 10 {
                xbb_full[pos] = match xbb_full[pos] {
                    b'A' => b'G',
                    b'G' => b'A',
                    b'C' => b'T',
                    _ => b'C',
                };
            }
        }
        let xbb_coding = std::str::from_utf8(&xbb_full[(COD_START - 1)..COD_END]).unwrap();
        assert_eq!(xbb_coding.len(), COD_LEN);

        let profiles_fasta = format!(
            ">BA.2\n{}\n>BJ.1\n{}\n>BA.5\n{}\n",
            String::from_utf8_lossy(&ba2),
            String::from_utf8_lossy(&bj1),
            String::from_utf8_lossy(&ba5),
        );
        let certified_tsv = format!(
            "# label\tcategory\tflags\tdag_parents\tn_seqs\tseq\nXBB\tCRF\tglobal_stem\t-\t50\t{}\n",
            xbb_coding
        );
        let compendium_entries = vec![CompendiumEntry {
            label: "XBB".to_string(),
            short_code: "XBB".to_string(),
            category: "CRF".to_string(),
            status: "Certified".to_string(),
            descriptive_name: "BA.2[BJ.1:S]".to_string(),
            structural_formula: "BA.2[BJ.1:S]".to_string(),
            backbone: "BA.2".to_string(),
            t_mrca: None,
            t_mrca_ci_low: None,
            t_mrca_ci_high: None,
            first_sampled_year: None,
            chronaeon_half_life_yr: None,
            hxb2_barcode: "BA.2:266-21562|BJ.1:21563-25384|BA.2:25385-29674".to_string(),
            scope: "full".to_string(),
        }];

        let dir = temp_pack_dir("sarscov2");
        std::fs::write(
            dir.join("manifest.json"),
            serde_json::to_string_pretty(&spec).unwrap(),
        )
        .unwrap();
        std::fs::write(dir.join("profiles.fasta"), profiles_fasta).unwrap();
        std::fs::write(dir.join("certified_profiles.tsv"), certified_tsv).unwrap();
        std::fs::write(
            dir.join("compendium.json"),
            serde_json::to_string_pretty(&compendium_entries).unwrap(),
        )
        .unwrap();

        let engine = EggmanEngine::from_pack_path(&dir).unwrap();
        let _ = std::fs::remove_dir_all(&dir);

        assert_eq!(engine.pathogen_name(), "SARS-CoV-2");
        assert_eq!(engine.reference_name(), "Wuhan-Hu-1");
        assert_eq!(engine.basal_trunk_count(), 3);
        assert_eq!(engine.certified_lineage_count(), 1);
        assert!(engine.spec().effective_max_query_len() > REF_LEN);

        let cfg = ScreeningConfig::default();

        // a + b. Unaligned >20 kb XBB query (~29.8 kb trimmed + 6-nt Spike deletion)
        let mut xbb_unaligned = Vec::with_capacity(REF_LEN);
        xbb_unaligned.extend_from_slice(&xbb_full[50..22000]);
        xbb_unaligned.extend_from_slice(&xbb_full[22006..(REF_LEN - 50)]);
        assert!(xbb_unaligned.len() > 20_000);
        let res_xbb = engine.classify_sequence(
            "SC2_XBB_Unaligned",
            std::str::from_utf8(&xbb_unaligned).unwrap(),
            &cfg,
        );
        assert_eq!(
            res_xbb.alignment_mode,
            "Auto-aligned (15-mer anchor + banded NW to Wuhan-Hu-1)"
        );
        assert_eq!(res_xbb.locus_scope, "full");
        assert_eq!(res_xbb.called_lineage, "XBB");
        assert_eq!(res_xbb.called_category, "CRF");
        assert_eq!(res_xbb.structural_formula, "BA.2[BJ.1:S]");

        // c. De novo BA.5 (1..21562) + BA.2 (21563..29903) -> URF BA.5[BA.2:S..N]
        let mut urf_full = ba5.clone();
        urf_full[(21563 - 1)..].copy_from_slice(&ba2[(21563 - 1)..]);
        let res_urf = engine.classify_sequence(
            "SC2_BA5_BA2_URF",
            std::str::from_utf8(&urf_full).unwrap(),
            &cfg,
        );
        assert!(!res_urf.is_clonal);
        assert_eq!(res_urf.called_lineage, "URF");
        assert_eq!(res_urf.structural_formula, "BA.5[BA.2:S..N]");

        // d. BA.2 with only 2 isolated BJ.1 SNPs in Spike -> pruned to BA.2(Clonal)
        let mut ba2_two_snps = ba2.clone();
        ba2_two_snps[22005] = bj1[22005];
        ba2_two_snps[22021] = bj1[22021];
        let res_pruned = engine.classify_sequence(
            "SC2_BA2_2SNPs",
            std::str::from_utf8(&ba2_two_snps).unwrap(),
            &cfg,
        );
        assert!(res_pruned.is_clonal);
        assert_eq!(res_pruned.called_lineage, "BA.2");
        assert_eq!(res_pruned.structural_formula, "BA.2(Clonal)");
    }

    #[test]
    fn embedded_builtin_packs_load_and_stitch_flu_segments() {
        let hiv = EggmanEngine::from_builtin_id("hiv1").expect("embedded hiv1");
        assert_eq!(hiv.basal_trunk_count(), 14);

        let iav = EggmanEngine::from_builtin_id("iav").expect("embedded iav");
        assert!(iav.spec().is_segmented);
        assert_eq!(iav.spec().landmarks.len(), 8);

        let sc2 = EggmanEngine::from_builtin_id("sc2").expect("embedded sc2");
        assert_eq!(sc2.reference_name(), "Wuhan-Hu-1");

        let sarbeco = EggmanEngine::from_builtin_id("sarbeco").expect("embedded sarbeco");
        assert_eq!(sarbeco.basal_trunk_count(), 5);

        // Extract circulating B3.13 profile from IAV_CERTIFIED_TSV and split into 8 out-of-order FASTA
        // records with 5'/3' UTR overhangs; also verify raw ChronAeon Tier-1 splice recovers exact
        // Tier-1 coordinates while Guard-4 blocks certifying a stemless 2021-EA1 chimera.
        let profs = parse_aligned_fasta_profiles(IAV_PROFILES_FASTA, iav.spec().ref_full_len)
            .expect("iav profiles");
        let ea1 = &profs["EA1"];
        let nam = &profs["NAm_B313"];
        let b313_line = IAV_CERTIFIED_TSV
            .lines()
            .find(|l| l.starts_with("B3.13\t"))
            .expect("B3.13 certified row");
        let b313_str = b313_line.split('\t').nth(5).expect("B3.13 seq column");
        let mut stemless_splice = nam.clone();
        for &lm_idx in &[2usize, 3, 5, 6] {
            let lm = &iav.spec().landmarks[lm_idx];
            for p in (lm.start - 1)..lm.end {
                stemless_splice[p] = ea1[p];
            }
        }
        let to_dna = |slice: &[i8]| -> String {
            slice
                .iter()
                .map(|&b| match b {
                    0 => 'A',
                    1 => 'C',
                    2 => 'G',
                    3 => 'T',
                    _ => 'N',
                })
                .collect()
        };
        let res_stemless = iav.classify_sequence(
            "Stemless_Tier1_Splice",
            &to_dna(&stemless_splice),
            &ScreeningConfig::default(),
        );
        assert_eq!(res_stemless.called_lineage, "URF");
        assert_eq!(res_stemless.structural_formula, "NAm_B313[EA1:PA+HA+NA+MP]");

        let order = [7usize, 3, 0, 5, 2, 6, 1, 4];
        let mut multi_fasta = String::new();
        for &idx in &order {
            let lm = &iav.spec().landmarks[idx];
            let seg_seq = &b313_str[(lm.start - 1)..lm.end];
            // Add realistic 5' UTR & 3' poly-U tail on HA & NA
            let padded = if idx == 3 || idx == 5 {
                format!("AGCAAAAGCAGGGG{}{}", seg_seq, "AACCCCTTTTT")
            } else {
                seg_seq.to_string()
            };
            multi_fasta.push_str(&format!(
                ">PP7551{:02}.1 Influenza A virus (A/dairy cow/Texas/24-008749-001/2024(H5N1)) segment {} ({}) complete cds\n{}\n",
                idx + 20,
                idx + 1,
                lm.short_name,
                padded
            ));
        }

        let (asm, res, _stitched) =
            iav.classify_fasta_with_segmentation(&multi_fasta, "auto", &ScreeningConfig::default());
        assert!(asm.was_multi_record_assembled);
        assert_eq!(asm.raw_record_count, 8);
        assert_eq!(asm.isolate_count, 1);
        assert_eq!(res.len(), 1);
        assert_eq!(res[0].called_lineage, "B3.13");
        assert_eq!(res[0].structural_formula, "NAm_B313[EA1:PA+HA+NA+MP]");

        // Test bare accession auto-merge (<= 8 distinct segments) vs duplicate HA collision guard
        let mut bare_8 = String::new();
        for idx in 0..8 {
            let lm = &iav.spec().landmarks[idx];
            bare_8.push_str(&format!(
                ">ACC_00{}\n{}\n",
                idx + 1,
                &b313_str[(lm.start - 1)..lm.end]
            ));
        }
        let (asm_bare, res_bare, _) =
            iav.classify_fasta_with_segmentation(&bare_8, "auto", &ScreeningConfig::default());
        assert_eq!(asm_bare.isolate_count, 1);
        assert_eq!(res_bare[0].called_lineage, "B3.13");

        let ha_lm = &iav.spec().landmarks[3];
        let dup_ha = format!(
            ">ISO_ONE\n{}\n>ISO_TWO\n{}\n",
            &b313_str[(ha_lm.start - 1)..ha_lm.end],
            to_dna(&ea1[(ha_lm.start - 1)..ha_lm.end])
        );
        let (asm_ha, _, _) =
            iav.classify_fasta_with_segmentation(&dup_ha, "auto", &ScreeningConfig::default());
        assert_eq!(asm_ha.isolate_count, 2);
        assert!(!asm_ha.was_multi_record_assembled);
    }

    #[test]
    fn trunk_spec_temporal_fields_backward_compat() {
        let legacy_trunk_json = r#"{"name": "B.1.1.7"}"#;
        let t: TrunkSpec = serde_json::from_str(legacy_trunk_json).unwrap();
        assert_eq!(t.name, "B.1.1.7");
        assert!(t.is_pure_basal);
        assert!(t.t_mrca.is_none());
        assert!(t.t_mrca_ci.is_none());
        assert!(t.first_sampled_year.is_none());
        assert!(t.chronaeon_half_life_yr.is_none());
        assert!(t.emergence_floor().is_none());

        let timed_json = r#"{
            "name": "LF.7",
            "t_mrca": 2024.55,
            "t_mrca_ci": [2024.42, 2024.68],
            "first_sampled_year": 2024.49,
            "chronaeon_half_life_yr": 0.25
        }"#;
        let t2: TrunkSpec = serde_json::from_str(timed_json).unwrap();
        assert_eq!(t2.t_mrca, Some(2024.55));
        assert_eq!(t2.t_mrca_ci, Some((2024.42, 2024.68)));
        assert_eq!(t2.first_sampled_year, Some(2024.49));
        assert_eq!(t2.chronaeon_half_life_yr, Some(0.25));
        assert!((t2.emergence_floor().unwrap() - 2024.42).abs() < 1e-9);

        let rt = serde_json::to_string(&t2).unwrap();
        let t3: TrunkSpec = serde_json::from_str(&rt).unwrap();
        assert_eq!(t3.emergence_floor(), t2.emergence_floor());

        assert_eq!(PathogenSpec::hiv1_group_m_default().effective_temporal_slack_yr(), 3.0);
        assert!(
            PathogenSpec::hiv1_group_m_default()
                .basal_trunks
                .iter()
                .all(|tr| tr.emergence_floor().is_some())
        );
    }

    #[test]
    fn header_decimal_year_parsing() {
        let y1 = parse_header_decimal_year("hCoV-19/USA/WA-1|2024-03-15|EPI_ISL_1").unwrap();
        assert!((y1 - 2024.202186).abs() < 1e-4, "got {}", y1);

        let y2 = parse_header_decimal_year("isolate/2022-11-01/spike").unwrap();
        assert!((y2 - 2022.832877).abs() < 1e-4, "got {}", y2);

        let y3 = parse_header_decimal_year("sc2|2021.45|BA.1").unwrap();
        assert!((y3 - 2021.45).abs() < 1e-9, "got {}", y3);

        let y4 = parse_header_decimal_year("sample_2023-07").unwrap();
        assert!((y4 - 2023.534247).abs() < 1e-4, "got {}", y4);

        let y5 = parse_header_decimal_year("clade_B|2023|ref").unwrap();
        assert!((2023.0..2024.0).contains(&y5), "got {}", y5);

        let y6 = parse_header_decimal_year("A/dairy cow/Texas/24-008749-001/2024(H5N1)").unwrap();
        assert!((2024.0..2025.0).contains(&y6), "got {}", y6);

        let y7 = parse_header_decimal_year("sample|2024-XX-XX|raw").unwrap();
        assert!((2024.0..2025.0).contains(&y7), "got {}", y7);

        // False-positive guards: taxids, accessions, strain serial numbers
        assert!(parse_header_decimal_year("2697049.2956128|PP755589.1|SRR20241234").is_none());
        assert!(parse_header_decimal_year("TX/2024-008749_sample").is_none());
    }

    #[test]
    fn temporal_guard_masks_future_sc2_donors_and_clades() {
        const REF_LEN: usize = 29_903;
        const COD_START: usize = 266;
        const COD_END: usize = 29_674;
        const COD_LEN: usize = COD_END - COD_START + 1;

        let timed_trunk = |name: &str, t_mrca: f64, lo: f64, first: f64| TrunkSpec {
            name: name.to_string(),
            is_pure_basal: true,
            penalty: 0.0,
            active_range: None,
            fallback_profile: None,
            display_alias: None,
            clonal_certified_label: None,
            t_mrca: Some(t_mrca),
            t_mrca_ci: Some((lo, t_mrca + 0.08)),
            first_sampled_year: Some(first),
            chronaeon_half_life_yr: Some(0.25),
        };

        let spec = PathogenSpec {
            pathogen_name: "SARS-CoV-2".to_string(),
            reference_name: "MN908947".to_string(),
            ref_full_len: REF_LEN,
            coding_start: COD_START,
            coding_end: COD_END,
            max_query_len: 35_000,
            min_mapped_coding_nt: 3_000,
            full_span_min_nt: 24_000,
            full_mapped_min_nt: 20_000,
            is_segmented: false,
            forbid_intra_segment_recomb: false,
            use_hiv1_legacy_gene_mapper: false,
            hypervariable_masks: Vec::new(),
            landmarks: vec![
                GeneLandmark { short_name: "ORF1ab".into(), parent_gene: None, parent_min_span: 0, start: 266, end: 21555 },
                GeneLandmark { short_name: "S".into(), parent_gene: None, parent_min_span: 0, start: 21563, end: 25384 },
                GeneLandmark { short_name: "3'-end".into(), parent_gene: None, parent_min_span: 0, start: 25385, end: 29674 },
            ],
            locus_scopes: vec![LocusScopeRule {
                name: "S".into(),
                min_start: 21_500,
                min_end: 0,
                max_end: 25_600,
                min_span_nt: 2_500,
            }],
            basal_trunks: vec![
                timed_trunk("BA.1", 2021.75, 2021.70, 2021.78),
                timed_trunk("BA.2", 2021.85, 2021.80, 2021.88),
                timed_trunk("BA.5", 2022.15, 2022.10, 2022.18),
                timed_trunk("KS.1.1", 2024.30, 2024.25, 2024.32),
                timed_trunk("LF.7", 2024.50, 2024.45, 2024.52),
                timed_trunk("LP.8.1", 2024.85, 2024.80, 2024.88),
            ],
            ghost_relabel_rule: None,
            align_founder_keys: Vec::new(),
            epi_founder_order: Vec::new(),
            epi_founder_penalty: 0.005,
            epi_certified_replacements: Vec::new(),
            pack_version: "test-temporal".to_string(),
            snapshot_date: "2026-11-17".to_string(),
            curation_sources: String::new(),
            certification_criteria: String::new(),
            default_temporal_slack_yr: Some(0.5),
            chronaeon_half_life_yr: Some(0.25),
        };

        let anc = deterministic_dna(REF_LEN, 0x2021_2025_5C02_0001);
        let ba1 = mutate_periodic(&anc, 24, 0, 1);
        let ba2 = mutate_periodic(&anc, 24, 4, 2);
        let ba5 = mutate_periodic(&anc, 24, 8, 3);
        let ks11 = mutate_periodic(&anc, 24, 12, 1);
        let lf7 = mutate_periodic(&anc, 24, 16, 2);
        let lp81 = mutate_periodic(&anc, 24, 20, 3);

        // Future Tier-2 recombinant clade XEC (2024.55) derived from KS.1.1 + LF.7,
        // PLUS a BA.1+BA.2 synthetic mosaic X_FUTURE whose parents are 2021 but clade itself emerged in 2024.6.
        let mut xec_full = ks11.clone();
        xec_full[(21563 - 1)..25384].copy_from_slice(&lf7[(21563 - 1)..25384]);
        for k in 0..25usize {
            let pos = COD_START + 250 + k * 900;
            if pos < COD_END && pos % 24 >= 22 {
                xec_full[pos] = if xec_full[pos] == b'A' { b'G' } else { b'A' };
            }
        }
        let mut xfut_full = ba1.clone();
        xfut_full[(21563 - 1)..25384].copy_from_slice(&ba2[(21563 - 1)..25384]);
        for k in 0..25usize {
            let pos = COD_START + 300 + k * 900;
            if pos < COD_END && pos % 24 >= 22 {
                xfut_full[pos] = if xfut_full[pos] == b'C' { b'T' } else { b'C' };
            }
        }

        let xec_cod = std::str::from_utf8(&xec_full[(COD_START - 1)..COD_END]).unwrap();
        let xfut_cod = std::str::from_utf8(&xfut_full[(COD_START - 1)..COD_END]).unwrap();
        assert_eq!(xec_cod.len(), COD_LEN);

        let profiles_fasta = format!(
            ">BA.1\n{}\n>BA.2\n{}\n>BA.5\n{}\n>KS.1.1\n{}\n>LF.7\n{}\n>LP.8.1\n{}\n",
            String::from_utf8_lossy(&ba1),
            String::from_utf8_lossy(&ba2),
            String::from_utf8_lossy(&ba5),
            String::from_utf8_lossy(&ks11),
            String::from_utf8_lossy(&lf7),
            String::from_utf8_lossy(&lp81),
        );
        let certified_tsv = format!(
            "# label\tcategory\tflags\tdag_parents\tn_seqs\tseq\nXEC\tCRF\tglobal_stem\t-\t50\t{}\nX_FUTURE\tCRF\tglobal_stem\t-\t50\t{}\n",
            xec_cod, xfut_cod
        );
        let compendium_entries = vec![
            CompendiumEntry {
                label: "XEC".to_string(),
                short_code: "XEC".to_string(),
                category: "CRF".to_string(),
                status: "Certified".to_string(),
                descriptive_name: "KS.1.1[LF.7:S]".to_string(),
                structural_formula: "KS.1.1[LF.7:S]".to_string(),
                backbone: "KS.1.1".to_string(),
                t_mrca: Some(2024.55),
                t_mrca_ci_low: Some(2024.50),
                t_mrca_ci_high: Some(2024.62),
                first_sampled_year: Some(2024.58),
                chronaeon_half_life_yr: Some(0.25),
                hxb2_barcode: "KS.1.1:266-21562|LF.7:21563-25384|KS.1.1:25385-29674".to_string(),
                scope: "full".to_string(),
            },
            CompendiumEntry {
                label: "X_FUTURE".to_string(),
                short_code: "X_FUTURE".to_string(),
                category: "CRF".to_string(),
                status: "Certified".to_string(),
                descriptive_name: "BA.1[BA.2:S]".to_string(),
                structural_formula: "BA.1[BA.2:S]".to_string(),
                backbone: "BA.1".to_string(),
                t_mrca: Some(2024.60),
                t_mrca_ci_low: Some(2024.55),
                t_mrca_ci_high: Some(2024.70),
                first_sampled_year: Some(2024.62),
                chronaeon_half_life_yr: Some(0.25),
                hxb2_barcode: "BA.1:266-21562|BA.2:21563-25384|BA.1:25385-29674".to_string(),
                scope: "full".to_string(),
            },
        ];

        let dir = temp_pack_dir("sc2_temporal");
        std::fs::write(dir.join("manifest.json"), serde_json::to_string_pretty(&spec).unwrap()).unwrap();
        std::fs::write(dir.join("profiles.fasta"), profiles_fasta).unwrap();
        std::fs::write(dir.join("certified_profiles.tsv"), certified_tsv).unwrap();
        std::fs::write(dir.join("compendium.json"), serde_json::to_string_pretty(&compendium_entries).unwrap()).unwrap();

        let engine = EggmanEngine::from_pack_path(&dir).unwrap();
        let cfg = ScreeningConfig::default();

        // 1. Default (unguarded) call recognizes LF.7 / LP.8.1 / XEC even with a 2021 header.
        let lf7_seq = String::from_utf8_lossy(&lf7);
        let lp81_seq = String::from_utf8_lossy(&lp81);
        let xec_seq = String::from_utf8_lossy(&xec_full);
        let xfut_seq = String::from_utf8_lossy(&xfut_full);

        let r_un = engine.classify_sequence("anachro|2021-12-01", &lf7_seq, &cfg);
        assert_eq!(r_un.called_lineage, "LF.7");
        assert!(r_un.masked_future_trunks.is_empty());

        // 2. With --temporal-guard on a 2021/2022 header, 2024-2025 donors are masked.
        let tg = TemporalGuardConfig {
            enabled: true,
            query_year: None,
            slack_yr: None,
        };
        let r_2021 = engine.classify_sequence_with_temporal("anachro|2021-12-01", &lf7_seq, &cfg, &tg);
        assert!(r_2021.query_year.unwrap() > 2021.9 && r_2021.query_year.unwrap() < 2022.0);
        for fut in ["KS.1.1", "LF.7", "LP.8.1"] {
            assert!(r_2021.masked_future_trunks.iter().any(|m| m == fut), "{:?}", r_2021.masked_future_trunks);
            assert_ne!(r_2021.called_lineage, fut);
            assert!(!r_2021.constituent_trunks.iter().any(|c| c == fut));
        }

        let r_2022_lp = engine.classify_sequence_with_temporal("early|2022-06-15", &lp81_seq, &cfg, &tg);
        assert_ne!(r_2022_lp.called_lineage, "LP.8.1");
        assert!(r_2022_lp.masked_future_trunks.contains(&"LP.8.1".to_string()));
        assert!(r_2022_lp.masked_future_trunks.contains(&"LF.7".to_string()));
        assert!(r_2022_lp.masked_future_trunks.contains(&"KS.1.1".to_string()));

        // Future Tier-2 XEC cannot be certified on 2022 query, nor matched in compendium.
        let r_2022_xec = engine.classify_sequence_with_temporal("early_xec|2022-06-15", &xec_seq, &cfg, &tg);
        assert_ne!(r_2022_xec.called_lineage, "XEC");
        assert!(!r_2022_xec.matches.iter().any(|m| m.label == "XEC"));

        // Future Tier-2 clade X_FUTURE whose parents are 2021 (BA.1 + BA.2) is blocked on 2022 query!
        let r_2022_xfut = engine.classify_sequence_with_temporal("early_xfut|2022-06-15", &xfut_seq, &cfg, &tg);
        assert_eq!(r_2022_xfut.called_lineage, "URF");
        assert!(!r_2022_xfut.matches.iter().any(|m| m.label == "X_FUTURE"));

        // 3. On a 2025 query, all trunks and certified clades remain active.
        let r_2025_lf = engine.classify_sequence_with_temporal("modern|2025-03-01", &lf7_seq, &cfg, &tg);
        assert!(r_2025_lf.masked_future_trunks.is_empty());
        assert_eq!(r_2025_lf.called_lineage, "LF.7");

        let r_2025_lp = engine.classify_sequence_with_temporal("modern|2025-03-01", &lp81_seq, &cfg, &tg);
        assert!(r_2025_lp.masked_future_trunks.is_empty());
        assert_eq!(r_2025_lp.called_lineage, "LP.8.1");

        let r_2025_xfut = engine.classify_sequence_with_temporal("modern_xfut|2025-03-01", &xfut_seq, &cfg, &tg);
        assert_eq!(r_2025_xfut.called_lineage, "X_FUTURE");

        // 4. Safe fallback if temporal threshold would leave < 2 trunks active (query_year = 1990).
        let tg_1990 = TemporalGuardConfig {
            enabled: true,
            query_year: Some(1990.0),
            slack_yr: Some(0.5),
        };
        let r_1990 = engine.classify_sequence_with_temporal("ancient", &lf7_seq, &cfg, &tg_1990);
        assert_eq!(r_1990.masked_future_trunks.len(), 4);
        assert_eq!(r_1990.query_year, Some(1990.0));
        assert!(!r_1990.segments.is_empty());

        let _ = std::fs::remove_dir_all(&dir);
    }
}

/// Parse an explicit or delimited decimal year from a FASTA header.
/// Supports `YYYY-MM-DD`, `YYYY/MM/DD`, `YYYY.frac`, `YYYY-MM`, `YYYY/MM`, `YYYY-XX-XX`, and delimited `YYYY`.
pub fn parse_header_decimal_year(header: &str) -> Option<f64> {
    let b = header.as_bytes();
    let n = b.len();
    let parse_u32 = |s: &[u8]| -> Option<u32> {
        if !s.iter().all(|c| c.is_ascii_digit()) {
            return None;
        }
        std::str::from_utf8(s).ok()?.parse::<u32>().ok()
    };
    let is_valid_yr = |y: u32| (1900..=2100).contains(&y);
    let left_ok = |i: usize| -> bool {
        i == 0 || (!b[i - 1].is_ascii_alphanumeric() && b[i - 1] != b'.')
    };

    // Pass 1: YYYY-MM-DD or YYYY/MM/DD
    if n >= 10 {
        for i in 0..=(n - 10) {
            if !left_ok(i) {
                continue;
            }
            let sep = b[i + 4];
            if (sep == b'-' || sep == b'/') && b[i + 7] == sep {
                if let (Some(y), Some(m), Some(d)) = (
                    parse_u32(&b[i..i + 4]),
                    parse_u32(&b[i + 5..i + 7]),
                    parse_u32(&b[i + 8..i + 10]),
                ) {
                    if let Some(dec) = ymd_to_decimal_year(y, m, d) {
                        let right_ok = i + 10 == n
                            || !b[i + 10].is_ascii_digit()
                                && b[i + 10] != b'.'
                                && !(b[i + 10] == b'-'
                                    && b.get(i + 11).is_some_and(|c| c.is_ascii_digit()));
                        if right_ok {
                            return Some(dec);
                        }
                    }
                }
            }
        }
    }

    // Pass 2: decimal year YYYY.fff (e.g. |2021.45|)
    if n >= 6 {
        for i in 0..=(n - 6) {
            if !left_ok(i) || b[i + 4] != b'.' || !b[i + 5].is_ascii_digit() {
                continue;
            }
            let Some(y) = parse_u32(&b[i..i + 4]) else {
                continue;
            };
            if !is_valid_yr(y) {
                continue;
            }
            let mut end = i + 5;
            while end < n && b[end].is_ascii_digit() && end - (i + 5) < 6 {
                end += 1;
            }
            let right_ok = end == n || (!b[end].is_ascii_alphanumeric() && b[end] != b'.');
            if right_ok {
                if let Ok(val) = std::str::from_utf8(&b[i..end]).unwrap_or("").parse::<f64>() {
                    if val.is_finite() && (1900.0..2101.0).contains(&val) {
                        return Some(val);
                    }
                }
            }
        }
    }

    // Pass 3: YYYY-MM or YYYY/MM (plus YYYY-XX-XX / YYYY-XX placeholders -> mid-year)
    if n >= 7 {
        for i in 0..=(n - 7) {
            if !left_ok(i) {
                continue;
            }
            let sep = b[i + 4];
            if sep != b'-' && sep != b'/' {
                continue;
            }
            let Some(y) = parse_u32(&b[i..i + 4]) else {
                continue;
            };
            if !is_valid_yr(y) {
                continue;
            }
            if (&b[i + 5..i + 7] == b"XX" || &b[i + 5..i + 7] == b"xx")
                && (i + 7 == n
                    || !b[i + 7].is_ascii_alphanumeric()
                    || (i + 10 <= n
                        && b[i + 7] == sep
                        && (&b[i + 8..i + 10] == b"XX" || &b[i + 8..i + 10] == b"xx")
                        && (i + 10 == n || !b[i + 10].is_ascii_alphanumeric())))
            {
                return Some(y as f64 + 0.5);
            }
            if let Some(m) = parse_u32(&b[i + 5..i + 7]) {
                let right_ok = i + 7 == n
                    || (!b[i + 7].is_ascii_alphanumeric()
                        && b[i + 7] != b'-'
                        && b[i + 7] != b'/'
                        && b[i + 7] != b'.');
                if right_ok {
                    if let Some(dec) = ymd_to_decimal_year(y, m, 15) {
                        return Some(dec);
                    }
                }
            }
        }
    }

    // Pass 4: delimited bare YYYY (e.g. |2023| or /2024(H5N1)) -> mid-year YYYY.5
    if n >= 4 {
        for i in 0..=(n - 4) {
            if !left_ok(i) {
                continue;
            }
            if let Some(y) = parse_u32(&b[i..i + 4]) {
                if !is_valid_yr(y) {
                    continue;
                }
                let right_ok = i + 4 == n
                    || (!b[i + 4].is_ascii_alphanumeric()
                        && b[i + 4] != b'.'
                        && !(b[i + 4] == b'-'
                            && b.get(i + 5).is_some_and(|c| c.is_ascii_alphanumeric())));
                if right_ok {
                    return Some(y as f64 + 0.5);
                }
            }
        }
    }

    None
}

fn is_leap_year(y: u32) -> bool {
    (y % 4 == 0 && y % 100 != 0) || (y % 400 == 0)
}

fn days_in_month(y: u32, m: u32) -> u32 {
    match m {
        1 | 3 | 5 | 7 | 8 | 10 | 12 => 31,
        4 | 6 | 9 | 11 => 30,
        2 => {
            if is_leap_year(y) {
                29
            } else {
                28
            }
        }
        _ => 0,
    }
}

fn ymd_to_decimal_year(y: u32, m: u32, d: u32) -> Option<f64> {
    if !(1900..=2100).contains(&y) || !(1..=12).contains(&m) {
        return None;
    }
    let dim = days_in_month(y, m);
    if d < 1 || d > dim {
        return None;
    }
    let mut doy0 = d - 1;
    for mm in 1..m {
        doy0 += days_in_month(y, mm);
    }
    let total_days = if is_leap_year(y) { 366.0 } else { 365.0 };
    Some(y as f64 + (doy0 as f64) / total_days)
}

/// Extract canonical Influenza isolate strain key or prefix before segment delimiter.
fn extract_flu_isolate_key(header: &str) -> Option<String> {
    let h = header.trim();
    if h.is_empty() {
        return None;
    }
    // 1. Check for canonical WHO Influenza A/B strain string: A/host/location/id/year(HxNy)
    for (idx, _) in h.match_indices('/') {
        if idx >= 1 {
            let prefix_char = h.as_bytes()[idx - 1];
            if prefix_char == b'A' || prefix_char == b'B' {
                let start = idx - 1;
                let prev_ok = start == 0
                    || matches!(
                        h.as_bytes()[start - 1],
                        b'(' | b'|' | b' ' | b'\t' | b'_' | b',' | b'['
                    );
                if prev_ok {
                    let sub = &h[start..];
                    // Require at least 3 slashes (A/loc/id/year or A/host/loc/id/year)
                    let enclosed_in_paren = start > 0 && h.as_bytes()[start - 1] == b'(';
                    let mut end_off = sub.len();
                    if enclosed_in_paren {
                        let mut depth = 1i32;
                        for (j, ch) in sub.char_indices() {
                            match ch {
                                '(' => depth += 1,
                                ')' => {
                                    depth -= 1;
                                    if depth == 0 {
                                        end_off = j;
                                        break;
                                    }
                                }
                                _ => {}
                            }
                        }
                    } else {
                        for (j, ch) in sub.char_indices() {
                            if ch == '|' || ch == '[' || ch == ']' {
                                end_off = j;
                                break;
                            }
                        }
                    }
                    let mut cand = sub[..end_off].trim();
                    if let Some(pos) = cand.to_ascii_lowercase().find(" segment ") {
                        cand = cand[..pos].trim();
                    }
                    if cand.matches('/').count() >= 3 && cand.len() >= 8 {
                        return Some(cand.to_string());
                    }
                }
            }
        }
    }

    // 2. Check pipe/underscore/slash delimited first token for trailing segment suffix (e.g., Sample1|PB2)
    let first_tok = h.split_whitespace().next().unwrap_or(h);
    let seg_tokens = [
        "PB2", "PB1", "PA", "HA", "NP", "NA", "MP", "M1", "M2", "NS", "NS1", "NEP", "SEG1",
        "SEG2", "SEG3", "SEG4", "SEG5", "SEG6", "SEG7", "SEG8",
    ];
    for delim in ['|', '_', '/', '-'] {
        if let Some(pos) = first_tok.rfind(delim) {
            let suffix = first_tok[pos + 1..].to_ascii_uppercase();
            if seg_tokens.contains(&suffix.as_str()) && pos >= 2 {
                return Some(first_tok[..pos].to_string());
            }
        }
    }
    None
}

/// Infer Flu A segment landmark index (0..7) strictly from delimited header tokens or `segment N`.
fn infer_flu_segment_from_header(header: &str, landmarks: &[GeneLandmark]) -> Option<usize> {
    let lower = header.to_ascii_lowercase();
    let find_lm = |name: &str| landmarks.iter().position(|lm| lm.short_name.eq_ignore_ascii_case(name));
    for (seg_num, name) in [
        ("segment 1", "PB2"),
        ("segment 2", "PB1"),
        ("segment 3", "PA"),
        ("segment 4", "HA"),
        ("segment 5", "NP"),
        ("segment 6", "NA"),
        ("segment 7", "MP"),
        ("segment 8", "NS"),
        ("seg 1", "PB2"),
        ("seg 2", "PB1"),
        ("seg 3", "PA"),
        ("seg 4", "HA"),
        ("seg 5", "NP"),
        ("seg 6", "NA"),
        ("seg 7", "MP"),
        ("seg 8", "NS"),
        ("seg1", "PB2"),
        ("seg2", "PB1"),
        ("seg3", "PA"),
        ("seg4", "HA"),
        ("seg5", "NP"),
        ("seg6", "NA"),
        ("seg7", "MP"),
        ("seg8", "NS"),
    ] {
        if lower.contains(seg_num) {
            return find_lm(name);
        }
    }
    let tokens: Vec<String> = lower
        .split(|c: char| !c.is_ascii_alphanumeric())
        .filter(|t| !t.is_empty())
        .map(|s| s.to_string())
        .collect();
    for tok in tokens.iter().rev() {
        let mapped = match tok.as_str() {
            "pb2" => Some("PB2"),
            "pb1" => Some("PB1"),
            "pa" => Some("PA"),
            "ha" | "hemagglutinin" | "haemagglutinin" => Some("HA"),
            "np" | "nucleoprotein" => Some("NP"),
            "na" | "neuraminidase" => Some("NA"),
            "mp" | "m1" | "m2" | "matrix" => Some("MP"),
            "ns" | "ns1" | "nep" | "nonstructural" => Some("NS"),
            _ => None,
        };
        if let Some(m) = mapped {
            if let Some(idx) = find_lm(m) {
                return Some(idx);
            }
        }
    }
    None
}

/// One-sided Fisher's exact test (`alternative = "greater"`) on a 2x2 contingency table `[[a, b], [c, d]]`.
pub fn fisher_exact_greater(a: usize, b: usize, c: usize, d: usize) -> f64 {
    let n = a + b + c + d;
    if n == 0 {
        return 1.0;
    }
    let r1 = a + b;
    let c1 = a + c;
    let max_a = r1.min(c1);

    let log_fact = |k: usize| -> f64 {
        if k <= 1 {
            0.0
        } else {
            // Stirling / Lanczos lgamma(k + 1)
            ln_gamma((k + 1) as f64)
        }
    };

    let log_denom = log_fact(n) - log_fact(r1) - log_fact(n - r1) - log_fact(c1) - log_fact(n - c1);
    let mut p_sum = 0.0f64;

    for x in a..=max_a {
        let y = r1 - x;
        let z = c1 - x;
        if n + x < r1 + c1 {
            continue;
        }
        let w = n + x - r1 - c1;
        let log_p = -log_denom - log_fact(x) - log_fact(y) - log_fact(z) - log_fact(w);
        p_sum += log_p.exp();
    }
    p_sum.clamp(0.0, 1.0)
}

/// One-sided exact binomial test upper tail `P(Bin(n, 0.5) >= k)`.
fn binom_upper_tail_half(k: usize, n: usize) -> f64 {
    if k == 0 {
        return 1.0;
    }
    if k > n || n == 0 {
        return 0.0;
    }
    let log_n_fac = ln_gamma((n + 1) as f64);
    let log_half_pow_n = -(n as f64) * std::f64::consts::LN_2;
    let mut sum = 0.0f64;
    for i in k..=n {
        let log_p = log_n_fac
            - ln_gamma((i + 1) as f64)
            - ln_gamma((n - i + 1) as f64)
            + log_half_pow_n;
        sum += log_p.exp();
    }
    sum.clamp(0.0, 1.0)
}

/// Lanczos approximation for `ln(Gamma(z))`.
fn ln_gamma(z: f64) -> f64 {
    const COEFFS: [f64; 9] = [
        0.99999999999980993,
        676.5203681218851,
        -1259.1392167224028,
        771.32342877765313,
        -176.61502916214059,
        12.507343278686905,
        -0.13857109526572012,
        9.9843695780195716e-6,
        1.5056327351493116e-7,
    ];
    let x = z - 1.0;
    let mut sum = COEFFS[0];
    for (i, &c) in COEFFS.iter().enumerate().skip(1) {
        sum += c / (x + i as f64);
    }
    let t = x + 7.5;
    0.5 * (2.0 * std::f64::consts::PI).ln() + (x + 0.5) * t.ln() - t + sum.ln()
}
