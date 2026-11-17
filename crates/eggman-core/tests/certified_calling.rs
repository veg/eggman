use eggman_core::*;

const PROFILES: &str = include_str!(concat!(env!("CARGO_MANIFEST_DIR"), "/data/profiles_9719.tsv"));
const FOUNDERS: &str = include_str!(concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/data/founder_profiles_9719.tsv"
));

fn load_tsv_seq(tsv: &str, name: &str) -> String {
    tsv.lines()
        .find_map(|l| {
            let (n, s) = l.split_once('\t')?;
            (n == name).then(|| s.trim().to_string())
        })
        .unwrap_or_else(|| panic!("Missing profile {} in TSV", name))
}

fn example(name: &str) -> String {
    std::fs::read_to_string(format!(
        "{}/../../docs/examples/{}.fasta",
        env!("CARGO_MANIFEST_DIR"),
        name
    ))
    .unwrap()
}

fn strip_gaps(seq: &str) -> String {
    seq.chars()
        .filter(|&c| matches!(c, 'A' | 'C' | 'G' | 'T' | 'U' | 'a' | 'c' | 'g' | 't' | 'u'))
        .collect()
}

#[test]
fn demos_unaligned_assign_certified_crfs_and_preserve_basal() {
    let engine = EggmanEngine::new();
    let cfg = ScreeningConfig::default();
    let recs = EggmanEngine::parse_fasta(&example("demos"));
    assert_eq!(recs.len(), 5);

    let results: Vec<SequenceResult> = recs
        .iter()
        .map(|(h, s)| engine.classify_sequence(h, s, &cfg))
        .collect();

    assert_eq!(results[0].called_lineage, "B");
    assert_eq!(results[0].called_category, "Pure_Subtype");
    assert_eq!(results[0].epidemiological_formula, "B(Clonal)");
    assert_eq!(results[0].basal_formula, "B(Clonal)");

    assert_eq!(results[1].called_lineage, "CRF01_AE");
    assert_eq!(results[1].called_category, "Certified_CRF");
    assert_eq!(results[1].epidemiological_formula, "CRF01_AE(Clonal)");
    assert!(results[1].basal_formula.contains("?E"));

    assert_eq!(results[2].called_lineage, "CRF02_AG");
    assert_eq!(results[2].called_category, "Certified_CRF");
    assert_eq!(results[2].epidemiological_formula, "CRF02_AG(Clonal)");
    assert!(results[2].basal_formula.contains("G0"));

    assert_eq!(results[3].called_lineage, "CRF06_cpx");
    assert_eq!(results[3].called_category, "Certified_CRF");
    assert_eq!(results[3].epidemiological_formula, "CRF06_cpx(Clonal)");

    assert_eq!(results[4].called_lineage, "CRF07_BC");
    assert_eq!(results[4].called_category, "Certified_CRF");
    assert_eq!(results[4].epidemiological_formula, "CRF07_BC(Clonal)");
    assert!(results[4].basal_formula.starts_with("C[B:"));
}

#[test]
fn crf01_02_and_subtype_g_examples_resolve_accurately() {
    let engine = EggmanEngine::new();
    let cfg = ScreeningConfig::default();

    for (hdr, seq) in EggmanEngine::parse_fasta(&example("crf01_02")) {
        let res = engine.classify_sequence(&hdr, &seq, &cfg);
        if hdr.contains("01_AE") {
            assert_eq!(res.called_lineage, "CRF01_AE", "failed on {}", hdr);
            assert_eq!(res.called_category, "Certified_CRF");
        } else if hdr.contains("02_AG") {
            assert_eq!(res.called_lineage, "CRF02_AG", "failed on {}", hdr);
            assert_eq!(res.called_category, "Certified_CRF");
        }
    }

    for (hdr, seq) in EggmanEngine::parse_fasta(&example("subtype_g")) {
        let res = engine.classify_sequence(&hdr, &seq, &cfg);
        assert_eq!(res.called_lineage, "G", "failed on {}", hdr);
        assert_eq!(res.called_category, "Ancestral_Recombinant_Subtype", "failed on {}", hdr);
    }
}

#[test]
fn pr_rt_partial_calls_crf01_ae_pure_b_and_crf01_b_urf() {
    let engine = EggmanEngine::new();
    let cfg = ScreeningConfig::default();

    // HXB2 PR-RT span: 2253..3554 (1-based HXB2) -> 0-based [2252..3554]
    let hxb2_full = load_tsv_seq(PROFILES, "HXB2");
    let cm240_full = load_tsv_seq(PROFILES, "CM240");
    let crf01_founder = load_tsv_seq(FOUNDERS, "CRF01_AE");

    let hxb2_pr_rt = strip_gaps(&hxb2_full[2252..3554]);
    let cm240_pr_rt = strip_gaps(&cm240_full[2252..3554]);

    let r_b = engine.classify_sequence("amplicon_1", &hxb2_pr_rt, &cfg);
    assert_eq!(r_b.locus_scope, "pol");
    assert_eq!(r_b.called_lineage, "B");
    assert_eq!(r_b.called_category, "Pure_Subtype");
    assert_eq!(r_b.epidemiological_formula, "B(Clonal)@pol");

    let r_01 = engine.classify_sequence("amplicon_2", &cm240_pr_rt, &cfg);
    assert_eq!(r_01.locus_scope, "pol");
    assert_eq!(r_01.called_lineage, "CRF01_AE");
    assert_eq!(r_01.called_category, "Certified_CRF");
    assert_eq!(r_01.epidemiological_formula, "CRF01_AE(Clonal)@pol");

    // Synthetic second-generation URF in PR-RT: 5' CRF01_AE (2253..2900) + 3' B (2901..3554)
    let chimera_pol = format!(
        "{}{}",
        strip_gaps(&crf01_founder[2252..2900]),
        strip_gaps(&hxb2_full[2900..3554])
    );
    let r_urf = engine.classify_sequence("amplicon_3", &chimera_pol, &cfg);
    assert_eq!(r_urf.locus_scope, "pol");
    assert_eq!(r_urf.called_lineage, "URF");
    assert_eq!(r_urf.called_category, "URF");
    assert!(
        r_urf.epidemiological_segments.len() >= 2,
        "expected >=2 epi segments, got {:?}",
        r_urf.epidemiological_segments
    );
    assert!(
        r_urf.epidemiological_formula.contains("CRF01_AE")
            && r_urf.epidemiological_formula.contains('B'),
        "expected CRF01_AE and B in epidemiological_formula, got {}",
        r_urf.epidemiological_formula
    );
}

#[test]
fn classifier_is_blind_to_adversarial_fasta_headers() {
    let engine = EggmanEngine::new();
    let cfg = ScreeningConfig::default();
    let recs = EggmanEngine::parse_fasta(&example("demos"));
    let (_, hxb2_seq) = &recs[0];
    let (_, cm240_seq) = &recs[1];

    let r1_orig = engine.classify_sequence("Ref.01_AE.TH.90.CM240", cm240_seq, &cfg);
    let r1_adv = engine.classify_sequence("Pure_B_Fake_Header|panel=P1A|truth=B", cm240_seq, &cfg);
    assert_eq!(r1_orig.called_lineage, r1_adv.called_lineage);
    assert_eq!(r1_orig.called_category, r1_adv.called_category);
    assert_eq!(r1_orig.epidemiological_formula, r1_adv.epidemiological_formula);
    assert_eq!(r1_orig.basal_formula, r1_adv.basal_formula);
    assert_eq!(r1_orig.hxb2_barcode, r1_adv.hxb2_barcode);

    let r0_orig = engine.classify_sequence("Ref.B.FR.83.HXB2", hxb2_seq, &cfg);
    let r0_adv = engine.classify_sequence("CRF01_AE_Fake_Header|truth=CRF01_AE", hxb2_seq, &cfg);
    assert_eq!(r0_orig.called_lineage, r0_adv.called_lineage);
    assert_eq!(r0_orig.called_category, r0_adv.called_category);
    assert_eq!(r0_orig.epidemiological_formula, r0_adv.epidemiological_formula);
}
