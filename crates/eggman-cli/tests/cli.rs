use std::io::{BufRead, BufReader};
use std::path::PathBuf;
use std::process::{Command, Stdio};

fn bin() -> Command {
    Command::new(env!("CARGO_BIN_EXE_eggman"))
}

fn demos() -> String {
    concat!(env!("CARGO_MANIFEST_DIR"), "/../../docs/examples/demos.fasta").to_string()
}

fn tmp(name: &str) -> PathBuf {
    let d = std::env::temp_dir().join(format!("eggman_cli_{}", std::process::id()));
    std::fs::create_dir_all(&d).unwrap();
    d.join(name)
}

#[test]
fn no_args_exit_2() {
    let o = bin().output().unwrap();
    assert_eq!(o.status.code(), Some(2));
    assert!(!o.stderr.is_empty());
}

#[test]
fn empty_file_exit_1() {
    let p = tmp("empty.fasta");
    std::fs::write(&p, "").unwrap();
    let o = bin().arg(&p).output().unwrap();
    assert_eq!(o.status.code(), Some(1));
    assert!(String::from_utf8_lossy(&o.stderr).contains("no FASTA records found"));
}

#[test]
fn nan_p_thresh_exit_2() {
    let o = bin().args([&demos(), "--p-thresh", "nan"]).output().unwrap();
    assert_eq!(o.status.code(), Some(2));
}

#[test]
fn unreadable_input_exit_2() {
    let o = bin().arg("/nonexistent/in.fasta").output().unwrap();
    assert_eq!(o.status.code(), Some(2));
}

#[test]
fn unwritable_out_exit_2() {
    let o = bin().args([&demos(), "--out", "/nonexistent/x.csv"]).output().unwrap();
    assert_eq!(o.status.code(), Some(2));
    assert!(!String::from_utf8_lossy(&o.stdout).contains("EGGMAN"));
}

#[test]
fn csv_quoting() {
    let p = tmp("q.fasta");
    std::fs::write(&p, ">=\"x\"\nACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n").unwrap();
    let out = tmp("q.csv");
    let o = bin().arg(&p).arg("--out").arg(&out).arg("-q").output().unwrap();
    assert_eq!(o.status.code(), Some(0), "{}", String::from_utf8_lossy(&o.stderr));
    let csv = std::fs::read_to_string(&out).unwrap();
    let line = csv.lines().nth(1).unwrap();
    assert!(line.starts_with("\"'=\"\"x\"\"\","), "{}", line);
}

#[test]
fn broken_pipe_no_panic() {
    let out = tmp("bp.json");
    let mut c = bin()
        .arg(demos())
        .arg("--json")
        .arg(&out)
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .unwrap();
    let mut rd = BufReader::new(c.stdout.take().unwrap());
    let mut line = String::new();
    rd.read_line(&mut line).unwrap();
    drop(rd);
    let st = c.wait().unwrap();
    assert_ne!(st.code(), Some(101));
    assert_eq!(st.code(), Some(0));
    assert!(out.exists());
}

#[test]
fn demos_five_results() {
    let out = tmp("demos.json");
    let o = bin().arg(demos()).arg("--json").arg(&out).arg("-q").output().unwrap();
    assert_eq!(o.status.code(), Some(0));
    let v: serde_json::Value = serde_json::from_str(&std::fs::read_to_string(&out).unwrap()).unwrap();
    assert_eq!(v.as_array().unwrap().len(), 5);
}

#[test]
fn export_default_pack_writes_hiv1_spec() {
    let out = tmp("default_pack_spec.json");
    let o = bin()
        .arg("--export-default-pack")
        .arg(&out)
        .output()
        .unwrap();
    assert_eq!(o.status.code(), Some(0));
    let raw = std::fs::read_to_string(&out).unwrap();
    let v: serde_json::Value = serde_json::from_str(&raw).unwrap();
    assert_eq!(v["pathogen_name"], "HIV-1 Group M");
    assert_eq!(v["reference_name"], "HXB2");
    assert_eq!(v["ref_full_len"], 9719);
    let _ = std::fs::remove_file(&out);
}

#[test]
fn cli_pack_flag_classifies_influenza_reassortant() {
    use eggman_core::{GeneLandmark, PathogenSpec, TrunkSpec};

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
    let trunk = |name: &str| TrunkSpec {
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
    };
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
        landmarks: segs
            .iter()
            .map(|&(name, start, end)| GeneLandmark {
                short_name: name.to_string(),
                parent_gene: None,
                parent_min_span: 0,
                start,
                end,
            })
            .collect(),
        locus_scopes: Vec::new(),
        basal_trunks: vec![trunk("H1N1pdm09"), trunk("H5N1")],
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

    const BASES: [u8; 4] = [b'A', b'C', b'G', b'T'];
    let mut x: u64 = 0x9E37_79B9_7F4A_7C15;
    let mut anc = Vec::with_capacity(total_len);
    for _ in 0..total_len {
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        anc.push(BASES[(x as usize) & 3]);
    }
    let mutate = |phase: usize, shift: usize| -> Vec<u8> {
        anc.iter()
            .enumerate()
            .map(|(i, &b)| {
                if i % 16 == phase {
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
    };
    let h1 = mutate(0, 1);
    let h5 = mutate(8, 2);

    let mut q62 = h1.clone();
    q62[(6916 - 1)..8693].copy_from_slice(&h5[(6916 - 1)..8693]);
    q62[(10259 - 1)..11671].copy_from_slice(&h5[(10259 - 1)..11671]);

    let pack_dir = std::env::temp_dir().join(format!("eggman_cli_flupack_{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&pack_dir);
    std::fs::create_dir_all(&pack_dir).unwrap();
    std::fs::write(
        pack_dir.join("manifest.json"),
        serde_json::to_string_pretty(&spec).unwrap(),
    )
    .unwrap();
    std::fs::write(
        pack_dir.join("profiles.fasta"),
        format!(
            ">H1N1pdm09\n{}\n>H5N1\n{}\n",
            String::from_utf8_lossy(&h1),
            String::from_utf8_lossy(&h5),
        ),
    )
    .unwrap();

    let query_fa = tmp("flu_query.fasta");
    let out_json = tmp("flu_out.json");
    std::fs::write(
        &query_fa,
        format!(">Flu6plus2\n{}\n", String::from_utf8_lossy(&q62)),
    )
    .unwrap();

    let o = bin()
        .arg(&query_fa)
        .arg("--pack")
        .arg(&pack_dir)
        .arg("--json")
        .arg(&out_json)
        .output()
        .unwrap();

    let _ = std::fs::remove_dir_all(&pack_dir);
    let _ = std::fs::remove_file(&query_fa);

    assert_eq!(o.status.code(), Some(0));
    let stderr = String::from_utf8_lossy(&o.stderr);
    assert!(stderr.contains("Influenza A"), "stderr was: {}", stderr);
    assert!(stderr.contains("FluA-8seg"), "stderr was: {}", stderr);

    let v: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(&out_json).unwrap()).unwrap();
    let _ = std::fs::remove_file(&out_json);
    let arr = v.as_array().unwrap();
    assert_eq!(arr.len(), 1);
    assert_eq!(arr[0]["called_lineage"], "URF");
    assert_eq!(arr[0]["structural_formula"], "H1N1pdm09[H5N1:HA+NA]");
    assert_eq!(arr[0]["pack_version"], "test");
    assert_eq!(arr[0]["snapshot_date"], "2026-11-16");
}

#[test]
fn export_pack_and_builtin_id_flag_work() {
    let p = tmp("exported_iav_pack.json");
    let o = bin()
        .arg("--pack")
        .arg("iav")
        .arg("--export-pack")
        .arg(&p)
        .output()
        .unwrap();
    assert_eq!(o.status.code(), Some(0));
    let v: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(&p).unwrap()).unwrap();
    let _ = std::fs::remove_file(&p);
    assert_eq!(v["spec"]["reference_name"], "FluA-2.3.4.4b-8seg");
    assert_eq!(v["spec"]["pack_version"], "2026.11-chronaeon-h5nx-v2");
    assert!(v["profiles_fasta"].as_str().unwrap().contains(">EA1"));
}

#[test]
fn cli_iav_stitches_unconcatenated_8segment_fasta() {
    let raw_iav = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../../docs/examples/iav_raw_8segments.fasta"
    );
    let out_json = tmp("iav_stitched_cli.json");
    let o = bin()
        .arg(raw_iav)
        .arg("--pack")
        .arg("iav")
        .arg("--json")
        .arg(&out_json)
        .arg("-q")
        .output()
        .unwrap();
    assert_eq!(
        o.status.code(),
        Some(0),
        "stderr: {}",
        String::from_utf8_lossy(&o.stderr)
    );
    let v: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(&out_json).unwrap()).unwrap();
    let _ = std::fs::remove_file(&out_json);
    let arr = v.as_array().unwrap();
    assert_eq!(arr.len(), 2);
    assert_eq!(arr[0]["called_lineage"], "B3.13");
    assert_eq!(arr[1]["called_lineage"], "D1.1");
    assert_eq!(arr[0]["pack_version"], "2026.11-chronaeon-h5nx-v2");
}

#[test]
fn cli_temporal_guard_and_query_year_override() {
    let out_json = tmp("demos_temporal.json");
    let o = bin()
        .arg(demos())
        .arg("--temporal-guard")
        .arg("--query-year")
        .arg("1985.0")
        .arg("--temporal-slack-yr")
        .arg("1.0")
        .arg("--json")
        .arg(&out_json)
        .arg("-q")
        .output()
        .unwrap();
    assert_eq!(o.status.code(), Some(0), "stderr: {}", String::from_utf8_lossy(&o.stderr));
    let v: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(&out_json).unwrap()).unwrap();
    let _ = std::fs::remove_file(&out_json);
    let arr = v.as_array().unwrap();
    assert_eq!(arr.len(), 5);
    for rec in arr {
        assert_eq!(rec["query_year"].as_f64().unwrap(), 1985.0);
        let masked: Vec<&str> = rec["masked_future_trunks"]
            .as_array()
            .unwrap()
            .iter()
            .filter_map(|x| x.as_str())
            .collect();
        // K (t_mrca=1996.0) exceeds 1985 + 1.0 = 1986.0 and must be masked
        assert!(masked.contains(&"K"), "expected K masked in {:?}", masked);
    }

    // Invalid query year must exit cleanly with code 2
    let bad = bin()
        .arg(demos())
        .arg("--temporal-guard")
        .arg("--query-year")
        .arg("1800.0")
        .arg("-q")
        .output()
        .unwrap();
    assert_eq!(bad.status.code(), Some(2));
}

