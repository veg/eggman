use clap::Parser;
use rayon::prelude::*;
use eggman_core::{
    EggmanEngine, PathogenProfilePack, PathogenSpec, ScreeningConfig, SequenceResult,
    TemporalFloorMode, TemporalGuardConfig,
};
use std::fs;
use std::io::{self, Write};
use std::process::exit;
use std::path::PathBuf;
use std::time::Instant;

#[derive(Parser, Debug)]
#[command(
    name = "eggman",
    version = "0.2.0",
    about = "EGGMAN: Exact Genomic Genealogy and Mosaic Ancestry Nomenclature across recombining and reassorting RNA viruses (--pack)"
)]
struct Args {
    /// Input FASTA file (aligned or unaligned, full-length or partial amplicon)
    #[arg(short, long)]
    fasta: Option<PathBuf>,

    /// Positional FASTA path if --fasta is omitted
    #[arg(index = 1)]
    input_pos: Option<PathBuf>,

    /// Custom pathogen profile pack (built-in id `iav`/`sc2`/`sarbeco`/`hiv1`, JSON bundle file, or directory containing manifest.json)
    #[arg(long)]
    pack: Option<PathBuf>,

    /// Export a self-contained PathogenProfilePack JSON bundle (from --pack <id_or_dir>, or default hiv1) to PATH and exit
    #[arg(long)]
    export_pack: Option<PathBuf>,

    /// Export the embedded HIV-1 Group M PathogenSpec template JSON to PATH and exit
    #[arg(long)]
    export_default_pack: Option<PathBuf>,

    /// Output CSV file path for tabular results
    #[arg(short, long)]
    out: Option<PathBuf>,

    /// Output JSON file path for full detailed segment/breakpoint report
    #[arg(long)]
    json: Option<PathBuf>,

    /// Maximum number of sequences to process (0 = all)
    #[arg(short, long, default_value_t = 0)]
    limit: usize,

    /// Template-switch penalty gamma (default 7.0 if omitted)
    #[arg(long)]
    gamma: Option<f64>,

    /// Minimum segment length in reference nucleotides (auto-scales to 120 full / 150 partial if omitted)
    #[arg(long)]
    min_span: Option<usize>,

    /// Minimum subtype-informative SNPs per segment
    #[arg(long, default_value_t = 5)]
    min_snps: usize,

    /// Fisher exact test p-value threshold
    #[arg(long, default_value_t = 1e-4)]
    p_thresh: f64,

    /// Disable masking of hypervariable loops
    #[arg(long, default_value_t = false)]
    no_loop_mask: bool,

    /// Optional custom Tier-2 certified profiles TSV path (overrides embedded HIV-1 table)
    #[arg(long, visible_alias = "certified-profiles")]
    profiles: Option<PathBuf>,

    /// Enable non-lookahead temporal donor admissibility guarding
    #[arg(long, default_value_t = false)]
    temporal_guard: bool,

    /// Explicit decimal query year override for --temporal-guard (default: parsed from FASTA header)
    #[arg(long)]
    query_year: Option<f64>,

    /// Temporal lookahead tolerance in years for --temporal-guard (default: pack setting in sampled mode, or 0.0 in submission mode)
    #[arg(long)]
    temporal_slack_yr: Option<f64>,

    /// Temporal floor reference mode for --temporal-guard (`sampled` or `submission`)
    #[arg(long, default_value = "sampled")]
    temporal_mode: String,

    /// Segmented virus FASTA grouping mode (`auto`, `single_isolate`, `per_record`)
    #[arg(long, default_value = "auto")]
    seg_mode: String,

    /// Quiet mode (suppress per-sequence console report)
    #[arg(short, long, default_value_t = false)]
    quiet: bool,
}

fn die(code: i32, msg: &str) -> ! {
    eprintln!("eggman: {}", msg);
    exit(code)
}

/// RFC-4180 quoted field; `guard` neutralises spreadsheet formula prefixes.
fn csv_quote(s: &str, guard: bool) -> String {
    let mut out = String::with_capacity(s.len() + 2);
    out.push('"');
    if guard && s.starts_with(['=', '+', '-', '@']) {
        out.push('\'');
    }
    out.push_str(&s.replace('"', "\"\""));
    out.push('"');
    out
}

fn main() {
    let args = Args::parse();

    if let Some(export_path) = &args.export_pack {
        let pack = if let Some(pack_arg) = &args.pack {
            PathogenProfilePack::from_path(pack_arg)
                .unwrap_or_else(|e| die(2, &format!("failed to load pathogen pack {:?}: {}", pack_arg, e)))
        } else {
            PathogenProfilePack::hiv1_group_m_default()
        };
        let json_str = serde_json::to_string_pretty(&pack)
            .unwrap_or_else(|e| die(2, &format!("failed to serialize PathogenProfilePack: {}", e)));
        fs::write(export_path, json_str + "\n")
            .unwrap_or_else(|e| die(2, &format!("failed to write pack bundle {:?}: {}", export_path, e)));
        eprintln!(
            "Exported PathogenProfilePack bundle ('{}', version '{}') to {:?}",
            pack.spec.pathogen_name, pack.spec.pack_version, export_path
        );
        return;
    }

    if let Some(export_path) = &args.export_default_pack {
        let json_str = serde_json::to_string_pretty(&PathogenSpec::hiv1_group_m_default())
            .unwrap_or_else(|e| die(2, &format!("failed to serialize default PathogenSpec: {}", e)));
        fs::write(export_path, json_str + "\n")
            .unwrap_or_else(|e| die(2, &format!("failed to write template {:?}: {}", export_path, e)));
        eprintln!(
            "Exported default HIV-1 Group M PathogenSpec manifest template to {:?}",
            export_path
        );
        return;
    }

    let fasta_path = match args.fasta.or(args.input_pos) {
        Some(p) => p,
        None => die(2, "please provide an input FASTA file via `<file.fasta>` or `--fasta <file.fasta>`"),
    };

    if args.pack.is_some() && args.profiles.is_some() {
        die(2, "cannot combine `--pack` with `--profiles` (specify certified profiles inside the pack)");
    }

    let engine_cfg = ScreeningConfig {
        gamma: args.gamma,
        min_span_nt: args.min_span,
        min_informative_snps: args.min_snps,
        p_fisher_thresh: args.p_thresh,
        mask_hypervariable_loops: !args.no_loop_mask,
    };
    if let Err(e) = engine_cfg.validate() {
        die(2, &format!("invalid configuration: {}", e));
    }

    let temporal_mode = match args.temporal_mode.trim().to_ascii_lowercase().as_str() {
        "sampled" => TemporalFloorMode::Sampled,
        "submission" => TemporalFloorMode::Submission,
        other => die(
            2,
            &format!(
                "invalid --temporal-mode '{}' (expected `sampled` or `submission`)",
                other
            ),
        ),
    };
    let temporal_cfg = TemporalGuardConfig {
        enabled: args.temporal_guard,
        query_year: args.query_year,
        slack_yr: args.temporal_slack_yr,
        mode: temporal_mode,
    };
    if let Err(e) = temporal_cfg.validate() {
        die(2, &format!("invalid temporal guard configuration: {}", e));
    }

    let t0 = Instant::now();
    let bytes = fs::read(&fasta_path)
        .unwrap_or_else(|e| die(2, &format!("failed to read FASTA file {:?}: {}", fasta_path, e)));
    let raw_text = String::from_utf8(bytes)
        .unwrap_or_else(|_| die(2, &format!("FASTA file {:?} is not valid UTF-8", fasta_path)));

    let mut csv_file = args.out.as_ref().map(|p| {
        fs::File::create(p).unwrap_or_else(|e| die(2, &format!("cannot create CSV {:?}: {}", p, e)))
    });
    let mut json_file = args.json.as_ref().map(|p| {
        fs::File::create(p).unwrap_or_else(|e| die(2, &format!("cannot create JSON {:?}: {}", p, e)))
    });

    let mut records = EggmanEngine::parse_fasta(&raw_text);
    if records.is_empty() {
        die(1, "no FASTA records found");
    }

    let engine = if let Some(pack_path) = &args.pack {
        EggmanEngine::from_pack_path(pack_path)
            .unwrap_or_else(|e| die(2, &format!("failed to load pathogen pack {:?}: {}", pack_path, e)))
    } else if let Some(p) = &args.profiles {
        let p_bytes = fs::read(p)
            .unwrap_or_else(|e| die(2, &format!("failed to read profiles TSV {:?}: {}", p, e)));
        let p_text = String::from_utf8(p_bytes)
            .unwrap_or_else(|_| die(2, &format!("profiles TSV {:?} is not valid UTF-8", p)));
        EggmanEngine::with_certified_profiles_tsv(&p_text)
    } else {
        EggmanEngine::new()
    };
    if !args.quiet {
        eprintln!(
            "[eggman] pathogen='{}' reference='{}' basal_trunks={} certified_lineages={}",
            engine.pathogen_name(),
            engine.reference_name(),
            engine.basal_trunk_count(),
            engine.certified_lineage_count()
        );
    }
    let cfg = engine_cfg;

    let (results, total_records): (Vec<SequenceResult>, usize) = if engine.spec().is_segmented
        && !args.seg_mode.eq_ignore_ascii_case("per_record")
        && records
            .iter()
            .any(|(_, s)| s.chars().filter(|&c| c != '|').count() < 8_000)
    {
        let (_, mut seg_results, _) = engine.classify_fasta_with_segmentation_temporal(
            &raw_text,
            &args.seg_mode,
            &cfg,
            &temporal_cfg,
        );
        let n_isolates = seg_results.len();
        if args.limit > 0 && seg_results.len() > args.limit {
            seg_results.truncate(args.limit);
        }
        (seg_results, n_isolates)
    } else {
        let n_recs = records.len();
        if args.limit > 0 && records.len() > args.limit {
            records.truncate(args.limit);
        }
        let res = records
            .par_iter()
            .map(|(hdr, seq)| engine.classify_sequence_with_temporal(hdr, seq, &cfg, &temporal_cfg))
            .collect();
        (res, n_recs)
    };

    let elapsed = t0.elapsed().as_secs_f64();

    if let Some(f) = csv_file.as_mut() {
        let mut csv_lines = Vec::with_capacity(results.len() + 1);
        csv_lines.push(
            "sequence_id,raw_length_nt,mapped_coding_nt,hxb2_span_start,hxb2_span_end,locus_scope,alignment_mode,called_lineage,called_category,call_identity_pct,call_stem_snps,call_p_value,epidemiological_formula,epidemiological_barcode,basal_formula,basal_barcode,descriptive_name,structural_formula,backbone_trunk,n_segments,n_breakpoints,hxb2_barcode,ci_bracketed_barcode,compendium_matches,query_year,masked_future_trunks".to_string(),
        );
        for r in &results {
            let matches_joined = r
                .matches
                .iter()
                .map(|m| format!("{} ({})", m.label, m.descriptive_name))
                .collect::<Vec<_>>()
                .join("; ");
            let qy_str = r
                .query_year
                .map(|y| format!("{:.4}", y))
                .unwrap_or_default();
            let masked_joined = r.masked_future_trunks.join(";");
            let (id_str, stem_str, p_str) = if let Some(sup) = &r.call_support {
                (
                    format!("{:.1}", sup.identity_pct),
                    format!("{}", sup.net_stem_snps),
                    format!("{:.2e}", sup.p_value),
                )
            } else {
                (String::new(), String::new(), String::new())
            };
            csv_lines.push(format!(
                "{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{}",
                csv_quote(&r.sequence_id, true),
                r.raw_length_nt,
                r.mapped_coding_nt,
                r.hxb2_span_start,
                r.hxb2_span_end,
                r.locus_scope,
                csv_quote(&r.alignment_mode, false),
                csv_quote(&r.called_lineage, false),
                csv_quote(&r.called_category, false),
                id_str,
                stem_str,
                p_str,
                csv_quote(&r.epidemiological_formula, false),
                csv_quote(&r.epidemiological_barcode, false),
                csv_quote(&r.basal_formula, false),
                csv_quote(&r.basal_barcode, false),
                csv_quote(&r.descriptive_name, false),
                csv_quote(&r.structural_formula, false),
                r.backbone_trunk,
                r.n_segments,
                r.n_breakpoints,
                csv_quote(&r.hxb2_barcode, false),
                csv_quote(&r.ci_bracketed_barcode, false),
                csv_quote(&matches_joined, true),
                qy_str,
                csv_quote(&masked_joined, false)
            ));
        }
        let path = args.out.as_ref().unwrap();
        f.write_all((csv_lines.join("\n") + "\n").as_bytes())
            .unwrap_or_else(|e| die(2, &format!("failed to write CSV {:?}: {}", path, e)));
    }

    if let Some(f) = json_file.as_mut() {
        let json_str = serde_json::to_string_pretty(&results)
            .unwrap_or_else(|e| die(2, &format!("failed to serialize JSON: {}", e)));
        let path = args.json.as_ref().unwrap();
        f.write_all(json_str.as_bytes())
            .unwrap_or_else(|e| die(2, &format!("failed to write JSON {:?}: {}", path, e)));
    }

    let mut out = String::new();
    if !args.quiet {
        let n_label = if results.len() < total_records {
            format!("{} of {} records", results.len(), total_records)
        } else {
            format!("{} sequences", results.len())
        };
        out.push_str("================================================================================================\n");
        out.push_str(&format!(
            "EGGMAN {} SUBTYPING & RECOMBINATION CLASSIFIER (ref: {}, N={}, {:.3}s)\n",
            engine.pathogen_name().to_uppercase(),
            engine.reference_name(),
            n_label,
            elapsed
        ));
        out.push_str("================================================================================================\n");

        for r in &results {
            out.push_str(&format!("\n> {}\n", r.sequence_id));
            out.push_str(&format!(
                "  Alignment / Scope   : {} | {} ({} {}..{}, {} mapped coding nt)\n",
                r.alignment_mode,
                r.locus_scope,
                engine.reference_name(),
                r.hxb2_span_start,
                r.hxb2_span_end,
                r.mapped_coding_nt
            ));
            if let Some(sup) = &r.call_support {
                out.push_str(&format!(
                    "  Called Lineage      : {} [{}] ({:.1}% id, {:+} stem SNPs, p={:.2e})\n",
                    r.called_lineage, r.called_category, sup.identity_pct, sup.net_stem_snps, sup.p_value
                ));
            } else {
                out.push_str(&format!(
                    "  Called Lineage      : {} [{}]\n",
                    r.called_lineage, r.called_category
                ));
            }
            out.push_str(&format!("  Epidemiological Map : {}\n", r.epidemiological_formula));
            out.push_str(&format!("  Epidemiological BC  : {}\n", r.epidemiological_barcode));
            out.push_str(&format!("  Basal Ancestral Map : {}\n", r.descriptive_name));
            out.push_str(&format!("  Basal Ref Barcode   : {}\n", r.hxb2_barcode));
            out.push_str(&format!("  CI-Bracketed Map    : {}\n", r.ci_bracketed_barcode));
            if r.query_year.is_some() || !r.masked_future_trunks.is_empty() {
                let qy = r
                    .query_year
                    .map(|y| format!("{:.4}", y))
                    .unwrap_or_else(|| "none".to_string());
                let masked = if r.masked_future_trunks.is_empty() {
                    "none".to_string()
                } else {
                    r.masked_future_trunks.join(",")
                };
                out.push_str(&format!("  Temporal Guard      : query_year={} | masked_trunks=[{}]\n", qy, masked));
            }
            let match_strs: Vec<String> = r
                .matches
                .iter()
                .take(4)
                .map(|m| format!("{} [{}]", m.label, m.match_type))
                .collect();
            out.push_str(&format!("  Compendium Match    : {}\n", match_strs.join("; ")));
            if !r.breakpoints.is_empty() {
                let bp_strs: Vec<String> = r
                    .breakpoints
                    .iter()
                    .map(|b| {
                        format!(
                            "{} {} ({}->{}, CI {}-{}, p={:.2e})",
                            engine.reference_name(),
                            b.bp_hxb2,
                            b.left_trunk,
                            b.right_trunk,
                            b.ci_hxb2_low,
                            b.ci_hxb2_high,
                            b.p_fisher
                        )
                    })
                    .collect();
                out.push_str(&format!("  Verified Junctions  : {}\n", bp_strs.join("; ")));
            }
        }
    }
    if let Some(p) = &args.out {
        out.push_str(&format!("\n[+] Wrote CSV results ({} sequences) to {:?}\n", results.len(), p));
    }
    if let Some(p) = &args.json {
        out.push_str(&format!("[+] Wrote JSON report to {:?}\n", p));
    }

    // A closed pipe (e.g. `| head`) is not an error: files are already written.
    let stdout = io::stdout();
    let mut lock = stdout.lock();
    if let Err(e) = lock.write_all(out.as_bytes()).and_then(|_| lock.flush()) {
        if e.kind() != io::ErrorKind::BrokenPipe {
            die(2, &format!("failed to write to stdout: {}", e));
        }
    }
}
