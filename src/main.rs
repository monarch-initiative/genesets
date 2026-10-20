use std::{
    collections::HashMap,
    fs::File,
    io::{self, Write},
    path::{Path, PathBuf},
};

use anyhow::{Context, Result, bail};
use clap::{Args, Parser, Subcommand};
use genesets_rs::{
    Correction, EnrichmentOptions, GeneUniverse, NamedSet, RawSet, SetFormat,
    compare::{
        DiffOptions, DiffSummaryContext, ResultFormat, compare_records, load_result_records,
        resolve_format, summarize_diff,
    },
    enrichment::{enrich_matrix, enrich_one},
    gsea::{Algorithm, GseaOptions, Ranking, Statistic, gsea_matrix},
    index::{
        TermGraph, build_flat_sets, build_gene_list_bitset, build_ontology_sets, build_term_graph,
        union_sets,
    },
    io::{
        default_set_name, read_annotation_pairs, read_closure_pairs, read_gene_list, read_name_map,
        read_name_pairs, read_rankings, read_sets,
    },
    mgsa::{BetaPrior, MgsaOptions, MgsaResult, mgsa_matrix},
    output::{
        OutputFormat, OutputOptions, write_diff_rows, write_gsea_rows, write_mgsa_rows,
        write_parquet_diff_rows, write_parquet_rows, write_rows,
    },
};
use serde::{Deserialize, Serialize};

#[derive(Parser, Debug)]
#[command(
    author,
    version,
    about = "Fast ontology-aware gene set enrichment analysis"
)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand, Debug)]
enum Commands {
    /// Enrich one query gene set against ontology terms or flat target sets.
    Enrich(EnrichArgs),
    /// Enrich many query sets against many target sets.
    Matrix(MatrixArgs),
    /// Run an enrich or matrix job from YAML.
    Run(RunArgs),
    /// Compare two result tables by significance crossing.
    Compare(CompareArgs),
    /// Bayesian model-based gene set analysis (MGSA), as in Ontologizer and wtfgenes.
    Mgsa(MgsaArgs),
    /// Preranked gene set enrichment analysis (GSEA) over scored gene lists.
    Gsea(GseaArgs),
}

#[derive(Args, Clone, Debug)]
struct TargetArgs {
    /// Gene-to-term pair table. Columns: gene_curie, term_curie.
    #[arg(long)]
    annotations: Option<PathBuf>,

    /// Flat target gene sets. Pairwise columns are set_id, gene_curie.
    #[arg(long)]
    target_sets: Option<PathBuf>,

    /// Format for --target-sets.
    #[arg(long, value_enum, default_value_t = SetFormat::Auto)]
    target_format: SetFormat,

    /// Optional term CURIE/name table.
    #[arg(long)]
    terms: Option<PathBuf>,

    /// Optional reflexive child-to-ancestor closure table.
    #[arg(long)]
    closure: Option<PathBuf>,
}

#[derive(Args, Clone, Debug)]
struct EnrichArgs {
    #[command(flatten)]
    targets: TargetArgs,

    /// Optional gene CURIE/name or symbol table.
    #[arg(long)]
    gene_names: Option<PathBuf>,

    /// Query gene set file.
    #[arg(long)]
    sample: PathBuf,

    /// Format for --sample.
    #[arg(long, value_enum, default_value_t = SetFormat::Auto)]
    sample_format: SetFormat,

    /// Select one set from a multi-set sample file.
    #[arg(long)]
    sample_set: Option<String>,

    /// Name to use when --sample is a one-column list.
    #[arg(long, default_value = "sample")]
    sample_name: String,

    /// Optional background gene list. Defaults to all genes in target sets.
    #[arg(long)]
    background: Option<PathBuf>,

    /// Suppress rows with overlap below this count.
    #[arg(long, default_value_t = 1)]
    min_overlap: u32,

    /// Suppress rows with raw p-value above this cutoff.
    #[arg(long, value_parser = parse_probability)]
    max_p_value: Option<f64>,

    /// Suppress rows with adjusted p-value above this cutoff.
    #[arg(long, value_parser = parse_probability)]
    max_p_adjust: Option<f64>,

    /// Multiple-testing correction.
    #[arg(long, value_enum, default_value_t = Correction::Bonferroni)]
    correction: Correction,

    /// Write output here instead of stdout. Required for parquet output.
    #[arg(short, long)]
    output: Option<PathBuf>,

    /// Output serialization format.
    #[arg(long, value_enum, default_value_t = OutputFormat::Tsv)]
    output_format: OutputFormat,

    /// Include semicolon-separated overlapping gene ids and names.
    #[arg(long)]
    overlap_genes: bool,

    /// Rayon worker threads. Defaults to Rayon auto-detection.
    #[arg(long)]
    threads: Option<usize>,
}

#[derive(Args, Clone, Debug)]
struct MatrixArgs {
    #[command(flatten)]
    targets: TargetArgs,

    /// Optional gene CURIE/name or symbol table.
    #[arg(long)]
    gene_names: Option<PathBuf>,

    /// Query gene sets. Omit with --queries-from-targets for term-vs-term runs.
    #[arg(long)]
    queries: Option<PathBuf>,

    /// Format for --queries.
    #[arg(long, value_enum, default_value_t = SetFormat::Auto)]
    query_format: SetFormat,

    /// Select one query set from a multi-set query file.
    #[arg(long)]
    query_set: Option<String>,

    /// Use the target sets themselves as queries.
    #[arg(long)]
    queries_from_targets: bool,

    /// Optional background gene list. Defaults to all genes in target sets.
    #[arg(long)]
    background: Option<PathBuf>,

    /// Suppress rows with overlap below this count.
    #[arg(long, default_value_t = 1)]
    min_overlap: u32,

    /// Suppress rows with raw p-value above this cutoff.
    #[arg(long, value_parser = parse_probability)]
    max_p_value: Option<f64>,

    /// Suppress rows with adjusted p-value above this cutoff.
    #[arg(long, value_parser = parse_probability)]
    max_p_adjust: Option<f64>,

    /// Multiple-testing correction.
    #[arg(long, value_enum, default_value_t = Correction::Bonferroni)]
    correction: Correction,

    /// Write output here instead of stdout. Required for parquet output.
    #[arg(short, long)]
    output: Option<PathBuf>,

    /// Output serialization format.
    #[arg(long, value_enum, default_value_t = OutputFormat::Tsv)]
    output_format: OutputFormat,

    /// Include semicolon-separated overlapping gene ids and names.
    #[arg(long)]
    overlap_genes: bool,

    /// Rayon worker threads. Defaults to Rayon auto-detection.
    #[arg(long)]
    threads: Option<usize>,
}

#[derive(Args, Clone, Debug)]
struct MgsaArgs {
    #[command(flatten)]
    targets: TargetArgs,

    /// Optional gene CURIE/name or symbol table.
    #[arg(long)]
    gene_names: Option<PathBuf>,

    /// Query gene set file. Every set in a multi-set file is analyzed
    /// independently unless --sample-set picks one.
    #[arg(long)]
    sample: PathBuf,

    /// Format for --sample.
    #[arg(long, value_enum, default_value_t = SetFormat::Auto)]
    sample_format: SetFormat,

    /// Select one set from a multi-set sample file.
    #[arg(long)]
    sample_set: Option<String>,

    /// Name to use when --sample is a one-column list.
    #[arg(long, default_value = "sample")]
    sample_name: String,

    /// Optional background gene list. Defaults to all genes in target sets.
    #[arg(long)]
    background: Option<PathBuf>,

    /// Recorded MCMC samples per modeled term.
    #[arg(long, default_value_t = 100)]
    samples_per_term: u32,

    /// Burn-in MCMC samples per modeled term.
    #[arg(long, default_value_t = 10)]
    burn_per_term: u32,

    /// Mode of the Beta prior on the probability that a term is active.
    #[arg(long, default_value_t = 0.5, value_parser = parse_probability)]
    term_prob: f64,

    /// Pseudocounts of the term prior. 0 gives a uniform prior.
    #[arg(long, default_value_t = 0.0)]
    term_count: f64,

    /// Mode of the Beta prior on the false-positive rate.
    #[arg(long, default_value_t = 0.5, value_parser = parse_probability)]
    false_pos_prob: f64,

    /// Pseudocounts of the false-positive prior.
    #[arg(long, default_value_t = 0.0)]
    false_pos_count: f64,

    /// Mode of the Beta prior on the false-negative rate.
    #[arg(long, default_value_t = 0.5, value_parser = parse_probability)]
    false_neg_prob: f64,

    /// Pseudocounts of the false-negative prior.
    #[arg(long, default_value_t = 0.0)]
    false_neg_count: f64,

    /// Relative rate of moves that toggle one term.
    #[arg(long, default_value_t = 1.0)]
    flip_rate: f64,

    /// Relative rate of moves that swap an active term for a parent, child,
    /// or sibling. Needs --closure; ignored for flat target sets.
    #[arg(long, default_value_t = 1.0)]
    step_rate: f64,

    /// Relative rate of moves that swap an active term for any other term.
    #[arg(long, default_value_t = 1.0)]
    jump_rate: f64,

    /// Random seed. Each query derives its own stream from it.
    #[arg(long, default_value_t = 123_456_789)]
    seed: u64,

    /// Suppress rows with posterior below this cutoff.
    #[arg(long, default_value_t = 0.0, value_parser = parse_probability)]
    min_posterior: f64,

    /// Write output here instead of stdout.
    #[arg(short, long)]
    output: Option<PathBuf>,

    /// Output serialization format. Parquet is not supported for MGSA yet.
    #[arg(long, value_enum, default_value_t = OutputFormat::Tsv)]
    output_format: OutputFormat,

    /// Include semicolon-separated overlapping gene ids and names.
    #[arg(long)]
    overlap_genes: bool,

    /// Optional YAML file with per-query sampler diagnostics and posterior
    /// rate estimates.
    #[arg(long)]
    summary_output: Option<PathBuf>,

    /// Rayon worker threads. Defaults to Rayon auto-detection.
    #[arg(long)]
    threads: Option<usize>,
}

#[derive(Args, Clone, Debug)]
struct GseaArgs {
    #[command(flatten)]
    targets: TargetArgs,

    /// Optional gene CURIE/name or symbol table.
    #[arg(long)]
    gene_names: Option<PathBuf>,

    /// Scored gene list: two columns (gene, score), as in a GSEA .rnk file,
    /// or three (ranking_id, gene, score) for several rankings at once.
    /// Every ranked gene counts toward the ranking, annotated or not.
    #[arg(long)]
    ranks: PathBuf,

    /// Name for a two-column ranking. Defaults to the file stem.
    #[arg(long)]
    ranks_name: Option<String>,

    /// Select one ranking from a three-column file.
    #[arg(long)]
    ranks_set: Option<String>,

    /// How each target is scored against the ranking.
    #[arg(long, value_enum, default_value_t = Statistic::Ks)]
    statistic: Statistic,

    /// How the ontology graph is used: score terms independently (classic),
    /// drop genes of significant descendants (elim), or drop a child's genes
    /// from its parent when the child is at least as significant (weight01).
    /// elim and weight01 need --closure.
    #[arg(long, value_enum, default_value_t = Algorithm::Classic)]
    algorithm: Algorithm,

    /// Raw p-value at or below which elim removes a term's genes from its
    /// ancestors.
    #[arg(long, default_value_t = 0.01, value_parser = parse_probability)]
    elim_cutoff: f64,

    /// Exponent on |score| for the running-sum step at set members.
    /// 0 is the unweighted Kolmogorov-Smirnov statistic; 1 is the GSEA default.
    #[arg(long, default_value_t = 1.0)]
    weight: f64,

    /// Random gene sets drawn for each KS null.
    #[arg(long, default_value_t = 1000)]
    permutations: u32,

    /// Skip targets with fewer ranked members.
    #[arg(long, default_value_t = 15)]
    min_size: u32,

    /// Skip targets with more ranked members.
    #[arg(long, default_value_t = 500)]
    max_size: u32,

    /// Random seed. Each ranking and set size derives its own stream from it.
    #[arg(long, default_value_t = 123_456_789)]
    seed: u64,

    /// Suppress rows with raw p-value above this cutoff.
    #[arg(long, value_parser = parse_probability)]
    max_p_value: Option<f64>,

    /// Suppress rows with Benjamini-Hochberg adjusted p-value above this cutoff.
    #[arg(long, value_parser = parse_probability)]
    max_p_adjust: Option<f64>,

    /// Write output here instead of stdout.
    #[arg(short, long)]
    output: Option<PathBuf>,

    /// Output serialization format. Parquet is not supported for GSEA yet.
    #[arg(long, value_enum, default_value_t = OutputFormat::Tsv)]
    output_format: OutputFormat,

    /// Include semicolon-separated leading-edge gene ids and names.
    #[arg(long)]
    leading_edge_genes: bool,

    /// Rayon worker threads. Defaults to Rayon auto-detection.
    #[arg(long)]
    threads: Option<usize>,
}

#[derive(Args, Clone, Debug)]
struct RunArgs {
    /// YAML configuration file. Relative paths inside it resolve from its directory.
    config: PathBuf,
}

#[derive(Args, Clone, Debug)]
struct CompareArgs {
    /// Left result table, usually the older ontology or annotation version.
    #[arg(long)]
    left: PathBuf,

    /// Right result table, usually the newer ontology or annotation version.
    #[arg(long)]
    right: PathBuf,

    /// Format for --left.
    #[arg(long, value_enum, default_value_t = ResultFormat::Auto)]
    left_format: ResultFormat,

    /// Format for --right.
    #[arg(long, value_enum, default_value_t = ResultFormat::Auto)]
    right_format: ResultFormat,

    /// Adjusted p-value threshold used for significance crossing.
    #[arg(long, default_value_t = 0.05, value_parser = parse_probability)]
    p_adjust_cutoff: f64,

    /// Emit only gained/lost threshold crossings, excluding shared significant pairs.
    #[arg(long)]
    crossings_only: bool,

    /// Write diff output here instead of stdout. Required for parquet output.
    #[arg(short, long)]
    output: Option<PathBuf>,

    /// Output serialization format.
    #[arg(long, value_enum, default_value_t = OutputFormat::Tsv)]
    output_format: OutputFormat,

    /// Optional YAML metadata summary output.
    #[arg(long)]
    metadata_output: Option<PathBuf>,
}

#[derive(Clone, Debug)]
enum TargetSource {
    Ontology {
        term_names: Vec<(String, String)>,
        closure: Vec<(String, String)>,
        annotations: Vec<(String, String)>,
    },
    Flat {
        term_names: Vec<(String, String)>,
        sets: Vec<RawSet>,
    },
}

impl TargetSource {
    fn load(args: &TargetArgs) -> Result<Self> {
        match (&args.annotations, &args.target_sets) {
            (Some(_), Some(_)) => bail!("use either --annotations or --target-sets, not both"),
            (None, None) => bail!("one of --annotations or --target-sets is required"),
            (Some(annotations), None) => {
                let term_names = read_optional_name_pairs(args.terms.as_deref())?;
                let closure = read_optional_closure(args.closure.as_deref())?;
                let annotations = read_annotation_pairs(annotations)?;
                Ok(Self::Ontology {
                    term_names,
                    closure,
                    annotations,
                })
            }
            (None, Some(target_sets)) => {
                if args.closure.is_some() {
                    bail!(
                        "--closure requires --annotations; flat --target-sets do not use ontology closure"
                    );
                }
                let term_names = read_optional_name_pairs(args.terms.as_deref())?;
                let sets = read_sets(
                    target_sets,
                    args.target_format,
                    &default_set_name(target_sets),
                )?;
                Ok(Self::Flat { term_names, sets })
            }
        }
    }

    fn add_genes_to_universe(&self, universe: &mut GeneUniverse) {
        match self {
            Self::Ontology { annotations, .. } => {
                for (gene, _) in annotations {
                    universe.add_gene(gene);
                }
            }
            Self::Flat { sets, .. } => {
                for set in sets {
                    universe.add_genes(&set.genes);
                }
            }
        }
    }

    fn build_graph(&self, sets: &[NamedSet]) -> Option<TermGraph> {
        match self {
            Self::Ontology { closure, .. } if !closure.is_empty() => {
                Some(build_term_graph(sets, closure))
            }
            _ => None,
        }
    }

    fn build_sets(&self, universe: &GeneUniverse) -> Vec<NamedSet> {
        match self {
            Self::Ontology {
                term_names,
                closure,
                annotations,
            } => build_ontology_sets(term_names, closure, annotations, universe),
            Self::Flat { term_names, sets } => {
                let term_name_map: HashMap<String, String> = term_names.iter().cloned().collect();
                let mut sets = build_flat_sets(sets, universe);
                for set in &mut sets {
                    if set.name.is_none() {
                        set.name = term_name_map.get(&set.id).cloned();
                    }
                }
                sets
            }
        }
    }
}

#[derive(Copy, Clone, Debug, Deserialize)]
#[serde(rename_all = "kebab-case")]
enum ConfigMode {
    Enrich,
    Matrix,
}

#[derive(Debug, Deserialize)]
struct RunConfig {
    mode: Option<ConfigMode>,
    ontology: Option<OntologyConfig>,
    input: InputConfig,
    background: Option<BackgroundConfig>,
    gene_names: Option<PathBuf>,
    output: Option<PathBuf>,
    output_format: Option<OutputFormat>,
    min_overlap: Option<u32>,
    max_p_value: Option<f64>,
    max_p_adjust: Option<f64>,
    correction: Option<Correction>,
    overlap_genes: Option<bool>,
    threads: Option<usize>,
}

impl RunConfig {
    fn resolve_paths(&mut self, base: &Path) {
        if let Some(ontology) = &mut self.ontology {
            ontology.resolve_paths(base);
        }
        self.input.resolve_paths(base);
        if let Some(background) = &mut self.background {
            background.resolve_paths(base);
        }
        resolve_optional_config_path(base, &mut self.gene_names);
        resolve_optional_config_path(base, &mut self.output);
    }
}

#[derive(Clone, Debug, Default, Deserialize)]
struct OntologyConfig {
    terms: Option<PathBuf>,
    closure: Option<PathBuf>,
    annotations: Option<PathBuf>,
    gene_names: Option<PathBuf>,
}

impl OntologyConfig {
    fn resolve_paths(&mut self, base: &Path) {
        resolve_optional_config_path(base, &mut self.terms);
        resolve_optional_config_path(base, &mut self.closure);
        resolve_optional_config_path(base, &mut self.annotations);
        resolve_optional_config_path(base, &mut self.gene_names);
    }
}

#[derive(Clone, Debug, Default, Deserialize)]
struct InputConfig {
    sample: Option<PathBuf>,
    sample_format: Option<SetFormat>,
    sample_set: Option<String>,
    sample_name: Option<String>,
    queries: Option<PathBuf>,
    query_format: Option<SetFormat>,
    query_set: Option<String>,
    queries_from_targets: Option<bool>,
    targets: Option<PathBuf>,
    target_format: Option<SetFormat>,
}

impl InputConfig {
    fn resolve_paths(&mut self, base: &Path) {
        resolve_optional_config_path(base, &mut self.sample);
        resolve_optional_config_path(base, &mut self.queries);
        resolve_optional_config_path(base, &mut self.targets);
    }
}

#[derive(Clone, Debug, Deserialize)]
#[serde(untagged)]
enum BackgroundConfig {
    Path(PathBuf),
    Object { file: PathBuf },
}

impl BackgroundConfig {
    fn path(&self) -> &Path {
        match self {
            Self::Path(path) => path,
            Self::Object { file } => file,
        }
    }
}

impl BackgroundConfig {
    fn resolve_paths(&mut self, base: &Path) {
        match self {
            Self::Path(path) => resolve_config_path(base, path),
            Self::Object { file } => resolve_config_path(base, file),
        }
    }
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    match cli.command {
        Commands::Enrich(args) => run_with_thread_pool(args.threads, || execute_enrich(args)),
        Commands::Matrix(args) => run_with_thread_pool(args.threads, || execute_matrix(args)),
        Commands::Run(args) => execute_run(args),
        Commands::Compare(args) => execute_compare(args),
        Commands::Mgsa(args) => run_with_thread_pool(args.threads, || execute_mgsa(args)),
        Commands::Gsea(args) => run_with_thread_pool(args.threads, || execute_gsea(args)),
    }
}

fn execute_run(args: RunArgs) -> Result<()> {
    let config_file = File::open(&args.config)
        .with_context(|| format!("failed to open config {}", args.config.display()))?;
    let mut config: RunConfig = serde_yaml::from_reader(config_file)
        .with_context(|| format!("failed to parse config {}", args.config.display()))?;
    let config_base = args.config.parent().unwrap_or_else(|| Path::new("."));
    config.resolve_paths(config_base);
    let threads = config.threads;

    run_with_thread_pool(threads, || {
        let mode = config.mode.unwrap_or_else(|| {
            if config.input.sample.is_some() {
                ConfigMode::Enrich
            } else {
                ConfigMode::Matrix
            }
        });

        match mode {
            ConfigMode::Enrich => execute_enrich(config_to_enrich_args(config)?),
            ConfigMode::Matrix => execute_matrix(config_to_matrix_args(config)?),
        }
    })
}

fn resolve_config_path(base: &Path, path: &mut PathBuf) {
    if path.is_relative() {
        *path = base.join(&path);
    }
}

fn resolve_optional_config_path(base: &Path, path: &mut Option<PathBuf>) {
    if let Some(path) = path {
        resolve_config_path(base, path);
    }
}

fn config_to_target_args(config: &RunConfig) -> TargetArgs {
    let ontology = config.ontology.clone().unwrap_or_default();
    TargetArgs {
        annotations: ontology.annotations,
        target_sets: config.input.targets.clone(),
        target_format: config.input.target_format.unwrap_or(SetFormat::Auto),
        terms: ontology.terms,
        closure: ontology.closure,
    }
}

fn config_gene_names(config: &RunConfig) -> Option<PathBuf> {
    config
        .ontology
        .as_ref()
        .and_then(|ontology| ontology.gene_names.clone())
        .or_else(|| config.gene_names.clone())
}

fn config_background(config: &RunConfig) -> Option<PathBuf> {
    config
        .background
        .as_ref()
        .map(|background| background.path().to_path_buf())
}

fn config_to_enrich_args(config: RunConfig) -> Result<EnrichArgs> {
    let sample = config
        .input
        .sample
        .clone()
        .context("YAML enrich mode requires input.sample")?;
    Ok(EnrichArgs {
        targets: config_to_target_args(&config),
        gene_names: config_gene_names(&config),
        sample,
        sample_format: config.input.sample_format.unwrap_or(SetFormat::Auto),
        sample_set: config.input.sample_set.clone(),
        sample_name: config
            .input
            .sample_name
            .clone()
            .unwrap_or_else(|| "sample".to_owned()),
        background: config_background(&config),
        min_overlap: config.min_overlap.unwrap_or(1),
        max_p_value: config.max_p_value,
        max_p_adjust: config.max_p_adjust,
        correction: config.correction.unwrap_or(Correction::Bonferroni),
        output: config.output.clone(),
        output_format: config.output_format.unwrap_or(OutputFormat::Tsv),
        overlap_genes: config.overlap_genes.unwrap_or(false),
        threads: None,
    })
}

fn config_to_matrix_args(config: RunConfig) -> Result<MatrixArgs> {
    Ok(MatrixArgs {
        targets: config_to_target_args(&config),
        gene_names: config_gene_names(&config),
        queries: config.input.queries.clone(),
        query_format: config.input.query_format.unwrap_or(SetFormat::Auto),
        query_set: config.input.query_set.clone(),
        queries_from_targets: config.input.queries_from_targets.unwrap_or(false),
        background: config_background(&config),
        min_overlap: config.min_overlap.unwrap_or(1),
        max_p_value: config.max_p_value,
        max_p_adjust: config.max_p_adjust,
        correction: config.correction.unwrap_or(Correction::Bonferroni),
        output: config.output.clone(),
        output_format: config.output_format.unwrap_or(OutputFormat::Tsv),
        overlap_genes: config.overlap_genes.unwrap_or(false),
        threads: None,
    })
}

fn execute_enrich(args: EnrichArgs) -> Result<()> {
    let target_source = TargetSource::load(&args.targets)?;
    let mut query_sets = read_sets(&args.sample, args.sample_format, &args.sample_name)
        .with_context(|| format!("failed to read sample sets from {}", args.sample.display()))?;
    query_sets = filter_raw_sets(query_sets, args.sample_set.as_deref(), "sample")?;
    let query_set = require_single_set(query_sets, "sample", "--sample-set")?;

    let background_genes = read_optional_gene_list(args.background.as_deref())?;
    let gene_names = read_optional_name_map(args.gene_names.as_deref())?;

    let mut universe = GeneUniverse::new();
    target_source.add_genes_to_universe(&mut universe);
    universe.add_genes(&query_set.genes);
    if let Some(background_genes) = &background_genes {
        universe.add_genes(background_genes);
    }
    universe.apply_gene_names(&gene_names);

    let target_sets = target_source.build_sets(&universe);
    ensure_non_empty_targets(&target_sets)?;
    let query = build_flat_sets(&[query_set], &universe)
        .into_iter()
        .next()
        .expect("one query set was provided");
    let background = build_background(background_genes.as_deref(), &target_sets, &universe);

    let options = EnrichmentOptions {
        min_overlap: args.min_overlap,
        correction: args.correction,
        include_overlap_genes: args.overlap_genes,
        max_p_value: args.max_p_value,
        max_p_adjust: args.max_p_adjust,
    };
    let rows = enrich_one(&query, &target_sets, &background, &options);
    let queries = std::slice::from_ref(&query);
    write_result_rows(
        args.output.as_deref(),
        &rows,
        queries,
        &target_sets,
        &universe,
        OutputOptions {
            format: args.output_format,
            include_overlap_genes: args.overlap_genes,
            correction: args.correction,
        },
    )
}

fn execute_matrix(args: MatrixArgs) -> Result<()> {
    let target_source = TargetSource::load(&args.targets)?;
    let mut query_sets = match (&args.queries, args.queries_from_targets) {
        (Some(_), true) => bail!("use either --queries or --queries-from-targets, not both"),
        (None, false) => bail!("--queries is required unless --queries-from-targets is set"),
        (Some(path), false) => read_sets(path, args.query_format, &default_set_name(path))
            .with_context(|| format!("failed to read query sets from {}", path.display()))?,
        (None, true) => Vec::new(),
    };
    if !args.queries_from_targets {
        query_sets = filter_raw_sets(query_sets, args.query_set.as_deref(), "query")?;
    }

    let background_genes = read_optional_gene_list(args.background.as_deref())?;
    let gene_names = read_optional_name_map(args.gene_names.as_deref())?;

    let mut universe = GeneUniverse::new();
    target_source.add_genes_to_universe(&mut universe);
    for query_set in &query_sets {
        universe.add_genes(&query_set.genes);
    }
    if let Some(background_genes) = &background_genes {
        universe.add_genes(background_genes);
    }
    universe.apply_gene_names(&gene_names);

    let target_sets = target_source.build_sets(&universe);
    ensure_non_empty_targets(&target_sets)?;
    let query_sets = if args.queries_from_targets {
        target_sets.clone()
    } else {
        build_flat_sets(&query_sets, &universe)
    };
    let query_sets = filter_named_sets(query_sets, args.query_set.as_deref(), "query")?;
    if query_sets.is_empty() {
        bail!("no query sets were loaded");
    }
    let background = build_background(background_genes.as_deref(), &target_sets, &universe);

    let options = EnrichmentOptions {
        min_overlap: args.min_overlap,
        correction: args.correction,
        include_overlap_genes: args.overlap_genes,
        max_p_value: args.max_p_value,
        max_p_adjust: args.max_p_adjust,
    };
    let rows = enrich_matrix(&query_sets, &target_sets, &background, &options);
    write_result_rows(
        args.output.as_deref(),
        &rows,
        &query_sets,
        &target_sets,
        &universe,
        OutputOptions {
            format: args.output_format,
            include_overlap_genes: args.overlap_genes,
            correction: args.correction,
        },
    )
}

#[derive(Serialize)]
struct MgsaSummary {
    query_id: String,
    query_size: u32,
    background_size: u32,
    modeled_terms: usize,
    samples: u64,
    burn: u64,
    acceptance_rate: f64,
    mean_active_terms: f64,
    term_rate: f64,
    false_pos_rate: f64,
    false_neg_rate: f64,
}

fn execute_mgsa(args: MgsaArgs) -> Result<()> {
    if args.output_format == OutputFormat::Parquet {
        bail!("mgsa supports --output-format tsv or null");
    }
    let rates = [args.flip_rate, args.step_rate, args.jump_rate];
    if rates.iter().any(|rate| *rate < 0.0) || rates.iter().all(|rate| *rate == 0.0) {
        bail!("move rates must be non-negative and not all zero");
    }
    for count in [args.term_count, args.false_pos_count, args.false_neg_count] {
        if count < 0.0 {
            bail!("prior pseudocounts must be non-negative");
        }
    }

    let target_source = TargetSource::load(&args.targets)?;
    let query_sets = read_sets(&args.sample, args.sample_format, &args.sample_name)
        .with_context(|| format!("failed to read sample sets from {}", args.sample.display()))?;
    let query_sets = filter_raw_sets(query_sets, args.sample_set.as_deref(), "sample")?;
    if query_sets.is_empty() {
        bail!("no sample set was loaded");
    }

    let background_genes = read_optional_gene_list(args.background.as_deref())?;
    let gene_names = read_optional_name_map(args.gene_names.as_deref())?;

    let mut universe = GeneUniverse::new();
    target_source.add_genes_to_universe(&mut universe);
    for query_set in &query_sets {
        universe.add_genes(&query_set.genes);
    }
    if let Some(background_genes) = &background_genes {
        universe.add_genes(background_genes);
    }
    universe.apply_gene_names(&gene_names);

    let target_sets = target_source.build_sets(&universe);
    ensure_non_empty_targets(&target_sets)?;
    let graph = target_source.build_graph(&target_sets);
    let queries = build_flat_sets(&query_sets, &universe);
    let background = build_background(background_genes.as_deref(), &target_sets, &universe);

    let options = MgsaOptions {
        term_prior: BetaPrior {
            mode: args.term_prob,
            count: args.term_count,
        },
        false_pos_prior: BetaPrior {
            mode: args.false_pos_prob,
            count: args.false_pos_count,
        },
        false_neg_prior: BetaPrior {
            mode: args.false_neg_prob,
            count: args.false_neg_count,
        },
        samples_per_term: args.samples_per_term,
        burn_per_term: args.burn_per_term,
        flip_rate: args.flip_rate,
        step_rate: args.step_rate,
        jump_rate: args.jump_rate,
        seed: args.seed,
        include_overlap_genes: args.overlap_genes,
    };
    let results = mgsa_matrix(
        &queries,
        &target_sets,
        &background,
        graph.as_ref(),
        &options,
    );

    if let Some(path) = &args.summary_output {
        let summaries: Vec<MgsaSummary> = results
            .iter()
            .map(|result| mgsa_summary(result, &queries))
            .collect();
        let file = File::create(path)
            .with_context(|| format!("failed to create summary {}", path.display()))?;
        serde_yaml::to_writer(file, &summaries)
            .with_context(|| format!("failed to write summary {}", path.display()))?;
    }
    if args.output_format == OutputFormat::Null {
        return Ok(());
    }
    write_output(args.output.as_deref(), |writer| {
        write_mgsa_rows(
            writer,
            &results,
            &queries,
            &target_sets,
            &universe,
            args.min_posterior,
            args.overlap_genes,
        )
    })
}

fn mgsa_summary(result: &MgsaResult, queries: &[NamedSet]) -> MgsaSummary {
    let stats = &result.stats;
    MgsaSummary {
        query_id: queries[result.query_index].id.clone(),
        query_size: result.query_size,
        background_size: result.background_size,
        modeled_terms: stats.modeled_terms,
        samples: stats.samples,
        burn: stats.burn,
        acceptance_rate: if stats.samples > 0 {
            stats.accepted_moves as f64 / stats.samples as f64
        } else {
            0.0
        },
        mean_active_terms: stats.mean_active_terms,
        term_rate: stats.term_rate,
        false_pos_rate: stats.false_pos_rate,
        false_neg_rate: stats.false_neg_rate,
    }
}

fn execute_gsea(args: GseaArgs) -> Result<()> {
    if args.output_format == OutputFormat::Parquet {
        bail!("gsea supports --output-format tsv or null");
    }
    if !args.weight.is_finite() || args.weight < 0.0 {
        bail!("--weight must be finite and non-negative");
    }
    if args.permutations == 0 {
        bail!("--permutations must be positive");
    }
    if args.min_size > args.max_size {
        bail!("--min-size must not exceed --max-size");
    }

    let target_source = TargetSource::load(&args.targets)?;
    let ranks_name = args
        .ranks_name
        .clone()
        .unwrap_or_else(|| default_set_name(&args.ranks));
    let mut ranked_lists = read_rankings(&args.ranks, &ranks_name)
        .with_context(|| format!("failed to read rankings from {}", args.ranks.display()))?;
    if let Some(wanted) = args.ranks_set.as_deref() {
        ranked_lists.retain(|list| list.id == wanted);
        if ranked_lists.is_empty() {
            bail!("no ranking named {wanted:?} was found");
        }
    }
    let gene_names = read_optional_name_map(args.gene_names.as_deref())?;

    let mut universe = GeneUniverse::new();
    target_source.add_genes_to_universe(&mut universe);
    for list in &ranked_lists {
        universe.add_genes(&list.genes);
    }
    universe.apply_gene_names(&gene_names);

    let target_sets = target_source.build_sets(&universe);
    ensure_non_empty_targets(&target_sets)?;
    let graph = target_source.build_graph(&target_sets);
    if args.algorithm != Algorithm::Classic && graph.is_none() {
        bail!("--algorithm elim and weight01 need --annotations with --closure");
    }
    let rankings = ranked_lists
        .iter()
        .map(|list| Ranking::new(list, &universe))
        .collect::<Result<Vec<_>>>()?;

    let options = GseaOptions {
        statistic: args.statistic,
        algorithm: args.algorithm,
        elim_cutoff: args.elim_cutoff,
        weight: args.weight,
        permutations: args.permutations,
        min_size: args.min_size,
        max_size: args.max_size,
        seed: args.seed,
        include_leading_edge: args.leading_edge_genes,
        max_p_value: args.max_p_value,
        max_p_adjust: args.max_p_adjust,
    };
    let rows = gsea_matrix(&rankings, &target_sets, graph.as_ref(), &options)?;
    if args.output_format == OutputFormat::Null {
        return Ok(());
    }
    write_output(args.output.as_deref(), |writer| {
        write_gsea_rows(
            writer,
            &rows,
            &rankings,
            &target_sets,
            &universe,
            args.leading_edge_genes,
        )
    })
}

fn execute_compare(args: CompareArgs) -> Result<()> {
    let left_format = resolve_format(&args.left, args.left_format)?;
    let right_format = resolve_format(&args.right, args.right_format)?;
    let left_records = load_result_records(&args.left, left_format)
        .with_context(|| format!("failed to load left results from {}", args.left.display()))?;
    let right_records = load_result_records(&args.right, right_format)
        .with_context(|| format!("failed to load right results from {}", args.right.display()))?;
    let options = DiffOptions {
        p_adjust_cutoff: args.p_adjust_cutoff,
        crossings_only: args.crossings_only,
    };
    let rows = compare_records(&left_records, &right_records, &options);
    let summary = summarize_diff(
        DiffSummaryContext {
            left_path: args.left.clone(),
            right_path: args.right.clone(),
            left_format,
            right_format,
            left_rows: left_records.len(),
            right_rows: right_records.len(),
        },
        &options,
        &rows,
    );
    if let Some(path) = &args.metadata_output {
        let file = File::create(path)
            .with_context(|| format!("failed to create metadata {}", path.display()))?;
        serde_yaml::to_writer(file, &summary)
            .with_context(|| format!("failed to write metadata {}", path.display()))?;
    }
    write_diff_result_rows(args.output.as_deref(), &rows, args.output_format)
}

fn read_optional_name_pairs(path: Option<&Path>) -> Result<Vec<(String, String)>> {
    path.map(read_name_pairs)
        .transpose()
        .map(Option::unwrap_or_default)
}

fn read_optional_name_map(path: Option<&Path>) -> Result<HashMap<String, String>> {
    path.map(read_name_map)
        .transpose()
        .map(Option::unwrap_or_default)
}

fn read_optional_closure(path: Option<&Path>) -> Result<Vec<(String, String)>> {
    path.map(read_closure_pairs)
        .transpose()
        .map(Option::unwrap_or_default)
}

fn read_optional_gene_list(path: Option<&Path>) -> Result<Option<Vec<String>>> {
    path.map(read_gene_list).transpose()
}

fn filter_raw_sets(sets: Vec<RawSet>, wanted: Option<&str>, label: &str) -> Result<Vec<RawSet>> {
    let Some(wanted) = wanted else {
        return Ok(sets);
    };
    let filtered: Vec<RawSet> = sets.into_iter().filter(|set| set.id == wanted).collect();
    if filtered.is_empty() {
        bail!("no {label} set named {wanted:?} was found");
    }
    Ok(filtered)
}

fn require_single_set(sets: Vec<RawSet>, label: &str, selector_arg: &str) -> Result<RawSet> {
    let mut sets = sets.into_iter();
    let Some(first) = sets.next() else {
        bail!("no {label} set was loaded");
    };
    if sets.next().is_some() {
        bail!("multiple {label} sets were loaded; choose one with {selector_arg}");
    }
    Ok(first)
}

fn ensure_non_empty_targets(target_sets: &[NamedSet]) -> Result<()> {
    if target_sets.is_empty() {
        bail!("no target sets were loaded");
    }
    Ok(())
}

fn build_background(
    explicit_background: Option<&[String]>,
    target_sets: &[NamedSet],
    universe: &GeneUniverse,
) -> genesets_rs::bitset::DenseBitSet {
    explicit_background.map_or_else(
        || union_sets(target_sets, universe.len()),
        |genes| build_gene_list_bitset(genes, universe),
    )
}

fn filter_named_sets(
    sets: Vec<NamedSet>,
    wanted: Option<&str>,
    label: &str,
) -> Result<Vec<NamedSet>> {
    let Some(wanted) = wanted else {
        return Ok(sets);
    };
    let filtered: Vec<NamedSet> = sets.into_iter().filter(|set| set.id == wanted).collect();
    if filtered.is_empty() {
        bail!("no {label} set named {wanted:?} was found");
    }
    Ok(filtered)
}

fn write_output<F>(path: Option<&Path>, write: F) -> Result<()>
where
    F: FnOnce(&mut dyn Write) -> Result<()>,
{
    if let Some(path) = path {
        let mut file = File::create(path)
            .with_context(|| format!("failed to create output {}", path.display()))?;
        write(&mut file)
    } else {
        let stdout = io::stdout();
        let mut handle = stdout.lock();
        write(&mut handle)
    }
}

fn write_result_rows(
    path: Option<&Path>,
    rows: &[genesets_rs::EnrichmentRow],
    queries: &[NamedSet],
    targets: &[NamedSet],
    universe: &GeneUniverse,
    options: OutputOptions,
) -> Result<()> {
    match options.format {
        OutputFormat::Parquet => {
            let path = path.context("parquet output requires --output")?;
            let file = File::create(path)
                .with_context(|| format!("failed to create output {}", path.display()))?;
            write_parquet_rows(file, rows, queries, targets, universe, options)
        }
        OutputFormat::Tsv | OutputFormat::Null => write_output(path, |writer| {
            write_rows(writer, rows, queries, targets, universe, options)
        }),
    }
}

fn write_diff_result_rows(
    path: Option<&Path>,
    rows: &[genesets_rs::compare::DiffRow],
    format: OutputFormat,
) -> Result<()> {
    match format {
        OutputFormat::Parquet => {
            let path = path.context("parquet output requires --output")?;
            let file = File::create(path)
                .with_context(|| format!("failed to create output {}", path.display()))?;
            write_parquet_diff_rows(file, rows)
        }
        OutputFormat::Tsv => write_output(path, |writer| write_diff_rows(writer, rows)),
        OutputFormat::Null => Ok(()),
    }
}

fn run_with_thread_pool<F, T>(threads: Option<usize>, run: F) -> Result<T>
where
    F: FnOnce() -> Result<T> + Send,
    T: Send,
{
    if let Some(threads) = threads {
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .context("failed to initialize Rayon thread pool")?;
        pool.install(run)
    } else {
        run()
    }
}

fn parse_probability(value: &str) -> Result<f64, String> {
    let parsed = value
        .parse::<f64>()
        .map_err(|error| format!("expected a number between 0 and 1: {error}"))?;
    if (0.0..=1.0).contains(&parsed) {
        Ok(parsed)
    } else {
        Err("expected a number between 0 and 1".to_owned())
    }
}
