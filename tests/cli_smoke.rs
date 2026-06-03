use std::{fs, path::Path, process::Command};

fn bin() -> &'static str {
    env!("CARGO_BIN_EXE_genesets-rs")
}

fn run(args: &[&str]) -> String {
    let output = Command::new(bin())
        .args(args)
        .output()
        .expect("failed to run genesets-rs");
    assert!(
        output.status.success(),
        "command failed\nstatus: {:?}\nstderr:\n{}\nstdout:\n{}",
        output.status.code(),
        String::from_utf8_lossy(&output.stderr),
        String::from_utf8_lossy(&output.stdout)
    );
    String::from_utf8(output.stdout).expect("stdout was not utf8")
}

fn run_from_dir(args: &[&str], cwd: &Path) -> String {
    let output = Command::new(bin())
        .args(args)
        .current_dir(cwd)
        .output()
        .expect("failed to run genesets-rs");
    assert!(
        output.status.success(),
        "command failed\nstatus: {:?}\nstderr:\n{}\nstdout:\n{}",
        output.status.code(),
        String::from_utf8_lossy(&output.stderr),
        String::from_utf8_lossy(&output.stdout)
    );
    String::from_utf8(output.stdout).expect("stdout was not utf8")
}

#[test]
fn enrich_example_outputs_expected_terms() {
    let stdout = run(&[
        "enrich",
        "--annotations",
        "examples/gene_terms.tsv",
        "--terms",
        "examples/terms.tsv",
        "--closure",
        "examples/closure.tsv",
        "--sample",
        "examples/sample.txt",
        "--background",
        "examples/background.txt",
        "--overlap-genes",
    ]);

    let lines: Vec<&str> = stdout.lines().collect();
    assert_eq!(lines.len(), 3);
    assert!(lines[0].starts_with("query_id\tquery_name\ttarget_id"));
    assert!(lines[1].contains("sample\t\tT:child\tChild process\t2\t3\t2\t5"));
    assert!(lines[1].contains("g1;g2"));
    assert!(lines[2].contains("sample\t\tT:parent\tParent process\t2\t3\t3\t5"));
}

#[test]
fn matrix_queries_from_targets_outputs_term_pairs() {
    let stdout = run(&[
        "matrix",
        "--annotations",
        "examples/gene_terms.tsv",
        "--terms",
        "examples/terms.tsv",
        "--closure",
        "examples/closure.tsv",
        "--queries-from-targets",
        "--background",
        "examples/background.txt",
    ]);

    let lines: Vec<&str> = stdout.lines().collect();
    assert_eq!(lines.len(), 6);
    assert!(stdout.contains("T:child\tChild process\tT:child\tChild process\t2\t2\t2\t5"));
    assert!(stdout.contains("T:parent\tParent process\tT:parent\tParent process\t3\t3\t3\t5"));
    assert!(stdout.contains("T:other\tOther process\tT:other\tOther process\t1\t1\t1\t5"));
}

#[test]
fn yaml_run_matches_enrich_command() {
    let stdout = run(&["run", "examples/enrich.yaml"]);

    assert!(stdout.contains("query_id\tquery_name\ttarget_id"));
    assert!(stdout.contains("sample\t\tT:child\tChild process\t2\t3\t2\t5"));
    assert!(stdout.contains("sample\t\tT:parent\tParent process\t2\t3\t3\t5"));
}

#[test]
fn yaml_paths_are_relative_to_config_file() {
    let config = Path::new(env!("CARGO_MANIFEST_DIR")).join("examples/enrich.yaml");
    let config_arg = config.to_str().expect("example config path was not utf8");
    let stdout = run_from_dir(&["run", config_arg], &std::env::temp_dir());

    assert!(stdout.contains("query_id\tquery_name\ttarget_id"));
    assert!(stdout.contains("sample\t\tT:child\tChild process\t2\t3\t2\t5"));
    assert!(stdout.contains("sample\t\tT:parent\tParent process\t2\t3\t3\t5"));
}

#[test]
fn enrich_can_filter_by_adjusted_p_value() {
    let stdout = run(&[
        "enrich",
        "--annotations",
        "examples/gene_terms.tsv",
        "--terms",
        "examples/terms.tsv",
        "--closure",
        "examples/closure.tsv",
        "--sample",
        "examples/sample.txt",
        "--background",
        "examples/background.txt",
        "--max-p-adjust",
        "0.95",
    ]);

    let lines: Vec<&str> = stdout.lines().collect();
    assert_eq!(lines.len(), 2);
    assert!(lines[1].contains("sample\t\tT:child\tChild process\t2\t3\t2\t5"));
}

#[test]
fn matrix_can_discard_output() {
    let stdout = run(&[
        "matrix",
        "--annotations",
        "examples/gene_terms.tsv",
        "--terms",
        "examples/terms.tsv",
        "--closure",
        "examples/closure.tsv",
        "--queries-from-targets",
        "--background",
        "examples/background.txt",
        "--output-format",
        "null",
    ]);

    assert!(stdout.is_empty());
}

#[test]
fn matrix_can_write_parquet_output() {
    let path = std::env::temp_dir().join(format!(
        "genesets-rs-{}-matrix-smoke.parquet",
        std::process::id()
    ));
    let _ = fs::remove_file(&path);

    let output = Command::new(bin())
        .args([
            "matrix",
            "--annotations",
            "examples/gene_terms.tsv",
            "--terms",
            "examples/terms.tsv",
            "--closure",
            "examples/closure.tsv",
            "--queries-from-targets",
            "--background",
            "examples/background.txt",
            "--output-format",
            "parquet",
            "--output",
        ])
        .arg(&path)
        .output()
        .expect("failed to run genesets-rs");
    assert!(
        output.status.success(),
        "command failed\nstatus: {:?}\nstderr:\n{}\nstdout:\n{}",
        output.status.code(),
        String::from_utf8_lossy(&output.stderr),
        String::from_utf8_lossy(&output.stdout)
    );
    assert!(output.stdout.is_empty());

    let bytes = fs::read(&path).expect("parquet output was not written");
    assert!(bytes.len() > 8);
    assert_eq!(&bytes[..4], b"PAR1");
    assert_eq!(&bytes[bytes.len() - 4..], b"PAR1");

    if command_exists("duckdb") {
        let sql = format!(
            "SELECT COUNT(*) FROM read_parquet('{}');",
            sql_string_literal(&path)
        );
        let duckdb = Command::new("duckdb")
            .args(["-csv", "-noheader", "-c", &sql])
            .output()
            .expect("failed to run duckdb");
        assert!(
            duckdb.status.success(),
            "duckdb failed\nstderr:\n{}\nstdout:\n{}",
            String::from_utf8_lossy(&duckdb.stderr),
            String::from_utf8_lossy(&duckdb.stdout)
        );
        assert_eq!(String::from_utf8_lossy(&duckdb.stdout).trim(), "5");
    }

    let _ = fs::remove_file(&path);
}

#[test]
fn compare_reports_threshold_crossings() {
    let prefix = std::env::temp_dir().join(format!("genesets-rs-{}", std::process::id()));
    let left = prefix.with_extension("left.tsv");
    let right = prefix.with_extension("right.tsv");
    let metadata = prefix.with_extension("compare.yaml");
    let _ = fs::remove_file(&left);
    let _ = fs::remove_file(&right);
    let _ = fs::remove_file(&metadata);

    fs::write(
        &left,
        concat!(
            "query_id\tquery_name\ttarget_id\ttarget_name\toverlap\tquery_size\ttarget_size\tbackground_size\tp_value\tp_adjust_bonferroni\n",
            "q1\tQuery 1\tt_lost\tLost term\t5\t10\t20\t100\t1e-6\t0.01\n",
            "q1\tQuery 1\tt_shared\tShared term\t6\t10\t30\t100\t1e-7\t0.02\n",
            "q2\tQuery 2\tt_gain\tGained term\t1\t8\t18\t100\t0.5\t0.9\n",
        ),
    )
    .unwrap();
    fs::write(
        &right,
        concat!(
            "query_id\tquery_name\ttarget_id\ttarget_name\toverlap\tquery_size\ttarget_size\tbackground_size\tp_value\tp_adjust_bonferroni\n",
            "q1\tQuery 1\tt_lost\tLost term\t2\t10\t19\t100\t0.2\t0.7\n",
            "q1\tQuery 1\tt_shared\tShared term\t7\t10\t33\t100\t1e-8\t0.001\n",
            "q2\tQuery 2\tt_gain\tGained term\t6\t8\t21\t100\t1e-5\t0.03\n",
        ),
    )
    .unwrap();

    let output = Command::new(bin())
        .args([
            "compare",
            "--left",
            left.to_str().unwrap(),
            "--right",
            right.to_str().unwrap(),
            "--p-adjust-cutoff",
            "0.05",
            "--metadata-output",
            metadata.to_str().unwrap(),
        ])
        .output()
        .expect("failed to run genesets-rs");
    assert!(
        output.status.success(),
        "command failed\nstatus: {:?}\nstderr:\n{}\nstdout:\n{}",
        output.status.code(),
        String::from_utf8_lossy(&output.stderr),
        String::from_utf8_lossy(&output.stdout)
    );
    let stdout = String::from_utf8(output.stdout).unwrap();
    assert!(stdout.contains("lost_significant\tq1\tt_lost"));
    assert!(stdout.contains("gained_significant\tq2\tt_gain"));
    assert!(stdout.contains("shared_significant\tq1\tt_shared"));

    let metadata = fs::read_to_string(&metadata).unwrap();
    assert!(metadata.contains("lost_significant: 1"));
    assert!(metadata.contains("gained_significant: 1"));
    assert!(metadata.contains("shared_significant: 1"));

    let _ = fs::remove_file(&left);
    let _ = fs::remove_file(&right);
    let _ = fs::remove_file(prefix.with_extension("compare.yaml"));
}

#[test]
fn compare_can_read_and_write_parquet() {
    let input = std::env::temp_dir().join(format!(
        "genesets-rs-{}-compare-input.parquet",
        std::process::id()
    ));
    let output_path = std::env::temp_dir().join(format!(
        "genesets-rs-{}-compare-output.parquet",
        std::process::id()
    ));
    let _ = fs::remove_file(&input);
    let _ = fs::remove_file(&output_path);

    let matrix = Command::new(bin())
        .args([
            "matrix",
            "--annotations",
            "examples/gene_terms.tsv",
            "--terms",
            "examples/terms.tsv",
            "--closure",
            "examples/closure.tsv",
            "--queries-from-targets",
            "--background",
            "examples/background.txt",
            "--output-format",
            "parquet",
            "--output",
        ])
        .arg(&input)
        .output()
        .expect("failed to run matrix");
    assert!(
        matrix.status.success(),
        "matrix failed\nstderr:\n{}\nstdout:\n{}",
        String::from_utf8_lossy(&matrix.stderr),
        String::from_utf8_lossy(&matrix.stdout)
    );

    let compare = Command::new(bin())
        .args([
            "compare",
            "--left",
            input.to_str().unwrap(),
            "--right",
            input.to_str().unwrap(),
            "--p-adjust-cutoff",
            "0.95",
            "--output-format",
            "parquet",
            "--output",
        ])
        .arg(&output_path)
        .output()
        .expect("failed to run compare");
    assert!(
        compare.status.success(),
        "compare failed\nstderr:\n{}\nstdout:\n{}",
        String::from_utf8_lossy(&compare.stderr),
        String::from_utf8_lossy(&compare.stdout)
    );
    let bytes = fs::read(&output_path).expect("parquet diff output was not written");
    assert_eq!(&bytes[..4], b"PAR1");
    assert_eq!(&bytes[bytes.len() - 4..], b"PAR1");

    if command_exists("duckdb") {
        let sql = format!(
            "SELECT COUNT(*) FROM read_parquet('{}');",
            sql_string_literal(&output_path)
        );
        let duckdb = Command::new("duckdb")
            .args(["-csv", "-noheader", "-c", &sql])
            .output()
            .expect("failed to run duckdb");
        assert!(
            duckdb.status.success(),
            "duckdb failed\nstderr:\n{}\nstdout:\n{}",
            String::from_utf8_lossy(&duckdb.stderr),
            String::from_utf8_lossy(&duckdb.stdout)
        );
        assert_eq!(String::from_utf8_lossy(&duckdb.stdout).trim(), "2");
    }

    let _ = fs::remove_file(&input);
    let _ = fs::remove_file(&output_path);
}

fn command_exists(command: &str) -> bool {
    Command::new(command)
        .arg("-version")
        .output()
        .is_ok_and(|output| output.status.success())
}

fn sql_string_literal(path: &Path) -> String {
    path.to_string_lossy().replace('\'', "''")
}
