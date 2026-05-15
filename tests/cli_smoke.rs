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

fn command_exists(command: &str) -> bool {
    Command::new(command)
        .arg("-version")
        .output()
        .is_ok_and(|output| output.status.success())
}

fn sql_string_literal(path: &Path) -> String {
    path.to_string_lossy().replace('\'', "''")
}
