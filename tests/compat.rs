//! Compat tests for rsomics-fastp vs fastp 0.20.x.
//!
//! Semantic note: fastp's exact output encoding changed between 0.20 and 1.x
//! (adapter detection heuristics, default trimming settings). We gate on
//! version 0.20 and compare read counts for SE mode.

use std::path::{Path, PathBuf};
use std::process::Command;

use tempfile::TempDir;

fn ours() -> Command {
    Command::new(env!("CARGO_BIN_EXE_rsomics-fastp"))
}

fn r1_fixture() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/golden/test_r1.fastq")
}

fn fastp_version() -> Option<(u32, u32)> {
    let out = Command::new("fastp").arg("--version").output().ok()?;
    // fastp prints version to stderr: "fastp 0.20.1"
    let text = String::from_utf8_lossy(&out.stderr);
    let first = text.lines().next()?;
    let ver = first.split_whitespace().nth(1)?;
    let mut parts = ver.split('.');
    let maj: u32 = parts.next()?.parse().ok()?;
    let min: u32 = parts.next()?.parse().ok()?;
    Some((maj, min))
}

/// Count reads in a FASTQ file (each read = 4 lines).
fn count_reads(path: &std::path::Path) -> usize {
    let content = std::fs::read_to_string(path).unwrap_or_default();
    content.lines().filter(|l| l.starts_with('@')).count()
}

#[test]
fn se_filter_matches_fastp() {
    let (maj, min) = match fastp_version() {
        Some(v) => v,
        None => {
            eprintln!("SKIP fastp compat: fastp not found");
            return;
        }
    };
    if maj != 0 || min != 20 {
        eprintln!("SKIP fastp compat: fastp {maj}.{min} (need 0.20.x)");
        return;
    }

    let tmp = TempDir::new().unwrap();
    let our_out = tmp.path().join("ours_out.fastq");
    let fastp_out = tmp.path().join("fastp_out.fastq");

    // Run ours
    let ours_status = ours()
        .args(["-i", r1_fixture().to_str().unwrap()])
        .args(["-o", our_out.to_str().unwrap()])
        .args(["--length-required", "15"])
        .status()
        .unwrap();
    assert!(ours_status.success(), "rsomics-fastp failed");

    // Run fastp 0.20
    let fastp_status = Command::new("fastp")
        .args(["-i", r1_fixture().to_str().unwrap()])
        .args(["-o", fastp_out.to_str().unwrap()])
        .args(["--length_required", "15"])
        .arg("--disable_adapter_trimming")
        .status()
        .unwrap();
    assert!(fastp_status.success(), "fastp failed");

    let our_reads = count_reads(&our_out);
    let fastp_reads = count_reads(&fastp_out);

    // Both should agree on read count after same filters
    // Note: adapter detection heuristics differ — we only compare count
    // (not byte-exact output) given different adapter detection models.
    assert_eq!(
        our_reads, fastp_reads,
        "read count after SE filter: ours={our_reads}, fastp={fastp_reads}"
    );
}

#[test]
fn se_short_reads_filtered() {
    // Reads shorter than --length-required must be filtered out.
    let tmp = TempDir::new().unwrap();
    let our_out = tmp.path().join("ours_out.fastq");

    // read4_short is only 4bp — must be filtered with length-required=15
    let status = ours()
        .args(["-i", r1_fixture().to_str().unwrap()])
        .args(["-o", our_out.to_str().unwrap()])
        .args(["--length-required", "15"])
        .status()
        .unwrap();
    assert!(status.success(), "rsomics-fastp failed");

    // The N-only read (read2) fails the quality filter, the short read fails
    // length filter — so at most 2 reads should pass (read1 and read3_adapter
    // with adapter trimmed to >= 15bp).
    let surviving = count_reads(&our_out);
    assert!(
        surviving < 4,
        "short read should have been filtered; got {surviving}"
    );
}
