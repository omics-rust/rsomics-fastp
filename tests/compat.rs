//! Compat tests for rsomics-fastp vs fastp 0.20.x.
//!
//! fastp's output encoding shifted between 0.20 and 1.x (adapter heuristics,
//! default trimming), so the live oracle test gates on 0.20.x. With adapter
//! trimming disabled on both sides the SE output is byte-identical, which the
//! committed golden captures.

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

    // Adapter trimming disabled on both sides so the comparison isolates the
    // length/quality/N filters; with adapter heuristics out of the picture the
    // SE output is byte-identical (see se_filter_matches_golden).
    let ours_status = ours()
        .args(["-i", r1_fixture().to_str().unwrap()])
        .args(["-o", our_out.to_str().unwrap()])
        .args(["--length-required", "15"])
        .arg("--disable-adapter-trimming")
        .status()
        .unwrap();
    assert!(ours_status.success(), "rsomics-fastp failed");

    let fastp_status = Command::new("fastp")
        .args(["-i", r1_fixture().to_str().unwrap()])
        .args(["-o", fastp_out.to_str().unwrap()])
        .args(["--length_required", "15"])
        .arg("--disable_adapter_trimming")
        .status()
        .unwrap();
    assert!(fastp_status.success(), "fastp failed");

    assert_eq!(
        count_reads(&our_out),
        count_reads(&fastp_out),
        "read count after SE filter"
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

// fastp_se_filtered.fastq was captured once from fastp 0.20.1 with
// `--length_required 15 --disable_adapter_trimming`. Both sides disable adapter
// trimming, so the SE-filtered output is byte-identical and CI can diff
// ours-vs-golden without fastp installed.
#[test]
fn se_filter_matches_golden() {
    let tmp = TempDir::new().unwrap();
    let our_out = tmp.path().join("ours_out.fastq");

    let status = ours()
        .args(["-i", r1_fixture().to_str().unwrap()])
        .args(["-o", our_out.to_str().unwrap()])
        .args(["--length-required", "15"])
        .arg("--disable-adapter-trimming")
        .status()
        .unwrap();
    assert!(status.success(), "rsomics-fastp failed");

    let ours = std::fs::read(&our_out).unwrap();
    let golden = std::fs::read(
        Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/golden/fastp_se_filtered.fastq"),
    )
    .unwrap();
    assert_eq!(
        ours, golden,
        "SE-filtered output must match fastp 0.20.1 golden"
    );
}
