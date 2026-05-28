use criterion::{Criterion, criterion_group, criterion_main};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::PathBuf;
use std::process::Command;

const N: usize = 200_000;
const READ_LEN: usize = 150;
const SEED: u64 = 0x00FA_5701;

fn ensure_fixture() -> PathBuf {
    let mut p = std::env::temp_dir();
    p.push(format!("rsomics-fastp-bench-{N}.fq"));
    if p.exists() {
        return p;
    }
    let mut w = BufWriter::new(File::create(&p).unwrap());
    let mut rng = SEED;
    let qual = [b'I'; READ_LEN];
    for i in 0..N {
        let mut seq = Vec::with_capacity(READ_LEN);
        for _ in 0..READ_LEN {
            rng = rng.wrapping_mul(6_364_136_223_846_793_005).wrapping_add(1);
            seq.push(b"ACGT"[((rng >> 33) & 3) as usize]);
        }
        writeln!(w, "@read{i}").unwrap();
        w.write_all(&seq).unwrap();
        w.write_all(b"\n+\n").unwrap();
        w.write_all(&qual).unwrap();
        w.write_all(b"\n").unwrap();
    }
    p
}

fn bench(c: &mut Criterion) {
    let fq = ensure_fixture();
    let ours = env!("CARGO_BIN_EXE_rsomics-fastp");
    let mut group = c.benchmark_group(format!("fastp/{N}reads"));
    group.sample_size(10);

    group.bench_function("rsomics-fastp-se", |b| {
        b.iter(|| {
            let out = Command::new(ours)
                .args(["-i", fq.to_str().unwrap()])
                .args(["-o", "/dev/null"])
                .output()
                .expect("run rsomics-fastp");
            assert!(out.status.success());
        });
    });

    if Command::new("fastp").arg("--version").output().is_ok() {
        group.bench_function("fastp-upstream-se", |b| {
            b.iter(|| {
                let out = Command::new("fastp")
                    .args(["-i", fq.to_str().unwrap()])
                    .args(["-o", "/dev/null"])
                    .args(["-j", "/dev/null"])
                    .args(["--thread", "1"])
                    .output()
                    .expect("run fastp");
                assert!(out.status.success());
            });
        });
    }

    group.finish();
}

criterion_group!(benches, bench);
criterion_main!(benches);
