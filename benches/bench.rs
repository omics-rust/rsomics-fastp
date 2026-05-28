use criterion::{Criterion, criterion_group, criterion_main};
use std::hint::black_box;
use std::path::PathBuf;
use std::process::Command;

fn bench_fastp(c: &mut Criterion) {
    let bin = env!("CARGO_BIN_EXE_rsomics-fastp");
    let manifest = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let fq = manifest.join("tests/golden/test_r1.fastq");
    c.bench_function("rsomics-fastp golden", |b| {
        b.iter(|| {
            let out = Command::new(black_box(bin))
                .args(["--in1", fq.to_str().unwrap(), "--out1", "/dev/null"])
                .output()
                .unwrap();
            assert!(out.status.success());
        });
    });
}

criterion_group!(benches, bench_fastp);
criterion_main!(benches);
