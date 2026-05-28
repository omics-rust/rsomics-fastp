use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::Path;

use flate2::Compression;
use flate2::write::GzEncoder;
use needletail::parse_fastx_file;

use crate::adapter::{AdapterConfig, detect_pe_overlap, trim_adapter_se};
use crate::filter::{FilterConfig, FilterResult, filter_read, update_stats};
use crate::poly::trim_poly_x;
use crate::report::{FilteringResult, ReadStats, Report, Summary};

#[derive(Debug, Clone)]
pub struct RunConfig {
    pub filter: FilterConfig,
    pub adapter: AdapterConfig,
    pub adapter_trimming: bool,
    pub trim_poly_g: bool,
    pub trim_poly_x: bool,
    pub poly_min_len: usize,
}

impl Default for RunConfig {
    fn default() -> Self {
        Self {
            filter: FilterConfig::default(),
            adapter: AdapterConfig::default(),
            adapter_trimming: true,
            trim_poly_g: false,
            trim_poly_x: false,
            poly_min_len: 10,
        }
    }
}

pub struct RunResult {
    pub before: ReadStats,
    pub after: ReadStats,
    pub filtering: FilteringResult,
}

fn open_writer(path: &Path) -> anyhow::Result<Box<dyn Write>> {
    let file = File::create(path)?;
    if path.extension().and_then(|e| e.to_str()) == Some("gz") {
        Ok(Box::new(GzEncoder::new(file, Compression::default())))
    } else {
        Ok(Box::new(BufWriter::new(file)))
    }
}

fn write_record(w: &mut dyn Write, name: &[u8], seq: &[u8], qual: &[u8]) -> io::Result<()> {
    w.write_all(b"@")?;
    w.write_all(name)?;
    w.write_all(b"\n")?;
    w.write_all(seq)?;
    w.write_all(b"\n+\n")?;
    w.write_all(qual)?;
    w.write_all(b"\n")
}

fn apply_poly_trim(seq: &[u8], len: usize, cfg: &RunConfig) -> usize {
    let mut t = len;
    if cfg.trim_poly_g {
        t = trim_poly_x(&seq[..t], b'G', cfg.poly_min_len);
    }
    if cfg.trim_poly_x {
        for base in [b'A', b'T', b'C', b'G'] {
            let u = trim_poly_x(&seq[..t], base, cfg.poly_min_len);
            if u < t {
                t = u;
            }
        }
    }
    t
}

pub fn run_se(in_path: &Path, out_path: &Path, cfg: &RunConfig) -> anyhow::Result<RunResult> {
    let mut writer = open_writer(out_path)?;
    let mut before = ReadStats::default();
    let mut after = ReadStats::default();
    let mut filtering = FilteringResult::default();

    let mut reader = parse_fastx_file(in_path)?;
    while let Some(record) = reader.next() {
        let rec = record?;
        let name = rec.id();
        let seq_cow = rec.seq();
        let raw_seq: &[u8] = &seq_cow;
        let raw_qual = rec.qual().unwrap_or(&[]);

        update_stats(&mut before, raw_seq, raw_qual);
        filtering.passed_filter_reads += 1;

        let mut t = raw_seq.len();
        t = apply_poly_trim(raw_seq, t, cfg);

        if cfg.adapter_trimming {
            t = trim_adapter_se(&raw_seq[..t], &cfg.adapter.adapter_r1, &cfg.adapter);
        }

        let seq = &raw_seq[..t];
        let qual = if raw_qual.is_empty() {
            raw_qual
        } else {
            &raw_qual[..t]
        };

        match filter_read(seq, qual, &cfg.filter) {
            FilterResult::Pass => {
                update_stats(&mut after, seq, qual);
                write_record(&mut *writer, name, seq, qual)?;
            }
            FilterResult::TooShort => {
                filtering.passed_filter_reads -= 1;
                filtering.too_short_reads += 1;
            }
            FilterResult::TooLong => {
                filtering.passed_filter_reads -= 1;
                filtering.too_long_reads += 1;
            }
            FilterResult::TooManyLowQualBases => {
                filtering.passed_filter_reads -= 1;
                filtering.low_quality_reads += 1;
            }
            #[allow(non_snake_case)]
            FilterResult::TooManyNBases => {
                filtering.passed_filter_reads -= 1;
                filtering.too_many_N_reads += 1;
            }
        }
    }

    Ok(RunResult {
        before,
        after,
        filtering,
    })
}

pub fn run_pe(
    in1: &Path,
    in2: &Path,
    out1: &Path,
    out2: &Path,
    cfg: &RunConfig,
) -> anyhow::Result<RunResult> {
    let mut w1 = open_writer(out1)?;
    let mut w2 = open_writer(out2)?;
    let mut before = ReadStats::default();
    let mut after = ReadStats::default();
    let mut filtering = FilteringResult::default();

    let mut r1_reader = parse_fastx_file(in1)?;
    let mut r2_reader = parse_fastx_file(in2)?;

    loop {
        let rec1 = r1_reader.next();
        let rec2 = r2_reader.next();
        match (rec1, rec2) {
            (None, None) => break,
            (Some(r1), Some(r2)) => {
                let r1 = r1?;
                let r2 = r2?;
                let name1 = r1.id();
                let name2 = r2.id();
                let seq1_cow = r1.seq();
                let seq2_cow = r2.seq();
                let raw1: &[u8] = &seq1_cow;
                let raw2: &[u8] = &seq2_cow;
                let q1 = r1.qual().unwrap_or(&[]);
                let q2 = r2.qual().unwrap_or(&[]);

                update_stats(&mut before, raw1, q1);
                update_stats(&mut before, raw2, q2);
                filtering.passed_filter_reads += 2;

                let mut t1 = raw1.len();
                let mut t2 = raw2.len();

                t1 = apply_poly_trim(raw1, t1, cfg);
                t2 = apply_poly_trim(raw2, t2, cfg);

                if cfg.adapter_trimming {
                    if let Some((s1, s2)) =
                        detect_pe_overlap(&raw1[..t1], &raw2[..t2], &cfg.adapter)
                    {
                        t1 = s1;
                        t2 = s2;
                    } else {
                        t1 = trim_adapter_se(&raw1[..t1], &cfg.adapter.adapter_r1, &cfg.adapter);
                        t2 = trim_adapter_se(&raw2[..t2], &cfg.adapter.adapter_r2, &cfg.adapter);
                    }
                }

                let seq1 = &raw1[..t1];
                let seq2 = &raw2[..t2];
                let qual1 = if q1.is_empty() { q1 } else { &q1[..t1] };
                let qual2 = if q2.is_empty() { q2 } else { &q2[..t2] };

                let res1 = filter_read(seq1, qual1, &cfg.filter);
                let res2 = filter_read(seq2, qual2, &cfg.filter);

                if res1 == FilterResult::Pass && res2 == FilterResult::Pass {
                    update_stats(&mut after, seq1, qual1);
                    update_stats(&mut after, seq2, qual2);
                    write_record(&mut *w1, name1, seq1, qual1)?;
                    write_record(&mut *w2, name2, seq2, qual2)?;
                } else {
                    filtering.passed_filter_reads -= 2;
                    for res in [&res1, &res2] {
                        match res {
                            FilterResult::TooShort => filtering.too_short_reads += 1,
                            FilterResult::TooLong => filtering.too_long_reads += 1,
                            FilterResult::TooManyLowQualBases => filtering.low_quality_reads += 1,
                            #[allow(non_snake_case)]
                            FilterResult::TooManyNBases => filtering.too_many_N_reads += 1,
                            FilterResult::Pass => {}
                        }
                    }
                }
            }
            _ => anyhow::bail!("paired FASTQ files have different record counts"),
        }
    }

    Ok(RunResult {
        before,
        after,
        filtering,
    })
}

pub fn write_report(result: &RunResult, json_path: &Path) -> anyhow::Result<()> {
    let report = Report {
        summary: Summary {
            before_filtering: result.before.clone(),
            after_filtering: result.after.clone(),
        },
        #[allow(non_snake_case)]
        filtering_result: FilteringResult {
            passed_filter_reads: result.filtering.passed_filter_reads,
            low_quality_reads: result.filtering.low_quality_reads,
            too_many_N_reads: result.filtering.too_many_N_reads,
            too_short_reads: result.filtering.too_short_reads,
            too_long_reads: result.filtering.too_long_reads,
        },
    };
    let json = serde_json::to_string_pretty(&report)?;
    std::fs::write(json_path, json)?;
    Ok(())
}
