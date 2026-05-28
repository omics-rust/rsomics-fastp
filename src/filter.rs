use crate::report::ReadStats;

#[derive(Debug, Clone)]
pub struct FilterConfig {
    pub min_qual: u8,
    pub max_low_qual_fraction: f64,
    pub max_n_bases: usize,
    pub min_length: usize,
    /// 0 means unlimited.
    pub max_length: usize,
}

impl Default for FilterConfig {
    fn default() -> Self {
        Self {
            min_qual: 15,
            max_low_qual_fraction: 0.4,
            max_n_bases: 5,
            min_length: 15,
            max_length: 0,
        }
    }
}

#[derive(Debug, PartialEq, Eq)]
pub enum FilterResult {
    Pass,
    TooShort,
    TooLong,
    TooManyLowQualBases,
    TooManyNBases,
}

pub fn filter_read(seq: &[u8], qual: &[u8], cfg: &FilterConfig) -> FilterResult {
    let len = seq.len();
    if len < cfg.min_length {
        return FilterResult::TooShort;
    }
    if cfg.max_length > 0 && len > cfg.max_length {
        return FilterResult::TooLong;
    }
    let mut n_count = 0usize;
    let mut low_qual_count = 0usize;
    for (&b, &q) in seq.iter().zip(qual.iter()) {
        if b == b'N' || b == b'n' {
            n_count += 1;
        }
        // FASTQ: qual bytes are Phred+33 ASCII values.
        if q < cfg.min_qual + 33 {
            low_qual_count += 1;
        }
    }
    if n_count > cfg.max_n_bases {
        return FilterResult::TooManyNBases;
    }
    let low_qual_frac = low_qual_count as f64 / len as f64;
    if low_qual_frac > cfg.max_low_qual_fraction {
        return FilterResult::TooManyLowQualBases;
    }
    FilterResult::Pass
}

pub fn update_stats(stats: &mut ReadStats, seq: &[u8], qual: &[u8]) {
    stats.total_reads += 1;
    stats.total_bases += seq.len() as u64;
    for &q in qual.iter() {
        let phred = q.saturating_sub(33);
        if phred >= 20 {
            stats.q20_bases += 1;
        }
        if phred >= 30 {
            stats.q30_bases += 1;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn hq(n: usize) -> Vec<u8> {
        vec![b'I'; n] // Phred 40 (well above threshold 15)
    }

    #[test]
    fn pass_high_qual_read() {
        let seq = b"ACGTACGTACGTACGT"; // 16 bases >= min_length 15
        let qual = hq(seq.len());
        assert_eq!(
            filter_read(seq, &qual, &FilterConfig::default()),
            FilterResult::Pass
        );
    }

    #[test]
    fn fail_too_short() {
        let seq = b"ACGTACGT"; // 8 bases < 15
        let qual = hq(seq.len());
        assert_eq!(
            filter_read(seq, &qual, &FilterConfig::default()),
            FilterResult::TooShort
        );
    }

    #[test]
    fn fail_too_long() {
        let seq = vec![b'A'; 100];
        let qual = hq(100);
        let cfg = FilterConfig {
            max_length: 50,
            ..FilterConfig::default()
        };
        assert_eq!(filter_read(&seq, &qual, &cfg), FilterResult::TooLong);
    }

    #[test]
    fn fail_low_qual() {
        // 7 of 16 bases below qual threshold => 43.75% > 40% limit
        let seq = b"ACGTACGTACGTACGT";
        let mut qual = hq(seq.len());
        for q in qual[..7].iter_mut() {
            *q = 33 + 10; // Phred 10, below threshold 15
        }
        assert_eq!(
            filter_read(seq, &qual, &FilterConfig::default()),
            FilterResult::TooManyLowQualBases
        );
    }

    #[test]
    fn fail_too_many_n() {
        // 6 N bases > max_n_bases 5
        let seq = b"NNNNNNACGTACGTAC"; // 16 bases, 6 N
        let qual = hq(seq.len());
        assert_eq!(
            filter_read(seq, &qual, &FilterConfig::default()),
            FilterResult::TooManyNBases
        );
    }
}
