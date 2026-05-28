/// Illumina TruSeq Read 1 adapter (3' end)
pub const ADAPTER_R1: &[u8] = b"AGATCGGAAGAGCACACGTCTGAACTCCAGTCA";
/// Illumina TruSeq Read 2 adapter (3' end)
pub const ADAPTER_R2: &[u8] = b"AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT";

#[derive(Debug, Clone)]
pub struct AdapterConfig {
    pub adapter_r1: Vec<u8>,
    pub adapter_r2: Vec<u8>,
    pub max_mismatch_rate: f64,
    pub min_overlap: usize,
}

impl Default for AdapterConfig {
    fn default() -> Self {
        Self {
            adapter_r1: ADAPTER_R1.to_vec(),
            adapter_r2: ADAPTER_R2.to_vec(),
            max_mismatch_rate: 0.1,
            min_overlap: 10,
        }
    }
}

/// Trim adapter from the 3' end of `seq`.
/// Returns the number of bases to keep from the 5' end.
pub fn trim_adapter_se(seq: &[u8], adapter: &[u8], cfg: &AdapterConfig) -> usize {
    let seq_len = seq.len();
    let adp_len = adapter.len();

    for overlap in (cfg.min_overlap..=seq_len.min(adp_len)).rev() {
        let seq_start = seq_len - overlap;
        let mismatches = seq[seq_start..]
            .iter()
            .zip(adapter[..overlap].iter())
            .filter(|&(a, b)| a != b)
            .count();
        let allowed = (overlap as f64 * cfg.max_mismatch_rate).ceil() as usize;
        if mismatches <= allowed {
            return seq_start;
        }
    }
    seq_len
}

fn rev_comp(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&b| match b {
            b'A' | b'a' => b'T',
            b'T' | b't' => b'A',
            b'G' | b'g' => b'C',
            b'C' | b'c' => b'G',
            other => other,
        })
        .collect()
}

/// Detect insert-based PE overlap. Returns `Some((r1_keep, r2_keep))` if found.
pub fn detect_pe_overlap(r1: &[u8], r2: &[u8], cfg: &AdapterConfig) -> Option<(usize, usize)> {
    let r2_rc = rev_comp(r2);
    let r1_len = r1.len();
    let r2_len = r2.len();

    for insert_size in (cfg.min_overlap..=r1_len.min(r2_len)).rev() {
        let mismatches = r1[..insert_size]
            .iter()
            .zip(r2_rc[..insert_size].iter())
            .filter(|&(a, b)| a != b)
            .count();
        let allowed = (insert_size as f64 * cfg.max_mismatch_rate).ceil() as usize;
        if mismatches <= allowed {
            return Some((insert_size, insert_size));
        }
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn trim_adapter_from_3_end() {
        let mut seq = b"ACGTACGTACGTACGTACGT".to_vec();
        seq.extend_from_slice(ADAPTER_R1);
        let cfg = AdapterConfig::default();
        let keep = trim_adapter_se(&seq, ADAPTER_R1, &cfg);
        assert_eq!(keep, 20);
    }

    #[test]
    fn no_adapter_untouched() {
        let seq = b"ACGTACGTACGTACGT";
        let cfg = AdapterConfig::default();
        assert_eq!(trim_adapter_se(seq, ADAPTER_R1, &cfg), seq.len());
    }
}
