/// Trim a poly-X run from the 3′ end of a sequence.
///
/// Algorithm matches fastp v0.20.1 `trimPolyG` / `trimPolyX`:
/// - Scans from the 3′ end, allowing up to 1 mismatch per 8 bases examined
///   (integer division) with a hard cap of 5 total mismatches.
/// - `run_start` tracks the most recently seen matching base (leftmost boundary
///   of the run); mismatches do not advance this boundary.
/// - Trims if at least `min_len` bases were examined when the scan ended.
///
/// `min_len`: minimum run length to trigger trimming (fastp default 10).
pub fn trim_poly_x(seq: &[u8], base: u8, min_len: usize) -> usize {
    const ALLOW_ONE_MISMATCH_FOR_EACH: usize = 8;
    const MAX_MISMATCH: usize = 5;

    let len = seq.len();
    if len < min_len {
        return len;
    }

    let base_upper = base.to_ascii_uppercase();
    let base_lower = base.to_ascii_lowercase();

    let mut mismatches = 0usize;
    let mut run_start = len; // most recently seen matching-base position
    let mut i = 0usize; // number of positions examined (from 3' end)

    while i < len {
        let pos = len - i - 1;
        let b = seq[pos];

        if b == base_upper || b == base_lower {
            run_start = pos;
        } else {
            mismatches += 1;
        }

        // fastp break condition: max absolute mismatches OR rate exceeded
        let allowed_by_rate = (i + 1) / ALLOW_ONE_MISMATCH_FOR_EACH;
        if mismatches > MAX_MISMATCH || (mismatches > allowed_by_rate && i + 1 >= min_len) {
            break;
        }

        i += 1;
    }

    // Trim if the scan covered at least `min_len` positions
    if i < min_len {
        return len;
    }
    run_start
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn no_poly_g_untouched() {
        // No G run at 3′ end — should not trim
        let seq = b"ACGTACGTACGT";
        assert_eq!(trim_poly_x(seq, b'G', 10), seq.len());
    }

    #[test]
    fn poly_g_trimmed() {
        // 11 G's at 3′ end. The G at pos 6 of the prefix is within the mismatch
        // budget and gets absorbed into the run (fastp v0.20.1 behaviour verified).
        // fastp trims to 6, leaving "ACGTAC".
        let seq = b"ACGTACGTGGGGGGGGGGG";
        let keep = trim_poly_x(seq, b'G', 10);
        assert_eq!(keep, 6);
        assert_eq!(&seq[..keep], b"ACGTAC");
    }

    #[test]
    fn poly_g_shorter_than_min_untouched() {
        // Only 5 G's — below min_len 10 — should not trim
        let seq = b"ACGTACGTGGGGG";
        assert_eq!(trim_poly_x(seq, b'G', 10), seq.len());
    }

    #[test]
    fn poly_g_with_one_mismatch_in_run() {
        // "AAAGGGCGGGGGGGGGGG" — 1 C mismatch inside a 16-char G-dominant run.
        // i=0..15: 1 mismatch in 15 bases = allowed=15/8=1 => not exceeded; continue.
        // At i=15: the next char is A (pos 2). mismatches=2, allowed=16/8=2. 2>2? no.
        // At i=16: A (pos 1). mismatches=3, allowed=17/8=2. 3>2 AND 17>=10 => break.
        // run_start = position of last G seen = 3 (the G at pos 3, 0-indexed).
        // Keep = 3 (trim everything from pos 3 onward).
        let seq = b"AAAGGGCGGGGGGGGGGG";
        let keep = trim_poly_x(seq, b'G', 10);
        assert!(keep <= 3, "expected keep<=3, got {keep}");
    }

    #[test]
    fn poly_x_c_trimmed() {
        // 12 C's at 3′ end
        let seq = b"ACGTACGTCCCCCCCCCCCC";
        let keep = trim_poly_x(seq, b'C', 10);
        assert_eq!(keep, 8);
        assert_eq!(&seq[..keep], b"ACGTACGT");
    }

    #[test]
    fn entire_sequence_is_poly_g() {
        let seq = b"GGGGGGGGGGGGGGGG";
        assert_eq!(trim_poly_x(seq, b'G', 10), 0);
    }

    #[test]
    fn short_seq_below_min_len() {
        let seq = b"GGGGG"; // 5 G's, min_len=10
        assert_eq!(trim_poly_x(seq, b'G', 10), seq.len());
    }
}
