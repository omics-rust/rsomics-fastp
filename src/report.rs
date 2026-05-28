use serde::{Deserialize, Serialize};

#[derive(Debug, Default, Clone, Serialize, Deserialize)]
pub struct ReadStats {
    pub total_reads: u64,
    pub total_bases: u64,
    pub q20_bases: u64,
    pub q30_bases: u64,
}

impl ReadStats {
    pub fn q20_rate(&self) -> f64 {
        if self.total_bases == 0 {
            0.0
        } else {
            self.q20_bases as f64 / self.total_bases as f64
        }
    }

    pub fn q30_rate(&self) -> f64 {
        if self.total_bases == 0 {
            0.0
        } else {
            self.q30_bases as f64 / self.total_bases as f64
        }
    }
}

#[allow(non_snake_case)]
#[derive(Debug, Default, Serialize, Deserialize)]
pub struct FilteringResult {
    pub passed_filter_reads: u64,
    pub low_quality_reads: u64,
    /// Field name matches fastp JSON schema.
    pub too_many_N_reads: u64,
    pub too_short_reads: u64,
    pub too_long_reads: u64,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct Report {
    pub summary: Summary,
    pub filtering_result: FilteringResult,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct Summary {
    pub before_filtering: ReadStats,
    pub after_filtering: ReadStats,
}
