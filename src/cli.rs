use std::path::PathBuf;

use clap::Parser;
use rsomics_common::{CommonFlags, Result, RsomicsError, Tool, ToolMeta};
use rsomics_fastp::adapter::AdapterConfig;
use rsomics_fastp::filter::FilterConfig;
use rsomics_fastp::run::{RunConfig, run_pe, run_se, write_report};

pub const META: ToolMeta = ToolMeta {
    name: env!("CARGO_PKG_NAME"),
    version: env!("CARGO_PKG_VERSION"),
};

#[derive(Parser, Debug)]
#[command(
    name = "rsomics-fastp",
    version,
    about = "FASTQ quality control and preprocessing",
    disable_help_flag = true
)]
pub struct Cli {
    #[command(flatten)]
    pub common: CommonFlags,

    /// Input FASTQ (R1 for PE).
    #[arg(short = 'i', long = "in1")]
    pub in1: PathBuf,

    /// Input FASTQ R2 (paired-end).
    #[arg(short = 'I', long = "in2")]
    pub in2: Option<PathBuf>,

    /// Output FASTQ R1.
    #[arg(short = 'o', long = "out1")]
    pub out1: PathBuf,

    /// Output FASTQ R2 (required for paired-end).
    #[arg(short = 'O', long = "out2")]
    pub out2: Option<PathBuf>,

    /// Write JSON report to this path (fastp-compatible schema).
    #[arg(short = 'j', long = "report-json")]
    pub report_json: Option<PathBuf>,

    /// Minimum base quality (Phred) to count as qualified.
    #[arg(long = "qualified-quality-phred", default_value_t = 15u8)]
    pub qualified_quality_phred: u8,

    /// Maximum fraction (%) of unqualified bases allowed per read.
    #[arg(long = "unqualified-percent-limit", default_value_t = 40u8)]
    pub unqualified_percent_limit: u8,

    /// Maximum number of N bases allowed per read.
    #[arg(long = "n-base-limit", default_value_t = 5usize)]
    pub n_base_limit: usize,

    /// Minimum read length after trimming.
    #[arg(long = "length-required", default_value_t = 15usize)]
    pub length_required: usize,

    /// Maximum read length (0 = unlimited).
    #[arg(long = "length-limit", default_value_t = 0usize)]
    pub length_limit: usize,

    /// Enable poly-G tail trimming (NextSeq/NovaSeq artifact).
    #[arg(long = "trim-poly-g")]
    pub trim_poly_g: bool,

    /// Enable poly-X tail trimming for all four bases.
    #[arg(long = "trim-poly-x")]
    pub trim_poly_x: bool,

    /// Minimum poly-G/X run length to trigger trimming.
    #[arg(long = "poly-g-min-len", default_value_t = 10usize)]
    pub poly_g_min_len: usize,

    /// Disable adapter trimming.
    #[arg(long = "disable-adapter-trimming")]
    pub disable_adapter_trimming: bool,

    /// Custom R1 adapter sequence (overrides Illumina TruSeq default).
    #[arg(long = "adapter-sequence")]
    pub adapter_sequence: Option<String>,

    /// Custom R2 adapter sequence (overrides Illumina TruSeq default).
    #[arg(long = "adapter-sequence-r2")]
    pub adapter_sequence_r2: Option<String>,
}

impl Tool for Cli {
    fn meta() -> ToolMeta {
        META
    }

    fn common(&self) -> &CommonFlags {
        &self.common
    }

    fn execute(self) -> Result<()> {
        let filter = FilterConfig {
            min_qual: self.qualified_quality_phred,
            max_low_qual_fraction: self.unqualified_percent_limit as f64 / 100.0,
            max_n_bases: self.n_base_limit,
            min_length: self.length_required,
            max_length: self.length_limit,
        };

        let mut adapter = AdapterConfig::default();
        if let Some(ref seq) = self.adapter_sequence {
            adapter.adapter_r1 = seq.as_bytes().to_vec();
        }
        if let Some(ref seq) = self.adapter_sequence_r2 {
            adapter.adapter_r2 = seq.as_bytes().to_vec();
        }

        let cfg = RunConfig {
            filter,
            adapter,
            adapter_trimming: !self.disable_adapter_trimming,
            trim_poly_g: self.trim_poly_g,
            trim_poly_x: self.trim_poly_x,
            poly_min_len: self.poly_g_min_len,
        };

        let result = if let (Some(in2), Some(out2)) = (&self.in2, &self.out2) {
            run_pe(&self.in1, in2, &self.out1, out2, &cfg)
                .map_err(|e| RsomicsError::UpstreamError(e.to_string()))?
        } else if self.in2.is_some() != self.out2.is_some() {
            return Err(RsomicsError::ConfigError(
                "--in2 and --out2 must be provided together".into(),
            ));
        } else {
            run_se(&self.in1, &self.out1, &cfg)
                .map_err(|e| RsomicsError::UpstreamError(e.to_string()))?
        };

        if let Some(ref json_path) = self.report_json {
            write_report(&result, json_path)
                .map_err(|e| RsomicsError::UpstreamError(e.to_string()))?;
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use clap::CommandFactory;

    #[test]
    fn cli_definition_is_valid() {
        super::Cli::command().debug_assert();
    }
}
