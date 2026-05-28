# rsomics-fastp

Fast FASTQ quality control and preprocessing — clean-room Rust reimplementation of fastp.

## Usage

```
rsomics-fastp -i in.fq -o out.fq
rsomics-fastp -i in_R1.fq -I in_R2.fq -o out_R1.fq -O out_R2.fq
rsomics-fastp -i in.fq -o out.fq --trim-poly-g --json report.json
```

## Key flags

| Flag | Default | Description |
|------|---------|-------------|
| `-i/--in1` | required | Input FASTQ (R1 for PE) |
| `-I/--in2` | — | Input FASTQ R2 (paired-end) |
| `-o/--out1` | required | Output FASTQ R1 |
| `-O/--out2` | — | Output FASTQ R2 |
| `-j/--json` | — | JSON report path |
| `--qualified-quality-phred` | 15 | Minimum Phred score for "qualified" base |
| `--unqualified-percent-limit` | 40 | Max % unqualified bases per read |
| `--n-base-limit` | 5 | Max N bases per read |
| `--length-required` | 15 | Minimum read length after trimming |
| `--length-limit` | 0 | Maximum read length (0 = unlimited) |
| `--trim-poly-g` | off | Trim poly-G tails (NextSeq/NovaSeq) |
| `--trim-poly-x` | off | Trim poly-X tails for all four bases |
| `--poly-g-min-len` | 10 | Minimum run length to trigger poly trimming |
| `--disable-adapter-trimming` | off | Disable Illumina adapter trimming |
| `--adapter-sequence` | TruSeq R1 | Custom R1 adapter |
| `--adapter-sequence-r2` | TruSeq R2 | Custom R2 adapter |

## Origin

This crate is an independent Rust reimplementation of `fastp` based on:
- Chen, S. et al. "fastp: an ultra-fast all-in-one FASTQ preprocessor."
  *Bioinformatics* 34(17):i884–i890 (2018). DOI: 10.1093/bioinformatics/bty560
- The fastp public format specification and observed black-box behaviour
- fastp v0.20.1 source (MIT license) — used as reference for algorithm details

Algorithm provenance: the poly-X trimming algorithm (mismatch budget and
firstGPos semantics) matches the behaviour of fastp v0.20.1 `trimPolyG`.
The adapter trimming uses a suffix-scan overlap approach consistent with
fastp's default logic.

License: MIT OR Apache-2.0.
Upstream credit: [fastp](https://github.com/OpenGene/fastp) (MIT license).
