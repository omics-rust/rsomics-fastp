# rsomics-fastp

Fast FASTQ quality control and preprocessing.

Reads one (SE) or two (PE) FASTQ files, applies adapter trimming, quality filtering,
length filtering, poly-G/X trimming, and N-base filtering. Writes passing reads to
output FASTQ file(s) and optionally a JSON summary report.

## Usage

```
rsomics-fastp -i in.fastq -o out.fastq [OPTIONS]
rsomics-fastp -i in_R1.fastq -I in_R2.fastq -o out_R1.fastq -O out_R2.fastq [OPTIONS]
```

## Origin

This crate reimplements the core filtering and adapter-trimming logic of
[fastp](https://github.com/OpenGene/fastp) (MIT license) based on:

- The fastp paper: Chen et al. (2018) *Bioinformatics* 34(17):i884–i890.
  DOI: 10.1093/bioinformatics/bty560
- The public file-format specification (FASTQ)
- Black-box behaviour testing against fastp 0.20.x

No source code from the upstream was used as reference during implementation
beyond the algorithm description in the paper.

License: MIT OR Apache-2.0.
Upstream credit: fastp <https://github.com/OpenGene/fastp> (MIT).
