# Demo: 4T1 Initial vs Tail Vein Lung Metastasis

A small, subsampled version of `experiments/4T1_Tail_Vein_vs_Lung_Met.toml` for quickly testing the pipeline.

## Contents

```
demo/
├── demo_4T1_Initial_vs_Tail_Vein.toml   # Config
├── fastq/                               # R1/R2 FASTQs, 100,000 read pairs per sample
└── results/                             # Expected output (abundances, fitness bootstraps, figure)
```

Samples: `Initial_Tail_Vein_Injection_A` and `initial_Tail_Vein_Injection_B` (references), and lung metastases 1–3.

## Running

From `SatSeq_Bulk_DNAseq_Analysis/`, with [uv](https://github.com/astral-sh/uv) and [mmseqs2](https://github.com/soedinglab/MMseqs2) installed:

```bash
uv run python src/satseq_analysis.py demo/demo_4T1_Initial_vs_Tail_Vein.toml
```

This takes about 2 minutes on a laptop. The mapping step is skipped if `results/` already contains an abundance file. Delete `demo/results/` to rerun everything from the reads.

## Expected output

![Demo fitness heatmap and ranked barplot](results/figures/fitness_barplot_data.xlsx.png)

Top: fitness for each variant (amino acid × position). Bottom: substitutions ranked by fitness, with bootstrap intervals.

## How the reads were subsampled

For each sample, lanes from both sequencing runs were concatenated. 100,000 read pairs were then drawn at random, with R1 and R2 kept paired:

```bash
paste <(zcat <R1 files> | paste - - - -) <(zcat <R2 files> | paste - - - -) \
  | shuf -n 100000 --random-source=<(yes 42) > tmp.tsv
cut -f1-4 tmp.tsv | tr '\t' '\n' | gzip -9 > <sample>_S1_L001_R1_001.fastq.gz
cut -f5-8 tmp.tsv | tr '\t' '\n' | gzip -9 > <sample>_S1_L001_R2_001.fastq.gz
```

## Notes

- Only R2 is used by the pipeline. R1 is included for completeness.
- With subsampled reads, many variants fall below the 8-count filter, so fitness estimates are sparse and noisier than in the full experiment. The demo shows that the pipeline runs end to end; it does not reproduce the full results.
