# SatSeq Processing and Manuscript Figures

This repository contains the processed datasets, analysis code, and figure-generation scripts for the SatSeq manuscript. It is organized into two top-level parts:

## Directory overview

```
SatSeq_submission_paper_figures/
├── datasets/
│   ├── FACS/
│   ├── bulk_fitness/
│   └── single_cell/
└── scripts_results_plots/
    ├── FACS_position_subs_modeling/
    ├── fitness_position_subs_modeling/
    ├── mSTING_Stability_predictions/
    ├── mutant_antigenicity/
    ├── cGAMP_sc/
    ├── tumor_sc_invivo/
    └── P13_K337_further/
        ├── RNAseq_bulk_time_points/
        ├── mass_spec_data_analysis/
        ├── mets_counts_survival/
        └── sc_analysis/
```

## Datasets

`datasets/` contains the common processed single-cell, FACS, and tumor fitness datasets, each in its own folder.

| Folder | Contents |
|---|---|
| `FACS/` | Collected per-variant FACS fitness results for all 4 readouts (STING, IFN, NFkB, ATF3), covering 3 replicates pooled across barcodes. |
| `bulk_fitness/` | Collected per-variant tumor fitness results matched across pairs of conditions (e.g., primary vs. mets) for the 4T1 and EO771 cell lines, covering all replicates pooled across barcodes. |
| `single_cell/` | MAST differential gene expression results for different conditions (e.g., cGAMP vs. H2O) and cell lines, plus the scripts to run the MAST jobs that generate these results from `adata_*.h5ad` files. |

The code and barcode-level data required to generate the `FACS/` and `bulk_fitness/` results are under `satseq_processing/SatSeq_Bulk_Analysis/`.

## Analysis scripts, results, and plots

`scripts_results_plots/` contains analysis-specific folders that use the datasets above and produce the final figures used in the manuscript.

| Folder | Description | Figures produced |
|---|---|---|
| `FACS_position_subs_modeling/` | Analysis and results for the hierarchical Bayesian mixed-effect model fit of FACS data, yielding per-position additive and interaction weights and amino acid substitution effect matrices. | Fig. 1f; Extended Data Figs. 2 and 4 (all panels) |
| `fitness_position_subs_modeling/` | Analysis and results for the hierarchical Bayesian mixed-effect model fit of tumor fitness data, yielding per-position additive and interaction weights and amino acid substitution effect matrices. | Fig. 3a, b, d; Extended Data Figs. 6a-l, 7d-g |
| `mSTING_Stability_predictions/` | In silico protein stability predictions for STING mutants from various models, plus the hierarchical Bayesian mixed-effect model fit that yields per-position additive and interaction weights and amino acid substitution effect matrices. | Extended Data Fig. 3 |
| `mutant_antigenicity/` | Runs MHCflurry to predict MHC binding of antigens derived from STING mutants. | Extended Data Figs. 6m-o |
| `cGAMP_sc/` | Analyzes FACS and MAST scores for cGAMP stimulation data to call and plot loss-of-function residues. | Fig. 2b, c, e, f; Extended Data Fig. 5a-e |
| `tumor_sc_invivo/` | Runs hotspot analysis of tumor single-cell data, and analyzes and plots hotspot module scores and gene correlations. Also contains structural and bond data for STING-TBK1-IRF3 complexes. | Fig. 4a-c, 5a; Extended Data Fig. 8 |
| `P13_K337_further/` | Additional data analysis for the P13 and K337 subclone experiments (see subfolders below). | See below |

### `P13_K337_further/` subfolders

| Subfolder | Description | Figures produced |
|---|---|---|
| `RNAseq_bulk_time_points/` | Bulk RNA sequencing time-series data and analysis. | Fig. 4d; Extended Data Fig. 9 |
| `mass_spec_data_analysis/` | Mass spectrometry data and analysis. | Fig. 5e; Extended Data Fig. 10a, b |
| `mets_counts_survival/` | Metastatic lesion count data and survival data, with their analysis and plotting. | Fig. 3e, f |
| `sc_analysis/` | Hotspot analysis of tumor and microenvironment single-cell data, including analysis and plotting of hotspot module scores, gene correlations, and immune response enrichment analysis (IREA). | Fig. 4e, h |

## Suggested workflow

1. Start from the barcode-level data and code in `../SatSeq_Bulk_RNASeq_Analysis/` to regenerate the pooled FACS and bulk tumor fitness results, or use the processed outputs directly from `datasets/FACS/` and `datasets/bulk_fitness/`.
2. For single-cell results, run the MAST scripts in `datasets/single_cell/` on the `adata_*.h5ad` files to regenerate the differential expression tables.
3. Run the analysis folders under `scripts_results_plots/` to reproduce the manuscript figures listed above.
