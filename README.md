# SatSeq code repository

Analysis code supporting the manuscript:

**SatSeq: An *in vivo* genotype to phenotype map of the metastatic functions of STING**  
Deyell, M. *et al.* (manuscript in preparation)

SatSeq links engineered STING variants to cellular phenotypes using sequencing barcodes. This repository contains PacBio barcode-map construction, assignment of variants to 10x single cells, bulk DNA sequencing analysis, and downstream manuscript analyses.

## Start here

| Goal | Starting point |
|---|---|
| Reproduce manuscript analyses from processed results | [Manuscript figures README](SatSeq_submission_paper_figures/README.md), including the figure-to-analysis table |
| Construct barcode-to-variant maps and assign variants to single cells | [PacBio and 10x README](SatSeq_10x_Variant_Calling/README.md) |
| Regenerate bulk barcode counts and enrichment/fitness estimates | [Bulk DNA sequencing README](SatSeq_Bulk_DNAseq_Analysis/README.md) |
| Explore additional bulk sequencing analyses | [Bulk analysis notebooks README](SatSeq_Bulk_DNAseq_Analysis/analysis/README.md) |

Each component has its own software requirements and working-directory assumptions. All code was tested/run on linux systems running Rocky 9.8.

## Workflow overview

```text
PacBio HiFi BAMs -> per-run barcode-to-variant maps -> combined cell-line maps
                                                             |
10x FASTQs -> Cell Ranger -> per-cell barcode consensus -------+
                                                             |
                                                 cell-to-variant assignments
                                                             |
                                                 single-cell analyses/figures

Bulk DNA FASTQs -> barcode counts -> enrichment/fitness -> models/figures
```

The biological cell lines are **4T1** and **EO771.lmb**. Historical scripts and filenames also use `EO771` for EO771.lmb. 


## Questions and maintenance

Please use [GitHub Issues](https://github.com/LaughneyLab/satseq_processing/issues) for reproducibility questions or bug reports, including the release/commit, command, software environment, and error message. Maintainer/contact: @mjhubisz.
