# PacBio barcode maps and 10x variant assignment

This directory contains the PacBio and 10x processing notebooks for the [SatSeq manuscript](../README.md). They build barcode-to-STING-variant maps and use those maps to label cells in single-cell RNA sequencing experiments. Historical `EO771` identifiers refer to EO771.lmb.

## Workflow and required inputs

Run notebooks with **this directory as the working directory**, in the order below. Sample data sets are provided in the `data` subdirectory;
full data sets are submitted to SRA/GEO linked to BioProject PRJNA1552726.

| Order | Notebook/script | Inputs | Outputs/use |
|---|---|---|---|
| 1 | `01_pacbio_parsing.ipynb` | One PacBio HiFi BAM, `data/STING.fa`, and the feature reference | Creates read/variant tables and barcode-to-mutation maps for a single pacbio run; set `bamfile`, `outdir`, and sample parameters for each run |
| 2 | `01b_combine_pacbio_maps.ipynb` | Per-run plasmid and cell-line barcode maps | Combined barcode-variant maps for 4T1 and EO771.lmb |
| 3 | `bash/run_cellranger.sh` | Matching gene-expression/custom-feature FASTQs, transcriptome reference, and feature reference | Cell Ranger count matrices, molecule information, and BAM output |
| 4 | `bash/bowtie_pileup.sh` | A capture's `outs/possorted_genome_bam.bam` and `outs/raw_feature_bc_matrix/barcodes.tsv.gz` | Barcode alignment BAM, consensus FASTA, and coverage table in that capture's `outs/` directory |
| 5 | `02_10x_variant_calling.ipynb` | Cell Ranger outputs, per-cell barcode consensus files, and the matching combined PacBio map | QC-filtered single-cell data with variant annotations, saved as `adata_{label}.h5ad` |




## Output annotations

The `02_10x_variant_calling.ipynb` notebook creates the fields below. 

| Field | Meaning in the supplied notebook |
|---|---|
| `obs['sample']` | Original capture identifier |
| `obs['variant_cr']` | Assignment from Cell Ranger custom-feature counts/reads after the notebook's criteria and reference-label corrections |
| `obs['consensus_barcode']` | Extracted combined barcode region from the per-cell consensus |
| `obs['variant_pacbio']` | Assignment supported by the PacBio barcode map |
| `obs['variant_substitution']` | Combined variant assignment before reduction to position-level labels |
| `obs['variant_pos']` | Position-level representation created by this notebook |
| `obs['haveVariant']` | Whether `variant` differs from the literal string `'None'` |
| `obsm['variant_counts']` / `obsm['variant_reads']` | Custom-feature UMI/read matrices separated from gene-expression columns |
| `layers['X_orig']` | raw UMI counts |
| `layers['X']` | normalized UMI counts  (each cell sum to 10000) |
| `layers['logX']` | log10-normalized UMI counts (log10(X+1) where X is normalized UMI) |


## Software environment

The notebooks contain version-reporting cells, which document all libraries used. They were all executed in a docker environment, which can be downloaded and started with the following commands. Run from your project directory:
```bash
docker pull biohpc/scrna2024:0.3
container=$((docker run --rm -it -d \
  -v $(pwd):/data \
  --workdir /data \
  -p 127.0.0.1:8888:8888 \
  --entrypoint jupyter \
  biohpc/scrna2024:0.3 \
  lab --ip=0.0.0.0 --port=8888 --no-browser --allow-root))
docker exec $container jupyter server list
```
The last command outputs information about the jupyter server, including a long token string. You can open a local web browser, navigate to [localhost:8888](localhost:8888), and enter this token to access the notebook.

Alignment of PacBio and 10x reads was done outside the docker environment, on Rocky 9.8 linux system, with:

| Tool | Evidence in the supplied code / version to confirm |
|---|---|
| Cell Ranger | 8.0.1 |
| Transcriptome | `refdata-gex-mm10-2020-A`|
| Bowtie2 | Required for barcode alignment; version: 2.5.1 |
| SAMtools | Required for BAM processing/consensus and compression; version 1.20 |

## Small example and validation

The `data/examples` directory contains small example data sets suitable for testing the jupyter notebooks.

For citation, license status, and questions, see the [main README](../README.md).
