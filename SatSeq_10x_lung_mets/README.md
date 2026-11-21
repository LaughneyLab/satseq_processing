# Lung Mets 10X analysis

This directory contains scripts needed to reproduce the analysis of the lung mets 10X data.

The 10x data here is regular scRNA-seq data; each sample has a different STING mutation (K337R, P13A, STOP, or REF). There are no custom barcodes in this 10x library. This analysis does not contain any novel methods.

## Workflow and required inputs

Run notebooks with **this directory as the working directory**. Sample data sets are provided in the `data` subdirectory;
full data sets are submitted to SRA/GEO linked to BioProject PRJNA1552726.

| Order | Notebook/script | Inputs | Outputs/use |
|---|---|---|---|
| 0 | `cellranger count` | RAW fastqs for each sample | `raw_feature_matrix.h5`, `molecule_info.h5`. Used cellranger-8.0.1 and mouse reference refdata-gex-mm10-2020A |
| 1 | `lung_mets_analysis.ipynb` | cellranger ouptut | filtered anndata object with cell type annotations, MILO results, MAST results |
| 1a | `cellAssign.R` | filtered anndata object, `data/lung_markers_with_4T1.csv`  | cell type assignment for each cell |
| 1b | `run_scDblFinder.R` | filtered anndata object | scDblFinder doublet probabilities per cell |
| 1c | `milo_analysis.ipynb` | filtered anndata object with neighborhhods | differential abundance results for pairiwse genotype comparisons |
| 1d | `MAST_wrapper.R` | filtered anndata object | differential expression for pairwise genotype comparisons |

The invocation of items `1a`, `1b`, `1c`, and `1d` above are described in the main notebook (item `1`), `lung_mets_analysis.ipynb`.

## Output annotations

The `lung_mets_analysis.ipynb` notebook creates the fields below. 

| Field | Meaning in the supplied notebook |
|---|---|
| `obs['cellAssign_0.8']` | Cell type assigned by cellAssign
| `obs['Macro Cell Type 4']` | Final cell types after refining cellAssign annotations |
| `obs['celltype_sub']` | Final cell type for each cell, including fine-grained cell types for T cell category |
| `layers['X_orig']` | XXXX |
| `layers['X']` | XXXX |
| `layers['logX']` | XXXX |


## Software environment

The notebooks contain version-reporting cells, which document all libraries used. 

The main notebook, `lung_mets_analysis.ipynb`, was executed in a docker environment, which can be downloaded and started with the following commands. Run from your project directory:
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
Relevant versions of python and packages in the main docker image include:

| Software | version |
| --- | --- |
|Python| 3.10.12|
|h5py| 3.10.0|
|matplotlib| 3.8.3|
|numpy| 1.26.4|
|pandas| 2.2.1|
|scanpy| 1.9.8|
|scipy| 1.12.0|
|scrublet| 0.2.3|
|seaborn| 0.13.2|
|textalloc| 1.1.8|
|scikit-learn| 1.4.1.post1|
| R | 4.3.2 |
| MAST | 1.28.0 |
| SingleCellExperiment | 1.24.0 |
| rhdf5 | 2.46.1 |
| Matrix | 1.6.3 |

The last command outputs information about the jupyter server, including a long token string. You can open a local web browser, navigate to [localhost:8888](localhost:8888), and enter this token to access the notebook.


The following three scripts were each run in separate docker environments. The main notebook, `lung_mets_analysis.ipynb`, includes instructions for launching these analyses. The `data` directory contains results for each analysis, so that you can also skip those sections and read the results:

| Main notebook section | code | docker image | Software versions in container |
| --- | --- | --- | --- | 
| Cell assign | `callCellAssign_gpu.sh` | `biohpc/scrna-gpu:0.0.1` (requires GPU) | R 4.0.2, cellassign 0.99.21, SingleCellExperiment 1.10.1, rhdf5 2.32.2 |
| scDblFinder | `run_scDblFinder.R` | `biohpc/scdblfinder:0.1` | R 4.6.0, scDblFinder 1.26.2, rhdf5 2.56.0 |
| MILO | `milo_analysis.ipynb` | `biohpc/milo:0.1` | R 4.4.2, miloR 2.0.0, zellkonverter 1.14.1, SingleCellExperiment 1.26.0, scater 1.32.1, Matrix 1.7.0, dplyr 1.1.4 |


For citation, license status, and questions, see the [main README](../README.md).
