# SatSeq processing code

This repository contains analysis code supporting the manuscript:

**SatSeq: An *in vivo* genotype to phenotype map of the metastatic functions of STING**  
Deyell, M. *et al.* (manuscript in preparation)

The notebooks in this repository describe the long read processing and variant calling steps used in the SatSeq workflow.

## Notebooks

### `01_pacbio_parsing.ipynb`

Processes PacBio sequencing data from mutant STING constructs containing a downstream SatSeq barcode. The notebook parses the sequencing reads and generates a barcode-to-mutation map linking SatSeq barcodes to STING variants.

### `01b_combine_pacbio_maps.ipynb`

Combines barcode-to-mutation maps generated from multiple PacBio sequencing runs from plasmid and 4T1/EO771.lmb cell lines. This notebook was used to generate the combined maps used in the study, but is not generally relevant for future satseq analysis.

### `02_10x_variant_calling.ipynb`

Processes 10x single-cell sequencing data from the SatSeq experiments. The PacBio-derived barcode-to-mutation map is used to assign STING mutations to individual cells.

## Repository structure

```text
.
+-- 01_pacbio_parsing.ipynb
+-- 01b_combine_pacbio_maps.ipynb
+-- 02_10x_variant_calling.ipynb
+-- bash/
+-- python/
+-- data/
```

- `bash/`  Bash scripts used by the analysis notebooks.
- `python/`  Python helper functions used by the notebooks.
- `data/`  Small supporting datasets required to run the notebooks.

## Software environment

The notebooks report the versions of Python and key Python packages used for each analysis. The docker image with the environment can be downloaded with:
```
docker pull biohpc/scrna2024:0.3
```

and a jupyter notebook running on port 8888 can be started in this environment with a command like:
```bash
CID=$(docker run -d biohpc/scrna2024 /workdir)
docker exec -it $CID jupyter server list
```

Where `/workdir` is your working directory. The first command starts the container, the second command will print the link to the jupyter notebook, including the access token.


## Citation

If you use this code, please cite the associated SatSeq manuscript. Citation information will be updated once the manuscript is available.


    
