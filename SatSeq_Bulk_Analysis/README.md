# SatSeq

Saturation mutagenesis sequencing analysis pipeline.

## Requirements

Install [uv](https://github.com/astral-sh/uv), create an environment from pyproject.toml:

```bash
uv sync
```

For mapping: Install [mmseqs2](https://github.com/soedinglab/MMseqs2)

## Project Structure

```
SatSeq/
├── src/
│   └── satseq_analysis.py    # Main pipeline script
├── experiments/              # TOML config files for each experiment
├── assets/                   # Barcode maps and reference files
│   ├── bcmap_4T1_*.txt.gz
│   ├── bcmap_EO771_*.txt.gz
│   └── SSBC_Full.csv
└── results/                  # Output directory (created automatically)
```

## Usage

Run a single experiment from the project root directory:

```bash
uv run python src/satseq_analysis.py experiments/<config>.toml
```

Run all experiments in a directory:

```bash
uv run python src/satseq_analysis.py experiments/
```

To overwrite existing bootstrap results, add the `--overwrite` flag:

```bash
uv run python src/satseq_analysis.py experiments/<config>.toml --overwrite
uv run python src/satseq_analysis.py experiments/ --overwrite
```

## Available Experiments

### 4T1 Cell Line - Passage Experiments

| Experiment | Config File | Description |
|------------|-------------|-------------|
| Initial vs MitogenHigh p5 | `experiments/4T1_Initial_vs_MitogenHigh_p5.toml` | Passage 1 vs Passage 5 (high mitogen) |
| Initial vs MitogenHigh p10 | `experiments/4T1_Initial_vs_MitogenHigh_p10.toml` | Passage 1 vs Passage 10 (high mitogen) |
| Initial vs MitogenHigh p20 | `experiments/4T1_Initial_vs_MitogenHigh_p20.toml` | Passage 1 vs Passage 20 (high mitogen) |
| Initial vs MitogenLow p5 | `experiments/4T1_Initial_vs_MitogenLow_p5.toml` | Passage 1 vs Passage 5 (low mitogen/2% FBS) |
| Initial vs MitogenLow p10 | `experiments/4T1_Initial_vs_MitogenLow_p10.toml` | Passage 1 vs Passage 10 (low mitogen/2% FBS) |
| Initial vs MitogenLow p20 | `experiments/4T1_Initial_vs_MitogenLow_p20.toml` | Passage 1 vs Passage 20 (low mitogen/2% FBS) |
| MitogenHigh p5 vs p10 | `experiments/4T1_MitogenHigh_p5_vs_p10.toml` | Passage 5 vs 10 (high mitogen, paired) |
| MitogenHigh p10 vs p20 | `experiments/4T1_MitogenHigh_p10_vs_p20.toml` | Passage 10 vs 20 (high mitogen, paired) |
| MitogenLow p5 vs p10 | `experiments/4T1_MitogenLow_p5_vs_p10.toml` | Passage 5 vs 10 (low mitogen, paired) |
| MitogenLow p10 vs p20 | `experiments/4T1_MitogenLow_p10_vs_p20.toml` | Passage 10 vs 20 (low mitogen, paired) |

### 4T1 Cell Line - STING FACS Experiments

| Experiment | Config File | Description |
|------------|-------------|-------------|
| STING FACS All vs Down | `experiments/4T1_STING_FACS_All_vs_Down.toml` | Total vs Down-sorted (paired) |
| STING FACS All vs Up | `experiments/4T1_STING_FACS_All_vs_Up.toml` | Total vs Up-sorted (paired) |

### 4T1 Cell Line - In Vivo Experiments

| Experiment | Config File | Description |
|------------|-------------|-------------|
| Mammary Fat Pad Initial vs Met | `experiments/4T1_Mammary_Fat_Pad_Initial_vs_Met.toml` | Initial vs Lung Metastasis |
| Mammary Fat Pad Initial vs Primary | `experiments/4T1_Mammary_Fat_Pad_Initial_vs_Primary.toml` | Initial vs Primary Tumor |
| Mammary Fat Pad Primary vs Met | `experiments/4T1_Mammary_Fat_Pad_Primary_vs_Met.toml` | Primary vs Metastasis (paired) |
| Tail Vein vs Lung Met | `experiments/4T1_Tail_Vein_vs_Lung_Met.toml` | Tail vein injection vs Lung Met |
| Ref Stop C147 K159 | `experiments/4T1_Ref_Stop_C147_K159_Initial_vs_Tail_Vein.toml` | Control experiment |

### EO771 Cell Line - Passage Experiments

| Experiment | Config File | Description |
|------------|-------------|-------------|
| Initial vs MitogenHigh p5 | `experiments/EO771_Initial_vs_MitogenHigh_p5.toml` | Passage 1 vs Passage 5 |
| Initial vs MitogenHigh p10 | `experiments/EO771_Initial_vs_MitogenHigh_p10.toml` | Passage 1 vs Passage 10 |
| Initial vs MitogenHigh p20 | `experiments/EO771_Initial_vs_MitogenHigh_p20.toml` | Passage 1 vs Passage 20 |
| MitogenHigh p5 vs p10 | `experiments/EO771_MitogenHigh_p5_vs_p10.toml` | Passage 5 vs 10 (paired) |
| MitogenHigh p10 vs p20 | `experiments/EO771_MitogenHigh_p10_vs_p20.toml` | Passage 10 vs 20 (paired) |

### EO771 Cell Line - In Vivo Experiments

| Experiment | Config File | Description |
|------------|-------------|-------------|
| Initial vs Tail Vein | `experiments/EO771_Initial_vs_Tail_Vein.toml` | Initial vs Tail vein injection |
| Mammary Fat Pad Initial vs Primary | `experiments/EO771_Mammary_Fat_Pad_Initial_vs_Primary.toml` | Initial vs Primary tumor |

## Example

Run the 4T1 Initial vs MitogenHigh p5 experiment:

```bash
uv run python src/satseq_analysis.py experiments/4T1_Initial_vs_MitogenHigh_p5.toml
```

Results will be saved to `results/4T1_Initial_vs_MitogenHigh_p5/`.

## Configuration

Each TOML config file specifies:
- `experiment_name`: Name for output files
- `output_dir`: Where results are saved (relative to project root)
- `barcodes_file`: Long-read barcode map (in `assets/`)
- `twist_file`: TWIST positional barcodes (in `assets/`)
- `data_dirs`: Directories containing FASTQ files
- `reference_names`: Reference sample prefixes
- `alt_names`: Treatment/alternate sample prefixes
- `dt`: Time unit for fitness calculation
- `paired`: Whether samples are paired (1:1 comparison)

## Analyses
Analysis and visualizations are in the `analysis/` directory.
