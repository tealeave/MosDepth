# Coverage and Uniformity Analysis using MosDepth

This project provides scripts to perform coverage and uniformity analysis on sequencing data using BAM files and a target BED file. It leverages MosDepth to calculate per-base coverage and then processes this data to determine threshold coverage, mean coverage per target, run uniformity metrics (including the 95% central fold difference), and optionally generate visualizations like correlation heatmaps and coverage histograms.

## Features

*   **MosDepth Integration**: Executes MosDepth to generate per-base coverage data from BAM files.
*   **Comprehensive Coverage Metrics**: Calculates mean coverage per target and coverage at specified thresholds (e.g., 10x, 20x).
*   **Uniformity Assessment**: Computes the central 95% fold difference, a metric indicating coverage uniformity (lower values are better).
*   **Data Aggregation**: Compiles results across multiple samples, providing both individual and average metrics.
*   **Visualization (Optional)**:
    *   Generates a pairwise correlation heatmap of coverage between samples (`--heatmap`).
    *   Generates a histogram of average per-base coverage across samples (`--histogram`).
*   **Batch Processing**: Designed to handle multiple samples efficiently.

## Prerequisites

1.  **Conda/Miniconda**: You need Conda to manage the Python environment and dependencies. You can install Miniconda from [here](https://docs.conda.io/en/latest/miniconda.html).
2.  **MosDepth**: The core coverage calculation tool, MosDepth, must be installed separately. Please follow the installation instructions on the [MosDepth GitHub Repository](https://github.com/brentp/mosdepth). Ensure the `mosdepth` executable is in your system's PATH.
3.  **Git**: Required for cloning repositories if needed (listed in `dependency.yml`).

## Installation

1.  **Clone the repository (if you haven't already):**
    ```bash
    git clone <repository-url>
    cd MosDepth
    ```
2.  **Create and activate the Conda environment:**
    This command uses the `dependency.yml` file to install the required Python packages (like pandas, numpy, matplotlib, seaborn, openpyxl, etc.) into an environment named `ambry`.
    ```bash
    conda env create -f dependency.yml
    conda activate ambry
    ```

## Scripts

*   **`get_coverage_info.py`**: The main script for running the analysis pipeline.
*   **`per_base_convert.py`**: (Purpose not detailed in original README - potentially a helper script).
*   **`MosDepth.job`**: (Purpose not detailed in original README - likely a template for job submission systems like SLURM/PBS).

## Usage: `get_coverage_info.py`

This script orchestrates the coverage analysis process.

```bash
python get_coverage_info.py --bam_list <path/to/your/bam_list.txt> \
                            --bed <path/to/your/target_regions.bed> \
                            [--heatmap] \
                            [--histogram] \
                            [--outdir <path/to/output_directory>]
```

### Input Files

*   **`--bam_list BAM_LIST`**: (Required) A plain text file containing the full paths to your input BAM files, one file path per line.
    Example `bam_list.txt`:
    ```
    /path/to/sample1.bam
    /path/to/sample2.bam
    /path/to/sample3.bam
    ```
*   **`--bed BED`**: (Required) A BED file defining the target regions for coverage analysis. Crucially, the **4th column** must contain unique identifiers for each region (e.g., `$gene_exon` format as mentioned in the original README). MosDepth uses this file to calculate coverage specifically over these regions.
    Example BED format (tab-separated):
    ```
    chr1  1000  1500  GeneA_exon1
    chr1  2500  3000  GeneA_exon2
    chr2  5000  5500  GeneB_exon1
    ```

### Optional Arguments

*   **`--heatmap`**: Generate a PNG file (`correlation_heatmap.png`) showing the pairwise correlation of mean target coverage between samples.
*   **`--histogram`**: Generate a PNG file (`Avg_PerBaseCovHist.png`) showing the distribution of average per-base coverage across all samples.
*   **`--outdir OUTDIR`**: Specify a directory where all output files will be saved. If not provided, output files are saved in the current working directory.
*   **`-h, --help`**: Show the help message and exit.

### Output Files

The script generates the following files (typically in the specified `--outdir` or the current directory):

*   **`PerTargetMeanCov.xlsx`**: An Excel file containing the mean coverage calculated for each target region (defined in the BED file) for every input sample. Includes an average across samples.
*   **`PerTargetThresholdCov.xlsx`**: An Excel file showing the breadth of coverage for each target region per sample, indicating the percentage of bases covered at predefined thresholds (e.g., >=1x, >=10x, >=20x). Includes averages.
*   **`Central95FoldDiff.xlsx`**: An Excel file containing the calculated central 95% fold difference for each sample, providing a measure of coverage uniformity.
*   **`correlation_heatmap.png`**: (Optional: if `--heatmap` is used) Pairwise correlation heatmap image.
*   **`Avg_PerBaseCovHist.png`**: (Optional: if `--histogram` is used) Average per-base coverage histogram image.
*   **MosDepth Outputs**: The script will also generate intermediate files from MosDepth runs (e.g., `.mosdepth.global.dist.txt`, `.mosdepth.region.dist.txt`, `.regions.bed.gz`, `.per-base.bed.gz`) within the output directory. These are used for the downstream calculations.
