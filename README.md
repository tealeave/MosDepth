# Coverage and Uniformity Analysis with MosDepth

This script is designed to perform coverage analysis from BAM files and a target BED file. It calculates per base coverage, threshold coverage, run uniformity, and 95% central fold difference. Additionally, it can generate a correlation heatmap and a per base coverage histogram.

## DEPENDENCIES
- `dependency.yml`
- [MosDepth GitHub Repository](https://github.com/brentp/mosdepth)

## Scripts

### `get_coverage_info.py`

#### Usage:
```
get_coverage_info.py [-h] --bam_list BAM_LIST --bed BED [--heatmap]
                            [--histogram] [--outdir OUTDIR]
```

#### Available Arguments:
- `-h, --help`: Show this help message and exit.
- `--bam_list BAM_LIST`: [Input] List of BAM files for coverage analysis.
- `--bed BED`: [Input] Target BED file where the 4th column has the `$gene_exon` info.
- `--heatmap`: [Optional] Generate coverage correlation heatmap.
- `--histogram`: [Optional] Generate per base coverage histogram.
- `--outdir OUTDIR`: [Optional] Path to the output directory.

#### Output:
- `PerTargetThresholdCov.xlsx`
- `PerTargetMeanCov.xlsx`
- `Central95FoldDiff.xlsx`
- `correlation_heatmap.png`
- `Avg_PerBaseCovHist.png`

## Features

1. **MosDepth Job Submission**: The script can either spawn child processes to run MosDepth or submit MosDepth jobs for batch processing, especially useful when dealing with a large number of samples.

2. **Coverage Compilation**: It compiles per target and threshold results, providing both individual sample data and averages across samples.

3. **Pairwise Correlation Heatmap**: If the `--heatmap` argument is provided, the script will plot a pairwise correlation heatmap.

4. **Central 95% Fold Difference Calculation**: The script calculates the central 95% fold difference, which is a metric reflecting the spread and coverage uniformity. The smaller the value, the better the uniformity.

5. **Per Base Coverage Histogram**: If the `--histogram` argument is provided, the script will plot a histogram of per base coverage.
