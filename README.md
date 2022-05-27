# Coverage and Uniformity analysis with MosDepth
this script calculates per base coverage, threshold coverage, run uniformity, 95% central fold diff, 
and plots correlation heatmap and per base coverage histogram

## DEPENDENCIES
`dependency.yml`<br />
`https://github.com/brentp/mosdepth`

## Scripts

```
get_coverage_info.py
    usage:
    get_coverage_info.py [-h] --bam_list BAM_LIST --bed BED [--heatmap]
                                [--histogram] [--outdir OUTDIR]

    available arguments:
    -h, --help           show this help message and exit
    --bam_list BAM_LIST  [Input] list of bam fils for cov analysis
    --bed BED            [Input] target bed file where 4th column has the $gene_exon info
    --heatmap            [Optional] genrate cov correlation heatmap
    --histogram          [Optional] genrate per base coverage histogram
    --outdir OUTDIR      [Optional] Path to output directory

    output: PerTargetThresholdCov.xlsx, PerTargetThresholdCov.xlsx, correlation_heatmap.png,
            Central95FoldDiff.xlsx, Avg_PerBaseCovHist.png


```