# This repo is a collection of per target coverage and threshold coverage scripts with a given list of BAM files and a target bed file

1. wrapper_for_mosdepth.sh: the driver script that calls different script and manage jobs
    - run it like: `bash wrapper_for_mosdepth.sh $your_bed_file`
    - line 18: modify to glob your target BAM files
    - line 35 and 39: give the full path of your python env
2. mosdepth.job: jod template file for qsubing individul samples
    - run it like: `qsub -v bam=$bam -v bed=$BED mosdepth.job`
    - line 17: give the full path of the mosdepth on your python env
3. compile_mosdepth_exonStat.py: compiles all the mosdepth results from `*.regions.bed.gz`, `*.thresholds.bed.gz`, 
    and `*.mosdepth.summary.txt` from mosdepth and detects any target with any sample below the coverage threshold
    - run it like: `python compile_mosdepth_exonStat.py $coverage_threshold`
    - output: `PerTargetThresholdCov.xlsx` and `PerTargetThresholdCov.xlsx`
4. corr_heatmap.py: uses PerTargetThresholdCov.xlsx to plot in-run sample cov correlation heatmap
    - run it like: `python corr_heatmap.py PerTargetThresholdCov.xlsx`
    - output: `correlation_heatmap.png`
5. main packages installed in python env:
    - pandas, xlsxwriter, xlrd, seaborn, matplotlib