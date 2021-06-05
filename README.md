# This repo is a collection of per target coverage and threshold coverage scripts

1. wrapper_for_mosdepth.sh: the driver script that calls different script and manage jobs
    - run it like: `bash wrapper_for_mosdepth.sh $your_bed_file`
    - line 16: modify to glob your target bam files
    - line 32: give the full path of the python env that you set up with the packages in "env_packages"
2. mosdepth.job: jod template file for qsubing individul samples
    - run it like: `qsub -v bam=$bam -v bed=$BED mosdepth.job`
    - line 17: give the full path of the python env that you set up with the packages in "env_packages"
3. compile_mosdepth_exonStat.py: compiles all the mosdepth results from `*.regions.bed.gz`, `*.thresholds.bed.gz`, and `*.mosdepth.summary.txt`
    - run it like: python compile_mosdepth_exonStat.py
    - output: PerTargetThresholdCov.xlsx and PerTargetThresholdCov.xlsx