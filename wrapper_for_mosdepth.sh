#!/bin/bash

# set -eux
WORK_DIR=$PWD

# Give this script a list of run id in "runs" for batch processing
# Give your annotated bed files with the target name as the fourth column as the positional arg
BED=$1


foo(){
         local run=$1
         mkdir ${run}
         cd ${run}

         ls /NGS/Nextseq/${run}/PANEL/PROJECT/Sample_*/*.bam | grep -v NEG > mose_bam_lst
         for bam in `cat mose_bam_lst`
            do
               qsub -v bam=$bam -v bed=${WORK_DIR}/$BED ${WORK_DIR}/mosdepth.job
            done

         # Wait for all the jobs complete
         num_running=`qstat -r | grep Full | grep "mosdepth" | wc -l`
         while [ ${num_running} -ne 0 ]
         do
         echo "Running ${num_running} jobs for mosdepth"
         sleep 30
         num_running=`qstat -r | grep Full | grep "mosdepth" | wc -l`
         done

         # Give the full path of your python env here, this will compile the mosedepth result
         /home/dlin/miniconda3/envs/mose/bin/python3.7 ${WORK_DIR}/compile_mosdepth_exonStat.py
         mv PerTargetMeanCov.xlsx ${run}_PerTargetMeanCov.xlsx
         mv PerTargetThresholdCov.xlsx ${run}_PerTargetThresholdCov.xlsx

         cp *.xlsx ${WORK_DIR}
         cd ${WORK_DIR}
}

for r in `cat runs`; do foo $r & done
