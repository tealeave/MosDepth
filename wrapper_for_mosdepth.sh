#!/bin/bash

# set -eux
WORK_DIR=$PWD

# Change line 18 to glob the bam files correctly
# Give your annotated bed files with cords like "chr1    10292381        10292497        Gene_RefSeqNum_ExonNum"
# Run it like "bash wrapper_for_mosdepth.sh $your_bed_file"
BED=$1
SEQUENCING_PLATFORM=


foo(){
         local run=$1
         mkdir ${run}
         cd ${run}

         ls /NGS/${SEQUENCING_PLATFORM}/${run}/Panel/Project/Sample*/*.bam | grep -v NEG > mose_bam_lst
         for bam in `cat mose_bam_lst`
            do
               qsub -v bam=$bam -v bed=${WORK_DIR}/$BED ${WORK_DIR}/mosdepth.job
            done

         # Wait for all the jobs to complete
         num_running=`qstat -r | grep Full | grep "mosdepth" | wc -l`
         while [ ${num_running} -ne 0 ]
         do
         echo "Running ${num_running} jobs for mosdepth"
         sleep 30
         num_running=`qstat -r | grep Full | grep "mosdepth" | wc -l`
         done

         # Give the full path of your python env here, this will compile the mosedepth result and detects targets with coverage below the threshold
         # In this case,  20x in ${run}_PerTargetMeanCov.xlsx
         /home/dlin/miniconda3/envs/mose/bin/python3.7 ${WORK_DIR}/compile_mosdepth_exonStat.py 20
         mv PerTargetMeanCov.xlsx ${run}_PerTargetMeanCov.xlsx
         mv PerTargetThresholdCov.xlsx ${run}_PerTargetThresholdCov.xlsx

         /home/dlin/miniconda3/envs/mose/bin/python3.7 ${WORK_DIR}/corr_heatmap.py ${run}_PerTargetMeanCov.xlsx
         mv correlation_heatmap.png ${run}_correlation_heatmap.png

         cp *.xlsx ${run}_correlation_heatmap.png ${WORK_DIR}
         cd ${WORK_DIR}
}

for r in `cat runs`; do foo $r & done
