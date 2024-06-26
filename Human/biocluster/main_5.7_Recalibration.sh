#!/bin/bash
set -e
analysis_dir=/group/bioinf/Users/Madina/pipeline_tested/human

echo "Started at `date`"

        sbatch -J 5.7_Recalibration $analysis_dir/biocluster/child_5.7_VQSR.sh $analysis_dir 


echo "Finished at `date`"
echo "__DONE__"


