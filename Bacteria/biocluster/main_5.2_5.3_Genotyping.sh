#!/bin/bash
set -e
analysis_dir=/group/bioinf/Data/Madina_test/pipeline/v1/

echo 'Started'
sbatch $analysis_dir/biocluster/child_5.2_5.3_Genotyping.sh $analysis_dir 

echo "Finished at `date`"
echo "__DONE__"

