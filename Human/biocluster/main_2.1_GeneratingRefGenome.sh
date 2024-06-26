#!/bin/bash
set -e
analysis_dir=/group/bioinf/Users/Madina/pipeline_tested/human

echo 'Started'

#chmod u+x /group/bioinf/Data/Madina_test/pipeline/v1/Human/jobs/child_2.1_GeneratingRefGenome.sh

      
echo "2.1. Generating reference genome"
sbatch -J 2.1_GeneratingRefGenome $analysis_dir/biocluster/child_2.1_GeneratingRefGenome.sh $analysis_dir 
   
echo "Finished at `date`"
echo "__DONE__"

