#!/bin/bash
set -e
analysis_dir=/group/bioinf/Data/Madina_test/pipeline/v1
input_vcf=$analysis_dir/5_VariantCalling/5.3_Genotyping/variants.vcf

echo "Started at `date`"

echo "Started variants filtering"
sbatch $analysis_dir/biocluster/child_5_Filtering.sh $analysis_dir $input_vcf
    
echo "Finished at `date`"
echo "__DONE__"


