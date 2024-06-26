#!/bin/bash
set -e
analysis_dir=/group/bioinf_biomarkers_rna/rna_us_analysis

echo "Started at `date`"

        sbatch -J 5.8_HardFiltering $analysis_dir/biocluster/child_5.8_HardFiltering.sh $analysis_dir 


echo "Finished at `date`"
echo "__DONE__"


