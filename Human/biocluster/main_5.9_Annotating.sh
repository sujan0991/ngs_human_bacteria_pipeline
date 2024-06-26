#!/bin/bash
set -e
analysis_dir=/group/bioinf_biomarkers_rna/rna_us_analysis

echo "Started at `date`"

        sbatch -J 5.9_Annotating $analysis_dir/biocluster/child_5.9_Annotating.sh $analysis_dir 


echo "Finished at `date`"
echo "__DONE__"


