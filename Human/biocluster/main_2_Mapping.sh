#!/bin/bash
set -e

analysis_dir=/group/bioinf/Users/Madina/pipeline_tested/human
processed_files=$analysis_dir/1_Preprocessing/1.2_CuratedData/1.2_Quality_trimmed                                       
paired_end="No" 
filename_extention=.fastq.gz # extension of the file to work on. Data curation works on .fastq.gz
readname=R

echo 'Started'

#chmod u+x /group/bioinf/Data/Madina_test/pipeline/v1/Human/jobs/child_2_Mapping.sh

      
echo "2.1. Generating reference genome"
for s in $(cat $analysis_dir/samplenames.txt)
    do
        samplename="$s"
        echo "Processing sample: $samplename"
        sbatch -J $samplename.2_mapping $analysis_dir/biocluster/child_2_Mapping.sh $analysis_dir $processed_files $paired_end $filename_extention $samplename $readname
    done

echo "Finished at `date`"
echo "__DONE__"

   
echo "Finished at `date`"
echo "__DONE__"

