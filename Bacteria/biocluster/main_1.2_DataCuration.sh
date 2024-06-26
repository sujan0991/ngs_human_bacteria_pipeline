#!/bin/bash
set -e

analysis_dir=/group/bioinf/Data/Madina_test/pipeline/v1
folder_nextstep=$analysis_dir/0_fastQ_data                          
paired_end="Yes" # "Yes" if reads are pair-ended
filename_extention=.fastq.gz
readname=R

echo 'Started'

#chmod u+x /group/bioinf/Data/Madina_test/pipeline/v1/biolcuster/child_1.2_DataCuration.sh

# Loop over each line in samplenames.txt
for s in $(cat $analysis_dir/samplenames.txt)
    do
        samplename="$s"
        echo "Processing sample: $samplename"
        sbatch -J 1.2_DataCuration.$samplename $analysis_dir/biocluster/child_1.2_DataCuration.sh $analysis_dir $folder_nextstep $paired_end $filename_extention $samplename $readname
    done

echo "Finished at `date`"
echo "__DONE__"
