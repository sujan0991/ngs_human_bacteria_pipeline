#!/bin/bash
analysis_dir=/group/bioinf/Users/Madina/pipeline_tested/human

folder_nextstep=$analysis_dir/0_fastQ_data_links # folder_nextstep is where the most recent .fastq.gz files will be stored. Includes subfolders "R1" and "R2" if reads are paired-end.
                                           # folder_nextstep will be updated as we move from one process to another; if process is skipped folder_nextstep is not updated
                                           # we begin where our .fastq.gz files ended up in 1.1_QualityCheckRaw: $analysis_dir/0_fastQ_data                                 
paired_end="No" # "Yes" if reads are pair-ended
filename_extention=.fastq.gz
readname=R

echo 'Started'

#chmod u+x /group/bioinf/Data/Madina_test/pipeline/child_1.2_DataCuration.sh

#!/bin/bash

# Loop over each line in samplenames.txt
for s in $(cat $analysis_dir/samplenames.txt)
    do
        samplename="$s"
        echo "Processing sample: $samplename"
        sbatch -J 1.2_DataCuration.$samplename $analysis_dir/biocluster/child_1.2_DataCuration.sh $analysis_dir $folder_nextstep $paired_end $filename_extention $samplename $readname
    done

echo "Finished at `date`"
echo "__DONE__"
