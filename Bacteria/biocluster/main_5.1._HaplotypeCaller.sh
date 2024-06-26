#!/bin/bash
analysis_dir=/group/bioinf/Users/Madina/pipeline_tested/bacteria 
bam_files=/group/bioinf/Users/Madina/pipeline_tested/bacteria/Picard_processed/3.3_markDuplicates                                         
filename_extention=.bam

echo 'Started'

#chmod u+x /group/bioinf/Data/Madina_test/pipeline/child_GATK_ecoli.sh

#!/bin/bash

# Loop over each line in samplenames.txt
for s in $(cat $bam_files/bam_files.txt)
    do
        samplename="$s"
        echo "Processing sample: $samplename"
        sbatch -J $samplename.HaplotypeCaller $analysis_dir/biocluster/child_5.1._HaplotypeCaller.sh $analysis_dir $bam_files $filename_extention $samplename
    done

echo "Finished at `date`"
echo "__DONE__"

