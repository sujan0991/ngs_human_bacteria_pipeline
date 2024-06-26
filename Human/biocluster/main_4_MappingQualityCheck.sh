#!/bin/bash
analysis_dir=/group/bioinf/Users/Madina/pipeline_tested/human
bam_files=/group/bioinf/Users/Madina/pipeline_tested/human/3_Picard_processed/3.3_markDuplicates
filename_extention=.bam

echo 'Started'

#chmod u+x /group/bioinf/Data/Madina_test/pipeline/child_III_Picard.sh

#!/bin/bash

find $bam_files -type f -name "*.bam" -exec basename {} \; | sed 's/\.[^.]*$//' > $bam_files/bam_files.txt
# Loop over each line in samplenames.txt
for s in $(cat $bam_files/bam_files.txt)
    do
        samplename="$s"
        echo "Processing sample: $samplename"
        sbatch -J $samplename.qualimap $analysis_dir/biocluster/child_4_MappingQualityCheck.sh $analysis_dir $bam_files $filename_extention $samplename
    done

echo "Finished at `date`"
echo "__DONE__"

