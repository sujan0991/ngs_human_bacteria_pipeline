#!/bin/bash
set -e 
analysis_dir=/group/bioinf/Users/Madina/pipeline_tested/human
sam_files=/group/bioinf/Users/Madina/pipeline_tested/human/2_mapping_STAR/2.2_Mapping/star_2pass/
filename_extention=.sam

echo 'Started'

#chmod u+x /group/bioinf/Data/Madina_test/pipeline/child_III_Picard.sh

#!/bin/bash

# Loop over each line in samplenames.txt
for s in $(cat $analysis_dir/sam_files.txt)
    do
        samplename="$s"
        echo "Processing sample: $samplename"
        sbatch -J $samplename.picard $analysis_dir/biocluster/child_3_Picard.sh $analysis_dir $sam_files $filename_extention $samplename
    done

echo "Finished at `date`"
echo "__DONE__"

