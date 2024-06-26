#!/bin/bash
analysis_dir=/group/bioinf/Users/Madina/pipeline_tested/bacteria
sam_files=$analysis_dir/2._AlignmentBacteria
filename_extention=.sam

echo 'Started'

#chmod u+x /group/bioinf/Data/Madina_test/pipeline/child_III_Picard.sh

#!/bin/bash

# Loop over each line in samplenames.txt
for s in $(cat $analysis_dir/samplenames.txt)
    do
        samplename="$s"
        echo "Processing sample: $samplename"
        sbatch -J 3_Picard.$samplename $analysis_dir/biocluster/child_3_Picard.sh $analysis_dir $sam_files $filename_extention $samplename
    done

echo "Finished at `date`"
echo "__DONE__"

