#!/bin/bash
set -e
analysis_dir=/group/bioinf_biomarkers_rna/rna_us_analysis 
bam_files=/group/bioinf_biomarkers_rna/rna_us_analysis/3_Picard_processed/3.3_markDuplicates                                        
filename_extention=.bam

echo 'Started'

#chmod u+x /group/bioinf/Data/Madina_test/pipeline/child_GATK_ecoli.sh

#!/bin/bash

# Loop over each line in samplenames.txt

for s in $(cat $bam_files/bam_files.txt)
    do
        samplename="$s"
        echo "Processing sample: $samplename"
        sbatch -J gatk_gvcf.$samplename $analysis_dir/biocluster/child_5.1_5.3_VariantCalling.sh $analysis_dir $bam_files $filename_extention $samplename
    done

echo "Finished at `date`"
echo "__DONE__"
