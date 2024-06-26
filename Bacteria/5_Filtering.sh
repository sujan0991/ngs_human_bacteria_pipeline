#!/usr/bin/env bash
set -e

# ******************************************************************************
#                   OVERVIEW
# ******************************************************************************
# I DATA CLEANING
# 1.1 Quality check of the raw data. FastQC and MultiQC
# 1.2 Data Curation
#        - excluding samples (optional)
#        - UMI extraction (optional; usually NOT needed for bacteria)
#        - cutting adapter sequences (optional)
#        - remove low quality reads (optional)
#        - removing poly-sequences (optional)
# 1.3 Quality check of the curated data

# II Alignment 
#       - DNA -> BWA
#       - RNA -> Eukaryotic -> STAR
#       - RNA -> Prokarytic -> BWA
# 2.1 Generating genome indexes files. (Needs to be done only once!) 
# 2.2 - 2.4 Alignment -> .sam files

# III Picard Processing
# 3.1 Sorting
# 3.2 Adding read group information
# 3.3 Marking duplicates and indexing

# IV Quality check of aligned bam files
# 4.1. Log.final.out to look at sumary stats of mapping 
#       mapping rate. >=75% is a good quality sample
#       multimapping. If multimapping is high, there are usually two issues: 1) sample was contaminated (-> check .fastq data quality) 2) alignment did not work (-> check mapping options)

# 4.2. Qualimap
# https://hbctraining.github.io/Intro-to-rnaseq-hpc-salmon/lessons/03_QC_STAR_and_Qualimap_run.html#qualimap

# V GATK Variant Calling
# 5.1. Call variants per sample - HaplotypeCaler
# 5.2. Consolidate gvcf-s into a single .vcf file
# 5.3. Genotyping 
# 5.4. Split SNPs and Indels
# ******************************************************************************

# ******************************************************************************
# In this script: STEP V. GATK Variant Calling: Germline short variant discovery (SNPs + Indels)
#                         
#                         https://gatk.broadinstitute.org/hc/en-us/articles/360035535932-Germline-short-variant-discovery-SNPs-Indels- 
# GATK Version: 4.2.5.0                                                                                                                              
# ******************************************************************************


# External shell arguments:
# ******************************************************************************
analysis_dir=$1 #/group/bioinf/Data/Madina_test/pipeline/v1 main analysis folder where all performed  steps of the pipeline are saved
input_vcf=$2 # $analysis_dir/5_VariantCalling/5.3_Genotyping/variants.vcf
# ******************************************************************************
# Internal arguments:
# *******************************************************************************
reference_files=/group/bioinf/Data/Madina_test/bacteria/reference_genomes/
reference_genome=$reference_files/reference.fna
memory=-Xmx40g
ploidy=1
maximum_nocall_fraction=0.15 # 0.2
sw=/group/bioinf/Users/Madina/sw/
softwares_GATK=/group/bioinf/Users/Madina/sw/gatk-4.2.5.0/gatk-package-4.2.5.0-local.jar # path to GATK .jar file
# ******************************************************************************
# ******************************************************************************
# ******************************************************************************

echo "Started Analysing at `date`"

# ******************************************************************************
# 5.4. SelectVariants: separating Indels from SNP
# ******************************************************************************
mkdir $analysis_dir/5_VariantCalling/5.4_SelectVariants/5.4.1_SNPs -p
mkdir $analysis_dir/5_VariantCalling/5.4_SelectVariants/5.4.2_Indels -p

java $memory -jar $softwares_GATK SelectVariants \
    --reference $reference_genome \
    --variant $input_vcf \
    --select-type-to-include SNP \
    --output $analysis_dir/5_VariantCalling/5.4_SelectVariants/5.4.1_SNPs/SNPs.vcf \
    --tmp-dir $TMPDIR > $TMPDIR/5.4.1_SNPs.txt 2>&1
    
mv $TMPDIR/5.4.1_SNPs.txt $analysis_dir/5_VariantCalling/5.4_SelectVariants/5.4.1_SNPs

java $memory -jar $softwares_GATK SelectVariants \
    --reference $reference_genome \
    --variant $input_vcf \
    --select-type-to-include INDEL \
    --output $analysis_dir/5_VariantCalling/5.4_SelectVariants/5.4.2_Indels/Indels.vcf \
    --tmp-dir $TMPDIR > $TMPDIR/5.4.2_Indels.txt 2>&1 
    
mv $TMPDIR/5.4.2_Indels.txt $analysis_dir/5_VariantCalling/5.4_SelectVariants/5.4.2_Indels

# ******************************************************************************
# 5.5. VariantFiltration: you might want to reconsider the filters based on your data!
# ******************************************************************************
mkdir $analysis_dir/5_VariantCalling/5.5_HardFiltering/_SNPs -p
mkdir $analysis_dir/5_VariantCalling/5.5_HardFiltering/_Indels -p

java $memory -jar $softwares_GATK VariantFiltration \
    --reference $reference_genome \
    --variant $analysis_dir/5_VariantCalling/5.4_SelectVariants/5.4.1_SNPs/SNPs.vcf \
    --filter-expression "QD < 2.0" --filter-name "QD2" \
    --filter-expression "QUAL < 30.0" --filter-name "QUAL30" \
    --filter-expression "FS > 60.0" --filter-name "FS60" \
    --output $analysis_dir/5_VariantCalling/5.5_HardFiltering/_SNPs/5.5_snps_filtered.vcf \
    --tmp-dir $TMPDIR > $TMPDIR/5.5_snps_filtered.txt 2>&1

mv $TMPDIR/5.5_snps_filtered.txt $analysis_dir/5_VariantCalling/5.5_HardFiltering/_SNPs

java $memory -jar $softwares_GATK VariantFiltration \
    --reference $reference_genome \
    --variant $analysis_dir/5_VariantCalling/5.4_SelectVariants/5.4.2_Indels/Indels.vcf \
    --filter-expression "QD < 2.0" --filter-name "QD2" \
    --filter-expression "QUAL < 30.0" --filter-name "QUAL30" \
    --filter-expression "FS > 200.0" --filter-name "FS200" \
    --output $analysis_dir/5_VariantCalling/5.5_HardFiltering/_Indels/5.5_indels_filtered.vcf \
    --tmp-dir $TMPDIR > $TMPDIR/5.5_indels_filtered.txt 2>&1 
    
mv $TMPDIR/5.5_indels_filtered.txt $analysis_dir/5_VariantCalling/5.5_HardFiltering/_Indels

# ******************************************************************************
# 5.6 Combining the filtered SNPs and INDELs with MergeVcfs
# This is now called as a picard function
# ******************************************************************************
mkdir $analysis_dir/5_VariantCalling/5.6_filtered -p

singularity exec -B /group /group/bioinf/Software/singularity/bacteria.sif java  $memory -Djava.io.tmpdir=$analysis_dir/5_VariantCalling/5.6_filtered/temp -jar $sw/picard.jar MergeVcfs \
    -I $analysis_dir/5_VariantCalling/5.5_HardFiltering/_SNPs/5.5_snps_filtered.vcf \
    -I $analysis_dir/5_VariantCalling/5.5_HardFiltering/_Indels/5.5_indels_filtered.vcf \
    -O $analysis_dir/5_VariantCalling/5.6_filtered/5.6_combined_filtered.vcf \
    > $analysis_dir/5_VariantCalling/5.6_filtered/5.6_combined_filtered.txt 2>&1 

# ******************************************************************************
# 5.7 Sorting the combined file with SortVcf
# This is now called as a picard function
# ******************************************************************************
mkdir $analysis_dir/5_VariantCalling/5.7_sorted -p

singularity exec -B /group /group/bioinf/Software/singularity/bacteria.sif java  $memory -Djava.io.tmpdir=$analysis_dir/5_VariantCalling/5.7_sorted/temp -jar $sw/picard.jar SortVcf \
    -I $analysis_dir/5_VariantCalling/5.6_filtered/5.6_combined_filtered.vcf \
    -O $analysis_dir/5_VariantCalling/5.7_sorted/5.7_sorted_combined_filtered.vcf \
    > $analysis_dir/5_VariantCalling/5.7_sorted/5.7_sorted_combined_filtered.txt 2>&1 

# ******************************************************************************
# 5.8 Applying the hard filteres with SelectVariants
# ******************************************************************************
mkdir $analysis_dir/5_VariantCalling/5.8_sorted_combined_filtered -p 

java $memory -jar $softwares_GATK SelectVariants \
    --reference $reference_genome \
    --exclude-filtered \
    --variant $analysis_dir/5_VariantCalling/5.7_sorted/5.7_sorted_combined_filtered.vcf \
    --output $analysis_dir/5_VariantCalling/5.8_sorted_combined_filtered/5.8_selected_sorted_combined_filtered.vcf \
    --tmp-dir $TMPDIR > $TMPDIR/5.8_selected_sorted_combined_filtered.txt 2>&1 
    
mv $TMPDIR/5.8_selected_sorted_combined_filtered.txt $analysis_dir/5_VariantCalling/5.8_sorted_combined_filtered


# ******************************************************************************
# 5.9 Filter low quality genotypes with VariantFiltration
# 5.9.1 Filtering based on low quality GQ with VariantFiltration
# ******************************************************************************
mkdir $analysis_dir/5_VariantCalling/5.9_ -p
mkdir $analysis_dir/5_VariantCalling/5.9_/5.9.1 -p
mkdir $analysis_dir/5_VariantCalling/5.9_/5.9.2.1 -p
mkdir $analysis_dir/5_VariantCalling/5.9_/5.9.2.2 -p
 
java $memory -jar $softwares_GATK VariantFiltration \
    --reference $reference_genome \
    --variant $analysis_dir/5_VariantCalling/5.8_sorted_combined_filtered/5.8_selected_sorted_combined_filtered.vcf \
    --filter-expression "GQ < 20" --filter-name "GQ20" \
    --output $analysis_dir/5_VariantCalling/5.9_/5.9.1/5.9.1_filtered_GQ20.vcf \
    --tmp-dir $TMPDIR > $TMPDIR/5.9.1_filtered_GQ20.txt 2>&1 
    
mv $TMPDIR/5.9.1_filtered_GQ20.txt $analysis_dir/5_VariantCalling/5.9_/5.9.1

# ******************************************************************************
# 5.9.2 Applying the low quality genotype filtere with SelectVariants
# 5.9.2.1 Changing the genotypes of the filtered variants to null or no call(./.)
# https://gatkforums.broadinstitute.org/gatk/discussion/2806/howto-apply-hard-filters-to-a-call-set
# https://gatkforums.broadinstitute.org/gatk/discussion/12350/
# https://gatkforums.broadinstitute.org/gatk/discussion/2806/howto-apply-hard-filters-to-a-call-set/p2
# ******************************************************************************
java $memory -jar $softwares_GATK SelectVariants \
    --reference $reference_genome \
    --set-filtered-gt-to-nocall \
    --variant $analysis_dir/5_VariantCalling/5.9_/5.9.1/5.9.1_filtered_GQ20.vcf \
    --output  $analysis_dir/5_VariantCalling/5.9_/5.9.2.1/5.9.2.1_filtered_GQ20_to_nocall_.vcf \
    --tmp-dir $TMPDIR > $TMPDIR/5.9.2.1_filtered_GQ20_to_nocall_.txt 2>&1 
    
mv $TMPDIR/5.9.2.1_filtered_GQ20_to_nocall_.txt $analysis_dir/5_VariantCalling/5.9_/5.9.2.1

# ******************************************************************************
# 5.9.2.2 Keeping only variants that are no call (./.) in maximum $maximum_nocall_fraction of samples.
# https://gatkforums.broadinstitute.org/gatk/discussion/9795/how-to-use-setfilteredgttonocall-and-maxnocallfraction-with-selectvariants
# ******************************************************************************
java $memory -jar $softwares_GATK SelectVariants \
    --reference $reference_genome \
    --max-nocall-fraction $maximum_nocall_fraction \
    --variant $analysis_dir/5_VariantCalling/5.9_/5.9.2.1/5.9.2.1_filtered_GQ20_to_nocall_.vcf \
    --output $analysis_dir/5_VariantCalling/5.9_/5.9.2.2/5.9.2.2_selected__filtered_GQ20_max_no_call_fraction.vcf \
    --tmp-dir $TMPDIR > $TMPDIR/5.9.2.2_selected__filtered_GQ20_max_no_call_fraction.txt 2>&1 
    
mv $TMPDIR/5.9.2.2_selected__filtered_GQ20_max_no_call_fraction.txt $analysis_dir/5_VariantCalling/5.9_/5.9.2.2

    
