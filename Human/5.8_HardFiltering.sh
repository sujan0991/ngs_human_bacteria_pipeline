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

# II Alignment -> .sam file
# 2.1 Generating genome indexes files. (Needs to be done only once!) 
# 2.2 - 2.4 Alignment 

# III Picard Processing -> .bam file
# 3.1 Sorting
# 3.2 Adding read group information
# 3.3 Marking duplicates and indexing

# IV Quality check of aligned bam files
# 4.1 Log.final.out to look at sumary stats of mapping 
#       mapping rate. >=75% is a good quality sample
#       multimapping. If multimapping is high, there are usually two issues: 1) sample was contaminated (-> check .fastq data quality) 2) alignment did not work (-> check mapping options)

# 4.2 Qualimap

# V Variant Calling
# 5.1 SplitNCigarReads
# 5.2.a BaseRecalibrator
# 5.2.b ApplyBQSR
# 5.3 HaplotypeCaller
# 5.4 GenomicsDBImport
# 5.5 GenotypeGVCFs
# 5.6 Merging per/chromosome VCFs into one
# 5.7 VQSR  
# 5.8 Hard filtering 
# ******************************************************************************

# ******************************************************************************
# In this script: Hard filtering
#  
# # https://gatk.broadinstitute.org/hc/en-us/articles/360035531112--How-to-Filter-variants-either-with-VQSR-or-by-hard-filtering                       
# https://gatk.broadinstitute.org/hc/en-us/articles/360035535932-Germline-short-variant-discovery-SNPs-Indels- 
# GATK Version: 4.2.5.0                                                                                                                              
# ******************************************************************************


# External shell arguments:
# ******************************************************************************
analysis_dir=$1 #/group/bioinf/Data/Madina_test/pipeline/v1 main analysis folder where all performed  steps of the pipeline are saved
# ******************************************************************************
# Internal arguments:
# *******************************************************************************
nullSites=0.2
reference_genome38=/group/bioinf/Data/Madina_test/pipeline/v1/Human/reference_genome/GATK_bundle/hg38/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta
resource_folder=/group/bioinf/Data/Madina_test/pipeline/v1/Human/reference_genome/GATK_bundle/hg38
memory=-Xmx30g 
sw=/group/bioinf/Users/Madina/sw/
softwares_GATK=/group/bioinf/Users/Madina/sw/gatk-4.2.5.0/gatk-package-4.2.5.0-local.jar # path to GATK .jar file
# ******************************************************************************
# ******************************************************************************

# Links

# population resource files are available at: https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/v0;tab=objects?prefix=&forceOnObjectsSortingFiltering=false

#   Joint germline variant calling
#   https://gatk.broadinstitute.org/hc/en-us/articles/360035890431-The-logic-of-joint-calling-for-germline-short-variants
#   https://gatk.broadinstitute.org/hc/en-us/articles/360035890411

#   VQSR
#   https://gatk.broadinstitute.org/hc/en-us/articles/360035531612?id=1259
#   https://gatk.broadinstitute.org/hc/en-us/articles/360035531112--How-to-Filter-variants-either-with-VQSR-or-by-hard-filtering

#   Hard filetring 
#   https://gatk.broadinstitute.org/hc/en-us/articles/360035890471
# ******************************************************************************

echo "Started Analysing at `date`"
# ******************************************************************************
#  Hard filtering
# ******************************************************************************
mkdir $analysis_dir/5_VariantCalling/5.8_filtering -p

# ******************************************************************************
# Hard filter SNPs and indels separately
# ******************************************************************************

# Subset to SNPs only 
java $memory -jar $softwares_GATK SelectVariants \
    --reference $reference_genome38 \
    --variant $analysis_dir/5_VariantCalling/5.7_recalibration/snp.indel.recalibrated.vcf \
    --select-type-to-include SNP \
    --output $analysis_dir/5_VariantCalling/5.8_filtering/SNPs.vcf \
    --tmp-dir $TMPDIR > $TMPDIR/5.8_SNPs.txt 2>&1

mv $TMPDIR/5.8_SNPs.txt $analysis_dir/5_VariantCalling/5.8_filtering

# Subset to indels
# to not lose mixed records set the variant type to "MIXED"
# note: MIXED records are treated as indels
java $memory -jar $softwares_GATK SelectVariants \
    --reference $reference_genome38 \
    --variant $analysis_dir/5_VariantCalling/5.7_recalibration/snp.indel.recalibrated.vcf \
    --select-type-to-include MIXED \
    --output $analysis_dir/5_VariantCalling/5.8_filtering/Indels.vcf \
    --tmp-dir $TMPDIR > $TMPDIR/5.8_Indels.txt 2>&1

mv $TMPDIR/5.8_Indels.txt $analysis_dir/5_VariantCalling/5.8_filtering

# ******************************************************************************
# Hard filtering SNPs
# https://gatk.broadinstitute.org/hc/en-us/articles/360036350452-VariantFiltration
# ******************************************************************************

# Filter low quality on SNPs
java $memory -jar $softwares_GATK VariantFiltration \
    -R $reference_genome38 \
    -V $analysis_dir/5_VariantCalling/5.8_filtering/SNPs.vcf \
    -G-filter "GQ < 20.0" -G-filter-name lowGQ  \
    -O $analysis_dir/5_VariantCalling/5.8_filtering/Gfiltered_SNPs.vcf \
    --tmp-dir $TMPDIR > $TMPDIR/Gfiltered_SNPs.txt 2>&1

mv $TMPDIR/Gfiltered_SNPs.txt $analysis_dir/5_VariantCalling/5.8_filtering

# SelectVariants1 on SNPs
java $memory -jar $softwares_GATK SelectVariants \
    -R $reference_genome38 \
    --exclude-filtered \
    --exclude-non-variants \
    --set-filtered-gt-to-nocall  \
    -V $analysis_dir/5_VariantCalling/5.8_filtering/Gfiltered_SNPs.vcf \
    -O $analysis_dir/5_VariantCalling/5.8_filtering/SelectVariants1_SNPs.vcf \
    --tmp-dir $TMPDIR > $TMPDIR/SelectVariants1_SNPs.txt 2>&1

mv $TMPDIR/SelectVariants1_SNPs.txt $analysis_dir/5_VariantCalling/5.8_filtering 

# SelectVariants2 on SNPs
java $memory -jar $softwares_GATK SelectVariants \
    -R $reference_genome38 \
    --max-nocall-fraction $nullSites \
    -V $analysis_dir/5_VariantCalling/5.8_filtering/SelectVariants1_SNPs.vcf \
    -O $analysis_dir/5_VariantCalling/5.8_filtering/SelectVariants2_SNPs_.$nullSites.vcf \
    --tmp-dir $TMPDIR > $TMPDIR/SelectVariants2_SNPs.txt 2>&1

mv $TMPDIR/SelectVariants2_SNPs.txt $analysis_dir/5_VariantCalling/5.8_filtering 

# ******************************************************************************
# Hard filtering Indels
# https://gatk.broadinstitute.org/hc/en-us/articles/360036350452-VariantFiltration
# ******************************************************************************
# Filter low quality on Indels
java $memory -jar $softwares_GATK VariantFiltration \
    -R $reference_genome38 \
    -V $analysis_dir/5_VariantCalling/5.8_filtering/Indels.vcf \
    -G-filter "GQ < 20.0" -G-filter-name lowGQ  \
    -O $analysis_dir/5_VariantCalling/5.8_filtering/Gfiltered_Indels.vcf \
    --tmp-dir $TMPDIR > $TMPDIR/Gfiltered_Indels.txt 2>&1

mv $TMPDIR/Gfiltered_Indels.txt $analysis_dir/5_VariantCalling/5.8_filtering

# SelectVariants1 on Indels
java $memory -jar $softwares_GATK SelectVariants \
    -R $reference_genome38 \
    --exclude-filtered \
    --exclude-non-variants \
    --set-filtered-gt-to-nocall  \
    -V $analysis_dir/5_VariantCalling/5.8_filtering/Gfiltered_Indels.vcf \
    -O $analysis_dir/5_VariantCalling/5.8_filtering/SelectVariants1_Indels.vcf \
    --tmp-dir $TMPDIR > $TMPDIR/SelectVariants1_Indels.txt 2>&1

mv $TMPDIR/SelectVariants1_Indels.txt $analysis_dir/5_VariantCalling/5.8_filtering 

# SelectVariants2 on Indels
java $memory -jar $softwares_GATK SelectVariants \
    -R $reference_genome38 \
    --max-nocall-fraction $nullSites \
    -V $analysis_dir/5_VariantCalling/5.8_filtering/SelectVariants1_Indels.vcf \
    -O $analysis_dir/5_VariantCalling/5.8_filtering/SelectVariants2_Indels_.$nullSites.vcf \
    --tmp-dir $TMPDIR > $TMPDIR/SelectVariants2_Indels.txt 2>&1

mv $TMPDIR/SelectVariants2_Indels.txt $analysis_dir/5_VariantCalling/5.8_filtering 

# ******************************************************************************
# Gather SNPs and Indels into one VCF file
# ******************************************************************************
singularity exec -B /group /group/bioinf/Software/singularity/bacteria.sif java  $memory -Djava.io.tmpdir=$analysis_dir/5_VariantCalling/5.8_filtering -jar $sw/picard.jar MergeVcfs \
    -I $analysis_dir/5_VariantCalling/5.8_filtering/SelectVariants2_SNPs_.$nullSites.vcf \
    -I $analysis_dir/5_VariantCalling/5.8_filtering/SelectVariants2_Indels_.$nullSites.vcf \
    -O $analysis_dir/5_VariantCalling/5.8_filtering/filtered.$nullSites.vcf \
    > $analysis_dir/5_VariantCalling/5.8_filtering/filtered.$nullSites.txt 2>&1 
    
# ******************************************************************************
# Evaluate the filtered files
# ******************************************************************************
singularity exec -B /group /group/bioinf/Software/singularity/bacteria.sif java  $memory -Djava.io.tmpdir=$analysis_dir/5_VariantCalling/5.8_filtering -jar $sw/picard.jar CollectVariantCallingMetrics \
    -I $analysis_dir/5_VariantCalling/5.8_filtering/filtered.$nullSites.vcf \
    -O $analysis_dir/5_VariantCalling/5.8_filtering/metrics \
    --DBSNP $resource_folder/resources_broad_hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf \
    > $analysis_dir/5_VariantCalling/5.8_filtering/filtering.metrics.txt 2>&1 





    

