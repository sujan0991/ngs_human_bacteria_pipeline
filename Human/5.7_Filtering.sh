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
# In this script: Variant Quality Score Recalibration (VQSR)
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
# If cohort is large, e.g. hundreds of unrelated samples first you need to filter variants (VariantFiltration) on ExcessHet
# Smaller cohorts do not need ExcessHet filtering
# ******************************************************************************

# ******************************************************************************
# Variant Quality Score Recalibration (VQSR): build the model
# done separately per variant type (SNPs and Indels)
# https://gatk.broadinstitute.org/hc/en-us/articles/360036510892-VariantRecalibrator
# https://gatk.broadinstitute.org/hc/en-us/articles/360035531112--How-to-Filter-variants-either-with-VQSR-or-by-hard-filtering
# ******************************************************************************
mkdir $analysis_dir/5_VariantCalling/5.7_recalibration -p

# VQSR for Indels
java $memory -jar $softwares_GATK VariantRecalibrator \
    -R $reference_genome38 \
    -V $analysis_dir/5_VariantCalling/5.6_Merging/5.6_MergedVariants.vcf \
    --resource:mills,known=false,training=true,truth=true,prior=12.0  $resource_folder/resources_broad_hg38_v0_Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
    --resource:axiomPoly,known=false,training=true,truth=false,prior=10 $resource_folder/resources_broad_hg38_v0_Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz \
    --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $resource_folder/resources_broad_hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf \
    -an QD -an ReadPosRankSum -an FS -an SOR -an DP \
    -mode INDEL  \
    -O $analysis_dir/5_VariantCalling/5.7_recalibration/INDEL_output.recal \
    --tranches-file $analysis_dir/5_VariantCalling/5.7_recalibration/INDEL_output.tranches \
    --rscript-file $analysis_dir/5_VariantCalling/5.7_recalibration/INDEL_output.plots.R \
    --tmp-dir $TMPDIR > $TMPDIR/INDELS_recalibration.txt 2>&1

mv $TMPDIR/INDELS_recalibration.txt $analysis_dir/5_VariantCalling/5.7_recalibration

# VQSR for SNPs
java $memory -jar $softwares_GATK VariantRecalibrator \
    -R $reference_genome38 \
    -V $analysis_dir/5_VariantCalling/5.6_Merging/5.6_MergedVariants.vcf \
    --resource:hapmap,known=false,training=true,truth=true,prior=15.0  $resource_folder/resources_broad_hg38_v0_hapmap_3.3.hg38.vcf.gz \
    --resource:omni,known=false,training=true,truth=true,prior=12.0  $resource_folder/resources_broad_hg38_v0_1000G_omni2.5.hg38.vcf.gz \
    --resource:1000G,known=false,training=true,truth=false,prior=10.0 $resource_folder/resources_broad_hg38_v0_1000G_phase1.snps.high_confidence.hg38.vcf.gz \
    --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $resource_folder/resources_broad_hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf \
    -an QD -an ReadPosRankSum -an FS -an SOR -an DP \
  	-mode SNP \
    -O $analysis_dir/5_VariantCalling/5.7_recalibration/SNP_output.recal \
    --tranches-file $analysis_dir/5_VariantCalling/5.7_recalibration/SNP_output.tranches \
    --rscript-file $analysis_dir/5_VariantCalling/5.7_recalibration/SNP_output.plots.R \
    --tmp-dir $TMPDIR > $TMPDIR/SNP_recalibration.txt 2>&1

mv $TMPDIR/SNP_recalibration.txt $analysis_dir/5_VariantCalling/5.7_recalibration
# ******************************************************************************
# Variant Quality Score Recalibration (VQSR): apply the model
# done separately per variant type (SNPs and Indels)
# https://gatk.broadinstitute.org/hc/en-us/articles/360036712971-ApplyVQSR
# ******************************************************************************
# ApplyVQSR for INDEL
java $memory -jar $softwares_GATK ApplyVQSR \
    -R $reference_genome38 \
    -V $analysis_dir/5_VariantCalling/5.7_recalibration/5.6_MergedVariants.vcf \
    --truth-sensitivity-filter-level 99.0 \
    --recal-file $analysis_dir/5_VariantCalling/5.7_recalibration/INDEL_output.recal \
    --tranches-file $analysis_dir/5_VariantCalling/5.7_recalibration/INDEL_output.tranches \
    -mode INDEL \
    -O $analysis_dir/5_VariantCalling/5.7_recalibration/indel.recalibrated.vcf \
    --tmp-dir $TMPDIR > $TMPDIR/recalibrated_final.txt 2>&1

mv $TMPDIR/recalibrated.txt $analysis_dir/5_VariantCalling/5.7_recalibration

# ApplyVQSR for SNPs
java $memory -jar $softwares_GATK ApplyVQSR \
    -R $reference_genome38 \
    -V $analysis_dir/5_VariantCalling/5.7_recalibration/indel.recalibrated.vcf \
    --truth-sensitivity-filter-level 99.0 \
    --recal-file $analysis_dir/5_VariantCalling/5.7_recalibration/SNP_output.recal \
    --tranches-file $analysis_dir/5_VariantCalling/5.7_recalibration/SNP_output.tranches \
    -mode SNP \
    -O $analysis_dir/5_VariantCalling/5.7_recalibration/snp.indel.recalibrated.vcf \
    --tmp-dir $TMPDIR > $TMPDIR/snp.indel.recalibrated.txt 2>&1

mv $TMPDIR/snp.indel.recalibrated.txt $analysis_dir/5_VariantCalling/5.7_recalibration


    

