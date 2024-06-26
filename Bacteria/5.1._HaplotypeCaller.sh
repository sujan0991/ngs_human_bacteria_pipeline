#!/usr/bin/env bash
set -e 

# ******************************************************************************
#                   OVERVIEW
# ******************************************************************************
# I Preprocessing
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
# 4.1. If alignment was done with STAR, Log.final.out to look at sumary stats of mapping 
#       mapping rate. >=75% is a good quality sample
#       multimapping. If multimapping is high, there are usually two issues: a) sample was contaminated (-> check .fastq data quality) b) alignment did not work (-> check mapping options)

# 4.2. Qualimap
# https://hbctraining.github.io/Intro-to-rnaseq-hpc-salmon/lessons/03_QC_STAR_and_Qualimap_run.html#qualimap

# V Variant Calling
# 5.1 HaplotypeCaller - .makes g.vcf file per each sample
# 5.2 ConsolidateGVCFs - combine g.vcf files into a single file
# 5.3 Genotyping - makes single .vcf file with all variants
# 5.4 Select variants - splits variants into SNPs and indels
# 5.5 - 5.9 Filtering, sorting, and combining
# ******************************************************************************

# ******************************************************************************
# In this script: STEP V. GATK Variant Calling: Germline short variant discovery (SNPs + Indels)
#                         5.1. HaplotypeCaller: per sample .g.vcf file
#                         https://gatk.broadinstitute.org/hc/en-us/articles/360035535932-Germline-short-variant-discovery-SNPs-Indels- 
# GATK Version: 4.2.5.0                                                                                                                              
# ******************************************************************************


# External shell arguments:
# ******************************************************************************
analysis_dir=$1 #/group/bioinf/Data/Madina_test/pipeline/v1 main analysis folder where all performed  steps of the pipeline are saved
bam_files=$2 # /group/bioinf/Data/Madina_test/pipeline/v1/Picard_processed/3.3_markDuplicates/                                          
filename_extention=$3 # .bam
samplename=$4 # name of the sample, e.g. DNA1
# ******************************************************************************
filename="${samplename}${filename_extention}"

# Internal arguments:
# *******************************************************************************
reference_files=/group/bioinf/Data/Madina_test/bacteria/reference_genomes/
reference_genome=$reference_files/reference.fna
memory=-Xmx40g
ploidy=1
format=GVCF
softwares_GATK=/group/bioinf/Users/Madina/sw/gatk-4.2.5.0/gatk-package-4.2.5.0-local.jar # path to GATK .jar file
# ******************************************************************************
# ******************************************************************************
# ******************************************************************************

echo "Started HaplotypeCaller at `date`"
# ******************************************************************************
# 5.1. Haplotypecaller: begin by calling variants per sample in order to produce a file in GVCF format      
# ******************************************************************************

mkdir $analysis_dir/5_VariantCalling/5.1_HaplotypeCaller/ -p

echo 'Step 5.1. Variant Calling: HaplotypeCaller' 
    
java $memory -jar $softwares_GATK HaplotypeCaller \
    -R $reference_genome \
    -I $bam_files/$filename \
    -ploidy $ploidy \
    --emit-ref-confidence $format \
    -O $analysis_dir/5_VariantCalling/5.1_HaplotypeCaller/$samplename.g.vcf \
    --tmp-dir $TMPDIR > $TMPDIR/$samplename.HaplotypeCaller.txt 2>&1

mv $TMPDIR/$samplename.HaplotypeCaller.txt $analysis_dir/5_VariantCalling/5.1_HaplotypeCaller




