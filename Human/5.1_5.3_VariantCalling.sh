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
# In this script: 5.1-5.3 Final output - per sample GVCF file
# ******************************************************************************

# The pipeline follows recommendations regarding GATK variant calling from here:
# https://gatk.broadinstitute.org/hc/en-us/articles/360035531192-RNAseq-short-variant-discovery-SNPs-Indels-
#  https://hpc.nih.gov/training/gatk_tutorial/

# For internal arguments you need to provide a reference genome file in the FASTA format and files with known sites 
# The known sites are used to mask out positions with known variation to avoid confusing real variation with errors. 
# Please make sure that alongside the vcf files with known sites you also save a .tbi file for the respective .vcf file
# Either dowload a .tbi file for the .vcf file or index yourself
# I used compatible files for germline variant calling from GATK bundle:
# https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/v0;tab=objects?prefix=&forceOnObjectsSortingFiltering=false

# ******************************************************************************
# External shell arguments:
# ******************************************************************************
analysis_dir=$1 #/group/bioinf/Data/Madina_test/pipeline/v1/Human main analysis folder where all performed  steps of the pipeline are saved
bam_files=$2 # /group/bioinf/Data/Madina_test/pipeline/v1/Human/Picard_processed/3.3_markDuplicates/                                          
filename_extention=$3 # .bam
samplename=$4 # name of the sample, e.g. DNA1
# ******************************************************************************
filename="${samplename}${filename_extention}" 

# Internal arguments:
# *******************************************************************************
reference_genome38=/group/bioinf/Data/Madina_test/pipeline/v1/Human/reference_genome/GATK_bundle/hg38/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta
known_dbsnp=/group/bioinf/Data/Madina_test/pipeline/v1/Human/reference_genome/GATK_bundle/hg38/resources_broad_hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf 
known_indels=/group/bioinf/Data/Madina_test/pipeline/v1/Human/reference_genome/GATK_bundle/hg38/resources_broad_hg38_v0_Homo_sapiens_assembly38.known_indels.vcf.gz
known_mills=/group/bioinf/Data/Madina_test/pipeline/v1/Human/reference_genome/GATK_bundle/hg38/resources_broad_hg38_v0_Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
memory=-Xmx30g 
ploidy=2 
softwares_GATK=/group/bioinf/Users/Madina/sw/gatk-4.2.5.0/gatk-package-4.2.5.0-local.jar # path to GATK .jar file
# ******************************************************************************
# ******************************************************************************
# ******************************************************************************

# Generate files if needed
if [ ! -s $reference_genome38.fai ]; then
    echo "Generating .fai file for reference"
    samtools faidx $reference_genome38
fi

filename38=$(echo $reference_genome38 | rev | cut -d"." -f2-  | rev)

if [ ! -s $filename38.dict ]; then
    echo "Generating .dict file for reference"
  java -jar $softwares_GATK CreateSequenceDictionary -R $reference_genome38
fi

# Directories for GATK output
mkdir $analysis_dir/5_VariantCalling/ -p
mkdir $analysis_dir/5_VariantCalling/5.1_SplitsCigar/ -p
mkdir $analysis_dir/5_VariantCalling/5.2_recal_data_table/ -p
mkdir $analysis_dir/5_VariantCalling/5.2_recal_reads/ -p
mkdir $analysis_dir/5_VariantCalling/5.3_BP_RNA_files/ -p

# ******************************************************************************
# SplitNCigarReads
#           - must be done before variant calling for RNAseq
#           - prevents introns from being icnluded in the alignment
# ******************************************************************************
echo "Step 5.1. Splitting spliced reads: SplitNCigarReads"
java $memory -jar $softwares_GATK SplitNCigarReads \
   -I $bam_files/$filename \
   -O $analysis_dir/5_VariantCalling/5.1_SplitsCigar/$filename \
   -R $reference_genome38 \
   --tmp-dir $TMPDIR > $TMPDIR/$samplename.SplitNCigarReads.txt 2>&1

mv $TMPDIR/$samplename.SplitNCigarReads.txt $analysis_dir/5_VariantCalling/5.1_SplitsCigar/
# ******************************************************************************
# Base Recalibration
#           -recalibrates base quality scores 
#           -initial statistics collection
#           -per-region statistics collection
#           -generation of recalibrated dataset
# ******************************************************************************
# 4.2.a Recalibrate base quality scores 
echo "Step 5.2.a. Base Recalibration" 
java $memory -jar $softwares_GATK BaseRecalibrator \
   -I $analysis_dir/5_VariantCalling/5.1_SplitsCigar/$filename \
   -R $reference_genome38 \
   -O $analysis_dir/5_VariantCalling/5.2_recal_data_table/$samplename.recall_data_table \
   --known-sites $known_dbsnp \
   --known-sites $known_indels \
   --known-sites $known_mills \
   --tmp-dir $TMPDIR > $TMPDIR/$samplename.BaseRecalibrator.txt 2>&1

mv $TMPDIR/$samplename.BaseRecalibrator.txt $analysis_dir/5_VariantCalling/5.2_recal_data_table/ 

# 4.2.b Apply recalibrated base quality scores 
echo "Step 5.2.b. ApplyBQSR"
java $memory -jar $softwares_GATK ApplyBQSR \
    -I $analysis_dir/5_VariantCalling/5.1_SplitsCigar/$filename \
    -R $reference_genome38 \
    -bqsr $analysis_dir/5_VariantCalling/5.2_recal_data_table/$samplename.recall_data_table \
    -O $analysis_dir/5_VariantCalling/5.2_recal_reads/$filename \
   --tmp-dir $TMPDIR > $TMPDIR/$samplename.ApplyBQSR.txt 2>&1

mv $TMPDIR/$samplename.ApplyBQSR.txt  $analysis_dir/5_VariantCalling/5.2_recal_reads/

# ******************************************************************************
# Variant Calling: HaplotypeCaller
#           -call potential variant sites per sample 
#           -our pipeline produces a GVCF file per sample
# ******************************************************************************
echo "Step 5.3. Variant Calling: HaplotypeCaller"    
java $memory -jar $softwares_GATK HaplotypeCaller \
    -R $reference_genome38 \
    -I $analysis_dir/5_VariantCalling/5.2_recal_reads/$filename \
    -ploidy $ploidy \
    --emit-ref-confidence GVCF --dont-use-soft-clipped-bases True \
    -O $analysis_dir/5_VariantCalling/5.3_BP_RNA_files/$samplename.g.vcf \
   --tmp-dir $TMPDIR > $TMPDIR/$samplename.HaplotypeCaller.txt 2>&1

mv $TMPDIR/$samplename.HaplotypeCaller.txt  $analysis_dir/5_VariantCalling/5.3_BP_RNA_files/




