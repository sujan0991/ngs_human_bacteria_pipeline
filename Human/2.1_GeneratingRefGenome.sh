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
# In this script: # 2.1 Generating genome indices files
#                       RNA-STAR version 2.7.9a
#                       Reference genome version: GRCh38.p13, see https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.39/ (RefSeq)
# ******************************************************************************

# ******************************************************************************
#                           !USER DEFINED HEADING!
# ******************************************************************************
# External shell arguments:
# ******************************************************************************
analysis_dir=$1 #/group/bioinf/Data/Madina_test/pipeline/v1/Human
# ******************************************************************************

# Internal arguments:
# *******************************************************************************
threads=4
reference_genome38=/group/bioinf/Data/Madina_test/pipeline/v1/Human/reference_genome/GATK_bundle/hg38/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta #FASTA file with reference genome (mandatory)
reference_gtf=/group/bioinf/Data/Madina_test/pipeline/v1/Human/reference_genome/older/hg38.gtf #GTF file with annptated transcripts (highly recommended)
temp_folder=/home/bioinf/maja488d/star
if [ -d "$temp_folder" ]; then
  rm -r $temp_folder
fi
# ******************************************************************************
# ******************************************************************************
# ******************************************************************************

#*******************************************************************************
# 2.1 Generating genome indexes files
#*******************************************************************************

# Important: generating indices files and mapping have to be done with same version of STAR

# Directory where the indices files will be stored. For human genome we need to have at least 100GB space on disk.
mkdir $analysis_dir/2_mapping_STAR/ -p
mkdir $analysis_dir/2_mapping_STAR/2.1_STAR_genomeDir/ -p
genomeDir=$analysis_dir/2_mapping_STAR/2.1_STAR_genomeDir


echo "Starting step 2.1 Generating genome indexes files"

#STAR manual https://raw.githubusercontent.com/alexdobin/STAR/master/doc/STARmanual.pdf

#   --genomeChrBinNbits: if working with a large (>5000 references), this helps to reduce RAM consumption
#                        recommended min(18,log2[max(GenomeLength/NumberOfReferences,ReadLength)])

#$softwares/STAR --runMode genomeGenerate \
 STAR --runMode genomeGenerate \
    --genomeDir $genomeDir \
    --genomeFastaFiles  $reference_genome38 \
    --sjdbGTFfile $reference_gtf \
    --sjdbOverhang 75 \
    --runThreadN $threads \
    --outTmpDir $temp_folder


echo "Finished step 2.1 Generating genome indexes files"

