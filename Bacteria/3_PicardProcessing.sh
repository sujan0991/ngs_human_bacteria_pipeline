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
# In this script: STEP 3. Picard Processing
#       - Adding read group information
#       - Sorting
#       - Marking duplicates and indexing
# ******************************************************************************

# ******************************************************************************
#                           !USER PROVIDED HEADING!
# ******************************************************************************

# External shell arguments:
# ******************************************************************************
analysis_dir=$1 #/group/bioinf/Data/Madina_test/pipeline/v1
sam_files=$2
filename_extention=$3 #.sam
samplename=$4 #DNA1
# ******************************************************************************

# Constructing the name of the file 
filename="${samplename}${filename_extention}"

# Internal arguments:
# ******************************************************************************
sw=/group/bioinf/Users/Madina/sw/
memory=-Xmx30g
# ******************************************************************************
# ******************************************************************************

# Picard tools require Java 1.17 or newer
# we store the Picard software in a singularity container
# singularity exec -B /group /group/bioinf/Software/singularity/bacteria.sif java
# -Xmx60g indicating the amount of memory to allocate to the Java virtual machine which will run the program. 60 gigabytes

# analysis_dir --- full path to the folder where all analysis files will be kept. They will be divided into different subfolders at each each step.
# sw --- folder with java based softwares, e.g. fastqc, picard, cutadapt
# sam_files -- path to where the aligned sam files were saved in the Mapping step
# ******************************************************************************

# ******************************************************************************
# LATEST VERSION PICARD UPDATES: https://github.com/broadinstitute/picard/wiki
#       - PICARD installation 
#           - Picard 3.0 requires Java 1.17 version
#           - Picard is now built using gradle

#       -PICARD syntax
#       -NOTE: Picard's command line syntax has changed
#       -wherever you used to do e.g. I=input.bam, you now do -I input.bam.
#       -see https://github.com/broadinstitute/picard/wiki/Command-Line-Syntax-Transition-For-Users-(Pre-Transition)
# ******************************************************************************

mkdir $analysis_dir/Picard_processed/3.1_sortSam/ -p
mkdir $analysis_dir/Picard_processed/3.2_AddOrReplaceReadGroups -p
mkdir $analysis_dir/Picard_processed/3.3_markDuplicates -p
mkdir $analysis_dir/Picard_processed/3.3_markDuplicates/deduppedBamSortedFiles -p
mkdir $analysis_dir/Picard_processed/temp -p

#       Note: it is important to explicitly set -CREATE_INDEX true in 3.3 MarkDuplicates. This is necessary for future variant calling.

     singularity exec -B /group /group/bioinf/Software/singularity/bacteria.sif java  $memory -Djava.io.tmpdir=$analysis_dir/Picard_processed/temp -jar $sw/picard.jar SortSam -I $sam_files/$filename -O $analysis_dir/Picard_processed/3.1_sortSam/$samplename.ordered.sam -SORT_ORDER coordinate -TMP_DIR $analysis_dir/Picard_processed/temp -MAX_RECORDS_IN_RAM 5000000 > $analysis_dir/Picard_processed/3.1_sortSam/3.1.picard_sorting.$samplename.txt 2>&1
    
     singularity exec -B /group /group/bioinf/Software/singularity/bacteria.sif java  $memory -Djava.io.tmpdir=$analysis_dir/Picard_processed/temp -jar $sw/picard.jar AddOrReplaceReadGroups -I $analysis_dir/Picard_processed/3.1_sortSam/$samplename.ordered.sam -O $analysis_dir/Picard_processed/3.2_AddOrReplaceReadGroups/$samplename.addReplace.bam -RGLB $filename -RGPL ILLUMINA -RGPU ILLUMINA -RGSM $filename -RGID $samplename -SORT_ORDER coordinate -COMPRESSION_LEVEL 5 -VALIDATION_STRINGENCY SILENT -MAX_RECORDS_IN_RAM 5000000 > $analysis_dir/Picard_processed/3.2_AddOrReplaceReadGroups/3.2.AddOrReplaceReadGroups.$samplename.txt 2>&1
     
     singularity exec -B /group /group/bioinf/Software/singularity/bacteria.sif java  $memory -Djava.io.tmpdir=$analysis_dir/Picard_processed/temp -jar  $sw/picard.jar MarkDuplicates -I $analysis_dir/Picard_processed/3.2_AddOrReplaceReadGroups/$samplename.addReplace.bam -O $analysis_dir/Picard_processed/3.3_markDuplicates/$samplename.marked.bam -CREATE_INDEX true -VALIDATION_STRINGENCY SILENT -ASSUME_SORT_ORDER coordinate -M $analysis_dir/Picard_processed/3.3_markDuplicates/deduppedBamSortedFiles/$samplename.output.metrics -MAX_RECORDS_IN_RAM 5000000 > $analysis_dir/Picard_processed/3.3_markDuplicates/3.3.mark_duplicates.$samplename.txt 2>&1
     
