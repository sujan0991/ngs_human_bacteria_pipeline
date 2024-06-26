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
# In this script: ALIGNMENT
#         - building the reference genome
#         - running aligner
# ******************************************************************************

# ******************************************************************************
#                           !USER DEFINED HEADING!
# ******************************************************************************

# External shell arguments:
# ******************************************************************************
analysis_dir=$1 #/group/bioinf/Data/Madina_test/pipeline/v1
fastQ_files_folder=$2 #$analysis_dir/1_Preprocessing/1.2_CuratedData/1.2_RemovePoly #See end of data curation: e.g., "Data curation complete. Files can be found in $analysis_dir/1_Preprocessing/1.2_CuratedData/1.2_RemovePoly"
paired_end=$3 #"Yes" # "Yes" if reads are pair-ended. Otherwise, provide whatever else, e.g. just an empty string "" or "No" or "fhgj"
filename_extention=$4 # extension of the file to work on. Here, BWA alignment works on .fastq.gz
samplename=$5 # name of the sample, e.g. DNA1
readname=$6 # R
# ******************************************************************************

# Constructing the name of the file
if [ "$paired_end" = "Yes" ]; then
    echo "Reads are paired-end" 
    filename="${samplename}_${readname}1${filename_extention}"
    filename2="${samplename}_${readname}2${filename_extention}"
        
else
    echo "Reads are single-end"
    filename="${samplename}_${readname}1${filename_extention}"
fi

# Internal arguments:
# ******************************************************************************
threads=4
indexing_algorithm=is
softwares=/group/bioinf/Users/Madina/sw/
reference_genome=/group/negin_ecoli/ACRAS/reference_Ecoli_K12_MG1655/reference.fna
# ******************************************************************************
# ******************************************************************************

# ******************************************************************************
# 2.1 Generating the reference genome
#       # -a: indexing algorithm. either bwtsw or is. Default: is. It is faster and simpler, but does not work with genomes larger than 2GB. 
#*******************************************************************************

# A: This creates a collection of files used by BWA to perform the alignment.
bwa index -a $indexing_algorithm $reference_genome

# B: This creates a file called reference.fa.fai, with one record per line for each of the contigs in the FASTA reference file. Each record is composed of the contig name, size, location, basesPerLine and bytesPerLine.
samtools faidx $reference_genome

#*******************************************************************************
# STEP 2.2 Aligning
#       BWA may produce multiple primary alignments for different part of a query sequence
#       -M: mark shorter split hits as secondary ---> compatible with Picard tools, e.g. mark duplicates

#       -k: minimal length of matches to start alignment. Default 20. Smaller values: easier to start alignment, but higher chance or multimapping and false positives
#           higher values: more stringent mapping
#       See https://bio-bwa.sourceforge.net/bwa.shtml
#*******************************************************************************
mkdir $analysis_dir/2._AlignmentBacteria -p

if [ "$paired_end" = "Yes" ]; then
    bwa mem -M -k 20 -t $threads $reference_genome $fastQ_files_folder/R1/$filename $fastQ_files_folder/R2/$filename2 > $analysis_dir/2._AlignmentBacteria/$samplename.aligning_log.txt 2>&1 -o $analysis_dir/2._AlignmentBacteria/$samplename.sam      
else               
    bwa mem -M -k 20 -t $threads $reference_genome $fastQ_files_folder/$filename > $analysis_dir/2._AlignmentBacteria/$samplename.aligning_log.txt 2>&1 -o $analysis_dir/2._AlignmentBacteria/$samplename.sam  
fi

#*******************************************************************************
# STEP 2.3 Check for quality of aligning (optional)
#*******************************************************************************

# check the number of mapped  reads
#$softwares/samtools view -F 4 -c $analysis_dir/$samplename.sam

# check the number of unmapped  reads
#$softwares/samtools view -f 4 -c $analysis_dir/$samplename.sam

# checking mapping quality scores 
#$softwares/samtools view -F 4 -h $analysis_dir/$samplename.sam | grep -v "^@" | cut -f 5 | sort | uniq -c


