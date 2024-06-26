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
# In this script: 1.3 QUALITY CONTROL of the CURATED DATA: generate FastQC and MultiQC reports
#         - run fastqc on all curated samples
#         - produce multiqc report for all curated samples
#         - fastqc and multiqc for all reverse and forward reads if data is paired-end
# ******************************************************************************

# ******************************************************************************
#                           !USER PROVIDED HEADING!
# ******************************************************************************

analysis_dir=/group/bioinf/Users/Madina/pipeline_tested/human

threads=1 # For FastQC: -t threads Specifies the number of files which can be processed simultaneously.  Each thread will be allocated 250MB of memory so you shouldn't run more threads than your available memory will cope with, and not more than 6 threads on a 32 bit machine

paired_end="No" # set to "Yes" if you want to additionally compare the quality of forward and reverse reads. Applies to paired-end data

# Input: curated .fastq.gz files. See from step 1.2: "Data curation complete. Files can be found in /group/bioinf/Data/Madina_test/pipeline/v1/1_Preprocessing/1.2_CuratedData/1.2_RemovePoly"
fastQ_files_folder=$analysis_dir/1_Preprocessing/1.2_CuratedData/1.2_Quality_trimmed
fastQ_files_R1=$fastQ_files_folder/R1
fastQ_files_R2=$fastQ_files_folder/R2

# analysis_dir --- full path to the folder where all analysis files will be kept. They will be divided into different subfolders at each each step.
# fastQ_files_folder --- where the raw .fastq data is kept. 
# fastQ_files_R1 --- for paired-end reads, subdirectory for the forward reads. Has path fastQ_files_folder/fastQ_files_R1 
# fastQ_files_R2 --- for paired-end reads, subdirectory for the reverse reads. Has path fastQ_files_folder/fastQ_files_R2
# threads --- number of cores/cpus the process will use. If not provided, the default value is 4.
# ******************************************************************************

# ******************************************************************************
# Curated data 
# ******************************************************************************

# create a vector with all samples
files=$(find "$fastQ_files_folder" -type f -name '*.fastq.gz' -print0 | tr '\0' ' ')


# ******************************************************************************
# Running FastQC
# Installation: https://raw.githubusercontent.com/s-andrews/FastQC/master/INSTALL.txt
# Syntax: https://home.cc.umanitoba.ca/~psgendb/doc/fastqc.help
# Documentation: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
# ******************************************************************************
mkdir $analysis_dir/1_Preprocessing/1.3_QualityCheck_Curated/ -p
mkdir $analysis_dir/1_Preprocessing/1.3_QualityCheck_Curated/1.3_fastqc_reports/ -p 
echo "Started at `date`"
echo "Running FastQC on all samples"

fastqc \
    --outdir $analysis_dir/1_Preprocessing/1.3_QualityCheck_Curated/1.3_fastqc_reports/ \
    --extract \
    -t $threads \
        $files > $analysis_dir/1_Preprocessing/1.3_QualityCheck_Curated/1.3_fastqc_log.txt 2>&1
        
# ******************************************************************************
# Runing multiqc on fastqc reports
# Multiqc aggregates multiple fastqc reports into a single .html report
# Inputs for multiqc are fastqc reports
# https://multiqc.info/docs/
# ******************************************************************************
mkdir $analysis_dir/1_Preprocessing/1.3_QualityCheck_Curated/1.3_multiqc_reports/ -p
# -m fastqc: because we're only reporting on FastQC, we specify this module to speed up multiqc
# --interactive: makes the .html file interactive
# -f: overwrites an existing report
echo "MultiQC: aggregating FastQC reports for all samples"
multiqc -m fastqc  $analysis_dir/1_Preprocessing/1.3_QualityCheck_Curated/1.3_fastqc_reports/ -o $analysis_dir/1_Preprocessing/1.3_QualityCheck_Curated/1.3_multiqc_reports/ --interactive -f > $analysis_dir/1_Preprocessing/1.3_QualityCheck_Curated/1.3_multiqc_log.txt 2>&1
echo "MultiQC report is in" $analysis_dir/1_Preprocessing/1.3_QualityCheck_Curated/1.3_multiqc_reports/ 

#*******************************************************************************
# In case the reads are paired-end, you might want to run quality checks only on forward reads
# or only on reverse reads
#*******************************************************************************
if [ "$paired_end" = "Yes" ]; then
    echo "Running quality control separately on forward and reverse reads"
    
# forward reads
    echo "First checking only the forward reads"
    mkdir $analysis_dir/1_Preprocessing/1.3_QualityCheck_Curated/1.3_FastqcFwdReads/ -p   
    files=$(find "$fastQ_files_folder/R1" -type f -name '*.fastq.gz' -print0 | tr '\0' ' ')
    
fastqc \
    --outdir $analysis_dir/1_Preprocessing/1.3_QualityCheck_Curated/1.3_FastqcFwdReads/ \
    --extract \
    -t $threads \
        $files > $analysis_dir/1_Preprocessing/1.3_QualityCheck_Curated/1.3_fastqc_fwd_reads_log.txt 2>&1
    mkdir $analysis_dir/1_Preprocessing/1.3_QualityCheck_Curated/1.3_MultiqcFwdReads/ -p

echo "MultiQC: aggregating FastQC reports for the forward reads"
multiqc -m fastqc  $analysis_dir/1_Preprocessing/1.3_QualityCheck_Curated/1.3_FastqcFwdReads/ -o  $analysis_dir/1_Preprocessing/1.3_QualityCheck_Curated/1.3_MultiqcFwdReads/ --interactive -f > $analysis_dir/1_Preprocessing/1.3_QualityCheck_Curated/1.3_multiqc_fwd_log.txt 2>&1
echo "MultiQC report for the forward reads (R1) is in" $analysis_dir/1_Preprocessing/1.3_QualityCheck_Curated/1.3_MultiqcFwdReads/

# reverse reads
    echo "Checking only the reverse reads"
    mkdir $analysis_dir/1_Preprocessing/1.3_QualityCheck_Curated/1.3_FastqcRevReads/ -p
   
    files=$(find "$fastQ_files_folder/R2" -type f -name '*.fastq.gz' -print0 | tr '\0' ' ')

      
fastqc \
    --outdir $analysis_dir/1_Preprocessing/1.3_QualityCheck_Curated/1.3_FastqcRevReads \
    --extract \
    -t $threads \
        $files > $analysis_dir/1_Preprocessing/1.3_QualityCheck_Curated/1.3_fastqc_rev_reads_log.txt 2>&1
    mkdir $analysis_dir/1_Preprocessing/1.3_QualityCheck_Curated/1.3_MultiqcRevReads/ -p

echo "MultiQC: aggregating FastQC reports for the reverse reads"
multiqc -m fastqc  $analysis_dir/1_Preprocessing/1.3_QualityCheck_Curated/1.3_FastqcRevReads/ -o $analysis_dir/1_Preprocessing/1.3_QualityCheck_Curated/1.3_MultiqcRevReads/ --interactive -f > $analysis_dir/1_Preprocessing/1.3_QualityCheck_Curated/1.3_multiqc_rev_log.txt 2>&1
echo "MultiQC report for the reverse reads (R2) is in" $analysis_dir/1_Preprocessing/1.3_QualityCheck_Curated/1.3_MultiqcRevReads/

else
    echo "Additional quality checks are not needed"
fi

echo "Finished at `date`"
echo "__DONE__"

    
    
    
    
    
    
    
    
