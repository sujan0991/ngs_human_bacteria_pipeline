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
# In this script: Quality check of bam files
# ******************************************************************************

# ******************************************************************************
#                           !USER DEFINED HEADING!
# ******************************************************************************
# External shell arguments:
analysis_dir=$1 #/group/bioinf/Data/Madina_test/pipeline/v1
bam_files=$2 #$analysis_dir/Picard_processed/3.3_markDuplicates
filename_extension=$3 #.bam
samplename=$4 
# ******************************************************************************
# Constructing the filename
filename="${samplename}${filename_extension}"
# ******************************************************************************

# ******************************************************************************
# QUALIMAP:
# 1. unset user display before running qualimap: type "unset DISPLAY" in the terminal. To set it back to the local system: DISPLAY=:0.0
# 2. make sure that folder with qualimap is appended to the PATH. Otherwise, qualimap command will not be found
export PATH=/group/bioinf/Users/Madina/sw/qualimap_v2.2.1:$PATH
# ******************************************************************************

# ******************************************************************************
# Qualimap
#       bamqc            Evaluate NGS mapping to a reference genome
#       rnaseq           Evaluate RNA-seq alignment data
#       counts           Counts data analysis (further RNA-seq data evaluation)
#       multi-bamqc      Compare QC reports from multiple NGS mappings
#       clustering       Cluster epigenomic signals
#       comp-counts      Compute feature counts
# Options
#       -outdir directory where qualimap reports will be stored passed to the 
#       If outdir is not provided, it will be created automatically in the same folder where BAM file is located.
#       -p Sequencing library protocol Can be either of the three: strand-specific-forward, strand-specific-reverse or non-strand-specific (default)
#       --run-bamqc Raw BAM files are provided as input. If this option is activated BAM QC process first will be run for each sample, then multi-sample analysis will be performed
#       documentation: http://qualimap.conesalab.org/doc_html/command_line.html#multi-sample-bam-qc   
# For speed
#       --java-mem-size -> The bigger the faster (the garbage collector will work less often).
#       -nw (number of windows) -> The bigger the slower with better resolution on the results and less Java memory needed. The smaller the faster with less resolution on the results and more Java memory needed. Every time the end of a window is reached, we launch a thread that computes the corresponding statistics for the window (please, read this for more info: http://qualimap.bioinfo.cipf.es/doc_html/analysis.html#advanced-parameters). Unfortunately there is no general rule to set an optimal value since it will depend on the data being analysed (peaks of coverage imply more memory and time needed for one window).
#       -nr (number of reads in the chunk) -> It controls how many reads are stored in RAM before they are computed. We do this to be able to keep the threads busy as much time as possible since the I/O from the disk is very time consuming. This should therefore be optimized with respect to the number of threads the Java memory size and the number of windows in order to minimize the time where threads are idle.
# ******************************************************************************
mkdir $analysis_dir/4_MappingQualityCheck/ -p
mkdir $analysis_dir/4_MappingQualityCheck/4_Qualimap_per_sample/$samplename -p 
mkdir $analysis_dir/4_MappingQualityCheck/4_Qualimap_aggregate/ -p

qualimap bamqc \
-bam $bam_files/$filename \
-outdir $analysis_dir/4_MappingQualityCheck/4_Qualimap_per_sample/$samplename \
--java-mem-size=16G \ > $analysis_dir/4_MappingQualityCheck/4_Qualimap_per_sample/$samplename/4_Qualimap_log.txt 2>&1

# aggregate individual qualimap reports into one
multiqc -m qualimap $analysis_dir/4_MappingQualityCheck/4_Qualimap_per_sample/ -o $analysis_dir/4_MappingQualityCheck/4_Qualimap_aggregate/ --interactive -f > $analysis_dir/4_MappingQualityCheck/4_Qualimap_aggregate/4_multiqc_Qualimap_log.txt 2>&1

