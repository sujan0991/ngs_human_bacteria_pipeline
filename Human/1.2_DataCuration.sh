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

#*******************************************************************************
#  In this script: 1.2 DATA CURATION
#        - excluding samples (optional)
#        - UMI extraction (optional; usually NOT needed for bacteria)
#        - cutting adapter sequences (optional)
#        - remove low quality reads (optional)
#        - removing poly-sequences (optional)
#*******************************************************************************


# ******************************************************************************
# ******************************************************************************
#                           !USER PROVIDED HEADING!
# ******************************************************************************
# ******************************************************************************

# External shell arguments:
# ******************************************************************************
analysis_dir=$1 # main analysis folder where all performed  steps of the pipeline are saved


folder_nextstep=$2 # folder_nextstep is where the most recent .fastq.gz files will be stored. Includes subfolders "R1" and "R2" if reads are paired-end.
                   # folder_nextstep will be updated as we move from one process to another; if process is skipped folder_nextstep is not updated
                   # we begin where our .fastq.gz files ended up in 1.1_QualityCheckRaw: $analysis_dir/0_fastQ_data
                                           
paired_end=$3 #"Yes" # "Yes" if reads are pair-ended. Otherwise, provide whatever else, e.g. just an empty string "" or "No" or "fhgj"
filename_extention=$4 # extension of the file to work on. Data curation works on .fastq.gz
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
# *******************************************************************************
filename_list_exclude=$analysis_dir/exclude.txt # If you need to exclude certain samples, provide a .txt file with sample names that need to be excluded    
 
softwares=/group/bioinf/Users/Madina/sw/ # Hint: you can add the folder with softwares to the path so that the computer knows where to look for when running a software (e.g. cutadapt, trimmomatic)

UMI_extraction="" # set to "Yes" is UMIs need to be extracted 
UMI_pattern="NNN"
UMI_pattern2="NNN" # if paired_end="Yes"
umi_method=string #alternative: regex

cut_adapters="Yes" # set to "Yes" if sequence(s) were not trimmed yet, otherwise the script will continue without trimming adapters
trimming_software="trimmomatic" # set to "cutadapt" or "trimmomatic". In paired-end mode, Trimmomatic uses the additional information contained in paired reads to better find adapters. 
#adapter_fwd="ACACTCTTTCCCTACACGACGCTCTTCCGATCT" # adapter sequence if cutadapt was chosen as a trimming software
#adapter_rev="GATCGGAAGAGCACACGTCTGAACTCCAGTCAC" # if paired_end="Yes", additionally provide the reverse adapter sequence
#illumina_adapters_file_PE=/group/bioinf/Data/Madina_test/pipeline/adapters/TruSeq2-PE.fa # for trimmomatic a file with adapters needs to be passed. PE
illumina_adapters_file_SE=/group/bioinf/Data/Madina_test/pipeline/adapters/TruSeq2-SE.fa

low_quality="Yes" # set to "Yes" if you need to trim low-quality reads

remove_poly="" # set to "Yes" if needed; otherwise this will be skipped
#poly_A="A{50}"
#poly_G="G{50}"
#poly_T="T{50}"
#poly_C="C{50}"
#wildcard="N{50}"

threads=4 # Cutadapt supports multicore and can processes reads on multiple cores, even while working on a single file. Speeds up reads processing. 
          # Trimmomatic is also a multithreaded tool.
# ******************************************************************************
# ******************************************************************************
# ******************************************************************************

# ******************************************************************************
# Check if the sample is in the filename list to be excluded 
# ******************************************************************************

if [ "$paired_end" = "Yes" ]; then
    if grep -q "$filename\|$filename2" "$filename_list_exclude"; then
        mkdir $analysis_dir/0_excluded_samples/ -p
        mv "$folder_nextstep/R1/$filename" "$analysis_dir/0_excluded_samples/$filename"
        mv "$folder_nextstep/R2/$filename2" "$analysis_dir/0_excluded_samples/$filename2"
        echo "Moved $filename and $filename2 to $analysis_dir/0_excluded_samples/"
    else
        echo "$filename or $filename2 not found in $filename_list_exclude"
    fi
else
    echo "Reads are single-end"
    if grep -q "$filename" "$filename_list_exclude"; then
        mkdir $analysis_dir/0_excluded_samples/ -p
        mv "$folder_nextstep/$filename" "$analysis_dir/0_excluded_samples/$filename"
        echo "Moved $filename to $analysis_dir/0_excluded_samples/"
    else
        echo "$filename not found in $filename_list_exclude"
    fi  
fi

echo "Data curation starts in folder" $folder_nextstep

# ******************************************************************************
# UMI-extraction
# ******************************************************************************

# UMI overview
#       https://umi-tools.readthedocs.io/en/latest/QUICK_START.html
#       UMI - short sequence at the start of a read that uniquely tags each fragment in a sample library. UMI labelling takes place before PCR amplification
#       UMIs need to be removed, but the fragment itself should be kept. Move UMIs to the sample name (extract command from UMI-tools)
#       filter out duplicate reads and PCR errors. so that only the unique reads are kept
#       reduces rates of false positive variant calls  

if [ "$UMI_extraction" = "Yes" ]; then
    echo "UMIs will be extracted."
    mkdir $analysis_dir/1_Preprocessing/1.2_CuratedData/1.2_UMI_extracted/ -p
               
    if [ "$paired_end" = "Yes" ]; then
        echo "Reads are paired-end" 
        mkdir $analysis_dir/1_Preprocessing/1.2_CuratedData/1.2_UMI_extracted/R1 -p
        mkdir $analysis_dir/1_Preprocessing/1.2_CuratedData/1.2_UMI_extracted/R2 -p
        /home/bioinf/maja488d/.local/bin/umi_tools extract \
        --bc-pattern $UMI_pattern \
        --bc-pattern2 $UMI_pattern2 \
        --extract-method $umi_method \
        --read2-in $folder_nextstep/R2/$filename2 \
        --read2-out $analysis_dir/1_Preprocessing/1.2_CuratedData/1.2_UMI_extracted/R2/$filename2 \
        -I $folder_nextstep/R1/$filename  \
        -S $analysis_dir/1_Preprocessing/1.2_CuratedData/1.2_UMI_extracted/R1/$filename > $analysis_dir/1_Preprocessing/1.2_CuratedData/1.2_UMI_extracted/1.2.umi_extract.$samplename.txt 2>&1  
    else
        echo "Reads are single-end"
        /home/bioinf/maja488d/.local/bin/umi_tools extract \
        --bc-pattern $UMI_pattern \
        --extract-method $umi_method \
        -I $folder_nextstep/$filename \
        -S $analysis_dir/1_Preprocessing/1.2_CuratedData/1.2_UMI_extracted/$filename > $analysis_dir/1_Preprocessing/1.2_CuratedData/1.2_UMI_extracted/1.2.umi_extract.$samplename.txt 2>&1     
    fi
    
    folder_nextstep=$analysis_dir/1_Preprocessing/1.2_CuratedData/1.2_UMI_extracted/
    echo "UMIs have been extracted"
    echo "Sample with extracted UMIs is in" $folder_nextstep
else 
    echo "UMI extraction is not needed"
fi

# ******************************************************************************
# Remove adapter sequences if it has not been done 
# ******************************************************************************

# Cutadapt documentation: https://cutadapt.readthedocs.io/en/stable/guide.html?highlight=cores#

# Adapter-search parameters
#       -e: allowed maximum error rate. number of errors/matching length. Default -e 0.1
#       -O: minimum overlap. Reduces the number of falsely trimmed bases. Default -O 3, i.e. the alignment algorithm requires that at least three bases of the adapter are aligned to the read. 
#       --no-indels: allow only for mismathces when searching adapters. If not enabled, error tolerance allows indels and mismatches
#       Adapter-search parameters can be set adapter-specific!
#       see more at https://cutadapt.readthedocs.io/en/stable/guide.html?highlight=cores#

# Speed-up trick. Applies to intermediate short-lived files
#        --compression-level=1
#        if not specified, cutadapt automatically detects that compression is needed (from extension .gz) and the default compression level is 6
#        setting compressio level to 1, will process files faster (no time spent on compressing output). BUT files will take more space memory-wise


if [ "$cut_adapters" = "Yes" ] && [ "$trimming_software" = "cutadapt" ]; then
    echo "Adapters will be removed with Cutadapt"
    mkdir $analysis_dir/1_Preprocessing/1.2_CuratedData/1.2_RemoveAdapters/ -p
    
    if [ "$paired_end" = "Yes" ]; then
        echo "Reads are paired-end"
        mkdir $analysis_dir/1_Preprocessing/1.2_CuratedData/1.2_RemoveAdapters/R1/ -p
        mkdir $analysis_dir/1_Preprocessing/1.2_CuratedData/1.2_RemoveAdapters/R2/ -p
        
        cutadapt \
        --cores=$threads \
        -a $adapter_fwd -A $adapter_rev \
        -e 0.1 \
        -O 3 \
        --discard-trimmed \
        -o $analysis_dir/1_Preprocessing/1.2_CuratedData/1.2_RemoveAdapters/R1/$filename \
        -p $analysis_dir/1_Preprocessing/1.2_CuratedData/1.2_RemoveAdapters/R2/$filename2 \
        $folder_nextstep/R1/$filename \
        $folder_nextstep/R2/$filename2 > $analysis_dir/1_Preprocessing/1.2_CuratedData/1.2_RemoveAdapters/1.2.cutadapt_log.$samplename.txt 2>&1   
               
    else 
        echo "Reads are single-end"
        cutadapt \
        --cores=$threads \
        -a $adapter_fwd \
        -e 0.1 \
        -O 3 \
        --discard-trimmed \
        -o $analysis_dir/1_Preprocessing/1.2_CuratedData/1.2_RemoveAdapters/$filename  \
        $folder_nextstep/$filename > $analysis_dir/1_Preprocessing/1.2_CuratedData/1.2_RemoveAdapters/1.2.cutadapt_log.$samplename.txt 2>&1     
        
    fi
    folder_nextstep=$analysis_dir/1_Preprocessing/1.2_CuratedData/1.2_RemoveAdapters
    echo "Adapter sequences have been removed"
    echo "Samples with removed adapters are in" $folder_nextstep
else 
    echo "Removing adapters with Cutadapt skipped"
fi

# ******************************************************************************
# Cutting adapters with trimmomatic
# ******************************************************************************

#       Trimmomatic can be better at trimming adapters for paired-end reads: palindrome mode
#       Adapter read-through situation: DNA fragment was shorter than the read and the sequences continued reading into the adapter
#       If a read contains the adapter only partially, it is hard to find and trim it. 
#       Trimmomatic in palindrome mode uses a second adapter identification strategy, specifically for adapter read-through,
#       which takes advantage of the added evidence available in paired-end data

#       Options
#           -trimlog: creates a log of all read trimmings wuth detailed info, incl. read name, surviving read length, trimmed amount from start etc
#           we do not save trimlog, because it is too large (several GB for human RNA)
#           instead we save a log of the process summary

if [ "$cut_adapters" = "Yes" ] && [ "$trimming_software" = "trimmomatic" ]; then
    echo "Adapters will be removed with Trimmomatic"
    mkdir $analysis_dir/1_Preprocessing/1.2_CuratedData/1.2_RemoveAdapters/ -p 
        
    if [ "$paired_end" = "Yes" ]; then
        echo "Reads are paired-end" 
          
        mkdir $analysis_dir/1_Preprocessing/1.2_CuratedData/1.2_RemoveAdapters/R1/ -p
        mkdir $analysis_dir/1_Preprocessing/1.2_CuratedData/1.2_RemoveAdapters/R2/ -p  
        mkdir $analysis_dir/1_Preprocessing/1.2_CuratedData/1.2_RemoveAdapters/unpaired -p
        
        java -jar $softwares/trimmomatic-0.39.jar \
        PE \
        -threads $threads \
        $folder_nextstep/R1/$filename $folder_nextstep/R2/$filename2 \
        $analysis_dir/1_Preprocessing/1.2_CuratedData/1.2_RemoveAdapters/R1/$filename \
        $analysis_dir/1_Preprocessing/1.2_CuratedData/1.2_RemoveAdapters/unpaired/$filename \
        $analysis_dir/1_Preprocessing/1.2_CuratedData/1.2_RemoveAdapters/R2/$filename2 \
        $analysis_dir/1_Preprocessing/1.2_CuratedData/1.2_RemoveAdapters/unpaired/$filename2 \
        ILLUMINACLIP:$illumina_adapters_file_PE:2:30:10:3:false > $analysis_dir/1_Preprocessing/1.2_CuratedData/1.2_RemoveAdapters/trimmomatic_adapters_log.$samplename.txt 2>&1  
        
  
   else 
        echo "Reads are single-end"
        java -jar $softwares/trimmomatic-0.39.jar \
        SE \
        -threads $threads \
        $folder_nextstep/$filename \
        $analysis_dir/1_Preprocessing/1.2_CuratedData/1.2_RemoveAdapters/$filename \
        ILLUMINACLIP:$illumina_adapters_file_SE:2:30:10 > $analysis_dir/1_Preprocessing/1.2_CuratedData/1.2_RemoveAdapters/trimmomatic_adapters_log.$samplename.txt 2>&1   
     
    fi
    folder_nextstep=$analysis_dir/1_Preprocessing/1.2_CuratedData/1.2_RemoveAdapters/
    echo "Adapters have been removed with trimmomatic"
    echo "Samples with removed adapters are in" $folder_nextstep
else
    echo "Trimming adapters with Trimmomatic skipped"
fi

# ******************************************************************************
# Removing Low Quality Reads
# ******************************************************************************
#   Adjust these flags to how you want to quality trim:
#       SLIDINGWINDOW: Performs a sliding window trimming approach. It starts scanning at the 5‟ end and clips the read once the average quality within the window falls below a threshold.
#       MAXINFO: An adaptive quality trimmer which balances read length and error rate to maximise the value of each read
#       LEADING: Cut bases off the start of a read, if below a threshold quality
#       TRAILING: Cut bases off the end of a read, if below a threshold quality
#       CROP: Cut the read to a specified length by removing bases from the end
#       HEADCROP: Cut the specified number of bases from the start of the read
#       MINLEN: Drop the read if it is below a specified length
#       AVGQUAL: Drop the read if the average quality is below the specified level
#       TOPHRED33: Convert quality scores to Phred-33
#       TOPHRED64: Convert quality scores to Phred-64 

#       Remark: if you want to trim adapters using trimmomatic  AND quality trim, then both steps can be done at once, i.e. ILLUMINACLIP and Quality related flags
#               The different processing steps occur in the order in which the steps are specified. It is recommended in most cases that adapter clipping (i.e. ILLUMINACLIP flag), if required, is done as early as possible

if [ "$low_quality" = "Yes" ]; then

    echo "Trimming low-quality reads with Trimmomatic"
    mkdir $analysis_dir/1_Preprocessing/1.2_CuratedData/1.2_Quality_trimmed/ -p 
   
    if [ "$paired_end" = "Yes" ]; then
        echo "Reads are paired-end" 
          
        mkdir $analysis_dir/1_Preprocessing/1.2_CuratedData/1.2_Quality_trimmed/R1/ -p
        mkdir $analysis_dir/1_Preprocessing/1.2_CuratedData/1.2_Quality_trimmed/R2/ -p  
        mkdir $analysis_dir/1_Preprocessing/1.2_CuratedData/1.2_Quality_trimmed/unpaired -p
        
        java -jar $softwares/trimmomatic-0.39.jar \
        PE \
        -threads $threads \
        $folder_nextstep/R1/$filename $folder_nextstep/R2/$filename2 \
        $analysis_dir/1_Preprocessing/1.2_CuratedData/1.2_Quality_trimmed/R1/$filename \
        $analysis_dir/1_Preprocessing/1.2_CuratedData/1.2_Quality_trimmed/unpaired/$filename \
        $analysis_dir/1_Preprocessing/1.2_CuratedData/1.2_Quality_trimmed/R2/$filename2 \
        $analysis_dir/1_Preprocessing/1.2_CuratedData/1.2_Quality_trimmed/unpaired/$filename2 \
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36 > $analysis_dir/1_Preprocessing/1.2_CuratedData/1.2_Quality_trimmed/quality_trimming_log.$samplename.txt 2>&1  
        
  
   else 
        echo "Reads are single-end"
        java -jar $softwares/trimmomatic-0.39.jar \
        SE \
        -threads $threads \
        $folder_nextstep/$filename \
        $analysis_dir/1_Preprocessing/1.2_CuratedData/1.2_Quality_trimmed/$filename \
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36 > $analysis_dir/1_Preprocessing/1.2_CuratedData/1.2_Quality_trimmed/quality_trimming_log.$samplename.txt 2>&1   
     
    fi
    folder_nextstep=$analysis_dir/1_Preprocessing/1.2_CuratedData/1.2_Quality_trimmed/
    echo "Low-quality reads have been removed"
    echo "Samples with removed low-quality reads are in" $folder_nextstep
else
    echo "Removing low quality sequences skipped"
fi

# ******************************************************************************
# Remove poly-sequences if needed
# ******************************************************************************
if [ "$remove_poly" = "Yes" ]; then
    echo "Removing poly-sequences with Cutadapt"
    mkdir $analysis_dir/1_Preprocessing/1.2_CuratedData/1.2_RemovePoly/ -p
    
    if [ "$paired_end" = "Yes" ]; then
        echo "Reads are paired-end"
        mkdir $analysis_dir/1_Preprocessing/1.2_CuratedData/1.2_RemovePoly/R1/ -p
        mkdir $analysis_dir/1_Preprocessing/1.2_CuratedData/1.2_RemovePoly/R2/ -p
        
        cutadapt \
        --cores=$threads \
        -b $poly_G -B $poly_G \
        --minimum-length=20 \
        -e 0.1 \
        -O 3 \
        --discard-trimmed \
        -o $analysis_dir/1_Preprocessing/1.2_CuratedData/1.2_RemovePoly/R1/$filename \
        -p $analysis_dir/1_Preprocessing/1.2_CuratedData/1.2_RemovePoly/R2/$filename2 \
        $folder_nextstep/R1/$filename \
        $folder_nextstep/R2/$filename2 > $analysis_dir/1_Preprocessing/1.2_CuratedData/1.2_RemovePoly/1.2.remove_poly.$samplename.txt 2>&1   
               
    else 
        echo "Reads are single-end"
        cutadapt \
        --cores=$threads \
        -b $poly_G \
        --minimum-length=20 \
        -e 0.1 \
        -O 3 \
        --discard-trimmed \
        -o $analysis_dir/1_Preprocessing/1.2_CuratedData/1.2_RemovePoly/$filename  \
        $folder_nextstep/$filename > $analysis_dir/1_Preprocessing/1.2_CuratedData/1.2_RemovePoly/1.2.remove_poly.$samplename.txt 2>&1     
        
    fi
    folder_nextstep=$analysis_dir/1_Preprocessing/1.2_CuratedData/1.2_RemovePoly 
    echo "Poly-sequences have been removed" 
    echo "Samples with removed poly-sequences are in" $folder_nextstep
else
    echo "Removing poly-sequences skipped"
fi

echo "Data curation complete. Files can be found in" $folder_nextstep


















