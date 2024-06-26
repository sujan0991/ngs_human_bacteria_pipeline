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
# In this script: 
#                2.2. First-pass STAR alignment
#                2.3. Splice junctions database generation
#                2.4. Second-pass STAR alignment
# ******************************************************************************

# ******************************************************************************
#                           !USER DEFINED HEADING!
# ******************************************************************************
# External shell arguments:
# ******************************************************************************
analysis_dir=$1 #/group/bioinf/Data/Madina_test/pipeline/v1/Human
processed_files=$2 #/group/bioinf/Data/Madina_test/pipeline/v1/Human/0_fastQ_data_links                                        
paired_end=$3 #"Yes" if reads are pair-ended
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
genomeDir=/group/bioinf/Data/Madina_test/pipeline/human/ega/2_mapping_STAR/2.1_STAR_genomeDir
reference_genome38=/group/bioinf/Data/Madina_test/pipeline/v1/Human/reference_genome/GATK_bundle/hg38/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta #FASTA file with reference genome (mandatory)
reference_gtf=/group/bioinf/Data/Madina_test/pipeline/v1/Human/reference_genome/older/hg38.gtf #GTF file with annptated transcripts (highly recommended)
threads=4
overhang=100
# ******************************************************************************
# STAR options:
#       --outFileNamePrefix: if not specified, all aligned .sam files will be named "Aligned.out.sam", which is not optimal. 
#                            therefore, we assign sample specific names. If this option is used, do not forget to adjust the --sjdbFileChrStartEnd option to a sample specific SJ.out.tab file!
#       --genomeChrBinNbits: if working with a large (>5000 references), this helps to reduce RAM consumption
#                            recommended min(18,log2[max(GenomeLength/NumberOfReferences,ReadLength)])
#       --outSAMtype: technically we can do first step of Picard processing - sorting by coordinate - in the mapper, by setting this option to BAM SortedByCoordinate

# ******************************************************************************
# Paired-end case
# ******************************************************************************
if [ "$paired_end" = "Yes" ]; then 

    echo "Reads are paired-end"
          
        temp_folder=/home/bioinf/maja488d/${samplename}.temp
        if [ -d "$temp_folder" ]; then
            rm -r $temp_folder
        fi    
        
        mkdir $analysis_dir/2_mapping_STAR/ -p
        mkdir $analysis_dir/2_mapping_STAR/2.2_Mapping/star_1pass/ -p
        mkdir $analysis_dir/2_mapping_STAR/2.2_Mapping/star_1pass/${samplename}/ -p
        cd $analysis_dir/2_mapping_STAR/2.2_Mapping/star_1pass/${samplename}/
        
# ******************************************************************************
# 2.2 Alignment 
# ******************************************************************************

        echo "Starting step 2.2: first-pass STAR alignment for sample" $sample

        STAR --runMode alignReads \
             --genomeDir  $genomeDir \
             --readFilesIn $processed_files/R1/$filename $processed_files/R2/$filename2 \
             --readFilesCommand zcat \
             --runThreadN $threads \
             --outFileNamePrefix $samplename \
             --outTmpDir $temp_folder

        echo "Finished step 2.2 for sample" ${samplename}
               
        if [ -d "$temp_folder" ]; then
            rm -r $temp_folder
        fi
        temp_folder=/home/bioinf/maja488d/${samplename}.temp

# ******************************************************************************
# 2.3 second-pass STAR 
# For the 2-pass STAR, a new index is then created using splice junction information contained in the file SJ.out.tab from the first pass
# *******************************************************************************

        mkdir $analysis_dir/2_mapping_STAR/2.1_STAR_genomeDir_pass/ -p
        mkdir $analysis_dir/2_mapping_STAR/2.1_STAR_genomeDir_pass/${samplename}/ -p

        echo "Starting step 2.3: Generating splice junctions database for sample" ${samplename}

       STAR  --runMode genomeGenerate \
             --genomeDir  $analysis_dir/2_mapping_STAR/2.1_STAR_genomeDir_pass/${samplename}/ \
             --genomeFastaFiles  $reference_genome38 \
             --sjdbFileChrStartEnd  $analysis_dir/2_mapping_STAR/2.2_Mapping/star_1pass/${samplename}/"${samplename}SJ.out.tab" \
             --sjdbOverhang 75 \
             --limitSjdbInsertNsj 2000000 \
             --runThreadN $threads \
             --outFileNamePrefix $samplename \
             --outTmpDir $temp_folder

        echo "Finished step 2.3 for sample" $samplename 
        
        if [ -d "$temp_folder" ]; then
            rm -r $temp_folder
        fi
        temp_folder=/home/bioinf/maja488d/${samplename}.temp

# ******************************************************************************
# 4 Alignment with 2-pass STAR 
# The resulting index is then used to produce the final alignments as follows
# ******************************************************************************
        mkdir  $analysis_dir/2_mapping_STAR/2.2_Mapping/star_2pass/ -p
        mkdir  $analysis_dir/2_mapping_STAR/2.2_Mapping/star_2pass/${samplename}/ -p
        cd  $analysis_dir/2_mapping_STAR/2.2_Mapping/star_2pass/${samplename}/

        echo "Starting step 2.4: second-pass STAR alignment for sample" ${samplename} 

        STAR --runMode alignReads \
             --genomeDir $analysis_dir/2_mapping_STAR/2.1_STAR_genomeDir_pass/${samplename}/ \
             --outSAMstrandField intronMotif \
             --readFilesIn $processed_files/R1/$filename $processed_files/R2/$filename2 \
             --readFilesCommand zcat \
             --runThreadN $threads \
             --outFileNamePrefix $samplename \
             --outTmpDir $temp_folder

        echo "Finished step 2.4 for sample "${samplename}

    rm -rf $analysis_dir/2_mapping_STAR/2.2_Mapping/star_1pass/*
    rm -rf  $analysis_dir/2_mapping_STAR/2.1_STAR_genomeDir_pass/*

# ******************************************************************************
# ******************************************************************************
else 
    echo "Reads are single-end"
    
        temp_folder=/home/bioinf/maja488d/${samplename}.temp
        if [ -d "$temp_folder" ]; then
            rm -r $temp_folder
        fi


        mkdir $analysis_dir/2_mapping_STAR/2.2_Mapping/star_1pass/ -p
        mkdir $analysis_dir/2_mapping_STAR/2.2_Mapping/star_1pass/$samplename/ -p
        cd $analysis_dir/2_mapping_STAR/2.2_Mapping/star_1pass/$samplename/
        
# ******************************************************************************
# 2.2 Alignment 
# ******************************************************************************


        echo "Starting step 2.2: first-pass STAR alignment for sample" $samplename

        STAR --runMode alignReads \
             --genomeDir  $genomeDir \
             --readFilesIn $processed_files/$filename \
             --readFilesCommand zcat \
             --runThreadN $threads \
             --outFileNamePrefix $samplename \
             --outTmpDir $temp_folder

        echo "Finished step 2.2 for sample" $samplename
    
            if [ -d "$temp_folder" ]; then
                rm -r $temp_folder
            fi
        temp_folder=/home/bioinf/maja488d/${samplename}.temp
        
# ******************************************************************************
# 2.3 second-pass STAR 
# For the 2-pass STAR, a new index is then created using splice junction information contained in the file SJ.out.tab from the first pass
# *******************************************************************************

        mkdir $analysis_dir/2_mapping_STAR/2.1_STAR_genomeDir_pass/ -p
        mkdir $analysis_dir/2_mapping_STAR/2.1_STAR_genomeDir_pass/$samplename/ -p

        echo "Starting step 2.3: Generating splice junctions database for sample" $samplename

        STAR --runMode genomeGenerate \
             --genomeDir  $analysis_dir/2_mapping_STAR/2.1_STAR_genomeDir_pass/$samplename \
             --genomeFastaFiles  $reference_genome38 \
             --sjdbFileChrStartEnd  $analysis_dir/2_mapping_STAR/2.2_Mapping/star_1pass/$samplename/"${samplename}SJ.out.tab" \
             --sjdbOverhang $overhang \
             --runThreadN $threads \
             --outFileNamePrefix $samplename \
             --outTmpDir $temp_folder

        echo "Finished step 2.3 for sample"  $samplename
    
        if [ -d "$temp_folder" ]; then
                rm -r $temp_folder
            fi
        temp_folder=/home/bioinf/maja488d/${samplename}.temp

# ******************************************************************************
# 4 Alignment with 2-pass STAR 
# The resulting index is then used to produce the final alignments as follows
# ******************************************************************************

        mkdir  $analysis_dir/2_mapping_STAR/2.2_Mapping/star_2pass/ -p
        mkdir  $analysis_dir/2_mapping_STAR/2.2_Mapping/star_2pass/$samplename/ -p
        cd  $analysis_dir/2_mapping_STAR/2.2_Mapping/star_2pass/$samplename

        echo "Starting step 2.4: second-pass STAR alignment for sample" $samplename   

        STAR \
                  --genomeDir $analysis_dir/2_mapping_STAR/2.1_STAR_genomeDir_pass/$samplename \
                  --readFilesIn $processed_files/$filename \
                  --readFilesCommand zcat \
                  --runThreadN $threads \
                  --outFileNamePrefix $samplename \
                  --outTmpDir $temp_folder

        echo "Finished step 2.4 for sample" $samplename
    
        rm -rf $analysis_dir/2_mapping_STAR/2.2_Mapping/star_1pass/$samplename
        rm -rf  $analysis_dir/2_mapping_STAR/2.1_STAR_genomeDir_pass/$samplename

fi

