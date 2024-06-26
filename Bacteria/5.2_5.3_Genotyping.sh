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
# 4.1. Log.final.out to look at sumary stats of mapping 
#       mapping rate. >=75% is a good quality sample
#       multimapping. If multimapping is high, there are usually two issues: 1) sample was contaminated (-> check .fastq data quality) 2) alignment did not work (-> check mapping options)

# 4.2. Qualimap
# https://hbctraining.github.io/Intro-to-rnaseq-hpc-salmon/lessons/03_QC_STAR_and_Qualimap_run.html#qualimap

# V GATK Variant Calling
# 5.1. Call variants per sample - HaplotypeCaler
# 5.2. Consolidate gvcf-s into a single .vcf file
# 5.3. Genotyping 
# 5.4. Split SNPs and Indels
# ******************************************************************************

# ******************************************************************************
# In this script: STEP V. GATK Variant Calling: Germline short variant discovery (SNPs + Indels)
#                         
# https://gatk.broadinstitute.org/hc/en-us/articles/360035535932-Germline-short-variant-discovery-SNPs-Indels- 
# GATK Version: 4.2.5.0                                                                                                                              
# ******************************************************************************


# External shell arguments:
# ******************************************************************************
analysis_dir=$1 # main analysis folder where all performed  steps of the pipeline are saved
# ******************************************************************************
# Internal arguments:
# *******************************************************************************
reference_files=/group/bioinf/Data/Madina_test/bacteria/reference_genomes/
reference_genome=$reference_files/reference.fna
memory=-Xmx40g
ploidy=1
sw=/group/bioinf/Users/Madina/sw/
softwares_GATK=/group/bioinf/Users/Madina/sw/gatk-4.2.5.0/gatk-package-4.2.5.0-local.jar # path to GATK .jar file
# ******************************************************************************
# ******************************************************************************
# ******************************************************************************

echo "Started Analysing at `date`"
# ******************************************************************************
# 5.2. Consolidate GVCFs: GenomicsDBImport
# --genomicsdb-workspace-path SHOULD POINT TO A NON-EXISTENT DIRECTORY  
# ******************************************************************************

# provide sample GVCFs in a map file 
# The sample map is a tab-delimited text file with sample_name--tab--path_to_sample_vcf per line.
find $analysis_dir/5_VariantCalling/5.1_HaplotypeCaller -type f -name "*.g.vcf" -printf "%f\t%p\n" > $analysis_dir/5_VariantCalling/5.1_HaplotypeCaller/gvcf_samples.txt    

# --genomicsdb-workspace-path SHOULD POINT TO A NON-EXISTENT DIRECTORY 
# -L One or more genomic intervals over which to operate. In the case of e.coli it is easy, because the organism has only once chromosome. So we work on that entire chromosome "NC_000913"
java $memory -jar $softwares_GATK GenomicsDBImport \
      --sample-name-map $analysis_dir/5_VariantCalling/5.1_HaplotypeCaller/gvcf_samples.txt \
      --genomicsdb-workspace-path $analysis_dir/5_VariantCalling/5.2_ConsolidateGVCFs \
      --reference $reference_genome \
      -L NC_000913:1-4641652 \
      --tmp-dir $TMPDIR > $TMPDIR/5.2_ConsolidateGVCFs.txt 2>&1

mv $TMPDIR/5.2_ConsolidateGVCFs.txt $analysis_dir/5_VariantCalling/5.2_ConsolidateGVCFs

# ******************************************************************************
# 5.3. Genotyping 
# ******************************************************************************
mkdir $analysis_dir/5_VariantCalling/5.3_Genotyping -p
java $memory -jar $softwares_GATK GenotypeGVCFs \
	--reference $reference_genome \
	-V gendb:////$analysis_dir/5_VariantCalling/5.2_ConsolidateGVCFs \
	--sample-ploidy $ploidy \
	-O $analysis_dir/5_VariantCalling/5.3_Genotyping/variants.vcf \
    --tmp-dir $TMPDIR > $TMPDIR/5.3_Genotyping.txt 2>&1

mv $TMPDIR/5.3_Genotyping.txt $analysis_dir/5_VariantCalling/5.3_Genotyping

