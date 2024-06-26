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
# ******************************************************************************

# ******************************************************************************
# In this script: STEP V. 5.4-5.6. Joint henotyping. Final output - vcf file for all samples for all genomic intervals
#                         
# https://gatk.broadinstitute.org/hc/en-us/articles/360035535932-Germline-short-variant-discovery-SNPs-Indels- 
# https://hpc.nih.gov/training/gatk_tutorial/genotype-gvcfs.html
# GATK Version: 4.2.5.0                                                                                                                              
# ******************************************************************************


# External shell arguments:
# ******************************************************************************
analysis_dir=$1 #/group/bioinf/Data/Madina_test/pipeline/v1 main analysis folder where all performed  steps of the pipeline are saved
chrom=$2 # list of chromosomes. chr 1 - chr 22 and chr X
# ******************************************************************************
# Internal arguments:
# *******************************************************************************
reference_genome38=/group/bioinf/Data/Madina_test/pipeline/v1/Human/reference_genome/GATK_bundle/hg38/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta
memory=-Xmx30g 
ploidy=2 
sw=/group/bioinf/Users/Madina/sw/
softwares_GATK=/group/bioinf/Users/Madina/sw/gatk-4.2.5.0/gatk-package-4.2.5.0-local.jar # path to GATK .jar file
# ******************************************************************************
# ******************************************************************************
# ******************************************************************************

echo "Started Analysing at `date`"
# ******************************************************************************
# 5.4. Consolidate GVCFs: GenomicsDBImport
#       --genomicsdb-workspace-path SHOULD POINT TO A NON-EXISTENT DIRECTORY  
#       To speedup, GenomicsDBImport was performed on each chromosome.
#       GVCFs are consolidated into a GenomicsDB datastore in order to improve scalability and speedup the next step: joint genotyping
# ******************************************************************************

find $analysis_dir/5_VariantCalling/5.3_BP_RNA_files/ -type f -name "*.g.vcf" -printf "%f\t%p\n" > $analysis_dir/5_VariantCalling/5.3_BP_RNA_files/gvcf_samples.txt  
mkdir $analysis_dir/5_VariantCalling/5.4_GenomicsDBImport -p
 
java $memory -jar $softwares_GATK GenomicsDBImport \
      --sample-name-map $analysis_dir/5_VariantCalling/5.3_BP_RNA_files/gvcf_samples.txt \
      --genomicsdb-workspace-path $analysis_dir/5_VariantCalling/5.4_GenomicsDBImport/$chrom \
      --reference $reference_genome38 \
      -L $chrom \
      --tmp-dir $TMPDIR > $TMPDIR/GenomicsDBImport.chr.$chrom.txt 2>&1

mv $TMPDIR/GenomicsDBImport.chr.$chrom.txt $analysis_dir/5_VariantCalling/5.4_GenomicsDBImport/$chrom


# ******************************************************************************
# 5.5. Joint genotyping per chromosome
#           GenotypeGVCFs uses the potential variants from the HaplotypeCaller and does the joint genotyping. 
#           It will look at the available information for each site from both variant and non-variant alleles across all samples, 
#           and will produce a VCF file containing only the sites that it found to be variant in at least one sample.
# ******************************************************************************
mkdir $analysis_dir/5_VariantCalling/5.5_Genotyping -p
mkdir $analysis_dir/5_VariantCalling/5.5_Genotyping/$chrom -p

java $memory -jar $softwares_GATK GenotypeGVCFs \
	--reference $reference_genome38 \
	-V gendb:////$analysis_dir/5_VariantCalling/5.4_GenomicsDBImport/$chrom \
	--sample-ploidy $ploidy \
	-O $analysis_dir/5_VariantCalling/5.5_Genotyping/$chrom/$chrom.variants.vcf \
    --tmp-dir $TMPDIR > $TMPDIR/Genotyping.chr.$chrom.txt 2>&1

mv $TMPDIR/Genotyping.chr.$chrom.txt $analysis_dir/5_VariantCalling/5.5_Genotyping/$chrom

# ******************************************************************************
# 5.6 Merging per chromosome vcfs into one file
# ****************************************************************************** 
vcf_chr=$analysis_dir/5_VariantCalling/5.5_Genotyping/
mkdir $analysis_dir/5_VariantCalling/5.6_Merging -p
#find $analysis_dir/5_GATK_gvcf/5.5_Genotyping -type f -name "*.vcf" -printf "%p\n" > $analysis_dir/5_GATK_gvcf/5.5_Genotyping/input_variant_files.list

singularity exec -B /group /group/bioinf/Software/singularity/bacteria.sif java  $memory -Djava.io.tmpdir=$analysis_dir/5_VariantCalling/5.6_Merging -jar $sw/picard.jar GatherVcfs \
    -I $vcf_chr/chr1/chr1.variants.vcf \
    -I $vcf_chr/chr2/chr2.variants.vcf \
    -I $vcf_chr/chr3/chr3.variants.vcf \
    -I $vcf_chr/chr4/chr4.variants.vcf \
    -I $vcf_chr/chr5/chr5.variants.vcf \
    -I $vcf_chr/chr6/chr6.variants.vcf \
    -I $vcf_chr/chr7/chr7.variants.vcf \
    -I $vcf_chr/chr8/chr8.variants.vcf \
    -I $vcf_chr/chr9/chr9.variants.vcf \
    -I $vcf_chr/chr10/chr10.variants.vcf \
    -I $vcf_chr/chr11/chr11.variants.vcf \
    -I $vcf_chr/chr12/chr12.variants.vcf \
    -I $vcf_chr/chr13/chr13.variants.vcf \
    -I $vcf_chr/chr14/chr14.variants.vcf \
    -I $vcf_chr/chr15/chr15.variants.vcf \
    -I $vcf_chr/chr16/chr16.variants.vcf \
    -I $vcf_chr/chr17/chr17.variants.vcf \
    -I $vcf_chr/chr18/chr18.variants.vcf \
    -I $vcf_chr/chr19/chr19.variants.vcf \
    -I $vcf_chr/chr20/chr20.variants.vcf \
    -I $vcf_chr/chr21/chr21.variants.vcf \
    -I $vcf_chr/chr22/chr22.variants.vcf \
    -I $vcf_chr/chrX/chrX.variants.vcf \
    -I $vcf_chr/chrY/chrY.variants.vcf \
    -O $analysis_dir/5_VariantCalling/5.6_Merging/5.6_MergedVariants.vcf \
    > $analysis_dir/5_VariantCalling/5.6_Merging/5.6_MergingVCFs.txt 2>&1 

