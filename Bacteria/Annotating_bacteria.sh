#!/usr/bin/env bash
set -e 

# External arguments
# ******************************************************************************
analysis_dir=$1
# ******************************************************************************

# Internal arguments:
# ******************************************************************************
sw=/group/bioinf/Users/Madina/sw/
memory=-Xmx30g
# ******************************************************************************
# ******************************************************************************

mkdir $analysis_dir/5_VariantCalling_Bacteria/ -p
mkdir $analysis_dir/5_VariantCalling_Bacteria/5.10_Annotation/ -p

# input: filtered vcf file from GATK
# to check if chromosome names are same

# check chromosome name in GATK vcf
# awk '{print $1}' $analysis_dir/ecoli.vcf | tail -n1 
# NC_000913 is the chromosome name in GATK vcf file 

# check chromosome name in snpeff database
# singularity exec -B /group /group/bioinf/Software/singularity/bacteria.sif java -Xmx4g -jar $sw/snpEff.jar -v Escherichia_coli
# 'Chromosome'

# rename chromosome in GATK vcf file and run snpeff
cat $analysis_dir/gatk_filtered.vcf | sed "s/NC_000913/Chromosome/" > $analysis_dir/ecoli_v1.vcf
singularity exec -B /group /group/bioinf/Software/singularity/bacteria.sif java $memory -jar $sw/snpEff.jar Escherichia_coli $analysis_dir/ecoli_v1.vcf > $analysis_dir/5_VariantCalling_Bacteria/5.10_Annotation/annotated.vcf


