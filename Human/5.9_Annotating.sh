#!/usr/bin/env bash
set -e 
analysis_dir=$1
download_data_sources="No" # "Yes" or nothing, i.e. empty string "", or whatever else if you do not need to download data sources
variants_type="germline" # "germline" or "somatic"

# If you need to download data sources 

if [ "$download_data_sources" = "Yes" ]; then
    if [ "$variants_type" = "germline" ]; then
    echo "Downloading data sources for germline variants" 
    java $memory -jar $softwares_GATK FuncotatorDataSourceDownloader \
     --germline \
     --validate-integrity \
     --extract-after-download       
    elif [ "$variants_type" = "somatic" ]; then
    echo "Downloading data sources for somatic variants"
    java $memory -jar $softwares_GATK FuncotatorDataSourceDownloader \
    --somatic \
    --validate-integrity \
    --extract-after-download
    else
    "No need to download data sources"
    fi
else
    echo "No need to download data sources"
fi 


# Internal arguments:
# *******************************************************************************
reference_genome38=/group/bioinf/Data/Madina_test/pipeline/v1/Human/reference_genome/GATK_bundle/hg38/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta
memory=-Xmx30g 
softwares_GATK=/group/bioinf/Users/Madina/sw/gatk-4.2.5.0/gatk-package-4.2.5.0-local.jar # path to GATK .jar file
data_sources=/group/bioinf_biomarkers_rna/rna_us_analysis/funcotator_dataSources.v1.7.20200521g
# ******************************************************************************
# ******************************************************************************

# mkdir $analysis_dir/5_VariantCalling/5.9_VariantAnnotator/ -p

java $memory -jar $softwares_GATK Funcotator \
    -R $reference_genome38 \
    -V $analysis_dir/5_VariantCalling/5.8_filtering/filtered.0.2.vcf \
    --ref-version hg38 \
    --data-sources-path $data_sources \
    --output $analysis_dir/5_VariantCalling/5.9_VariantAnnotator/annotated.final.vcf \
    --output-file-format VCF \
    --tmp-dir $TMPDIR > $TMPDIR/5.9_Annotating.txt 2>&1

mv $TMPDIR/5.9_Annotating.txt $analysis_dir/5_VariantCalling/5.9_VariantAnnotator/


