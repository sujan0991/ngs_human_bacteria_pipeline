#!/usr/bin/Rscript 
print("Start gene expression analysis")
# ******************************************************************************
# Load packages 
# ******************************************************************************
library(GenomicAlignments)
library(GenomicFeatures)
library(BiocParallel)
library(DESeq2, quietly = TRUE)
library(ape, warn.conflicts = FALSE)
library(PCAtools)
library(biomaRt)
library(scatterplot3d)
library(BatchQC)
library(readxl)
library(preprocessCore)
library(gplots)
library(pheatmap)
# ******************************************************************************
# ******************************************************************************

# duplicated reads and gene expression
# Use UMIs can help, see Klepikova et al., 2017
# Swati Parekh, 2016 "The impact of amplification on differential expression analyses by RNA-seq". Conclusion: removing duplicates can worsen the power and the False Discovery Rate (FDR) for differential gene expression. 
# Some suggest though removal of duplicates (Dozgomorov 2015)
# https://support.bioconductor.org/p/83485/
# https://www.biostars.org/p/55648/
# https://www.biostars.org/p/14283/

# several transcripts - one gene. what to do? 
# https://bioinformatics.stackexchange.com/questions/5281/how-to-deal-with-duplicate-genes-having-different-expression-values
# https://www.biostars.org/p/244850/

# ******************************************************************************
# Overview
# 1. Read quantification -> matrix with raw not normalized reads counts
# 2. Annotating 
# 3. Summing transcripts from same genes
# 4. Quality control: investigating for the presence of batch effect
# 5. Differential gene expression analysis using DESeq2
# 6. DESeq2: visualizing results
# 7. PCAtools: everything Principal Components Analysis
# ******************************************************************************

# ******************************************************************************
# USER-PROVIDED HEADING
# ******************************************************************************
regions <- "exons" # do you want to retrieve transcripts only from exons or from all regions. Provide either "exons" or "all"
gtfFile <- "/group/bioinf/Data/Madina_test/pipeline/v1/Human/reference_genome/older/hg38.gtf"
ncores <- 12 # Number of workers (BiocParallel package). Defaults to all cores available as determined by detectCores
# ******************************************************************************
# ******************************************************************************

# ******************************************************************************
# Input: Bam files
# yieldSize: how many reads should be processed at a time
# ******************************************************************************
bamFileList <- c(
  "/group/bioinf_biomarkers_rna/rna_us_analysis/3_Picard_processed/3.3_markDuplicates/BRPC-17-587.marked.bam",
  "/group/bioinf_biomarkers_rna/rna_us_analysis/3_Picard_processed/3.3_markDuplicates/BRPC-20-387.marked.bam",
  "/group/bioinf_biomarkers_rna/rna_us_analysis/3_Picard_processed/3.3_markDuplicates/BRPC-21-363-1.marked.bam",
  "/group/bioinf_biomarkers_rna/rna_us_analysis/3_Picard_processed/3.3_markDuplicates/BRPC-21-363-2.marked.bam",
  "/group/bioinf_biomarkers_rna/rna_us_analysis/3_Picard_processed/3.3_markDuplicates/BRPC-21-412.marked.bam",
  "/group/bioinf_biomarkers_rna/rna_us_analysis/3_Picard_processed/3.3_markDuplicates/BRPC-22-145.marked.bam",
  "/group/bioinf_biomarkers_rna/rna_us_analysis/3_Picard_processed/3.3_markDuplicates/BRPC-22-287.marked.bam",
 "/group/bioinf_biomarkers_rna/rna_us_analysis/3_Picard_processed/3.3_markDuplicates/BRPC-23-015.marked.bam",
 "/group/bioinf_biomarkers_rna/rna_us_analysis/3_Picard_processed/3.3_markDuplicates/BRPC-23-067.marked.bam", 
 "/group/bioinf_biomarkers_rna/rna_us_analysis/3_Picard_processed/3.3_markDuplicates/BRPC-23-096.marked.bam", 
 "/group/bioinf_biomarkers_rna/rna_us_analysis/3_Picard_processed/3.3_markDuplicates/BRPC-23-104-1.marked.bam", 
 "/group/bioinf_biomarkers_rna/rna_us_analysis/3_Picard_processed/3.3_markDuplicates/BRPC-23-104-2.marked.bam", 
 "/group/bioinf_biomarkers_rna/rna_us_analysis/3_Picard_processed/3.3_markDuplicates/BRPC-23-113.marked.bam"
)
bfl <- BamFileList(bamFileList, yieldSize = 50000, index = character())
# ******************************************************************************
# ******************************************************************************

# ******************************************************************************
# Read the gene model from a GTF file
# ******************************************************************************
# Read the transcriptome annotation file
txdb <- makeTxDbFromGFF(gtfFile, format="gtf", organism="Homo sapiens")

# Retrieve transcripts from the transcriptome annotation
if (regions=="exons") {
print("Get the exons grouped by gene:")
eByg <- exonsBy(txdb, by = c("gene"))
} else if (regions=="all") {
print("All genes")
eByg <- genes(txdb)
} else {
print("Please provide trancriptomic regions of interest")
}
# ******************************************************************************
# ******************************************************************************

# ******************************************************************************
# Initialize the BiocParallel backend
# ******************************************************************************
multicoreParam <- MulticoreParam(workers = ncores)
register(multicoreParam)
registered()
# ******************************************************************************
# ******************************************************************************

# ******************************************************************************
# Preparing the count matrix: reads quantification  within exonic  regions
# Raw, not normalized count values. This is important for downstreat DeSeq2 analysis
# ******************************************************************************
# https://bioconductor.org/packages/devel/bioc/vignettes/SummarizedExperiment/inst/doc/SummarizedExperiment.html
# summarizeOverlaps: calculate the overlap between genomic features in aligned bam sample and eByg
# mode:  "Union", "IntersectionStrict", or "IntersectionNotEmpty"
# "Union" : (Default) Reads that overlap any portion of exactly one feature are counted. Reads that overlap multiple features are discarded. This is the most conservative of the 3 modes. 
# inter.feature = TRUE: reads mapping to multiple features are dropped (i.e., not counted). 
# produced data class: RangedSummarizedExperiment
# assumption: experiment was not strand-specific
counteByg <- bplapply(
    bfl, 
    function(x) summarizeOverlaps(
    eByg,
    x,
    mode = "Union", 
    ignore.strand = TRUE, 
    inter.feature = TRUE,
    singleEnd = FALSE, 
    BPPARAM = multicoreParam))
    
countDFeByg <- sapply(seq(along = counteByg), function(x) assays(counteByg[[x]])$counts)
rownames(countDFeByg) <- names(rowRanges(counteByg[[1]]))
colnames(countDFeByg) <- names(bfl)
read_counts <- as.data.frame(countDFeByg)
write.csv(read_counts, file="Gene_expression_analysis/countDFeByg.csv")

# ******************************************************************************
# Annotating
# https://bioconductor.org/packages/release/bioc/vignettes/biomaRt/inst/doc/accessing_ensembl.html
# ******************************************************************************
# annotate transcript ids to gene names
mart <- useMart("ensembl","hsapiens_gene_ensembl")
read_counts <- read.csv(file = "Gene_expression_analysis/countDFeByg.csv", header=TRUE, stringsAsFactors=FALSE)

annotations <- getBM(attributes = c('ensembl_transcript_id_version', 
                            'ensembl_gene_id', 
                            'external_transcript_name',
                            'external_gene_name'),
             filters = 'ensembl_transcript_id_version', 
             values = read_counts[,1],
             mart = mart)
write.csv(annotations, file="Gene_expression_analysis/annotations.csv")
names(read_counts)[names(read_counts) == 'X'] <- 'ensembl_transcript_id_version'
# annotated raw reads
annotated_rr <- merge(annotations,read_counts,by="ensembl_transcript_id_version")
annotated_rr <- annotated_rr[, -c(1,3,4)]
length(unique(annotated_rr$ensembl_gene_id)) #69019 gene ids

counts_sum <- aggregate(annotated_rr[,2:14], by=list(annotated_rr$ensembl_gene_id), "sum")
names(counts_sum)[names(counts_sum) == 'Group.1'] <- 'ensembl_gene_id'
rownames(counts_sum) <- counts_sum$ensembl_gene_id
counts_sum$ensembl_gene_id <- NULL
write.csv(counts_sum, file="Gene_expression_analysis/annotated_raw_counts.csv")

# ******************************************************************************
# Check for the presence of batch effect
# ******************************************************************************
# metadata
metadata <- read_excel("Gene_expression_analysis/duke_peterallen_pdac_2023-07-09.xlsx")
metadata$samplenames <- c(
  "BRPC.17.587.marked.bam",
  "BRPC.20.387.marked.bam",
  "BRPC.21.363.1.marked.bam",
  "BRPC.21.363.2.marked.bam",
  "BRPC.21.412.marked.bam",
  "BRPC.22.145.marked.bam",
  "BRPC.22.287.marked.bam",
  "BRPC.23.015.marked.bam",
  "BRPC.23.067.marked.bam",
  "BRPC.23.096.marked.bam",
  "BRPC.23.104.1.marked.bam",
  "BRPC.23.104.2.marked.bam",
  "BRPC.23.113.marked.bam"
)
metadata$stage <-as.character(metadata$stage)
rownames(metadata) <- metadata$samplenames
metadata$samplenames <- NULL

# simulated data
# quantile normalization. input needs to be a numeric matrix
read_counts_norm <- normalize.quantiles(data.matrix(counts_sum))
read_counts_norm <- data.frame(read_counts_norm)
batch<-metadata$Lab
condition <- metadata$stage
batchQC(read_counts_norm, batch=batch, condition=condition,
        report_file="batch_report.html", report_dir=".",
        report_option_binary="111111111",
        view_report=FALSE, interactive=TRUE)

# ******************************************************************************
# Differential expression analysis with DeSeq2
# ******************************************************************************
colData <- data.frame(
  sampleName = rownames(metadata),
  condition_sample = factor(metadata$sample),
  condition_stage = factor(metadata$stage)
)

# pass as input count matrix per gene level
dds <- DESeqDataSetFromMatrix(countData = counts_sum,
                              colData = colData,
                              design = ~condition_stage + condition_sample)
 
ddsHTSeq <- DESeq(dds)
resHTSeq <- results(ddsHTSeq)
# Sort this table by p-value (smaller p-values on top), and save it to a file so that we can later import it into Excel.
orderedRes <- resHTSeq[ order(resHTSeq$padj), ]
write.csv(as.data.frame(orderedRes), file="Gene_expression_analysis/DESeq2.csv")

# Normalization
# There are several ways to normalize data
# 1
normCounts <- counts(ddsHTSeq, normalized = TRUE)
write.csv(as.data.frame(normCounts), file="Gene_expression_analysis/Normalized_counts.csv")

# 2 varianceStabilizingTransformation
#blind=TRUE should be used for comparing samples in an manner unbiased by prior information on samples, for example to perform sample QA (quality assurance). blind=FALSE should be used for transforming data for downstream analysis, where the full #use of the design information should be made. blind=FALSE will skip re-estimation of the dispersion trend, if this has already been calculated. If many of genes have large differences in counts due to the experimental design, it is important to set #blind=FALSE for downstream analysis.
vst <- assay(vst(ddsHTSeq, blind=FALSE))
# 3rlogTransformation
#This function transforms the count data to the log2 scale in a way which minimizes differences between samples for rows with small counts, and which normalizes with respect to library size. The rlog transformation produces a similar variance #stabilizing effect as varianceStabilizingTransformation, though rlog is more robust in the case when the size factors vary widely. The transformation is useful when checking for outliers or as input for machine learning techniques such as #clustering or linear discriminant analysis. rlog takes as input a DESeqDataSet and returns a RangedSummarizedExperiment object. 
rld <- rlogTransformation(ddsHTSeq, blind=FALSE)
# DeSeq2 also takes as input a SummarizedExperiment object 
# Create a SummarizedExperiment object: count matrix + experiment metadata
#countDF_se <- SummarizedExperiment::SummarizedExperiment(assays = countDFeByg,colData = colData)
#dds <- DESeqDataSet(countDF_se, design = ~condition_stage + condition_sample)

# ******************************************************************************
# Visualizing DeSeq2 results: https://gtpb.github.io/ADER18F/pages/tutorial1.html
# ******************************************************************************

# Hierarchical clustering of samples based on their gene expression data
d <- cor(assay(rlog(dds)), method = "spearman") # pearson: linear; larger datasets
hc <- hclust(dist(1 - d))
png("Gene_expression_analysis/sample_tree.png")
plot.phylo(as.phylo(hc), type = "p", edge.col = "blue", edge.width = 2, show.node.label = TRUE, no.margin = TRUE)
dev.off()

png("Gene_expression_analysis/dispersion.png")
plotDispEsts(ddsHTSeq)
dev.off()  

# inspect the distribution of p-values
png("Gene_expression_analysis/p_values_distribution.png")
hist(resHTSeq$pvalue, breaks=0:50/50, xlab="p value", main="Histogram of nominal p values")
dev.off()

# MA plot displays the relationship between a genes’ mean expression and its fold-change between experimental conditions,
png("Gene_expression_analysis/MA.png")
plotMA(resHTSeq)
dev.off()  

# Volcano plot displays the relationship between fold-change and evidence of differential expression (represented as -log p-value)
png("Gene_expression_analysis/volcano_plot.png")
plot(resHTSeq$log2FoldChange, -log10(resHTSeq$pvalue), xlab="log2 Fold-change", ylab="-log P-value", pch=20, cex=0.5)
points(resHTSeq$log2FoldChange[ resHTSeq$padj<0.05 ], -log10(resHTSeq$pvalue[ resHTSeq$padj<0.05 ]), col="red", pch=20, cex=0.5)
abline(v=0, h=-log10(0.05), lty="dashed", col="grey")
dev.off()

dists <- dist(t(as.data.frame(normCounts)))
# headmap of distances
png("Gene_expression_analysis/heatmap.png")
heatmap(as.matrix(dists), main="Clustering of euclidean distances", scale="none") 
dev.off()

# there are no significant genes, nothing to plot
diffgenes <- rownames(resHTSeq)[ which(resHTSeq$padj < 0.05) ]
diffcounts <- normCounts[ diffgenes, ]
png("Gene_expression_analysis/heatmap.2.png")
heatmap.2(diffcounts, 
          labRow = "", 
          trace = "none", density.info = "none",
          scale = "row",
          distfun = function(x) as.dist(1 - cor(t(x))))
dev.off()

# select the 20 most differentially expressed genes
select <- row.names(orderedRes[1:20, ])
# transform the counts to log10
log10_normCounts <- log10(normCounts + 1)
# get the values for the selected genes
values <- log10_normCounts[ select, ]
png("Gene_expression_analysis/pheatmap.png")
pheatmap(values,
         scale = "none", 
         cluster_rows = FALSE, 
         cluster_cols = FALSE,
         fontsize_row = 8,
         annotation_names_col = FALSE,
         gaps_col = c(3,6),
         display_numbers = TRUE,
         number_format = "%.2f",         
         height=12,
         width=6)
dev.off()
# ******************************************************************************

# ******************************************************************************
# PCAtools: everything Principal Components Analysis
# ******************************************************************************

# rownames of metadata should coincide with colnames of count matrix
normalized_counts <- read.csv(file = "Gene_expression_analysis/Normalized_counts.csv", header = TRUE)
rownames(normalized_counts) <- normalized_counts$X
normalized_counts$X <- NULL
metadata$samplenames <- c(
   "BRPC.17.587.marked.bam",
   "BRPC.20.387.marked.bam",
   "BRPC.21.363.1.marked.bam",
   "BRPC.21.363.2.marked.bam",
   "BRPC.21.412.marked.bam",
   "BRPC.22.145.marked.bam",
   "BRPC.22.287.marked.bam",
   "BRPC.23.015.marked.bam",
   "BRPC.23.067.marked.bam",
   "BRPC.23.096.marked.bam",
   "BRPC.23.104.1.marked.bam",
   "BRPC.23.104.2.marked.bam",
   "BRPC.23.113.marked.bam"
 )
rownames(metadata) <- metadata$samplenames
metadata$samplenames <- NULL
metadata$Lab <- c(rep("A",6), rep("B",7))
p <- pca(normalized_counts, metadata=metadata, removeVar=0.1) 
screeplot <- screeplot(p, axisLabSize=8, titleLabSize = 12)
ggsave("screeplot_pca.png", screeplot, dpi = 700)
biplot <- biplot(p, showLoadings=TRUE, labSize=3, pointSize=5, sizeLoadingsNames=1, colby = 'stage') # Logical,
ggsave("Gene_expression_analysis/biplot_pca.png", biplot, dpi = 1200)



project.pca <- prcomp(t(normalized_counts))
summary(project.pca)
#Determine the proportion of variance of each component
#Proportion of variance equals (PC stdev^2) / (sum all PCs stdev^2)
project.pca.proportionvariances <- ((project.pca$sdev^2) / (sum(project.pca$sdev^2)))*100
# scree plot
png("Gene_expression_analysis/screeplot2.png")
barplot(project.pca.proportionvariances, cex.names=1, xlab=paste("Principal component (PC), 1-", length(project.pca$sdev)), ylab="Proportion of variation (%)", main="Scree plot", ylim=c(0,100))
dev.off()

# pairs plot
png("Gene_expression_analysis/pairs_plot_PC1_PC6.png")
par(cex=1.0, cex.axis=0.8, cex.main=0.8)
pairs(project.pca$x[,1:6], col="black", main="Principal components analysis bi-plot\nPCs 1-6", pch=16)
dev.off()
png("Gene_expression_analysis/pairs_plot_PC7_PC13.png")
par(cex=1.0, cex.axis=0.8, cex.main=0.8)
pairs(project.pca$x[,7:13], col="black", main="Principal components analysis bi-plot\nPCs 7-13", pch=16)
dev.off()

#biplots
png("Gene_expression_analysis/pairwise_pca.png")
par(mar=c(4,4,4,4), mfrow=c(1,3), cex=1.0, cex.main=0.8, cex.axis=0.8)

#Plots scatter plot for PC 1 and 2
plot(project.pca$x, type="n", main="Principal components analysis bi-plot", xlab=paste("PC1, ", round(project.pca.proportionvariances[1], 2), "%"), ylab=paste("PC2, ", round(project.pca.proportionvariances[2], 2), "%"))
points(project.pca$x, col="black", pch=16, cex=1)

#Plots scatter plot for PC 1 and 3
plot(project.pca$x[,1], project.pca$x[,3], type="n", main="Principal components analysis bi-plot", xlab=paste("PC1, ", round(project.pca.proportionvariances[1], 2), "%"), ylab=paste("PC3, ", round(project.pca.proportionvariances[3], 2), "%"))
points(project.pca$x[,1], project.pca$x[,3], col="black", pch=16, cex=1)

#Plots scatter plot for PC 2 and 3
plot(project.pca$x[,2], project.pca$x[,3], type="n", main="Principal components analysis bi-plot", xlab=paste("PC2, ", round(project.pca.proportionvariances[2], 2), "%"), ylab=paste("PC3, ", round(project.pca.proportionvariances[3], 2), "%"))
points(project.pca$x[,2], project.pca$x[,3], col="black", pch=16, cex=1)
dev.off()

cat("End gene expresion analysis")








