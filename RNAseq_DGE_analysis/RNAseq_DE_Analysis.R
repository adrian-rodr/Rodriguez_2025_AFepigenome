


####################################################
####################################################
############ Re-analysis of GE data ################
############ (Alison et al., 2019) #################
####################################################
####################################################


library(dplyr)
library(ggplot2)
library(DESeq2)
library(Rsubread)
library(biomaRt)
library(ggvenn)
library(ggsignif)



STAR_BAM_files <- c("GSM3666571_SR_LA_Homo_sapiens.GRCh38_Male_SRR8714853-72_7_RNA-Seq.bam",
                    "GSM3666572_SR_RA_Homo_sapiens.GRCh38_Male_SRR8714873-92_7_RNA-Seq.bam",
                    "GSM3666573_SR_LA_Homo_sapiens.GRCh38_Male_SRR8714893-12_12b_RNA-Seq.bam",
                    "GSM3666574_SR_RA_Homo_sapiens.GRCh38_Male_SRR8714913-32_12b_RNA-Seq.bam",
                    "GSM3666575_SR_LA_Homo_sapiens.GRCh38_Male_SRR8714933-52_20_RNA-Seq.bam",
                    "GSM3666576_SR_RA_Homo_sapiens.GRCh38_Male_SRR8714953-72_20_RNA-Seq.bam",
                    "GSM3666577_SR_LA_Homo_sapiens.GRCh38_Male_SRR8714973-92_27_RNA-Seq.bam",
                    "GSM3666578_SR_RA_Homo_sapiens.GRCh38_Male_SRR8714993-5012_27_RNA-Seq.bam",
                    "GSM3666579_SR_LA_Homo_sapiens.GRCh38_Male_SRR8715013-32_30_RNA-Seq.bam",
                    "GSM3666580_SR_RA_Homo_sapiens.GRCh38_Male_SRR8715033-52_30_RNA-Seq.bam",
                    "GSM3666581_AF_LA_Homo_sapiens.GRCh38_Male_SRR8715053-72_5_RNA-Seq.bam",
                    "GSM3666582_AF_RA_Homo_sapiens.GRCh38_Male_SRR8715073-82_5_RNA-Seq.bam",
                    "GSM3666583_AF_LA_Homo_sapiens.GRCh38_Male_SRR8715083-112_14_RNA-Seq.bam",
                    "GSM3666584_AF_RA_Homo_sapiens.GRCh38_Male_SRR8715113-32_14_RNA-Seq.bam",
                    "GSM3666585_AF_LA_Homo_sapiens.GRCh38_Male_SRR8715133-52_25_RNA-Seq.bam",
                    "GSM3666586_AF_RA_Homo_sapiens.GRCh38_Male_SRR8715153-72_25_RNA-Seq.bam",
                    "GSM3666587_AF_LA_Homo_sapiens.GRCh38_Male_SRR8715173-92_28_RNA-Seq.bam",
                    "GSM3666588_AF_RA_Homo_sapiens.GRCh38_Male_SRR8715193-212_28_RNA-Seq.bam",
                    "GSM3666589_AF_LA_Homo_sapiens.GRCh38_Male_SRR8715213-32_31_RNA-Seq.bam",
                    "GSM3666590_AF_RA_Homo_sapiens.GRCh38_Male_SRR8715233-52_31_RNA-Seq.bam"
)


count_matrix <- featureCounts(files=STAR_BAM_files, annot.ext ="Homo_sapiens.GRCh38.108.gtf", isGTFAnnotationFile = TRUE, isPairedEnd = TRUE, requireBothEndsMapped = TRUE, countMultiMappingReads = FALSE, allowMultiOverlap = FALSE)

save(count_matrix, file="count_matrix.RData")


setwd("/Rodriguez_2025_AFepigenome/RNAseq_DGE_analysis/")

load("count_matrix.RData")
class(count_matrix)

# For DESeq2, the matrix of read counts can be directly provided from the "counts" element in the list output.




#########################################################################
#########################################################################
#########################              ##################################
#########################   AF vs SR   ##################################
#########################              ##################################
#########################################################################
#########################################################################


# Construct DESeqDataSet object 

coldata <- read.delim("coldata.txt")

dds_STAR_featureCounts <- DESeqDataSetFromMatrix(countData = count_matrix$counts,
                                                 colData = coldata,
                                                 design = ~ Condition + Tissue)



#########################

# Use bioMart to retrieve human genome genes metadata (Ensembl ID, gene symbol and biotype) to incorporate info to rowData of se object

library(biomaRt)
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl") #This will connect us to the most recent version of the Human Genes BioMart
# Create vactor with all possible values that can be given to 'biotype' filter
all_biotypes_vector <- listFilterOptions(ensembl, filter='biotype')
gene_annotation <- getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol','gene_biotype'),
                         filters = 'biotype', #We're filtering for all biotypes of genes (protein coding, pseudogenes, etc.)
                         values = all_biotypes_vector, 
                         mart = ensembl)
# Once we have a data frame with all the metadata for genes, we can incorporate this columns in the rowData data frame of our DESeqDataSet object
# SOURCE: https://support.bioconductor.org/p/62374/


#To do this, you first have to make sure that these new columns line up row for row with the SummarizedExperiment object
all(rownames(dds_STAR_featureCounts) == gene_annotation$ensembl_gene_id)

#If this is not the case, you need to reorder the gene_annotation data.frame using the match() function, so that the Ensembl IDs match up, then check again that everything is lined up
gene_annotation <- gene_annotation[match(rownames(dds_STAR_featureCounts), gene_annotation$ensembl_gene_id),]
dim(gene_annotation) #62,703
dim(dds_STAR_featureCounts) #62,703
gene_annotation
all(rownames(dds_STAR_featureCounts) == gene_annotation$ensembl_gene_id)  #This is now TRUE

sum(is.na(gene_annotation$ensembl_gene_id))
sum(duplicated(gene_annotation$ensembl_gene_id))
# If the sums above are both 0, then you proceed to add to the metadata columns of the rowData with the following command:
mcols(dds_STAR_featureCounts) <- cbind(mcols(dds_STAR_featureCounts), gene_annotation)
dim(dds_STAR_featureCounts)
rowData(dds_STAR_featureCounts)


#### We're going to run the analyses filtering for protein coding genes first only:
unique(rowData(dds_STAR_featureCounts)$gene_biotype) #remove duplicate elements/rows
dds_STAR_featureCounts_protein_coding <- subset(dds_STAR_featureCounts,rowData(dds_STAR_featureCounts)$gene_biotype=="protein_coding")
dim(dds_STAR_featureCounts_protein_coding) #20,036 genes

#Pre-filtering rows with very small counts
# There are many rows (genes) with zero or very small counts --> Remove those which have a total count of less than 5
dds_STAR_featureCounts_protein_coding_filtered <- dds_STAR_featureCounts_protein_coding[ rowSums(assay(dds_STAR_featureCounts_protein_coding)) >= 5, ]
dim(dds_STAR_featureCounts_protein_coding_filtered) #18,602 genes

#Specify the reference level for each variable
#Set SR as the reference level
colData(dds_STAR_featureCounts_protein_coding_filtered)$Condition <- factor(colData(dds_STAR_featureCounts_protein_coding_filtered)$Condition)
colData(dds_STAR_featureCounts_protein_coding_filtered)$Condition <- relevel(colData(dds_STAR_featureCounts_protein_coding_filtered)$Condition, ref="SR")
colData(dds_STAR_featureCounts_protein_coding_filtered)$Condition

#Set LA as the reference level
colData(dds_STAR_featureCounts_protein_coding_filtered)$Tissue <- factor(colData(dds_STAR_featureCounts_protein_coding_filtered)$Tissue)
colData(dds_STAR_featureCounts_protein_coding_filtered)$Tissue <- relevel(colData(dds_STAR_featureCounts_protein_coding_filtered)$Tissue, ref="LA")
colData(dds_STAR_featureCounts_protein_coding_filtered)$Tissue
colData(dds_STAR_featureCounts_protein_coding_filtered)$Tissue


###########################################################################################################


#### DESeq2

colnames(dds_STAR_featureCounts_protein_coding_filtered) <- colData(dds_STAR_featureCounts_protein_coding_filtered)$names


### EXPLORATORY ANALYSIS AND VISUALIZATION ###

# Visual exploration of sample relationships --> Transformation of the counts for computing distances or making plots

vsd_STAR_featureCounts_protein_coding_filtered <- DESeq2::vst(dds_STAR_featureCounts_protein_coding_filtered, blind=FALSE)
class(vsd_STAR_featureCounts_protein_coding_filtered)

# PCA

DESeq2::plotPCA(vsd_STAR_featureCounts_protein_coding_filtered, intgroup = "Tissue")
DESeq2::plotPCA(vsd_STAR_featureCounts_protein_coding_filtered, intgroup = "Condition")


### DIFFERENTIAL EXPRESSION ANALYSIS ###

dds_STAR_featureCounts_protein_coding_filtered <- DESeq2::DESeq(dds_STAR_featureCounts_protein_coding_filtered, minReplicatesForReplace=Inf)


### RA vs LA
res_RAvsLA_STAR_featureCounts <- DESeq2::results(dds_STAR_featureCounts_protein_coding_filtered, contrast=c("Tissue","RA","LA"), cooksCutoff=FALSE, independentFiltering=FALSE)

head(res_RAvsLA_STAR_featureCounts)
## The following code outputs the set of RA vs LA DEGs (as shown in 📎 Supplementary Table 4).
write.csv(res_RAvsLA_STAR_featureCounts,"res_RAvsLA_STAR_featureCounts2.csv") 
summary(res_RAvsLA_STAR_featureCounts)

# 📊 Figure S4b
res_RAvsLA_STAR_featureCounts_csv <- read.csv("res_RAvsLA_STAR_featureCounts2.csv")
res_RAvsLA_STAR_featureCounts_csv <- res_RAvsLA_STAR_featureCounts_csv[res_RAvsLA_STAR_featureCounts_csv$X != 'ENSG00000277836', ] # remove outlier for visualization improvement
res_RAvsLA_STAR_featureCounts_csv_mut = mutate(res_RAvsLA_STAR_featureCounts_csv, UPORDOWN=ifelse(res_RAvsLA_STAR_featureCounts_csv$padj>0.05,"non",ifelse(res_RAvsLA_STAR_featureCounts_csv$log2FoldChange < 0,"LA","RA")))
res_RAvsLA_STAR_featureCounts_volcanoplot <- ggplot(res_RAvsLA_STAR_featureCounts_csv_mut, aes(log2FoldChange, -log(padj,10)))  +
  geom_point(aes(col = UPORDOWN), size = 0.8) +
  xlab(expression("log"[2]*"FC")) + 
  ylab(expression("-log"[10]*"P-adj")) +
  scale_color_manual(values = c("#D881BE", "lightgrey", "#575689")) +
  guides(colour = guide_legend(override.aes = list(size=1.5))) + theme(panel.background = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
res_RAvsLA_STAR_featureCounts_volcanoplot + theme(
  axis.title.x = element_text(size = 15), # Increase x-axis label size
  axis.title.y = element_text(size = 15), # Increase y-axis label size
  axis.text = element_text(size = 12)     # Increase axis text size
)

### AF vs SR
res_AFvsSR_STAR_featureCounts <- DESeq2::results(dds_STAR_featureCounts_protein_coding_filtered, contrast=c("Condition","AF","SR"), cooksCutoff=FALSE, independentFiltering=FALSE)
head(res_AFvsSR_STAR_featureCounts)
## The following code outputs the set of AF vs SR DEGs (as shown in 📎 Supplementary Table 4).
write.csv(res_AFvsSR_STAR_featureCounts,"res_AFvsSR_STAR_featureCounts2.csv")
summary(res_AFvsSR_STAR_featureCounts)


# 📊 Figure S4b
res_AFvsSR_STAR_featureCounts_csv <- read.csv("res_AFvsSR_STAR_featureCounts2.csv")
res_AFvsSR_STAR_featureCounts_csv <- res_AFvsSR_STAR_featureCounts_csv[res_AFvsSR_STAR_featureCounts_csv$X != 'ENSG00000277836', ]
res_AFvsSR_STAR_featureCounts_csv_mut = mutate(res_AFvsSR_STAR_featureCounts_csv, UPORDOWN=ifelse(res_AFvsSR_STAR_featureCounts_csv$padj>0.05,"non",ifelse(res_AFvsSR_STAR_featureCounts_csv$log2FoldChange < 0,"SR","AF")))
res_AFvsSR_STAR_featureCounts_volcanoplot <- ggplot(res_AFvsSR_STAR_featureCounts_csv_mut, aes(log2FoldChange, -log(padj,10)))  +
  geom_point(aes(col = UPORDOWN), size = 0.8) +
  xlab(expression("log"[2]*"FC")) + 
  ylab(expression("-log"[10]*"P-adj")) +
  scale_color_manual(values = c("#F4997F", "lightgrey", "#43B9D1")) +  
  guides(colour = guide_legend(override.aes = list(size=1.5))) + theme(panel.background = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
res_AFvsSR_STAR_featureCounts_volcanoplot + theme(
  axis.title.x = element_text(size = 15), # Increase x-axis label size
  axis.title.y = element_text(size = 15), # Increase y-axis label size
  axis.text = element_text(size = 12)     # Increase axis text size
)



## Genes of interest - expression boxplots

SLC43A3_expr <- read.table("SLC43A3_for_Boxplot.txt", sep="\t",header=TRUE) #File used to define sample type
counts <- counts(dds_STAR_featureCounts_protein_coding_filtered, normalized = TRUE)



### NPPB 
# 📊 Figures S4g and 4
counts_df_NPPB <- counts['ENSG00000120937',]
names(counts_df_NPPB)
sample_type <- SLC43A3_expr$Sample_type
disease_status <- c("SR","SR","SR","SR","SR","SR","SR","SR","SR","SR","AF","AF","AF","AF","AF","AF","AF","AF","AF","AF")
anatomical_side <- c("LA","RA","LA","RA","LA","RA","LA","RA","LA","RA")
X<-unname(counts_df_NPPB)
NPPB_expr <- data.frame(Sample = names(counts_df_NPPB),
                        Sample_type = sample_type,
                        Disease_status = disease_status,
                        Anat_side = anatomical_side,
                        X = X)
NPPB_expr_expr_boxplot_doxplot <- ggplot(NPPB_expr, aes(x = Disease_status, y = X, fill=Disease_status))  +
  geom_boxplot(colour = "black", width=0.3) +
  geom_dotplot(binaxis = "y", stackdir='center', dotsize = 1.3) + scale_fill_manual(values = c("#F4997F","#43B9D1")) +
  scale_y_continuous(name = "Normalized reads") + 
  scale_x_discrete(name = "Group") + 
  theme_classic() + ggtitle("NPPB expression") + 
  theme(plot.title = element_text(hjust = 0.5,face="bold",size =18)) 
NPPB_expr_expr_boxplot_doxplot +
  theme(axis.text.x = element_text(size = 14), axis.title.y = element_text(size = 14)) +
  geom_signif(comparisons = list(c("AF", "SR")), map_signif_level=FALSE,annotations = c(" ") )


### VASH1
# 📊 Figure 4
counts_df_VASH1 <- counts['ENSG00000071246',]
names(counts_df_VASH1)
sample_type <- SLC43A3_expr$Sample_type
disease_status <- c("SR","SR","SR","SR","SR","SR","SR","SR","SR","SR","AF","AF","AF","AF","AF","AF","AF","AF","AF","AF")
anatomical_side <- c("LA","RA","LA","RA","LA","RA","LA","RA","LA","RA")
X<-unname(counts_df_VASH1)
VASH1_expr <- data.frame(Sample = names(counts_df_VASH1),
                         Sample_type = sample_type,
                         Disease_status = disease_status,
                         Anat_side = anatomical_side,
                         X = X)
VASH1_expr_boxplot_doxplot <- ggplot(VASH1_expr, aes(x = Disease_status, y = X, fill=Disease_status))  +
  geom_boxplot(colour = "black", width=0.3) +
  geom_dotplot(binaxis = "y", stackdir='center', dotsize = 1.3) + scale_fill_manual(values = c("#F4997F","#43B9D1")) +
  scale_y_continuous(name = "Normalized reads") + 
  scale_x_discrete(name = "Group") + 
  theme_classic() + ggtitle("VASH1 expression") + 
  theme(plot.title = element_text(hjust = 0.5,face="bold",size =18)) 
VASH1_expr_boxplot_doxplot +
  theme(axis.text.x = element_text(size = 14), axis.title.y = element_text(size = 14)) +
  geom_signif(comparisons = list(c("AF", "SR")), map_signif_level=FALSE,annotations = c(" ") )


### ANGPTL2
# 📊 Figure 4
counts_df_ANGPTL2 <- counts['ENSG00000136859',]
names(counts_df_ANGPTL2)
sample_type <- SLC43A3_expr$Sample_type
disease_status <- c("SR","SR","SR","SR","SR","SR","SR","SR","SR","SR","AF","AF","AF","AF","AF","AF","AF","AF","AF","AF")
anatomical_side <- c("LA","RA","LA","RA","LA","RA","LA","RA","LA","RA")
X<-unname(counts_df_ANGPTL2)
ANGPTL2_expr <- data.frame(Sample = names(counts_df_ANGPTL2),
                           Sample_type = sample_type,
                           Disease_status = disease_status,
                           Anat_side = anatomical_side,
                           X = X)

ANGPTL2_expr_boxplot_doxplot <- ggplot(ANGPTL2_expr, aes(x = Disease_status, y = X, fill=Disease_status))  +
  geom_boxplot(colour = "black", width=0.3) +
  geom_dotplot(binaxis = "y", stackdir='center', dotsize = 1.3) + scale_fill_manual(values = c("#F4997F","#43B9D1")) +
  scale_y_continuous(name = "Normalized reads") + 
  scale_x_discrete(name = "Group") + 
  theme_classic() + ggtitle("ANGPTL2 expression") + 
  theme(plot.title = element_text(hjust = 0.5,face="bold",size =18)) 
ANGPTL2_expr_boxplot_doxplot +
  theme(axis.text.x = element_text(size = 14), axis.title.y = element_text(size = 14)) +
  geom_signif(comparisons = list(c("AF", "SR")), map_signif_level=FALSE,annotations = c(" ") )


### STRN
# 📊 Figure 4
counts_df_STRN <- counts['ENSG00000115808',]
names(counts_df_STRN)
sample_type <- SLC43A3_expr$Sample_type
disease_status <- c("SR","SR","SR","SR","SR","SR","SR","SR","SR","SR","AF","AF","AF","AF","AF","AF","AF","AF","AF","AF")
anatomical_side <- c("LA","RA","LA","RA","LA","RA","LA","RA","LA","RA")
X<-unname(counts_df_STRN)
STRN_expr <- data.frame(Sample = names(counts_df_STRN),
                        Sample_type = sample_type,
                        Disease_status = disease_status,
                        Anat_side = anatomical_side,
                        X = X)
STRN_expr_boxplot_doxplot <- ggplot(STRN_expr, aes(x = Disease_status, y = X, fill=Disease_status))  +
  geom_boxplot(colour = "black", width=0.3) +
  geom_dotplot(binaxis = "y", stackdir='center', dotsize = 1.3) + scale_fill_manual(values = c("#F4997F","#43B9D1")) +
  scale_y_continuous(name = "Normalized reads") + 
  scale_x_discrete(name = "Group") + 
  theme_classic() + ggtitle("STRN expression") + 
  theme(plot.title = element_text(hjust = 0.5,face="bold",size =18)) 
STRN_expr_boxplot_doxplot +
  theme(axis.text.x = element_text(size = 14), axis.title.y = element_text(size = 14))+
  geom_signif(comparisons = list(c("AF", "SR")), map_signif_level=FALSE,annotations = c(" ") )


### LRRC32
# 📊 Figure 4
counts_df_LRRC32 <- counts['ENSG00000137507',]
names(counts_df_LRRC32)
sample_type <- SLC43A3_expr$Sample_type
disease_status <- c("SR","SR","SR","SR","SR","SR","SR","SR","SR","SR","AF","AF","AF","AF","AF","AF","AF","AF","AF","AF")
anatomical_side <- c("LA","RA","LA","RA","LA","RA","LA","RA","LA","RA")
X<-unname(counts_df_LRRC32)
LRRC32_expr <- data.frame(Sample = names(counts_df_LRRC32),
                          Sample_type = sample_type,
                          Disease_status = disease_status,
                          Anat_side = anatomical_side,
                          X = X)
LRRC32_expr_boxplot_doxplot <- ggplot(LRRC32_expr, aes(x = Disease_status, y = X, fill=Disease_status))  +
  geom_boxplot(colour = "black", width=0.3) +
  geom_dotplot(binaxis = "y", stackdir='center', dotsize = 0.9) + scale_fill_manual(values = c("#F4997F","#43B9D1")) +
  scale_y_continuous(name = "Normalized reads") + 
  scale_x_discrete(name = "Group") + 
  theme_classic() + ggtitle("LRRC32 expression") + 
  theme(plot.title = element_text(hjust = 0.5,face="bold",size =18)) 
LRRC32_expr_boxplot_doxplot +
  theme(axis.text.x = element_text(size = 14), axis.title.y = element_text(size = 14))+
  geom_signif(comparisons = list(c("AF", "SR")), map_signif_level=FALSE,annotations = c(" ") )

# 📊 Figure 3c
LRRC32_expr_boxplot_doxplot <- ggplot(LRRC32_expr, aes(x = Sample_type, y = X, fill = Sample_type))  +
  geom_boxplot(colour = "black", width=0.3, fill = "#F2FDFD")  +
  geom_dotplot(binaxis = "y", stackdir='center', dotsize = 1.5) + scale_fill_manual(values = c("#C96274","#E1D08C","#7A96B8","#9CC6B2")) +
  scale_y_continuous(name = "Normalized reads") + 
  scale_x_discrete(name = "Group") + 
  theme_classic() + ggtitle(NULL) + 
  theme(plot.title = element_text(hjust = 0.5,face="bold",size =18), axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12), axis.title.y = element_text(size = 14)) 
LRRC32_expr_boxplot_doxplot
LRRC32_expr_boxplot_doxplot + 
  geom_signif(comparisons = list(c("AFLA", "AFRA")), map_signif_level=FALSE, y_position = 7700, annotations = c("n.s."), textsize = 5 )+
  geom_signif(comparisons = list(c("AFLA", "SRLA")), map_signif_level=FALSE,  y_position = 8000, annotations = c("***"), textsize = 5 ) +
  geom_signif(comparisons = list(c("AFLA", "SRRA")), map_signif_level=FALSE,  y_position = 8300, annotations = c("n.s"), textsize = 5 ) 



### SLC43A3
# 📊 Figures 4 and 6c
counts_df_SLC43A3 <- counts['ENSG00000134802',]
names(counts_df_SLC43A3)
sample_type <- SLC43A3_expr$Sample_type
disease_status <- c("SR","SR","SR","SR","SR","SR","SR","SR","SR","SR","AF","AF","AF","AF","AF","AF","AF","AF","AF","AF")
anatomical_side <- c("LA","RA","LA","RA","LA","RA","LA","RA","LA","RA")
X<-unname(counts_df_SLC43A3)
SLC43A3_expr <- data.frame(Sample = names(counts_df_SLC43A3),
                           Sample_type = sample_type,
                           Disease_status = disease_status,
                           Anat_side = anatomical_side,
                           X = X)

SLC43A3_expr_boxplot_doxplot <- ggplot(SLC43A3_expr, aes(x = Disease_status, y = X, fill=Disease_status))  +
  geom_boxplot(colour = "black", width=0.3) +
  geom_dotplot(binaxis = "y", stackdir='center', dotsize = 0.9) + scale_fill_manual(values = c("#F4997F","#43B9D1")) +
  scale_y_continuous(name = "Normalized reads") + 
  scale_x_discrete(name = "Group") + 
  theme_classic() + ggtitle("SLC43A3 expression") + 
  theme(plot.title = element_text(hjust = 0.5,face="bold",size =18)) 
SLC43A3_expr_boxplot_doxplot +
  theme(axis.text.x = element_text(size = 14), axis.title.y = element_text(size = 14))+
  geom_signif(comparisons = list(c("AF", "SR")), map_signif_level=FALSE,annotations = c(" ") )

### FOXS1
counts_df_FOXS1 <- counts['ENSG00000179772',] #ENSG00000179772
names(counts_df_FOXS1)
sample_type <- SLC43A3_expr$Sample_type
disease_status <- c("SR","SR","SR","SR","SR","SR","SR","SR","SR","SR","AF","AF","AF","AF","AF","AF","AF","AF","AF","AF")
anatomical_side <- c("LA","RA","LA","RA","LA","RA","LA","RA","LA","RA")
X<-unname(counts_df_FOXS1)
FOXS1_expr <- data.frame(Sample = names(counts_df_FOXS1),
                         Sample_type = sample_type,
                         Disease_status = disease_status,
                         Anat_side = anatomical_side,
                         X = X)

FOXS1_expr_boxplot_doxplot <- ggplot(FOXS1_expr, aes(x = Disease_status, y = X, fill=Disease_status))  +
  geom_boxplot(colour = "black", width=0.3) +
  geom_dotplot(binaxis = "y", stackdir='center', dotsize = 1.3) + scale_fill_manual(values = c("#F4997F","#43B9D1")) +
  scale_y_continuous(name = "Normalized reads") + 
  scale_x_discrete(name = "Group") + 
  theme_classic() + ggtitle("FOXS1 expression") + 
  theme(plot.title = element_text(hjust = 0.5,face="bold",size =18)) 
FOXS1_expr_boxplot_doxplot +
  theme(axis.text.x = element_text(size = 14), axis.title.y = element_text(size = 14))+
  geom_signif(comparisons = list(c("AF", "SR")), map_signif_level=FALSE,annotations = c(" ") )

### SCX
# 📊 Figure 4
counts_df_SCX <- counts['ENSG00000260428',]
names(counts_df_SCX)
sample_type <- SLC43A3_expr$Sample_type
disease_status <- c("SR","SR","SR","SR","SR","SR","SR","SR","SR","SR","AF","AF","AF","AF","AF","AF","AF","AF","AF","AF")
anatomical_side <- c("LA","RA","LA","RA","LA","RA","LA","RA","LA","RA")
X<-unname(counts_df_SCX)
SCX_expr <- data.frame(Sample = names(counts_df_SCX),
                       Sample_type = sample_type,
                       Disease_status = disease_status,
                       Anat_side = anatomical_side,
                       X = X)

SCX_expr_boxplot_doxplot <- ggplot(SCX_expr, aes(x = Disease_status, y = X, fill=Disease_status))  +
  geom_boxplot(colour = "black", width=0.3) +
  geom_dotplot(binaxis = "y", stackdir='center', dotsize = 1.3) + scale_fill_manual(values = c("#F4997F","#43B9D1")) +
  scale_y_continuous(name = "Normalized reads") + 
  scale_x_discrete(name = "Group") + 
  theme_classic() + ggtitle("SCX expression") + 
  theme(plot.title = element_text(hjust = 0.5,face="bold",size =18)) 
SCX_expr_boxplot_doxplot +
  theme(axis.text.x = element_text(size = 14), axis.title.y = element_text(size = 14))+
  geom_signif(comparisons = list(c("AF", "SR")), map_signif_level=FALSE,annotations = c(" ") )


# 📊 Figure 3d
SCX_expr_boxplot_doxplot <- ggplot(SCX_expr, aes(x = Sample_type, y = X, fill = Sample_type))  +
  geom_boxplot(colour = "black", width=0.3, fill = "#F2FDFD")  +
  geom_dotplot(binaxis = "y", stackdir='center', dotsize = 1.5) + scale_fill_manual(values = c("#C96274","#E1D08C","#7A96B8","#9CC6B2")) +
  scale_y_continuous(name = "Normalized reads") + 
  scale_x_discrete(name = "Group") + 
  theme_classic() + ggtitle(NULL) + 
  theme(plot.title = element_text(hjust = 0.5,face="bold",size =18), axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12), axis.title.y = element_text(size = 14)) 
SCX_expr_boxplot_doxplot
SCX_expr_boxplot_doxplot + 
  geom_signif(comparisons = list(c("AFRA", "AFLA")), map_signif_level=FALSE, y_position = 275, annotations = c("n.s."), textsize = 5 )+
  geom_signif(comparisons = list(c("AFRA", "SRLA")), map_signif_level=FALSE,  y_position = 255, annotations = c("***"), textsize = 5 ) +
  geom_signif(comparisons = list(c("AFRA", "SRRA")), map_signif_level=FALSE,  y_position = 275, annotations = c("*"), textsize = 5 ) 


### MT1X
# 📊 Figure 4
counts_df_MT1X <- counts['ENSG00000187193',]
names(counts_df_MT1X)
sample_type <- SLC43A3_expr$Sample_type
disease_status <- c("SR","SR","SR","SR","SR","SR","SR","SR","SR","SR","AF","AF","AF","AF","AF","AF","AF","AF","AF","AF")
anatomical_side <- c("LA","RA","LA","RA","LA","RA","LA","RA","LA","RA")
X<-unname(counts_df_MT1X)
MT1X_expr <- data.frame(Sample = names(counts_df_MT1X),
                        Sample_type = sample_type,
                        Disease_status = disease_status,
                        Anat_side = anatomical_side,
                        X = X)
MT1X_expr_expr_boxplot_doxplot <- ggplot(MT1X_expr, aes(x = Disease_status, y = X, fill=Disease_status))  +
  geom_boxplot(colour = "black", width=0.3) +
  geom_dotplot(binaxis = "y", stackdir='center', dotsize = 1.3) + scale_fill_manual(values = c("#F4997F","#43B9D1")) +
  scale_y_continuous(name = "Normalized reads") + 
  scale_x_discrete(name = "Group") + 
  theme_classic() + ggtitle("MT1X expression") + 
  theme(plot.title = element_text(hjust = 0.5,face="bold",size =18)) 
MT1X_expr_expr_boxplot_doxplot +
  theme(axis.text.x = element_text(size = 14), axis.title.y = element_text(size = 14))+
  geom_signif(comparisons = list(c("AF", "SR")), map_signif_level=FALSE,annotations = c(" ") )

### RGS6
# 📊 Figure 4
counts_df_RGS6 <- counts['ENSG00000182732',]
names(counts_df_RGS6)
sample_type <- SLC43A3_expr$Sample_type
disease_status <- c("SR","SR","SR","SR","SR","SR","SR","SR","SR","SR","AF","AF","AF","AF","AF","AF","AF","AF","AF","AF")
anatomical_side <- c("LA","RA","LA","RA","LA","RA","LA","RA","LA","RA")
X<-unname(counts_df_RGS6)
RGS6_expr <- data.frame(Sample = names(counts_df_RGS6),
                        Sample_type = sample_type,
                        Disease_status = disease_status,
                        Anat_side = anatomical_side,
                        X = X)
RGS6_expr_expr_boxplot_doxplot <- ggplot(RGS6_expr, aes(x = Disease_status, y = X, fill=Disease_status))  +
  geom_boxplot(colour = "black", width=0.3) +
  geom_dotplot(binaxis = "y", stackdir='center', dotsize = 1.3) + scale_fill_manual(values = c("#F4997F","#43B9D1")) +
  scale_y_continuous(name = "Normalized reads") + 
  scale_x_discrete(name = "Group") + 
  theme_classic() + ggtitle("RGS6 expression") + 
  theme(plot.title = element_text(hjust = 0.5,face="bold",size =18)) 
RGS6_expr_expr_boxplot_doxplot+
  geom_signif(comparisons = list(c("AF", "SR")), map_signif_level=FALSE,annotations = c(" ") )


### ANGPT1
# 📊 Figure 4
counts_df_ANGPT1 <- counts['ENSG00000154188',]
names(counts_df_ANGPT1)
sample_type <- SLC43A3_expr$Sample_type
disease_status <- c("SR","SR","SR","SR","SR","SR","SR","SR","SR","SR","AF","AF","AF","AF","AF","AF","AF","AF","AF","AF")
anatomical_side <- c("LA","RA","LA","RA","LA","RA","LA","RA","LA","RA")
X<-unname(counts_df_ANGPT1)
ANGPT1_expr <- data.frame(Sample = names(counts_df_ANGPT1),
                          Sample_type = sample_type,
                          Disease_status = disease_status,
                          Anat_side = anatomical_side,
                          X = X)
ANGPT1_expr_expr_boxplot_doxplot <- ggplot(ANGPT1_expr, aes(x = Disease_status, y = X, fill=Disease_status))  +
  geom_boxplot(colour = "black", width=0.3) +
  geom_dotplot(binaxis = "y", stackdir='center', dotsize = 1) + scale_fill_manual(values = c("#F4997F","#43B9D1")) +
  scale_y_continuous(name = "Normalized reads") + 
  scale_x_discrete(name = "Group") + 
  theme_classic() + ggtitle("ANGPT1 expression") + 
  theme(plot.title = element_text(hjust = 0.5,face="bold",size =18)) 
ANGPT1_expr_expr_boxplot_doxplot

### BLM
# 📊 Figure 4
counts_df_BLM <- counts['ENSG00000197299',]
names(counts_df_BLM)
sample_type <- SLC43A3_expr$Sample_type
disease_status <- c("SR","SR","SR","SR","SR","SR","SR","SR","SR","SR","AF","AF","AF","AF","AF","AF","AF","AF","AF","AF")
anatomical_side <- c("LA","RA","LA","RA","LA","RA","LA","RA","LA","RA")
X<-unname(counts_df_BLM)
BLM_expr <- data.frame(Sample = names(counts_df_BLM),
                         Sample_type = sample_type,
                         Disease_status = disease_status,
                         Anat_side = anatomical_side,
                         X = X)
BLM_expr_expr_boxplot_doxplot <- ggplot(BLM_expr, aes(x = Disease_status, y = X, fill=Disease_status))  +
  geom_boxplot(colour = "black", width=0.3) +
  geom_dotplot(binaxis = "y", stackdir='center', dotsize = 1.3) + scale_fill_manual(values = c("#F4997F","#43B9D1")) +
  scale_y_continuous(name = "Normalized reads") + 
  scale_x_discrete(name = "Group") + 
  theme_classic() + ggtitle("BLM expression") + 
  theme(plot.title = element_text(hjust = 0.5,face="bold",size =18)) 
BLM_expr_expr_boxplot_doxplot+
  geom_signif(comparisons = list(c("AF", "SR")), map_signif_level=FALSE,annotations = c(" ") )


### NR2F1
# 📊 Figure 4
counts_df_NR2F1 <- counts['ENSG00000175745',]
names(counts_df_NR2F1)
sample_type <- SLC43A3_expr$Sample_type
disease_status <- c("SR","SR","SR","SR","SR","SR","SR","SR","SR","SR","AF","AF","AF","AF","AF","AF","AF","AF","AF","AF")
anatomical_side <- c("LA","RA","LA","RA","LA","RA","LA","RA","LA","RA")
X<-unname(counts_df_NR2F1)
NR2F1_expr <- data.frame(Sample = names(counts_df_NR2F1),
                       Sample_type = sample_type,
                       Disease_status = disease_status,
                       Anat_side = anatomical_side,
                       X = X)
NR2F1_expr_expr_boxplot_doxplot <- ggplot(NR2F1_expr, aes(x = Disease_status, y = X, fill=Disease_status))  +
  geom_boxplot(colour = "black", width=0.3) +
  geom_dotplot(binaxis = "y", stackdir='center', dotsize = 1.3) + scale_fill_manual(values = c("#F4997F","#43B9D1")) +
  scale_y_continuous(name = "Normalized reads") + 
  scale_x_discrete(name = "Group") + 
  theme_classic() + ggtitle("NR2F1 expression") + 
  theme(plot.title = element_text(hjust = 0.5,face="bold",size =18)) 
NR2F1_expr_expr_boxplot_doxplot+
  geom_signif(comparisons = list(c("AF", "SR")), map_signif_level=FALSE,annotations = c(" ") )


### MYLK2
# 📊 Figure 6d
counts_df_MYLK2 <- counts['ENSG00000101306',]
names(counts_df_MYLK2)
sample_type <- SLC43A3_expr$Sample_type
disease_status <- c("SR","SR","SR","SR","SR","SR","SR","SR","SR","SR","AF","AF","AF","AF","AF","AF","AF","AF","AF","AF")
anatomical_side <- c("LA","RA","LA","RA","LA","RA","LA","RA","LA","RA")
X<-unname(counts_df_MYLK2)
MYLK2_expr <- data.frame(Sample = names(counts_df_MYLK2),
                         Sample_type = sample_type,
                         Disease_status = disease_status,
                         Anat_side = anatomical_side,
                         X = X)
MYLK2_expr_expr_boxplot_doxplot <- ggplot(MYLK2_expr, aes(x = Disease_status, y = X, fill=Disease_status))  +
  geom_boxplot(colour = "black", width=0.3) +
  geom_dotplot(binaxis = "y", stackdir='center', dotsize = 1.3) + scale_fill_manual(values = c("#F4997F","#43B9D1")) +
  scale_y_continuous(name = "Normalized reads") + 
  scale_x_discrete(name = "Group") + 
  theme_classic() + ggtitle("MYLK2 expression") + 
  theme(plot.title = element_text(hjust = 0.5,face="bold",size =18)) 
MYLK2_expr_expr_boxplot_doxplot+
  theme(axis.text.x = element_text(size = 14), axis.title.y = element_text(size = 14))+
  geom_signif(comparisons = list(c("AF", "SR")), map_signif_level=FALSE,annotations = c(" ") )


### LRRC4B
# 📊 Figures 7b and S7b
counts_df_LRRC4B <- counts['ENSG00000131409',]
names(counts_df_LRRC4B)
sample_type <- SLC43A3_expr$Sample_type
disease_status <- c("SR","SR","SR","SR","SR","SR","SR","SR","SR","SR","AF","AF","AF","AF","AF","AF","AF","AF","AF","AF")
anatomical_side <- c("LA","RA","LA","RA","LA","RA","LA","RA","LA","RA")
X<-unname(counts_df_LRRC4B)
LRRC4B_expr <- data.frame(Sample = names(counts_df_LRRC4B),
                          Sample_type = sample_type,
                          Disease_status = disease_status,
                          Anat_side = anatomical_side,
                          X = X)
LRRC4B_expr_boxplot_doxplot <- ggplot(LRRC4B_expr, aes(x = Disease_status, y = X, fill=Disease_status))  +
  geom_boxplot(colour = "black", width=0.3) +
  geom_dotplot(binaxis = "y", stackdir='center', dotsize = 1.2) + scale_fill_manual(values = c("#F4997F","#43B9D1")) +
  scale_y_continuous(name = "Normalized reads") + 
  scale_x_discrete(name = "Group") + 
  theme_classic() + ggtitle("LRRC4B expression") + 
  theme(plot.title = element_text(hjust = 0.5,face="bold",size =18)) 
LRRC4B_expr_boxplot_doxplot +
  theme(axis.text.x = element_text(size = 14), axis.title.y = element_text(size = 14))+
  geom_signif(comparisons = list(c("AF", "SR")), map_signif_level=FALSE,annotations = c(" ") )


### ZNF727
# 📊 Figure 7d
counts_df_ZNF727 <- counts['ENSG00000164107',]
names(counts_df_ZNF727)
sample_type <- SLC43A3_expr$Sample_type
disease_status <- c("SR","SR","SR","SR","SR","SR","SR","SR","SR","SR","AF","AF","AF","AF","AF","AF","AF","AF","AF","AF")
anatomical_side <- c("LA","RA","LA","RA","LA","RA","LA","RA","LA","RA")
X<-unname(counts_df_ZNF727)
ZNF727_expr <- data.frame(Sample = names(counts_df_ZNF727),
                          Sample_type = sample_type,
                          Disease_status = disease_status,
                          Anat_side = anatomical_side,
                          X = X)
ZNF727_expr_expr_boxplot_doxplot <- ggplot(ZNF727_expr, aes(x = Disease_status, y = X, fill=Disease_status))  +
  geom_boxplot(colour = "black", width=0.3) +
  geom_dotplot(binaxis = "y", stackdir='center', dotsize = 1.5) + scale_fill_manual(values = c("#F4997F","#43B9D1")) +
  scale_y_continuous(name = "Normalized reads") + 
  scale_x_discrete(name = "Group") + 
  theme_classic() + ggtitle("ZNF727 expression") + 
  theme(plot.title = element_text(hjust = 0.5,face="bold",size =18)) 

ZNF727_expr_expr_boxplot_doxplot+
  geom_signif(comparisons = list(c("AF", "SR")), map_signif_level=FALSE,annotations = c(" ") )


### COL23A1
# 📊 Figure S4f
counts_df_COL23A1 <- counts['ENSG00000050767',]
names(counts_df_COL23A1)
sample_type <- SLC43A3_expr$Sample_type
disease_status <- c("SR","SR","SR","SR","SR","SR","SR","SR","SR","SR","AF","AF","AF","AF","AF","AF","AF","AF","AF","AF")
anatomical_side <- c("LA","RA","LA","RA","LA","RA","LA","RA","LA","RA")
X<-unname(counts_df_COL23A1)
COL23A1_expr <- data.frame(Sample = names(counts_df_COL23A1),
                           Sample_type = sample_type,
                           Disease_status = disease_status,
                           Anat_side = anatomical_side,
                           X = X)
COL23A1_expr_boxplot_doxplot <- ggplot(COL23A1_expr, aes(x = Disease_status, y = X, fill=Disease_status))  +
  geom_boxplot(colour = "black", width=0.3, fill = "#F2FDFD")  +
  geom_dotplot(binaxis = "y", stackdir='center', dotsize = 1.3) + scale_fill_manual(values = c("#F4997F","#43B9D1")) +
  scale_y_continuous(name = "Normalized reads") + 
  scale_x_discrete(name = "Group") + 
  theme_classic() + ggtitle("COL23A1 expression") + 
  theme(plot.title = element_text(hjust = 0.5,face="bold",size =18)) 
COL23A1_expr_boxplot_doxplot +
  theme(axis.text.x = element_text(size = 14), axis.title.y = element_text(size = 14)) +
  geom_signif(comparisons = list(c("AF", "SR")), map_signif_level=FALSE,annotations = c("*") )


### SUMF1
# 📊 Figure S7a
counts_df_SUMF1 <- counts['ENSG00000144455',]
names(counts_df_SUMF1)
sample_type <- SLC43A3_expr$Sample_type
disease_status <- c("SR","SR","SR","SR","SR","SR","SR","SR","SR","SR","AF","AF","AF","AF","AF","AF","AF","AF","AF","AF")
anatomical_side <- c("LA","RA","LA","RA","LA","RA","LA","RA","LA","RA")
X<-unname(counts_df_SUMF1)
SUMF1_expr <- data.frame(Sample = names(counts_df_SUMF1),
                          Sample_type = sample_type,
                          Disease_status = disease_status,
                          Anat_side = anatomical_side,
                          X = X)
SUMF1_expr_expr_boxplot_doxplot <- ggplot(MNT_expr, aes(x = Disease_status, y = X, fill=Disease_status))  +
  geom_boxplot(colour = "black", width=0.3) +
  geom_dotplot(binaxis = "y", stackdir='center', dotsize = 1.5) + scale_fill_manual(values = c("#F4997F","#43B9D1")) +
  scale_y_continuous(name = "Normalized reads") + 
  scale_x_discrete(name = "Group") + 
  theme_classic() + ggtitle("SUMF1 expression") + 
  theme(plot.title = element_text(hjust = 0.5,face="bold",size =18)) 
SUMF1_expr_expr_boxplot_doxplot+
  geom_signif(comparisons = list(c("AF", "SR")), map_signif_level=FALSE,annotations = c(" ") )


### MNT
# 📊 Figure S7a
counts_df_MNT <- counts['ENSG00000070444',]
names(counts_df_MNT)
sample_type <- SLC43A3_expr$Sample_type
disease_status <- c("SR","SR","SR","SR","SR","SR","SR","SR","SR","SR","AF","AF","AF","AF","AF","AF","AF","AF","AF","AF")
anatomical_side <- c("LA","RA","LA","RA","LA","RA","LA","RA","LA","RA")
X<-unname(counts_df_SUMF1)
MNT_expr <- data.frame(Sample = names(counts_df_MNT),
                         Sample_type = sample_type,
                         Disease_status = disease_status,
                         Anat_side = anatomical_side,
                         X = X)
MNT_expr_expr_boxplot_doxplot <- ggplot(MNT_expr, aes(x = Disease_status, y = X, fill=Disease_status))  +
  geom_boxplot(colour = "black", width=0.3) +
  geom_dotplot(binaxis = "y", stackdir='center', dotsize = 1.5) + scale_fill_manual(values = c("#F4997F","#43B9D1")) +
  scale_y_continuous(name = "Normalized reads") + 
  scale_x_discrete(name = "Group") + 
  theme_classic() + ggtitle("MNT expression") + 
  theme(plot.title = element_text(hjust = 0.5,face="bold",size =18)) 
MNT_expr_expr_boxplot_doxplot+
  geom_signif(comparisons = list(c("AF", "SR")), map_signif_level=FALSE,annotations = c(" ") )


### SYT3
# 📊 Figure S7b
counts_df_SYT3 <- counts['ENSG00000070444',]
names(counts_df_SYT3)
sample_type <- SLC43A3_expr$Sample_type
disease_status <- c("SR","SR","SR","SR","SR","SR","SR","SR","SR","SR","AF","AF","AF","AF","AF","AF","AF","AF","AF","AF")
anatomical_side <- c("LA","RA","LA","RA","LA","RA","LA","RA","LA","RA")
X<-unname(counts_df_SUMF1)
SYT3_expr <- data.frame(Sample = names(counts_df_MNT),
                       Sample_type = sample_type,
                       Disease_status = disease_status,
                       Anat_side = anatomical_side,
                       X = X)
SYT3_expr_expr_boxplot_doxplot <- ggplot(SYT3_expr, aes(x = Disease_status, y = X, fill=Disease_status))  +
  geom_boxplot(colour = "black", width=0.3) +
  geom_dotplot(binaxis = "y", stackdir='center', dotsize = 1.5) + scale_fill_manual(values = c("#F4997F","#43B9D1")) +
  scale_y_continuous(name = "Normalized reads") + 
  scale_x_discrete(name = "Group") + 
  theme_classic() + ggtitle("SYT3 expression") + 
  theme(plot.title = element_text(hjust = 0.5,face="bold",size =18)) 
SYT3_expr_expr_boxplot_doxplot+
  geom_signif(comparisons = list(c("AF", "SR")), map_signif_level=FALSE,annotations = c(" ") )


### SHANK1
# 📊 Figure S7b
counts_df_SHANK1 <- counts['ENSG00000006468',]
names(counts_df_SHANK1)
sample_type <- SLC43A3_expr$Sample_type
disease_status <- c("SR","SR","SR","SR","SR","SR","SR","SR","SR","SR","AF","AF","AF","AF","AF","AF","AF","AF","AF","AF")
anatomical_side <- c("LA","RA","LA","RA","LA","RA","LA","RA","LA","RA")
X<-unname(counts_df_SHANK1)
SHANK1_expr <- data.frame(Sample = names(counts_df_SHANK1),
                          Sample_type = sample_type,
                          Disease_status = disease_status,
                          Anat_side = anatomical_side,
                          X = X)
SHANK1_expr_boxplot_doxplot <- ggplot(SHANK1_expr, aes(x = Disease_status, y = X, fill=Disease_status))  +
  geom_boxplot(colour = "black", width=0.3) +
  geom_dotplot(binaxis = "y", stackdir='center', dotsize = 1.2) + scale_fill_manual(values = c("#F4997F","#43B9D1")) +
  scale_y_continuous(name = "Normalized reads") + 
  scale_x_discrete(name = "Group") + 
  theme_classic() + ggtitle("SHANK1 expression") + 
  theme(plot.title = element_text(hjust = 0.5,face="bold",size =18)) 
SHANK1_expr_boxplot_doxplot +
  theme(axis.text.x = element_text(size = 14), axis.title.y = element_text(size = 14))+
  geom_signif(comparisons = list(c("AF", "SR")), map_signif_level=FALSE,annotations = c(" ") )


### SETMAR
# 📊 Figure S7a
counts_df_SETMAR <- counts['ENSG00000170364',]
names(counts_df_SETMAR)
sample_type <- SLC43A3_expr$Sample_type
disease_status <- c("SR","SR","SR","SR","SR","SR","SR","SR","SR","SR","AF","AF","AF","AF","AF","AF","AF","AF","AF","AF")
anatomical_side <- c("LA","RA","LA","RA","LA","RA","LA","RA","LA","RA")
X<-unname(counts_df_SETMAR)
SETMAR_expr <- data.frame(Sample = names(counts_df_SETMAR),
                          Sample_type = sample_type,
                          Disease_status = disease_status,
                          Anat_side = anatomical_side,
                          X = X)
SETMAR_expr_boxplot_doxplot <- ggplot(SETMAR_expr, aes(x = Disease_status, y = X, fill=Disease_status))  +
  geom_boxplot(colour = "black", width=0.3) +
  geom_dotplot(binaxis = "y", stackdir='center', dotsize = 1.2) + scale_fill_manual(values = c("#F4997F","#43B9D1")) +
  scale_y_continuous(name = "Normalized reads") + 
  scale_x_discrete(name = "Group") + 
  theme_classic() + ggtitle("SETMAR expression") + 
  theme(plot.title = element_text(hjust = 0.5,face="bold",size =18)) 
SETMAR_expr_boxplot_doxplot +
  theme(axis.text.x = element_text(size = 14), axis.title.y = element_text(size = 14))+
  geom_signif(comparisons = list(c("AF", "SR")), map_signif_level=FALSE,annotations = c(" ") )



### WNT5A
# 📊 Figure 3e
counts_df_WNT5A <- counts['ENSG00000114251',]
names(counts_df_WNT5A)
sample_type <- SLC43A3_expr$Sample_type
disease_status <- c("SR","SR","SR","SR","SR","SR","SR","SR","SR","SR","AF","AF","AF","AF","AF","AF","AF","AF","AF","AF")
anatomical_side <- c("LA","RA","LA","RA","LA","RA","LA","RA","LA","RA")
X<-unname(counts_df_WNT5A)
WNT5A_expr <- data.frame(Sample = names(counts_df_WNT5A),
                          Sample_type = sample_type,
                          Disease_status = disease_status,
                          Anat_side = anatomical_side,
                          X = X)
WNT5A_expr_boxplot_doxplot <- ggplot(WNT5A_expr, aes(x = Sample_type, y = X, fill = Sample_type))  +
  geom_boxplot(colour = "black", width=0.3, fill = "#F2FDFD")  +
  geom_dotplot(binaxis = "y", stackdir='center', dotsize = 2.3) + scale_fill_manual(values = c("#C96274","#E1D08C","#7A96B8","#9CC6B2")) +
  scale_y_continuous(name = "Normalized reads") + 
  scale_x_discrete(name = "Group") + 
  theme_classic() + ggtitle(NULL) + 
  theme(plot.title = element_text(hjust = 0.5,face="bold",size =17), axis.text.x = element_text(size = 17), axis.text.y = element_text(size = 12), axis.title.y = element_text(size = 17)) 
WNT5A_expr_boxplot_doxplot
WNT5A_expr_boxplot_doxplot + 
  geom_signif(comparisons = list(c("SRRA", "SRLA")), map_signif_level=FALSE, y_position = 2900, annotations = c("***"), textsize = 8 )+
  geom_signif(comparisons = list(c("SRRA", "AFRA")), map_signif_level=FALSE, y_position = 3200, annotations = c("**"), textsize = 8 ) +
  geom_signif(comparisons = list(c("SRRA", "AFLA")), map_signif_level=FALSE, y_position = 3500, annotations = c("***"), textsize = 8 ) 



################################################################################
################################################################################
#########################                     ##################################
#########################   4-group analyses  ##################################
#########################                     ##################################
################################################################################
################################################################################

## AFRA vs SRRA
## AFLA vs SRLA

coldata_4groups <- read.delim("coldata_4groups.txt")

dds_STAR_featureCounts_protein_coding_filtered_4groups <- dds_STAR_featureCounts_protein_coding_filtered
colData(dds_STAR_featureCounts_protein_coding_filtered_4groups) <- DataFrame(coldata_4groups)
colData(dds_STAR_featureCounts_protein_coding_filtered_4groups)

#Set SRRA as the reference level
colData(dds_STAR_featureCounts_protein_coding_filtered_4groups)$ConSide <- factor(colData(dds_STAR_featureCounts_protein_coding_filtered_4groups)$ConSide)
colData(dds_STAR_featureCounts_protein_coding_filtered_4groups)$ConSide <- relevel(colData(dds_STAR_featureCounts_protein_coding_filtered_4groups)$ConSide, ref="SRRA")
colData(dds_STAR_featureCounts_protein_coding_filtered_4groups)$ConSide


dds_STAR_featureCounts_protein_coding_filtered_4groups
dds_STAR_featureCounts_protein_coding_filtered_4groups <- DESeqDataSet(dds_STAR_featureCounts_protein_coding_filtered_4groups, design = ~ ConSide)
colnames(dds_STAR_featureCounts_protein_coding_filtered_4groups) <- colData(dds_STAR_featureCounts_protein_coding_filtered_4groups)$names



### EXPLORATORY ANALYSIS AND VISUALIZATION ###

# PCA
# 📊 Figure S4a

PCA_4groups <- DESeq2::plotPCA(vsd_STAR_featureCounts_protein_coding_filtered_4groups, intgroup = "ConSide")
PCA_4groups <- PCA_4groups + scale_color_manual(values=c("#91C7B1", "#D85B74", "#E3D081", "#7298BC"))
PCA_4groups+
  theme(
    axis.title.x = element_text(size = 15), 
    axis.title.y = element_text(size = 15), 
    axis.text = element_text(size = 12)     
  )

# 
PCA_4groups <- PCA_4groups +
  geom_point(size = 5) +  
  theme_minimal() +       
  theme(panel.background = element_blank(),  
        plot.background = element_blank(),   
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),  
        axis.line = element_line())


### DIFFERENTIAL EXPRESSION ANALYSIS ###

dds_STAR_featureCounts_protein_coding_filtered_4groups <- DESeq2::DESeq(dds_STAR_featureCounts_protein_coding_filtered_4groups, minReplicatesForReplace=Inf)



# SR-RA VS AF-RA
## The following code outputs the set of SR-RA vs AF-RA DEGs (as shown in 📎 Supplementary Table 4).
res_SRRAvsAFRA_STAR_featureCounts <- DESeq2::results(dds_STAR_featureCounts_protein_coding_filtered_4groups, contrast=c("ConSide","SRRA","AFRA"), cooksCutoff=FALSE, independentFiltering=FALSE)
head(res_SRRAvsAFRA_STAR_featureCounts)
write.csv(res_SRRAvsAFRA_STAR_featureCounts,"res_SRRAvsAFRA_STAR_featureCounts.csv")
summary(res_SRRAvsAFRA_STAR_featureCounts)

res_SRRAvsAFRA_STAR_featureCounts_csv <- read.csv("res_SRRAvsAFRA_STAR_featureCounts.csv")
res_SRRAvsAFRA_STAR_featureCounts_csv <- res_SRRAvsAFRA_STAR_featureCounts_csv[res_SRRAvsAFRA_STAR_featureCounts_csv$X != 'ENSG00000277836', ] # remove outliers for visualization
res_SRRAvsAFRA_STAR_featureCounts_csv <- res_SRRAvsAFRA_STAR_featureCounts_csv[res_SRRAvsAFRA_STAR_featureCounts_csv$X != 'ENSG00000131737', ]

res_SRRAvsAFRA_STAR_featureCounts_csv_mut = mutate(res_SRRAvsAFRA_STAR_featureCounts_csv, UPORDOWN=ifelse(res_SRRAvsAFRA_STAR_featureCounts_csv$padj>0.05,"non",ifelse(res_SRRAvsAFRA_STAR_featureCounts_csv$log2FoldChange < 0,"AF-RA","SR-RA")))
res_SRRAvsAFRA_STAR_featureCounts_volcanoplot <- ggplot(res_SRRAvsAFRA_STAR_featureCounts_csv_mut, aes(log2FoldChange, -log(padj,10)))  +
  geom_point(aes(col = UPORDOWN), size = 0.8) +
  xlab(expression("log"[2]*"FC")) + 
  ylab(expression("-log"[10]*"P-adj")) +
  scale_color_manual(values = c("#DFCD7F", "lightgrey", "#91C8B2")) +
  guides(colour = guide_legend(override.aes = list(size=1.5))) + theme(panel.background = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
res_SRRAvsAFRA_STAR_featureCounts_volcanoplot + theme(
  axis.title.x = element_text(size = 15), 
  axis.title.y = element_text(size = 15), 
  axis.text = element_text(size = 12)     
) # 📊 Figure S4c


# SR-LA VS AF-LA
## The following code outputs the set of SR-LA vs AF-LA DEGs (as shown in 📎 Supplementary Table 4).
res_SRLAvsAFLA_STAR_featureCounts <- DESeq2::results(dds_STAR_featureCounts_protein_coding_filtered_4groups, contrast=c("ConSide","SRLA","AFLA"), cooksCutoff=FALSE, independentFiltering=FALSE)
head(res_SRLAvsAFLA_STAR_featureCounts)
write.csv(res_SRLAvsAFLA_STAR_featureCounts,"res_SRLAvsAFLA_STAR_featureCounts.csv")
summary(res_SRLAvsAFLA_STAR_featureCounts)

res_SRLAvsAFLA_STAR_featureCounts_csv <- read.csv("res_SRLAvsAFLA_STAR_featureCounts.csv")
res_SRLAvsAFLA_STAR_featureCounts_csv <- res_SRLAvsAFLA_STAR_featureCounts_csv[res_SRLAvsAFLA_STAR_featureCounts_csv$X != 'ENSG00000277836', ]  # remove outliers for visualization
res_SRLAvsAFLA_STAR_featureCounts_csv <- res_SRLAvsAFLA_STAR_featureCounts_csv[res_SRLAvsAFLA_STAR_featureCounts_csv$X != 'ENSG00000131737', ]

res_SRLAvsAFLA_STAR_featureCounts_csv_mut = mutate(res_SRLAvsAFLA_STAR_featureCounts_csv, UPORDOWN=ifelse(res_SRLAvsAFLA_STAR_featureCounts_csv$padj>0.05,"non",ifelse(res_SRLAvsAFLA_STAR_featureCounts_csv$log2FoldChange < 0,"AF-LA","SR-LA")))
res_SRLAvsAFLA_STAR_featureCounts_volcanoplot <- ggplot(res_SRLAvsAFLA_STAR_featureCounts_csv_mut, aes(log2FoldChange, -log(padj,10)))  +
  geom_point(aes(col = UPORDOWN), size = 0.8) +
  xlab(expression("log"[2]*"FC")) + 
  ylab(expression("-log"[10]*"P-adj")) +
  scale_color_manual(values = c("#D85C74", "lightgrey", "#7297BD")) +
  guides(colour = guide_legend(override.aes = list(size=1.5))) + theme(panel.background = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
res_SRLAvsAFLA_STAR_featureCounts_volcanoplot + theme(
  axis.title.x = element_text(size = 15), 
  axis.title.y = element_text(size = 15), 
  axis.text = element_text(size = 12)     
) # 📊 Figure S4c



## Overlap of genes found to be differentially expressed in each anatomical side 

res_SRRAvsAFRA_STAR_featureCounts_csv_sign <- res_SRRAvsAFRA_STAR_featureCounts_csv[res_SRRAvsAFRA_STAR_featureCounts_csv$padj < 0.05, ]

res_SRLAvsAFLA_STAR_featureCounts_csv_sign <- res_SRLAvsAFLA_STAR_featureCounts_csv[res_SRLAvsAFLA_STAR_featureCounts_csv$padj < 0.05, ]

intersectt <- function (x, y) 
{
  if (is.null(x) || is.null(y)) 
    return(NULL)
  u <- as.vector(x)
  v <- as.vector(y)
  c(u[!duplicated(unclass(u)) & (match(u, v, 0L) > 0L)], v[numeric()])
}

Overlap_of_DEGs_betweenLAandRA <- intersectt(res_SRLAvsAFLA_STAR_featureCounts_csv_sign$X,res_SRRAvsAFRA_STAR_featureCounts_csv_sign$X)  
Overlap_of_DEGs_betweenLAandRA_df <- data.frame(Common_genes = Overlap_of_DEGs_betweenLAandRA)

# 📊 Figure S4c - Venn Diagram
Overlap_of_DEGs_betweenLAandRA_x <- list(
  A = res_SRLAvsAFLA_STAR_featureCounts_csv_sign$X, 
  B = res_SRRAvsAFRA_STAR_featureCounts_csv_sign$X)
Overlap_of_DEGs_betweenLAandRA_VennDiagram <- ggvenn(
  Overlap_of_DEGs_betweenLAandRA_x, 
  fill_color = c("#D882BD", "#585688" ),
  stroke_size = 0.5, set_name_size = 4, text_size = 6
)
Overlap_of_DEGs_betweenLAandRA_VennDiagram



### VISUALIZATION OF TFs EXPRESSION IN ATRIAL SAMPLES - Related to Figure 7 (TF enrichment analyses in DMRs) ###


### TPM representation

load("count_matrix.RData")
if (!all(rownames(count_matrix$counts) %in% count_matrix$annotation$GeneID)) {
  stop("Some genes in count_matrix are not found in annotation!")
}
missing_genes <- setdiff(rownames(count_matrix$counts), count_matrix$annotation$GeneID)
missing_genes
# All genes are matching 

# Ensure gene lengths are in kilobases
count_matrix$annotation$Length_kb <- count_matrix$annotation$Length / 1000

# Merge gene lengths into count_matrix based on gene_id
count_matrix_reads <- count_matrix$counts
count_matrix_reads <- as.data.frame(count_matrix_reads)
count_matrix_reads$gene_id <- rownames(count_matrix_reads)
merged_df <- merge(count_matrix_reads, count_matrix$annotation, by.x = "gene_id", by.y = "GeneID")

# Drop gene_id for calculations
rownames(merged_df) <- merged_df$gene_id
merged_df$gene_id <- NULL

# Calculate TPM
calculate_tpm <- function(counts, gene_lengths_kb) {
  # Normalize counts by gene length
  rate <- counts / gene_lengths_kb
  # Scale by total rate per sample
  tpm <- t(t(rate) / colSums(rate)) * 1e6
  return(tpm)
}

# Apply the function
gene_lengths_kb <- merged_df$Length_kb  # Extract gene lengths
raw_counts <- merged_df[, -c(26,25,24,23,22,21)]  # Exclude annotation columns

raw_counts/gene_lengths_kb

tpm_matrix <- calculate_tpm(raw_counts, gene_lengths_kb)

### Represent heatmap for DMR TF binding enrichment (selection from HOMER de novo TF enrichment analysis)

selected_genes_Hypermeth <- c("ENSG00000185551","ENSG00000137090","ENSG00000006468","ENSG00000169957","ENSG00000083817","ENSG00000152433","ENSG00000166823","ENSG00000164920",
                              "ENSG00000155090","ENSG00000175387","ENSG00000173275","ENSG00000136870") 

selected_genes_Hypometh <- c("ENSG00000196628","ENSG00000204366","ENSG00000140262","ENSG00000131061",
                             "ENSG00000117595", "ENSG00000128652", "ENSG00000147421", "ENSG00000129194", "ENSG00000111046", "ENSG00000177494", "ENSG00000102145", "ENSG00000116017")

coldata <- read.delim("coldata.txt")

subset_tpm_Hypermeth <- tpm_matrix[rownames(tpm_matrix) %in% selected_genes_Hypermeth, ]
subset_tpm_Hypometh <- tpm_matrix[rownames(tpm_matrix) %in% selected_genes_Hypometh, ]

disease_samples <- c("GSM3666581_AF_LA_Homo_sapiens.GRCh38_Male_SRR8715053-72_5_RNA-Seq.bam","GSM3666582_AF_RA_Homo_sapiens.GRCh38_Male_SRR8715073-82_5_RNA-Seq.bam",
                     "GSM3666583_AF_LA_Homo_sapiens.GRCh38_Male_SRR8715083-112_14_RNA-Seq.bam", "GSM3666584_AF_RA_Homo_sapiens.GRCh38_Male_SRR8715113-32_14_RNA-Seq.bam",
                     "GSM3666585_AF_LA_Homo_sapiens.GRCh38_Male_SRR8715133-52_25_RNA-Seq.bam","GSM3666586_AF_RA_Homo_sapiens.GRCh38_Male_SRR8715153-72_25_RNA-Seq.bam",
                     "GSM3666587_AF_LA_Homo_sapiens.GRCh38_Male_SRR8715173-92_28_RNA-Seq.bam","GSM3666588_AF_RA_Homo_sapiens.GRCh38_Male_SRR8715193-212_28_RNA-Seq.bam",
                     "GSM3666589_AF_LA_Homo_sapiens.GRCh38_Male_SRR8715213-32_31_RNA-Seq.bam","GSM3666590_AF_RA_Homo_sapiens.GRCh38_Male_SRR8715233-52_31_RNA-Seq.bam")
control_samples <- c("GSM3666571_SR_LA_Homo_sapiens.GRCh38_Male_SRR8714853-72_7_RNA-Seq.bam","GSM3666572_SR_RA_Homo_sapiens.GRCh38_Male_SRR8714873-92_7_RNA-Seq.bam",
                     "GSM3666573_SR_LA_Homo_sapiens.GRCh38_Male_SRR8714893-12_12b_RNA-Seq.bam", "GSM3666574_SR_RA_Homo_sapiens.GRCh38_Male_SRR8714913-32_12b_RNA-Seq.bam",
                     "GSM3666575_SR_LA_Homo_sapiens.GRCh38_Male_SRR8714933-52_20_RNA-Seq.bam","GSM3666576_SR_RA_Homo_sapiens.GRCh38_Male_SRR8714953-72_20_RNA-Seq.bam",
                     "GSM3666577_SR_LA_Homo_sapiens.GRCh38_Male_SRR8714973-92_27_RNA-Seq.bam","GSM3666578_SR_RA_Homo_sapiens.GRCh38_Male_SRR8714993-5012_27_RNA-Seq.bam",
                     "GSM3666579_SR_LA_Homo_sapiens.GRCh38_Male_SRR8715013-32_30_RNA-Seq.bam","GSM3666580_SR_RA_Homo_sapiens.GRCh38_Male_SRR8715033-52_30_RNA-Seq.bam")

median_disease_hyper <- apply(subset_tpm_Hypermeth[, disease_samples], 1, median, na.rm = TRUE)
median_control_hyper <- apply(subset_tpm_Hypermeth[, control_samples], 1, median, na.rm = TRUE)

median_disease_hypo <- apply(subset_tpm_Hypometh[, disease_samples], 1, median, na.rm = TRUE)
median_control_hypo <- apply(subset_tpm_Hypometh[, control_samples], 1, median, na.rm = TRUE)

median_matrix_hyper <- cbind(Control = median_control_hyper, Disease = median_disease_hyper)
median_matrix_hypo <- cbind(Control = median_control_hypo, Disease = median_disease_hypo)

log_median_matrix_hyper <- log2(median_matrix_hyper + 1)
log_median_matrix_hypo <- log2(median_matrix_hypo + 1)

row_order_hyper <- c("ENSG00000185551","ENSG00000137090","ENSG00000006468","ENSG00000169957","ENSG00000083817","ENSG00000152433","ENSG00000166823","ENSG00000164920",
                     "ENSG00000155090","ENSG00000175387","ENSG00000173275","ENSG00000136870")

row_order_hypo <- c("ENSG00000196628","ENSG00000204366","ENSG00000140262","ENSG00000131061",
                    "ENSG00000117595", "ENSG00000128652", "ENSG00000147421", "ENSG00000129194", "ENSG00000111046", "ENSG00000177494", "ENSG00000102145", "ENSG00000116017")

col_order <- c("Control","Disease")

log_median_matrix_ordered_hyper <- log_median_matrix_hyper[row_order_hyper,col_order]
log_median_matrix_ordered_hypo <- log_median_matrix_hypo[row_order_hypo,col_order]

# 📊 Figure 7a - Hypermethylated DMRs 
pheatmap(log_median_matrix_ordered_hyper,
         cluster_rows = F,
         cluster_cols = F,
         scale =  "none",
         main = "Median TPM (log2) in Disease vs Control",
         color = c("#7FFFD4","#76EEC6","#66CDAA","#458B74"))

# 📊 Figure 7a - Hypomethylated DMRs 
pheatmap(log_median_matrix_ordered_hypo,
         cluster_rows = F,
         cluster_cols = F,
         breaks = seq(0,6),
         scale =  "none",
         main = "Median TPM (log2) in Disease vs Control",
         color = long_palette)


### Represent heatmap for selection of AF-related TFs
# ZFHX3, PRRX1, NKX2-5, PITX2, TBX5, TBX3, SOX5, HAND2, NEURL, SHOX2, ETV1, CREB2

selected_genes <- c("ENSG00000140836", "ENSG00000116132", "ENSG00000183072", "ENSG00000164093", "ENSG00000089225", "ENSG00000135111", "ENSG00000134532", "ENSG00000164107",
                    "ENSG00000107954","ENSG00000168779","ENSG00000006468","ENSG00000128272")

subset_tpm <- tpm_matrix[rownames(tpm_matrix) %in% selected_genes, ]

LA_samples <- c("GSM3666581_AF_LA_Homo_sapiens.GRCh38_Male_SRR8715053-72_5_RNA-Seq.bam", "GSM3666583_AF_LA_Homo_sapiens.GRCh38_Male_SRR8715083-112_14_RNA-Seq.bam",
                "GSM3666585_AF_LA_Homo_sapiens.GRCh38_Male_SRR8715133-52_25_RNA-Seq.bam", "GSM3666587_AF_LA_Homo_sapiens.GRCh38_Male_SRR8715173-92_28_RNA-Seq.bam",
                "GSM3666589_AF_LA_Homo_sapiens.GRCh38_Male_SRR8715213-32_31_RNA-Seq.bam", "GSM3666571_SR_LA_Homo_sapiens.GRCh38_Male_SRR8714853-72_7_RNA-Seq.bam",
                "GSM3666573_SR_LA_Homo_sapiens.GRCh38_Male_SRR8714893-12_12b_RNA-Seq.bam", "GSM3666575_SR_LA_Homo_sapiens.GRCh38_Male_SRR8714933-52_20_RNA-Seq.bam",
                "GSM3666577_SR_LA_Homo_sapiens.GRCh38_Male_SRR8714973-92_27_RNA-Seq.bam", "GSM3666579_SR_LA_Homo_sapiens.GRCh38_Male_SRR8715013-32_30_RNA-Seq.bam")

RA_samples <- c("GSM3666582_AF_RA_Homo_sapiens.GRCh38_Male_SRR8715073-82_5_RNA-Seq.bam", "GSM3666584_AF_RA_Homo_sapiens.GRCh38_Male_SRR8715113-32_14_RNA-Seq.bam",
                "GSM3666586_AF_RA_Homo_sapiens.GRCh38_Male_SRR8715153-72_25_RNA-Seq.bam", "GSM3666588_AF_RA_Homo_sapiens.GRCh38_Male_SRR8715193-212_28_RNA-Seq.bam",
                "GSM3666590_AF_RA_Homo_sapiens.GRCh38_Male_SRR8715233-52_31_RNA-Seq.bam", "GSM3666572_SR_RA_Homo_sapiens.GRCh38_Male_SRR8714873-92_7_RNA-Seq.bam",
                "GSM3666574_SR_RA_Homo_sapiens.GRCh38_Male_SRR8714913-32_12b_RNA-Seq.bam", "GSM3666576_SR_RA_Homo_sapiens.GRCh38_Male_SRR8714953-72_20_RNA-Seq.bam",
                "GSM3666578_SR_RA_Homo_sapiens.GRCh38_Male_SRR8714993-5012_27_RNA-Seq.bam", "GSM3666580_SR_RA_Homo_sapiens.GRCh38_Male_SRR8715033-52_30_RNA-Seq.bam")

## AF vs SR

median_disease <- apply(subset_tpm[, disease_samples], 1, median, na.rm = TRUE)
median_control <- apply(subset_tpm[, control_samples], 1, median, na.rm = TRUE)

median_matrix <- cbind(Control = median_control, Disease = median_disease)

log_median_matrix <- log2(median_matrix + 1)


row_order <- c("ENSG00000140836", "ENSG00000116132", "ENSG00000183072", "ENSG00000164093", "ENSG00000089225", "ENSG00000135111", "ENSG00000134532", "ENSG00000164107",
               "ENSG00000107954","ENSG00000168779","ENSG00000006468","ENSG00000128272")
col_order <- c("Control","Disease")
log_median_matrix_ordered <- log_median_matrix[row_order,col_order]

# 📊 Figure 7c - AF vs SR
pheatmap(log_median_matrix_ordered,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         scale =  "none",
         main = "Median TPM (log2) in Disease vs Control",
         color = c("#7FFFD4","#76EEC6","#66CDAA","#458B74") )


## RA vs LA

median_LA <- apply(subset_tpm[, LA_samples], 1, median, na.rm = TRUE)
median_RA <- apply(subset_tpm[, RA_samples], 1, median, na.rm = TRUE)

median_matrix <- cbind(LA = median_LA, RA = median_RA)

log_median_matrix <- log2(median_matrix + 1)

row_order <- c("ENSG00000140836", "ENSG00000116132", "ENSG00000183072", "ENSG00000164093", "ENSG00000089225", "ENSG00000135111", "ENSG00000134532", "ENSG00000164107",
               "ENSG00000107954","ENSG00000168779","ENSG00000006468","ENSG00000128272")

col_order <- c("LA","RA")
log_median_matrix_ordered <- log_median_matrix[row_order,col_order]

# 📊 Figure 7c - RA vs LA
pheatmap(log_median_matrix_ordered,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         scale =  "none",
         main = "Median TPM (log2) in LA vs RA",
         color = c("#ecdff5","#C7A7E3","#9C6BC9","#6B3FA0") )


