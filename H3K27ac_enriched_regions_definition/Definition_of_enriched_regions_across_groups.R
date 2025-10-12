

library(DiffBind)
library(dplyr)
library(ggplot2)
library(reshape2)
library(Repitools)
library(ChIPpeakAnno)
library(nVennR)
library(clusterProfiler)




############################################################################
##### Exploration of peaks and correlation across individual samples #######
############################################################################

# 📊 Figures S1b and S1c

setwd("/Rodriguez_2025_AFepigenome/H3K27ac_enriched_regions_definition/")

AF_dataset_samples_details_Nov23 <- read.table("AF_dataset_libraries_peaks_metadata_2.txt", sep="\t",header=TRUE)

Peaks_called_per_sample_plot <- ggplot(AF_dataset_samples_details_Nov23, aes(x = factor(Sample,levels = c("25-AF-LA","28-AF-LA","47-AF-LA","31-AF-LA","2231-AF-LA","313-AF-LA","25-AF-RA","2231-AF-RA","31-AF-RA","313-AF-RA","18-AF-RA","28-AF-RA","12b-SR-LA","45-SR-LA","20-SR-LA","27-SR-LA","12b-SR-RA","27-SR-RA","33-SR-RA","20-SR-RA")), y = called_peaks, fill=Subgroup)) + geom_bar(stat="identity",color="black",size=0.3) + scale_y_continuous("Peaks called per sample") + scale_x_discrete("Atrial sample") + theme_classic() + theme(axis.text.x = element_text(angle=45, hjust = 1)) + scale_fill_manual(values = c("#D85B74","#E3D081","#7298BC","#91C7B1")) + geom_text(aes(label=called_peaks), position=position_dodge(width=0.9), vjust=-0.25, size=3.7)
Peaks_called_per_sample_plot # 📊 Figure S1b

# Calculate average per group 
df_avg <- AF_dataset_samples_details_Nov23 %>%
  group_by(Subgroup) %>%
  summarise(avg_called_peaks = mean(called_peaks, na.rm = TRUE))
print(df_avg)


#Explore peaks length per sample

setwd("/Rodriguez_2025_AFepigenome/ChIP-seq_basic_processing/")
vl0029_25_AF_RA <- as.data.frame(read.table("25_AF_LA_peaks.narrowPeak",header = FALSE, sep="\t",stringsAsFactors=FALSE))
vl0214_2231_AF_RA <- as.data.frame(read.table("2231_AF_RA_peaks.narrowPeak",header = FALSE, sep="\t",stringsAsFactors=FALSE))
vl0258_31_AF_RA <- as.data.frame(read.table("31_AF_RA_peaks.narrowPeak",header = FALSE, sep="\t",stringsAsFactors=FALSE))
vl0260_313_AF_RA <- as.data.frame(read.table("313_AF_RA_peaks.narrowPeak",header = FALSE, sep="\t",stringsAsFactors=FALSE))
vl0326_18_AF_RA <- as.data.frame(read.table("18_AF_RA_peaks.narrowPeak",header = FALSE, sep="\t",stringsAsFactors=FALSE))
vl0033_28_AF_RA <- as.data.frame(read.table("28_AF_RA_peaks.narrowPeak",header = FALSE, sep="\t",stringsAsFactors=FALSE))
vl0028_25_AF_LA <- as.data.frame(read.table("25_AF_LA_peaks.narrowPeak",header = FALSE, sep="\t",stringsAsFactors=FALSE))
vl0032_28_AF_LA <- as.data.frame(read.table("28_AF_LA_peaks.narrowPeak",header = FALSE, sep="\t",stringsAsFactors=FALSE))
vl0212_47_AF_LA <- as.data.frame(read.table("47_AF_LA_peaks.narrowPeak",header = FALSE, sep="\t",stringsAsFactors=FALSE))
vl0215_31_AF_LA <- as.data.frame(read.table("31_AF_LA_peaks.narrowPeak",header = FALSE, sep="\t",stringsAsFactors=FALSE))
vl0256_2231_AF_LA <- as.data.frame(read.table("2231_AF_LA_peaks.narrowPeak",header = FALSE, sep="\t",stringsAsFactors=FALSE))
vl0257_313_AF_LA <- as.data.frame(read.table("313_AF_LA_peaks.narrowPeak",header = FALSE, sep="\t",stringsAsFactors=FALSE))
vl0327_18_AF_LA <- as.data.frame(read.table("18_AF_LA_peaks.narrowPeak",header = FALSE, sep="\t",stringsAsFactors=FALSE))
vl0027_12b_SR_RA <- as.data.frame(read.table("12b_SR_RA_peaks.narrowPeak",header = FALSE, sep="\t",stringsAsFactors=FALSE))
vl0031_27_SR_RA <- as.data.frame(read.table("27_SR_RA_peaks.narrowPeak",header = FALSE, sep="\t",stringsAsFactors=FALSE))
vl0218_33_SR_RA <- as.data.frame(read.table("33_SR_RA_peaks.narrowPeak",header = FALSE, sep="\t",stringsAsFactors=FALSE))
vl0261_20_SR_RA <- as.data.frame(read.table("20_SR_RA_peaks.narrowPeak",header = FALSE, sep="\t",stringsAsFactors=FALSE))
vl0026_12b_SR_LA <- as.data.frame(read.table("12b_SR_LA_peaks.narrowPeak",header = FALSE, sep="\t",stringsAsFactors=FALSE))
vl0213_45_SR_LA <- as.data.frame(read.table("45_SR_LA_peaks.narrowPeak",header = FALSE, sep="\t",stringsAsFactors=FALSE))
vl0217_20_SR_LA <- as.data.frame(read.table("20_SR_LA_peaks.narrowPeak",header = FALSE, sep="\t",stringsAsFactors=FALSE))
vl0328_27_SR_LA <- as.data.frame(read.table("27_SR_LA_peaks.narrowPeak",header = FALSE, sep="\t",stringsAsFactors=FALSE))

vl0029_25_AF_RA_length <- vl0029_25_AF_RA$V3 - vl0029_25_AF_RA$V2
vl0214_2231_AF_RA_length <- vl0214_2231_AF_RA$V3 - vl0214_2231_AF_RA$V2
vl0258_31_AF_RA_length <- vl0258_31_AF_RA$V3 - vl0258_31_AF_RA$V2
vl0260_313_AF_RA_length <- vl0260_313_AF_RA$V3 - vl0260_313_AF_RA$V2
vl0326_18_AF_RA_length <- vl0326_18_AF_RA$V3 - vl0326_18_AF_RA$V2
vl0033_28_AF_RA_length <- vl0033_28_AF_RA$V3 - vl0033_28_AF_RA$V2
vl0028_25_AF_LA_length <- vl0028_25_AF_LA$V3 - vl0028_25_AF_LA$V2
vl0032_28_AF_LA_length <- vl0032_28_AF_LA$V3 - vl0032_28_AF_LA$V2
vl0212_47_AF_LA_length <- vl0212_47_AF_LA$V3 - vl0212_47_AF_LA$V2
vl0215_31_AF_LA_length <- vl0215_31_AF_LA$V3 - vl0215_31_AF_LA$V2
vl0256_2231_AF_LA_length <- vl0256_2231_AF_LA$V3 - vl0256_2231_AF_LA$V2
vl0257_313_AF_LA_length <- vl0257_313_AF_LA$V3 - vl0257_313_AF_LA$V2
vl0327_18_AF_LA_length <- vl0327_18_AF_LA$V3 - vl0327_18_AF_LA$V2
vl0027_12b_SR_RA_length <- vl0027_12b_SR_RA$V3 - vl0027_12b_SR_RA$V2
vl0031_27_SR_RA_length <- vl0031_27_SR_RA$V3 - vl0031_27_SR_RA$V2
vl0218_33_SR_RA_length <- vl0218_33_SR_RA$V3 - vl0218_33_SR_RA$V2
vl0261_20_SR_RA_length <- vl0261_20_SR_RA$V3 - vl0261_20_SR_RA$V2
vl0026_12b_SR_LA_length <- vl0026_12b_SR_LA$V3 - vl0026_12b_SR_LA$V2
vl0213_45_SR_LA_length <- vl0213_45_SR_LA$V3 - vl0213_45_SR_LA$V2
vl0217_20_SR_LA_length <- vl0217_20_SR_LA$V3 - vl0217_20_SR_LA$V2
vl0328_27_SR_LA_length <- vl0328_27_SR_LA$V3 - vl0328_27_SR_LA$V2

vl0029_25_AF_RA_length <- melt(vl0029_25_AF_RA_length)
vl0214_2231_AF_RA_length <- melt(vl0214_2231_AF_RA_length)
vl0258_31_AF_RA_length <- melt(vl0258_31_AF_RA_length)
vl0260_313_AF_RA_length <- melt(vl0260_313_AF_RA_length)
vl0326_18_AF_RA_length <- melt(vl0326_18_AF_RA_length)
vl0033_28_AF_RA_length <- melt(vl0033_28_AF_RA_length)
vl0028_25_AF_LA_length <- melt(vl0028_25_AF_LA_length)
vl0032_28_AF_LA_length <- melt(vl0032_28_AF_LA_length)
vl0212_47_AF_LA_length <- melt(vl0212_47_AF_LA_length)
vl0215_31_AF_LA_length <- melt(vl0215_31_AF_LA_length)
vl0256_2231_AF_LA_length <- melt(vl0256_2231_AF_LA_length)
vl0257_313_AF_LA_length <- melt(vl0257_313_AF_LA_length)
vl0327_18_AF_LA_length <- melt(vl0327_18_AF_LA_length)
vl0027_12b_SR_RA_length <- melt(vl0027_12b_SR_RA_length)
vl0031_27_SR_RA_length<- melt(vl0031_27_SR_RA_length)
vl0218_33_SR_RA_length<- melt(vl0218_33_SR_RA_length)
vl0261_20_SR_RA_length<- melt(vl0261_20_SR_RA_length)
vl0026_12b_SR_LA_length<- melt(vl0026_12b_SR_LA_length)
vl0213_45_SR_LA_length<- melt(vl0213_45_SR_LA_length)
vl0217_20_SR_LA_length<- melt(vl0217_20_SR_LA_length)
vl0328_27_SR_LA_length<- melt(vl0328_27_SR_LA_length)

vl0029_25_AF_RA_length$sample <- "vl0029_25_AF_RA"
vl0214_2231_AF_RA_length$sample <- "vl0214_2231_AF_RA"
vl0258_31_AF_RA_length$sample <- "vl0258_31_AF_RA"
vl0260_313_AF_RA_length$sample <- "vl0260_313_AF_RA"
vl0326_18_AF_RA_length$sample <- "vl0326_18_AF_RA"
vl0033_28_AF_RA_length$sample <- "vl0033_28_AF_RA"
vl0028_25_AF_LA_length$sample <- "vl0028_25_AF_LA"
vl0032_28_AF_LA_length$sample <- "vl0032_28_AF_LA"
vl0212_47_AF_LA_length$sample <- "vl0212_47_AF_LA"
vl0215_31_AF_LA_length$sample <- "vl0215_31_AF_LA"
vl0256_2231_AF_LA_length$sample <- "vl0256_2231_AF_LA"
vl0257_313_AF_LA_length$sample <- "vl0257_313_AF_LA"
vl0327_18_AF_LA_length$sample <- "vl0327_18_AF_LA"
vl0027_12b_SR_RA_length$sample <- "vl0027_12b_SR_RA"
vl0031_27_SR_RA_length$sample <- "vl0031_27_SR_RA"
vl0218_33_SR_RA_length$sample<- "vl0218_33_SR_RA"
vl0261_20_SR_RA_length$sample<- "vl0261_20_SR_RA"
vl0026_12b_SR_LA_length$sample <- "vl0026_12b_SR_LA"
vl0213_45_SR_LA_length$sample <- "vl0213_45_SR_LA"
vl0217_20_SR_LA_length$sample <- "vl0217_20_SR_LA"
vl0328_27_SR_LA_length$sample <- "vl0328_27_SR_LA"

vl0029_25_AF_RA_length$sampleid <- "25-AF-RA"
vl0214_2231_AF_RA_length$sampleid <- "2231-AF-RA"
vl0258_31_AF_RA_length$sampleid <- "31-AF-RA"
vl0260_313_AF_RA_length$sampleid <- "313-AF-RA"
vl0326_18_AF_RA_length$sampleid <- "18-AF-RA"
vl0033_28_AF_RA_length$sampleid <- "28-AF-RA"
vl0028_25_AF_LA_length$sampleid <- "25-AF-LA"
vl0032_28_AF_LA_length$sampleid <- "28-AF-LA"
vl0212_47_AF_LA_length$sampleid <- "47-AF-LA"
vl0215_31_AF_LA_length$sampleid <- "31-AF-LA"
vl0256_2231_AF_LA_length$sampleid <- "2231-AF-LA"
vl0257_313_AF_LA_length$sampleid <- "313-AF-LA"
vl0327_18_AF_LA_length$sampleid <- "18-AF-LA"
vl0027_12b_SR_RA_length$sampleid <- "12b-SR-RA"
vl0031_27_SR_RA_length$sampleid <- "27-SR-RA"
vl0218_33_SR_RA_length$sampleid<- "33-SR-RA"
vl0261_20_SR_RA_length$sampleid<- "20-SR-RA"
vl0026_12b_SR_LA_length$sampleid <- "12b-SR-LA"
vl0213_45_SR_LA_length$sampleid <- "45-SR-LA"
vl0217_20_SR_LA_length$sampleid <- "20-SR-LA"
vl0328_27_SR_LA_length$sampleid <- "27-SR-LA"

mean(vl0029_25_AF_RA_length$value)
mean(vl0214_2231_AF_RA_length$value)
mean(vl0258_31_AF_RA_length$value)
mean(vl0260_313_AF_RA_length$value)
mean(vl0326_18_AF_RA_length$value)
mean(vl0033_28_AF_RA_length$value)
mean(vl0028_25_AF_LA_length$value)
mean(vl0032_28_AF_LA_length$value)
mean(vl0212_47_AF_LA_length$value)
mean(vl0215_31_AF_LA_length$value)
mean(vl0256_2231_AF_LA_length$value)
mean(vl0257_313_AF_LA_length$value)
mean(vl0327_18_AF_LA_length$value)
mean(vl0027_12b_SR_RA_length$value)
mean(vl0031_27_SR_RA_length$value)
mean(vl0218_33_SR_RA_length$value)
mean(vl0261_20_SR_RA_length$value)
mean(vl0026_12b_SR_LA_length$value)
mean(vl0213_45_SR_LA_length$value)
mean(vl0217_20_SR_LA_length$value)
mean(vl0328_27_SR_LA_length$value)

vl0029_25_AF_RA_length$class <- "AF-RA"
vl0214_2231_AF_RA_length$class <- "AF-RA"
vl0258_31_AF_RA_length$class <- "AF-RA"
vl0260_313_AF_RA_length$class <- "AF-RA"
vl0326_18_AF_RA_length$class <- "AF-RA"
vl0033_28_AF_RA_length$class <- "AF-RA"
vl0028_25_AF_LA_length$class <- "AF-LA"
vl0032_28_AF_LA_length$class <- "AF-LA"
vl0212_47_AF_LA_length$class <- "AF-LA"
vl0215_31_AF_LA_length$class <- "AF-LA"
vl0256_2231_AF_LA_length$class <- "AF-LA"
vl0257_313_AF_LA_length$class <- "AF-LA"
vl0327_18_AF_LA_length$class <- "AF-LA"
vl0027_12b_SR_RA_length$class <- "SR-RA"
vl0031_27_SR_RA_length$class <- "SR-RA"
vl0218_33_SR_RA_length$class<- "SR-RA"
vl0261_20_SR_RA_length$class<- "SR-RA"
vl0026_12b_SR_LA_length$class <- "SR-LA"
vl0213_45_SR_LA_length$class <- "SR-LA"
vl0217_20_SR_LA_length$class <- "SR-LA"
vl0328_27_SR_LA_length$class <- "SR-LA"

peaksets_lengths <- rbind(vl0029_25_AF_RA_length,vl0214_2231_AF_RA_length,vl0260_313_AF_RA_length,vl0258_31_AF_RA_length,
                          vl0326_18_AF_RA_length,vl0033_28_AF_RA_length,vl0028_25_AF_LA_length,vl0032_28_AF_LA_length,
                          vl0212_47_AF_LA_length,vl0215_31_AF_LA_length,vl0256_2231_AF_LA_length,vl0257_313_AF_LA_length,
                          vl0027_12b_SR_RA_length,vl0031_27_SR_RA_length,vl0218_33_SR_RA_length,
                          vl0261_20_SR_RA_length,vl0026_12b_SR_LA_length,vl0213_45_SR_LA_length,vl0217_20_SR_LA_length,
                          vl0328_27_SR_LA_length)

mean(peaksets_lengths$value)

level_order <- c('vl0029_25_AF_RA','vl0214_2231_AF_RA','vl0260_313_AF_RA','vl0258_31_AF_RA','vl0326_18_AF_RA','vl0033_28_AF_RA','vl0028_25_AF_LA','vl0032_28_AF_LA',
                 'vl0212_47_AF_LA','vl0215_31_AF_LA','vl0256_2231_AF_LA','vl0257_313_AF_LA','vl0027_12b_SR_RA','vl0031_27_SR_RA','vl0218_33_SR_RA','vl0261_20_SR_RA',
                 'vl0026_12b_SR_LA','vl0213_45_SR_LA','vl0217_20_SR_LA','vl0328_27_SR_LA')

## Explore average peakset lengths in each group:
df_avg_length <- peaksets_lengths %>%
  group_by(class) %>%
  summarise(avg_peak_length = mean(value, na.rm = TRUE))
print(df_avg_length)

level_order <- c("25-AF-LA","28-AF-LA","47-AF-LA","31-AF-LA","2231-AF-LA","313-AF-LA","25-AF-RA","2231-AF-RA","31-AF-RA","313-AF-RA","18-AF-RA","28-AF-RA","12b-SR-LA","45-SR-LA","20-SR-LA","27-SR-LA","12b-SR-RA","27-SR-RA","33-SR-RA","20-SR-RA")
sample_peakset_lengths_boxplot_nooutliers <- ggplot(peaksets_lengths, aes(x=factor(sampleid, level = level_order), y=value)) + geom_boxplot(aes(fill=class),outlier.shape=NA) + coord_cartesian(ylim = c(0, 3500)) + scale_y_continuous("Peak length (bp)") + scale_x_discrete("Atrial sample") + theme_classic() + theme(axis.text.x = element_text(angle=45, hjust = 1) )  + scale_fill_manual(values=c('#D85B74', '#E3D081', '#7298BC','#91C7B1'))                                                                                                                                                                                                                                                                                                                                                                                                                                                                       
sample_peakset_lengths_boxplot_nooutliers # 📊 Figure S1c


#Explore correlation across samples and distribution in multidimensional space


setwd("/Rodriguez_2025_AFepigenome/H3K27ac_enriched_regions_definition/")


load("DbaCountMatrix_AF_genrich_updated.RData") #Upload DBA object using count matrix with Genrich peaks as reference (improves visualisation)
matrix

#PCA plot (subgroups)
dba.plotPCA(matrix,DBA_TISSUE:DBA_CONDITION,vColors = c("#7298BC","#91C7B1","#D85B74","#E3D081"),label = DBA_ID)
dba.plotPCA(matrix,DBA_TISSUE:DBA_CONDITION,vColors = c("#7298BC","#91C7B1","#D85B74","#E3D081")) # 📊 Figure 1a


#Correlation heatmap

load("/Rodriguez_2025_AFepigenome/ChIP-seq_basic_processing/AF_final_dataset_H3K27ac_Dec23_dba.RData") #Upload DBA object as generated by DiffBind (consensus peakset across samples used as reference) - as generated in ChIP-seq_basic_processing/ChIPseq_DiffEnrich_Analysis_DiffBind_TCseq.R script 
AF_final_dataset_H3K27ac_Dec23_dba <- dba(AF_final_dataset_H3K27ac_Dec23_dba)
### H3K27ac correlation plot --> heatmap based on count scores (corrected by read density)
debug(dba.plotHeatmap)
dba.plotHeatmap(AF_final_dataset_H3K27ac_Dec23_dba) #save "res" object as "correlation_values_AF_H3K27ac_dataset.csv" containing correlation values across samples to local environment"
#write.csv(res,file="corr_values_AF_H3K27ac_dataset_2.csv")
undebug(dba.plotHeatmap)

corr_values_AF_H3K27ac_dataset <- read.csv("corr_values_AF_H3K27ac_dataset_2.csv")
library(reshape2)
corr_values_AF_H3K27ac_dataset_melted <- melt(corr_values_AF_H3K27ac_dataset)
head(corr_values_AF_H3K27ac_dataset_melted)
peaksets_orderforheatmaps <- c("vl0028","vl0032","vl0212","vl0215","vl0256","vl0257","vl0026","vl0213","vl0217","vl0328","vl0029","vl0033","vl0214","vl0258","vl0260","vl0326","vl0027","vl0031","vl0218","vl0261")
library(ggplot2)

peaksets_orderforheatmaps <- c("vl0212","vl0215","vl0028","vl0032","vl0256","vl0257","vl0213","vl0217","vl0026","vl0328", "vl0326","vl0029","vl0033","vl0214","vl0258","vl0260","vl0027","vl0031","vl0261","vl0218")
corr_values_heatmap <- ggplot(data = corr_values_AF_H3K27ac_dataset_melted, aes(x=factor(X,level=peaksets_orderforheatmaps), y=factor(variable, level=peaksets_orderforheatmaps), fill=value)) +
  geom_tile() + scale_fill_gradientn(colours=c("#FFFFCC", "#C7E9B4", "#2F8EBF", "#225EA8", "#0C2C84")) + theme(axis.text=element_text(size=10),axis.text.x=element_text(angle = -90, hjust = 0))
corr_values_heatmap # 📊 Figure 1b





###############################################################################################
##### Definition of enriched regions across sample groups based on pairwise comparisons #######
###############################################################################################



# Explore pairwise comparisons peaksets 

## AFLA

setwd("/Rodriguez_2025_AFepigenome/H3K27ac_enriched_regions_definition/AF_LA_enriched_peaks")
myFiles <- list.files(pattern="*.bed")

#AFLA
diff_12b_SR_RA_vs_2231_AF_LA_c2.0_cond2.bed <- as.data.frame(read.table("diff_12b_SR_RA_vs_2231_AF_LA_c2.0_cond2.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_12b_SR_RA_vs_25_AF_LA_c2.0_cond2.bed <- as.data.frame(read.table("diff_12b_SR_RA_vs_25_AF_LA_c2.0_cond2.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_12b_SR_RA_vs_28_AF_LA_c2.0_cond2.bed <- as.data.frame(read.table("diff_12b_SR_RA_vs_28_AF_LA_c2.0_cond2.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_12b_SR_RA_vs_31_AF_LA_c2.0_cond2.bed <- as.data.frame(read.table("diff_12b_SR_RA_vs_31_AF_LA_c2.0_cond2.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_12b_SR_RA_vs_313_AF_LA_c2.0_cond2.bed <- as.data.frame(read.table("diff_12b_SR_RA_vs_313_AF_LA_c2.0_cond2.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_12b_SR_RA_vs_47_AF_LA_c2.0_cond2.bed <- as.data.frame(read.table("diff_12b_SR_RA_vs_47_AF_LA_c2.0_cond2.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_18_AF_RA_vs_2231_AF_LA_c2.0_cond2.bed <- as.data.frame(read.table("diff_18_AF_RA_vs_2231_AF_LA_c2.0_cond2.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_18_AF_RA_vs_25_AF_LA_c2.0_cond2.bed <- as.data.frame(read.table("diff_18_AF_RA_vs_25_AF_LA_c2.0_cond2.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_18_AF_RA_vs_28_AF_LA_c2.0_cond2.bed <- as.data.frame(read.table("diff_18_AF_RA_vs_28_AF_LA_c2.0_cond2.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_18_AF_RA_vs_31_AF_LA_c2.0_cond2.bed <- as.data.frame(read.table("diff_18_AF_RA_vs_31_AF_LA_c2.0_cond2.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_18_AF_RA_vs_313_AF_LA_c2.0_cond2.bed <- as.data.frame(read.table("diff_18_AF_RA_vs_313_AF_LA_c2.0_cond2.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_18_AF_RA_vs_47_AF_LA_c2.0_cond2.bed <- as.data.frame(read.table("diff_18_AF_RA_vs_47_AF_LA_c2.0_cond2.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_20_SR_RA_vs_2231_AF_LA_c2.0_cond2.bed <- as.data.frame(read.table("diff_20_SR_RA_vs_2231_AF_LA_c2.0_cond2.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_20_SR_RA_vs_25_AF_LA_c2.0_cond2.bed <- as.data.frame(read.table("diff_20_SR_RA_vs_25_AF_LA_c2.0_cond2.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_20_SR_RA_vs_28_AF_LA_c2.0_cond2.bed <- as.data.frame(read.table("diff_20_SR_RA_vs_28_AF_LA_c2.0_cond2.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_20_SR_RA_vs_31_AF_LA_c2.0_cond2.bed <- as.data.frame(read.table("diff_20_SR_RA_vs_31_AF_LA_c2.0_cond2.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_20_SR_RA_vs_313_AF_LA_c2.0_cond2.bed <- as.data.frame(read.table("diff_20_SR_RA_vs_313_AF_LA_c2.0_cond2.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_20_SR_RA_vs_47_AF_LA_c2.0_cond2.bed <- as.data.frame(read.table("diff_20_SR_RA_vs_47_AF_LA_c2.0_cond2.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_2231_AF_LA_vs_20_SR_LA_c2.0_cond1.bed <- as.data.frame(read.table("diff_2231_AF_LA_vs_20_SR_LA_c2.0_cond1.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_2231_AF_LA_vs_27_SR_LA_c2.0_cond1.bed <- as.data.frame(read.table("diff_2231_AF_LA_vs_27_SR_LA_c2.0_cond1.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_2231_AF_LA_vs_45_SR_LA_c2.0_cond1.bed <- as.data.frame(read.table("diff_2231_AF_LA_vs_45_SR_LA_c2.0_cond1.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_2231_AF_RA_vs_2231_AF_LA_c2.0_cond2.bed <- as.data.frame(read.table("diff_2231_AF_RA_vs_2231_AF_LA_c2.0_cond2.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_2231_AF_RA_vs_25_AF_LA_c2.0_cond2.bed <- as.data.frame(read.table("diff_2231_AF_RA_vs_25_AF_LA_c2.0_cond2.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_2231_AF_RA_vs_28_AF_LA_c2.0_cond2.bed <- as.data.frame(read.table("diff_2231_AF_RA_vs_28_AF_LA_c2.0_cond2.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_2231_AF_RA_vs_31_AF_LA_c2.0_cond2.bed <- as.data.frame(read.table("diff_2231_AF_RA_vs_31_AF_LA_c2.0_cond2.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_2231_AF_RA_vs_313_AF_LA_c2.0_cond2.bed <- as.data.frame(read.table("diff_2231_AF_RA_vs_313_AF_LA_c2.0_cond2.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_2231_AF_RA_vs_47_AF_LA_c2.0_cond2.bed <- as.data.frame(read.table("diff_2231_AF_RA_vs_47_AF_LA_c2.0_cond2.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_2231_AF_LA_vs_12b_SR_LA_c2.0_cond1.bed <- as.data.frame(read.table("diff_2231_AF_LA_vs_12b_SR_LA_c2.0_cond1.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_25_AF_LA_vs_12b_SR_LA_c2.0_cond1.bed <- as.data.frame(read.table("diff_25_AF_LA_vs_12b_SR_LA_c2.0_cond1.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_25_AF_LA_vs_27_SR_LA_c2.0_cond1.bed <- as.data.frame(read.table("diff_25_AF_LA_vs_27_SR_LA_c2.0_cond1.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_25_AF_LA_vs_20_SR_LA_c2.0_cond1.bed <- as.data.frame(read.table("diff_25_AF_LA_vs_20_SR_LA_c2.0_cond1.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_25_AF_LA_vs_45_SR_LA_c2.0_cond1.bed <- as.data.frame(read.table("diff_25_AF_LA_vs_45_SR_LA_c2.0_cond1.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_25_AF_RA_vs_2231_AF_LA_c2.0_cond2.bed <- as.data.frame(read.table("diff_25_AF_RA_vs_2231_AF_LA_c2.0_cond2.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_25_AF_RA_vs_25_AF_LA_c2.0_cond2.bed <- as.data.frame(read.table("diff_25_AF_RA_vs_25_AF_LA_c2.0_cond2.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_25_AF_RA_vs_28_AF_LA_c2.0_cond2.bed <- as.data.frame(read.table("diff_25_AF_RA_vs_28_AF_LA_c2.0_cond2.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_25_AF_RA_vs_31_AF_LA_c2.0_cond2.bed <- as.data.frame(read.table("diff_12b_SR_RA_vs_25_AF_LA_c2.0_cond2.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_25_AF_RA_vs_313_AF_LA_c2.0_cond2.bed <- as.data.frame(read.table("diff_25_AF_RA_vs_313_AF_LA_c2.0_cond2.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_25_AF_RA_vs_47_AF_LA_c2.0_cond2.bed <- as.data.frame(read.table("diff_25_AF_RA_vs_47_AF_LA_c2.0_cond2.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_27_SR_RA_vs_2231_AF_LA_c2.0_cond2.bed <- as.data.frame(read.table("diff_27_SR_RA_vs_2231_AF_LA_c2.0_cond2.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_27_SR_RA_vs_25_AF_LA_c2.0_cond2.bed <- as.data.frame(read.table("diff_27_SR_RA_vs_25_AF_LA_c2.0_cond2.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_27_SR_RA_vs_28_AF_LA_c2.0_cond2.bed <- as.data.frame(read.table("diff_27_SR_RA_vs_28_AF_LA_c2.0_cond2.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_27_SR_RA_vs_31_AF_LA_c2.0_cond2.bed <- as.data.frame(read.table("diff_27_SR_RA_vs_31_AF_LA_c2.0_cond2.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_27_SR_RA_vs_313_AF_LA_c2.0_cond2.bed <- as.data.frame(read.table("diff_27_SR_RA_vs_313_AF_LA_c2.0_cond2.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_27_SR_RA_vs_47_AF_LA_c2.0_cond2.bed <- as.data.frame(read.table("diff_27_SR_RA_vs_47_AF_LA_c2.0_cond2.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_28_AF_LA_vs_12b_SR_LA_c2.0_cond1.bed <- as.data.frame(read.table("diff_28_AF_LA_vs_12b_SR_LA_c2.0_cond1.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_28_AF_LA_vs_20_SR_LA_c2.0_cond1.bed <- as.data.frame(read.table("diff_28_AF_LA_vs_20_SR_LA_c2.0_cond1.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_28_AF_LA_vs_27_SR_LA_c2.0_cond1.bed <- as.data.frame(read.table("diff_28_AF_LA_vs_27_SR_LA_c2.0_cond1.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_28_AF_LA_vs_45_SR_LA_c2.0_cond1.bed <- as.data.frame(read.table("diff_28_AF_LA_vs_45_SR_LA_c2.0_cond1.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_28_AF_RA_vs_2231_AF_LA_c2.0_cond2.bed <- as.data.frame(read.table("diff_28_AF_RA_vs_2231_AF_LA_c2.0_cond2.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_28_AF_RA_vs_25_AF_LA_c2.0_cond2.bed <- as.data.frame(read.table("diff_28_AF_RA_vs_25_AF_LA_c2.0_cond2.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_28_AF_RA_vs_28_AF_LA_c2.0_cond2.bed <- as.data.frame(read.table("diff_28_AF_RA_vs_28_AF_LA_c2.0_cond2.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_28_AF_RA_vs_31_AF_LA_c2.0_cond2.bed <- as.data.frame(read.table("diff_28_AF_RA_vs_31_AF_LA_c2.0_cond2.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_28_AF_RA_vs_313_AF_LA_c2.0_cond2.bed <- as.data.frame(read.table("diff_28_AF_RA_vs_313_AF_LA_c2.0_cond2.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_31_AF_LA_vs_12b_SR_LA_c2.0_cond1.bed <- as.data.frame(read.table("diff_31_AF_LA_vs_12b_SR_LA_c2.0_cond1.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_31_AF_LA_vs_20_SR_LA_c2.0_cond1.bed <- as.data.frame(read.table("diff_31_AF_LA_vs_20_SR_LA_c2.0_cond1.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_31_AF_LA_vs_27_SR_LA_c2.0_cond1.bed <- as.data.frame(read.table("diff_31_AF_LA_vs_27_SR_LA_c2.0_cond1.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_31_AF_LA_vs_45_SR_LA_c2.0_cond1.bed <- as.data.frame(read.table("diff_31_AF_LA_vs_45_SR_LA_c2.0_cond1.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_31_AF_RA_vs_2231_AF_LA_c2.0_cond2.bed <- as.data.frame(read.table("diff_31_AF_RA_vs_2231_AF_LA_c2.0_cond2.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_31_AF_RA_vs_28_AF_LA_c2.0_cond2.bed <- as.data.frame(read.table("diff_31_AF_RA_vs_28_AF_LA_c2.0_cond2.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_31_AF_RA_vs_31_AF_LA_c2.0_cond2.bed <- as.data.frame(read.table("diff_31_AF_RA_vs_31_AF_LA_c2.0_cond2.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_31_AF_RA_vs_313_AF_LA_c2.0_cond2.bed <- as.data.frame(read.table("diff_31_AF_RA_vs_313_AF_LA_c2.0_cond2.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_31_AF_RA_vs_47_AF_LA_c2.0_cond2.bed <- as.data.frame(read.table("diff_31_AF_RA_vs_47_AF_LA_c2.0_cond2.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_313_AF_LA_vs_12b_SR_LA_c2.0_cond1.bed <- as.data.frame(read.table("diff_313_AF_LA_vs_12b_SR_LA_c2.0_cond1.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_313_AF_LA_vs_20_SR_LA_c2.0_cond1.bed <- as.data.frame(read.table("diff_313_AF_LA_vs_20_SR_LA_c2.0_cond1.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_313_AF_LA_vs_27_SR_LA_c2.0_cond1.bed <- as.data.frame(read.table("diff_313_AF_LA_vs_27_SR_LA_c2.0_cond1.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_313_AF_LA_vs_45_SR_LA_c2.0_cond1.bed <- as.data.frame(read.table("diff_313_AF_LA_vs_45_SR_LA_c2.0_cond1.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_313_AF_RA_vs_2231_AF_LA_c2.0_cond2.bed <- as.data.frame(read.table("diff_313_AF_RA_vs_2231_AF_LA_c2.0_cond2.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_313_AF_RA_vs_25_AF_LA_c2.0_cond2.bed <- as.data.frame(read.table("diff_313_AF_RA_vs_25_AF_LA_c2.0_cond2.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_313_AF_RA_vs_28_AF_LA_c2.0_cond2.bed <- as.data.frame(read.table("diff_313_AF_RA_vs_28_AF_LA_c2.0_cond2.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_313_AF_RA_vs_31_AF_LA_c2.0_cond2.bed <- as.data.frame(read.table("diff_313_AF_RA_vs_31_AF_LA_c2.0_cond2.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_313_AF_RA_vs_313_AF_LA_c2.0_cond2.bed <- as.data.frame(read.table("diff_313_AF_RA_vs_313_AF_LA_c2.0_cond2.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_313_AF_RA_vs_47_AF_LA_c2.0_cond2.bed <- as.data.frame(read.table("diff_313_AF_RA_vs_47_AF_LA_c2.0_cond2.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_33_SR_RA_vs_2231_AF_LA_c2.0_cond2.bed <- as.data.frame(read.table("diff_33_SR_RA_vs_2231_AF_LA_c2.0_cond2.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_33_SR_RA_vs_25_AF_LA_c2.0_cond2.bed <- as.data.frame(read.table("diff_33_SR_RA_vs_25_AF_LA_c2.0_cond2.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_33_SR_RA_vs_28_AF_LA_c2.0_cond2.bed <- as.data.frame(read.table("diff_33_SR_RA_vs_28_AF_LA_c2.0_cond2.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_33_SR_RA_vs_31_AF_LA_c2.0_cond2.bed <- as.data.frame(read.table("diff_33_SR_RA_vs_31_AF_LA_c2.0_cond2.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_33_SR_RA_vs_313_AF_LA_c2.0_cond2.bed <- as.data.frame(read.table("diff_33_SR_RA_vs_313_AF_LA_c2.0_cond2.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_33_SR_RA_vs_47_AF_LA_c2.0_cond2.bed <- as.data.frame(read.table("diff_33_SR_RA_vs_47_AF_LA_c2.0_cond2.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_47_AF_LA_vs_12b_SR_LA_c2.0_cond1.bed <- as.data.frame(read.table("diff_47_AF_LA_vs_12b_SR_LA_c2.0_cond1.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_47_AF_LA_vs_27_SR_LA_c2.0_cond1.bed <- as.data.frame(read.table("diff_47_AF_LA_vs_27_SR_LA_c2.0_cond1.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_47_AF_LA_vs_45_SR_LA_c2.0_cond1.bed <- as.data.frame(read.table("diff_47_AF_LA_vs_45_SR_LA_c2.0_cond1.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_47_AF_LA_vs_20_SR_LA_c2.0_cond1.bed <- as.data.frame(read.table("diff_47_AF_LA_vs_20_SR_LA_c2.0_cond1.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_28_AF_RA_vs_47_AF_LA_c2.0_cond2.bed <- as.data.frame(read.table("diff_28_AF_RA_vs_47_AF_LA_c2.0_cond2.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_31_AF_RA_vs_25_AF_LA_c2.0_cond2.bed <- as.data.frame(read.table("diff_31_AF_RA_vs_25_AF_LA_c2.0_cond2.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))

diff_12b_SR_RA_vs_2231_AF_LA_c2.0_cond2.bed_length <- diff_12b_SR_RA_vs_2231_AF_LA_c2.0_cond2.bed$V3 - diff_12b_SR_RA_vs_2231_AF_LA_c2.0_cond2.bed$V2
diff_12b_SR_RA_vs_25_AF_LA_c2.0_cond2.bed_length <- diff_12b_SR_RA_vs_25_AF_LA_c2.0_cond2.bed$V3 - diff_12b_SR_RA_vs_25_AF_LA_c2.0_cond2.bed$V2
diff_12b_SR_RA_vs_28_AF_LA_c2.0_cond2.bed_length <- diff_12b_SR_RA_vs_28_AF_LA_c2.0_cond2.bed$V3 - diff_12b_SR_RA_vs_28_AF_LA_c2.0_cond2.bed$V2
diff_12b_SR_RA_vs_31_AF_LA_c2.0_cond2.bed_length<- diff_12b_SR_RA_vs_31_AF_LA_c2.0_cond2.bed$V3 - diff_12b_SR_RA_vs_31_AF_LA_c2.0_cond2.bed$V2
diff_12b_SR_RA_vs_313_AF_LA_c2.0_cond2.bed_length<- diff_12b_SR_RA_vs_313_AF_LA_c2.0_cond2.bed$V3 - diff_12b_SR_RA_vs_313_AF_LA_c2.0_cond2.bed$V2
diff_12b_SR_RA_vs_47_AF_LA_c2.0_cond2.bed_length<- diff_12b_SR_RA_vs_47_AF_LA_c2.0_cond2.bed$V3 - diff_12b_SR_RA_vs_47_AF_LA_c2.0_cond2.bed$V2
diff_18_AF_RA_vs_2231_AF_LA_c2.0_cond2.bed_length<- diff_18_AF_RA_vs_2231_AF_LA_c2.0_cond2.bed$V3 - diff_18_AF_RA_vs_2231_AF_LA_c2.0_cond2.bed$V2
diff_18_AF_RA_vs_25_AF_LA_c2.0_cond2.bed_length<- diff_18_AF_RA_vs_25_AF_LA_c2.0_cond2.bed$V3 - diff_18_AF_RA_vs_25_AF_LA_c2.0_cond2.bed$V2
diff_18_AF_RA_vs_28_AF_LA_c2.0_cond2.bed_length<- diff_18_AF_RA_vs_28_AF_LA_c2.0_cond2.bed$V3 - diff_18_AF_RA_vs_28_AF_LA_c2.0_cond2.bed$V2
diff_18_AF_RA_vs_31_AF_LA_c2.0_cond2.bed_length<- diff_18_AF_RA_vs_31_AF_LA_c2.0_cond2.bed$V3 - diff_18_AF_RA_vs_31_AF_LA_c2.0_cond2.bed$V2
diff_18_AF_RA_vs_313_AF_LA_c2.0_cond2.bed_length<- diff_18_AF_RA_vs_31_AF_LA_c2.0_cond2.bed$V3 - diff_18_AF_RA_vs_31_AF_LA_c2.0_cond2.bed$V2
diff_18_AF_RA_vs_47_AF_LA_c2.0_cond2.bed_length<- diff_18_AF_RA_vs_47_AF_LA_c2.0_cond2.bed$V3 - diff_18_AF_RA_vs_47_AF_LA_c2.0_cond2.bed$V2
diff_20_SR_RA_vs_2231_AF_LA_c2.0_cond2.bed_length<- diff_20_SR_RA_vs_2231_AF_LA_c2.0_cond2.bed$V3 - diff_20_SR_RA_vs_2231_AF_LA_c2.0_cond2.bed$V2
diff_20_SR_RA_vs_25_AF_LA_c2.0_cond2.bed_length<- diff_20_SR_RA_vs_25_AF_LA_c2.0_cond2.bed$V3 - diff_20_SR_RA_vs_25_AF_LA_c2.0_cond2.bed$V2
diff_20_SR_RA_vs_28_AF_LA_c2.0_cond2.bed_length<- diff_20_SR_RA_vs_28_AF_LA_c2.0_cond2.bed$V3 - diff_20_SR_RA_vs_28_AF_LA_c2.0_cond2.bed$V2
diff_20_SR_RA_vs_31_AF_LA_c2.0_cond2.bed_length<- diff_20_SR_RA_vs_31_AF_LA_c2.0_cond2.bed$V3 - diff_20_SR_RA_vs_31_AF_LA_c2.0_cond2.bed$V2
diff_20_SR_RA_vs_313_AF_LA_c2.0_cond2.bed_length<- diff_20_SR_RA_vs_313_AF_LA_c2.0_cond2.bed$V3 - diff_20_SR_RA_vs_313_AF_LA_c2.0_cond2.bed$V2
diff_20_SR_RA_vs_47_AF_LA_c2.0_cond2.bed_length<- diff_20_SR_RA_vs_47_AF_LA_c2.0_cond2.bed$V3 - diff_20_SR_RA_vs_47_AF_LA_c2.0_cond2.bed$V2
diff_2231_AF_LA_vs_12b_SR_LA_c2.0_cond1.bed_length<- diff_2231_AF_LA_vs_12b_SR_LA_c2.0_cond1.bed$V3 - diff_2231_AF_LA_vs_12b_SR_LA_c2.0_cond1.bed$V2
diff_2231_AF_LA_vs_20_SR_LA_c2.0_cond1.bed_length<- diff_2231_AF_LA_vs_20_SR_LA_c2.0_cond1.bed$V3 - diff_2231_AF_LA_vs_20_SR_LA_c2.0_cond1.bed$V2
diff_2231_AF_LA_vs_27_SR_LA_c2.0_cond1.bed_length<- diff_2231_AF_LA_vs_27_SR_LA_c2.0_cond1.bed$V3 - diff_2231_AF_LA_vs_27_SR_LA_c2.0_cond1.bed$V2
diff_2231_AF_LA_vs_45_SR_LA_c2.0_cond1.bed_length<- diff_2231_AF_LA_vs_45_SR_LA_c2.0_cond1.bed$V3 - diff_2231_AF_LA_vs_45_SR_LA_c2.0_cond1.bed$V2
diff_2231_AF_RA_vs_2231_AF_LA_c2.0_cond2.bed_length<- diff_2231_AF_RA_vs_2231_AF_LA_c2.0_cond2.bed$V3 - diff_2231_AF_RA_vs_2231_AF_LA_c2.0_cond2.bed$V2
diff_2231_AF_RA_vs_25_AF_LA_c2.0_cond2.bed_length<- diff_2231_AF_RA_vs_25_AF_LA_c2.0_cond2.bed$V3 - diff_2231_AF_RA_vs_25_AF_LA_c2.0_cond2.bed$V2
diff_2231_AF_RA_vs_28_AF_LA_c2.0_cond2.bed_length<- diff_2231_AF_RA_vs_28_AF_LA_c2.0_cond2.bed$V3 - diff_2231_AF_RA_vs_28_AF_LA_c2.0_cond2.bed$V2
diff_2231_AF_RA_vs_31_AF_LA_c2.0_cond2.bed_length<- diff_2231_AF_RA_vs_31_AF_LA_c2.0_cond2.bed$V3 - diff_2231_AF_RA_vs_31_AF_LA_c2.0_cond2.bed$V2
diff_2231_AF_RA_vs_313_AF_LA_c2.0_cond2.bed_length<- diff_2231_AF_RA_vs_313_AF_LA_c2.0_cond2.bed$V3 - diff_2231_AF_RA_vs_313_AF_LA_c2.0_cond2.bed$V2
diff_2231_AF_RA_vs_47_AF_LA_c2.0_cond2.bed_length<- diff_2231_AF_RA_vs_47_AF_LA_c2.0_cond2.bed$V3 - diff_2231_AF_RA_vs_47_AF_LA_c2.0_cond2.bed$V2
diff_25_AF_LA_vs_12b_SR_LA_c2.0_cond1.bed_length<- diff_25_AF_LA_vs_12b_SR_LA_c2.0_cond1.bed$V3 - diff_25_AF_LA_vs_12b_SR_LA_c2.0_cond1.bed$V2
diff_25_AF_LA_vs_20_SR_LA_c2.0_cond1.bed_length<- diff_25_AF_LA_vs_20_SR_LA_c2.0_cond1.bed$V3 - diff_25_AF_LA_vs_20_SR_LA_c2.0_cond1.bed$V2
diff_25_AF_LA_vs_27_SR_LA_c2.0_cond1.bed_length<- diff_25_AF_LA_vs_27_SR_LA_c2.0_cond1.bed$V3 - diff_25_AF_LA_vs_27_SR_LA_c2.0_cond1.bed$V2
diff_25_AF_LA_vs_45_SR_LA_c2.0_cond1.bed_length<- diff_25_AF_LA_vs_45_SR_LA_c2.0_cond1.bed$V3 - diff_25_AF_LA_vs_45_SR_LA_c2.0_cond1.bed$V2
diff_25_AF_RA_vs_2231_AF_LA_c2.0_cond2.bed_length<- diff_25_AF_RA_vs_2231_AF_LA_c2.0_cond2.bed$V3 - diff_25_AF_RA_vs_2231_AF_LA_c2.0_cond2.bed$V2
diff_25_AF_RA_vs_25_AF_LA_c2.0_cond2.bed_length<- diff_25_AF_RA_vs_25_AF_LA_c2.0_cond2.bed$V3 - diff_25_AF_RA_vs_25_AF_LA_c2.0_cond2.bed$V2
diff_25_AF_RA_vs_28_AF_LA_c2.0_cond2.bed_length<- diff_25_AF_RA_vs_28_AF_LA_c2.0_cond2.bed$V3 - diff_25_AF_RA_vs_28_AF_LA_c2.0_cond2.bed$V2
diff_25_AF_RA_vs_31_AF_LA_c2.0_cond2.bed_length<- diff_25_AF_RA_vs_31_AF_LA_c2.0_cond2.bed$V3 - diff_25_AF_RA_vs_31_AF_LA_c2.0_cond2.bed$V2
diff_25_AF_RA_vs_313_AF_LA_c2.0_cond2.bed_length<- diff_25_AF_RA_vs_313_AF_LA_c2.0_cond2.bed$V3 - diff_25_AF_RA_vs_313_AF_LA_c2.0_cond2.bed$V2
diff_25_AF_RA_vs_47_AF_LA_c2.0_cond2.bed_length<- diff_25_AF_RA_vs_47_AF_LA_c2.0_cond2.bed$V3 - diff_25_AF_RA_vs_47_AF_LA_c2.0_cond2.bed$V2
diff_27_SR_RA_vs_2231_AF_LA_c2.0_cond2.bed_length<- diff_27_SR_RA_vs_2231_AF_LA_c2.0_cond2.bed$V3 - diff_27_SR_RA_vs_2231_AF_LA_c2.0_cond2.bed$V2
diff_27_SR_RA_vs_25_AF_LA_c2.0_cond2.bed_length<- diff_27_SR_RA_vs_25_AF_LA_c2.0_cond2.bed$V3 - diff_27_SR_RA_vs_25_AF_LA_c2.0_cond2.bed$V2
diff_27_SR_RA_vs_28_AF_LA_c2.0_cond2.bed_length<- diff_27_SR_RA_vs_28_AF_LA_c2.0_cond2.bed$V3 - diff_27_SR_RA_vs_28_AF_LA_c2.0_cond2.bed$V2
diff_27_SR_RA_vs_31_AF_LA_c2.0_cond2.bed_length<- diff_27_SR_RA_vs_31_AF_LA_c2.0_cond2.bed$V3 - diff_27_SR_RA_vs_31_AF_LA_c2.0_cond2.bed$V2
diff_27_SR_RA_vs_313_AF_LA_c2.0_cond2.bed_length<- diff_27_SR_RA_vs_313_AF_LA_c2.0_cond2.bed$V3 - diff_27_SR_RA_vs_313_AF_LA_c2.0_cond2.bed$V2
diff_27_SR_RA_vs_47_AF_LA_c2.0_cond2.bed_length<- diff_27_SR_RA_vs_47_AF_LA_c2.0_cond2.bed$V3 - diff_27_SR_RA_vs_47_AF_LA_c2.0_cond2.bed$V2
diff_28_AF_LA_vs_12b_SR_LA_c2.0_cond1.bed_length<- diff_28_AF_LA_vs_12b_SR_LA_c2.0_cond1.bed$V3 - diff_28_AF_LA_vs_12b_SR_LA_c2.0_cond1.bed$V2
diff_28_AF_LA_vs_20_SR_LA_c2.0_cond1.bed_length<- diff_28_AF_LA_vs_20_SR_LA_c2.0_cond1.bed$V3 - diff_28_AF_LA_vs_20_SR_LA_c2.0_cond1.bed$V2
diff_28_AF_LA_vs_27_SR_LA_c2.0_cond1.bed_length<- diff_28_AF_LA_vs_27_SR_LA_c2.0_cond1.bed$V3 - diff_28_AF_LA_vs_27_SR_LA_c2.0_cond1.bed$V2
diff_28_AF_LA_vs_45_SR_LA_c2.0_cond1.bed_length<- diff_28_AF_LA_vs_45_SR_LA_c2.0_cond1.bed$V3 - diff_28_AF_LA_vs_45_SR_LA_c2.0_cond1.bed$V2
diff_28_AF_RA_vs_2231_AF_LA_c2.0_cond2.bed_length<- diff_28_AF_RA_vs_2231_AF_LA_c2.0_cond2.bed$V3 - diff_28_AF_RA_vs_2231_AF_LA_c2.0_cond2.bed$V2
diff_28_AF_RA_vs_25_AF_LA_c2.0_cond2.bed_length<- diff_28_AF_RA_vs_25_AF_LA_c2.0_cond2.bed$V3 - diff_28_AF_RA_vs_25_AF_LA_c2.0_cond2.bed$V2
diff_28_AF_RA_vs_28_AF_LA_c2.0_cond2.bed_length<- diff_28_AF_RA_vs_28_AF_LA_c2.0_cond2.bed$V3 - diff_28_AF_RA_vs_28_AF_LA_c2.0_cond2.bed$V2
diff_28_AF_RA_vs_31_AF_LA_c2.0_cond2.bed_length<- diff_28_AF_RA_vs_31_AF_LA_c2.0_cond2.bed$V3 - diff_28_AF_RA_vs_31_AF_LA_c2.0_cond2.bed$V2
diff_28_AF_RA_vs_313_AF_LA_c2.0_cond2.bed_length<- diff_28_AF_RA_vs_313_AF_LA_c2.0_cond2.bed$V3 - diff_28_AF_RA_vs_313_AF_LA_c2.0_cond2.bed$V2
diff_28_AF_RA_vs_47_AF_LA_c2.0_cond2.bed_length<- diff_28_AF_RA_vs_47_AF_LA_c2.0_cond2.bed$V3 - diff_28_AF_RA_vs_47_AF_LA_c2.0_cond2.bed$V2
diff_31_AF_LA_vs_12b_SR_LA_c2.0_cond1.bed_length<- diff_31_AF_LA_vs_12b_SR_LA_c2.0_cond1.bed$V3 - diff_31_AF_LA_vs_12b_SR_LA_c2.0_cond1.bed$V2
diff_31_AF_LA_vs_20_SR_LA_c2.0_cond1.bed_length<- diff_31_AF_LA_vs_20_SR_LA_c2.0_cond1.bed$V3 - diff_31_AF_LA_vs_20_SR_LA_c2.0_cond1.bed$V2
diff_31_AF_LA_vs_27_SR_LA_c2.0_cond1.bed_length<- diff_31_AF_LA_vs_27_SR_LA_c2.0_cond1.bed$V3 - diff_31_AF_LA_vs_27_SR_LA_c2.0_cond1.bed$V2
diff_31_AF_LA_vs_45_SR_LA_c2.0_cond1.bed_length<- diff_31_AF_LA_vs_45_SR_LA_c2.0_cond1.bed$V3 - diff_31_AF_LA_vs_45_SR_LA_c2.0_cond1.bed$V2
diff_31_AF_RA_vs_2231_AF_LA_c2.0_cond2.bed_length<- diff_31_AF_RA_vs_2231_AF_LA_c2.0_cond2.bed$V3 - diff_31_AF_RA_vs_2231_AF_LA_c2.0_cond2.bed$V2
diff_31_AF_RA_vs_25_AF_LA_c2.0_cond2.bed_length<- diff_31_AF_RA_vs_25_AF_LA_c2.0_cond2.bed$V3 - diff_31_AF_RA_vs_25_AF_LA_c2.0_cond2.bed$V2
diff_31_AF_RA_vs_28_AF_LA_c2.0_cond2.bed_length<- diff_31_AF_RA_vs_28_AF_LA_c2.0_cond2.bed$V3 - diff_31_AF_RA_vs_28_AF_LA_c2.0_cond2.bed$V2
diff_31_AF_RA_vs_31_AF_LA_c2.0_cond2.bed_length<- diff_31_AF_RA_vs_31_AF_LA_c2.0_cond2.bed$V3 - diff_31_AF_RA_vs_31_AF_LA_c2.0_cond2.bed$V2
diff_31_AF_RA_vs_313_AF_LA_c2.0_cond2.bed_length<- diff_31_AF_RA_vs_313_AF_LA_c2.0_cond2.bed$V3 - diff_31_AF_RA_vs_313_AF_LA_c2.0_cond2.bed$V2
diff_31_AF_RA_vs_47_AF_LA_c2.0_cond2.bed_length<- diff_31_AF_RA_vs_47_AF_LA_c2.0_cond2.bed$V3 - diff_31_AF_RA_vs_47_AF_LA_c2.0_cond2.bed$V2
diff_313_AF_LA_vs_12b_SR_LA_c2.0_cond1.bed_length<- diff_313_AF_LA_vs_12b_SR_LA_c2.0_cond1.bed$V3 - diff_313_AF_LA_vs_12b_SR_LA_c2.0_cond1.bed$V2
diff_313_AF_LA_vs_20_SR_LA_c2.0_cond1.bed_length<- diff_313_AF_LA_vs_20_SR_LA_c2.0_cond1.bed$V3 - diff_313_AF_LA_vs_20_SR_LA_c2.0_cond1.bed$V2
diff_313_AF_LA_vs_27_SR_LA_c2.0_cond1.bed_length<- diff_313_AF_LA_vs_27_SR_LA_c2.0_cond1.bed$V3 - diff_313_AF_LA_vs_27_SR_LA_c2.0_cond1.bed$V2
diff_313_AF_LA_vs_45_SR_LA_c2.0_cond1.bed_length<- diff_313_AF_LA_vs_45_SR_LA_c2.0_cond1.bed$V3 - diff_313_AF_LA_vs_45_SR_LA_c2.0_cond1.bed$V2
diff_313_AF_RA_vs_2231_AF_LA_c2.0_cond2.bed_length<- diff_313_AF_RA_vs_2231_AF_LA_c2.0_cond2.bed$V3 - diff_313_AF_RA_vs_2231_AF_LA_c2.0_cond2.bed$V2
diff_313_AF_RA_vs_25_AF_LA_c2.0_cond2.bed_length<- diff_313_AF_RA_vs_25_AF_LA_c2.0_cond2.bed$V3 - diff_313_AF_RA_vs_25_AF_LA_c2.0_cond2.bed$V2
diff_313_AF_RA_vs_28_AF_LA_c2.0_cond2.bed_length<- diff_313_AF_RA_vs_28_AF_LA_c2.0_cond2.bed$V3 - diff_313_AF_RA_vs_28_AF_LA_c2.0_cond2.bed$V2
diff_313_AF_RA_vs_31_AF_LA_c2.0_cond2.bed_length<- diff_313_AF_RA_vs_31_AF_LA_c2.0_cond2.bed$V3 - diff_313_AF_RA_vs_31_AF_LA_c2.0_cond2.bed$V2
diff_313_AF_RA_vs_313_AF_LA_c2.0_cond2.bed_length<- diff_313_AF_RA_vs_313_AF_LA_c2.0_cond2.bed$V3 - diff_313_AF_RA_vs_313_AF_LA_c2.0_cond2.bed$V2
diff_313_AF_RA_vs_47_AF_LA_c2.0_cond2.bed_length<- diff_313_AF_RA_vs_47_AF_LA_c2.0_cond2.bed$V3 - diff_313_AF_RA_vs_47_AF_LA_c2.0_cond2.bed$V2
diff_33_SR_RA_vs_2231_AF_LA_c2.0_cond2.bed_length<- diff_33_SR_RA_vs_2231_AF_LA_c2.0_cond2.bed$V3 - diff_33_SR_RA_vs_2231_AF_LA_c2.0_cond2.bed$V2
diff_33_SR_RA_vs_25_AF_LA_c2.0_cond2.bed_length<- diff_33_SR_RA_vs_25_AF_LA_c2.0_cond2.bed$V3 - diff_33_SR_RA_vs_25_AF_LA_c2.0_cond2.bed$V2
diff_33_SR_RA_vs_28_AF_LA_c2.0_cond2.bed_length<- diff_33_SR_RA_vs_28_AF_LA_c2.0_cond2.bed$V3 - diff_33_SR_RA_vs_28_AF_LA_c2.0_cond2.bed$V2
diff_33_SR_RA_vs_31_AF_LA_c2.0_cond2.bed_length<- diff_33_SR_RA_vs_31_AF_LA_c2.0_cond2.bed$V3 - diff_33_SR_RA_vs_31_AF_LA_c2.0_cond2.bed$V2
diff_33_SR_RA_vs_313_AF_LA_c2.0_cond2.bed_length<- diff_33_SR_RA_vs_313_AF_LA_c2.0_cond2.bed$V3 - diff_33_SR_RA_vs_313_AF_LA_c2.0_cond2.bed$V2
diff_33_SR_RA_vs_47_AF_LA_c2.0_cond2.bed_length<- diff_33_SR_RA_vs_47_AF_LA_c2.0_cond2.bed$V3 - diff_33_SR_RA_vs_47_AF_LA_c2.0_cond2.bed$V2
diff_47_AF_LA_vs_12b_SR_LA_c2.0_cond1.bed_length<- diff_47_AF_LA_vs_12b_SR_LA_c2.0_cond1.bed$V3 - diff_47_AF_LA_vs_12b_SR_LA_c2.0_cond1.bed$V2
diff_47_AF_LA_vs_20_SR_LA_c2.0_cond1.bed_length<- diff_47_AF_LA_vs_20_SR_LA_c2.0_cond1.bed$V3 - diff_47_AF_LA_vs_20_SR_LA_c2.0_cond1.bed$V2
diff_47_AF_LA_vs_27_SR_LA_c2.0_cond1.bed_length<- diff_47_AF_LA_vs_27_SR_LA_c2.0_cond1.bed$V3 - diff_47_AF_LA_vs_27_SR_LA_c2.0_cond1.bed$V2
diff_47_AF_LA_vs_45_SR_LA_c2.0_cond1.bed_length<- diff_47_AF_LA_vs_45_SR_LA_c2.0_cond1.bed$V3 - diff_47_AF_LA_vs_45_SR_LA_c2.0_cond1.bed$V2



diff_12b_SR_RA_vs_2231_AF_LA_c2.0_cond2.bed_length <- melt(diff_12b_SR_RA_vs_2231_AF_LA_c2.0_cond2.bed_length)
diff_12b_SR_RA_vs_25_AF_LA_c2.0_cond2.bed_length <- melt(diff_12b_SR_RA_vs_25_AF_LA_c2.0_cond2.bed_length)
diff_12b_SR_RA_vs_28_AF_LA_c2.0_cond2.bed_length <- melt(diff_12b_SR_RA_vs_28_AF_LA_c2.0_cond2.bed_length)
diff_12b_SR_RA_vs_31_AF_LA_c2.0_cond2.bed_length<- melt(diff_12b_SR_RA_vs_31_AF_LA_c2.0_cond2.bed_length)
diff_12b_SR_RA_vs_313_AF_LA_c2.0_cond2.bed_length<- melt(diff_12b_SR_RA_vs_313_AF_LA_c2.0_cond2.bed_length)
diff_12b_SR_RA_vs_47_AF_LA_c2.0_cond2.bed_length<- melt(diff_12b_SR_RA_vs_47_AF_LA_c2.0_cond2.bed_length)
diff_18_AF_RA_vs_2231_AF_LA_c2.0_cond2.bed_length<- melt(diff_18_AF_RA_vs_2231_AF_LA_c2.0_cond2.bed_length)
diff_18_AF_RA_vs_25_AF_LA_c2.0_cond2.bed_length<- melt(diff_18_AF_RA_vs_25_AF_LA_c2.0_cond2.bed_length)
diff_18_AF_RA_vs_28_AF_LA_c2.0_cond2.bed_length<- melt(diff_18_AF_RA_vs_28_AF_LA_c2.0_cond2.bed_length)
diff_18_AF_RA_vs_31_AF_LA_c2.0_cond2.bed_length<- melt(diff_18_AF_RA_vs_31_AF_LA_c2.0_cond2.bed_length)
diff_18_AF_RA_vs_313_AF_LA_c2.0_cond2.bed_length<- melt(diff_18_AF_RA_vs_313_AF_LA_c2.0_cond2.bed_length)
diff_18_AF_RA_vs_47_AF_LA_c2.0_cond2.bed_length<- melt(diff_18_AF_RA_vs_47_AF_LA_c2.0_cond2.bed_length)
diff_20_SR_RA_vs_2231_AF_LA_c2.0_cond2.bed_length<- melt(diff_20_SR_RA_vs_2231_AF_LA_c2.0_cond2.bed_length)
diff_20_SR_RA_vs_25_AF_LA_c2.0_cond2.bed_length<- melt(diff_20_SR_RA_vs_25_AF_LA_c2.0_cond2.bed_length)
diff_20_SR_RA_vs_28_AF_LA_c2.0_cond2.bed_length<- melt(diff_20_SR_RA_vs_28_AF_LA_c2.0_cond2.bed_length)
diff_20_SR_RA_vs_31_AF_LA_c2.0_cond2.bed_length<- melt(diff_20_SR_RA_vs_31_AF_LA_c2.0_cond2.bed_length)
diff_20_SR_RA_vs_313_AF_LA_c2.0_cond2.bed_length<- melt(diff_20_SR_RA_vs_313_AF_LA_c2.0_cond2.bed_length)
diff_20_SR_RA_vs_47_AF_LA_c2.0_cond2.bed_length<- melt(diff_20_SR_RA_vs_47_AF_LA_c2.0_cond2.bed_length)
diff_2231_AF_LA_vs_12b_SR_LA_c2.0_cond1.bed_length<- melt(diff_2231_AF_LA_vs_12b_SR_LA_c2.0_cond1.bed_length)
diff_2231_AF_LA_vs_20_SR_LA_c2.0_cond1.bed_length<- melt(diff_2231_AF_LA_vs_20_SR_LA_c2.0_cond1.bed_length)
diff_2231_AF_LA_vs_27_SR_LA_c2.0_cond1.bed_length<- melt(diff_2231_AF_LA_vs_27_SR_LA_c2.0_cond1.bed_length)
diff_2231_AF_LA_vs_45_SR_LA_c2.0_cond1.bed_length<- melt(diff_2231_AF_LA_vs_45_SR_LA_c2.0_cond1.bed_length)
diff_2231_AF_RA_vs_2231_AF_LA_c2.0_cond2.bed_length<- melt(diff_2231_AF_RA_vs_2231_AF_LA_c2.0_cond2.bed_length)
diff_2231_AF_RA_vs_25_AF_LA_c2.0_cond2.bed_length<- melt(diff_2231_AF_RA_vs_25_AF_LA_c2.0_cond2.bed_length)
diff_2231_AF_RA_vs_28_AF_LA_c2.0_cond2.bed_length<- melt(diff_2231_AF_RA_vs_28_AF_LA_c2.0_cond2.bed_length)
diff_2231_AF_RA_vs_31_AF_LA_c2.0_cond2.bed_length<- melt(diff_2231_AF_RA_vs_28_AF_LA_c2.0_cond2.bed_length)
diff_2231_AF_RA_vs_31_AF_LA_c2.0_cond2.bed_length = subset(diff_2231_AF_RA_vs_31_AF_LA_c2.0_cond2.bed_length, select = -c(variable) )
diff_2231_AF_RA_vs_313_AF_LA_c2.0_cond2.bed_length<- melt(diff_2231_AF_RA_vs_313_AF_LA_c2.0_cond2.bed_length)
diff_2231_AF_RA_vs_47_AF_LA_c2.0_cond2.bed_length<- melt(diff_2231_AF_RA_vs_47_AF_LA_c2.0_cond2.bed_length)
diff_25_AF_LA_vs_12b_SR_LA_c2.0_cond1.bed_length<- melt(diff_25_AF_LA_vs_12b_SR_LA_c2.0_cond1.bed_length)
diff_25_AF_LA_vs_20_SR_LA_c2.0_cond1.bed_length<-melt(diff_25_AF_LA_vs_20_SR_LA_c2.0_cond1.bed_length)
diff_25_AF_LA_vs_27_SR_LA_c2.0_cond1.bed_length<- melt(diff_25_AF_LA_vs_27_SR_LA_c2.0_cond1.bed_length)
diff_25_AF_LA_vs_45_SR_LA_c2.0_cond1.bed_length<- melt(diff_25_AF_LA_vs_45_SR_LA_c2.0_cond1.bed_length)
diff_25_AF_RA_vs_2231_AF_LA_c2.0_cond2.bed_length<- melt(diff_25_AF_RA_vs_2231_AF_LA_c2.0_cond2.bed_length)
diff_25_AF_RA_vs_25_AF_LA_c2.0_cond2.bed_length<- melt(diff_25_AF_RA_vs_25_AF_LA_c2.0_cond2.bed_length)
diff_25_AF_RA_vs_28_AF_LA_c2.0_cond2.bed_length<- melt(diff_25_AF_RA_vs_28_AF_LA_c2.0_cond2.bed_length)
diff_25_AF_RA_vs_31_AF_LA_c2.0_cond2.bed_length<- melt(diff_25_AF_RA_vs_31_AF_LA_c2.0_cond2.bed_length)
diff_25_AF_RA_vs_313_AF_LA_c2.0_cond2.bed_length<- melt(diff_25_AF_RA_vs_313_AF_LA_c2.0_cond2.bed_length)
diff_25_AF_RA_vs_47_AF_LA_c2.0_cond2.bed_length<- melt(diff_25_AF_RA_vs_47_AF_LA_c2.0_cond2.bed_length)
diff_27_SR_RA_vs_2231_AF_LA_c2.0_cond2.bed_length<- melt(diff_27_SR_RA_vs_2231_AF_LA_c2.0_cond2.bed_length)
diff_27_SR_RA_vs_25_AF_LA_c2.0_cond2.bed_length<- melt(diff_27_SR_RA_vs_25_AF_LA_c2.0_cond2.bed_length)
diff_27_SR_RA_vs_28_AF_LA_c2.0_cond2.bed_length<- melt(diff_27_SR_RA_vs_28_AF_LA_c2.0_cond2.bed_length)
diff_27_SR_RA_vs_31_AF_LA_c2.0_cond2.bed_length<- melt(diff_27_SR_RA_vs_31_AF_LA_c2.0_cond2.bed_length)
diff_27_SR_RA_vs_313_AF_LA_c2.0_cond2.bed_length<- melt(diff_27_SR_RA_vs_313_AF_LA_c2.0_cond2.bed_length)
diff_27_SR_RA_vs_47_AF_LA_c2.0_cond2.bed_length<- melt(diff_27_SR_RA_vs_47_AF_LA_c2.0_cond2.bed_length)
diff_28_AF_LA_vs_12b_SR_LA_c2.0_cond1.bed_length<- melt(diff_28_AF_LA_vs_12b_SR_LA_c2.0_cond1.bed_length)
diff_28_AF_LA_vs_20_SR_LA_c2.0_cond1.bed_length<- melt(diff_28_AF_LA_vs_20_SR_LA_c2.0_cond1.bed_length)
diff_28_AF_LA_vs_27_SR_LA_c2.0_cond1.bed_length<- melt(diff_28_AF_LA_vs_27_SR_LA_c2.0_cond1.bed_length)
diff_28_AF_LA_vs_45_SR_LA_c2.0_cond1.bed_length<- melt(diff_28_AF_LA_vs_45_SR_LA_c2.0_cond1.bed_length)
diff_28_AF_RA_vs_2231_AF_LA_c2.0_cond2.bed_length<- melt(diff_28_AF_RA_vs_2231_AF_LA_c2.0_cond2.bed_length)
diff_28_AF_RA_vs_25_AF_LA_c2.0_cond2.bed_length<- melt(diff_28_AF_RA_vs_25_AF_LA_c2.0_cond2.bed_length)
diff_28_AF_RA_vs_28_AF_LA_c2.0_cond2.bed_length<- melt(diff_28_AF_RA_vs_28_AF_LA_c2.0_cond2.bed_length)
diff_28_AF_RA_vs_31_AF_LA_c2.0_cond2.bed_length<- melt(diff_28_AF_RA_vs_31_AF_LA_c2.0_cond2.bed_length)
diff_28_AF_RA_vs_313_AF_LA_c2.0_cond2.bed_length<- melt(diff_28_AF_RA_vs_313_AF_LA_c2.0_cond2.bed_length)
diff_28_AF_RA_vs_47_AF_LA_c2.0_cond2.bed_length<- melt(diff_28_AF_RA_vs_47_AF_LA_c2.0_cond2.bed_length)
diff_31_AF_LA_vs_12b_SR_LA_c2.0_cond1.bed_length<- melt(diff_31_AF_LA_vs_12b_SR_LA_c2.0_cond1.bed_length)
diff_31_AF_LA_vs_20_SR_LA_c2.0_cond1.bed_length<-melt(diff_31_AF_LA_vs_20_SR_LA_c2.0_cond1.bed_length)
diff_31_AF_LA_vs_27_SR_LA_c2.0_cond1.bed_length<-melt(diff_31_AF_LA_vs_27_SR_LA_c2.0_cond1.bed_length)
diff_31_AF_LA_vs_45_SR_LA_c2.0_cond1.bed_length<- melt(diff_31_AF_LA_vs_45_SR_LA_c2.0_cond1.bed_length)
diff_31_AF_RA_vs_2231_AF_LA_c2.0_cond2.bed_length<- melt(diff_31_AF_RA_vs_2231_AF_LA_c2.0_cond2.bed_length)
diff_31_AF_RA_vs_25_AF_LA_c2.0_cond2.bed_length<-melt(diff_31_AF_RA_vs_25_AF_LA_c2.0_cond2.bed_length)
diff_31_AF_RA_vs_28_AF_LA_c2.0_cond2.bed_length<- melt(diff_31_AF_RA_vs_28_AF_LA_c2.0_cond2.bed_length)
diff_31_AF_RA_vs_31_AF_LA_c2.0_cond2.bed_length<- melt(diff_31_AF_RA_vs_31_AF_LA_c2.0_cond2.bed_length)
diff_31_AF_RA_vs_313_AF_LA_c2.0_cond2.bed_length<- melt(diff_31_AF_RA_vs_313_AF_LA_c2.0_cond2.bed_length)
diff_31_AF_RA_vs_47_AF_LA_c2.0_cond2.bed_length<- melt(diff_31_AF_RA_vs_47_AF_LA_c2.0_cond2.bed_length)
diff_313_AF_LA_vs_12b_SR_LA_c2.0_cond1.bed_length<- melt(diff_313_AF_LA_vs_12b_SR_LA_c2.0_cond1.bed_length)
diff_313_AF_LA_vs_20_SR_LA_c2.0_cond1.bed_length<- melt(diff_313_AF_LA_vs_20_SR_LA_c2.0_cond1.bed_length)
diff_313_AF_LA_vs_27_SR_LA_c2.0_cond1.bed_length<- melt(diff_313_AF_LA_vs_27_SR_LA_c2.0_cond1.bed_length)
diff_313_AF_LA_vs_45_SR_LA_c2.0_cond1.bed_length<- melt(diff_313_AF_LA_vs_45_SR_LA_c2.0_cond1.bed_length)
diff_313_AF_RA_vs_2231_AF_LA_c2.0_cond2.bed_length<- melt(diff_313_AF_RA_vs_2231_AF_LA_c2.0_cond2.bed_length)
diff_313_AF_RA_vs_25_AF_LA_c2.0_cond2.bed_length<- melt(diff_313_AF_RA_vs_25_AF_LA_c2.0_cond2.bed_length)
diff_313_AF_RA_vs_28_AF_LA_c2.0_cond2.bed_length<- melt(diff_313_AF_RA_vs_28_AF_LA_c2.0_cond2.bed_length)
diff_313_AF_RA_vs_31_AF_LA_c2.0_cond2.bed_length<- melt(diff_313_AF_RA_vs_31_AF_LA_c2.0_cond2.bed_length)
diff_313_AF_RA_vs_313_AF_LA_c2.0_cond2.bed_length<- melt(diff_313_AF_RA_vs_313_AF_LA_c2.0_cond2.bed_length)
diff_313_AF_RA_vs_47_AF_LA_c2.0_cond2.bed_length<-melt(diff_313_AF_RA_vs_47_AF_LA_c2.0_cond2.bed_length)
diff_33_SR_RA_vs_2231_AF_LA_c2.0_cond2.bed_length<- melt(diff_33_SR_RA_vs_2231_AF_LA_c2.0_cond2.bed_length)
diff_33_SR_RA_vs_25_AF_LA_c2.0_cond2.bed_length<- melt(diff_33_SR_RA_vs_25_AF_LA_c2.0_cond2.bed_length)
diff_33_SR_RA_vs_28_AF_LA_c2.0_cond2.bed_length<- melt(diff_33_SR_RA_vs_28_AF_LA_c2.0_cond2.bed_length)
diff_33_SR_RA_vs_31_AF_LA_c2.0_cond2.bed_length<- melt(diff_33_SR_RA_vs_31_AF_LA_c2.0_cond2.bed_length)
diff_33_SR_RA_vs_313_AF_LA_c2.0_cond2.bed_length<- melt(diff_33_SR_RA_vs_313_AF_LA_c2.0_cond2.bed_length)
diff_33_SR_RA_vs_47_AF_LA_c2.0_cond2.bed_length<- melt(diff_33_SR_RA_vs_47_AF_LA_c2.0_cond2.bed_length)
diff_47_AF_LA_vs_12b_SR_LA_c2.0_cond1.bed_length<- melt(diff_47_AF_LA_vs_12b_SR_LA_c2.0_cond1.bed_length)
diff_47_AF_LA_vs_20_SR_LA_c2.0_cond1.bed_length<- melt(diff_47_AF_LA_vs_20_SR_LA_c2.0_cond1.bed_length)
diff_47_AF_LA_vs_27_SR_LA_c2.0_cond1.bed_length<- melt(diff_47_AF_LA_vs_27_SR_LA_c2.0_cond1.bed_length)
diff_47_AF_LA_vs_45_SR_LA_c2.0_cond1.bed_length<- melt(diff_47_AF_LA_vs_45_SR_LA_c2.0_cond1.bed_length)


diff_12b_SR_RA_vs_2231_AF_LA_c2.0_cond2.bed_length$group <- "AF_LA"
diff_12b_SR_RA_vs_25_AF_LA_c2.0_cond2.bed_length$group <- "AF_LA"
diff_12b_SR_RA_vs_28_AF_LA_c2.0_cond2.bed_length$group <- "AF_LA"
diff_12b_SR_RA_vs_31_AF_LA_c2.0_cond2.bed_length$group <- "AF_LA"
diff_12b_SR_RA_vs_313_AF_LA_c2.0_cond2.bed_length$group <- "AF_LA"
diff_12b_SR_RA_vs_47_AF_LA_c2.0_cond2.bed_length$group <- "AF_LA"
diff_18_AF_RA_vs_2231_AF_LA_c2.0_cond2.bed_length$group <- "AF_LA"
diff_18_AF_RA_vs_25_AF_LA_c2.0_cond2.bed_length$group <- "AF_LA"
diff_18_AF_RA_vs_28_AF_LA_c2.0_cond2.bed_length$group <- "AF_LA"
diff_18_AF_RA_vs_31_AF_LA_c2.0_cond2.bed_length$group <- "AF_LA"
diff_18_AF_RA_vs_313_AF_LA_c2.0_cond2.bed_length$group <- "AF_LA"
diff_18_AF_RA_vs_47_AF_LA_c2.0_cond2.bed_length$group <- "AF_LA"
diff_20_SR_RA_vs_2231_AF_LA_c2.0_cond2.bed_length$group <- "AF_LA"
diff_20_SR_RA_vs_25_AF_LA_c2.0_cond2.bed_length$group <- "AF_LA"
diff_20_SR_RA_vs_28_AF_LA_c2.0_cond2.bed_length$group <- "AF_LA"
diff_20_SR_RA_vs_31_AF_LA_c2.0_cond2.bed_length$group <- "AF_LA"
diff_20_SR_RA_vs_313_AF_LA_c2.0_cond2.bed_length$group <- "AF_LA"
diff_20_SR_RA_vs_47_AF_LA_c2.0_cond2.bed_length$group <- "AF_LA"
diff_2231_AF_LA_vs_12b_SR_LA_c2.0_cond1.bed_length$group <- "AF_LA"
diff_2231_AF_LA_vs_20_SR_LA_c2.0_cond1.bed_length$group <- "AF_LA"
diff_2231_AF_LA_vs_27_SR_LA_c2.0_cond1.bed_length$group <- "AF_LA"
diff_2231_AF_LA_vs_45_SR_LA_c2.0_cond1.bed_length$group <- "AF_LA"
diff_2231_AF_RA_vs_2231_AF_LA_c2.0_cond2.bed_length$group <- "AF_LA"
diff_2231_AF_RA_vs_25_AF_LA_c2.0_cond2.bed_length$group <- "AF_LA"
diff_2231_AF_RA_vs_28_AF_LA_c2.0_cond2.bed_length$group <- "AF_LA"
diff_2231_AF_RA_vs_31_AF_LA_c2.0_cond2.bed_length$group <- "AF_LA"
diff_2231_AF_RA_vs_313_AF_LA_c2.0_cond2.bed_length$group <- "AF_LA"
diff_2231_AF_RA_vs_47_AF_LA_c2.0_cond2.bed_length$group <- "AF_LA"
diff_25_AF_LA_vs_12b_SR_LA_c2.0_cond1.bed_length$group <- "AF_LA"
diff_25_AF_LA_vs_20_SR_LA_c2.0_cond1.bed_length$group <- "AF_LA"
diff_25_AF_LA_vs_27_SR_LA_c2.0_cond1.bed_length$group <- "AF_LA"
diff_25_AF_LA_vs_45_SR_LA_c2.0_cond1.bed_length$group <- "AF_LA"
diff_25_AF_RA_vs_2231_AF_LA_c2.0_cond2.bed_length$group <- "AF_LA"
diff_25_AF_RA_vs_25_AF_LA_c2.0_cond2.bed_length$group <- "AF_LA"
diff_25_AF_RA_vs_28_AF_LA_c2.0_cond2.bed_length$group <- "AF_LA"
diff_25_AF_RA_vs_31_AF_LA_c2.0_cond2.bed_length$group <- "AF_LA"
diff_25_AF_RA_vs_313_AF_LA_c2.0_cond2.bed_length$group <- "AF_LA"
diff_25_AF_RA_vs_47_AF_LA_c2.0_cond2.bed_length$group <- "AF_LA"
diff_27_SR_RA_vs_2231_AF_LA_c2.0_cond2.bed_length$group <- "AF_LA"
diff_27_SR_RA_vs_25_AF_LA_c2.0_cond2.bed_length$group <- "AF_LA"
diff_27_SR_RA_vs_28_AF_LA_c2.0_cond2.bed_length$group <- "AF_LA"
diff_27_SR_RA_vs_31_AF_LA_c2.0_cond2.bed_length$group <- "AF_LA"
diff_27_SR_RA_vs_313_AF_LA_c2.0_cond2.bed_length$group <- "AF_LA"
diff_27_SR_RA_vs_47_AF_LA_c2.0_cond2.bed_length$group <- "AF_LA"
diff_28_AF_LA_vs_12b_SR_LA_c2.0_cond1.bed_length$group <- "AF_LA"
diff_28_AF_LA_vs_20_SR_LA_c2.0_cond1.bed_length$group <- "AF_LA"
diff_28_AF_LA_vs_27_SR_LA_c2.0_cond1.bed_length$group <- "AF_LA"
diff_28_AF_LA_vs_45_SR_LA_c2.0_cond1.bed_length$group <- "AF_LA"
diff_28_AF_RA_vs_2231_AF_LA_c2.0_cond2.bed_length$group <- "AF_LA"
diff_28_AF_RA_vs_25_AF_LA_c2.0_cond2.bed_length$group <- "AF_LA"
diff_28_AF_RA_vs_28_AF_LA_c2.0_cond2.bed_length$group <- "AF_LA"
diff_28_AF_RA_vs_31_AF_LA_c2.0_cond2.bed_length$group <- "AF_LA"
diff_28_AF_RA_vs_313_AF_LA_c2.0_cond2.bed_length$group <- "AF_LA"
diff_28_AF_RA_vs_47_AF_LA_c2.0_cond2.bed_length$group <- "AF_LA"
diff_31_AF_LA_vs_12b_SR_LA_c2.0_cond1.bed_length$group <- "AF_LA"
diff_31_AF_LA_vs_20_SR_LA_c2.0_cond1.bed_length$group <- "AF_LA"
diff_31_AF_LA_vs_27_SR_LA_c2.0_cond1.bed_length$group <- "AF_LA"
diff_31_AF_LA_vs_45_SR_LA_c2.0_cond1.bed_length$group <- "AF_LA"
diff_31_AF_RA_vs_2231_AF_LA_c2.0_cond2.bed_length$group <- "AF_LA"
diff_31_AF_RA_vs_25_AF_LA_c2.0_cond2.bed_length$group <- "AF_LA"
diff_31_AF_RA_vs_28_AF_LA_c2.0_cond2.bed_length$group <- "AF_LA"
diff_31_AF_RA_vs_31_AF_LA_c2.0_cond2.bed_length$group <- "AF_LA"
diff_31_AF_RA_vs_313_AF_LA_c2.0_cond2.bed_length$group <- "AF_LA"
diff_31_AF_RA_vs_47_AF_LA_c2.0_cond2.bed_length$group <- "AF_LA"
diff_313_AF_LA_vs_12b_SR_LA_c2.0_cond1.bed_length$group <- "AF_LA"
diff_313_AF_LA_vs_20_SR_LA_c2.0_cond1.bed_length$group <- "AF_LA"
diff_313_AF_LA_vs_27_SR_LA_c2.0_cond1.bed_length$group <- "AF_LA"
diff_313_AF_LA_vs_45_SR_LA_c2.0_cond1.bed_length$group <- "AF_LA"
diff_313_AF_RA_vs_2231_AF_LA_c2.0_cond2.bed_length$group <- "AF_LA"
diff_313_AF_RA_vs_25_AF_LA_c2.0_cond2.bed_length$group <- "AF_LA"
diff_313_AF_RA_vs_28_AF_LA_c2.0_cond2.bed_length$group <- "AF_LA"
diff_313_AF_RA_vs_31_AF_LA_c2.0_cond2.bed_length$group <- "AF_LA"
diff_313_AF_RA_vs_313_AF_LA_c2.0_cond2.bed_length$group <- "AF_LA"
diff_313_AF_RA_vs_47_AF_LA_c2.0_cond2.bed_length$group <- "AF_LA"
diff_33_SR_RA_vs_2231_AF_LA_c2.0_cond2.bed_length$group <- "AF_LA"
diff_33_SR_RA_vs_25_AF_LA_c2.0_cond2.bed_length$group <- "AF_LA"
diff_33_SR_RA_vs_28_AF_LA_c2.0_cond2.bed_length$group <- "AF_LA"
diff_33_SR_RA_vs_31_AF_LA_c2.0_cond2.bed_length$group <- "AF_LA"
diff_33_SR_RA_vs_313_AF_LA_c2.0_cond2.bed_length$group <- "AF_LA"
diff_33_SR_RA_vs_47_AF_LA_c2.0_cond2.bed_length$group <- "AF_LA"
diff_47_AF_LA_vs_12b_SR_LA_c2.0_cond1.bed_length$group <- "AF_LA"
diff_47_AF_LA_vs_20_SR_LA_c2.0_cond1.bed_length$group <- "AF_LA"
diff_47_AF_LA_vs_27_SR_LA_c2.0_cond1.bed_length$group <- "AF_LA"
diff_47_AF_LA_vs_45_SR_LA_c2.0_cond1.bed_length$group <- "AF_LA"

diff_12b_SR_RA_vs_2231_AF_LA_c2.0_cond2.bed_length$sample <- "diff_12b_SR_RA_vs_2231_AF_LA_c2.0_cond2.bed"
diff_12b_SR_RA_vs_25_AF_LA_c2.0_cond2.bed_length$sample <- "diff_12b_SR_RA_vs_25_AF_LA_c2.0_cond2.bed"
diff_12b_SR_RA_vs_28_AF_LA_c2.0_cond2.bed_length$sample <- "diff_12b_SR_RA_vs_28_AF_LA_c2.0_cond2.bed"
diff_12b_SR_RA_vs_31_AF_LA_c2.0_cond2.bed_length$sample <- "diff_12b_SR_RA_vs_31_AF_LA_c2.0_cond2.bed"
diff_12b_SR_RA_vs_313_AF_LA_c2.0_cond2.bed_length$sample <- "diff_12b_SR_RA_vs_313_AF_LA_c2.0_cond2.bed"
diff_12b_SR_RA_vs_47_AF_LA_c2.0_cond2.bed_length$sample <- "diff_12b_SR_RA_vs_47_AF_LA_c2.0_cond2.bed"
diff_18_AF_RA_vs_2231_AF_LA_c2.0_cond2.bed_length$sample <- "diff_18_AF_RA_vs_2231_AF_LA_c2.0_cond2.bed"
diff_18_AF_RA_vs_25_AF_LA_c2.0_cond2.bed_length$sample <- "diff_18_AF_RA_vs_25_AF_LA_c2.0_cond2.bed"
diff_18_AF_RA_vs_28_AF_LA_c2.0_cond2.bed_length$sample <- "diff_18_AF_RA_vs_28_AF_LA_c2.0_cond2.bed"
diff_18_AF_RA_vs_31_AF_LA_c2.0_cond2.bed_length$sample <- "diff_18_AF_RA_vs_31_AF_LA_c2.0_cond2.bed"
diff_18_AF_RA_vs_313_AF_LA_c2.0_cond2.bed_length$sample <- "diff_18_AF_RA_vs_313_AF_LA_c2.0_cond2.bed"
diff_18_AF_RA_vs_47_AF_LA_c2.0_cond2.bed_length$sample <- "diff_18_AF_RA_vs_47_AF_LA_c2.0_cond2.bed"
diff_20_SR_RA_vs_2231_AF_LA_c2.0_cond2.bed_length$sample <- "diff_20_SR_RA_vs_2231_AF_LA_c2.0_cond2.bed"
diff_20_SR_RA_vs_25_AF_LA_c2.0_cond2.bed_length$sample <- "diff_20_SR_RA_vs_25_AF_LA_c2.0_cond2.bed"
diff_20_SR_RA_vs_28_AF_LA_c2.0_cond2.bed_length$sample <- "diff_20_SR_RA_vs_28_AF_LA_c2.0_cond2.bed"
diff_20_SR_RA_vs_31_AF_LA_c2.0_cond2.bed_length$sample <- "diff_20_SR_RA_vs_31_AF_LA_c2.0_cond2.bed"
diff_20_SR_RA_vs_313_AF_LA_c2.0_cond2.bed_length$sample <- "diff_20_SR_RA_vs_313_AF_LA_c2.0_cond2.bed"
diff_20_SR_RA_vs_47_AF_LA_c2.0_cond2.bed_length$sample <- "diff_20_SR_RA_vs_47_AF_LA_c2.0_cond2.bed"
diff_2231_AF_LA_vs_12b_SR_LA_c2.0_cond1.bed_length$sample <- "diff_2231_AF_LA_vs_12b_SR_LA_c2.0_cond1.bed"
diff_2231_AF_LA_vs_20_SR_LA_c2.0_cond1.bed_length$sample <- "diff_2231_AF_LA_vs_20_SR_LA_c2.0_cond1.bed"
diff_2231_AF_LA_vs_27_SR_LA_c2.0_cond1.bed_length$sample <- "diff_2231_AF_LA_vs_27_SR_LA_c2.0_cond1.bed"
diff_2231_AF_LA_vs_45_SR_LA_c2.0_cond1.bed_length$sample <- "diff_2231_AF_LA_vs_45_SR_LA_c2.0_cond1.bed"
diff_2231_AF_RA_vs_2231_AF_LA_c2.0_cond2.bed_length$sample <- "diff_2231_AF_RA_vs_2231_AF_LA_c2.0_cond2.bed"
diff_2231_AF_RA_vs_25_AF_LA_c2.0_cond2.bed_length$sample <- "diff_2231_AF_RA_vs_25_AF_LA_c2.0_cond2.bed"
diff_2231_AF_RA_vs_28_AF_LA_c2.0_cond2.bed_length$sample <- "diff_2231_AF_RA_vs_28_AF_LA_c2.0_cond2.bed"
diff_2231_AF_RA_vs_31_AF_LA_c2.0_cond2.bed_length$sample <- "diff_2231_AF_RA_vs_31_AF_LA_c2.0_cond2.bed"
diff_2231_AF_RA_vs_313_AF_LA_c2.0_cond2.bed_length$sample <- "diff_2231_AF_RA_vs_313_AF_LA_c2.0_cond2.bed"
diff_2231_AF_RA_vs_47_AF_LA_c2.0_cond2.bed_length$sample <- "diff_2231_AF_RA_vs_47_AF_LA_c2.0_cond2.bed"
diff_25_AF_LA_vs_12b_SR_LA_c2.0_cond1.bed_length$sample <- "diff_25_AF_LA_vs_12b_SR_LA_c2.0_cond1.bed"
diff_25_AF_LA_vs_20_SR_LA_c2.0_cond1.bed_length$sample <- "diff_25_AF_LA_vs_20_SR_LA_c2.0_cond1.bed"
diff_25_AF_LA_vs_27_SR_LA_c2.0_cond1.bed_length$sample <- "diff_25_AF_LA_vs_27_SR_LA_c2.0_cond1.bed"
diff_25_AF_LA_vs_45_SR_LA_c2.0_cond1.bed_length$sample <- "diff_25_AF_LA_vs_45_SR_LA_c2.0_cond1.bed"
diff_25_AF_RA_vs_2231_AF_LA_c2.0_cond2.bed_length$sample <- "diff_25_AF_RA_vs_2231_AF_LA_c2.0_cond2.bed"
diff_25_AF_RA_vs_25_AF_LA_c2.0_cond2.bed_length$sample <- "diff_25_AF_RA_vs_25_AF_LA_c2.0_cond2.bed"
diff_25_AF_RA_vs_28_AF_LA_c2.0_cond2.bed_length$sample <- "diff_25_AF_RA_vs_28_AF_LA_c2.0_cond2.bed"
diff_25_AF_RA_vs_31_AF_LA_c2.0_cond2.bed_length$sample <- "diff_25_AF_RA_vs_31_AF_LA_c2.0_cond2.bed"
diff_25_AF_RA_vs_313_AF_LA_c2.0_cond2.bed_length$sample <- "diff_25_AF_RA_vs_313_AF_LA_c2.0_cond2.bed"
diff_25_AF_RA_vs_47_AF_LA_c2.0_cond2.bed_length$sample <- "diff_25_AF_RA_vs_47_AF_LA_c2.0_cond2.bed"
diff_27_SR_RA_vs_2231_AF_LA_c2.0_cond2.bed_length$sample <- "diff_27_SR_RA_vs_2231_AF_LA_c2.0_cond2.bed"
diff_27_SR_RA_vs_25_AF_LA_c2.0_cond2.bed_length$sample <- "diff_27_SR_RA_vs_25_AF_LA_c2.0_cond2.bed"
diff_27_SR_RA_vs_28_AF_LA_c2.0_cond2.bed_length$sample <- "diff_27_SR_RA_vs_28_AF_LA_c2.0_cond2.bed"
diff_27_SR_RA_vs_31_AF_LA_c2.0_cond2.bed_length$sample <- "diff_27_SR_RA_vs_31_AF_LA_c2.0_cond2.bed"
diff_27_SR_RA_vs_313_AF_LA_c2.0_cond2.bed_length$sample <- "diff_27_SR_RA_vs_313_AF_LA_c2.0_cond2.bed"
diff_27_SR_RA_vs_47_AF_LA_c2.0_cond2.bed_length$sample <- "diff_27_SR_RA_vs_47_AF_LA_c2.0_cond2.bed"
diff_28_AF_LA_vs_12b_SR_LA_c2.0_cond1.bed_length$sample <- "diff_28_AF_LA_vs_12b_SR_LA_c2.0_cond1.bed"
diff_28_AF_LA_vs_20_SR_LA_c2.0_cond1.bed_length$sample <- "diff_28_AF_LA_vs_20_SR_LA_c2.0_cond1.bed"
diff_28_AF_LA_vs_27_SR_LA_c2.0_cond1.bed_length$sample <- "diff_28_AF_LA_vs_27_SR_LA_c2.0_cond1.bed"
diff_28_AF_LA_vs_45_SR_LA_c2.0_cond1.bed_length$sample <- "diff_28_AF_LA_vs_45_SR_LA_c2.0_cond1.bed"
diff_28_AF_RA_vs_2231_AF_LA_c2.0_cond2.bed_length$sample <- "diff_28_AF_RA_vs_2231_AF_LA_c2.0_cond2.bed"
diff_28_AF_RA_vs_25_AF_LA_c2.0_cond2.bed_length$sample <- "diff_28_AF_RA_vs_25_AF_LA_c2.0_cond2.bed"
diff_28_AF_RA_vs_28_AF_LA_c2.0_cond2.bed_length$sample <- "diff_28_AF_RA_vs_28_AF_LA_c2.0_cond2.bed"
diff_28_AF_RA_vs_31_AF_LA_c2.0_cond2.bed_length$sample <- "diff_28_AF_RA_vs_31_AF_LA_c2.0_cond2.bed"
diff_28_AF_RA_vs_313_AF_LA_c2.0_cond2.bed_length$sample <- "diff_28_AF_RA_vs_313_AF_LA_c2.0_cond2.bed"
diff_28_AF_RA_vs_47_AF_LA_c2.0_cond2.bed_length$sample <- "diff_28_AF_RA_vs_47_AF_LA_c2.0_cond2.bed"
diff_31_AF_LA_vs_12b_SR_LA_c2.0_cond1.bed_length$sample <- "diff_31_AF_LA_vs_12b_SR_LA_c2.0_cond1.bed"
diff_31_AF_LA_vs_20_SR_LA_c2.0_cond1.bed_length$sample <- "diff_31_AF_LA_vs_20_SR_LA_c2.0_cond1.bed"
diff_31_AF_LA_vs_27_SR_LA_c2.0_cond1.bed_length$sample <- "diff_31_AF_LA_vs_27_SR_LA_c2.0_cond1.bed"
diff_31_AF_LA_vs_45_SR_LA_c2.0_cond1.bed_length$sample <- "diff_31_AF_LA_vs_45_SR_LA_c2.0_cond1.bed"
diff_31_AF_RA_vs_2231_AF_LA_c2.0_cond2.bed_length$sample <- "diff_31_AF_RA_vs_2231_AF_LA_c2.0_cond2.bed"
diff_31_AF_RA_vs_25_AF_LA_c2.0_cond2.bed_length$sample <- "diff_31_AF_RA_vs_25_AF_LA_c2.0_cond2.bed"
diff_31_AF_RA_vs_28_AF_LA_c2.0_cond2.bed_length$sample <- "diff_31_AF_RA_vs_28_AF_LA_c2.0_cond2.bed"
diff_31_AF_RA_vs_31_AF_LA_c2.0_cond2.bed_length$sample <- "diff_31_AF_RA_vs_31_AF_LA_c2.0_cond2.bed"
diff_31_AF_RA_vs_313_AF_LA_c2.0_cond2.bed_length$sample <- "diff_31_AF_RA_vs_313_AF_LA_c2.0_cond2.bed"
diff_31_AF_RA_vs_47_AF_LA_c2.0_cond2.bed_length$sample <- "diff_31_AF_RA_vs_47_AF_LA_c2.0_cond2.bed"
diff_313_AF_LA_vs_12b_SR_LA_c2.0_cond1.bed_length$sample <- "diff_313_AF_LA_vs_12b_SR_LA_c2.0_cond1.bed"
diff_313_AF_LA_vs_20_SR_LA_c2.0_cond1.bed_length$sample <- "diff_313_AF_LA_vs_20_SR_LA_c2.0_cond1.bed"
diff_313_AF_LA_vs_27_SR_LA_c2.0_cond1.bed_length$sample <- "diff_313_AF_LA_vs_27_SR_LA_c2.0_cond1.bed"
diff_313_AF_LA_vs_45_SR_LA_c2.0_cond1.bed_length$sample <- "diff_313_AF_LA_vs_45_SR_LA_c2.0_cond1.bed"
diff_313_AF_RA_vs_2231_AF_LA_c2.0_cond2.bed_length$sample <- "diff_313_AF_RA_vs_2231_AF_LA_c2.0_cond2.bed"
diff_313_AF_RA_vs_25_AF_LA_c2.0_cond2.bed_length$sample <- "diff_313_AF_RA_vs_25_AF_LA_c2.0_cond2.bed"
diff_313_AF_RA_vs_28_AF_LA_c2.0_cond2.bed_length$sample <- "diff_313_AF_RA_vs_28_AF_LA_c2.0_cond2.bed"
diff_313_AF_RA_vs_31_AF_LA_c2.0_cond2.bed_length$sample <- "diff_313_AF_RA_vs_31_AF_LA_c2.0_cond2.bed"
diff_313_AF_RA_vs_313_AF_LA_c2.0_cond2.bed_length$sample <- "diff_313_AF_RA_vs_313_AF_LA_c2.0_cond2.bed"
diff_313_AF_RA_vs_47_AF_LA_c2.0_cond2.bed_length$sample <- "diff_313_AF_RA_vs_47_AF_LA_c2.0_cond2.bed"
diff_33_SR_RA_vs_2231_AF_LA_c2.0_cond2.bed_length$sample <- "diff_33_SR_RA_vs_2231_AF_LA_c2.0_cond2.bed"
diff_33_SR_RA_vs_25_AF_LA_c2.0_cond2.bed_length$sample <- "diff_33_SR_RA_vs_25_AF_LA_c2.0_cond2.bed"
diff_33_SR_RA_vs_28_AF_LA_c2.0_cond2.bed_length$sample <- "diff_33_SR_RA_vs_28_AF_LA_c2.0_cond2.bed"
diff_33_SR_RA_vs_31_AF_LA_c2.0_cond2.bed_length$sample <- "diff_33_SR_RA_vs_31_AF_LA_c2.0_cond2.bed"
diff_33_SR_RA_vs_313_AF_LA_c2.0_cond2.bed_length$sample <- "diff_33_SR_RA_vs_313_AF_LA_c2.0_cond2.bed"
diff_33_SR_RA_vs_47_AF_LA_c2.0_cond2.bed_length$sample <- "diff_33_SR_RA_vs_47_AF_LA_c2.0_cond2.bed"
diff_47_AF_LA_vs_12b_SR_LA_c2.0_cond1.bed_length$sample <- "diff_47_AF_LA_vs_12b_SR_LA_c2.0_cond1.bed"
diff_47_AF_LA_vs_20_SR_LA_c2.0_cond1.bed_length$sample <- "diff_47_AF_LA_vs_20_SR_LA_c2.0_cond1.bed"
diff_47_AF_LA_vs_27_SR_LA_c2.0_cond1.bed_length$sample <- "diff_47_AF_LA_vs_27_SR_LA_c2.0_cond1.bed"
diff_47_AF_LA_vs_45_SR_LA_c2.0_cond1.bed_length$sample <- "diff_47_AF_LA_vs_45_SR_LA_c2.0_cond1.bed"


peaksets_lengths_AF_LA_pairwise <- rbind(diff_12b_SR_RA_vs_2231_AF_LA_c2.0_cond2.bed_length,
                                         diff_12b_SR_RA_vs_25_AF_LA_c2.0_cond2.bed_length,
                                         diff_12b_SR_RA_vs_28_AF_LA_c2.0_cond2.bed_length,
                                         diff_12b_SR_RA_vs_31_AF_LA_c2.0_cond2.bed_length,
                                         diff_12b_SR_RA_vs_313_AF_LA_c2.0_cond2.bed_length,
                                         diff_12b_SR_RA_vs_47_AF_LA_c2.0_cond2.bed_length,
                                         diff_18_AF_RA_vs_2231_AF_LA_c2.0_cond2.bed_length,
                                         diff_18_AF_RA_vs_25_AF_LA_c2.0_cond2.bed_length,
                                         diff_18_AF_RA_vs_28_AF_LA_c2.0_cond2.bed_length,
                                         diff_18_AF_RA_vs_31_AF_LA_c2.0_cond2.bed_length,
                                         diff_18_AF_RA_vs_313_AF_LA_c2.0_cond2.bed_length,
                                         diff_18_AF_RA_vs_47_AF_LA_c2.0_cond2.bed_length,
                                         diff_20_SR_RA_vs_2231_AF_LA_c2.0_cond2.bed_length,
                                         diff_20_SR_RA_vs_25_AF_LA_c2.0_cond2.bed_length,
                                         diff_20_SR_RA_vs_28_AF_LA_c2.0_cond2.bed_length,
                                         diff_20_SR_RA_vs_31_AF_LA_c2.0_cond2.bed_length,
                                         diff_20_SR_RA_vs_313_AF_LA_c2.0_cond2.bed_length,
                                         diff_20_SR_RA_vs_47_AF_LA_c2.0_cond2.bed_length,
                                         diff_2231_AF_LA_vs_12b_SR_LA_c2.0_cond1.bed_length,
                                         diff_2231_AF_LA_vs_20_SR_LA_c2.0_cond1.bed_length,
                                         diff_2231_AF_LA_vs_27_SR_LA_c2.0_cond1.bed_length,
                                         diff_2231_AF_LA_vs_45_SR_LA_c2.0_cond1.bed_length,
                                         diff_2231_AF_RA_vs_2231_AF_LA_c2.0_cond2.bed_length,
                                         diff_2231_AF_RA_vs_25_AF_LA_c2.0_cond2.bed_length,
                                         diff_2231_AF_RA_vs_28_AF_LA_c2.0_cond2.bed_length,
                                         diff_2231_AF_RA_vs_31_AF_LA_c2.0_cond2.bed_length,
                                         diff_2231_AF_RA_vs_313_AF_LA_c2.0_cond2.bed_length,
                                         diff_2231_AF_RA_vs_47_AF_LA_c2.0_cond2.bed_length,
                                         diff_25_AF_LA_vs_12b_SR_LA_c2.0_cond1.bed_length,
                                         diff_25_AF_LA_vs_20_SR_LA_c2.0_cond1.bed_length,
                                         diff_25_AF_LA_vs_27_SR_LA_c2.0_cond1.bed_length,
                                         diff_25_AF_LA_vs_45_SR_LA_c2.0_cond1.bed_length,
                                         diff_25_AF_RA_vs_2231_AF_LA_c2.0_cond2.bed_length,
                                         diff_25_AF_RA_vs_25_AF_LA_c2.0_cond2.bed_length,
                                         diff_25_AF_RA_vs_28_AF_LA_c2.0_cond2.bed_length,
                                         diff_25_AF_RA_vs_31_AF_LA_c2.0_cond2.bed_length,
                                         diff_25_AF_RA_vs_313_AF_LA_c2.0_cond2.bed_length,
                                         diff_25_AF_RA_vs_47_AF_LA_c2.0_cond2.bed_length,
                                         diff_27_SR_RA_vs_2231_AF_LA_c2.0_cond2.bed_length,
                                         diff_27_SR_RA_vs_25_AF_LA_c2.0_cond2.bed_length,
                                         diff_27_SR_RA_vs_28_AF_LA_c2.0_cond2.bed_length,
                                         diff_27_SR_RA_vs_31_AF_LA_c2.0_cond2.bed_length,
                                         diff_27_SR_RA_vs_313_AF_LA_c2.0_cond2.bed_length,
                                         diff_27_SR_RA_vs_47_AF_LA_c2.0_cond2.bed_length,
                                         diff_28_AF_LA_vs_12b_SR_LA_c2.0_cond1.bed_length,
                                         diff_28_AF_LA_vs_20_SR_LA_c2.0_cond1.bed_length,
                                         diff_28_AF_LA_vs_27_SR_LA_c2.0_cond1.bed_length,
                                         diff_28_AF_LA_vs_45_SR_LA_c2.0_cond1.bed_length,
                                         diff_28_AF_RA_vs_2231_AF_LA_c2.0_cond2.bed_length,
                                         diff_28_AF_RA_vs_25_AF_LA_c2.0_cond2.bed_length,
                                         diff_28_AF_RA_vs_28_AF_LA_c2.0_cond2.bed_length,
                                         diff_28_AF_RA_vs_31_AF_LA_c2.0_cond2.bed_length,
                                         diff_28_AF_RA_vs_313_AF_LA_c2.0_cond2.bed_length,
                                         diff_28_AF_RA_vs_47_AF_LA_c2.0_cond2.bed_length,
                                         diff_31_AF_LA_vs_12b_SR_LA_c2.0_cond1.bed_length,
                                         diff_31_AF_LA_vs_20_SR_LA_c2.0_cond1.bed_length,
                                         diff_31_AF_LA_vs_27_SR_LA_c2.0_cond1.bed_length,
                                         diff_31_AF_LA_vs_45_SR_LA_c2.0_cond1.bed_length,
                                         diff_31_AF_RA_vs_2231_AF_LA_c2.0_cond2.bed_length,
                                         diff_31_AF_RA_vs_25_AF_LA_c2.0_cond2.bed_length,
                                         diff_31_AF_RA_vs_28_AF_LA_c2.0_cond2.bed_length,
                                         diff_31_AF_RA_vs_31_AF_LA_c2.0_cond2.bed_length,
                                         diff_31_AF_RA_vs_313_AF_LA_c2.0_cond2.bed_length,
                                         diff_31_AF_RA_vs_47_AF_LA_c2.0_cond2.bed_length,
                                         diff_313_AF_LA_vs_12b_SR_LA_c2.0_cond1.bed_length,
                                         diff_313_AF_LA_vs_20_SR_LA_c2.0_cond1.bed_length,
                                         diff_313_AF_LA_vs_27_SR_LA_c2.0_cond1.bed_length,
                                         diff_313_AF_LA_vs_45_SR_LA_c2.0_cond1.bed_length,
                                         diff_313_AF_RA_vs_2231_AF_LA_c2.0_cond2.bed_length,
                                         diff_313_AF_RA_vs_25_AF_LA_c2.0_cond2.bed_length,
                                         diff_313_AF_RA_vs_28_AF_LA_c2.0_cond2.bed_length,
                                         diff_313_AF_RA_vs_31_AF_LA_c2.0_cond2.bed_length,
                                         diff_313_AF_RA_vs_313_AF_LA_c2.0_cond2.bed_length,
                                         diff_313_AF_RA_vs_47_AF_LA_c2.0_cond2.bed_length,
                                         diff_33_SR_RA_vs_2231_AF_LA_c2.0_cond2.bed_length,
                                         diff_33_SR_RA_vs_25_AF_LA_c2.0_cond2.bed_length,
                                         diff_33_SR_RA_vs_28_AF_LA_c2.0_cond2.bed_length,
                                         diff_33_SR_RA_vs_31_AF_LA_c2.0_cond2.bed_length,
                                         diff_33_SR_RA_vs_313_AF_LA_c2.0_cond2.bed_length,
                                         diff_33_SR_RA_vs_47_AF_LA_c2.0_cond2.bed_length,
                                         diff_47_AF_LA_vs_12b_SR_LA_c2.0_cond1.bed_length,
                                         diff_47_AF_LA_vs_20_SR_LA_c2.0_cond1.bed_length,
                                         diff_47_AF_LA_vs_27_SR_LA_c2.0_cond1.bed_length,
                                         diff_47_AF_LA_vs_45_SR_LA_c2.0_cond1.bed_length)


sample_peakset_lengths_boxplot <- ggplot(peaksets_lengths_AF_LA_pairwise, aes(x=sample, y=value)) + geom_boxplot(aes(fill=group),outlier.size = 0.05) + scale_y_continuous("Peak length (bp)") + scale_x_discrete("enriched_peakset_comparison") + theme_classic() + theme(axis.text.x = element_text(angle=45, hjust = 1) )  + scale_fill_manual(values=c('#D85B74')) + theme(axis.text.x = element_blank())                                                                                                                                                                                                                                                                                                                                                                                                                                                                   
sample_peakset_lengths_boxplot
sample_peakset_lengths_boxplot_nooutliers <- ggplot(peaksets_lengths_AF_LA_pairwise, aes(x=sample, y=value)) + geom_boxplot(aes(fill=group),outlier.shape=NA) + coord_cartesian(ylim = c(0, 2500))+ scale_y_continuous("Peak length (bp)") + scale_x_discrete("enriched_peakset_comparison") + theme_classic() + theme(axis.text.x = element_text(angle=45, hjust = 1) )  + scale_fill_manual(values=c('#D85B74'))    
sample_peakset_lengths_boxplot_nooutliers
sample_peakset_lengths_boxplot_nooutliers_noxaxis <- ggplot(peaksets_lengths_AF_LA_pairwise, aes(x=sample, y=value)) + geom_boxplot(aes(fill=group),outlier.shape=NA) + coord_cartesian(ylim = c(300, 2200))+ scale_y_continuous("Peak length (bp)") + scale_x_discrete("enriched_peakset_comparison") + theme_classic() + theme(axis.text.x = element_blank())  + scale_fill_manual(values=c('#D85B74'))    
sample_peakset_lengths_boxplot_nooutliers_noxaxis


#AF-RA

setwd("/Rodriguez_2025_AFepigenome/H3K27ac_enriched_regions_definition/AF_RA_enriched_peaks")
myFiles <- list.files(pattern="*.bed")

diff_18_AF_RA_vs_12b_SR_LA_c2.0_cond1.bed<- as.data.frame(read.table("diff_18_AF_RA_vs_12b_SR_LA_c2.0_cond1.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_18_AF_RA_vs_12b_SR_RA_c2.0_cond1.bed<- as.data.frame(read.table("diff_18_AF_RA_vs_12b_SR_RA_c2.0_cond1.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_18_AF_RA_vs_20_SR_LA_c2.0_cond1.bed<- as.data.frame(read.table("diff_18_AF_RA_vs_20_SR_LA_c2.0_cond1.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_18_AF_RA_vs_20_SR_RA_c2.0_cond1.bed<- as.data.frame(read.table("diff_18_AF_RA_vs_20_SR_RA_c2.0_cond1.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_18_AF_RA_vs_2231_AF_LA_c2.0_cond1.bed<- as.data.frame(read.table("diff_18_AF_RA_vs_2231_AF_LA_c2.0_cond1.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_18_AF_RA_vs_25_AF_LA_c2.0_cond1.bed<- as.data.frame(read.table("diff_18_AF_RA_vs_25_AF_LA_c2.0_cond1.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_18_AF_RA_vs_27_SR_LA_c2.0_cond1.bed<- as.data.frame(read.table("diff_18_AF_RA_vs_27_SR_LA_c2.0_cond1.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_18_AF_RA_vs_27_SR_RA_c2.0_cond1.bed<- as.data.frame(read.table("diff_18_AF_RA_vs_27_SR_RA_c2.0_cond1.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_18_AF_RA_vs_28_AF_LA_c2.0_cond1.bed<- as.data.frame(read.table("diff_18_AF_RA_vs_28_AF_LA_c2.0_cond1.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_18_AF_RA_vs_31_AF_LA_c2.0_cond1.bed<- as.data.frame(read.table("diff_18_AF_RA_vs_31_AF_LA_c2.0_cond1.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_18_AF_RA_vs_313_AF_LA_c2.0_cond1.bed<- as.data.frame(read.table("diff_18_AF_RA_vs_313_AF_LA_c2.0_cond1.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_18_AF_RA_vs_33_SR_RA_c2.0_cond1.bed<- as.data.frame(read.table("diff_18_AF_RA_vs_33_SR_RA_c2.0_cond1.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_18_AF_RA_vs_45_SR_LA_c2.0_cond1.bed<- as.data.frame(read.table("diff_18_AF_RA_vs_45_SR_LA_c2.0_cond1.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_18_AF_RA_vs_47_AF_LA_c2.0_cond1.bed<- as.data.frame(read.table("diff_18_AF_RA_vs_47_AF_LA_c2.0_cond1.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_2231_AF_RA_vs_12b_SR_LA_c2.0_cond1.bed<- as.data.frame(read.table("diff_2231_AF_RA_vs_12b_SR_LA_c2.0_cond1.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_2231_AF_RA_vs_12b_SR_RA_c2.0_cond1.bed<- as.data.frame(read.table("diff_2231_AF_RA_vs_12b_SR_RA_c2.0_cond1.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_2231_AF_RA_vs_20_SR_LA_c2.0_cond1.bed<- as.data.frame(read.table("diff_2231_AF_RA_vs_20_SR_LA_c2.0_cond1.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_2231_AF_RA_vs_20_SR_RA_c2.0_cond1.bed<- as.data.frame(read.table("diff_2231_AF_RA_vs_20_SR_RA_c2.0_cond1.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_2231_AF_RA_vs_2231_AF_LA_c2.0_cond1.bed<- as.data.frame(read.table("diff_2231_AF_RA_vs_2231_AF_LA_c2.0_cond1.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_2231_AF_RA_vs_25_AF_LA_c2.0_cond1.bed<- as.data.frame(read.table("diff_2231_AF_RA_vs_25_AF_LA_c2.0_cond1.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_2231_AF_RA_vs_27_SR_LA_c2.0_cond1.bed<- as.data.frame(read.table("diff_2231_AF_RA_vs_27_SR_LA_c2.0_cond1.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_2231_AF_RA_vs_27_SR_RA_c2.0_cond1.bed<- as.data.frame(read.table("diff_2231_AF_RA_vs_27_SR_RA_c2.0_cond1.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_2231_AF_RA_vs_28_AF_LA_c2.0_cond1.bed<- as.data.frame(read.table("diff_2231_AF_RA_vs_28_AF_LA_c2.0_cond1.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_2231_AF_RA_vs_31_AF_LA_c2.0_cond1.bed<- as.data.frame(read.table("diff_2231_AF_RA_vs_31_AF_LA_c2.0_cond1.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_2231_AF_RA_vs_313_AF_LA_c2.0_cond1.bed<- as.data.frame(read.table("diff_2231_AF_RA_vs_313_AF_LA_c2.0_cond1.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_2231_AF_RA_vs_33_SR_RA_c2.0_cond1.bed<- as.data.frame(read.table("diff_2231_AF_RA_vs_33_SR_RA_c2.0_cond1.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_2231_AF_RA_vs_45_SR_LA_c2.0_cond1.bed<- as.data.frame(read.table("diff_2231_AF_RA_vs_45_SR_LA_c2.0_cond1.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_2231_AF_RA_vs_47_AF_LA_c2.0_cond1.bed<- as.data.frame(read.table("diff_2231_AF_RA_vs_47_AF_LA_c2.0_cond1.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_25_AF_RA_vs_12b_SR_LA_c2.0_cond1.bed<- as.data.frame(read.table("diff_25_AF_RA_vs_12b_SR_LA_c2.0_cond1.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_25_AF_RA_vs_12b_SR_RA_c2.0_cond1.bed<- as.data.frame(read.table("diff_25_AF_RA_vs_12b_SR_RA_c2.0_cond1.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_25_AF_RA_vs_20_SR_LA_c2.0_cond1.bed<- as.data.frame(read.table("diff_25_AF_RA_vs_20_SR_LA_c2.0_cond1.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_25_AF_RA_vs_20_SR_RA_c2.0_cond1.bed<- as.data.frame(read.table("diff_25_AF_RA_vs_20_SR_RA_c2.0_cond1.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_25_AF_RA_vs_2231_AF_LA_c2.0_cond1.bed<- as.data.frame(read.table("diff_25_AF_RA_vs_2231_AF_LA_c2.0_cond1.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_25_AF_RA_vs_25_AF_LA_c2.0_cond1.bed<- as.data.frame(read.table("diff_25_AF_RA_vs_25_AF_LA_c2.0_cond1.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_25_AF_RA_vs_27_SR_LA_c2.0_cond1.bed<- as.data.frame(read.table("diff_25_AF_RA_vs_27_SR_LA_c2.0_cond1.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_25_AF_RA_vs_27_SR_RA_c2.0_cond1.bed<- as.data.frame(read.table("diff_25_AF_RA_vs_27_SR_RA_c2.0_cond1.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_25_AF_RA_vs_28_AF_LA_c2.0_cond1.bed<- as.data.frame(read.table("diff_25_AF_RA_vs_28_AF_LA_c2.0_cond1.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_25_AF_RA_vs_31_AF_LA_c2.0_cond1.bed<- as.data.frame(read.table("diff_25_AF_RA_vs_31_AF_LA_c2.0_cond1.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_25_AF_RA_vs_313_AF_LA_c2.0_cond1.bed<- as.data.frame(read.table("diff_25_AF_RA_vs_313_AF_LA_c2.0_cond1.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_25_AF_RA_vs_33_SR_RA_c2.0_cond1.bed<- as.data.frame(read.table("diff_25_AF_RA_vs_33_SR_RA_c2.0_cond1.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_25_AF_RA_vs_45_SR_LA_c2.0_cond1.bed<- as.data.frame(read.table("diff_25_AF_RA_vs_45_SR_LA_c2.0_cond1.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_25_AF_RA_vs_47_AF_LA_c2.0_cond1.bed<- as.data.frame(read.table("diff_25_AF_RA_vs_47_AF_LA_c2.0_cond1.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_28_AF_RA_vs_12b_SR_LA_c2.0_cond1.bed<- as.data.frame(read.table("diff_28_AF_RA_vs_12b_SR_LA_c2.0_cond1.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_28_AF_RA_vs_12b_SR_RA_c2.0_cond1.bed<- as.data.frame(read.table("diff_28_AF_RA_vs_12b_SR_RA_c2.0_cond1.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_28_AF_RA_vs_20_SR_LA_c2.0_cond1.bed<- as.data.frame(read.table("diff_28_AF_RA_vs_20_SR_LA_c2.0_cond1.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_28_AF_RA_vs_20_SR_RA_c2.0_cond1.bed<- as.data.frame(read.table("diff_28_AF_RA_vs_20_SR_RA_c2.0_cond1.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_28_AF_RA_vs_2231_AF_LA_c2.0_cond1.bed<- as.data.frame(read.table("diff_28_AF_RA_vs_2231_AF_LA_c2.0_cond1.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_28_AF_RA_vs_25_AF_LA_c2.0_cond1.bed<- as.data.frame(read.table("diff_28_AF_RA_vs_25_AF_LA_c2.0_cond1.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_28_AF_RA_vs_27_SR_LA_c2.0_cond1.bed<- as.data.frame(read.table("diff_28_AF_RA_vs_27_SR_LA_c2.0_cond1.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_28_AF_RA_vs_27_SR_RA_c2.0_cond1.bed<- as.data.frame(read.table("diff_28_AF_RA_vs_27_SR_RA_c2.0_cond1.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_28_AF_RA_vs_28_AF_LA_c2.0_cond1.bed<- as.data.frame(read.table("diff_28_AF_RA_vs_28_AF_LA_c2.0_cond1.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_28_AF_RA_vs_31_AF_LA_c2.0_cond1.bed<- as.data.frame(read.table("diff_28_AF_RA_vs_31_AF_LA_c2.0_cond1.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_28_AF_RA_vs_313_AF_LA_c2.0_cond1.bed<- as.data.frame(read.table("diff_28_AF_RA_vs_313_AF_LA_c2.0_cond1.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_28_AF_RA_vs_33_SR_RA_c2.0_cond1.bed<- as.data.frame(read.table("diff_28_AF_RA_vs_33_SR_RA_c2.0_cond1.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_28_AF_RA_vs_45_SR_LA_c2.0_cond1.bed<- as.data.frame(read.table("diff_28_AF_RA_vs_45_SR_LA_c2.0_cond1.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_28_AF_RA_vs_47_AF_LA_c2.0_cond1.bed<- as.data.frame(read.table("diff_28_AF_RA_vs_47_AF_LA_c2.0_cond1.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_31_AF_RA_vs_12b_SR_LA_c2.0_cond1.bed<- as.data.frame(read.table("diff_31_AF_RA_vs_12b_SR_LA_c2.0_cond1.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_31_AF_RA_vs_12b_SR_RA_c2.0_cond1.bed<- as.data.frame(read.table("diff_31_AF_RA_vs_12b_SR_RA_c2.0_cond1.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_31_AF_RA_vs_20_SR_LA_c2.0_cond1.bed<- as.data.frame(read.table("diff_31_AF_RA_vs_20_SR_LA_c2.0_cond1.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_31_AF_RA_vs_20_SR_RA_c2.0_cond1.bed<- as.data.frame(read.table("diff_31_AF_RA_vs_20_SR_RA_c2.0_cond1.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_31_AF_RA_vs_2231_AF_LA_c2.0_cond1.bed<- as.data.frame(read.table("diff_31_AF_RA_vs_2231_AF_LA_c2.0_cond1.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_31_AF_RA_vs_25_AF_LA_c2.0_cond1.bed<- as.data.frame(read.table("diff_31_AF_RA_vs_25_AF_LA_c2.0_cond1.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_31_AF_RA_vs_27_SR_LA_c2.0_cond1.bed<- as.data.frame(read.table("diff_31_AF_RA_vs_27_SR_LA_c2.0_cond1.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_31_AF_RA_vs_27_SR_RA_c2.0_cond1.bed<- as.data.frame(read.table("diff_31_AF_RA_vs_27_SR_RA_c2.0_cond1.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_31_AF_RA_vs_28_AF_LA_c2.0_cond1.bed<- as.data.frame(read.table("diff_31_AF_RA_vs_28_AF_LA_c2.0_cond1.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_31_AF_RA_vs_31_AF_LA_c2.0_cond1.bed<- as.data.frame(read.table("diff_31_AF_RA_vs_31_AF_LA_c2.0_cond1.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_31_AF_RA_vs_313_AF_LA_c2.0_cond1.bed<- as.data.frame(read.table("diff_31_AF_RA_vs_313_AF_LA_c2.0_cond1.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_31_AF_RA_vs_33_SR_RA_c2.0_cond1.bed<- as.data.frame(read.table("diff_31_AF_RA_vs_33_SR_RA_c2.0_cond1.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_31_AF_RA_vs_45_SR_LA_c2.0_cond1.bed<- as.data.frame(read.table("diff_31_AF_RA_vs_45_SR_LA_c2.0_cond1.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_31_AF_RA_vs_47_AF_LA_c2.0_cond1.bed<- as.data.frame(read.table("diff_31_AF_RA_vs_47_AF_LA_c2.0_cond1.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_313_AF_RA_vs_12b_SR_LA_c2.0_cond1.bed<- as.data.frame(read.table("diff_313_AF_RA_vs_12b_SR_LA_c2.0_cond1.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_313_AF_RA_vs_12b_SR_RA_c2.0_cond1.bed<- as.data.frame(read.table("diff_313_AF_RA_vs_12b_SR_RA_c2.0_cond1.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_313_AF_RA_vs_20_SR_LA_c2.0_cond1.bed<- as.data.frame(read.table("diff_313_AF_RA_vs_20_SR_LA_c2.0_cond1.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_313_AF_RA_vs_20_SR_RA_c2.0_cond1.bed<- as.data.frame(read.table("diff_313_AF_RA_vs_20_SR_RA_c2.0_cond1.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_313_AF_RA_vs_2231_AF_LA_c2.0_cond1.bed<- as.data.frame(read.table("diff_313_AF_RA_vs_2231_AF_LA_c2.0_cond1.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_313_AF_RA_vs_25_AF_LA_c2.0_cond1.bed<- as.data.frame(read.table("diff_313_AF_RA_vs_25_AF_LA_c2.0_cond1.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_313_AF_RA_vs_27_SR_LA_c2.0_cond1.bed<- as.data.frame(read.table("diff_313_AF_RA_vs_27_SR_LA_c2.0_cond1.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_313_AF_RA_vs_27_SR_RA_c2.0_cond1.bed<- as.data.frame(read.table("diff_313_AF_RA_vs_27_SR_RA_c2.0_cond1.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_313_AF_RA_vs_28_AF_LA_c2.0_cond1.bed<- as.data.frame(read.table("diff_313_AF_RA_vs_28_AF_LA_c2.0_cond1.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_313_AF_RA_vs_31_AF_LA_c2.0_cond1.bed<- as.data.frame(read.table("diff_313_AF_RA_vs_31_AF_LA_c2.0_cond1.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_313_AF_RA_vs_313_AF_LA_c2.0_cond1.bed<- as.data.frame(read.table("diff_313_AF_RA_vs_313_AF_LA_c2.0_cond1.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_313_AF_RA_vs_33_SR_RA_c2.0_cond1.bed<- as.data.frame(read.table("diff_313_AF_RA_vs_33_SR_RA_c2.0_cond1.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_313_AF_RA_vs_45_SR_LA_c2.0_cond1.bed<- as.data.frame(read.table("diff_313_AF_RA_vs_45_SR_LA_c2.0_cond1.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_313_AF_RA_vs_47_AF_LA_c2.0_cond1.bed<- as.data.frame(read.table("diff_313_AF_RA_vs_47_AF_LA_c2.0_cond1.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))

diff_18_AF_RA_vs_12b_SR_LA_c2.0_cond1.bed_length<- diff_18_AF_RA_vs_12b_SR_LA_c2.0_cond1.bed$V3 - diff_18_AF_RA_vs_12b_SR_LA_c2.0_cond1.bed$V2
diff_18_AF_RA_vs_12b_SR_RA_c2.0_cond1.bed_length<- diff_18_AF_RA_vs_12b_SR_RA_c2.0_cond1.bed$V3 - diff_18_AF_RA_vs_12b_SR_RA_c2.0_cond1.bed$V2
diff_18_AF_RA_vs_20_SR_LA_c2.0_cond1.bed_length<- diff_18_AF_RA_vs_20_SR_LA_c2.0_cond1.bed$V3 - diff_18_AF_RA_vs_20_SR_LA_c2.0_cond1.bed$V2
diff_18_AF_RA_vs_20_SR_RA_c2.0_cond1.bed_length<- diff_18_AF_RA_vs_20_SR_RA_c2.0_cond1.bed$V3 - diff_18_AF_RA_vs_20_SR_RA_c2.0_cond1.bed$V2
diff_18_AF_RA_vs_2231_AF_LA_c2.0_cond1.bed_length<- diff_18_AF_RA_vs_2231_AF_LA_c2.0_cond1.bed$V3 - diff_18_AF_RA_vs_2231_AF_LA_c2.0_cond1.bed$V2
diff_18_AF_RA_vs_25_AF_LA_c2.0_cond1.bed_length<- diff_18_AF_RA_vs_25_AF_LA_c2.0_cond1.bed$V3 - diff_18_AF_RA_vs_25_AF_LA_c2.0_cond1.bed$V2
diff_18_AF_RA_vs_27_SR_LA_c2.0_cond1.bed_length<- diff_18_AF_RA_vs_27_SR_LA_c2.0_cond1.bed$V3 - diff_18_AF_RA_vs_27_SR_LA_c2.0_cond1.bed$V2
diff_18_AF_RA_vs_27_SR_RA_c2.0_cond1.bed_length<- diff_18_AF_RA_vs_27_SR_RA_c2.0_cond1.bed$V3 - diff_18_AF_RA_vs_27_SR_RA_c2.0_cond1.bed$V2
diff_18_AF_RA_vs_28_AF_LA_c2.0_cond1.bed_length<- diff_18_AF_RA_vs_28_AF_LA_c2.0_cond1.bed$V3 - diff_18_AF_RA_vs_28_AF_LA_c2.0_cond1.bed$V2
diff_18_AF_RA_vs_31_AF_LA_c2.0_cond1.bed_length<- diff_18_AF_RA_vs_31_AF_LA_c2.0_cond1.bed$V3 - diff_18_AF_RA_vs_31_AF_LA_c2.0_cond1.bed$V2
diff_18_AF_RA_vs_313_AF_LA_c2.0_cond1.bed_length<- diff_18_AF_RA_vs_313_AF_LA_c2.0_cond1.bed$V3 - diff_18_AF_RA_vs_313_AF_LA_c2.0_cond1.bed$V2
diff_18_AF_RA_vs_33_SR_RA_c2.0_cond1.bed_length<- diff_18_AF_RA_vs_33_SR_RA_c2.0_cond1.bed$V3 - diff_18_AF_RA_vs_33_SR_RA_c2.0_cond1.bed$V2
diff_18_AF_RA_vs_45_SR_LA_c2.0_cond1.bed_length<- diff_18_AF_RA_vs_45_SR_LA_c2.0_cond1.bed$V3 - diff_18_AF_RA_vs_45_SR_LA_c2.0_cond1.bed$V2
diff_18_AF_RA_vs_47_AF_LA_c2.0_cond1.bed_length<- diff_18_AF_RA_vs_47_AF_LA_c2.0_cond1.bed$V3 - diff_18_AF_RA_vs_47_AF_LA_c2.0_cond1.bed$V2
diff_2231_AF_RA_vs_12b_SR_LA_c2.0_cond1.bed_length<- diff_2231_AF_RA_vs_12b_SR_LA_c2.0_cond1.bed$V3 - diff_2231_AF_RA_vs_12b_SR_LA_c2.0_cond1.bed$V2
diff_2231_AF_RA_vs_12b_SR_RA_c2.0_cond1.bed_length<- diff_2231_AF_RA_vs_12b_SR_RA_c2.0_cond1.bed$V3 - diff_2231_AF_RA_vs_12b_SR_RA_c2.0_cond1.bed$V2
diff_2231_AF_RA_vs_20_SR_LA_c2.0_cond1.bed_length<- diff_2231_AF_RA_vs_20_SR_LA_c2.0_cond1.bed$V3 - diff_2231_AF_RA_vs_20_SR_LA_c2.0_cond1.bed$V2
diff_2231_AF_RA_vs_20_SR_RA_c2.0_cond1.bed_length<- diff_2231_AF_RA_vs_20_SR_RA_c2.0_cond1.bed$V3 - diff_2231_AF_RA_vs_20_SR_RA_c2.0_cond1.bed$V2
diff_2231_AF_RA_vs_2231_AF_LA_c2.0_cond1.bed_length<- diff_2231_AF_RA_vs_2231_AF_LA_c2.0_cond1.bed$V3 - diff_2231_AF_RA_vs_2231_AF_LA_c2.0_cond1.bed$V2
diff_2231_AF_RA_vs_25_AF_LA_c2.0_cond1.bed_length<- diff_2231_AF_RA_vs_25_AF_LA_c2.0_cond1.bed$V3 - diff_2231_AF_RA_vs_25_AF_LA_c2.0_cond1.bed$V2
diff_2231_AF_RA_vs_27_SR_LA_c2.0_cond1.bed_length<- diff_2231_AF_RA_vs_27_SR_LA_c2.0_cond1.bed$V3 - diff_2231_AF_RA_vs_27_SR_LA_c2.0_cond1.bed$V2
diff_2231_AF_RA_vs_27_SR_RA_c2.0_cond1.bed_length<- diff_2231_AF_RA_vs_27_SR_RA_c2.0_cond1.bed$V3 - diff_2231_AF_RA_vs_27_SR_RA_c2.0_cond1.bed$V2
diff_2231_AF_RA_vs_28_AF_LA_c2.0_cond1.bed_length<- diff_2231_AF_RA_vs_28_AF_LA_c2.0_cond1.bed$V3 - diff_2231_AF_RA_vs_28_AF_LA_c2.0_cond1.bed$V2
diff_2231_AF_RA_vs_31_AF_LA_c2.0_cond1.bed_length<- diff_2231_AF_RA_vs_31_AF_LA_c2.0_cond1.bed$V3 - diff_2231_AF_RA_vs_31_AF_LA_c2.0_cond1.bed$V2
diff_2231_AF_RA_vs_313_AF_LA_c2.0_cond1.bed_length<- diff_2231_AF_RA_vs_313_AF_LA_c2.0_cond1.bed$V3 - diff_2231_AF_RA_vs_313_AF_LA_c2.0_cond1.bed$V2
diff_2231_AF_RA_vs_33_SR_RA_c2.0_cond1.bed_length<- diff_2231_AF_RA_vs_33_SR_RA_c2.0_cond1.bed$V3 - diff_2231_AF_RA_vs_33_SR_RA_c2.0_cond1.bed$V2
diff_2231_AF_RA_vs_45_SR_LA_c2.0_cond1.bed_length<- diff_2231_AF_RA_vs_45_SR_LA_c2.0_cond1.bed$V3 - diff_2231_AF_RA_vs_45_SR_LA_c2.0_cond1.bed$V2
diff_2231_AF_RA_vs_47_AF_LA_c2.0_cond1.bed_length<- diff_2231_AF_RA_vs_47_AF_LA_c2.0_cond1.bed$V3 - diff_2231_AF_RA_vs_47_AF_LA_c2.0_cond1.bed$V2
diff_25_AF_RA_vs_12b_SR_LA_c2.0_cond1.bed_length<- diff_25_AF_RA_vs_12b_SR_LA_c2.0_cond1.bed$V3 - diff_25_AF_RA_vs_12b_SR_LA_c2.0_cond1.bed$V2
diff_25_AF_RA_vs_12b_SR_RA_c2.0_cond1.bed_length<- diff_25_AF_RA_vs_12b_SR_RA_c2.0_cond1.bed$V3 - diff_25_AF_RA_vs_12b_SR_RA_c2.0_cond1.bed$V2
diff_25_AF_RA_vs_20_SR_LA_c2.0_cond1.bed_length<- diff_25_AF_RA_vs_20_SR_LA_c2.0_cond1.bed$V3 - diff_25_AF_RA_vs_20_SR_LA_c2.0_cond1.bed$V2
diff_25_AF_RA_vs_20_SR_RA_c2.0_cond1.bed_length<- diff_25_AF_RA_vs_20_SR_RA_c2.0_cond1.bed$V3 - diff_25_AF_RA_vs_20_SR_RA_c2.0_cond1.bed$V2
diff_25_AF_RA_vs_2231_AF_LA_c2.0_cond1.bed_length<- diff_25_AF_RA_vs_2231_AF_LA_c2.0_cond1.bed$V3 - diff_25_AF_RA_vs_2231_AF_LA_c2.0_cond1.bed$V2
diff_25_AF_RA_vs_25_AF_LA_c2.0_cond1.bed_length<- diff_25_AF_RA_vs_25_AF_LA_c2.0_cond1.bed$V3 - diff_25_AF_RA_vs_25_AF_LA_c2.0_cond1.bed$V2
diff_25_AF_RA_vs_27_SR_LA_c2.0_cond1.bed_length<- diff_25_AF_RA_vs_27_SR_LA_c2.0_cond1.bed$V3 - diff_25_AF_RA_vs_27_SR_LA_c2.0_cond1.bed$V2
diff_25_AF_RA_vs_27_SR_RA_c2.0_cond1.bed_length<- diff_25_AF_RA_vs_27_SR_RA_c2.0_cond1.bed$V3 - diff_25_AF_RA_vs_27_SR_RA_c2.0_cond1.bed$V2
diff_25_AF_RA_vs_28_AF_LA_c2.0_cond1.bed_length<- diff_25_AF_RA_vs_28_AF_LA_c2.0_cond1.bed$V3 - diff_25_AF_RA_vs_28_AF_LA_c2.0_cond1.bed$V2
diff_25_AF_RA_vs_31_AF_LA_c2.0_cond1.bed_length<- diff_25_AF_RA_vs_31_AF_LA_c2.0_cond1.bed$V3 - diff_25_AF_RA_vs_31_AF_LA_c2.0_cond1.bed$V2
diff_25_AF_RA_vs_313_AF_LA_c2.0_cond1.bed_length<- diff_25_AF_RA_vs_313_AF_LA_c2.0_cond1.bed$V3 - diff_25_AF_RA_vs_313_AF_LA_c2.0_cond1.bed$V2
diff_25_AF_RA_vs_33_SR_RA_c2.0_cond1.bed_length<- diff_25_AF_RA_vs_33_SR_RA_c2.0_cond1.bed$V3 - diff_25_AF_RA_vs_33_SR_RA_c2.0_cond1.bed$V2
diff_25_AF_RA_vs_45_SR_LA_c2.0_cond1.bed_length<- diff_25_AF_RA_vs_45_SR_LA_c2.0_cond1.bed$V3 - diff_25_AF_RA_vs_45_SR_LA_c2.0_cond1.bed$V2
diff_25_AF_RA_vs_47_AF_LA_c2.0_cond1.bed_length<- diff_25_AF_RA_vs_47_AF_LA_c2.0_cond1.bed$V3 - diff_25_AF_RA_vs_47_AF_LA_c2.0_cond1.bed$V2
diff_28_AF_RA_vs_12b_SR_LA_c2.0_cond1.bed_length<- diff_28_AF_RA_vs_12b_SR_LA_c2.0_cond1.bed$V3 - diff_28_AF_RA_vs_12b_SR_LA_c2.0_cond1.bed$V2
diff_28_AF_RA_vs_12b_SR_RA_c2.0_cond1.bed_length<- diff_28_AF_RA_vs_12b_SR_RA_c2.0_cond1.bed$V3 - diff_28_AF_RA_vs_12b_SR_RA_c2.0_cond1.bed$V2
diff_28_AF_RA_vs_20_SR_LA_c2.0_cond1.bed_length<- diff_28_AF_RA_vs_20_SR_LA_c2.0_cond1.bed$V3 - diff_28_AF_RA_vs_20_SR_LA_c2.0_cond1.bed$V2
diff_28_AF_RA_vs_20_SR_RA_c2.0_cond1.bed_length<- diff_28_AF_RA_vs_20_SR_RA_c2.0_cond1.bed$V3 - diff_28_AF_RA_vs_20_SR_RA_c2.0_cond1.bed$V2
diff_28_AF_RA_vs_2231_AF_LA_c2.0_cond1.bed_length<- diff_28_AF_RA_vs_2231_AF_LA_c2.0_cond1.bed$V3 - diff_28_AF_RA_vs_2231_AF_LA_c2.0_cond1.bed$V2
diff_28_AF_RA_vs_25_AF_LA_c2.0_cond1.bed_length<- diff_28_AF_RA_vs_25_AF_LA_c2.0_cond1.bed$V3 - diff_28_AF_RA_vs_25_AF_LA_c2.0_cond1.bed$V2
diff_28_AF_RA_vs_27_SR_LA_c2.0_cond1.bed_length<- diff_28_AF_RA_vs_27_SR_LA_c2.0_cond1.bed$V3 - diff_28_AF_RA_vs_27_SR_LA_c2.0_cond1.bed$V2
diff_28_AF_RA_vs_27_SR_RA_c2.0_cond1.bed_length<- diff_28_AF_RA_vs_27_SR_RA_c2.0_cond1.bed$V3 - diff_28_AF_RA_vs_27_SR_RA_c2.0_cond1.bed$V2
diff_28_AF_RA_vs_28_AF_LA_c2.0_cond1.bed_length<- diff_28_AF_RA_vs_28_AF_LA_c2.0_cond1.bed$V3 - diff_28_AF_RA_vs_28_AF_LA_c2.0_cond1.bed$V2
diff_28_AF_RA_vs_31_AF_LA_c2.0_cond1.bed_length<- diff_28_AF_RA_vs_31_AF_LA_c2.0_cond1.bed$V3 - diff_28_AF_RA_vs_31_AF_LA_c2.0_cond1.bed$V2
diff_28_AF_RA_vs_313_AF_LA_c2.0_cond1.bed_length<- diff_28_AF_RA_vs_313_AF_LA_c2.0_cond1.bed$V3 - diff_28_AF_RA_vs_313_AF_LA_c2.0_cond1.bed$V2
diff_28_AF_RA_vs_33_SR_RA_c2.0_cond1.bed_length<- diff_28_AF_RA_vs_33_SR_RA_c2.0_cond1.bed$V3 - diff_28_AF_RA_vs_33_SR_RA_c2.0_cond1.bed$V2
diff_28_AF_RA_vs_45_SR_LA_c2.0_cond1.bed_length<- diff_28_AF_RA_vs_45_SR_LA_c2.0_cond1.bed$V3 - diff_28_AF_RA_vs_45_SR_LA_c2.0_cond1.bed$V2
diff_28_AF_RA_vs_47_AF_LA_c2.0_cond1.bed_length<- diff_28_AF_RA_vs_47_AF_LA_c2.0_cond1.bed$V3 - diff_28_AF_RA_vs_47_AF_LA_c2.0_cond1.bed$V2
diff_31_AF_RA_vs_12b_SR_LA_c2.0_cond1.bed_length<- diff_31_AF_RA_vs_12b_SR_LA_c2.0_cond1.bed$V3 - diff_31_AF_RA_vs_12b_SR_LA_c2.0_cond1.bed$V2
diff_31_AF_RA_vs_12b_SR_RA_c2.0_cond1.bed_length<- diff_31_AF_RA_vs_12b_SR_RA_c2.0_cond1.bed$V3 - diff_31_AF_RA_vs_12b_SR_RA_c2.0_cond1.bed$V2
diff_31_AF_RA_vs_20_SR_LA_c2.0_cond1.bed_length<- diff_31_AF_RA_vs_20_SR_LA_c2.0_cond1.bed$V3 - diff_31_AF_RA_vs_20_SR_LA_c2.0_cond1.bed$V2
diff_31_AF_RA_vs_20_SR_RA_c2.0_cond1.bed_length<- diff_31_AF_RA_vs_20_SR_RA_c2.0_cond1.bed$V3 - diff_31_AF_RA_vs_20_SR_RA_c2.0_cond1.bed$V2
diff_31_AF_RA_vs_2231_AF_LA_c2.0_cond1.bed_length<- diff_31_AF_RA_vs_2231_AF_LA_c2.0_cond1.bed$V3 - diff_31_AF_RA_vs_2231_AF_LA_c2.0_cond1.bed$V2
diff_31_AF_RA_vs_25_AF_LA_c2.0_cond1.bed_length<- diff_31_AF_RA_vs_25_AF_LA_c2.0_cond1.bed$V3 - diff_31_AF_RA_vs_25_AF_LA_c2.0_cond1.bed$V2
diff_31_AF_RA_vs_27_SR_LA_c2.0_cond1.bed_length<- diff_31_AF_RA_vs_27_SR_LA_c2.0_cond1.bed$V3 - diff_31_AF_RA_vs_27_SR_LA_c2.0_cond1.bed$V2
diff_31_AF_RA_vs_27_SR_RA_c2.0_cond1.bed_length<- diff_31_AF_RA_vs_27_SR_RA_c2.0_cond1.bed$V3 - diff_31_AF_RA_vs_27_SR_RA_c2.0_cond1.bed$V2
diff_31_AF_RA_vs_28_AF_LA_c2.0_cond1.bed_length<- diff_31_AF_RA_vs_28_AF_LA_c2.0_cond1.bed$V3 - diff_31_AF_RA_vs_28_AF_LA_c2.0_cond1.bed$V2
diff_31_AF_RA_vs_31_AF_LA_c2.0_cond1.bed_length<- diff_31_AF_RA_vs_31_AF_LA_c2.0_cond1.bed$V3 - diff_31_AF_RA_vs_31_AF_LA_c2.0_cond1.bed$V2
diff_31_AF_RA_vs_313_AF_LA_c2.0_cond1.bed_length<- diff_31_AF_RA_vs_313_AF_LA_c2.0_cond1.bed$V3 - diff_31_AF_RA_vs_313_AF_LA_c2.0_cond1.bed$V2
diff_31_AF_RA_vs_33_SR_RA_c2.0_cond1.bed_length<- diff_31_AF_RA_vs_33_SR_RA_c2.0_cond1.bed$V3 - diff_31_AF_RA_vs_33_SR_RA_c2.0_cond1.bed$V2
diff_31_AF_RA_vs_45_SR_LA_c2.0_cond1.bed_length<- diff_31_AF_RA_vs_45_SR_LA_c2.0_cond1.bed$V3 - diff_31_AF_RA_vs_45_SR_LA_c2.0_cond1.bed$V2
diff_31_AF_RA_vs_47_AF_LA_c2.0_cond1.bed_length<- diff_31_AF_RA_vs_47_AF_LA_c2.0_cond1.bed$V3 - diff_31_AF_RA_vs_47_AF_LA_c2.0_cond1.bed$V2
diff_313_AF_RA_vs_12b_SR_LA_c2.0_cond1.bed_length<- diff_313_AF_RA_vs_12b_SR_LA_c2.0_cond1.bed$V3 - diff_313_AF_RA_vs_12b_SR_LA_c2.0_cond1.bed$V2
diff_313_AF_RA_vs_12b_SR_RA_c2.0_cond1.bed_length<- diff_313_AF_RA_vs_12b_SR_RA_c2.0_cond1.bed$V3 - diff_313_AF_RA_vs_12b_SR_RA_c2.0_cond1.bed$V2
diff_313_AF_RA_vs_20_SR_LA_c2.0_cond1.bed_length<- diff_313_AF_RA_vs_20_SR_LA_c2.0_cond1.bed$V3 - diff_313_AF_RA_vs_20_SR_LA_c2.0_cond1.bed$V2
diff_313_AF_RA_vs_20_SR_RA_c2.0_cond1.bed_length<- diff_313_AF_RA_vs_20_SR_RA_c2.0_cond1.bed$V3 - diff_313_AF_RA_vs_20_SR_RA_c2.0_cond1.bed$V2
diff_313_AF_RA_vs_2231_AF_LA_c2.0_cond1.bed_length<- diff_313_AF_RA_vs_2231_AF_LA_c2.0_cond1.bed$V3 - diff_313_AF_RA_vs_2231_AF_LA_c2.0_cond1.bed$V2
diff_313_AF_RA_vs_25_AF_LA_c2.0_cond1.bed_length<- diff_313_AF_RA_vs_25_AF_LA_c2.0_cond1.bed$V3 - diff_313_AF_RA_vs_25_AF_LA_c2.0_cond1.bed$V2
diff_313_AF_RA_vs_27_SR_LA_c2.0_cond1.bed_length<- diff_313_AF_RA_vs_27_SR_LA_c2.0_cond1.bed$V3 - diff_313_AF_RA_vs_27_SR_LA_c2.0_cond1.bed$V2
diff_313_AF_RA_vs_27_SR_RA_c2.0_cond1.bed_length<- diff_313_AF_RA_vs_27_SR_RA_c2.0_cond1.bed$V3 - diff_313_AF_RA_vs_27_SR_RA_c2.0_cond1.bed$V2
diff_313_AF_RA_vs_28_AF_LA_c2.0_cond1.bed_length<- diff_313_AF_RA_vs_28_AF_LA_c2.0_cond1.bed$V3 - diff_313_AF_RA_vs_28_AF_LA_c2.0_cond1.bed$V2
diff_313_AF_RA_vs_31_AF_LA_c2.0_cond1.bed_length<- diff_313_AF_RA_vs_31_AF_LA_c2.0_cond1.bed$V3 - diff_313_AF_RA_vs_31_AF_LA_c2.0_cond1.bed$V2
diff_313_AF_RA_vs_313_AF_LA_c2.0_cond1.bed_length<- diff_313_AF_RA_vs_313_AF_LA_c2.0_cond1.bed$V3 - diff_313_AF_RA_vs_313_AF_LA_c2.0_cond1.bed$V2
diff_313_AF_RA_vs_33_SR_RA_c2.0_cond1.bed_length<- diff_313_AF_RA_vs_33_SR_RA_c2.0_cond1.bed$V3 - diff_313_AF_RA_vs_33_SR_RA_c2.0_cond1.bed$V2
diff_313_AF_RA_vs_45_SR_LA_c2.0_cond1.bed_length<- diff_313_AF_RA_vs_45_SR_LA_c2.0_cond1.bed$V3 - diff_313_AF_RA_vs_45_SR_LA_c2.0_cond1.bed$V2
diff_313_AF_RA_vs_47_AF_LA_c2.0_cond1.bed_length<- diff_313_AF_RA_vs_45_SR_LA_c2.0_cond1.bed$V3 - diff_313_AF_RA_vs_45_SR_LA_c2.0_cond1.bed$V2

diff_18_AF_RA_vs_12b_SR_LA_c2.0_cond1.bed_length<-melt(diff_18_AF_RA_vs_12b_SR_LA_c2.0_cond1.bed_length)
diff_18_AF_RA_vs_12b_SR_RA_c2.0_cond1.bed_length<-melt(diff_18_AF_RA_vs_12b_SR_RA_c2.0_cond1.bed_length)
diff_18_AF_RA_vs_20_SR_LA_c2.0_cond1.bed_length<-melt(diff_18_AF_RA_vs_20_SR_LA_c2.0_cond1.bed_length)
diff_18_AF_RA_vs_20_SR_RA_c2.0_cond1.bed_length<-melt(diff_18_AF_RA_vs_20_SR_RA_c2.0_cond1.bed_length)
diff_18_AF_RA_vs_2231_AF_LA_c2.0_cond1.bed_length<-melt(diff_18_AF_RA_vs_2231_AF_LA_c2.0_cond1.bed_length)
diff_18_AF_RA_vs_25_AF_LA_c2.0_cond1.bed_length<-melt(diff_18_AF_RA_vs_25_AF_LA_c2.0_cond1.bed_length)
diff_18_AF_RA_vs_27_SR_LA_c2.0_cond1.bed_length<-melt(diff_18_AF_RA_vs_27_SR_LA_c2.0_cond1.bed_length)
diff_18_AF_RA_vs_27_SR_RA_c2.0_cond1.bed_length<-melt(diff_18_AF_RA_vs_27_SR_RA_c2.0_cond1.bed_length)
diff_18_AF_RA_vs_28_AF_LA_c2.0_cond1.bed_length<-melt(diff_18_AF_RA_vs_28_AF_LA_c2.0_cond1.bed_length)
diff_18_AF_RA_vs_31_AF_LA_c2.0_cond1.bed_length<-melt(diff_18_AF_RA_vs_31_AF_LA_c2.0_cond1.bed_length)
diff_18_AF_RA_vs_313_AF_LA_c2.0_cond1.bed_length<-melt(diff_18_AF_RA_vs_313_AF_LA_c2.0_cond1.bed_length)
diff_18_AF_RA_vs_33_SR_RA_c2.0_cond1.bed_length<-melt(diff_18_AF_RA_vs_33_SR_RA_c2.0_cond1.bed_length)
diff_18_AF_RA_vs_45_SR_LA_c2.0_cond1.bed_length<-melt(diff_18_AF_RA_vs_45_SR_LA_c2.0_cond1.bed_length)
diff_18_AF_RA_vs_47_AF_LA_c2.0_cond1.bed_length<-melt(diff_18_AF_RA_vs_47_AF_LA_c2.0_cond1.bed_length)
diff_2231_AF_RA_vs_12b_SR_LA_c2.0_cond1.bed_length<-melt(diff_2231_AF_RA_vs_12b_SR_LA_c2.0_cond1.bed_length)
diff_2231_AF_RA_vs_12b_SR_RA_c2.0_cond1.bed_length<-melt(diff_2231_AF_RA_vs_12b_SR_RA_c2.0_cond1.bed_length)
diff_2231_AF_RA_vs_20_SR_LA_c2.0_cond1.bed_length<-melt(diff_2231_AF_RA_vs_20_SR_LA_c2.0_cond1.bed_length)
diff_2231_AF_RA_vs_20_SR_RA_c2.0_cond1.bed_length<-melt(diff_2231_AF_RA_vs_20_SR_RA_c2.0_cond1.bed_length)
diff_2231_AF_RA_vs_2231_AF_LA_c2.0_cond1.bed_length<-melt(diff_2231_AF_RA_vs_2231_AF_LA_c2.0_cond1.bed_length)
diff_2231_AF_RA_vs_25_AF_LA_c2.0_cond1.bed_length<-melt(diff_2231_AF_RA_vs_25_AF_LA_c2.0_cond1.bed_length)
diff_2231_AF_RA_vs_27_SR_LA_c2.0_cond1.bed_length<-melt(diff_2231_AF_RA_vs_27_SR_LA_c2.0_cond1.bed_length)
diff_2231_AF_RA_vs_27_SR_RA_c2.0_cond1.bed_length<-melt(diff_2231_AF_RA_vs_27_SR_RA_c2.0_cond1.bed_length)
diff_2231_AF_RA_vs_28_AF_LA_c2.0_cond1.bed_length<-melt(diff_2231_AF_RA_vs_28_AF_LA_c2.0_cond1.bed_length)
diff_2231_AF_RA_vs_31_AF_LA_c2.0_cond1.bed_length<-melt(diff_2231_AF_RA_vs_31_AF_LA_c2.0_cond1.bed_length)
diff_2231_AF_RA_vs_313_AF_LA_c2.0_cond1.bed_length<-melt(diff_2231_AF_RA_vs_313_AF_LA_c2.0_cond1.bed_length)
diff_2231_AF_RA_vs_33_SR_RA_c2.0_cond1.bed_length<-melt(diff_2231_AF_RA_vs_33_SR_RA_c2.0_cond1.bed_length)
diff_2231_AF_RA_vs_45_SR_LA_c2.0_cond1.bed_length<-melt(diff_2231_AF_RA_vs_45_SR_LA_c2.0_cond1.bed_length)
diff_2231_AF_RA_vs_47_AF_LA_c2.0_cond1.bed_length<-melt(diff_2231_AF_RA_vs_47_AF_LA_c2.0_cond1.bed_length)
diff_25_AF_RA_vs_12b_SR_LA_c2.0_cond1.bed_length<-melt(diff_25_AF_RA_vs_12b_SR_LA_c2.0_cond1.bed_length)
diff_25_AF_RA_vs_12b_SR_RA_c2.0_cond1.bed_length<-melt(diff_25_AF_RA_vs_12b_SR_RA_c2.0_cond1.bed_length)
diff_25_AF_RA_vs_20_SR_LA_c2.0_cond1.bed_length<-melt(diff_25_AF_RA_vs_20_SR_LA_c2.0_cond1.bed_length)
diff_25_AF_RA_vs_20_SR_RA_c2.0_cond1.bed_length<-melt(diff_25_AF_RA_vs_20_SR_RA_c2.0_cond1.bed_length)
diff_25_AF_RA_vs_2231_AF_LA_c2.0_cond1.bed_length<-melt(diff_25_AF_RA_vs_2231_AF_LA_c2.0_cond1.bed_length)
diff_25_AF_RA_vs_25_AF_LA_c2.0_cond1.bed_length<-melt(diff_25_AF_RA_vs_25_AF_LA_c2.0_cond1.bed_length)
diff_25_AF_RA_vs_27_SR_LA_c2.0_cond1.bed_length<-melt(diff_25_AF_RA_vs_27_SR_LA_c2.0_cond1.bed_length)
diff_25_AF_RA_vs_27_SR_RA_c2.0_cond1.bed_length<-melt(diff_25_AF_RA_vs_27_SR_RA_c2.0_cond1.bed_length)
diff_25_AF_RA_vs_28_AF_LA_c2.0_cond1.bed_length<-melt(diff_25_AF_RA_vs_28_AF_LA_c2.0_cond1.bed_length)
diff_25_AF_RA_vs_31_AF_LA_c2.0_cond1.bed_length<-melt(diff_25_AF_RA_vs_31_AF_LA_c2.0_cond1.bed_length)
diff_25_AF_RA_vs_313_AF_LA_c2.0_cond1.bed_length<-melt(diff_25_AF_RA_vs_313_AF_LA_c2.0_cond1.bed_length)
diff_25_AF_RA_vs_33_SR_RA_c2.0_cond1.bed_length<-melt(diff_25_AF_RA_vs_33_SR_RA_c2.0_cond1.bed_length)
diff_25_AF_RA_vs_45_SR_LA_c2.0_cond1.bed_length<-melt(diff_25_AF_RA_vs_45_SR_LA_c2.0_cond1.bed_length)
diff_25_AF_RA_vs_47_AF_LA_c2.0_cond1.bed_length<-melt(diff_25_AF_RA_vs_47_AF_LA_c2.0_cond1.bed_length)
diff_28_AF_RA_vs_12b_SR_LA_c2.0_cond1.bed_length<-melt(diff_28_AF_RA_vs_12b_SR_LA_c2.0_cond1.bed_length)
diff_28_AF_RA_vs_12b_SR_RA_c2.0_cond1.bed_length<-melt(diff_28_AF_RA_vs_12b_SR_RA_c2.0_cond1.bed_length)
diff_28_AF_RA_vs_20_SR_LA_c2.0_cond1.bed_length<-melt(diff_28_AF_RA_vs_20_SR_LA_c2.0_cond1.bed_length)
diff_28_AF_RA_vs_20_SR_RA_c2.0_cond1.bed_length<-melt(diff_28_AF_RA_vs_20_SR_RA_c2.0_cond1.bed_length)
diff_28_AF_RA_vs_2231_AF_LA_c2.0_cond1.bed_length<-melt(diff_28_AF_RA_vs_2231_AF_LA_c2.0_cond1.bed_length)
diff_28_AF_RA_vs_25_AF_LA_c2.0_cond1.bed_length<-melt(diff_28_AF_RA_vs_25_AF_LA_c2.0_cond1.bed_length)
diff_28_AF_RA_vs_27_SR_LA_c2.0_cond1.bed_length<-melt(diff_28_AF_RA_vs_27_SR_LA_c2.0_cond1.bed_length)
diff_28_AF_RA_vs_27_SR_RA_c2.0_cond1.bed_length<-melt(diff_28_AF_RA_vs_27_SR_RA_c2.0_cond1.bed_length)
diff_28_AF_RA_vs_28_AF_LA_c2.0_cond1.bed_length<-melt(diff_28_AF_RA_vs_28_AF_LA_c2.0_cond1.bed_length)
diff_28_AF_RA_vs_31_AF_LA_c2.0_cond1.bed_length<-melt(diff_28_AF_RA_vs_31_AF_LA_c2.0_cond1.bed_length)
diff_28_AF_RA_vs_313_AF_LA_c2.0_cond1.bed_length<-melt(diff_28_AF_RA_vs_313_AF_LA_c2.0_cond1.bed_length)
diff_28_AF_RA_vs_33_SR_RA_c2.0_cond1.bed_length<-melt(diff_28_AF_RA_vs_33_SR_RA_c2.0_cond1.bed_length)
diff_28_AF_RA_vs_45_SR_LA_c2.0_cond1.bed_length<-melt(diff_28_AF_RA_vs_45_SR_LA_c2.0_cond1.bed_length)
diff_28_AF_RA_vs_47_AF_LA_c2.0_cond1.bed_length<-melt(diff_28_AF_RA_vs_47_AF_LA_c2.0_cond1.bed_length)
diff_31_AF_RA_vs_12b_SR_LA_c2.0_cond1.bed_length<-melt(diff_31_AF_RA_vs_12b_SR_LA_c2.0_cond1.bed_length)
diff_31_AF_RA_vs_12b_SR_RA_c2.0_cond1.bed_length<-melt(diff_31_AF_RA_vs_12b_SR_RA_c2.0_cond1.bed_length)
diff_31_AF_RA_vs_20_SR_LA_c2.0_cond1.bed_length<-melt(diff_31_AF_RA_vs_20_SR_LA_c2.0_cond1.bed_length)
diff_31_AF_RA_vs_20_SR_RA_c2.0_cond1.bed_length<-melt(diff_31_AF_RA_vs_20_SR_RA_c2.0_cond1.bed_length)
diff_31_AF_RA_vs_2231_AF_LA_c2.0_cond1.bed_length<-melt(diff_31_AF_RA_vs_2231_AF_LA_c2.0_cond1.bed_length)
diff_31_AF_RA_vs_25_AF_LA_c2.0_cond1.bed_length<-melt(diff_31_AF_RA_vs_25_AF_LA_c2.0_cond1.bed_length)
diff_31_AF_RA_vs_27_SR_LA_c2.0_cond1.bed_length<-melt(diff_31_AF_RA_vs_27_SR_LA_c2.0_cond1.bed_length)
diff_31_AF_RA_vs_27_SR_RA_c2.0_cond1.bed_length<-melt(diff_31_AF_RA_vs_27_SR_RA_c2.0_cond1.bed_length)
diff_31_AF_RA_vs_28_AF_LA_c2.0_cond1.bed_length<-melt(diff_31_AF_RA_vs_28_AF_LA_c2.0_cond1.bed_length)
diff_31_AF_RA_vs_31_AF_LA_c2.0_cond1.bed_length<-melt(diff_31_AF_RA_vs_31_AF_LA_c2.0_cond1.bed_length)
diff_31_AF_RA_vs_313_AF_LA_c2.0_cond1.bed_length<-melt(diff_31_AF_RA_vs_313_AF_LA_c2.0_cond1.bed_length)
diff_31_AF_RA_vs_33_SR_RA_c2.0_cond1.bed_length<-melt(diff_31_AF_RA_vs_33_SR_RA_c2.0_cond1.bed_length)
diff_31_AF_RA_vs_45_SR_LA_c2.0_cond1.bed_length<-melt(diff_31_AF_RA_vs_45_SR_LA_c2.0_cond1.bed_length)
diff_31_AF_RA_vs_47_AF_LA_c2.0_cond1.bed_length<-melt(diff_31_AF_RA_vs_47_AF_LA_c2.0_cond1.bed_length)
diff_313_AF_RA_vs_12b_SR_LA_c2.0_cond1.bed_length<-melt(diff_313_AF_RA_vs_12b_SR_LA_c2.0_cond1.bed_length)
diff_313_AF_RA_vs_12b_SR_RA_c2.0_cond1.bed_length<-melt(diff_313_AF_RA_vs_12b_SR_RA_c2.0_cond1.bed_length)
diff_313_AF_RA_vs_20_SR_LA_c2.0_cond1.bed_length<-melt(diff_313_AF_RA_vs_20_SR_LA_c2.0_cond1.bed_length)
diff_313_AF_RA_vs_20_SR_RA_c2.0_cond1.bed_length<-melt(diff_313_AF_RA_vs_20_SR_RA_c2.0_cond1.bed_length)
diff_313_AF_RA_vs_2231_AF_LA_c2.0_cond1.bed_length<-melt(diff_313_AF_RA_vs_2231_AF_LA_c2.0_cond1.bed_length)
diff_313_AF_RA_vs_25_AF_LA_c2.0_cond1.bed_length<-melt(diff_313_AF_RA_vs_25_AF_LA_c2.0_cond1.bed_length)
diff_313_AF_RA_vs_27_SR_LA_c2.0_cond1.bed_length<-melt(diff_313_AF_RA_vs_27_SR_LA_c2.0_cond1.bed_length)
diff_313_AF_RA_vs_27_SR_RA_c2.0_cond1.bed_length<-melt(diff_313_AF_RA_vs_27_SR_RA_c2.0_cond1.bed_length)
diff_313_AF_RA_vs_28_AF_LA_c2.0_cond1.bed_length<-melt(diff_313_AF_RA_vs_28_AF_LA_c2.0_cond1.bed_length)
diff_313_AF_RA_vs_31_AF_LA_c2.0_cond1.bed_length<-melt(diff_313_AF_RA_vs_31_AF_LA_c2.0_cond1.bed_length)
diff_313_AF_RA_vs_313_AF_LA_c2.0_cond1.bed_length<-melt(diff_313_AF_RA_vs_313_AF_LA_c2.0_cond1.bed_length)
diff_313_AF_RA_vs_33_SR_RA_c2.0_cond1.bed_length<-melt(diff_313_AF_RA_vs_33_SR_RA_c2.0_cond1.bed_length)
diff_313_AF_RA_vs_45_SR_LA_c2.0_cond1.bed_length<-melt(diff_313_AF_RA_vs_45_SR_LA_c2.0_cond1.bed_length)
diff_313_AF_RA_vs_47_AF_LA_c2.0_cond1.bed_length<-melt(diff_313_AF_RA_vs_47_AF_LA_c2.0_cond1.bed_length)

diff_18_AF_RA_vs_12b_SR_LA_c2.0_cond1.bed_length$group <- "AF_RA"
diff_18_AF_RA_vs_12b_SR_RA_c2.0_cond1.bed_length$group <- "AF_RA"
diff_18_AF_RA_vs_20_SR_LA_c2.0_cond1.bed_length$group <- "AF_RA"
diff_18_AF_RA_vs_20_SR_RA_c2.0_cond1.bed_length$group <- "AF_RA"
diff_18_AF_RA_vs_2231_AF_LA_c2.0_cond1.bed_length$group <- "AF_RA"
diff_18_AF_RA_vs_25_AF_LA_c2.0_cond1.bed_length$group <- "AF_RA"
diff_18_AF_RA_vs_27_SR_LA_c2.0_cond1.bed_length$group <- "AF_RA"
diff_18_AF_RA_vs_27_SR_RA_c2.0_cond1.bed_length$group <- "AF_RA"
diff_18_AF_RA_vs_28_AF_LA_c2.0_cond1.bed_length$group <- "AF_RA"
diff_18_AF_RA_vs_31_AF_LA_c2.0_cond1.bed_length$group <- "AF_RA"
diff_18_AF_RA_vs_313_AF_LA_c2.0_cond1.bed_length$group <- "AF_RA"
diff_18_AF_RA_vs_33_SR_RA_c2.0_cond1.bed_length$group <- "AF_RA"
diff_18_AF_RA_vs_45_SR_LA_c2.0_cond1.bed_length$group <- "AF_RA"
diff_18_AF_RA_vs_47_AF_LA_c2.0_cond1.bed_length$group <- "AF_RA"
diff_2231_AF_RA_vs_12b_SR_LA_c2.0_cond1.bed_length$group <- "AF_RA"
diff_2231_AF_RA_vs_12b_SR_RA_c2.0_cond1.bed_length$group <- "AF_RA"
diff_2231_AF_RA_vs_20_SR_LA_c2.0_cond1.bed_length$group <- "AF_RA"
diff_2231_AF_RA_vs_20_SR_RA_c2.0_cond1.bed_length$group <- "AF_RA"
diff_2231_AF_RA_vs_2231_AF_LA_c2.0_cond1.bed_length$group <- "AF_RA"
diff_2231_AF_RA_vs_25_AF_LA_c2.0_cond1.bed_length$group <- "AF_RA"
diff_2231_AF_RA_vs_27_SR_LA_c2.0_cond1.bed_length$group <- "AF_RA"
diff_2231_AF_RA_vs_27_SR_RA_c2.0_cond1.bed_length$group <- "AF_RA"
diff_2231_AF_RA_vs_28_AF_LA_c2.0_cond1.bed_length$group <- "AF_RA"
diff_2231_AF_RA_vs_31_AF_LA_c2.0_cond1.bed_length$group <- "AF_RA"
diff_2231_AF_RA_vs_313_AF_LA_c2.0_cond1.bed_length$group <- "AF_RA"
diff_2231_AF_RA_vs_33_SR_RA_c2.0_cond1.bed_length$group <- "AF_RA"
diff_2231_AF_RA_vs_45_SR_LA_c2.0_cond1.bed_length$group <- "AF_RA"
diff_2231_AF_RA_vs_47_AF_LA_c2.0_cond1.bed_length$group <- "AF_RA"
diff_25_AF_RA_vs_12b_SR_LA_c2.0_cond1.bed_length$group <- "AF_RA"
diff_25_AF_RA_vs_12b_SR_RA_c2.0_cond1.bed_length$group <- "AF_RA"
diff_25_AF_RA_vs_20_SR_LA_c2.0_cond1.bed_length$group <- "AF_RA"
diff_25_AF_RA_vs_20_SR_RA_c2.0_cond1.bed_length$group <- "AF_RA"
diff_25_AF_RA_vs_2231_AF_LA_c2.0_cond1.bed_length$group <- "AF_RA"
diff_25_AF_RA_vs_25_AF_LA_c2.0_cond1.bed_length$group <- "AF_RA"
diff_25_AF_RA_vs_27_SR_LA_c2.0_cond1.bed_length$group <- "AF_RA"
diff_25_AF_RA_vs_27_SR_RA_c2.0_cond1.bed_length$group <- "AF_RA"
diff_25_AF_RA_vs_28_AF_LA_c2.0_cond1.bed_length$group <- "AF_RA"
diff_25_AF_RA_vs_31_AF_LA_c2.0_cond1.bed_length$group <- "AF_RA"
diff_25_AF_RA_vs_313_AF_LA_c2.0_cond1.bed_length$group <- "AF_RA"
diff_25_AF_RA_vs_33_SR_RA_c2.0_cond1.bed_length$group <- "AF_RA"
diff_25_AF_RA_vs_45_SR_LA_c2.0_cond1.bed_length$group <- "AF_RA"
diff_25_AF_RA_vs_47_AF_LA_c2.0_cond1.bed_length$group <- "AF_RA"
diff_28_AF_RA_vs_12b_SR_LA_c2.0_cond1.bed_length$group <- "AF_RA"
diff_28_AF_RA_vs_12b_SR_RA_c2.0_cond1.bed_length$group <- "AF_RA"
diff_28_AF_RA_vs_20_SR_LA_c2.0_cond1.bed_length$group <- "AF_RA"
diff_28_AF_RA_vs_20_SR_RA_c2.0_cond1.bed_length$group <- "AF_RA"
diff_28_AF_RA_vs_2231_AF_LA_c2.0_cond1.bed_length$group <- "AF_RA"
diff_28_AF_RA_vs_25_AF_LA_c2.0_cond1.bed_length$group <- "AF_RA"
diff_28_AF_RA_vs_27_SR_LA_c2.0_cond1.bed_length$group <- "AF_RA"
diff_28_AF_RA_vs_27_SR_RA_c2.0_cond1.bed_length$group <- "AF_RA"
diff_28_AF_RA_vs_28_AF_LA_c2.0_cond1.bed_length$group <- "AF_RA"
diff_28_AF_RA_vs_31_AF_LA_c2.0_cond1.bed_length$group <- "AF_RA"
diff_28_AF_RA_vs_313_AF_LA_c2.0_cond1.bed_length$group <- "AF_RA"
diff_28_AF_RA_vs_33_SR_RA_c2.0_cond1.bed_length$group <- "AF_RA"
diff_28_AF_RA_vs_45_SR_LA_c2.0_cond1.bed_length$group <- "AF_RA"
diff_28_AF_RA_vs_47_AF_LA_c2.0_cond1.bed_length$group <- "AF_RA"
diff_31_AF_RA_vs_12b_SR_LA_c2.0_cond1.bed_length$group <- "AF_RA"
diff_31_AF_RA_vs_12b_SR_RA_c2.0_cond1.bed_length$group <- "AF_RA"
diff_31_AF_RA_vs_20_SR_LA_c2.0_cond1.bed_length$group <- "AF_RA"
diff_31_AF_RA_vs_20_SR_RA_c2.0_cond1.bed_length$group <- "AF_RA"
diff_31_AF_RA_vs_2231_AF_LA_c2.0_cond1.bed_length$group <- "AF_RA"
diff_31_AF_RA_vs_25_AF_LA_c2.0_cond1.bed_length$group <- "AF_RA"
diff_31_AF_RA_vs_27_SR_LA_c2.0_cond1.bed_length$group <- "AF_RA"
diff_31_AF_RA_vs_27_SR_RA_c2.0_cond1.bed_length$group <- "AF_RA"
diff_31_AF_RA_vs_28_AF_LA_c2.0_cond1.bed_length$group <- "AF_RA"
diff_31_AF_RA_vs_31_AF_LA_c2.0_cond1.bed_length$group <- "AF_RA"
diff_31_AF_RA_vs_313_AF_LA_c2.0_cond1.bed_length$group <- "AF_RA"
diff_31_AF_RA_vs_33_SR_RA_c2.0_cond1.bed_length$group <- "AF_RA"
diff_31_AF_RA_vs_45_SR_LA_c2.0_cond1.bed_length$group <- "AF_RA"
diff_31_AF_RA_vs_47_AF_LA_c2.0_cond1.bed_length$group <- "AF_RA"
diff_313_AF_RA_vs_12b_SR_LA_c2.0_cond1.bed_length$group <- "AF_RA"
diff_313_AF_RA_vs_12b_SR_RA_c2.0_cond1.bed_length$group <- "AF_RA"
diff_313_AF_RA_vs_20_SR_LA_c2.0_cond1.bed_length$group <- "AF_RA"
diff_313_AF_RA_vs_20_SR_RA_c2.0_cond1.bed_length$group <- "AF_RA"
diff_313_AF_RA_vs_2231_AF_LA_c2.0_cond1.bed_length$group <- "AF_RA"
diff_313_AF_RA_vs_25_AF_LA_c2.0_cond1.bed_length$group <- "AF_RA"
diff_313_AF_RA_vs_27_SR_LA_c2.0_cond1.bed_length$group <- "AF_RA"
diff_313_AF_RA_vs_27_SR_RA_c2.0_cond1.bed_length$group <- "AF_RA"
diff_313_AF_RA_vs_28_AF_LA_c2.0_cond1.bed_length$group <- "AF_RA"
diff_313_AF_RA_vs_31_AF_LA_c2.0_cond1.bed_length$group <- "AF_RA"
diff_313_AF_RA_vs_313_AF_LA_c2.0_cond1.bed_length$group <- "AF_RA"
diff_313_AF_RA_vs_33_SR_RA_c2.0_cond1.bed_length$group <- "AF_RA"
diff_313_AF_RA_vs_45_SR_LA_c2.0_cond1.bed_length$group <- "AF_RA"
diff_313_AF_RA_vs_47_AF_LA_c2.0_cond1.bed_length$group <- "AF_RA"

diff_18_AF_RA_vs_12b_SR_LA_c2.0_cond1.bed_length$sample <- "diff_18_AF_RA_vs_12b_SR_LA_c2.0_cond1.bed"
diff_18_AF_RA_vs_12b_SR_RA_c2.0_cond1.bed_length$sample <- "diff_18_AF_RA_vs_12b_SR_RA_c2.0_cond1.bed"
diff_18_AF_RA_vs_20_SR_LA_c2.0_cond1.bed_length$sample <- "diff_18_AF_RA_vs_20_SR_LA_c2.0_cond1.bed"
diff_18_AF_RA_vs_20_SR_RA_c2.0_cond1.bed_length$sample <- "diff_18_AF_RA_vs_20_SR_RA_c2.0_cond1.bed"
diff_18_AF_RA_vs_2231_AF_LA_c2.0_cond1.bed_length$sample <- "diff_18_AF_RA_vs_2231_AF_LA_c2.0_cond1.bed"
diff_18_AF_RA_vs_25_AF_LA_c2.0_cond1.bed_length$sample <- "diff_18_AF_RA_vs_25_AF_LA_c2.0_cond1.bed"
diff_18_AF_RA_vs_27_SR_LA_c2.0_cond1.bed_length$sample <- "diff_18_AF_RA_vs_27_SR_LA_c2.0_cond1.bed"
diff_18_AF_RA_vs_27_SR_RA_c2.0_cond1.bed_length$sample <- "diff_18_AF_RA_vs_27_SR_RA_c2.0_cond1.bed"
diff_18_AF_RA_vs_28_AF_LA_c2.0_cond1.bed_length$sample <- "diff_18_AF_RA_vs_28_AF_LA_c2.0_cond1.bed"
diff_18_AF_RA_vs_31_AF_LA_c2.0_cond1.bed_length$sample <- "diff_18_AF_RA_vs_31_AF_LA_c2.0_cond1.bed"
diff_18_AF_RA_vs_313_AF_LA_c2.0_cond1.bed_length$sample <- "diff_18_AF_RA_vs_313_AF_LA_c2.0_cond1.bed"
diff_18_AF_RA_vs_33_SR_RA_c2.0_cond1.bed_length$sample <- "diff_18_AF_RA_vs_33_SR_RA_c2.0_cond1.bed"
diff_18_AF_RA_vs_45_SR_LA_c2.0_cond1.bed_length$sample <- "diff_18_AF_RA_vs_45_SR_LA_c2.0_cond1.bed"
diff_18_AF_RA_vs_47_AF_LA_c2.0_cond1.bed_length$sample <- "diff_18_AF_RA_vs_47_AF_LA_c2.0_cond1.bed"
diff_2231_AF_RA_vs_12b_SR_LA_c2.0_cond1.bed_length$sample <- "diff_2231_AF_RA_vs_12b_SR_LA_c2.0_cond1.bed"
diff_2231_AF_RA_vs_12b_SR_RA_c2.0_cond1.bed_length$sample <- "diff_2231_AF_RA_vs_12b_SR_RA_c2.0_cond1.bed"
diff_2231_AF_RA_vs_20_SR_LA_c2.0_cond1.bed_length$sample <- "diff_2231_AF_RA_vs_20_SR_LA_c2.0_cond1.bed"
diff_2231_AF_RA_vs_20_SR_RA_c2.0_cond1.bed_length$sample <- "diff_2231_AF_RA_vs_20_SR_RA_c2.0_cond1.bed"
diff_2231_AF_RA_vs_2231_AF_LA_c2.0_cond1.bed_length$sample <- "diff_2231_AF_RA_vs_2231_AF_LA_c2.0_cond1.bed"
diff_2231_AF_RA_vs_25_AF_LA_c2.0_cond1.bed_length$sample <- "diff_2231_AF_RA_vs_25_AF_LA_c2.0_cond1.bed"
diff_2231_AF_RA_vs_27_SR_LA_c2.0_cond1.bed_length$sample <- "diff_2231_AF_RA_vs_27_SR_LA_c2.0_cond1.bed"
diff_2231_AF_RA_vs_27_SR_RA_c2.0_cond1.bed_length$sample <- "diff_2231_AF_RA_vs_27_SR_RA_c2.0_cond1.bed"
diff_2231_AF_RA_vs_28_AF_LA_c2.0_cond1.bed_length$sample <- "diff_2231_AF_RA_vs_28_AF_LA_c2.0_cond1.bed"
diff_2231_AF_RA_vs_31_AF_LA_c2.0_cond1.bed_length$sample <- "diff_2231_AF_RA_vs_31_AF_LA_c2.0_cond1.bed"
diff_2231_AF_RA_vs_313_AF_LA_c2.0_cond1.bed_length$sample <- "diff_2231_AF_RA_vs_313_AF_LA_c2.0_cond1.bed"
diff_2231_AF_RA_vs_33_SR_RA_c2.0_cond1.bed_length$sample <- "diff_2231_AF_RA_vs_33_SR_RA_c2.0_cond1.bed"
diff_2231_AF_RA_vs_45_SR_LA_c2.0_cond1.bed_length$sample <- "diff_2231_AF_RA_vs_45_SR_LA_c2.0_cond1.bed"
diff_2231_AF_RA_vs_47_AF_LA_c2.0_cond1.bed_length$sample <- "diff_2231_AF_RA_vs_47_AF_LA_c2.0_cond1.bed"
diff_25_AF_RA_vs_12b_SR_LA_c2.0_cond1.bed_length$sample <- "diff_25_AF_RA_vs_12b_SR_LA_c2.0_cond1.bed"
diff_25_AF_RA_vs_12b_SR_RA_c2.0_cond1.bed_length$sample <- "diff_25_AF_RA_vs_12b_SR_RA_c2.0_cond1.bed"
diff_25_AF_RA_vs_20_SR_LA_c2.0_cond1.bed_length$sample <- "diff_25_AF_RA_vs_20_SR_LA_c2.0_cond1.bed"
diff_25_AF_RA_vs_20_SR_RA_c2.0_cond1.bed_length$sample <- "diff_25_AF_RA_vs_20_SR_RA_c2.0_cond1.bed"
diff_25_AF_RA_vs_2231_AF_LA_c2.0_cond1.bed_length$sample <- "diff_25_AF_RA_vs_2231_AF_LA_c2.0_cond1.bed"
diff_25_AF_RA_vs_25_AF_LA_c2.0_cond1.bed_length$sample <- "diff_25_AF_RA_vs_25_AF_LA_c2.0_cond1.bed"
diff_25_AF_RA_vs_27_SR_LA_c2.0_cond1.bed_length$sample <- "diff_25_AF_RA_vs_27_SR_LA_c2.0_cond1.bed"
diff_25_AF_RA_vs_27_SR_RA_c2.0_cond1.bed_length$sample <- "diff_25_AF_RA_vs_27_SR_RA_c2.0_cond1.bed"
diff_25_AF_RA_vs_28_AF_LA_c2.0_cond1.bed_length$sample <- "diff_25_AF_RA_vs_28_AF_LA_c2.0_cond1.bed"
diff_25_AF_RA_vs_31_AF_LA_c2.0_cond1.bed_length$sample <- "diff_25_AF_RA_vs_31_AF_LA_c2.0_cond1.bed"
diff_25_AF_RA_vs_313_AF_LA_c2.0_cond1.bed_length$sample <- "diff_25_AF_RA_vs_313_AF_LA_c2.0_cond1.bed"
diff_25_AF_RA_vs_33_SR_RA_c2.0_cond1.bed_length$sample <- "diff_25_AF_RA_vs_33_SR_RA_c2.0_cond1.bed"
diff_25_AF_RA_vs_45_SR_LA_c2.0_cond1.bed_length$sample <- "diff_25_AF_RA_vs_45_SR_LA_c2.0_cond1.bed"
diff_25_AF_RA_vs_47_AF_LA_c2.0_cond1.bed_length$sample <- "diff_25_AF_RA_vs_47_AF_LA_c2.0_cond1.bed"
diff_28_AF_RA_vs_12b_SR_LA_c2.0_cond1.bed_length$sample <- "diff_28_AF_RA_vs_12b_SR_LA_c2.0_cond1.bed"
diff_28_AF_RA_vs_12b_SR_RA_c2.0_cond1.bed_length$sample <- "diff_28_AF_RA_vs_12b_SR_RA_c2.0_cond1.bed"
diff_28_AF_RA_vs_20_SR_LA_c2.0_cond1.bed_length$sample <- "diff_28_AF_RA_vs_20_SR_LA_c2.0_cond1.bed"
diff_28_AF_RA_vs_20_SR_RA_c2.0_cond1.bed_length$sample <- "diff_28_AF_RA_vs_20_SR_RA_c2.0_cond1.bed"
diff_28_AF_RA_vs_2231_AF_LA_c2.0_cond1.bed_length$sample <- "diff_28_AF_RA_vs_2231_AF_LA_c2.0_cond1.bed"
diff_28_AF_RA_vs_25_AF_LA_c2.0_cond1.bed_length$sample <- "diff_28_AF_RA_vs_25_AF_LA_c2.0_cond1.bed"
diff_28_AF_RA_vs_27_SR_LA_c2.0_cond1.bed_length$sample <- "diff_28_AF_RA_vs_27_SR_LA_c2.0_cond1.bed"
diff_28_AF_RA_vs_27_SR_RA_c2.0_cond1.bed_length$sample <- "diff_28_AF_RA_vs_27_SR_RA_c2.0_cond1.bed"
diff_28_AF_RA_vs_28_AF_LA_c2.0_cond1.bed_length$sample <- "diff_28_AF_RA_vs_28_AF_LA_c2.0_cond1.bed"
diff_28_AF_RA_vs_31_AF_LA_c2.0_cond1.bed_length$sample <- "diff_28_AF_RA_vs_31_AF_LA_c2.0_cond1.bed"
diff_28_AF_RA_vs_313_AF_LA_c2.0_cond1.bed_length$sample <- "diff_28_AF_RA_vs_313_AF_LA_c2.0_cond1.bed"
diff_28_AF_RA_vs_33_SR_RA_c2.0_cond1.bed_length$sample <- "diff_28_AF_RA_vs_33_SR_RA_c2.0_cond1.bed"
diff_28_AF_RA_vs_45_SR_LA_c2.0_cond1.bed_length$sample <- "diff_28_AF_RA_vs_45_SR_LA_c2.0_cond1.bed"
diff_28_AF_RA_vs_47_AF_LA_c2.0_cond1.bed_length$sample <- "diff_28_AF_RA_vs_47_AF_LA_c2.0_cond1.bed"
diff_31_AF_RA_vs_12b_SR_LA_c2.0_cond1.bed_length$sample <- "diff_31_AF_RA_vs_12b_SR_LA_c2.0_cond1.bed"
diff_31_AF_RA_vs_12b_SR_RA_c2.0_cond1.bed_length$sample <- "diff_31_AF_RA_vs_12b_SR_RA_c2.0_cond1.bed"
diff_31_AF_RA_vs_20_SR_LA_c2.0_cond1.bed_length$sample <- "diff_31_AF_RA_vs_20_SR_LA_c2.0_cond1.bed"
diff_31_AF_RA_vs_20_SR_RA_c2.0_cond1.bed_length$sample <- "diff_31_AF_RA_vs_20_SR_RA_c2.0_cond1.bed"
diff_31_AF_RA_vs_2231_AF_LA_c2.0_cond1.bed_length$sample <- "diff_31_AF_RA_vs_2231_AF_LA_c2.0_cond1.bed"
diff_31_AF_RA_vs_25_AF_LA_c2.0_cond1.bed_length$sample <- "diff_31_AF_RA_vs_25_AF_LA_c2.0_cond1.bed"
diff_31_AF_RA_vs_27_SR_LA_c2.0_cond1.bed_length$sample <- "diff_31_AF_RA_vs_27_SR_LA_c2.0_cond1.bed"
diff_31_AF_RA_vs_27_SR_RA_c2.0_cond1.bed_length$sample <- "diff_31_AF_RA_vs_27_SR_RA_c2.0_cond1.bed"
diff_31_AF_RA_vs_28_AF_LA_c2.0_cond1.bed_length$sample <- "diff_31_AF_RA_vs_28_AF_LA_c2.0_cond1.bed"
diff_31_AF_RA_vs_31_AF_LA_c2.0_cond1.bed_length$sample <- "diff_31_AF_RA_vs_31_AF_LA_c2.0_cond1.bed"
diff_31_AF_RA_vs_313_AF_LA_c2.0_cond1.bed_length$sample <- "diff_31_AF_RA_vs_313_AF_LA_c2.0_cond1.bed"
diff_31_AF_RA_vs_33_SR_RA_c2.0_cond1.bed_length$sample <- "diff_31_AF_RA_vs_33_SR_RA_c2.0_cond1.bed"
diff_31_AF_RA_vs_45_SR_LA_c2.0_cond1.bed_length$sample <- "diff_31_AF_RA_vs_45_SR_LA_c2.0_cond1.bed"
diff_31_AF_RA_vs_47_AF_LA_c2.0_cond1.bed_length$sample <- "diff_31_AF_RA_vs_47_AF_LA_c2.0_cond1.bed"
diff_313_AF_RA_vs_12b_SR_LA_c2.0_cond1.bed_length$sample <- "diff_313_AF_RA_vs_12b_SR_LA_c2.0_cond1.bed"
diff_313_AF_RA_vs_12b_SR_RA_c2.0_cond1.bed_length$sample <- "diff_313_AF_RA_vs_12b_SR_RA_c2.0_cond1.bed"
diff_313_AF_RA_vs_20_SR_LA_c2.0_cond1.bed_length$sample <- "diff_313_AF_RA_vs_20_SR_LA_c2.0_cond1.bed"
diff_313_AF_RA_vs_20_SR_RA_c2.0_cond1.bed_length$sample <- "diff_313_AF_RA_vs_20_SR_RA_c2.0_cond1.bed"
diff_313_AF_RA_vs_2231_AF_LA_c2.0_cond1.bed_length$sample <- "diff_313_AF_RA_vs_2231_AF_LA_c2.0_cond1.bed"
diff_313_AF_RA_vs_25_AF_LA_c2.0_cond1.bed_length$sample <- "diff_313_AF_RA_vs_25_AF_LA_c2.0_cond1.bed"
diff_313_AF_RA_vs_27_SR_LA_c2.0_cond1.bed_length$sample <- "diff_313_AF_RA_vs_27_SR_LA_c2.0_cond1.bed"
diff_313_AF_RA_vs_27_SR_RA_c2.0_cond1.bed_length$sample <- "diff_313_AF_RA_vs_27_SR_RA_c2.0_cond1.bed"
diff_313_AF_RA_vs_28_AF_LA_c2.0_cond1.bed_length$sample <- "diff_313_AF_RA_vs_28_AF_LA_c2.0_cond1.bed"
diff_313_AF_RA_vs_31_AF_LA_c2.0_cond1.bed_length$sample <- "diff_313_AF_RA_vs_31_AF_LA_c2.0_cond1.bed"
diff_313_AF_RA_vs_313_AF_LA_c2.0_cond1.bed_length$sample <- "diff_313_AF_RA_vs_313_AF_LA_c2.0_cond1.bed"
diff_313_AF_RA_vs_33_SR_RA_c2.0_cond1.bed_length$sample <- "diff_313_AF_RA_vs_33_SR_RA_c2.0_cond1.bed"
diff_313_AF_RA_vs_45_SR_LA_c2.0_cond1.bed_length$sample <- "diff_313_AF_RA_vs_45_SR_LA_c2.0_cond1.bed"
diff_313_AF_RA_vs_47_AF_LA_c2.0_cond1.bed_length$sample <- "diff_313_AF_RA_vs_47_AF_LA_c2.0_cond1.bed"

peaksets_lengths_AF_RA_pairwise <- rbind(diff_18_AF_RA_vs_12b_SR_LA_c2.0_cond1.bed_length,
                                         diff_18_AF_RA_vs_12b_SR_RA_c2.0_cond1.bed_length,
                                         diff_18_AF_RA_vs_20_SR_LA_c2.0_cond1.bed_length,
                                         diff_18_AF_RA_vs_20_SR_RA_c2.0_cond1.bed_length,
                                         diff_18_AF_RA_vs_2231_AF_LA_c2.0_cond1.bed_length,
                                         diff_18_AF_RA_vs_25_AF_LA_c2.0_cond1.bed_length,
                                         diff_18_AF_RA_vs_27_SR_LA_c2.0_cond1.bed_length,
                                         diff_18_AF_RA_vs_27_SR_RA_c2.0_cond1.bed_length,
                                         diff_18_AF_RA_vs_28_AF_LA_c2.0_cond1.bed_length,
                                         diff_18_AF_RA_vs_31_AF_LA_c2.0_cond1.bed_length,
                                         diff_18_AF_RA_vs_313_AF_LA_c2.0_cond1.bed_length,
                                         diff_18_AF_RA_vs_33_SR_RA_c2.0_cond1.bed_length,
                                         diff_18_AF_RA_vs_45_SR_LA_c2.0_cond1.bed_length,
                                         diff_18_AF_RA_vs_47_AF_LA_c2.0_cond1.bed_length,
                                         diff_2231_AF_RA_vs_12b_SR_LA_c2.0_cond1.bed_length,
                                         diff_2231_AF_RA_vs_12b_SR_RA_c2.0_cond1.bed_length,
                                         diff_2231_AF_RA_vs_20_SR_LA_c2.0_cond1.bed_length,
                                         diff_2231_AF_RA_vs_20_SR_RA_c2.0_cond1.bed_length,
                                         diff_2231_AF_RA_vs_2231_AF_LA_c2.0_cond1.bed_length,
                                         diff_2231_AF_RA_vs_25_AF_LA_c2.0_cond1.bed_length,
                                         diff_2231_AF_RA_vs_27_SR_LA_c2.0_cond1.bed_length,
                                         diff_2231_AF_RA_vs_27_SR_RA_c2.0_cond1.bed_length,
                                         diff_2231_AF_RA_vs_28_AF_LA_c2.0_cond1.bed_length,
                                         diff_2231_AF_RA_vs_31_AF_LA_c2.0_cond1.bed_length,
                                         diff_2231_AF_RA_vs_313_AF_LA_c2.0_cond1.bed_length,
                                         diff_2231_AF_RA_vs_33_SR_RA_c2.0_cond1.bed_length,
                                         diff_2231_AF_RA_vs_45_SR_LA_c2.0_cond1.bed_length,
                                         diff_2231_AF_RA_vs_47_AF_LA_c2.0_cond1.bed_length,
                                         diff_25_AF_RA_vs_12b_SR_LA_c2.0_cond1.bed_length,
                                         diff_25_AF_RA_vs_12b_SR_RA_c2.0_cond1.bed_length,
                                         diff_25_AF_RA_vs_20_SR_LA_c2.0_cond1.bed_length,
                                         diff_25_AF_RA_vs_20_SR_RA_c2.0_cond1.bed_length,
                                         diff_25_AF_RA_vs_2231_AF_LA_c2.0_cond1.bed_length,
                                         diff_25_AF_RA_vs_25_AF_LA_c2.0_cond1.bed_length,
                                         diff_25_AF_RA_vs_27_SR_LA_c2.0_cond1.bed_length,
                                         diff_25_AF_RA_vs_27_SR_RA_c2.0_cond1.bed_length,
                                         diff_25_AF_RA_vs_28_AF_LA_c2.0_cond1.bed_length,
                                         diff_25_AF_RA_vs_31_AF_LA_c2.0_cond1.bed_length,
                                         diff_25_AF_RA_vs_313_AF_LA_c2.0_cond1.bed_length,
                                         diff_25_AF_RA_vs_33_SR_RA_c2.0_cond1.bed_length,
                                         diff_25_AF_RA_vs_45_SR_LA_c2.0_cond1.bed_length,
                                         diff_25_AF_RA_vs_47_AF_LA_c2.0_cond1.bed_length,
                                         diff_28_AF_RA_vs_12b_SR_LA_c2.0_cond1.bed_length,
                                         diff_28_AF_RA_vs_12b_SR_RA_c2.0_cond1.bed_length,
                                         diff_28_AF_RA_vs_20_SR_LA_c2.0_cond1.bed_length,
                                         diff_28_AF_RA_vs_20_SR_RA_c2.0_cond1.bed_length,
                                         diff_28_AF_RA_vs_2231_AF_LA_c2.0_cond1.bed_length,
                                         diff_28_AF_RA_vs_25_AF_LA_c2.0_cond1.bed_length,
                                         diff_28_AF_RA_vs_27_SR_LA_c2.0_cond1.bed_length,
                                         diff_28_AF_RA_vs_27_SR_RA_c2.0_cond1.bed_length,
                                         diff_28_AF_RA_vs_28_AF_LA_c2.0_cond1.bed_length,
                                         diff_28_AF_RA_vs_31_AF_LA_c2.0_cond1.bed_length,
                                         diff_28_AF_RA_vs_313_AF_LA_c2.0_cond1.bed_length,
                                         diff_28_AF_RA_vs_33_SR_RA_c2.0_cond1.bed_length,
                                         diff_28_AF_RA_vs_45_SR_LA_c2.0_cond1.bed_length,
                                         diff_28_AF_RA_vs_47_AF_LA_c2.0_cond1.bed_length,
                                         diff_31_AF_RA_vs_12b_SR_LA_c2.0_cond1.bed_length,
                                         diff_31_AF_RA_vs_12b_SR_RA_c2.0_cond1.bed_length,
                                         diff_31_AF_RA_vs_20_SR_LA_c2.0_cond1.bed_length,
                                         diff_31_AF_RA_vs_20_SR_RA_c2.0_cond1.bed_length,
                                         diff_31_AF_RA_vs_2231_AF_LA_c2.0_cond1.bed_length,
                                         diff_31_AF_RA_vs_25_AF_LA_c2.0_cond1.bed_length,
                                         diff_31_AF_RA_vs_27_SR_LA_c2.0_cond1.bed_length,
                                         diff_31_AF_RA_vs_27_SR_RA_c2.0_cond1.bed_length,
                                         diff_31_AF_RA_vs_28_AF_LA_c2.0_cond1.bed_length,
                                         diff_31_AF_RA_vs_31_AF_LA_c2.0_cond1.bed_length,
                                         diff_31_AF_RA_vs_313_AF_LA_c2.0_cond1.bed_length,
                                         diff_31_AF_RA_vs_33_SR_RA_c2.0_cond1.bed_length,
                                         diff_31_AF_RA_vs_45_SR_LA_c2.0_cond1.bed_length,
                                         diff_31_AF_RA_vs_47_AF_LA_c2.0_cond1.bed_length,
                                         diff_313_AF_RA_vs_12b_SR_LA_c2.0_cond1.bed_length,
                                         diff_313_AF_RA_vs_12b_SR_RA_c2.0_cond1.bed_length,
                                         diff_313_AF_RA_vs_20_SR_LA_c2.0_cond1.bed_length,
                                         diff_313_AF_RA_vs_20_SR_RA_c2.0_cond1.bed_length,
                                         diff_313_AF_RA_vs_2231_AF_LA_c2.0_cond1.bed_length,
                                         diff_313_AF_RA_vs_25_AF_LA_c2.0_cond1.bed_length,
                                         diff_313_AF_RA_vs_27_SR_LA_c2.0_cond1.bed_length,
                                         diff_313_AF_RA_vs_27_SR_RA_c2.0_cond1.bed_length,
                                         diff_313_AF_RA_vs_28_AF_LA_c2.0_cond1.bed_length,
                                         diff_313_AF_RA_vs_31_AF_LA_c2.0_cond1.bed_length,
                                         diff_313_AF_RA_vs_313_AF_LA_c2.0_cond1.bed_length,
                                         diff_313_AF_RA_vs_33_SR_RA_c2.0_cond1.bed_length,
                                         diff_313_AF_RA_vs_45_SR_LA_c2.0_cond1.bed_length,
                                         diff_313_AF_RA_vs_47_AF_LA_c2.0_cond1.bed_length)

sample_peakset_lengths_boxplot <- ggplot(peaksets_lengths_AF_RA_pairwise, aes(x=sample, y=value)) + geom_boxplot(aes(fill=group),outlier.size = 0.05) + scale_y_continuous("Peak length (bp)") + scale_x_discrete("enriched_peakset_comparison") + theme_classic() + theme(axis.text.x = element_text(angle=45, hjust = 1) )  + scale_fill_manual(values=c('#E3D081'))  + theme(axis.text.x = element_blank())                                                                                                                                                                                                                                                                                                                                                                                                                                                                     

sample_peakset_lengths_boxplot
sample_peakset_lengths_boxplot_nooutliers <- ggplot(peaksets_lengths_AF_RA_pairwise, aes(x=sample, y=value)) + geom_boxplot(aes(fill=group),outlier.shape=NA) + coord_cartesian(ylim = c(0, 1750))+ scale_y_continuous("Peak length (bp)") + scale_x_discrete("enriched_peakset_comparison") + theme_classic() + theme(axis.text.x = element_text(angle=45, hjust = 1) ) + scale_fill_manual(values=c('#E3D081'))        
sample_peakset_lengths_boxplot_nooutliers
sample_peakset_lengths_boxplot_nooutliers_noxaxis <- ggplot(peaksets_lengths_AF_RA_pairwise, aes(x=sample, y=value)) + geom_boxplot(aes(fill=group),outlier.shape=NA) + coord_cartesian(ylim = c(300, 1750))+ scale_y_continuous("Peak length (bp)") + scale_x_discrete("enriched_peakset_comparison") + theme_classic() + theme(axis.text.x = element_blank() ) + scale_fill_manual(values=c('#E3D081'))        
sample_peakset_lengths_boxplot_nooutliers_noxaxis


#SR-LA

setwd("/Rodriguez_2025_AFepigenome/H3K27ac_enriched_regions_definition/SR_LA_enriched_peaks")
myFiles <- list.files(pattern="*.bed")


diff_12b_SR_RA_vs_12b_SR_LA_c2.0_cond2.bed<- as.data.frame(read.table("diff_12b_SR_RA_vs_12b_SR_LA_c2.0_cond2.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_12b_SR_RA_vs_20_SR_LA_c2.0_cond2.bed<- as.data.frame(read.table("diff_12b_SR_RA_vs_20_SR_LA_c2.0_cond2.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_12b_SR_RA_vs_27_SR_LA_c2.0_cond2.bed<- as.data.frame(read.table("diff_12b_SR_RA_vs_27_SR_LA_c2.0_cond2.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_12b_SR_RA_vs_45_SR_LA_c2.0_cond2.bed<- as.data.frame(read.table("diff_12b_SR_RA_vs_45_SR_LA_c2.0_cond2.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_18_AF_RA_vs_12b_SR_LA_c2.0_cond2.bed<- as.data.frame(read.table("diff_18_AF_RA_vs_12b_SR_LA_c2.0_cond2.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_18_AF_RA_vs_20_SR_LA_c2.0_cond2.bed<- as.data.frame(read.table("diff_18_AF_RA_vs_20_SR_LA_c2.0_cond2.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_18_AF_RA_vs_27_SR_LA_c2.0_cond2.bed<- as.data.frame(read.table("diff_18_AF_RA_vs_27_SR_LA_c2.0_cond2.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_18_AF_RA_vs_45_SR_LA_c2.0_cond2.bed<- as.data.frame(read.table("diff_18_AF_RA_vs_45_SR_LA_c2.0_cond2.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_20_SR_RA_vs_12b_SR_LA_c2.0_cond2.bed<- as.data.frame(read.table("diff_20_SR_RA_vs_12b_SR_LA_c2.0_cond2.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_20_SR_RA_vs_20_SR_LA_c2.0_cond2.bed<- as.data.frame(read.table("diff_20_SR_RA_vs_20_SR_LA_c2.0_cond2.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_20_SR_RA_vs_27_SR_LA_c2.0_cond2.bed<- as.data.frame(read.table("diff_20_SR_RA_vs_27_SR_LA_c2.0_cond2.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_20_SR_RA_vs_45_SR_LA_c2.0_cond2.bed<- as.data.frame(read.table("diff_20_SR_RA_vs_45_SR_LA_c2.0_cond2.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_2231_AF_LA_vs_12b_SR_LA_c2.0_cond2.bed<- as.data.frame(read.table("diff_2231_AF_LA_vs_12b_SR_LA_c2.0_cond2.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_2231_AF_LA_vs_20_SR_LA_c2.0_cond2.bed<- as.data.frame(read.table("diff_2231_AF_LA_vs_20_SR_LA_c2.0_cond2.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_2231_AF_LA_vs_27_SR_LA_c2.0_cond2.bed<- as.data.frame(read.table("diff_2231_AF_LA_vs_27_SR_LA_c2.0_cond2.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_2231_AF_LA_vs_45_SR_LA_c2.0_cond2.bed<- as.data.frame(read.table("diff_2231_AF_LA_vs_45_SR_LA_c2.0_cond2.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_2231_AF_RA_vs_12b_SR_LA_c2.0_cond2.bed<- as.data.frame(read.table("diff_2231_AF_RA_vs_12b_SR_LA_c2.0_cond2.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_2231_AF_RA_vs_20_SR_LA_c2.0_cond2.bed<- as.data.frame(read.table("diff_2231_AF_RA_vs_20_SR_LA_c2.0_cond2.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_2231_AF_RA_vs_27_SR_LA_c2.0_cond2.bed<- as.data.frame(read.table("diff_2231_AF_RA_vs_27_SR_LA_c2.0_cond2.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_2231_AF_RA_vs_45_SR_LA_c2.0_cond2.bed<- as.data.frame(read.table("diff_2231_AF_RA_vs_45_SR_LA_c2.0_cond2.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_25_AF_LA_vs_12b_SR_LA_c2.0_cond2.bed<- as.data.frame(read.table("diff_25_AF_LA_vs_12b_SR_LA_c2.0_cond2.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_25_AF_LA_vs_20_SR_LA_c2.0_cond2.bed<- as.data.frame(read.table("diff_25_AF_LA_vs_20_SR_LA_c2.0_cond2.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_25_AF_LA_vs_27_SR_LA_c2.0_cond2.bed<- as.data.frame(read.table("diff_25_AF_LA_vs_27_SR_LA_c2.0_cond2.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_25_AF_LA_vs_45_SR_LA_c2.0_cond2.bed<- as.data.frame(read.table("diff_25_AF_LA_vs_45_SR_LA_c2.0_cond2.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_25_AF_RA_vs_12b_SR_LA_c2.0_cond2.bed<- as.data.frame(read.table("diff_25_AF_RA_vs_12b_SR_LA_c2.0_cond2.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_25_AF_RA_vs_20_SR_LA_c2.0_cond2.bed<- as.data.frame(read.table("diff_25_AF_RA_vs_20_SR_LA_c2.0_cond2.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_25_AF_RA_vs_27_SR_LA_c2.0_cond2.bed<- as.data.frame(read.table("diff_25_AF_RA_vs_27_SR_LA_c2.0_cond2.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_25_AF_RA_vs_45_SR_LA_c2.0_cond2.bed<- as.data.frame(read.table("diff_25_AF_RA_vs_45_SR_LA_c2.0_cond2.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_27_SR_RA_vs_12b_SR_LA_c2.0_cond2.bed<- as.data.frame(read.table("diff_27_SR_RA_vs_12b_SR_LA_c2.0_cond2.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_27_SR_RA_vs_20_SR_LA_c2.0_cond2.bed<- as.data.frame(read.table("diff_27_SR_RA_vs_20_SR_LA_c2.0_cond2.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_27_SR_RA_vs_27_SR_LA_c2.0_cond2.bed<- as.data.frame(read.table("diff_27_SR_RA_vs_27_SR_LA_c2.0_cond2.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_27_SR_RA_vs_45_SR_LA_c2.0_cond2.bed<- as.data.frame(read.table("diff_27_SR_RA_vs_45_SR_LA_c2.0_cond2.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_28_AF_LA_vs_12b_SR_LA_c2.0_cond2.bed<- as.data.frame(read.table("diff_28_AF_LA_vs_12b_SR_LA_c2.0_cond2.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_28_AF_LA_vs_20_SR_LA_c2.0_cond2.bed<- as.data.frame(read.table("diff_28_AF_LA_vs_20_SR_LA_c2.0_cond2.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_28_AF_LA_vs_27_SR_LA_c2.0_cond2.bed<- as.data.frame(read.table("diff_28_AF_LA_vs_27_SR_LA_c2.0_cond2.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_28_AF_LA_vs_45_SR_LA_c2.0_cond2.bed<- as.data.frame(read.table("diff_28_AF_LA_vs_45_SR_LA_c2.0_cond2.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_28_AF_RA_vs_12b_SR_LA_c2.0_cond2.bed<- as.data.frame(read.table("diff_28_AF_RA_vs_12b_SR_LA_c2.0_cond2.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_28_AF_RA_vs_20_SR_LA_c2.0_cond2.bed<- as.data.frame(read.table("diff_28_AF_RA_vs_20_SR_LA_c2.0_cond2.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_28_AF_RA_vs_27_SR_LA_c2.0_cond2.bed<- as.data.frame(read.table("diff_28_AF_RA_vs_27_SR_LA_c2.0_cond2.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_28_AF_RA_vs_45_SR_LA_c2.0_cond2.bed<- as.data.frame(read.table("diff_28_AF_RA_vs_45_SR_LA_c2.0_cond2.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_31_AF_LA_vs_12b_SR_LA_c2.0_cond2.bed<- as.data.frame(read.table("diff_31_AF_LA_vs_12b_SR_LA_c2.0_cond2.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_31_AF_LA_vs_20_SR_LA_c2.0_cond2.bed<- as.data.frame(read.table("diff_31_AF_LA_vs_20_SR_LA_c2.0_cond2.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_31_AF_LA_vs_27_SR_LA_c2.0_cond2.bed<- as.data.frame(read.table("diff_31_AF_LA_vs_27_SR_LA_c2.0_cond2.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_31_AF_LA_vs_45_SR_LA_c2.0_cond2.bed<- as.data.frame(read.table("diff_31_AF_LA_vs_45_SR_LA_c2.0_cond2.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_31_AF_RA_vs_12b_SR_LA_c2.0_cond2.bed<- as.data.frame(read.table("diff_31_AF_RA_vs_12b_SR_LA_c2.0_cond2.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_31_AF_RA_vs_20_SR_LA_c2.0_cond2.bed<- as.data.frame(read.table("diff_31_AF_RA_vs_20_SR_LA_c2.0_cond2.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_31_AF_RA_vs_27_SR_LA_c2.0_cond2.bed<- as.data.frame(read.table("diff_31_AF_RA_vs_27_SR_LA_c2.0_cond2.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_31_AF_RA_vs_45_SR_LA_c2.0_cond2.bed<- as.data.frame(read.table("diff_31_AF_RA_vs_45_SR_LA_c2.0_cond2.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_313_AF_LA_vs_12b_SR_LA_c2.0_cond2.bed<- as.data.frame(read.table("diff_313_AF_LA_vs_12b_SR_LA_c2.0_cond2.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_313_AF_LA_vs_20_SR_LA_c2.0_cond2.bed<- as.data.frame(read.table("diff_313_AF_LA_vs_20_SR_LA_c2.0_cond2.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_313_AF_LA_vs_27_SR_LA_c2.0_cond2.bed<- as.data.frame(read.table("diff_313_AF_LA_vs_27_SR_LA_c2.0_cond2.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_313_AF_LA_vs_45_SR_LA_c2.0_cond2.bed<- as.data.frame(read.table("diff_313_AF_LA_vs_45_SR_LA_c2.0_cond2.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_313_AF_RA_vs_12b_SR_LA_c2.0_cond2.bed<- as.data.frame(read.table("diff_313_AF_RA_vs_12b_SR_LA_c2.0_cond2.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_313_AF_RA_vs_20_SR_LA_c2.0_cond2.bed<- as.data.frame(read.table("diff_313_AF_RA_vs_20_SR_LA_c2.0_cond2.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_313_AF_RA_vs_27_SR_LA_c2.0_cond2.bed<- as.data.frame(read.table("diff_313_AF_RA_vs_27_SR_LA_c2.0_cond2.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_313_AF_RA_vs_45_SR_LA_c2.0_cond2.bed<- as.data.frame(read.table("diff_313_AF_RA_vs_45_SR_LA_c2.0_cond2.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_33_SR_RA_vs_12b_SR_LA_c2.0_cond2.bed<- as.data.frame(read.table("diff_33_SR_RA_vs_12b_SR_LA_c2.0_cond2.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_33_SR_RA_vs_20_SR_LA_c2.0_cond2.bed<- as.data.frame(read.table("diff_33_SR_RA_vs_20_SR_LA_c2.0_cond2.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_33_SR_RA_vs_27_SR_LA_c2.0_cond2.bed<- as.data.frame(read.table("diff_33_SR_RA_vs_27_SR_LA_c2.0_cond2.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_33_SR_RA_vs_45_SR_LA_c2.0_cond2.bed<- as.data.frame(read.table("diff_33_SR_RA_vs_45_SR_LA_c2.0_cond2.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_47_AF_LA_vs_12b_SR_LA_c2.0_cond2.bed<- as.data.frame(read.table("diff_47_AF_LA_vs_12b_SR_LA_c2.0_cond2.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_47_AF_LA_vs_20_SR_LA_c2.0_cond2.bed<- as.data.frame(read.table("diff_47_AF_LA_vs_20_SR_LA_c2.0_cond2.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_47_AF_LA_vs_27_SR_LA_c2.0_cond2.bed<- as.data.frame(read.table("diff_47_AF_LA_vs_27_SR_LA_c2.0_cond2.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_47_AF_LA_vs_45_SR_LA_c2.0_cond2.bed<- as.data.frame(read.table("diff_47_AF_LA_vs_45_SR_LA_c2.0_cond2.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))


diff_12b_SR_RA_vs_12b_SR_LA_c2.0_cond2.bed_length<- diff_12b_SR_RA_vs_12b_SR_LA_c2.0_cond2.bed$V3 - diff_12b_SR_RA_vs_12b_SR_LA_c2.0_cond2.bed$V2
diff_12b_SR_RA_vs_20_SR_LA_c2.0_cond2.bed_length<- diff_12b_SR_RA_vs_20_SR_LA_c2.0_cond2.bed$V3 - diff_12b_SR_RA_vs_20_SR_LA_c2.0_cond2.bed$V2
diff_12b_SR_RA_vs_27_SR_LA_c2.0_cond2.bed_length<- diff_12b_SR_RA_vs_27_SR_LA_c2.0_cond2.bed$V3 - diff_12b_SR_RA_vs_27_SR_LA_c2.0_cond2.bed$V2
diff_12b_SR_RA_vs_45_SR_LA_c2.0_cond2.bed_length<- diff_12b_SR_RA_vs_45_SR_LA_c2.0_cond2.bed$V3 - diff_12b_SR_RA_vs_45_SR_LA_c2.0_cond2.bed$V2
diff_18_AF_RA_vs_12b_SR_LA_c2.0_cond2.bed_length<- diff_18_AF_RA_vs_12b_SR_LA_c2.0_cond2.bed$V3 - diff_18_AF_RA_vs_12b_SR_LA_c2.0_cond2.bed$V2
diff_18_AF_RA_vs_20_SR_LA_c2.0_cond2.bed_length<- diff_18_AF_RA_vs_20_SR_LA_c2.0_cond2.bed$V3 - diff_18_AF_RA_vs_20_SR_LA_c2.0_cond2.bed$V2
diff_18_AF_RA_vs_27_SR_LA_c2.0_cond2.bed_length<- diff_18_AF_RA_vs_27_SR_LA_c2.0_cond2.bed$V3 - diff_18_AF_RA_vs_27_SR_LA_c2.0_cond2.bed$V2
diff_18_AF_RA_vs_45_SR_LA_c2.0_cond2.bed_length<- diff_18_AF_RA_vs_45_SR_LA_c2.0_cond2.bed$V3 - diff_18_AF_RA_vs_45_SR_LA_c2.0_cond2.bed$V2
diff_20_SR_RA_vs_12b_SR_LA_c2.0_cond2.bed_length<- diff_20_SR_RA_vs_12b_SR_LA_c2.0_cond2.bed$V3 - diff_20_SR_RA_vs_12b_SR_LA_c2.0_cond2.bed$V2
diff_20_SR_RA_vs_20_SR_LA_c2.0_cond2.bed_length<- diff_20_SR_RA_vs_20_SR_LA_c2.0_cond2.bed$V3 - diff_20_SR_RA_vs_20_SR_LA_c2.0_cond2.bed$V2
diff_20_SR_RA_vs_27_SR_LA_c2.0_cond2.bed_length<- diff_20_SR_RA_vs_27_SR_LA_c2.0_cond2.bed$V3 - diff_20_SR_RA_vs_27_SR_LA_c2.0_cond2.bed$V2
diff_20_SR_RA_vs_45_SR_LA_c2.0_cond2.bed_length<- diff_20_SR_RA_vs_45_SR_LA_c2.0_cond2.bed$V3 - diff_20_SR_RA_vs_45_SR_LA_c2.0_cond2.bed$V2
diff_2231_AF_LA_vs_12b_SR_LA_c2.0_cond2.bed_length<- diff_2231_AF_LA_vs_12b_SR_LA_c2.0_cond2.bed$V3 - diff_2231_AF_LA_vs_12b_SR_LA_c2.0_cond2.bed$V2
diff_2231_AF_LA_vs_20_SR_LA_c2.0_cond2.bed_length<- diff_2231_AF_LA_vs_20_SR_LA_c2.0_cond2.bed$V3 - diff_2231_AF_LA_vs_20_SR_LA_c2.0_cond2.bed$V2
diff_2231_AF_LA_vs_27_SR_LA_c2.0_cond2.bed_length<- diff_2231_AF_LA_vs_27_SR_LA_c2.0_cond2.bed$V3 - diff_2231_AF_LA_vs_27_SR_LA_c2.0_cond2.bed$V2
diff_2231_AF_LA_vs_45_SR_LA_c2.0_cond2.bed_length<- diff_2231_AF_LA_vs_45_SR_LA_c2.0_cond2.bed$V3 - diff_2231_AF_LA_vs_45_SR_LA_c2.0_cond2.bed$V2
diff_2231_AF_RA_vs_12b_SR_LA_c2.0_cond2.bed_length<- diff_2231_AF_RA_vs_12b_SR_LA_c2.0_cond2.bed$V3 - diff_2231_AF_RA_vs_12b_SR_LA_c2.0_cond2.bed$V2
diff_2231_AF_RA_vs_20_SR_LA_c2.0_cond2.bed_length<- diff_2231_AF_RA_vs_20_SR_LA_c2.0_cond2.bed$V3 - diff_2231_AF_RA_vs_20_SR_LA_c2.0_cond2.bed$V2
diff_2231_AF_RA_vs_27_SR_LA_c2.0_cond2.bed_length<- diff_2231_AF_RA_vs_27_SR_LA_c2.0_cond2.bed$V3 - diff_2231_AF_RA_vs_27_SR_LA_c2.0_cond2.bed$V2
diff_2231_AF_RA_vs_45_SR_LA_c2.0_cond2.bed_length<- diff_2231_AF_RA_vs_45_SR_LA_c2.0_cond2.bed$V3 - diff_2231_AF_RA_vs_45_SR_LA_c2.0_cond2.bed$V2
diff_25_AF_LA_vs_12b_SR_LA_c2.0_cond2.bed_length<- diff_25_AF_LA_vs_12b_SR_LA_c2.0_cond2.bed$V3 - diff_25_AF_LA_vs_12b_SR_LA_c2.0_cond2.bed$V2
diff_25_AF_LA_vs_20_SR_LA_c2.0_cond2.bed_length<- diff_25_AF_LA_vs_20_SR_LA_c2.0_cond2.bed$V3 - diff_25_AF_LA_vs_20_SR_LA_c2.0_cond2.bed$V2
diff_25_AF_LA_vs_27_SR_LA_c2.0_cond2.bed_length<- diff_25_AF_LA_vs_27_SR_LA_c2.0_cond2.bed$V3 - diff_25_AF_LA_vs_27_SR_LA_c2.0_cond2.bed$V2
diff_25_AF_LA_vs_45_SR_LA_c2.0_cond2.bed_length<- diff_25_AF_LA_vs_45_SR_LA_c2.0_cond2.bed$V3 - diff_25_AF_LA_vs_45_SR_LA_c2.0_cond2.bed$V2
diff_25_AF_RA_vs_12b_SR_LA_c2.0_cond2.bed_length<- diff_25_AF_RA_vs_12b_SR_LA_c2.0_cond2.bed$V3 - diff_25_AF_RA_vs_12b_SR_LA_c2.0_cond2.bed$V2
diff_25_AF_RA_vs_20_SR_LA_c2.0_cond2.bed_length<- diff_25_AF_RA_vs_20_SR_LA_c2.0_cond2.bed$V3 - diff_25_AF_RA_vs_20_SR_LA_c2.0_cond2.bed$V2
diff_25_AF_RA_vs_27_SR_LA_c2.0_cond2.bed_length<- diff_25_AF_RA_vs_27_SR_LA_c2.0_cond2.bed$V3 - diff_25_AF_RA_vs_27_SR_LA_c2.0_cond2.bed$V2
diff_25_AF_RA_vs_45_SR_LA_c2.0_cond2.bed_length<- diff_25_AF_RA_vs_45_SR_LA_c2.0_cond2.bed$V3 - diff_25_AF_RA_vs_45_SR_LA_c2.0_cond2.bed$V2
diff_27_SR_RA_vs_12b_SR_LA_c2.0_cond2.bed_length<- diff_27_SR_RA_vs_12b_SR_LA_c2.0_cond2.bed$V3 - diff_27_SR_RA_vs_12b_SR_LA_c2.0_cond2.bed$V2
diff_27_SR_RA_vs_20_SR_LA_c2.0_cond2.bed_length<- diff_27_SR_RA_vs_20_SR_LA_c2.0_cond2.bed$V3 - diff_27_SR_RA_vs_20_SR_LA_c2.0_cond2.bed$V2
diff_27_SR_RA_vs_27_SR_LA_c2.0_cond2.bed_length<- diff_27_SR_RA_vs_27_SR_LA_c2.0_cond2.bed$V3 - diff_27_SR_RA_vs_27_SR_LA_c2.0_cond2.bed$V2
diff_27_SR_RA_vs_45_SR_LA_c2.0_cond2.bed_length<- diff_27_SR_RA_vs_45_SR_LA_c2.0_cond2.bed$V3 - diff_27_SR_RA_vs_45_SR_LA_c2.0_cond2.bed$V2
diff_28_AF_LA_vs_12b_SR_LA_c2.0_cond2.bed_length<- diff_28_AF_LA_vs_12b_SR_LA_c2.0_cond2.bed$V3 - diff_28_AF_LA_vs_12b_SR_LA_c2.0_cond2.bed$V2
diff_28_AF_LA_vs_20_SR_LA_c2.0_cond2.bed_length<- diff_28_AF_LA_vs_20_SR_LA_c2.0_cond2.bed$V3 - diff_28_AF_LA_vs_20_SR_LA_c2.0_cond2.bed$V2
diff_28_AF_LA_vs_27_SR_LA_c2.0_cond2.bed_length<- diff_28_AF_LA_vs_27_SR_LA_c2.0_cond2.bed$V3 - diff_28_AF_LA_vs_27_SR_LA_c2.0_cond2.bed$V2
diff_28_AF_LA_vs_45_SR_LA_c2.0_cond2.bed_length<- diff_28_AF_LA_vs_45_SR_LA_c2.0_cond2.bed$V3 - diff_28_AF_LA_vs_45_SR_LA_c2.0_cond2.bed$V2
diff_28_AF_RA_vs_12b_SR_LA_c2.0_cond2.bed_length<- diff_28_AF_RA_vs_12b_SR_LA_c2.0_cond2.bed$V3 - diff_28_AF_RA_vs_12b_SR_LA_c2.0_cond2.bed$V2
diff_28_AF_RA_vs_20_SR_LA_c2.0_cond2.bed_length<- diff_28_AF_RA_vs_20_SR_LA_c2.0_cond2.bed$V3 - diff_28_AF_RA_vs_20_SR_LA_c2.0_cond2.bed$V2
diff_28_AF_RA_vs_27_SR_LA_c2.0_cond2.bed_length<- diff_28_AF_RA_vs_27_SR_LA_c2.0_cond2.bed$V3 - diff_28_AF_RA_vs_27_SR_LA_c2.0_cond2.bed$V2
diff_28_AF_RA_vs_45_SR_LA_c2.0_cond2.bed_length<- diff_28_AF_RA_vs_45_SR_LA_c2.0_cond2.bed$V3 - diff_28_AF_RA_vs_45_SR_LA_c2.0_cond2.bed$V2
diff_31_AF_LA_vs_12b_SR_LA_c2.0_cond2.bed_length<- diff_31_AF_LA_vs_12b_SR_LA_c2.0_cond2.bed$V3 - diff_31_AF_LA_vs_12b_SR_LA_c2.0_cond2.bed$V2
diff_31_AF_LA_vs_20_SR_LA_c2.0_cond2.bed_length<- diff_31_AF_LA_vs_20_SR_LA_c2.0_cond2.bed$V3 - diff_31_AF_LA_vs_20_SR_LA_c2.0_cond2.bed$V2
diff_31_AF_LA_vs_27_SR_LA_c2.0_cond2.bed_length<- diff_31_AF_LA_vs_27_SR_LA_c2.0_cond2.bed$V3 - diff_31_AF_LA_vs_27_SR_LA_c2.0_cond2.bed$V2
diff_31_AF_LA_vs_45_SR_LA_c2.0_cond2.bed_length<- diff_31_AF_LA_vs_45_SR_LA_c2.0_cond2.bed$V3 - diff_31_AF_LA_vs_45_SR_LA_c2.0_cond2.bed$V2
diff_31_AF_RA_vs_12b_SR_LA_c2.0_cond2.bed_length<- diff_31_AF_RA_vs_12b_SR_LA_c2.0_cond2.bed$V3 - diff_31_AF_RA_vs_12b_SR_LA_c2.0_cond2.bed$V2
diff_31_AF_RA_vs_20_SR_LA_c2.0_cond2.bed_length<- diff_31_AF_RA_vs_20_SR_LA_c2.0_cond2.bed$V3 - diff_31_AF_RA_vs_20_SR_LA_c2.0_cond2.bed$V2
diff_31_AF_RA_vs_27_SR_LA_c2.0_cond2.bed_length<- diff_31_AF_RA_vs_27_SR_LA_c2.0_cond2.bed$V3 - diff_31_AF_RA_vs_27_SR_LA_c2.0_cond2.bed$V2
diff_31_AF_RA_vs_45_SR_LA_c2.0_cond2.bed_length<- diff_31_AF_RA_vs_45_SR_LA_c2.0_cond2.bed$V3 - diff_31_AF_RA_vs_45_SR_LA_c2.0_cond2.bed$V2
diff_313_AF_LA_vs_12b_SR_LA_c2.0_cond2.bed_length<- diff_313_AF_LA_vs_12b_SR_LA_c2.0_cond2.bed$V3 - diff_313_AF_LA_vs_12b_SR_LA_c2.0_cond2.bed$V2
diff_313_AF_LA_vs_20_SR_LA_c2.0_cond2.bed_length<- diff_313_AF_LA_vs_20_SR_LA_c2.0_cond2.bed$V3 - diff_313_AF_LA_vs_20_SR_LA_c2.0_cond2.bed$V2
diff_313_AF_LA_vs_27_SR_LA_c2.0_cond2.bed_length<- diff_313_AF_LA_vs_27_SR_LA_c2.0_cond2.bed$V3 - diff_313_AF_LA_vs_27_SR_LA_c2.0_cond2.bed$V2
diff_313_AF_LA_vs_45_SR_LA_c2.0_cond2.bed_length<- diff_313_AF_LA_vs_45_SR_LA_c2.0_cond2.bed$V3 - diff_313_AF_LA_vs_45_SR_LA_c2.0_cond2.bed$V2
diff_313_AF_RA_vs_12b_SR_LA_c2.0_cond2.bed_length<- diff_313_AF_RA_vs_12b_SR_LA_c2.0_cond2.bed$V3 - diff_313_AF_RA_vs_12b_SR_LA_c2.0_cond2.bed$V2
diff_313_AF_RA_vs_20_SR_LA_c2.0_cond2.bed_length<- diff_313_AF_RA_vs_20_SR_LA_c2.0_cond2.bed$V3 - diff_313_AF_RA_vs_20_SR_LA_c2.0_cond2.bed$V2
diff_313_AF_RA_vs_27_SR_LA_c2.0_cond2.bed_length<- diff_313_AF_RA_vs_20_SR_LA_c2.0_cond2.bed$V3 - diff_313_AF_RA_vs_20_SR_LA_c2.0_cond2.bed$V2
diff_313_AF_RA_vs_45_SR_LA_c2.0_cond2.bed_length<- diff_313_AF_RA_vs_45_SR_LA_c2.0_cond2.bed$V3 - diff_313_AF_RA_vs_45_SR_LA_c2.0_cond2.bed$V2
diff_33_SR_RA_vs_12b_SR_LA_c2.0_cond2.bed_length<- diff_33_SR_RA_vs_12b_SR_LA_c2.0_cond2.bed$V3 - diff_33_SR_RA_vs_12b_SR_LA_c2.0_cond2.bed$V2
diff_33_SR_RA_vs_20_SR_LA_c2.0_cond2.bed_length<- diff_33_SR_RA_vs_20_SR_LA_c2.0_cond2.bed$V3 - diff_33_SR_RA_vs_20_SR_LA_c2.0_cond2.bed$V2
diff_33_SR_RA_vs_27_SR_LA_c2.0_cond2.bed_length<- diff_33_SR_RA_vs_27_SR_LA_c2.0_cond2.bed$V3 - diff_33_SR_RA_vs_27_SR_LA_c2.0_cond2.bed$V2
diff_33_SR_RA_vs_45_SR_LA_c2.0_cond2.bed_length<- diff_33_SR_RA_vs_45_SR_LA_c2.0_cond2.bed$V3 - diff_33_SR_RA_vs_45_SR_LA_c2.0_cond2.bed$V2
diff_47_AF_LA_vs_12b_SR_LA_c2.0_cond2.bed_length<- diff_47_AF_LA_vs_12b_SR_LA_c2.0_cond2.bed$V3 - diff_47_AF_LA_vs_12b_SR_LA_c2.0_cond2.bed$V2
diff_47_AF_LA_vs_20_SR_LA_c2.0_cond2.bed_length<- diff_47_AF_LA_vs_20_SR_LA_c2.0_cond2.bed$V3 - diff_47_AF_LA_vs_20_SR_LA_c2.0_cond2.bed$V2
diff_47_AF_LA_vs_27_SR_LA_c2.0_cond2.bed_length<- diff_47_AF_LA_vs_27_SR_LA_c2.0_cond2.bed$V3 - diff_47_AF_LA_vs_27_SR_LA_c2.0_cond2.bed$V2
diff_47_AF_LA_vs_45_SR_LA_c2.0_cond2.bed_length<- diff_47_AF_LA_vs_45_SR_LA_c2.0_cond2.bed$V3 - diff_47_AF_LA_vs_45_SR_LA_c2.0_cond2.bed$V2



diff_12b_SR_RA_vs_12b_SR_LA_c2.0_cond2.bed_length<- melt(diff_12b_SR_RA_vs_12b_SR_LA_c2.0_cond2.bed_length)
diff_12b_SR_RA_vs_20_SR_LA_c2.0_cond2.bed_length<- melt(diff_12b_SR_RA_vs_20_SR_LA_c2.0_cond2.bed_length)
diff_12b_SR_RA_vs_27_SR_LA_c2.0_cond2.bed_length<- melt(diff_12b_SR_RA_vs_27_SR_LA_c2.0_cond2.bed_length)
diff_12b_SR_RA_vs_45_SR_LA_c2.0_cond2.bed_length<- melt(diff_12b_SR_RA_vs_45_SR_LA_c2.0_cond2.bed_length)
diff_18_AF_RA_vs_12b_SR_LA_c2.0_cond2.bed_length<- melt(diff_18_AF_RA_vs_12b_SR_LA_c2.0_cond2.bed_length)
diff_18_AF_RA_vs_20_SR_LA_c2.0_cond2.bed_length<- melt(diff_18_AF_RA_vs_20_SR_LA_c2.0_cond2.bed_length)
diff_18_AF_RA_vs_27_SR_LA_c2.0_cond2.bed_length<- melt(diff_18_AF_RA_vs_27_SR_LA_c2.0_cond2.bed_length)
diff_18_AF_RA_vs_45_SR_LA_c2.0_cond2.bed_length<- melt(diff_18_AF_RA_vs_45_SR_LA_c2.0_cond2.bed_length)
diff_20_SR_RA_vs_12b_SR_LA_c2.0_cond2.bed_length<- melt(diff_20_SR_RA_vs_12b_SR_LA_c2.0_cond2.bed_length)
diff_20_SR_RA_vs_20_SR_LA_c2.0_cond2.bed_length<- melt(diff_20_SR_RA_vs_20_SR_LA_c2.0_cond2.bed_length)
diff_20_SR_RA_vs_27_SR_LA_c2.0_cond2.bed_length<- melt(diff_20_SR_RA_vs_27_SR_LA_c2.0_cond2.bed_length)
diff_20_SR_RA_vs_45_SR_LA_c2.0_cond2.bed_length<- melt(diff_20_SR_RA_vs_45_SR_LA_c2.0_cond2.bed_length)
diff_2231_AF_LA_vs_12b_SR_LA_c2.0_cond2.bed_length<- melt(diff_2231_AF_LA_vs_12b_SR_LA_c2.0_cond2.bed_length)
diff_2231_AF_LA_vs_20_SR_LA_c2.0_cond2.bed_length<- melt(diff_2231_AF_LA_vs_20_SR_LA_c2.0_cond2.bed_length)
diff_2231_AF_LA_vs_27_SR_LA_c2.0_cond2.bed_length<- melt(diff_2231_AF_LA_vs_27_SR_LA_c2.0_cond2.bed_length)
diff_2231_AF_LA_vs_45_SR_LA_c2.0_cond2.bed_length<- melt(diff_2231_AF_LA_vs_45_SR_LA_c2.0_cond2.bed_length)
diff_2231_AF_RA_vs_12b_SR_LA_c2.0_cond2.bed_length<- melt(diff_2231_AF_RA_vs_12b_SR_LA_c2.0_cond2.bed_length)
diff_2231_AF_RA_vs_20_SR_LA_c2.0_cond2.bed_length<- melt(diff_2231_AF_RA_vs_20_SR_LA_c2.0_cond2.bed_length)
diff_2231_AF_RA_vs_27_SR_LA_c2.0_cond2.bed_length<- melt(diff_2231_AF_RA_vs_27_SR_LA_c2.0_cond2.bed_length)
diff_2231_AF_RA_vs_45_SR_LA_c2.0_cond2.bed_length<- melt(diff_2231_AF_RA_vs_45_SR_LA_c2.0_cond2.bed_length)
diff_25_AF_LA_vs_12b_SR_LA_c2.0_cond2.bed_length<- melt(diff_25_AF_LA_vs_12b_SR_LA_c2.0_cond2.bed_length)
diff_25_AF_LA_vs_20_SR_LA_c2.0_cond2.bed_length<- melt(diff_25_AF_LA_vs_20_SR_LA_c2.0_cond2.bed_length)
diff_25_AF_LA_vs_27_SR_LA_c2.0_cond2.bed_length<- melt(diff_25_AF_LA_vs_27_SR_LA_c2.0_cond2.bed_length)
diff_25_AF_LA_vs_45_SR_LA_c2.0_cond2.bed_length<- melt(diff_25_AF_LA_vs_45_SR_LA_c2.0_cond2.bed_length)
diff_25_AF_RA_vs_12b_SR_LA_c2.0_cond2.bed_length<- melt(diff_25_AF_RA_vs_12b_SR_LA_c2.0_cond2.bed_length)
diff_25_AF_RA_vs_20_SR_LA_c2.0_cond2.bed_length<- melt(diff_25_AF_RA_vs_20_SR_LA_c2.0_cond2.bed_length)
diff_25_AF_RA_vs_27_SR_LA_c2.0_cond2.bed_length<- melt(diff_25_AF_RA_vs_27_SR_LA_c2.0_cond2.bed_length)
diff_25_AF_RA_vs_45_SR_LA_c2.0_cond2.bed_length<- melt(diff_25_AF_RA_vs_45_SR_LA_c2.0_cond2.bed_length)
diff_27_SR_RA_vs_12b_SR_LA_c2.0_cond2.bed_length<- melt(diff_27_SR_RA_vs_12b_SR_LA_c2.0_cond2.bed_length)
diff_27_SR_RA_vs_20_SR_LA_c2.0_cond2.bed_length<- melt(diff_27_SR_RA_vs_20_SR_LA_c2.0_cond2.bed_length)
diff_27_SR_RA_vs_27_SR_LA_c2.0_cond2.bed_length<- melt(diff_27_SR_RA_vs_27_SR_LA_c2.0_cond2.bed_length)
diff_27_SR_RA_vs_45_SR_LA_c2.0_cond2.bed_length<- melt(diff_27_SR_RA_vs_45_SR_LA_c2.0_cond2.bed_length)
diff_28_AF_LA_vs_12b_SR_LA_c2.0_cond2.bed_length<- melt(diff_28_AF_LA_vs_12b_SR_LA_c2.0_cond2.bed_length)
diff_28_AF_LA_vs_20_SR_LA_c2.0_cond2.bed_length<- melt(diff_28_AF_LA_vs_20_SR_LA_c2.0_cond2.bed_length)
diff_28_AF_LA_vs_27_SR_LA_c2.0_cond2.bed_length<- melt(diff_28_AF_LA_vs_27_SR_LA_c2.0_cond2.bed_length)
diff_28_AF_LA_vs_45_SR_LA_c2.0_cond2.bed_length<- melt(diff_28_AF_LA_vs_45_SR_LA_c2.0_cond2.bed_length)
diff_28_AF_RA_vs_12b_SR_LA_c2.0_cond2.bed_length<- melt(diff_28_AF_RA_vs_12b_SR_LA_c2.0_cond2.bed_length)
diff_28_AF_RA_vs_20_SR_LA_c2.0_cond2.bed_length<- melt(diff_28_AF_RA_vs_20_SR_LA_c2.0_cond2.bed_length)
diff_28_AF_RA_vs_27_SR_LA_c2.0_cond2.bed_length<- melt(diff_28_AF_RA_vs_27_SR_LA_c2.0_cond2.bed_length)
diff_28_AF_RA_vs_45_SR_LA_c2.0_cond2.bed_length<- melt(diff_28_AF_RA_vs_45_SR_LA_c2.0_cond2.bed_length)
diff_31_AF_LA_vs_12b_SR_LA_c2.0_cond2.bed_length<- melt(diff_31_AF_LA_vs_12b_SR_LA_c2.0_cond2.bed_length)
diff_31_AF_LA_vs_20_SR_LA_c2.0_cond2.bed_length<- melt(diff_31_AF_LA_vs_20_SR_LA_c2.0_cond2.bed_length)
diff_31_AF_LA_vs_27_SR_LA_c2.0_cond2.bed_length<- melt(diff_31_AF_LA_vs_27_SR_LA_c2.0_cond2.bed_length)
diff_31_AF_LA_vs_45_SR_LA_c2.0_cond2.bed_length<- melt(diff_31_AF_LA_vs_45_SR_LA_c2.0_cond2.bed_length)
diff_31_AF_RA_vs_12b_SR_LA_c2.0_cond2.bed_length<- melt(diff_31_AF_RA_vs_12b_SR_LA_c2.0_cond2.bed_length)
diff_31_AF_RA_vs_20_SR_LA_c2.0_cond2.bed_length<- melt(diff_31_AF_RA_vs_20_SR_LA_c2.0_cond2.bed_length)
diff_31_AF_RA_vs_27_SR_LA_c2.0_cond2.bed_length<- melt(diff_31_AF_RA_vs_27_SR_LA_c2.0_cond2.bed_length)
diff_31_AF_RA_vs_45_SR_LA_c2.0_cond2.bed_length<- melt(diff_31_AF_RA_vs_45_SR_LA_c2.0_cond2.bed_length)
diff_313_AF_LA_vs_12b_SR_LA_c2.0_cond2.bed_length<- melt(diff_313_AF_LA_vs_12b_SR_LA_c2.0_cond2.bed_length)
diff_313_AF_LA_vs_20_SR_LA_c2.0_cond2.bed_length<- melt(diff_313_AF_LA_vs_20_SR_LA_c2.0_cond2.bed_length)
diff_313_AF_LA_vs_27_SR_LA_c2.0_cond2.bed_length<- melt(diff_313_AF_LA_vs_27_SR_LA_c2.0_cond2.bed_length)
diff_313_AF_LA_vs_45_SR_LA_c2.0_cond2.bed_length<- melt(diff_313_AF_LA_vs_45_SR_LA_c2.0_cond2.bed_length)
diff_313_AF_RA_vs_12b_SR_LA_c2.0_cond2.bed_length<- melt(diff_313_AF_RA_vs_12b_SR_LA_c2.0_cond2.bed_length)
diff_313_AF_RA_vs_20_SR_LA_c2.0_cond2.bed_length<- melt(diff_313_AF_RA_vs_20_SR_LA_c2.0_cond2.bed_length)
diff_313_AF_RA_vs_27_SR_LA_c2.0_cond2.bed_length<- melt(diff_313_AF_RA_vs_27_SR_LA_c2.0_cond2.bed_length)
diff_313_AF_RA_vs_45_SR_LA_c2.0_cond2.bed_length<- melt(diff_313_AF_RA_vs_45_SR_LA_c2.0_cond2.bed_length)
diff_33_SR_RA_vs_12b_SR_LA_c2.0_cond2.bed_length<- melt(diff_33_SR_RA_vs_12b_SR_LA_c2.0_cond2.bed_length)
diff_33_SR_RA_vs_20_SR_LA_c2.0_cond2.bed_length<- melt(diff_33_SR_RA_vs_20_SR_LA_c2.0_cond2.bed_length)
diff_33_SR_RA_vs_27_SR_LA_c2.0_cond2.bed_length<- melt(diff_33_SR_RA_vs_27_SR_LA_c2.0_cond2.bed_length)
diff_33_SR_RA_vs_45_SR_LA_c2.0_cond2.bed_length<- melt(diff_33_SR_RA_vs_45_SR_LA_c2.0_cond2.bed_length)
diff_47_AF_LA_vs_12b_SR_LA_c2.0_cond2.bed_length<- melt(diff_47_AF_LA_vs_12b_SR_LA_c2.0_cond2.bed_length)
diff_47_AF_LA_vs_20_SR_LA_c2.0_cond2.bed_length<- melt(diff_47_AF_LA_vs_20_SR_LA_c2.0_cond2.bed_length)
diff_47_AF_LA_vs_27_SR_LA_c2.0_cond2.bed_length<- melt(diff_47_AF_LA_vs_27_SR_LA_c2.0_cond2.bed_length)
diff_47_AF_LA_vs_45_SR_LA_c2.0_cond2.bed_length<- melt(diff_47_AF_LA_vs_45_SR_LA_c2.0_cond2.bed_length)

diff_12b_SR_RA_vs_12b_SR_LA_c2.0_cond2.bed_length$group <- "SR_LA"
diff_12b_SR_RA_vs_20_SR_LA_c2.0_cond2.bed_length$group <- "SR_LA"
diff_12b_SR_RA_vs_27_SR_LA_c2.0_cond2.bed_length$group <- "SR_LA"
diff_12b_SR_RA_vs_45_SR_LA_c2.0_cond2.bed_length$group <- "SR_LA"
diff_18_AF_RA_vs_12b_SR_LA_c2.0_cond2.bed_length$group <- "SR_LA"
diff_18_AF_RA_vs_20_SR_LA_c2.0_cond2.bed_length$group <- "SR_LA"
diff_18_AF_RA_vs_27_SR_LA_c2.0_cond2.bed_length$group <- "SR_LA"
diff_18_AF_RA_vs_45_SR_LA_c2.0_cond2.bed_length$group <- "SR_LA"
diff_20_SR_RA_vs_12b_SR_LA_c2.0_cond2.bed_length$group <- "SR_LA"
diff_20_SR_RA_vs_20_SR_LA_c2.0_cond2.bed_length$group <- "SR_LA"
diff_20_SR_RA_vs_27_SR_LA_c2.0_cond2.bed_length$group <- "SR_LA"
diff_20_SR_RA_vs_45_SR_LA_c2.0_cond2.bed_length$group <- "SR_LA"
diff_2231_AF_LA_vs_12b_SR_LA_c2.0_cond2.bed_length$group <- "SR_LA"
diff_2231_AF_LA_vs_20_SR_LA_c2.0_cond2.bed_length$group <- "SR_LA"
diff_2231_AF_LA_vs_27_SR_LA_c2.0_cond2.bed_length$group <- "SR_LA"
diff_2231_AF_LA_vs_45_SR_LA_c2.0_cond2.bed_length$group <- "SR_LA"
diff_2231_AF_RA_vs_12b_SR_LA_c2.0_cond2.bed_length$group <- "SR_LA"
diff_2231_AF_RA_vs_20_SR_LA_c2.0_cond2.bed_length$group <- "SR_LA"
diff_2231_AF_RA_vs_27_SR_LA_c2.0_cond2.bed_length$group <- "SR_LA"
diff_2231_AF_RA_vs_45_SR_LA_c2.0_cond2.bed_length$group <- "SR_LA"
diff_25_AF_LA_vs_12b_SR_LA_c2.0_cond2.bed_length$group <- "SR_LA"
diff_25_AF_LA_vs_20_SR_LA_c2.0_cond2.bed_length$group <- "SR_LA"
diff_25_AF_LA_vs_27_SR_LA_c2.0_cond2.bed_length$group <- "SR_LA"
diff_25_AF_LA_vs_45_SR_LA_c2.0_cond2.bed_length$group <- "SR_LA"
diff_25_AF_RA_vs_12b_SR_LA_c2.0_cond2.bed_length$group <- "SR_LA"
diff_25_AF_RA_vs_20_SR_LA_c2.0_cond2.bed_length$group <- "SR_LA"
diff_25_AF_RA_vs_27_SR_LA_c2.0_cond2.bed_length$group <- "SR_LA"
diff_25_AF_RA_vs_45_SR_LA_c2.0_cond2.bed_length$group <- "SR_LA"
diff_27_SR_RA_vs_12b_SR_LA_c2.0_cond2.bed_length$group <- "SR_LA"
diff_27_SR_RA_vs_20_SR_LA_c2.0_cond2.bed_length$group <- "SR_LA"
diff_27_SR_RA_vs_27_SR_LA_c2.0_cond2.bed_length$group <- "SR_LA"
diff_27_SR_RA_vs_45_SR_LA_c2.0_cond2.bed_length$group <- "SR_LA"
diff_28_AF_LA_vs_12b_SR_LA_c2.0_cond2.bed_length$group <- "SR_LA"
diff_28_AF_LA_vs_20_SR_LA_c2.0_cond2.bed_length$group <- "SR_LA"
diff_28_AF_LA_vs_27_SR_LA_c2.0_cond2.bed_length$group <- "SR_LA"
diff_28_AF_LA_vs_45_SR_LA_c2.0_cond2.bed_length$group <- "SR_LA"
diff_28_AF_RA_vs_12b_SR_LA_c2.0_cond2.bed_length$group <- "SR_LA"
diff_28_AF_RA_vs_20_SR_LA_c2.0_cond2.bed_length$group <- "SR_LA"
diff_28_AF_RA_vs_27_SR_LA_c2.0_cond2.bed_length$group <- "SR_LA"
diff_28_AF_RA_vs_45_SR_LA_c2.0_cond2.bed_length$group <- "SR_LA"
diff_31_AF_LA_vs_12b_SR_LA_c2.0_cond2.bed_length$group <- "SR_LA"
diff_31_AF_LA_vs_20_SR_LA_c2.0_cond2.bed_length$group <- "SR_LA"
diff_31_AF_LA_vs_27_SR_LA_c2.0_cond2.bed_length$group <- "SR_LA"
diff_31_AF_LA_vs_45_SR_LA_c2.0_cond2.bed_length$group <- "SR_LA"
diff_31_AF_RA_vs_12b_SR_LA_c2.0_cond2.bed_length$group <- "SR_LA"
diff_31_AF_RA_vs_20_SR_LA_c2.0_cond2.bed_length$group <- "SR_LA"
diff_31_AF_RA_vs_27_SR_LA_c2.0_cond2.bed_length$group <- "SR_LA"
diff_31_AF_RA_vs_45_SR_LA_c2.0_cond2.bed_length$group <- "SR_LA"
diff_313_AF_LA_vs_12b_SR_LA_c2.0_cond2.bed_length$group <- "SR_LA"
diff_313_AF_LA_vs_20_SR_LA_c2.0_cond2.bed_length$group <- "SR_LA"
diff_313_AF_LA_vs_27_SR_LA_c2.0_cond2.bed_length$group <- "SR_LA"
diff_313_AF_LA_vs_45_SR_LA_c2.0_cond2.bed_length$group <- "SR_LA"
diff_313_AF_RA_vs_12b_SR_LA_c2.0_cond2.bed_length$group <- "SR_LA"
diff_313_AF_RA_vs_20_SR_LA_c2.0_cond2.bed_length$group <- "SR_LA"
diff_313_AF_RA_vs_27_SR_LA_c2.0_cond2.bed_length$group <- "SR_LA"
diff_313_AF_RA_vs_45_SR_LA_c2.0_cond2.bed_length$group <- "SR_LA"
diff_33_SR_RA_vs_12b_SR_LA_c2.0_cond2.bed_length$group <- "SR_LA"
diff_33_SR_RA_vs_20_SR_LA_c2.0_cond2.bed_length$group <- "SR_LA"
diff_33_SR_RA_vs_27_SR_LA_c2.0_cond2.bed_length$group <- "SR_LA"
diff_33_SR_RA_vs_45_SR_LA_c2.0_cond2.bed_length$group <- "SR_LA"
diff_47_AF_LA_vs_12b_SR_LA_c2.0_cond2.bed_length$group <- "SR_LA"
diff_47_AF_LA_vs_20_SR_LA_c2.0_cond2.bed_length$group <- "SR_LA"
diff_47_AF_LA_vs_27_SR_LA_c2.0_cond2.bed_length$group <- "SR_LA"
diff_47_AF_LA_vs_45_SR_LA_c2.0_cond2.bed_length$group <- "SR_LA"

diff_12b_SR_RA_vs_12b_SR_LA_c2.0_cond2.bed_length$sample <- "diff_12b_SR_RA_vs_12b_SR_LA_c2.0_cond2.bed"
diff_12b_SR_RA_vs_20_SR_LA_c2.0_cond2.bed_length$sample <- "diff_12b_SR_RA_vs_20_SR_LA_c2.0_cond2.bed"
diff_12b_SR_RA_vs_27_SR_LA_c2.0_cond2.bed_length$sample <- "diff_12b_SR_RA_vs_27_SR_LA_c2.0_cond2.bed"
diff_12b_SR_RA_vs_45_SR_LA_c2.0_cond2.bed_length$sample <- "diff_12b_SR_RA_vs_45_SR_LA_c2.0_cond2.bed"
diff_18_AF_RA_vs_12b_SR_LA_c2.0_cond2.bed_length$sample <- "diff_18_AF_RA_vs_12b_SR_LA_c2.0_cond2.bed"
diff_18_AF_RA_vs_20_SR_LA_c2.0_cond2.bed_length$sample <- "diff_18_AF_RA_vs_20_SR_LA_c2.0_cond2.bed"
diff_18_AF_RA_vs_27_SR_LA_c2.0_cond2.bed_length$sample <- "diff_18_AF_RA_vs_27_SR_LA_c2.0_cond2.bed"
diff_18_AF_RA_vs_45_SR_LA_c2.0_cond2.bed_length$sample <- "diff_18_AF_RA_vs_45_SR_LA_c2.0_cond2.bed"
diff_20_SR_RA_vs_12b_SR_LA_c2.0_cond2.bed_length$sample <- "diff_20_SR_RA_vs_12b_SR_LA_c2.0_cond2.bed"
diff_20_SR_RA_vs_20_SR_LA_c2.0_cond2.bed_length$sample <- "diff_20_SR_RA_vs_20_SR_LA_c2.0_cond2.bed"
diff_20_SR_RA_vs_27_SR_LA_c2.0_cond2.bed_length$sample <- "diff_20_SR_RA_vs_27_SR_LA_c2.0_cond2.bed"
diff_20_SR_RA_vs_45_SR_LA_c2.0_cond2.bed_length$sample <- "diff_20_SR_RA_vs_45_SR_LA_c2.0_cond2.bed"
diff_2231_AF_LA_vs_12b_SR_LA_c2.0_cond2.bed_length$sample <- "diff_2231_AF_LA_vs_12b_SR_LA_c2.0_cond2.bed"
diff_2231_AF_LA_vs_20_SR_LA_c2.0_cond2.bed_length$sample <- "diff_2231_AF_LA_vs_20_SR_LA_c2.0_cond2.bed"
diff_2231_AF_LA_vs_27_SR_LA_c2.0_cond2.bed_length$sample <- "diff_2231_AF_LA_vs_27_SR_LA_c2.0_cond2.bed"
diff_2231_AF_LA_vs_45_SR_LA_c2.0_cond2.bed_length$sample <- "diff_2231_AF_RA_vs_45_SR_LA_c2.0_cond2.bed"
diff_2231_AF_RA_vs_12b_SR_LA_c2.0_cond2.bed_length$sample <- "diff_25_AF_LA_vs_12b_SR_LA_c2.0_cond2.bed"
diff_2231_AF_RA_vs_20_SR_LA_c2.0_cond2.bed_length$sample <- "diff_25_AF_LA_vs_20_SR_LA_c2.0_cond2.bed"
diff_2231_AF_RA_vs_27_SR_LA_c2.0_cond2.bed_length$sample <- "diff_25_AF_LA_vs_27_SR_LA_c2.0_cond2.bed"
diff_2231_AF_RA_vs_45_SR_LA_c2.0_cond2.bed_length$sample <- "diff_25_AF_LA_vs_45_SR_LA_c2.0_cond2.bed"
diff_25_AF_LA_vs_12b_SR_LA_c2.0_cond2.bed_length$sample <- "diff_25_AF_RA_vs_12b_SR_LA_c2.0_cond2.bed"
diff_25_AF_LA_vs_20_SR_LA_c2.0_cond2.bed_length$sample <- "diff_25_AF_RA_vs_20_SR_LA_c2.0_cond2.bed"
diff_25_AF_LA_vs_27_SR_LA_c2.0_cond2.bed_length$sample <- "diff_25_AF_RA_vs_27_SR_LA_c2.0_cond2.bed"
diff_25_AF_LA_vs_45_SR_LA_c2.0_cond2.bed_length$sample <- "diff_25_AF_LA_vs_45_SR_LA_c2.0_cond2.bed"
diff_25_AF_RA_vs_12b_SR_LA_c2.0_cond2.bed_length$sample <- "diff_25_AF_RA_vs_12b_SR_LA_c2.0_cond2.bed"
diff_25_AF_RA_vs_20_SR_LA_c2.0_cond2.bed_length$sample <- "diff_25_AF_RA_vs_20_SR_LA_c2.0_cond2.bed"
diff_25_AF_RA_vs_27_SR_LA_c2.0_cond2.bed_length$sample <- "diff_25_AF_RA_vs_27_SR_LA_c2.0_cond2.bed"
diff_25_AF_RA_vs_45_SR_LA_c2.0_cond2.bed_length$sample <- "diff_25_AF_RA_vs_45_SR_LA_c2.0_cond2.bed"
diff_27_SR_RA_vs_12b_SR_LA_c2.0_cond2.bed_length$sample <- "diff_27_SR_RA_vs_12b_SR_LA_c2.0_cond2.bed"
diff_27_SR_RA_vs_20_SR_LA_c2.0_cond2.bed_length$sample <- "diff_27_SR_RA_vs_20_SR_LA_c2.0_cond2.bed"
diff_27_SR_RA_vs_27_SR_LA_c2.0_cond2.bed_length$sample <- "diff_27_SR_RA_vs_27_SR_LA_c2.0_cond2.bed"
diff_27_SR_RA_vs_45_SR_LA_c2.0_cond2.bed_length$sample <- "diff_27_SR_RA_vs_45_SR_LA_c2.0_cond2.bed"
diff_28_AF_LA_vs_12b_SR_LA_c2.0_cond2.bed_length$sample <- "diff_28_AF_LA_vs_12b_SR_LA_c2.0_cond2.bed"
diff_28_AF_LA_vs_20_SR_LA_c2.0_cond2.bed_length$sample <- "diff_28_AF_LA_vs_20_SR_LA_c2.0_cond2.bed"
diff_28_AF_LA_vs_27_SR_LA_c2.0_cond2.bed_length$sample <- "diff_28_AF_LA_vs_27_SR_LA_c2.0_cond2.bed"
diff_28_AF_LA_vs_45_SR_LA_c2.0_cond2.bed_length$sample <- "diff_28_AF_LA_vs_45_SR_LA_c2.0_cond2.bed"
diff_28_AF_RA_vs_12b_SR_LA_c2.0_cond2.bed_length$sample <- "diff_28_AF_RA_vs_12b_SR_LA_c2.0_cond2.bed"
diff_28_AF_RA_vs_20_SR_LA_c2.0_cond2.bed_length$sample <- "diff_28_AF_RA_vs_20_SR_LA_c2.0_cond2.bed"
diff_28_AF_RA_vs_27_SR_LA_c2.0_cond2.bed_length$sample <- "diff_28_AF_RA_vs_27_SR_LA_c2.0_cond2.bed"
diff_28_AF_RA_vs_45_SR_LA_c2.0_cond2.bed_length$sample <- "diff_28_AF_RA_vs_45_SR_LA_c2.0_cond2.bed"
diff_31_AF_LA_vs_12b_SR_LA_c2.0_cond2.bed_length$sample <- "diff_31_AF_LA_vs_12b_SR_LA_c2.0_cond2.bed"
diff_31_AF_LA_vs_20_SR_LA_c2.0_cond2.bed_length$sample <- "diff_31_AF_LA_vs_20_SR_LA_c2.0_cond2.bed"
diff_31_AF_LA_vs_27_SR_LA_c2.0_cond2.bed_length$sample <- "diff_31_AF_LA_vs_27_SR_LA_c2.0_cond2.bed"
diff_31_AF_LA_vs_45_SR_LA_c2.0_cond2.bed_length$sample <- "diff_31_AF_LA_vs_45_SR_LA_c2.0_cond2.bed"
diff_31_AF_RA_vs_12b_SR_LA_c2.0_cond2.bed_length$sample <- "diff_31_AF_RA_vs_12b_SR_LA_c2.0_cond2.bed"
diff_31_AF_RA_vs_20_SR_LA_c2.0_cond2.bed_length$sample <- "diff_31_AF_RA_vs_20_SR_LA_c2.0_cond2.bed"
diff_31_AF_RA_vs_27_SR_LA_c2.0_cond2.bed_length$sample <- "diff_31_AF_RA_vs_27_SR_LA_c2.0_cond2.bed"
diff_31_AF_RA_vs_45_SR_LA_c2.0_cond2.bed_length$sample <- "diff_31_AF_RA_vs_45_SR_LA_c2.0_cond2.bed"
diff_313_AF_LA_vs_12b_SR_LA_c2.0_cond2.bed_length$sample <- "diff_313_AF_LA_vs_12b_SR_LA_c2.0_cond2.bed"
diff_313_AF_LA_vs_20_SR_LA_c2.0_cond2.bed_length$sample <- "diff_313_AF_LA_vs_20_SR_LA_c2.0_cond2.bed"
diff_313_AF_LA_vs_27_SR_LA_c2.0_cond2.bed_length$sample <- "diff_313_AF_LA_vs_27_SR_LA_c2.0_cond2.bed"
diff_313_AF_LA_vs_45_SR_LA_c2.0_cond2.bed_length$sample <- "diff_313_AF_LA_vs_45_SR_LA_c2.0_cond2.bed"
diff_313_AF_RA_vs_12b_SR_LA_c2.0_cond2.bed_length$sample <- "diff_313_AF_RA_vs_12b_SR_LA_c2.0_cond2.bed"
diff_313_AF_RA_vs_20_SR_LA_c2.0_cond2.bed_length$sample <- "diff_313_AF_RA_vs_20_SR_LA_c2.0_cond2.bed"
diff_313_AF_RA_vs_27_SR_LA_c2.0_cond2.bed_length$sample <- "diff_313_AF_RA_vs_27_SR_LA_c2.0_cond2.bed"
diff_313_AF_RA_vs_45_SR_LA_c2.0_cond2.bed_length$sample <- "diff_313_AF_RA_vs_45_SR_LA_c2.0_cond2.bed"
diff_33_SR_RA_vs_12b_SR_LA_c2.0_cond2.bed_length$sample <- "diff_33_SR_RA_vs_12b_SR_LA_c2.0_cond2.bed"
diff_33_SR_RA_vs_20_SR_LA_c2.0_cond2.bed_length$sample <- "diff_33_SR_RA_vs_20_SR_LA_c2.0_cond2.bed"
diff_33_SR_RA_vs_27_SR_LA_c2.0_cond2.bed_length$sample <- "diff_33_SR_RA_vs_27_SR_LA_c2.0_cond2.bed"
diff_33_SR_RA_vs_45_SR_LA_c2.0_cond2.bed_length$sample <- "diff_33_SR_RA_vs_45_SR_LA_c2.0_cond2.bed"
diff_47_AF_LA_vs_12b_SR_LA_c2.0_cond2.bed_length$sample <- "diff_47_AF_LA_vs_12b_SR_LA_c2.0_cond2.bed"
diff_47_AF_LA_vs_20_SR_LA_c2.0_cond2.bed_length$sample <- "diff_47_AF_LA_vs_20_SR_LA_c2.0_cond2.bed"
diff_47_AF_LA_vs_27_SR_LA_c2.0_cond2.bed_length$sample <- "diff_47_AF_LA_vs_27_SR_LA_c2.0_cond2.bed"
diff_47_AF_LA_vs_45_SR_LA_c2.0_cond2.bed_length$sample <- "diff_47_AF_LA_vs_45_SR_LA_c2.0_cond2.bed"


peaksets_lengths_SR_LA_pairwise <- rbind(diff_12b_SR_RA_vs_12b_SR_LA_c2.0_cond2.bed_length,
                                         diff_12b_SR_RA_vs_20_SR_LA_c2.0_cond2.bed_length,
                                         diff_12b_SR_RA_vs_27_SR_LA_c2.0_cond2.bed_length,
                                         diff_12b_SR_RA_vs_45_SR_LA_c2.0_cond2.bed_length,
                                         diff_18_AF_RA_vs_12b_SR_LA_c2.0_cond2.bed_length,
                                         diff_18_AF_RA_vs_20_SR_LA_c2.0_cond2.bed_length,
                                         diff_18_AF_RA_vs_27_SR_LA_c2.0_cond2.bed_length,
                                         diff_18_AF_RA_vs_45_SR_LA_c2.0_cond2.bed_length,
                                         diff_20_SR_RA_vs_12b_SR_LA_c2.0_cond2.bed_length,
                                         diff_20_SR_RA_vs_20_SR_LA_c2.0_cond2.bed_length,
                                         diff_20_SR_RA_vs_27_SR_LA_c2.0_cond2.bed_length,
                                         diff_20_SR_RA_vs_45_SR_LA_c2.0_cond2.bed_length,
                                         diff_2231_AF_LA_vs_12b_SR_LA_c2.0_cond2.bed_length,
                                         diff_2231_AF_LA_vs_20_SR_LA_c2.0_cond2.bed_length,
                                         diff_2231_AF_LA_vs_27_SR_LA_c2.0_cond2.bed_length,
                                         diff_2231_AF_LA_vs_45_SR_LA_c2.0_cond2.bed_length,
                                         diff_2231_AF_RA_vs_12b_SR_LA_c2.0_cond2.bed_length,
                                         diff_2231_AF_RA_vs_20_SR_LA_c2.0_cond2.bed_length,
                                         diff_2231_AF_RA_vs_27_SR_LA_c2.0_cond2.bed_length,
                                         diff_2231_AF_RA_vs_45_SR_LA_c2.0_cond2.bed_length,
                                         diff_25_AF_LA_vs_12b_SR_LA_c2.0_cond2.bed_length,
                                         diff_25_AF_LA_vs_20_SR_LA_c2.0_cond2.bed_length,
                                         diff_25_AF_LA_vs_27_SR_LA_c2.0_cond2.bed_length,
                                         diff_25_AF_LA_vs_45_SR_LA_c2.0_cond2.bed_length,
                                         diff_25_AF_RA_vs_12b_SR_LA_c2.0_cond2.bed_length,
                                         diff_25_AF_RA_vs_20_SR_LA_c2.0_cond2.bed_length,
                                         diff_25_AF_RA_vs_27_SR_LA_c2.0_cond2.bed_length,
                                         diff_25_AF_RA_vs_45_SR_LA_c2.0_cond2.bed_length,
                                         diff_27_SR_RA_vs_12b_SR_LA_c2.0_cond2.bed_length,
                                         diff_27_SR_RA_vs_20_SR_LA_c2.0_cond2.bed_length,
                                         diff_27_SR_RA_vs_27_SR_LA_c2.0_cond2.bed_length,
                                         diff_27_SR_RA_vs_45_SR_LA_c2.0_cond2.bed_length,
                                         diff_28_AF_LA_vs_12b_SR_LA_c2.0_cond2.bed_length,
                                         diff_28_AF_LA_vs_20_SR_LA_c2.0_cond2.bed_length,
                                         diff_28_AF_LA_vs_27_SR_LA_c2.0_cond2.bed_length,
                                         diff_28_AF_LA_vs_45_SR_LA_c2.0_cond2.bed_length,
                                         diff_28_AF_RA_vs_12b_SR_LA_c2.0_cond2.bed_length,
                                         diff_28_AF_RA_vs_20_SR_LA_c2.0_cond2.bed_length,
                                         diff_28_AF_RA_vs_27_SR_LA_c2.0_cond2.bed_length,
                                         diff_28_AF_RA_vs_45_SR_LA_c2.0_cond2.bed_length,
                                         diff_31_AF_LA_vs_12b_SR_LA_c2.0_cond2.bed_length,
                                         diff_31_AF_LA_vs_20_SR_LA_c2.0_cond2.bed_length,
                                         diff_31_AF_LA_vs_27_SR_LA_c2.0_cond2.bed_length,
                                         diff_31_AF_LA_vs_45_SR_LA_c2.0_cond2.bed_length,
                                         diff_31_AF_RA_vs_12b_SR_LA_c2.0_cond2.bed_length,
                                         diff_31_AF_RA_vs_20_SR_LA_c2.0_cond2.bed_length,
                                         diff_31_AF_RA_vs_27_SR_LA_c2.0_cond2.bed_length,
                                         diff_31_AF_RA_vs_45_SR_LA_c2.0_cond2.bed_length,
                                         diff_313_AF_LA_vs_12b_SR_LA_c2.0_cond2.bed_length,
                                         diff_313_AF_LA_vs_20_SR_LA_c2.0_cond2.bed_length,
                                         diff_313_AF_LA_vs_27_SR_LA_c2.0_cond2.bed_length,
                                         diff_313_AF_LA_vs_45_SR_LA_c2.0_cond2.bed_length,
                                         diff_313_AF_RA_vs_12b_SR_LA_c2.0_cond2.bed_length,
                                         diff_313_AF_RA_vs_20_SR_LA_c2.0_cond2.bed_length,
                                         diff_313_AF_RA_vs_27_SR_LA_c2.0_cond2.bed_length,
                                         diff_313_AF_RA_vs_45_SR_LA_c2.0_cond2.bed_length,
                                         diff_33_SR_RA_vs_12b_SR_LA_c2.0_cond2.bed_length,
                                         diff_33_SR_RA_vs_20_SR_LA_c2.0_cond2.bed_length,
                                         diff_33_SR_RA_vs_27_SR_LA_c2.0_cond2.bed_length,
                                         diff_33_SR_RA_vs_45_SR_LA_c2.0_cond2.bed_length,
                                         diff_47_AF_LA_vs_12b_SR_LA_c2.0_cond2.bed_length,
                                         diff_47_AF_LA_vs_20_SR_LA_c2.0_cond2.bed_length,
                                         diff_47_AF_LA_vs_27_SR_LA_c2.0_cond2.bed_length,
                                         diff_47_AF_LA_vs_45_SR_LA_c2.0_cond2.bed_length)

sample_peakset_lengths_boxplot <- ggplot(peaksets_lengths_SR_LA_pairwise, aes(x=sample, y=value)) + geom_boxplot(aes(fill=group),outlier.size = 0.05) + scale_y_continuous("Peak length (bp)") + scale_x_discrete("enriched_peakset_comparison") + theme_classic() + theme(axis.text.x = element_text(angle=45, hjust = 1) )  + scale_fill_manual(values=c('#7298BC')) + theme(axis.text.x = element_blank())                                                                                                                                                                                                                                                                                                                                                                                                                                                                      

sample_peakset_lengths_boxplot
sample_peakset_lengths_boxplot_nooutliers <- ggplot(peaksets_lengths_SR_LA_pairwise, aes(x=sample, y=value)) + geom_boxplot(aes(fill=group),outlier.shape=NA) + coord_cartesian(ylim = c(0, 1500))+ scale_y_continuous("Peak length (bp)") + scale_x_discrete("enriched_peakset_comparison") + theme_classic() + theme(axis.text.x = element_text(angle=45, hjust = 1) ) + scale_fill_manual(values=c('#7298BC'))        
sample_peakset_lengths_boxplot_nooutliers
sample_peakset_lengths_boxplot_nooutliers <- ggplot(peaksets_lengths_SR_LA_pairwise, aes(x=sample, y=value)) + geom_boxplot(aes(fill=group),outlier.shape=NA) + coord_cartesian(ylim = c(300, 1500))+ scale_y_continuous("Peak length (bp)") + scale_x_discrete("enriched_peakset_comparison") + theme_classic() + theme(axis.text.x = element_blank() ) + scale_fill_manual(values=c('#7298BC'))        
sample_peakset_lengths_boxplot_nooutliers



#SR-RA

setwd("/Rodriguez_2025_AFepigenome/H3K27ac_enriched_regions_definition/SR_RA_enriched_peaks")
myFiles <- list.files(pattern="*.bed")



diff_12b_SR_RA_vs_12b_SR_LA_c2.0_cond1.bed<-as.data.frame(read.table("diff_12b_SR_RA_vs_12b_SR_LA_c2.0_cond1.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_12b_SR_RA_vs_20_SR_LA_c2.0_cond1.bed<-as.data.frame(read.table("diff_12b_SR_RA_vs_20_SR_LA_c2.0_cond1.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_12b_SR_RA_vs_2231_AF_LA_c2.0_cond1.bed<-as.data.frame(read.table("diff_12b_SR_RA_vs_2231_AF_LA_c2.0_cond1.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_12b_SR_RA_vs_25_AF_LA_c2.0_cond1.bed<-as.data.frame(read.table("diff_12b_SR_RA_vs_25_AF_LA_c2.0_cond1.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_12b_SR_RA_vs_27_SR_LA_c2.0_cond1.bed<-as.data.frame(read.table("diff_12b_SR_RA_vs_27_SR_LA_c2.0_cond1.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_12b_SR_RA_vs_28_AF_LA_c2.0_cond1.bed<-as.data.frame(read.table("diff_12b_SR_RA_vs_28_AF_LA_c2.0_cond1.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_12b_SR_RA_vs_31_AF_LA_c2.0_cond1.bed<-as.data.frame(read.table("diff_12b_SR_RA_vs_31_AF_LA_c2.0_cond1.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_12b_SR_RA_vs_313_AF_LA_c2.0_cond1.bed<-as.data.frame(read.table("diff_12b_SR_RA_vs_313_AF_LA_c2.0_cond1.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_12b_SR_RA_vs_45_SR_LA_c2.0_cond1.bed<-as.data.frame(read.table("diff_12b_SR_RA_vs_45_SR_LA_c2.0_cond1.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_12b_SR_RA_vs_47_AF_LA_c2.0_cond1.bed<-as.data.frame(read.table("diff_12b_SR_RA_vs_47_AF_LA_c2.0_cond1.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_18_AF_RA_vs_12b_SR_RA_c2.0_cond2.bed<-as.data.frame(read.table("diff_18_AF_RA_vs_12b_SR_RA_c2.0_cond2.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_18_AF_RA_vs_20_SR_RA_c2.0_cond2.bed<-as.data.frame(read.table("diff_18_AF_RA_vs_20_SR_RA_c2.0_cond2.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_18_AF_RA_vs_27_SR_RA_c2.0_cond2.bed<-as.data.frame(read.table("diff_18_AF_RA_vs_27_SR_RA_c2.0_cond2.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_18_AF_RA_vs_33_SR_RA_c2.0_cond2.bed<-as.data.frame(read.table("diff_18_AF_RA_vs_33_SR_RA_c2.0_cond2.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_20_SR_RA_vs_12b_SR_LA_c2.0_cond1.bed<-as.data.frame(read.table("diff_20_SR_RA_vs_12b_SR_LA_c2.0_cond1.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_20_SR_RA_vs_20_SR_LA_c2.0_cond1.bed<-as.data.frame(read.table("diff_20_SR_RA_vs_20_SR_LA_c2.0_cond1.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_20_SR_RA_vs_2231_AF_LA_c2.0_cond1.bed<-as.data.frame(read.table("diff_20_SR_RA_vs_2231_AF_LA_c2.0_cond1.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_20_SR_RA_vs_25_AF_LA_c2.0_cond1.bed<-as.data.frame(read.table("diff_20_SR_RA_vs_25_AF_LA_c2.0_cond1.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_20_SR_RA_vs_27_SR_LA_c2.0_cond1.bed<-as.data.frame(read.table("diff_20_SR_RA_vs_27_SR_LA_c2.0_cond1.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_20_SR_RA_vs_28_AF_LA_c2.0_cond1.bed<-as.data.frame(read.table("diff_20_SR_RA_vs_28_AF_LA_c2.0_cond1.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_20_SR_RA_vs_31_AF_LA_c2.0_cond1.bed<-as.data.frame(read.table("diff_20_SR_RA_vs_31_AF_LA_c2.0_cond1.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_20_SR_RA_vs_313_AF_LA_c2.0_cond1.bed<-as.data.frame(read.table("diff_20_SR_RA_vs_313_AF_LA_c2.0_cond1.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_20_SR_RA_vs_45_SR_LA_c2.0_cond1.bed<-as.data.frame(read.table("diff_20_SR_RA_vs_45_SR_LA_c2.0_cond1.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_20_SR_RA_vs_47_AF_LA_c2.0_cond1.bed<-as.data.frame(read.table("diff_20_SR_RA_vs_47_AF_LA_c2.0_cond1.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_2231_AF_RA_vs_12b_SR_RA_c2.0_cond2.bed<-as.data.frame(read.table("diff_2231_AF_RA_vs_12b_SR_RA_c2.0_cond2.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_2231_AF_RA_vs_20_SR_RA_c2.0_cond2.bed<-as.data.frame(read.table("diff_2231_AF_RA_vs_20_SR_RA_c2.0_cond2.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_2231_AF_RA_vs_27_SR_RA_c2.0_cond2.bed<-as.data.frame(read.table("diff_2231_AF_RA_vs_27_SR_RA_c2.0_cond2.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_2231_AF_RA_vs_33_SR_RA_c2.0_cond2.bed<-as.data.frame(read.table("diff_2231_AF_RA_vs_33_SR_RA_c2.0_cond2.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_25_AF_RA_vs_12b_SR_RA_c2.0_cond2.bed<-as.data.frame(read.table("diff_25_AF_RA_vs_12b_SR_RA_c2.0_cond2.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_25_AF_RA_vs_20_SR_RA_c2.0_cond2.bed<-as.data.frame(read.table("diff_25_AF_RA_vs_20_SR_RA_c2.0_cond2.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_25_AF_RA_vs_27_SR_RA_c2.0_cond2.bed<-as.data.frame(read.table("diff_25_AF_RA_vs_27_SR_RA_c2.0_cond2.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_25_AF_RA_vs_33_SR_RA_c2.0_cond2.bed<-as.data.frame(read.table("diff_25_AF_RA_vs_33_SR_RA_c2.0_cond2.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_27_SR_RA_vs_12b_SR_LA_c2.0_cond1.bed<-as.data.frame(read.table("diff_27_SR_RA_vs_12b_SR_LA_c2.0_cond1.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_27_SR_RA_vs_20_SR_LA_c2.0_cond1.bed<-as.data.frame(read.table("diff_27_SR_RA_vs_20_SR_LA_c2.0_cond1.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_27_SR_RA_vs_2231_AF_LA_c2.0_cond1.bed<-as.data.frame(read.table("diff_27_SR_RA_vs_2231_AF_LA_c2.0_cond1.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_27_SR_RA_vs_25_AF_LA_c2.0_cond1.bed<-as.data.frame(read.table("diff_27_SR_RA_vs_25_AF_LA_c2.0_cond1.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_27_SR_RA_vs_27_SR_LA_c2.0_cond1.bed<-as.data.frame(read.table("diff_27_SR_RA_vs_27_SR_LA_c2.0_cond1.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_27_SR_RA_vs_28_AF_LA_c2.0_cond1.bed<-as.data.frame(read.table("diff_27_SR_RA_vs_28_AF_LA_c2.0_cond1.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_27_SR_RA_vs_31_AF_LA_c2.0_cond1.bed<-as.data.frame(read.table("diff_27_SR_RA_vs_31_AF_LA_c2.0_cond1.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_27_SR_RA_vs_313_AF_LA_c2.0_cond1.bed<-as.data.frame(read.table("diff_27_SR_RA_vs_313_AF_LA_c2.0_cond1.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_27_SR_RA_vs_45_SR_LA_c2.0_cond1.bed<-as.data.frame(read.table("diff_27_SR_RA_vs_45_SR_LA_c2.0_cond1.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_27_SR_RA_vs_47_AF_LA_c2.0_cond1.bed<-as.data.frame(read.table("diff_27_SR_RA_vs_47_AF_LA_c2.0_cond1.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_28_AF_RA_vs_12b_SR_RA_c2.0_cond2.bed<-as.data.frame(read.table("diff_28_AF_RA_vs_12b_SR_RA_c2.0_cond2.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_28_AF_RA_vs_20_SR_RA_c2.0_cond2.bed<-as.data.frame(read.table("diff_28_AF_RA_vs_20_SR_RA_c2.0_cond2.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_28_AF_RA_vs_27_SR_RA_c2.0_cond2.bed<-as.data.frame(read.table("diff_28_AF_RA_vs_27_SR_RA_c2.0_cond2.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_28_AF_RA_vs_33_SR_RA_c2.0_cond2.bed<-as.data.frame(read.table("diff_28_AF_RA_vs_33_SR_RA_c2.0_cond2.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_31_AF_RA_vs_12b_SR_RA_c2.0_cond2.bed<-as.data.frame(read.table("diff_31_AF_RA_vs_12b_SR_RA_c2.0_cond2.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_31_AF_RA_vs_20_SR_RA_c2.0_cond2.bed<-as.data.frame(read.table("diff_31_AF_RA_vs_20_SR_RA_c2.0_cond2.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_31_AF_RA_vs_27_SR_RA_c2.0_cond2.bed<-as.data.frame(read.table("diff_31_AF_RA_vs_27_SR_RA_c2.0_cond2.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_31_AF_RA_vs_33_SR_RA_c2.0_cond2.bed<-as.data.frame(read.table("diff_31_AF_RA_vs_33_SR_RA_c2.0_cond2.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_313_AF_RA_vs_12b_SR_RA_c2.0_cond2.bed<-as.data.frame(read.table("diff_313_AF_RA_vs_12b_SR_RA_c2.0_cond2.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_313_AF_RA_vs_20_SR_RA_c2.0_cond2.bed<-as.data.frame(read.table("diff_313_AF_RA_vs_20_SR_RA_c2.0_cond2.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_313_AF_RA_vs_27_SR_RA_c2.0_cond2.bed<-as.data.frame(read.table("diff_313_AF_RA_vs_27_SR_RA_c2.0_cond2.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_313_AF_RA_vs_33_SR_RA_c2.0_cond2.bed<-as.data.frame(read.table("diff_313_AF_RA_vs_33_SR_RA_c2.0_cond2.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_33_SR_RA_vs_12b_SR_LA_c2.0_cond1.bed<-as.data.frame(read.table("diff_33_SR_RA_vs_12b_SR_LA_c2.0_cond1.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_33_SR_RA_vs_20_SR_LA_c2.0_cond1.bed<-as.data.frame(read.table("diff_33_SR_RA_vs_20_SR_LA_c2.0_cond1.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_33_SR_RA_vs_2231_AF_LA_c2.0_cond1.bed<-as.data.frame(read.table("diff_33_SR_RA_vs_2231_AF_LA_c2.0_cond1.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_33_SR_RA_vs_25_AF_LA_c2.0_cond1.bed<-as.data.frame(read.table("diff_33_SR_RA_vs_25_AF_LA_c2.0_cond1.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_33_SR_RA_vs_27_SR_LA_c2.0_cond1.bed<-as.data.frame(read.table("diff_33_SR_RA_vs_27_SR_LA_c2.0_cond1.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_33_SR_RA_vs_28_AF_LA_c2.0_cond1.bed<-as.data.frame(read.table("diff_33_SR_RA_vs_28_AF_LA_c2.0_cond1.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_33_SR_RA_vs_31_AF_LA_c2.0_cond1.bed<-as.data.frame(read.table("diff_33_SR_RA_vs_31_AF_LA_c2.0_cond1.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_33_SR_RA_vs_313_AF_LA_c2.0_cond1.bed<-as.data.frame(read.table("diff_33_SR_RA_vs_313_AF_LA_c2.0_cond1.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_33_SR_RA_vs_45_SR_LA_c2.0_cond1.bed<-as.data.frame(read.table("diff_33_SR_RA_vs_45_SR_LA_c2.0_cond1.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))
diff_33_SR_RA_vs_47_AF_LA_c2.0_cond1.bed<-as.data.frame(read.table("diff_33_SR_RA_vs_47_AF_LA_c2.0_cond1.bed",skip=1, header = FALSE, sep="\t",stringsAsFactors=FALSE))


diff_12b_SR_RA_vs_12b_SR_LA_c2.0_cond1.bed_length<- diff_12b_SR_RA_vs_12b_SR_LA_c2.0_cond1.bed$V3 - diff_12b_SR_RA_vs_12b_SR_LA_c2.0_cond1.bed$V2
diff_12b_SR_RA_vs_20_SR_LA_c2.0_cond1.bed_length<- diff_12b_SR_RA_vs_20_SR_LA_c2.0_cond1.bed$V3 - diff_12b_SR_RA_vs_20_SR_LA_c2.0_cond1.bed$V2
diff_12b_SR_RA_vs_2231_AF_LA_c2.0_cond1.bed_length<- diff_12b_SR_RA_vs_2231_AF_LA_c2.0_cond1.bed$V3 - diff_12b_SR_RA_vs_2231_AF_LA_c2.0_cond1.bed$V2
diff_12b_SR_RA_vs_25_AF_LA_c2.0_cond1.bed_length<- diff_12b_SR_RA_vs_25_AF_LA_c2.0_cond1.bed$V3 - diff_12b_SR_RA_vs_25_AF_LA_c2.0_cond1.bed$V2
diff_12b_SR_RA_vs_27_SR_LA_c2.0_cond1.bed_length<- diff_12b_SR_RA_vs_27_SR_LA_c2.0_cond1.bed$V3 - diff_12b_SR_RA_vs_27_SR_LA_c2.0_cond1.bed$V2
diff_12b_SR_RA_vs_28_AF_LA_c2.0_cond1.bed_length<- diff_12b_SR_RA_vs_28_AF_LA_c2.0_cond1.bed$V3 - diff_12b_SR_RA_vs_28_AF_LA_c2.0_cond1.bed$V2
diff_12b_SR_RA_vs_31_AF_LA_c2.0_cond1.bed_length<- diff_12b_SR_RA_vs_31_AF_LA_c2.0_cond1.bed$V3 - diff_12b_SR_RA_vs_31_AF_LA_c2.0_cond1.bed$V2
diff_12b_SR_RA_vs_313_AF_LA_c2.0_cond1.bed_length<- diff_12b_SR_RA_vs_313_AF_LA_c2.0_cond1.bed$V3 - diff_12b_SR_RA_vs_313_AF_LA_c2.0_cond1.bed$V2
diff_12b_SR_RA_vs_45_SR_LA_c2.0_cond1.bed_length<- diff_12b_SR_RA_vs_45_SR_LA_c2.0_cond1.bed$V3 - diff_12b_SR_RA_vs_45_SR_LA_c2.0_cond1.bed$V2
diff_12b_SR_RA_vs_47_AF_LA_c2.0_cond1.bed_length<- diff_12b_SR_RA_vs_47_AF_LA_c2.0_cond1.bed$V3 - diff_12b_SR_RA_vs_47_AF_LA_c2.0_cond1.bed$V2
diff_18_AF_RA_vs_12b_SR_RA_c2.0_cond2.bed_length<- diff_18_AF_RA_vs_12b_SR_RA_c2.0_cond2.bed$V3 - diff_18_AF_RA_vs_12b_SR_RA_c2.0_cond2.bed$V2
diff_18_AF_RA_vs_20_SR_RA_c2.0_cond2.bed_length<- diff_18_AF_RA_vs_20_SR_RA_c2.0_cond2.bed$V3 - diff_18_AF_RA_vs_20_SR_RA_c2.0_cond2.bed$V2
diff_18_AF_RA_vs_27_SR_RA_c2.0_cond2.bed_length<- diff_18_AF_RA_vs_27_SR_RA_c2.0_cond2.bed$V3 - diff_18_AF_RA_vs_27_SR_RA_c2.0_cond2.bed$V2
diff_18_AF_RA_vs_33_SR_RA_c2.0_cond2.bed_length<- diff_18_AF_RA_vs_33_SR_RA_c2.0_cond2.bed$V3 - diff_18_AF_RA_vs_33_SR_RA_c2.0_cond2.bed$V2
diff_20_SR_RA_vs_12b_SR_LA_c2.0_cond1.bed_length<- diff_20_SR_RA_vs_12b_SR_LA_c2.0_cond1.bed$V3 - diff_20_SR_RA_vs_12b_SR_LA_c2.0_cond1.bed$V2
diff_20_SR_RA_vs_20_SR_LA_c2.0_cond1.bed_length<- diff_20_SR_RA_vs_20_SR_LA_c2.0_cond1.bed$V3 - diff_20_SR_RA_vs_20_SR_LA_c2.0_cond1.bed$V2
diff_20_SR_RA_vs_2231_AF_LA_c2.0_cond1.bed_length<- diff_20_SR_RA_vs_2231_AF_LA_c2.0_cond1.bed$V3 - diff_20_SR_RA_vs_2231_AF_LA_c2.0_cond1.bed$V2
diff_20_SR_RA_vs_25_AF_LA_c2.0_cond1.bed_length<- diff_20_SR_RA_vs_25_AF_LA_c2.0_cond1.bed$V3 - diff_20_SR_RA_vs_25_AF_LA_c2.0_cond1.bed$V2
diff_20_SR_RA_vs_27_SR_LA_c2.0_cond1.bed_length<- diff_20_SR_RA_vs_27_SR_LA_c2.0_cond1.bed$V3 - diff_20_SR_RA_vs_27_SR_LA_c2.0_cond1.bed$V2
diff_20_SR_RA_vs_28_AF_LA_c2.0_cond1.bed_length<- diff_20_SR_RA_vs_28_AF_LA_c2.0_cond1.bed$V3 - diff_20_SR_RA_vs_28_AF_LA_c2.0_cond1.bed$V2
diff_20_SR_RA_vs_31_AF_LA_c2.0_cond1.bed_length<- diff_20_SR_RA_vs_31_AF_LA_c2.0_cond1.bed$V3 - diff_20_SR_RA_vs_31_AF_LA_c2.0_cond1.bed$V2
diff_20_SR_RA_vs_313_AF_LA_c2.0_cond1.bed_length<- diff_20_SR_RA_vs_313_AF_LA_c2.0_cond1.bed$V3 - diff_20_SR_RA_vs_313_AF_LA_c2.0_cond1.bed$V2
diff_20_SR_RA_vs_45_SR_LA_c2.0_cond1.bed_length<- diff_20_SR_RA_vs_45_SR_LA_c2.0_cond1.bed$V3 - diff_20_SR_RA_vs_45_SR_LA_c2.0_cond1.bed$V2
diff_20_SR_RA_vs_47_AF_LA_c2.0_cond1.bed_length<- diff_20_SR_RA_vs_47_AF_LA_c2.0_cond1.bed$V3 - diff_20_SR_RA_vs_47_AF_LA_c2.0_cond1.bed$V2
diff_2231_AF_RA_vs_12b_SR_RA_c2.0_cond2.bed_length<- diff_2231_AF_RA_vs_12b_SR_RA_c2.0_cond2.bed$V3 - diff_2231_AF_RA_vs_12b_SR_RA_c2.0_cond2.bed$V2
diff_2231_AF_RA_vs_20_SR_RA_c2.0_cond2.bed_length<- diff_2231_AF_RA_vs_20_SR_RA_c2.0_cond2.bed$V3 - diff_2231_AF_RA_vs_20_SR_RA_c2.0_cond2.bed$V2
diff_2231_AF_RA_vs_27_SR_RA_c2.0_cond2.bed_length<- diff_2231_AF_RA_vs_27_SR_RA_c2.0_cond2.bed$V3 - diff_2231_AF_RA_vs_27_SR_RA_c2.0_cond2.bed$V2
diff_2231_AF_RA_vs_33_SR_RA_c2.0_cond2.bed_length<- diff_2231_AF_RA_vs_33_SR_RA_c2.0_cond2.bed$V3 - diff_2231_AF_RA_vs_33_SR_RA_c2.0_cond2.bed$V2
diff_25_AF_RA_vs_12b_SR_RA_c2.0_cond2.bed_length<- diff_25_AF_RA_vs_12b_SR_RA_c2.0_cond2.bed$V3 - diff_25_AF_RA_vs_12b_SR_RA_c2.0_cond2.bed$V2
diff_25_AF_RA_vs_20_SR_RA_c2.0_cond2.bed_length<- diff_25_AF_RA_vs_20_SR_RA_c2.0_cond2.bed$V3 - diff_25_AF_RA_vs_20_SR_RA_c2.0_cond2.bed$V2
diff_25_AF_RA_vs_27_SR_RA_c2.0_cond2.bed_length<- diff_25_AF_RA_vs_27_SR_RA_c2.0_cond2.bed$V3 - diff_25_AF_RA_vs_27_SR_RA_c2.0_cond2.bed$V2
diff_25_AF_RA_vs_33_SR_RA_c2.0_cond2.bed_length<- diff_25_AF_RA_vs_33_SR_RA_c2.0_cond2.bed$V3 - diff_25_AF_RA_vs_33_SR_RA_c2.0_cond2.bed$V2
diff_27_SR_RA_vs_12b_SR_LA_c2.0_cond1.bed_length<- diff_27_SR_RA_vs_12b_SR_LA_c2.0_cond1.bed$V3 - diff_27_SR_RA_vs_12b_SR_LA_c2.0_cond1.bed$V2
diff_27_SR_RA_vs_20_SR_LA_c2.0_cond1.bed_length<- diff_27_SR_RA_vs_20_SR_LA_c2.0_cond1.bed$V3 - diff_27_SR_RA_vs_20_SR_LA_c2.0_cond1.bed$V2
diff_27_SR_RA_vs_2231_AF_LA_c2.0_cond1.bed_length<- diff_27_SR_RA_vs_2231_AF_LA_c2.0_cond1.bed$V3 - diff_27_SR_RA_vs_2231_AF_LA_c2.0_cond1.bed$V2
diff_27_SR_RA_vs_25_AF_LA_c2.0_cond1.bed_length<- diff_27_SR_RA_vs_25_AF_LA_c2.0_cond1.bed$V3 - diff_27_SR_RA_vs_25_AF_LA_c2.0_cond1.bed$V2
diff_27_SR_RA_vs_27_SR_LA_c2.0_cond1.bed_length<- diff_27_SR_RA_vs_27_SR_LA_c2.0_cond1.bed$V3 - diff_27_SR_RA_vs_27_SR_LA_c2.0_cond1.bed$V2
diff_27_SR_RA_vs_28_AF_LA_c2.0_cond1.bed_length<- diff_27_SR_RA_vs_28_AF_LA_c2.0_cond1.bed$V3 - diff_27_SR_RA_vs_28_AF_LA_c2.0_cond1.bed$V2
diff_27_SR_RA_vs_31_AF_LA_c2.0_cond1.bed_length<- diff_27_SR_RA_vs_31_AF_LA_c2.0_cond1.bed$V3 - diff_27_SR_RA_vs_31_AF_LA_c2.0_cond1.bed$V2
diff_27_SR_RA_vs_313_AF_LA_c2.0_cond1.bed_length<- diff_27_SR_RA_vs_313_AF_LA_c2.0_cond1.bed$V3 - diff_27_SR_RA_vs_313_AF_LA_c2.0_cond1.bed$V2
diff_27_SR_RA_vs_45_SR_LA_c2.0_cond1.bed_length<- diff_27_SR_RA_vs_45_SR_LA_c2.0_cond1.bed$V3 - diff_27_SR_RA_vs_45_SR_LA_c2.0_cond1.bed$V2
diff_27_SR_RA_vs_47_AF_LA_c2.0_cond1.bed_length<- diff_27_SR_RA_vs_47_AF_LA_c2.0_cond1.bed$V3 - diff_27_SR_RA_vs_47_AF_LA_c2.0_cond1.bed$V2
diff_28_AF_RA_vs_12b_SR_RA_c2.0_cond2.bed_length<- diff_28_AF_RA_vs_12b_SR_RA_c2.0_cond2.bed$V3 - diff_28_AF_RA_vs_12b_SR_RA_c2.0_cond2.bed$V2
diff_28_AF_RA_vs_20_SR_RA_c2.0_cond2.bed_length<- diff_28_AF_RA_vs_20_SR_RA_c2.0_cond2.bed$V3 - diff_28_AF_RA_vs_20_SR_RA_c2.0_cond2.bed$V2
diff_28_AF_RA_vs_27_SR_RA_c2.0_cond2.bed_length<- diff_28_AF_RA_vs_27_SR_RA_c2.0_cond2.bed$V3 - diff_28_AF_RA_vs_27_SR_RA_c2.0_cond2.bed$V2
diff_28_AF_RA_vs_33_SR_RA_c2.0_cond2.bed_length<- diff_28_AF_RA_vs_33_SR_RA_c2.0_cond2.bed$V3 - diff_28_AF_RA_vs_33_SR_RA_c2.0_cond2.bed$V2
diff_31_AF_RA_vs_12b_SR_RA_c2.0_cond2.bed_length<- diff_31_AF_RA_vs_12b_SR_RA_c2.0_cond2.bed$V3 - diff_31_AF_RA_vs_12b_SR_RA_c2.0_cond2.bed$V2
diff_31_AF_RA_vs_20_SR_RA_c2.0_cond2.bed_length<- diff_31_AF_RA_vs_20_SR_RA_c2.0_cond2.bed$V3 - diff_31_AF_RA_vs_20_SR_RA_c2.0_cond2.bed$V2
diff_31_AF_RA_vs_27_SR_RA_c2.0_cond2.bed_length<- diff_31_AF_RA_vs_27_SR_RA_c2.0_cond2.bed$V3 - diff_31_AF_RA_vs_27_SR_RA_c2.0_cond2.bed$V2
diff_31_AF_RA_vs_33_SR_RA_c2.0_cond2.bed_length<- diff_31_AF_RA_vs_33_SR_RA_c2.0_cond2.bed$V3 - diff_31_AF_RA_vs_33_SR_RA_c2.0_cond2.bed$V2
diff_313_AF_RA_vs_12b_SR_RA_c2.0_cond2.bed_length<- diff_313_AF_RA_vs_12b_SR_RA_c2.0_cond2.bed$V3 - diff_313_AF_RA_vs_12b_SR_RA_c2.0_cond2.bed$V2
diff_313_AF_RA_vs_20_SR_RA_c2.0_cond2.bed_length<- diff_313_AF_RA_vs_20_SR_RA_c2.0_cond2.bed$V3 - diff_313_AF_RA_vs_20_SR_RA_c2.0_cond2.bed$V2
diff_313_AF_RA_vs_27_SR_RA_c2.0_cond2.bed_length<- diff_313_AF_RA_vs_27_SR_RA_c2.0_cond2.bed$V3 - diff_313_AF_RA_vs_27_SR_RA_c2.0_cond2.bed$V2
diff_313_AF_RA_vs_33_SR_RA_c2.0_cond2.bed_length<- diff_313_AF_RA_vs_33_SR_RA_c2.0_cond2.bed$V3 - diff_313_AF_RA_vs_33_SR_RA_c2.0_cond2.bed$V2
diff_33_SR_RA_vs_12b_SR_LA_c2.0_cond1.bed_length<- diff_33_SR_RA_vs_12b_SR_LA_c2.0_cond1.bed$V3 - diff_33_SR_RA_vs_12b_SR_LA_c2.0_cond1.bed$V2
diff_33_SR_RA_vs_20_SR_LA_c2.0_cond1.bed_length<- diff_33_SR_RA_vs_20_SR_LA_c2.0_cond1.bed$V3 - diff_33_SR_RA_vs_20_SR_LA_c2.0_cond1.bed$V2
diff_33_SR_RA_vs_2231_AF_LA_c2.0_cond1.bed_length<- diff_33_SR_RA_vs_2231_AF_LA_c2.0_cond1.bed$V3 - diff_33_SR_RA_vs_2231_AF_LA_c2.0_cond1.bed$V2
diff_33_SR_RA_vs_25_AF_LA_c2.0_cond1.bed_length<- diff_33_SR_RA_vs_25_AF_LA_c2.0_cond1.bed$V3 - diff_33_SR_RA_vs_25_AF_LA_c2.0_cond1.bed$V2
diff_33_SR_RA_vs_27_SR_LA_c2.0_cond1.bed_length<- diff_33_SR_RA_vs_27_SR_LA_c2.0_cond1.bed$V3 - diff_33_SR_RA_vs_27_SR_LA_c2.0_cond1.bed$V2
diff_33_SR_RA_vs_28_AF_LA_c2.0_cond1.bed_length<- diff_33_SR_RA_vs_28_AF_LA_c2.0_cond1.bed$V3 - diff_33_SR_RA_vs_28_AF_LA_c2.0_cond1.bed$V2
diff_33_SR_RA_vs_31_AF_LA_c2.0_cond1.bed_length<- diff_33_SR_RA_vs_31_AF_LA_c2.0_cond1.bed$V3 - diff_33_SR_RA_vs_31_AF_LA_c2.0_cond1.bed$V2
diff_33_SR_RA_vs_313_AF_LA_c2.0_cond1.bed_length<- diff_33_SR_RA_vs_313_AF_LA_c2.0_cond1.bed$V3 - diff_33_SR_RA_vs_313_AF_LA_c2.0_cond1.bed$V2
diff_33_SR_RA_vs_45_SR_LA_c2.0_cond1.bed_length<- diff_33_SR_RA_vs_45_SR_LA_c2.0_cond1.bed$V3 - diff_33_SR_RA_vs_45_SR_LA_c2.0_cond1.bed$V2
diff_33_SR_RA_vs_47_AF_LA_c2.0_cond1.bed_length<- diff_33_SR_RA_vs_47_AF_LA_c2.0_cond1.bed$V3 - diff_33_SR_RA_vs_47_AF_LA_c2.0_cond1.bed$V2


diff_12b_SR_RA_vs_12b_SR_LA_c2.0_cond1.bed_length <- melt(diff_12b_SR_RA_vs_12b_SR_LA_c2.0_cond1.bed_length)
diff_12b_SR_RA_vs_20_SR_LA_c2.0_cond1.bed_length <- melt(diff_12b_SR_RA_vs_20_SR_LA_c2.0_cond1.bed_length)
diff_12b_SR_RA_vs_2231_AF_LA_c2.0_cond1.bed_length  <- melt(diff_12b_SR_RA_vs_2231_AF_LA_c2.0_cond1.bed_length)
diff_12b_SR_RA_vs_25_AF_LA_c2.0_cond1.bed_length  <- melt(diff_12b_SR_RA_vs_25_AF_LA_c2.0_cond1.bed_length)
diff_12b_SR_RA_vs_27_SR_LA_c2.0_cond1.bed_length  <- melt(diff_12b_SR_RA_vs_27_SR_LA_c2.0_cond1.bed_length)
diff_12b_SR_RA_vs_28_AF_LA_c2.0_cond1.bed_length  <- melt(diff_12b_SR_RA_vs_28_AF_LA_c2.0_cond1.bed_length) 
diff_12b_SR_RA_vs_31_AF_LA_c2.0_cond1.bed_length  <- melt(diff_12b_SR_RA_vs_31_AF_LA_c2.0_cond1.bed_length)
diff_12b_SR_RA_vs_313_AF_LA_c2.0_cond1.bed_length  <- melt(diff_12b_SR_RA_vs_313_AF_LA_c2.0_cond1.bed_length)
diff_12b_SR_RA_vs_45_SR_LA_c2.0_cond1.bed_length  <- melt(diff_12b_SR_RA_vs_45_SR_LA_c2.0_cond1.bed_length)
diff_12b_SR_RA_vs_47_AF_LA_c2.0_cond1.bed_length  <- melt(diff_12b_SR_RA_vs_47_AF_LA_c2.0_cond1.bed_length)
diff_18_AF_RA_vs_12b_SR_RA_c2.0_cond2.bed_length  <- melt(diff_18_AF_RA_vs_12b_SR_RA_c2.0_cond2.bed_length)
diff_18_AF_RA_vs_20_SR_RA_c2.0_cond2.bed_length  <- melt(diff_18_AF_RA_vs_20_SR_RA_c2.0_cond2.bed_length)
diff_18_AF_RA_vs_27_SR_RA_c2.0_cond2.bed_length  <- melt(diff_18_AF_RA_vs_27_SR_RA_c2.0_cond2.bed_length)
diff_18_AF_RA_vs_33_SR_RA_c2.0_cond2.bed_length  <- melt(diff_18_AF_RA_vs_33_SR_RA_c2.0_cond2.bed_length)
diff_20_SR_RA_vs_12b_SR_LA_c2.0_cond1.bed_length  <- melt(diff_20_SR_RA_vs_12b_SR_LA_c2.0_cond1.bed_length)
diff_20_SR_RA_vs_20_SR_LA_c2.0_cond1.bed_length  <- melt(diff_20_SR_RA_vs_20_SR_LA_c2.0_cond1.bed_length)
diff_20_SR_RA_vs_2231_AF_LA_c2.0_cond1.bed_length  <- melt(diff_20_SR_RA_vs_2231_AF_LA_c2.0_cond1.bed_length)
diff_20_SR_RA_vs_25_AF_LA_c2.0_cond1.bed_length  <- melt(diff_20_SR_RA_vs_25_AF_LA_c2.0_cond1.bed_length)
diff_20_SR_RA_vs_27_SR_LA_c2.0_cond1.bed_length  <- melt(diff_20_SR_RA_vs_27_SR_LA_c2.0_cond1.bed_length)
diff_20_SR_RA_vs_28_AF_LA_c2.0_cond1.bed_length  <- melt(diff_20_SR_RA_vs_28_AF_LA_c2.0_cond1.bed_length)
diff_20_SR_RA_vs_31_AF_LA_c2.0_cond1.bed_length  <- melt(diff_20_SR_RA_vs_31_AF_LA_c2.0_cond1.bed_length)
diff_20_SR_RA_vs_313_AF_LA_c2.0_cond1.bed_length  <- melt(diff_20_SR_RA_vs_313_AF_LA_c2.0_cond1.bed_length)
diff_20_SR_RA_vs_45_SR_LA_c2.0_cond1.bed_length  <- melt(diff_20_SR_RA_vs_45_SR_LA_c2.0_cond1.bed_length)
diff_20_SR_RA_vs_47_AF_LA_c2.0_cond1.bed_length  <- melt(diff_20_SR_RA_vs_47_AF_LA_c2.0_cond1.bed_length)
diff_2231_AF_RA_vs_12b_SR_RA_c2.0_cond2.bed_length  <- melt(diff_2231_AF_RA_vs_12b_SR_RA_c2.0_cond2.bed_length)
diff_2231_AF_RA_vs_20_SR_RA_c2.0_cond2.bed_length  <- melt(diff_2231_AF_RA_vs_20_SR_RA_c2.0_cond2.bed_length)
diff_2231_AF_RA_vs_27_SR_RA_c2.0_cond2.bed_length  <- melt(diff_2231_AF_RA_vs_27_SR_RA_c2.0_cond2.bed_length)
diff_2231_AF_RA_vs_33_SR_RA_c2.0_cond2.bed_length  <- melt(diff_2231_AF_RA_vs_33_SR_RA_c2.0_cond2.bed_length)
diff_25_AF_RA_vs_12b_SR_RA_c2.0_cond2.bed_length  <- melt(diff_25_AF_RA_vs_12b_SR_RA_c2.0_cond2.bed_length)
diff_25_AF_RA_vs_20_SR_RA_c2.0_cond2.bed_length  <- melt(diff_25_AF_RA_vs_20_SR_RA_c2.0_cond2.bed_length)
diff_25_AF_RA_vs_27_SR_RA_c2.0_cond2.bed_length  <- melt(diff_25_AF_RA_vs_27_SR_RA_c2.0_cond2.bed_length)
diff_25_AF_RA_vs_33_SR_RA_c2.0_cond2.bed_length <- melt(diff_25_AF_RA_vs_33_SR_RA_c2.0_cond2.bed_length)
diff_27_SR_RA_vs_12b_SR_LA_c2.0_cond1.bed_length <- melt(diff_27_SR_RA_vs_12b_SR_LA_c2.0_cond1.bed_length)
diff_27_SR_RA_vs_20_SR_LA_c2.0_cond1.bed_length  <- melt(diff_27_SR_RA_vs_20_SR_LA_c2.0_cond1.bed_length)
diff_27_SR_RA_vs_2231_AF_LA_c2.0_cond1.bed_length  <- melt(diff_27_SR_RA_vs_2231_AF_LA_c2.0_cond1.bed_length)
diff_27_SR_RA_vs_25_AF_LA_c2.0_cond1.bed_length  <- melt(diff_27_SR_RA_vs_25_AF_LA_c2.0_cond1.bed_length)
diff_27_SR_RA_vs_27_SR_LA_c2.0_cond1.bed_length  <- melt(diff_27_SR_RA_vs_27_SR_LA_c2.0_cond1.bed_length)
diff_27_SR_RA_vs_28_AF_LA_c2.0_cond1.bed_length  <- melt(diff_27_SR_RA_vs_28_AF_LA_c2.0_cond1.bed_length)
diff_27_SR_RA_vs_31_AF_LA_c2.0_cond1.bed_length  <- melt(diff_27_SR_RA_vs_31_AF_LA_c2.0_cond1.bed_length)
diff_27_SR_RA_vs_313_AF_LA_c2.0_cond1.bed_length  <- melt(diff_27_SR_RA_vs_313_AF_LA_c2.0_cond1.bed_length)
diff_27_SR_RA_vs_45_SR_LA_c2.0_cond1.bed_length  <- melt(diff_27_SR_RA_vs_45_SR_LA_c2.0_cond1.bed_length)
diff_27_SR_RA_vs_47_AF_LA_c2.0_cond1.bed_length  <- melt(diff_27_SR_RA_vs_47_AF_LA_c2.0_cond1.bed_length)
diff_28_AF_RA_vs_12b_SR_RA_c2.0_cond2.bed_length  <- melt(diff_28_AF_RA_vs_12b_SR_RA_c2.0_cond2.bed_length)
diff_28_AF_RA_vs_20_SR_RA_c2.0_cond2.bed_length  <- melt(diff_28_AF_RA_vs_20_SR_RA_c2.0_cond2.bed_length)
diff_28_AF_RA_vs_27_SR_RA_c2.0_cond2.bed_length  <- melt(diff_28_AF_RA_vs_27_SR_RA_c2.0_cond2.bed_length)
diff_28_AF_RA_vs_33_SR_RA_c2.0_cond2.bed_length  <- melt(diff_28_AF_RA_vs_33_SR_RA_c2.0_cond2.bed_length)
diff_31_AF_RA_vs_12b_SR_RA_c2.0_cond2.bed_length  <- melt(diff_31_AF_RA_vs_12b_SR_RA_c2.0_cond2.bed_length)
diff_31_AF_RA_vs_20_SR_RA_c2.0_cond2.bed_length  <- melt(diff_31_AF_RA_vs_20_SR_RA_c2.0_cond2.bed_length)
diff_31_AF_RA_vs_27_SR_RA_c2.0_cond2.bed_length  <- melt(diff_31_AF_RA_vs_27_SR_RA_c2.0_cond2.bed_length)
diff_31_AF_RA_vs_33_SR_RA_c2.0_cond2.bed_length  <- melt(diff_31_AF_RA_vs_33_SR_RA_c2.0_cond2.bed_length)
diff_313_AF_RA_vs_12b_SR_RA_c2.0_cond2.bed_length  <- melt(diff_313_AF_RA_vs_12b_SR_RA_c2.0_cond2.bed_length)
diff_313_AF_RA_vs_20_SR_RA_c2.0_cond2.bed_length  <- melt(diff_313_AF_RA_vs_20_SR_RA_c2.0_cond2.bed_length)
diff_313_AF_RA_vs_27_SR_RA_c2.0_cond2.bed_length  <- melt(diff_313_AF_RA_vs_27_SR_RA_c2.0_cond2.bed_length)
diff_313_AF_RA_vs_33_SR_RA_c2.0_cond2.bed_length  <- melt(diff_313_AF_RA_vs_33_SR_RA_c2.0_cond2.bed_length)
diff_33_SR_RA_vs_12b_SR_LA_c2.0_cond1.bed_length  <- melt(diff_33_SR_RA_vs_12b_SR_LA_c2.0_cond1.bed_length)
diff_33_SR_RA_vs_20_SR_LA_c2.0_cond1.bed_length  <- melt(diff_33_SR_RA_vs_20_SR_LA_c2.0_cond1.bed_length)
diff_33_SR_RA_vs_2231_AF_LA_c2.0_cond1.bed_length  <- melt(diff_33_SR_RA_vs_2231_AF_LA_c2.0_cond1.bed_length)
diff_33_SR_RA_vs_25_AF_LA_c2.0_cond1.bed_length  <- melt(diff_33_SR_RA_vs_25_AF_LA_c2.0_cond1.bed_length)
diff_33_SR_RA_vs_27_SR_LA_c2.0_cond1.bed_length  <- melt(diff_33_SR_RA_vs_27_SR_LA_c2.0_cond1.bed_length)
diff_33_SR_RA_vs_28_AF_LA_c2.0_cond1.bed_length <- melt(diff_33_SR_RA_vs_28_AF_LA_c2.0_cond1.bed_length)
diff_33_SR_RA_vs_31_AF_LA_c2.0_cond1.bed_length  <- melt(diff_33_SR_RA_vs_31_AF_LA_c2.0_cond1.bed_length)
diff_33_SR_RA_vs_313_AF_LA_c2.0_cond1.bed_length  <- melt(diff_33_SR_RA_vs_313_AF_LA_c2.0_cond1.bed_length)
diff_33_SR_RA_vs_45_SR_LA_c2.0_cond1.bed_length  <- melt(diff_33_SR_RA_vs_45_SR_LA_c2.0_cond1.bed_length)
diff_33_SR_RA_vs_47_AF_LA_c2.0_cond1.bed_length  <- melt(diff_33_SR_RA_vs_47_AF_LA_c2.0_cond1.bed_length)


diff_12b_SR_RA_vs_12b_SR_LA_c2.0_cond1.bed_length$group <- "SR_RA"
diff_12b_SR_RA_vs_20_SR_LA_c2.0_cond1.bed_length$group <- "SR_RA"
diff_12b_SR_RA_vs_2231_AF_LA_c2.0_cond1.bed_length$group <- "SR_RA"
diff_12b_SR_RA_vs_25_AF_LA_c2.0_cond1.bed_length$group <- "SR_RA"
diff_12b_SR_RA_vs_27_SR_LA_c2.0_cond1.bed_length$group <- "SR_RA" 
diff_12b_SR_RA_vs_28_AF_LA_c2.0_cond1.bed_length$group <- "SR_RA"
diff_12b_SR_RA_vs_31_AF_LA_c2.0_cond1.bed_length$group <- "SR_RA" 
diff_12b_SR_RA_vs_313_AF_LA_c2.0_cond1.bed_length$group <- "SR_RA"
diff_12b_SR_RA_vs_45_SR_LA_c2.0_cond1.bed_length$group <- "SR_RA" 
diff_12b_SR_RA_vs_47_AF_LA_c2.0_cond1.bed_length$group <- "SR_RA"
diff_18_AF_RA_vs_12b_SR_RA_c2.0_cond2.bed_length$group <- "SR_RA"
diff_18_AF_RA_vs_20_SR_RA_c2.0_cond2.bed_length$group <- "SR_RA"
diff_18_AF_RA_vs_27_SR_RA_c2.0_cond2.bed_length$group <- "SR_RA"
diff_18_AF_RA_vs_33_SR_RA_c2.0_cond2.bed_length$group <- "SR_RA"
diff_20_SR_RA_vs_12b_SR_LA_c2.0_cond1.bed_length$group <- "SR_RA"
diff_20_SR_RA_vs_20_SR_LA_c2.0_cond1.bed_length$group <- "SR_RA"
diff_20_SR_RA_vs_2231_AF_LA_c2.0_cond1.bed_length$group <- "SR_RA"
diff_20_SR_RA_vs_25_AF_LA_c2.0_cond1.bed_length$group <- "SR_RA"
diff_20_SR_RA_vs_27_SR_LA_c2.0_cond1.bed_length$group <- "SR_RA"
diff_20_SR_RA_vs_28_AF_LA_c2.0_cond1.bed_length$group <- "SR_RA"
diff_20_SR_RA_vs_31_AF_LA_c2.0_cond1.bed_length$group <- "SR_RA" 
diff_20_SR_RA_vs_313_AF_LA_c2.0_cond1.bed_length$group <- "SR_RA" 
diff_20_SR_RA_vs_45_SR_LA_c2.0_cond1.bed_length$group <- "SR_RA" 
diff_20_SR_RA_vs_47_AF_LA_c2.0_cond1.bed_length$group <- "SR_RA" 
diff_2231_AF_RA_vs_12b_SR_RA_c2.0_cond2.bed_length$group <- "SR_RA" 
diff_2231_AF_RA_vs_20_SR_RA_c2.0_cond2.bed_length$group <- "SR_RA" 
diff_2231_AF_RA_vs_27_SR_RA_c2.0_cond2.bed_length$group <- "SR_RA" 
diff_2231_AF_RA_vs_33_SR_RA_c2.0_cond2.bed_length$group <- "SR_RA" 
diff_25_AF_RA_vs_12b_SR_RA_c2.0_cond2.bed_length$group <- "SR_RA" 
diff_25_AF_RA_vs_20_SR_RA_c2.0_cond2.bed_length$group <- "SR_RA"
diff_25_AF_RA_vs_27_SR_RA_c2.0_cond2.bed_length$group <- "SR_RA" 
diff_25_AF_RA_vs_33_SR_RA_c2.0_cond2.bed_length$group <- "SR_RA" 
diff_27_SR_RA_vs_12b_SR_LA_c2.0_cond1.bed_length$group <- "SR_RA" 
diff_27_SR_RA_vs_20_SR_LA_c2.0_cond1.bed_length$group <- "SR_RA" 
diff_27_SR_RA_vs_2231_AF_LA_c2.0_cond1.bed_length$group <- "SR_RA" 
diff_27_SR_RA_vs_25_AF_LA_c2.0_cond1.bed_length$group <- "SR_RA"
diff_27_SR_RA_vs_27_SR_LA_c2.0_cond1.bed_length$group <- "SR_RA" 
diff_27_SR_RA_vs_28_AF_LA_c2.0_cond1.bed_length$group <- "SR_RA" 
diff_27_SR_RA_vs_31_AF_LA_c2.0_cond1.bed_length$group <- "SR_RA" 
diff_27_SR_RA_vs_313_AF_LA_c2.0_cond1.bed_length$group <- "SR_RA" 
diff_27_SR_RA_vs_45_SR_LA_c2.0_cond1.bed_length$group <- "SR_RA" 
diff_27_SR_RA_vs_47_AF_LA_c2.0_cond1.bed_length$group <- "SR_RA" 
diff_28_AF_RA_vs_12b_SR_RA_c2.0_cond2.bed_length$group <- "SR_RA" 
diff_28_AF_RA_vs_20_SR_RA_c2.0_cond2.bed_length$group <- "SR_RA" 
diff_28_AF_RA_vs_27_SR_RA_c2.0_cond2.bed_length$group <- "SR_RA" 
diff_28_AF_RA_vs_33_SR_RA_c2.0_cond2.bed_length$group <- "SR_RA" 
diff_31_AF_RA_vs_12b_SR_RA_c2.0_cond2.bed_length$group <- "SR_RA" 
diff_31_AF_RA_vs_20_SR_RA_c2.0_cond2.bed_length$group <- "SR_RA" 
diff_31_AF_RA_vs_27_SR_RA_c2.0_cond2.bed_length$group <- "SR_RA" 
diff_31_AF_RA_vs_33_SR_RA_c2.0_cond2.bed_length$group <- "SR_RA" 
diff_313_AF_RA_vs_12b_SR_RA_c2.0_cond2.bed_length$group <- "SR_RA" 
diff_313_AF_RA_vs_20_SR_RA_c2.0_cond2.bed_length$group <- "SR_RA" 
diff_313_AF_RA_vs_27_SR_RA_c2.0_cond2.bed_length$group <- "SR_RA"
diff_313_AF_RA_vs_33_SR_RA_c2.0_cond2.bed_length$group <- "SR_RA" 
diff_33_SR_RA_vs_12b_SR_LA_c2.0_cond1.bed_length$group <- "SR_RA" 
diff_33_SR_RA_vs_20_SR_LA_c2.0_cond1.bed_length$group <- "SR_RA" 
diff_33_SR_RA_vs_2231_AF_LA_c2.0_cond1.bed_length$group <- "SR_RA" 
diff_33_SR_RA_vs_25_AF_LA_c2.0_cond1.bed_length$group <- "SR_RA" 
diff_33_SR_RA_vs_27_SR_LA_c2.0_cond1.bed_length$group <- "SR_RA" 
diff_33_SR_RA_vs_28_AF_LA_c2.0_cond1.bed_length$group <- "SR_RA"
diff_33_SR_RA_vs_31_AF_LA_c2.0_cond1.bed_length$group <- "SR_RA" 
diff_33_SR_RA_vs_313_AF_LA_c2.0_cond1.bed_length$group <- "SR_RA" 
diff_33_SR_RA_vs_45_SR_LA_c2.0_cond1.bed_length$group <- "SR_RA" 
diff_33_SR_RA_vs_47_AF_LA_c2.0_cond1.bed_length$group <- "SR_RA" 

diff_12b_SR_RA_vs_12b_SR_LA_c2.0_cond1.bed_length$sample <- "diff_12b_SR_RA_vs_12b_SR_LA_c2.0_cond1.bed" 
diff_12b_SR_RA_vs_20_SR_LA_c2.0_cond1.bed_length$sample <- "diff_12b_SR_RA_vs_20_SR_LA_c2.0_cond1.bed" 
diff_12b_SR_RA_vs_2231_AF_LA_c2.0_cond1.bed_length$sample <- "diff_12b_SR_RA_vs_2231_AF_LA_c2.0_cond1.bed"  
diff_12b_SR_RA_vs_25_AF_LA_c2.0_cond1.bed_length$sample <- "diff_12b_SR_RA_vs_25_AF_LA_c2.0_cond1.bed"  
diff_12b_SR_RA_vs_27_SR_LA_c2.0_cond1.bed_length$sample <- "diff_12b_SR_RA_vs_27_SR_LA_c2.0_cond1.bed"  
diff_12b_SR_RA_vs_28_AF_LA_c2.0_cond1.bed_length$sample <- "diff_12b_SR_RA_vs_28_AF_LA_c2.0_cond1.bed"  
diff_12b_SR_RA_vs_31_AF_LA_c2.0_cond1.bed_length$sample <- "diff_12b_SR_RA_vs_31_AF_LA_c2.0_cond1.bed"  
diff_12b_SR_RA_vs_313_AF_LA_c2.0_cond1.bed_length$sample <- "diff_12b_SR_RA_vs_313_AF_LA_c2.0_cond1.bed"  
diff_12b_SR_RA_vs_45_SR_LA_c2.0_cond1.bed_length$sample <- "diff_12b_SR_RA_vs_45_SR_LA_c2.0_cond1.bed"  
diff_12b_SR_RA_vs_47_AF_LA_c2.0_cond1.bed_length$sample <- "diff_12b_SR_RA_vs_47_AF_LA_c2.0_cond1.bed"  
diff_18_AF_RA_vs_12b_SR_RA_c2.0_cond2.bed_length$sample <- "diff_18_AF_RA_vs_12b_SR_RA_c2.0_cond2.bed"  
diff_18_AF_RA_vs_20_SR_RA_c2.0_cond2.bed_length$sample <- "diff_18_AF_RA_vs_20_SR_RA_c2.0_cond2.bed"  
diff_18_AF_RA_vs_27_SR_RA_c2.0_cond2.bed_length$sample <- "diff_18_AF_RA_vs_27_SR_RA_c2.0_cond2.bed"  
diff_18_AF_RA_vs_33_SR_RA_c2.0_cond2.bed_length$sample <- "diff_18_AF_RA_vs_33_SR_RA_c2.0_cond2.bed"  
diff_20_SR_RA_vs_12b_SR_LA_c2.0_cond1.bed_length$sample <- "diff_20_SR_RA_vs_12b_SR_LA_c2.0_cond1.bed"  
diff_20_SR_RA_vs_20_SR_LA_c2.0_cond1.bed_length$sample <- "diff_20_SR_RA_vs_20_SR_LA_c2.0_cond1.bed"  
diff_20_SR_RA_vs_2231_AF_LA_c2.0_cond1.bed_length$sample <- "diff_20_SR_RA_vs_2231_AF_LA_c2.0_cond1.bed"  
diff_20_SR_RA_vs_25_AF_LA_c2.0_cond1.bed_length$sample <- "diff_20_SR_RA_vs_25_AF_LA_c2.0_cond1.bed"  
diff_20_SR_RA_vs_27_SR_LA_c2.0_cond1.bed_length$sample <- "diff_20_SR_RA_vs_27_SR_LA_c2.0_cond1.bed"  
diff_20_SR_RA_vs_28_AF_LA_c2.0_cond1.bed_length$sample <- "diff_20_SR_RA_vs_28_AF_LA_c2.0_cond1.bed"  
diff_20_SR_RA_vs_31_AF_LA_c2.0_cond1.bed_length$sample <- "diff_20_SR_RA_vs_31_AF_LA_c2.0_cond1.bed"  
diff_20_SR_RA_vs_313_AF_LA_c2.0_cond1.bed_length$sample <- "diff_20_SR_RA_vs_313_AF_LA_c2.0_cond1.bed"  
diff_20_SR_RA_vs_45_SR_LA_c2.0_cond1.bed_length$sample <- "diff_20_SR_RA_vs_45_SR_LA_c2.0_cond1.bed"  
diff_20_SR_RA_vs_47_AF_LA_c2.0_cond1.bed_length$sample <- "diff_20_SR_RA_vs_47_AF_LA_c2.0_cond1.bed"  
diff_2231_AF_RA_vs_12b_SR_RA_c2.0_cond2.bed_length$sample <- "diff_2231_AF_RA_vs_12b_SR_RA_c2.0_cond2.bed"  
diff_2231_AF_RA_vs_20_SR_RA_c2.0_cond2.bed_length$sample <- "diff_2231_AF_RA_vs_20_SR_RA_c2.0_cond2.bed"  
diff_2231_AF_RA_vs_27_SR_RA_c2.0_cond2.bed_length$sample <- "diff_2231_AF_RA_vs_27_SR_RA_c2.0_cond2.bed"  
diff_2231_AF_RA_vs_33_SR_RA_c2.0_cond2.bed_length$sample <- "diff_2231_AF_RA_vs_33_SR_RA_c2.0_cond2.bed"  
diff_25_AF_RA_vs_12b_SR_RA_c2.0_cond2.bed_length$sample <- "diff_25_AF_RA_vs_12b_SR_RA_c2.0_cond2.bed"  
diff_25_AF_RA_vs_20_SR_RA_c2.0_cond2.bed_length$sample <- "diff_25_AF_RA_vs_20_SR_RA_c2.0_cond2.bed"  
diff_25_AF_RA_vs_27_SR_RA_c2.0_cond2.bed_length$sample <- "diff_25_AF_RA_vs_27_SR_RA_c2.0_cond2.bed"  
diff_25_AF_RA_vs_33_SR_RA_c2.0_cond2.bed_length$sample <- "diff_25_AF_RA_vs_33_SR_RA_c2.0_cond2.bed"  
diff_27_SR_RA_vs_12b_SR_LA_c2.0_cond1.bed_length$sample <- "diff_27_SR_RA_vs_12b_SR_LA_c2.0_cond1.bed"  
diff_27_SR_RA_vs_20_SR_LA_c2.0_cond1.bed_length$sample <- "diff_27_SR_RA_vs_20_SR_LA_c2.0_cond1.bed" 
diff_27_SR_RA_vs_2231_AF_LA_c2.0_cond1.bed_length$sample <- "diff_27_SR_RA_vs_2231_AF_LA_c2.0_cond1.bed" 
diff_27_SR_RA_vs_25_AF_LA_c2.0_cond1.bed_length$sample <- "diff_27_SR_RA_vs_25_AF_LA_c2.0_cond1.bed"  
diff_27_SR_RA_vs_27_SR_LA_c2.0_cond1.bed_length$sample <- "diff_27_SR_RA_vs_27_SR_LA_c2.0_cond1.bed"  
diff_27_SR_RA_vs_28_AF_LA_c2.0_cond1.bed_length$sample <- "diff_27_SR_RA_vs_28_AF_LA_c2.0_cond1.bed"  
diff_27_SR_RA_vs_31_AF_LA_c2.0_cond1.bed_length$sample <- "diff_27_SR_RA_vs_31_AF_LA_c2.0_cond1.bed"  
diff_27_SR_RA_vs_313_AF_LA_c2.0_cond1.bed_length$sample <- "diff_27_SR_RA_vs_313_AF_LA_c2.0_cond1.bed"  
diff_27_SR_RA_vs_45_SR_LA_c2.0_cond1.bed_length$sample <- "diff_27_SR_RA_vs_45_SR_LA_c2.0_cond1.bed"  
diff_27_SR_RA_vs_47_AF_LA_c2.0_cond1.bed_length$sample <- "diff_27_SR_RA_vs_47_AF_LA_c2.0_cond1.bed"  
diff_28_AF_RA_vs_12b_SR_RA_c2.0_cond2.bed_length$sample <- "diff_28_AF_RA_vs_12b_SR_RA_c2.0_cond2.bed"  
diff_28_AF_RA_vs_20_SR_RA_c2.0_cond2.bed_length$sample <- "diff_28_AF_RA_vs_20_SR_RA_c2.0_cond2.bed"  
diff_28_AF_RA_vs_27_SR_RA_c2.0_cond2.bed_length$sample <- "diff_28_AF_RA_vs_27_SR_RA_c2.0_cond2.bed"  
diff_28_AF_RA_vs_33_SR_RA_c2.0_cond2.bed_length$sample <- "diff_28_AF_RA_vs_33_SR_RA_c2.0_cond2.bed"  
diff_31_AF_RA_vs_12b_SR_RA_c2.0_cond2.bed_length$sample <- "diff_31_AF_RA_vs_12b_SR_RA_c2.0_cond2.bed"  
diff_31_AF_RA_vs_20_SR_RA_c2.0_cond2.bed_length$sample <- "diff_31_AF_RA_vs_20_SR_RA_c2.0_cond2.bed"  
diff_31_AF_RA_vs_27_SR_RA_c2.0_cond2.bed_length$sample <- "diff_31_AF_RA_vs_27_SR_RA_c2.0_cond2.bed"  
diff_31_AF_RA_vs_33_SR_RA_c2.0_cond2.bed_length$sample <- "diff_31_AF_RA_vs_33_SR_RA_c2.0_cond2.bed"  
diff_313_AF_RA_vs_12b_SR_RA_c2.0_cond2.bed_length$sample <- "diff_313_AF_RA_vs_12b_SR_RA_c2.0_cond2.bed"  
diff_313_AF_RA_vs_20_SR_RA_c2.0_cond2.bed_length$sample <- "diff_313_AF_RA_vs_20_SR_RA_c2.0_cond2.bed"  
diff_313_AF_RA_vs_27_SR_RA_c2.0_cond2.bed_length$sample <- "diff_313_AF_RA_vs_27_SR_RA_c2.0_cond2.bed"  
diff_313_AF_RA_vs_33_SR_RA_c2.0_cond2.bed_length$sample <- "diff_313_AF_RA_vs_33_SR_RA_c2.0_cond2.bed"  
diff_33_SR_RA_vs_12b_SR_LA_c2.0_cond1.bed_length$sample <- "diff_33_SR_RA_vs_12b_SR_LA_c2.0_cond1.bed"  
diff_33_SR_RA_vs_20_SR_LA_c2.0_cond1.bed_length$sample <- "diff_33_SR_RA_vs_20_SR_LA_c2.0_cond1.bed"  
diff_33_SR_RA_vs_2231_AF_LA_c2.0_cond1.bed_length$sample <- "diff_33_SR_RA_vs_2231_AF_LA_c2.0_cond1.bed"  
diff_33_SR_RA_vs_25_AF_LA_c2.0_cond1.bed_length$sample <- "diff_33_SR_RA_vs_25_AF_LA_c2.0_cond1.bed"  
diff_33_SR_RA_vs_27_SR_LA_c2.0_cond1.bed_length$sample <- "diff_33_SR_RA_vs_27_SR_LA_c2.0_cond1.bed"  
diff_33_SR_RA_vs_28_AF_LA_c2.0_cond1.bed_length$sample <- "diff_33_SR_RA_vs_28_AF_LA_c2.0_cond1.bed" 
diff_33_SR_RA_vs_31_AF_LA_c2.0_cond1.bed_length$sample <- "diff_33_SR_RA_vs_31_AF_LA_c2.0_cond1.bed"  
diff_33_SR_RA_vs_313_AF_LA_c2.0_cond1.bed_length$sample <- "diff_33_SR_RA_vs_313_AF_LA_c2.0_cond1.bed"  
diff_33_SR_RA_vs_45_SR_LA_c2.0_cond1.bed_length$sample <- "diff_33_SR_RA_vs_45_SR_LA_c2.0_cond1.bed"  
diff_33_SR_RA_vs_47_AF_LA_c2.0_cond1.bed_length$sample <- "diff_33_SR_RA_vs_47_AF_LA_c2.0_cond1.bed"  

peaksets_lengths_SR_RA_pairwise <- rbind(diff_12b_SR_RA_vs_12b_SR_LA_c2.0_cond1.bed_length,
                                         diff_12b_SR_RA_vs_20_SR_LA_c2.0_cond1.bed_length,
                                         diff_12b_SR_RA_vs_2231_AF_LA_c2.0_cond1.bed_length, 
                                         diff_12b_SR_RA_vs_25_AF_LA_c2.0_cond1.bed_length,
                                         diff_12b_SR_RA_vs_27_SR_LA_c2.0_cond1.bed_length, 
                                         diff_12b_SR_RA_vs_28_AF_LA_c2.0_cond1.bed_length, 
                                         diff_12b_SR_RA_vs_31_AF_LA_c2.0_cond1.bed_length, 
                                         diff_12b_SR_RA_vs_313_AF_LA_c2.0_cond1.bed_length, 
                                         diff_12b_SR_RA_vs_45_SR_LA_c2.0_cond1.bed_length, 
                                         diff_12b_SR_RA_vs_47_AF_LA_c2.0_cond1.bed_length, 
                                         diff_18_AF_RA_vs_12b_SR_RA_c2.0_cond2.bed_length ,
                                         diff_18_AF_RA_vs_20_SR_RA_c2.0_cond2.bed_length ,
                                         diff_18_AF_RA_vs_27_SR_RA_c2.0_cond2.bed_length ,
                                         diff_18_AF_RA_vs_33_SR_RA_c2.0_cond2.bed_length ,
                                         diff_20_SR_RA_vs_12b_SR_LA_c2.0_cond1.bed_length ,
                                         diff_20_SR_RA_vs_20_SR_LA_c2.0_cond1.bed_length ,
                                         diff_20_SR_RA_vs_2231_AF_LA_c2.0_cond1.bed_length, 
                                         diff_20_SR_RA_vs_25_AF_LA_c2.0_cond1.bed_length ,
                                         diff_20_SR_RA_vs_27_SR_LA_c2.0_cond1.bed_length ,
                                         diff_20_SR_RA_vs_28_AF_LA_c2.0_cond1.bed_length ,
                                         diff_20_SR_RA_vs_31_AF_LA_c2.0_cond1.bed_length ,
                                         diff_20_SR_RA_vs_313_AF_LA_c2.0_cond1.bed_length ,
                                         diff_20_SR_RA_vs_45_SR_LA_c2.0_cond1.bed_length ,
                                         diff_20_SR_RA_vs_47_AF_LA_c2.0_cond1.bed_length ,
                                         diff_2231_AF_RA_vs_12b_SR_RA_c2.0_cond2.bed_length, 
                                         diff_2231_AF_RA_vs_20_SR_RA_c2.0_cond2.bed_length ,
                                         diff_2231_AF_RA_vs_27_SR_RA_c2.0_cond2.bed_length ,
                                         diff_2231_AF_RA_vs_33_SR_RA_c2.0_cond2.bed_length ,
                                         diff_25_AF_RA_vs_12b_SR_RA_c2.0_cond2.bed_length ,
                                         diff_25_AF_RA_vs_20_SR_RA_c2.0_cond2.bed_length ,
                                         diff_25_AF_RA_vs_27_SR_RA_c2.0_cond2.bed_length ,
                                         diff_25_AF_RA_vs_33_SR_RA_c2.0_cond2.bed_length ,
                                         diff_27_SR_RA_vs_12b_SR_LA_c2.0_cond1.bed_length ,
                                         diff_27_SR_RA_vs_20_SR_LA_c2.0_cond1.bed_length ,
                                         diff_27_SR_RA_vs_2231_AF_LA_c2.0_cond1.bed_length, 
                                         diff_27_SR_RA_vs_25_AF_LA_c2.0_cond1.bed_length ,
                                         diff_27_SR_RA_vs_27_SR_LA_c2.0_cond1.bed_length ,
                                         diff_27_SR_RA_vs_28_AF_LA_c2.0_cond1.bed_length ,
                                         diff_27_SR_RA_vs_31_AF_LA_c2.0_cond1.bed_length ,
                                         diff_27_SR_RA_vs_313_AF_LA_c2.0_cond1.bed_length ,
                                         diff_27_SR_RA_vs_45_SR_LA_c2.0_cond1.bed_length ,
                                         diff_27_SR_RA_vs_47_AF_LA_c2.0_cond1.bed_length ,
                                         diff_28_AF_RA_vs_12b_SR_RA_c2.0_cond2.bed_length ,
                                         diff_28_AF_RA_vs_20_SR_RA_c2.0_cond2.bed_length ,
                                         diff_28_AF_RA_vs_27_SR_RA_c2.0_cond2.bed_length ,
                                         diff_28_AF_RA_vs_33_SR_RA_c2.0_cond2.bed_length ,
                                         diff_31_AF_RA_vs_12b_SR_RA_c2.0_cond2.bed_length ,
                                         diff_31_AF_RA_vs_20_SR_RA_c2.0_cond2.bed_length ,
                                         diff_31_AF_RA_vs_27_SR_RA_c2.0_cond2.bed_length ,
                                         diff_31_AF_RA_vs_33_SR_RA_c2.0_cond2.bed_length ,
                                         diff_313_AF_RA_vs_12b_SR_RA_c2.0_cond2.bed_length ,
                                         diff_313_AF_RA_vs_20_SR_RA_c2.0_cond2.bed_length ,
                                         diff_313_AF_RA_vs_27_SR_RA_c2.0_cond2.bed_length ,
                                         diff_313_AF_RA_vs_33_SR_RA_c2.0_cond2.bed_length ,
                                         diff_33_SR_RA_vs_12b_SR_LA_c2.0_cond1.bed_length ,
                                         diff_33_SR_RA_vs_20_SR_LA_c2.0_cond1.bed_length ,
                                         diff_33_SR_RA_vs_2231_AF_LA_c2.0_cond1.bed_length, 
                                         diff_33_SR_RA_vs_25_AF_LA_c2.0_cond1.bed_length ,
                                         diff_33_SR_RA_vs_27_SR_LA_c2.0_cond1.bed_length ,
                                         diff_33_SR_RA_vs_28_AF_LA_c2.0_cond1.bed_length,
                                         diff_33_SR_RA_vs_31_AF_LA_c2.0_cond1.bed_length ,
                                         diff_33_SR_RA_vs_313_AF_LA_c2.0_cond1.bed_length ,
                                         diff_33_SR_RA_vs_45_SR_LA_c2.0_cond1.bed_length ,
                                         diff_33_SR_RA_vs_47_AF_LA_c2.0_cond1.bed_length )

sample_peakset_lengths_boxplot <- ggplot(peaksets_lengths_SR_RA_pairwise, aes(x=sample, y=value)) + geom_boxplot(aes(fill=group),outlier.size = 0.05) + scale_y_continuous("Peak length (bp)") + scale_x_discrete("enriched_peakset_comparison") + theme_classic() + theme(axis.text.x = element_text(angle=45, hjust = 1) )  + scale_fill_manual(values=c('#91C7B1')) + theme(axis.text.x = element_blank())                                                                                                                                                                                                                                                                                                                                                                                                                                                                      
sample_peakset_lengths_boxplot
sample_peakset_lengths_boxplot_nooutliers <- ggplot(peaksets_lengths_SR_RA_pairwise, aes(x=sample, y=value)) + geom_boxplot(aes(fill=group),outlier.shape=NA) + coord_cartesian(ylim = c(0, 1000))+ scale_y_continuous("Peak length (bp)") + scale_x_discrete("enriched_peakset_comparison") + theme_classic() + theme(axis.text.x = element_text(angle=45, hjust = 1) ) + scale_fill_manual(values=c('#91C7B1'))        
sample_peakset_lengths_boxplot_nooutliers
sample_peakset_lengths_boxplot_nooutliers <- ggplot(peaksets_lengths_SR_RA_pairwise, aes(x=sample, y=value)) + geom_boxplot(aes(fill=group),outlier.shape=NA) + coord_cartesian(ylim = c(0, 2500))+ scale_y_continuous("Peak length (bp)") + scale_x_discrete("enriched_peakset_comparison") + theme_classic() + theme(axis.text.x = element_text(angle=45, hjust = 1) ) + scale_fill_manual(values=c('#91C7B1'))        
sample_peakset_lengths_boxplot_nooutliers
sample_peakset_lengths_boxplot_nooutliers_noxaxis <- ggplot(peaksets_lengths_SR_RA_pairwise, aes(x=sample, y=value)) + geom_boxplot(aes(fill=group),outlier.shape=NA) + coord_cartesian(ylim = c(300, 1030))+ scale_y_continuous("Peak length (bp)") + scale_x_discrete("enriched_peakset_comparison") + theme_classic() + theme(axis.text.x = element_blank() ) + scale_fill_manual(values=c('#91C7B1'))        
sample_peakset_lengths_boxplot_nooutliers_noxaxis



#################################################################################################################################################################

### Definition of H3K27ac-enriched regions in each sample group
# 1) Collate pairwise comparisons for each category
# 2) Filter peaksets with Genrich peaks (greylisted regions have been previously removed in these)
# 3) Overlap across categories and retain unique peaks only 

#################################################################################################################################################################


## 1. Generate candidate peakset for each of the four categories (Aim for roughly 1,000 peaks in each peakset)



#AF-RA
setwd("/Rodriguez_2025_AFepigenome/H3K27ac_enriched_regions_definition/AF_RA_enriched_peaks")
AF_RA_enriched_files <- list.files(pattern = "bed")
length(AF_RA_enriched_files) #84
#Prepare empty df file with just header
DBA_AF_RA_enriched <- setNames(data.frame(matrix(ncol = 9, nrow = 84)), c("SampleID", "Tissue", "Factor", "Condition", "bamReads", "ControlID", "bamControl", "PeakCaller","Peaks"))
DBA_AF_RA_enriched$SampleID <- AF_RA_enriched_files
DBA_AF_RA_enriched$Peaks <- AF_RA_enriched_files
DBA_AF_RA_enriched
#Convert to csv to complete editing
write.csv(DBA_AF_RA_enriched,"DBA_AF_RA_enriched.csv")
## csv file requires manual editing to fill in columns
    # Tissue: copy values from the "SampleID" column
    # Factor: set to "H3K27ac" for all rows
    # Condition: copy values from the "SampleID" column
    # ControlID: copy values from the "SampleID" column
    # PeakCaller:  set to "bed" for all rows
DBA_AF_RA_enriched_consensus_minOverlap35 <- dba(sampleSheet="DBA_AF_RA_enriched.csv",minOverlap = 35) #1,233 peaks
DBA_AF_RA_enriched_consensus_minOverlap35_peaks <- dba.peakset(DBA_AF_RA_enriched_consensus_minOverlap35, bRetrieve=TRUE)
DBA_AF_RA_enriched_consensus_olap.rate <- dba.overlap(DBA_AF_RA_enriched_consensus_minOverlap35,mode=DBA_OLAP_RATE)
plot(DBA_AF_RA_enriched_consensus_olap.rate,type='b',ylab='# peaks',xlab='Overlap at least this many peaksets')


#AF-LA
setwd("/Rodriguez_2025_AFepigenome/H3K27ac_enriched_regions_definition/AF_LA_enriched_peaks")
AF_LA_enriched_files <- list.files(pattern = "bed")
length(AF_LA_enriched_files) #84
#Prepare empty df file with just header
DBA_AF_LA_enriched <- setNames(data.frame(matrix(ncol = 9, nrow = 84)), c("SampleID", "Tissue", "Factor", "Condition", "bamReads", "ControlID", "bamControl", "PeakCaller","Peaks"))
DBA_AF_LA_enriched$SampleID <- AF_LA_enriched_files
DBA_AF_LA_enriched$Peaks <- AF_LA_enriched_files
DBA_AF_LA_enriched
#Convert to csv to complete editing
write.csv(DBA_AF_LA_enriched,"DBA_AF_LA_enriched.csv")
## csv file requires manual editing to fill in columns
    # Tissue: copy values from the "SampleID" column
    # Factor: set to "H3K27ac" for all rows
    # Condition: copy values from the "SampleID" column
    # ControlID: copy values from the "SampleID" column
    # PeakCaller:  set to "bed" for all rows
DBA_AF_LA_enriched_consensus_minOverlap42 <- dba(sampleSheet="DBA_AF_LA_enriched.csv",minOverlap = 42) #1,296 peaks
DBA_AF_LA_enriched_consensus_minOverlap42_peaks <- dba.peakset(DBA_AF_LA_enriched_consensus_minOverlap42, bRetrieve=TRUE)
DBA_AF_LA_enriched_consensus_olap.rate <- dba.overlap(DBA_AF_LA_enriched_consensus_minOverlap42,mode=DBA_OLAP_RATE)
plot(DBA_AF_LA_enriched_consensus_olap.rate,type='b',ylab='# peaks',xlab='Overlap at least this many peaksets')


#SR-RA
setwd("/Rodriguez_2025_AFepigenome/H3K27ac_enriched_regions_definition/SR_RA_enriched_peaks")
SR_RA_enriched_files <- list.files(pattern = "bed")
length(SR_RA_enriched_files) #64
#Prepare empty df file with just header
DBA_SR_RA_enriched <- setNames(data.frame(matrix(ncol = 9, nrow = 64)), c("SampleID", "Tissue", "Factor", "Condition", "bamReads", "ControlID", "bamControl", "PeakCaller","Peaks"))
DBA_SR_RA_enriched$SampleID <- SR_RA_enriched_files
DBA_SR_RA_enriched$Peaks <- SR_RA_enriched_files
DBA_SR_RA_enriched
#Convert to csv to complete editing
write.csv(DBA_SR_RA_enriched,"DBA_SR_RA_enriched.csv")
## csv file requires manual editing to fill in columns
    # Tissue: copy values from the "SampleID" column
    # Factor: set to "H3K27ac" for all rows
    # Condition: copy values from the "SampleID" column
    # ControlID: copy values from the "SampleID" column
    # PeakCaller:  set to "bed" for all rows
DBA_SR_RA_enriched_consensus_minOverlap10 <- dba(sampleSheet="DBA_SR_RA_enriched.csv",minOverlap = 10) #1,311 peaks
DBA_SR_RA_enriched_consensus_minOverlap10_peaks <- dba.peakset(DBA_SR_RA_enriched_consensus_minOverlap10, bRetrieve=TRUE)
DBA_SR_RA_enriched_consensus_olap.rate <- dba.overlap(DBA_SR_RA_enriched_consensus_minOverlap10,mode=DBA_OLAP_RATE)
plot(DBA_SR_RA_enriched_consensus_olap.rate,type='b',ylab='# peaks',xlab='Overlap at least this many peaksets')


#SR-LA
setwd("/Rodriguez_2025_AFepigenome/H3K27ac_enriched_regions_definition/SR_LA_enriched_peaks")
SR_LA_enriched_files <- list.files(pattern = "bed")
length(SR_LA_enriched_files) #64
#Prepare empty df file with just header
DBA_SR_LA_enriched <- setNames(data.frame(matrix(ncol = 9, nrow = 64)), c("SampleID", "Tissue", "Factor", "Condition", "bamReads", "ControlID", "bamControl", "PeakCaller","Peaks"))
DBA_SR_LA_enriched$SampleID <- SR_LA_enriched_files
DBA_SR_LA_enriched$Peaks <- SR_LA_enriched_files
DBA_SR_LA_enriched
#Convert to csv to complete editing
write.csv(DBA_SR_LA_enriched,"DBA_SR_LA_enriched.csv")
## csv file requires manual editing to fill in columns
    # Tissue: copy values from the "SampleID" column
    # Factor: set to "H3K27ac" for all rows
    # Condition: copy values from the "SampleID" column
    # ControlID: copy values from the "SampleID" column
    # PeakCaller:  set to "bed" for all rows
DBA_SR_LA_enriched_consensus_minOverlap26 <- dba(sampleSheet="DBA_SR_LA_enriched.csv",minOverlap = 26) #1,284 peaks
DBA_SR_LA_enriched_consensus_minOverlap26_peaks <- dba.peakset(DBA_SR_LA_enriched_consensus_minOverlap26, bRetrieve=TRUE)
DBA_SR_LA_enriched_consensus_olap.rate <- dba.overlap(DBA_SR_LA_enriched_consensus_minOverlap26,mode=DBA_OLAP_RATE)
plot(DBA_SR_LA_enriched_consensus_olap.rate,type='b',ylab='# peaks',xlab='Overlap at least this many peaksets')


#AF
setwd("/Rodriguez_2025_AFepigenome/H3K27ac_enriched_regions_definition/AF_enriched_peaks")
AF_enriched_files <- list.files(pattern = "bed")
length(AF_enriched_files) #96
#Prepare empty df file with just header
DBA_AF_enriched <- setNames(data.frame(matrix(ncol = 9, nrow = 96)), c("SampleID", "Tissue", "Factor", "Condition", "bamReads", "ControlID", "bamControl", "PeakCaller","Peaks"))
DBA_AF_enriched$SampleID <- AF_enriched_files
DBA_AF_enriched$Peaks <- AF_enriched_files
AF_enriched_files
#Convert to csv to complete editing
## csv file requires manual editing to fill in columns
    # Tissue: copy values from the "SampleID" column
    # Factor: set to "H3K27ac" for all rows
    # Condition: copy values from the "SampleID" column
    # ControlID: copy values from the "SampleID" column
    # PeakCaller:  set to "bed" for all rows
DBA_AF_enriched_consensus_minOverlap48 <- dba(sampleSheet="DBA_AF_enriched.csv",minOverlap = 48) #1,005 peaks
DBA_AF_enriched_consensus_minOverlap48_peaks <- dba.peakset(DBA_AF_enriched_consensus_minOverlap48, bRetrieve=TRUE)
DBA_AF_enriched_consensus_olap.rate <- dba.overlap(DBA_AF_enriched_consensus_minOverlap48,mode=DBA_OLAP_RATE)
plot(DBA_AF_enriched_consensus_olap.rate,type='b',ylab='# peaks',xlab='Overlap at least this many peaksets')


#SR
setwd("/Rodriguez_2025_AFepigenome/H3K27ac_enriched_regions_definition/SR_enriched_peaks")
SR_enriched_files <- list.files(pattern = "bed")
length(SR_enriched_files) #96
#Prepare empty df file with just header
DBA_SR_enriched <- setNames(data.frame(matrix(ncol = 9, nrow = 96)), c("SampleID", "Tissue", "Factor", "Condition", "bamReads", "ControlID", "bamControl", "PeakCaller","Peaks"))
DBA_SR_enriched$SampleID <- SR_enriched_files
DBA_SR_enriched$Peaks <- SR_enriched_files
SR_enriched_files
#Convert to csv to complete editing
## csv file requires manual editing to fill in columns
    # Tissue: copy values from the "SampleID" column
    # Factor: set to "H3K27ac" for all rows
    # Condition: copy values from the "SampleID" column
    # ControlID: copy values from the "SampleID" column
    # PeakCaller:  set to "bed" for all rows
write.csv(DBA_SR_enriched,"DBA_SR_enriched.csv")
DBA_SR_enriched_consensus_minOverlap24 <- dba(sampleSheet="DBA_SR_enriched.csv",minOverlap = 24) #1,112 peaks
DBA_SR_enriched_consensus_minOverlap24_peaks <- dba.peakset(DBA_SR_enriched_consensus_minOverlap24, bRetrieve=TRUE)
DB_SR_enriched_consensus_olap.rate <- dba.overlap(DBA_SR_enriched_consensus_minOverlap24,mode=DBA_OLAP_RATE)
plot(DBA_SR_enriched_consensus_olap.rate,type='b',ylab='# peaks',xlab='Overlap at least this many peaksets')



#################################################################################################################################################################


## 2. Filter with Genrich peaks 

#### Re-define peaksets but correcting for peaks redundancy using Genrich peaks (to subset larger peaks into smaller ones) --> using bedtools intersect
### Greylisted regions have been removed from the Genrich peaks



setwd("/Rodriguez_2025_AFepigenome/H3K27ac_enriched_regions_definition/")

#AF-RA-enriched
DBA_AF_RA_enriched_consensus_minOverlap35_peaks
granges(DBA_AF_RA_enriched_consensus_minOverlap35_peaks)
DBA_AF_RA_enriched_consensus_minOverlap35_peaks_df <- annoGR2DF(DBA_AF_RA_enriched_consensus_minOverlap35_peaks)
DBA_AF_RA_enriched_consensus_minOverlap35_peaks_df <- subset(DBA_AF_RA_enriched_consensus_minOverlap35_peaks_df, select = c(chr, start, end))
write.table(DBA_AF_RA_enriched_consensus_minOverlap35_peaks_df, "DBA_AF_RA_enriched_consensus.bed",row.names = F,col.names = F, sep="\t", quote=FALSE)


#AF-LA-enriched
DBA_AF_LA_enriched_consensus_minOverlap42_peaks
granges(DBA_AF_LA_enriched_consensus_minOverlap42_peaks)
DBA_AF_LA_enriched_consensus_minOverlap42_peaks_df <- annoGR2DF(DBA_AF_LA_enriched_consensus_minOverlap42_peaks)
DBA_AF_LA_enriched_consensus_minOverlap42_peaks_df <- subset(DBA_AF_LA_enriched_consensus_minOverlap42_peaks_df, select = c(chr, start, end))
write.table(DBA_AF_LA_enriched_consensus_minOverlap42_peaks_df, "DBA_AF_LA_enriched_consensus.bed",row.names = F,col.names = F, sep="\t", quote=FALSE)


#SR-RA-enriched
DBA_SR_RA_enriched_consensus_minOverlap10_peaks
granges(DBA_SR_RA_enriched_consensus_minOverlap10_peaks)
DBA_SR_RA_enriched_consensus_minOverlap10_peaks_df <- annoGR2DF(DBA_SR_RA_enriched_consensus_minOverlap10_peaks)
DBA_SR_RA_enriched_consensus_minOverlap10_peaks_df <- subset(DBA_SR_RA_enriched_consensus_minOverlap10_peaks_df, select = c(chr, start, end))
write.table(DBA_SR_RA_enriched_consensus_minOverlap10_peaks_df, "DBA_SR_RA_enriched_consensus.bed",row.names = F,col.names = F, sep="\t", quote=FALSE)


#SR-LA-enriched
DBA_SR_LA_enriched_consensus_minOverlap26_peaks
granges(DBA_SR_LA_enriched_consensus_minOverlap26_peaks)
DBA_SR_LA_enriched_consensus_minOverlap26_peaks_df <- annoGR2DF(DBA_SR_LA_enriched_consensus_minOverlap26_peaks)
DBA_SR_LA_enriched_consensus_minOverlap26_peaks_df <- subset(DBA_SR_LA_enriched_consensus_minOverlap26_peaks_df, select = c(chr, start, end))
write.table(DBA_SR_LA_enriched_consensus_minOverlap26_peaks_df, "DBA_SR_LA_enriched_consensus.bed",row.names = F,col.names = F, sep="\t", quote=FALSE)



#AF-enriched
DBA_AF_enriched_consensus_minOverlap48_peaks
granges(DBA_AF_enriched_consensus_minOverlap48_peaks)
DBA_AF_enriched_consensus_minOverlap48_peaks_df <- annoGR2DF(DBA_AF_enriched_consensus_minOverlap48_peaks)
DBA_AF_enriched_consensus_minOverlap48_peaks_df <- subset(DBA_AF_enriched_consensus_minOverlap48_peaks_df, select = c(chr, start, end))
write.table(DBA_AF_enriched_consensus_minOverlap48_peaks_df, "DBA_AF_enriched_consensus.bed",row.names = F,col.names = F, sep="\t", quote=FALSE)

#SR-enriched
DBA_SR_enriched_consensus_minOverlap24_peaks
granges(DBA_SR_enriched_consensus_minOverlap24_peaks)
DBA_SR_enriched_consensus_minOverlap24_peaks_df <- annoGR2DF(DBA_SR_enriched_consensus_minOverlap24_peaks)
DBA_SR_enriched_consensus_minOverlap24_peaks_df <- subset(DBA_SR_enriched_consensus_minOverlap24_peaks_df, select = c(chr, start, end))
write.table(DBA_SR_enriched_consensus_minOverlap24_peaks_df, "DBA_SR_enriched_consensus.bed",row.names = F,col.names = F, sep="\t", quote=FALSE)


# NOTE:
# - These commands should be run in a UNIX shell (e.g., on the cluster)
# - Make sure bedtools is loaded (module load bedtools)
# - Output files will be used as input in the next R step.

# RUN IN TERMINAL:

# module load bedtools 
# bedtools intersect -a DBA_AF_LA_enriched_consensus.bed -b AF_LA_Genrich_greylistExcluded_g10.bed > DBA_AF_LA_enriched_consensus_Genrich_filtered.bed
# bedtools intersect -a DBA_AF_RA_enriched_consensus.bed -b AF_RA_Genrich_greylistExcluded_g10.bed > DBA_AF_RA_enriched_consensus_Genrich_filtered.bed
# bedtools intersect -a DBA_SR_LA_enriched_consensus.bed -b SR_LA_Genrich_greylistExcluded_g10.bed > DBA_SR_LA_enriched_consensus_Genrich_filtered.bed
# bedtools intersect -a DBA_SR_RA_enriched_consensus.bed -b SR_RA_Genrich_greylistExcluded_g10.bed > DBA_SR_RA_enriched_consensus_Genrich_filtered.bed

# bedtools intersect -a DBA_AF_enriched_consensus.bed -b AF_Genrich_greylistExlcuded_g10.bed > DBA_AF_enriched_consensus_Genrich_filtered.bed
# bedtools intersect -a DBA_SR_enriched_consensus.bed -b SR_Genrich_greylistExcluded_g10.bed > DBA_SR_enriched_consensus_Genrich_filtered.bed


#################################################################################################################################################################


## 2. Explore overlap across peaksets 
## Get rid of overlapping peaks to define final peaksets in each category


setwd("/Rodriguez_2025_AFepigenome/H3K27ac_enriched_regions_definition/")


#AF-RA-enriched
DBA_AF_RA_enriched_consensus_Genrich_filtered <- as.data.frame(read.table("DBA_AF_RA_enriched_consensus_Genrich_filtered.bed",header = FALSE, sep="\t",stringsAsFactors=FALSE))
colnames(DBA_AF_RA_enriched_consensus_Genrich_filtered) <- c("chr","start","end")
DBA_AF_RA_enriched_consensus_Genrich_filtered_gr <- makeGRangesFromDataFrame(DBA_AF_RA_enriched_consensus_Genrich_filtered)

#AF-LA-enriched
DBA_AF_LA_enriched_consensus_Genrich_filtered <- as.data.frame(read.table("DBA_AF_LA_enriched_consensus_Genrich_filtered.bed",header = FALSE, sep="\t",stringsAsFactors=FALSE))
colnames(DBA_AF_LA_enriched_consensus_Genrich_filtered) <- c("chr","start","end")
DBA_AF_LA_enriched_consensus_Genrich_filtered_gr <- makeGRangesFromDataFrame(DBA_AF_LA_enriched_consensus_Genrich_filtered)


#SR-RA-enriched
DBA_SR_RA_enriched_consensus_Genrich_filtered <- as.data.frame(read.table("DBA_SR_RA_enriched_consensus_Genrich_filtered.bed",header = FALSE, sep="\t",stringsAsFactors=FALSE))
colnames(DBA_SR_RA_enriched_consensus_Genrich_filtered) <- c("chr","start","end")
DBA_SR_RA_enriched_consensus_Genrich_filtered_gr <- makeGRangesFromDataFrame(DBA_SR_RA_enriched_consensus_Genrich_filtered)


#SR-LA-enriched
DBA_SR_LA_enriched_consensus_Genrich_filtered <- as.data.frame(read.table("DBA_SR_LA_enriched_consensus_Genrich_filtered.bed",header = FALSE, sep="\t",stringsAsFactors=FALSE))
colnames(DBA_SR_LA_enriched_consensus_Genrich_filtered) <- c("chr","start","end")
DBA_SR_LA_enriched_consensus_Genrich_filtered_gr <- makeGRangesFromDataFrame(DBA_SR_LA_enriched_consensus_Genrich_filtered)


#AF-enriched
DBA_AF_enriched_consensus_Genrich_filtered <- as.data.frame(read.table("DBA_AF_enriched_consensus_Genrich_filtered.bed",header = FALSE, sep="\t",stringsAsFactors=FALSE))
colnames(DBA_AF_enriched_consensus_Genrich_filtered) <- c("chr","start","end")
DBA_AF_enriched_consensus_Genrich_filtered_gr <- makeGRangesFromDataFrame(DBA_AF_enriched_consensus_Genrich_filtered)


#SR-enriched
DBA_SR_enriched_consensus_Genrich_filtered <- as.data.frame(read.table("DBA_SR_enriched_consensus_Genrich_filtered.bed",header = FALSE, sep="\t",stringsAsFactors=FALSE))
colnames(DBA_SR_enriched_consensus_Genrich_filtered) <- c("chr","start","end")
DBA_SR_enriched_consensus_Genrich_filtered_gr <- makeGRangesFromDataFrame(DBA_SR_enriched_consensus_Genrich_filtered)


library(ChIPpeakAnno)


# Illustrate Venn diagram to explore overlapping regions
makeVennDiagram(list(DBA_AF_RA_enriched_consensus_Genrich_filtered_gr, DBA_AF_LA_enriched_consensus_Genrich_filtered_gr,DBA_SR_RA_enriched_consensus_Genrich_filtered_gr,DBA_SR_LA_enriched_consensus_Genrich_filtered_gr), NameOfPeaks=c("AF-RA-enriched", "AF-LA-enriched","SR-RA-enriched","SR-LA-enriched"),
                totalTest=2000,scaled=FALSE, euler.d=FALSE, 
                fill=c("#E3D081", "#D85B74","#91C7B1","#7298BC"), # circle fill color
                col=c("black", "black","black","black"), #circle border color
                cat.col=c("black", "black","black","black"))

makeVennDiagram(list(DBA_AF_enriched_consensus_Genrich_filtered_gr, DBA_SR_enriched_consensus_Genrich_filtered_gr), NameOfPeaks=c("SR-enriched", "AF-enriched"),
                totalTest=2000,scaled=FALSE, euler.d=FALSE, 
                fill=c("#F4997F","#43B9D1"), # circle fill color
                col=c("black", "black"), #circle border color
                cat.col=c("black", "black"))



# Area-proportional Venn diagrams

#devtools::install_github("vqf/nVennR")
library(nVennR)

myV <- createVennObj(nSets = 4, sNames = c('AF-RA', 'AF-LA', 'SR-RA', 'SR-LA'), 
                     sSizes = c(0, 955, 1027, 53, 920, 281, 46, 34, 1205, 65, 88, 14, 101, 80, 34, 43))
myV <- plotVenn(nVennObj = myV, outFile="VennDiagram_samplegroups_peaksets_updatedAug24.svg")
# 📊 Figure 2a
showSVG(nVennObj = myV, setColor= c('#E3D081', '#D85B74', '#91C7B1', '#7298BC'), outFile="VennDiagram_samplegroups_peaksets_updatedAug24.svg")


setwd("/Users/adrianrodriguez/Downloads")
library(nVennR)
myV <- createVennObj(nSets = 2, sNames = c('AF', 'SR'), 
                     sSizes = c(0, 998, 855, 264))
myV <- plotVenn(nVennObj = myV, outFile="VennDiagram_samplegroups_peaksets_updated-AFSR_Aug24.svg")
showSVG(nVennObj = myV, setColor= c('#F4997F', '#43B9D1'), outFile="VennDiagram_samplegroups_peaksets_updated-AFSR_Aug24.svg")


# Retain unique peaks only 
Overlap_across_4_sets <- findOverlapsOfPeaks(DBA_AF_RA_enriched_consensus_Genrich_filtered_gr, DBA_AF_LA_enriched_consensus_Genrich_filtered_gr,DBA_SR_RA_enriched_consensus_Genrich_filtered_gr,DBA_SR_LA_enriched_consensus_Genrich_filtered_gr)

All_unique_peaks <- Overlap_across_4_sets$uniquePeaks #4,107
library(Repitools)
All_unique_peaks_df <- annoGR2DF(granges(All_unique_peaks))
rownames(All_unique_peaks_df)


Overlap_across_AFandSR_sets <- findOverlapsOfPeaks(DBA_AF_enriched_consensus_Genrich_filtered_gr,DBA_SR_enriched_consensus_Genrich_filtered_gr)
library(Repitools)
All_unique_peaks_AFandSR <- Overlap_across_AFandSR_sets$uniquePeaks #1,853
All_unique_peaks_AFandSR_df <- annoGR2DF(granges(All_unique_peaks_AFandSR))
rownames(All_unique_peaks_AFandSR_df)


## The following code outputs the final sets of regions defined as H3K27ac-enriched for each group: AF, SR, AF-LA, AF-RA, SR-LA and SR-RA (as shown in 📎 Supplementary Table 3).

#AF unique peaks
All_unique_peaks_AFandSR_df
AF_enriched_rownumbers <-  grep("^DBA_AF_enriched_consensus",rownames(All_unique_peaks_AFandSR_df))
AF_enriched_nonoverlapping_df <- All_unique_peaks_AFandSR_df[AF_enriched_rownumbers,]
dim(AF_enriched_nonoverlapping_df) #855
AF_enriched_nonoverlapping_gr <- makeGRangesFromDataFrame(AF_enriched_nonoverlapping_df)
setwd("/Users/adrianrodriguez/Desktop/PhD/ChIP-Seq_data_analysis/AnalysisNov23/Analyses_UPDATED_Aug24/c_2") 
AF_enriched_nonoverlapping_df$peak_name <- rownames(AF_enriched_nonoverlapping_df)
AF_enriched_nonoverlapping_df <- AF_enriched_nonoverlapping_df[c(-4)]
write.table(AF_enriched_nonoverlapping_df, "AF_enriched_nonoverlappingwithallsets_GenrichGreylistedFiltered.bed",row.names = F,col.names = F, sep="\t", quote=FALSE)


#SR unique peaks
All_unique_peaks_AFandSR_df
SR_enriched_rownumbers <-  grep("^DBA_SR_enriched_consensus",rownames(All_unique_peaks_AFandSR_df))
SR_enriched_nonoverlapping_df <- All_unique_peaks_AFandSR_df[SR_enriched_rownumbers,]
dim(SR_enriched_nonoverlapping_df) #921
SR_enriched_nonoverlapping_gr <- makeGRangesFromDataFrame(SR_enriched_nonoverlapping_df)
setwd("/Users/adrianrodriguez/Desktop/PhD/ChIP-Seq_data_analysis/AnalysisNov23/Analyses_UPDATED_Aug24/c_2") 
SR_enriched_nonoverlapping_df$peak_name <- rownames(SR_enriched_nonoverlapping_df)
SR_enriched_nonoverlapping_df <- SR_enriched_nonoverlapping_df[c(-4)]
write.table(SR_enriched_nonoverlapping_df, "SR_enriched_nonoverlappingwithallsets_GenrichGreylistedFiltered.bed",row.names = F,col.names = F, sep="\t", quote=FALSE)


#AF_RA unique peaks
#Extract from All_unique_peaks_df those rows corresponding to AF-RA-enriched unique 
All_unique_peaks_df
AF_RA_enriched_rownumbers <- grep("^DBA_AF_RA_enriched_consensus",rownames(All_unique_peaks_df))
AF_RA_enriched_nonoverlappingwithallsets_df <- All_unique_peaks_df[AF_RA_enriched_rownumbers,]
dim(AF_RA_enriched_nonoverlappingwithallsets_df) #1,205
AF_RA_enriched_nonoverlappingwithallsets_gr <- makeGRangesFromDataFrame(AF_RA_enriched_nonoverlappingwithallsets_df)
setwd("/Users/adrianrodriguez/Desktop/PhD/ChIP-Seq_data_analysis/AnalysisNov23/Analyses_UPDATED_Aug24/c_2") 
AF_RA_enriched_nonoverlappingwithallsets_df$peak_name <- rownames(AF_RA_enriched_nonoverlappingwithallsets_df)
AF_RA_enriched_nonoverlappingwithallsets_df <- AF_RA_enriched_nonoverlappingwithallsets_df[c(-4)]
write.table(AF_RA_enriched_nonoverlappingwithallsets_df, "AF_RA_enriched_nonoverlappingwithallsets_GenrichGreylistedFiltered.bed",row.names = F,col.names = F, sep="\t", quote=FALSE)

#AF_LA unique peaks
#Extract from All_unique_peaks_df those rows corresponding to AF-LA-enriched unique
All_unique_peaks_df
AF_LA_enriched_rownumbers <- grep("^DBA_AF_LA_enriched_consensus",rownames(All_unique_peaks_df))
AF_LA_enriched_nonoverlappingwithallsets_df <- All_unique_peaks_df[AF_LA_enriched_rownumbers,]
dim(AF_LA_enriched_nonoverlappingwithallsets_df)
AF_LA_enriched_nonoverlappingwithallsets_gr <- makeGRangesFromDataFrame(AF_LA_enriched_nonoverlappingwithallsets_df)
setwd("/Users/adrianrodriguez/Desktop/PhD/ChIP-Seq_data_analysis/AnalysisNov23/Analyses_UPDATED_Aug24/c_2") 
AF_LA_enriched_nonoverlappingwithallsets_df$peak_name <- rownames(AF_LA_enriched_nonoverlappingwithallsets_df)
AF_LA_enriched_nonoverlappingwithallsets_df <- AF_LA_enriched_nonoverlappingwithallsets_df[c(-4)]
write.table(AF_LA_enriched_nonoverlappingwithallsets_df, "AF_LA_enriched_nonoverlappingwithallsets_GenrichGreylistedFiltered.bed",row.names = F,col.names = F, sep="\t", quote=FALSE)

#SR_RA unique peaks
#Extract from All_unique_peaks_df those rows corresponding to SR-enriched unique 
All_unique_peaks_df
SR_RA_enriched_rownumbers <- grep("^DBA_SR_RA_enriched_consensus",rownames(All_unique_peaks_df))
SR_RA_enriched_nonoverlappingwithallsets_df <- All_unique_peaks_df[SR_RA_enriched_rownumbers,]
dim(SR_RA_enriched_nonoverlappingwithallsets_df)
SR_RA_enriched_nonoverlappingwithallsets_gr <- makeGRangesFromDataFrame(SR_RA_enriched_nonoverlappingwithallsets_df)
setwd("/Users/adrianrodriguez/Desktop/PhD/ChIP-Seq_data_analysis/AnalysisNov23/Analyses_UPDATED_Aug24/c_2") 
SR_RA_enriched_nonoverlappingwithallsets_df$peak_name <- rownames(SR_RA_enriched_nonoverlappingwithallsets_df)
SR_RA_enriched_nonoverlappingwithallsets_df <- SR_RA_enriched_nonoverlappingwithallsets_df[c(-4)]
write.table(SR_RA_enriched_nonoverlappingwithallsets_df, "SR_RA_enriched_nonoverlappingwithallsets_GenrichGreylistedFiltered.bed",row.names = F,col.names = F, sep="\t", quote=FALSE)


#SR_LA unique peaks
#Extract from All_unique_peaks_df those rows corresponding to AF-enriched unique 
All_unique_peaks_df
SR_LA_enriched_rownumbers <- grep("^DBA_SR_LA_enriched_consensus",rownames(All_unique_peaks_df))
SR_LA_enriched_nonoverlappingwithallsets_df <- All_unique_peaks_df[SR_LA_enriched_rownumbers,]
dim(SR_LA_enriched_nonoverlappingwithallsets_df)
SR_LA_enriched_nonoverlappingwithallsets_df
rownames(SR_LA_enriched_nonoverlappingwithallsets_df)
SR_LA_enriched_nonoverlappingwithallsets_df$peak_name <- rownames(SR_LA_enriched_nonoverlappingwithallsets_df)
SR_LA_enriched_nonoverlappingwithallsets_df <- SR_LA_enriched_nonoverlappingwithallsets_df[c(-4)]
write.table(SR_LA_enriched_nonoverlappingwithallsets_df, "SR_LA_enriched_nonoverlappingwithallsets_GenrichGreylistedFiltered.bed",row.names = F,col.names = F, sep="\t", quote=FALSE)


### Generate bed file for ALL unique peaks to use as background for downstream analyses

All_unique_peaks_df 
dim(All_unique_peaks_df)
rownames(All_unique_peaks_df)
All_unique_peaks_df$peak_name <- rownames(All_unique_peaks_df)
All_unique_peaks_df <- All_unique_peaks_df[c(-4)]
write.table(All_unique_peaks_df, "1All_unique_peaks.bed_GenrichGreylistedFiltered",row.names = F,col.names = F, sep="\t", quote=FALSE)


All_unique_peaks_AFandSR_df 
dim(All_unique_peaks_AFandSR_df)
rownames(All_unique_peaks_AFandSR_df)
All_unique_peaks_AFandSR_df$peak_name <- rownames(All_unique_peaks_AFandSR_df)
All_unique_peaks_AFandSR_df <- All_unique_peaks_AFandSR_df[c(-4)]
write.table(All_unique_peaks_AFandSR_df, "All_unique_peaks_AFandSR_df_GenrichGreylistedFiltered.bed",row.names = F,col.names = F, sep="\t", quote=FALSE)



#################################################################################################################################################################


# Explore average peakset length 


setwd("/Rodriguez_2025_AFepigenome/H3K27ac_enriched_regions_definition/")


AF_LA_enriched_nonoverlappingwithallsets_df <- as.data.frame(read.table("AF_LA_enriched_nonoverlappingwithallsets_GenrichGreylistedFiltered.bed",header = FALSE, sep="\t",stringsAsFactors=FALSE))
AF_RA_enriched_nonoverlappingwithallsets_df <- as.data.frame(read.table("AF_RA_enriched_nonoverlappingwithallsets_GenrichGreylistedFiltered.bed",header = FALSE, sep="\t",stringsAsFactors=FALSE))
SR_LA_enriched_nonoverlappingwithallsets_df <- as.data.frame(read.table("SR_LA_enriched_nonoverlappingwithallsets_GenrichGreylistedFiltered.bed",header = FALSE, sep="\t",stringsAsFactors=FALSE))
SR_RA_enriched_nonoverlappingwithallsets_df <- as.data.frame(read.table("SR_RA_enriched_nonoverlappingwithallsets_GenrichGreylistedFiltered.bed",header = FALSE, sep="\t",stringsAsFactors=FALSE))

AF_RA_peaks_length <- AF_RA_enriched_nonoverlappingwithallsets_df$V3 - AF_RA_enriched_nonoverlappingwithallsets_df$V2
AF_LA_peaks_length <- AF_LA_enriched_nonoverlappingwithallsets_df$V3 - AF_LA_enriched_nonoverlappingwithallsets_df$V2
SR_RA_peaks_length <- SR_RA_enriched_nonoverlappingwithallsets_df$V3 - SR_RA_enriched_nonoverlappingwithallsets_df$V2
SR_LA_peaks_length <- SR_LA_enriched_nonoverlappingwithallsets_df$V3 - SR_LA_enriched_nonoverlappingwithallsets_df$V2


AF_RA_peaks_length <- melt(AF_RA_peaks_length)
AF_RA_peaks_length$class <- "AF_RA"

AF_LA_peaks_length <- melt(AF_LA_peaks_length)
AF_LA_peaks_length$class <- "AF_LA"

SR_RA_peaks_length <- melt(SR_RA_peaks_length)
SR_RA_peaks_length$class <- "SR_RA"

SR_LA_peaks_length <- melt(SR_LA_peaks_length)
SR_LA_peaks_length$class <- "SR_LA"

peaksets_lengths <- rbind(AF_RA_peaks_length,AF_LA_peaks_length,SR_RA_peaks_length,SR_LA_peaks_length)
#📊 Figures S2b 
ggplot(peaksets_lengths, aes(x=class, y=value)) + geom_boxplot(aes(fill=class)) + scale_y_continuous("Peak length (bp)") + scale_x_discrete("Subgroup") + theme_classic() + theme(axis.text.x = element_text(hjust = 0.5, size = 12), axis.title.y = element_text(size = 14)) + scale_fill_manual(values = c("#D85B74", "#E3D081","#7298BC","#91C7B1"))                                                                                                                                                                                                                                                                                                                                                                                                                                                       

#Calculate average length
df_avg <- peaksets_lengths %>%
  group_by(class) %>%
  summarise(avg_value = mean(value, na.rm = TRUE))
df_avg



#################################################################################################################################################################


# Explore reads coverage across defined sets of enriched regions

setwd("/Rodriguez_2025_AFepigenome/H3K27ac_enriched_regions_definition/")


#Assign names of samples corresponding to each category 
AFLA <- c("vl0028","vl0032","vl0212","vl0215","vl0256","vl0257")
AFRA <- c("vl0029","vl0033","vl0260","vl0214","vl0258","vl0326")
SRLA <- c("vl0026","vl0213","vl0217","vl0328")
SRRA <- c("vl0027","vl0031","vl0218","vl0261")

AF <- c("vl0028","vl0032","vl0212","vl0215","vl0256","vl0257","vl0029","vl0033","vl0260","vl0214","vl0258","vl0326")
SR <- c("vl0026","vl0213","vl0217","vl0328","vl0027","vl0031","vl0218","vl0261")

#Establish order of column samples 
colorder3 <- c(AFLA, AFRA, SRLA, SRRA) 
colorder4 <- c(AF, SR) 


data <- read.delim("genrich_c2_merged_fpkm_normalized.csv", header=T, sep="\t")


# Reorder by AFLA, AFRA, SRLA, SRRA
for (i in 1:length(colorder3)){
  lib <- colorder3[i]
  id <- grep(lib, colnames(data))
  if (i==1){so <- id}
  else{so <- c(so, id)}
}
datao <- data[,so]


#Subset rows corresponding to Genrich peaksets and enriched peaksets 

datao$peak_ID <- rownames(datao)

datao_enrichedpeaksets <- datao[grep("enriched_consensus_Genrich_filtered_gr", datao$peak_ID), ]
datao_Genrichpeaksets <- datao[grep("_peak_", datao$peak_ID), ]

datao_enrichedpeaksets <- datao_enrichedpeaksets[, -21]
datao_Genrichpeaksets <- datao_Genrichpeaksets[, -21]


#Represent as boxplots for each group of interest

library(stringr)

#Subset data for each of the four peaksets defined: AFLA, AFRA, SRLA, SRRA

idspl <- str_split(rownames(datao_enrichedpeaksets), "_")
for(i in 1:length(idspl)){
  if(i==1){ids <- idspl[[i]][1]}
  else{ids <- c(ids, idspl[[i]][1]) }
}
data2 <- cbind (datao_enrichedpeaksets, ids)

####################
### AF-LA peakset ##
####################

subset <- data2[data2$ids=="AFLA",]
boxplot(subset[,-21], las=2)#By library

#Aggregate data across same sample type: AFLA, AFRA, SRLA, SRRA


#AF-LA 
id <- match(AFLA, colnames(subset))
set <- subset[,id]
for (i in 1:dim(set)[2]){
  if(i==1){ag <- set[,i]}
  else{ag <- c(ag, set[,i])}
}

c1 <- ag
i1 <- length(c1)
l1 <- rep("AFLA",i1)

#AF-RA
id <- match(AFRA, colnames(subset))
set <- subset[,id]
for (i in 1:dim(set)[2]){
  if(i==1){ag <- set[,i]}
  else{ag <- c(ag, set[,i])}
}

c2 <- ag
i2 <- length(c2)
l2 <- rep("AFRA",i2)

#SR-LA
id <- match(SRLA, colnames(subset))
set <- subset[,id]
for (i in 1:dim(set)[2]){
  if(i==1){ag <- set[,i]}
  else{ag <- c(ag, set[,i])}
}

c3 <- ag 
i3 <- length(c3)
l3 <- rep("SRLA",i3)

#SR-RA
id <- match(SRRA, colnames(subset))
set <- subset[,id]
for (i in 1:dim(set)[2]){
  if(i==1){ag <- set[,i]}
  else{ag <- c(ag, set[,i])}
}

c4 <- ag 
i4 <- length(c4)
l4 <- rep("SRRA",i4)


df <- as.data.frame(cbind(c(c1,c2,c3,c4),c(l1,l2,l3,l4)))
colnames(df) <- c("fpkm", "sample_type")
df$sample_type <- as.factor(df$sample_type)
df$fpkm <- as.numeric(df$fpkm)
boxplot(fpkm ~ sample_type, data = df, main ="AF-LA")

#stats
wilcox.test(df$fpkm[df$sample_type=="AFLA"], df$fpkm[df$sample_type=="AFRA"], alternative="greater") # p-value < 2.2e-16
wilcox.test(df$fpkm[df$sample_type=="AFLA"], df$fpkm[df$sample_type=="SRLA"], alternative="greater") # p-value < 1.243e-07
wilcox.test(df$fpkm[df$sample_type=="AFLA"], df$fpkm[df$sample_type=="SRRA"], alternative="greater") # p-value < 2.2e-16
p_values <- c(2.2e-16, 1.243e-07, 2.2e-16)
p_adjusted <- p.adjust(p_values, method = "bonferroni")

#📊 Figure 2b (top left - AFLA)
ggplot(df, aes(x=sample_type, y=fpkm, fill=sample_type)) + 
  geom_violin(colour = "black", width=0.6) + scale_y_continuous(name = "FPKM") + 
  scale_fill_manual(values=c("#D85B74","#E3D081","#7298BC","#91C7B1")) +
  theme_classic() + geom_boxplot(colour = "black", fill="#F2FDFD", width=0.2,outlier.shape = NA) +
  labs(title="AF-LA enriched peakset") + 
  theme(axis.text.x = element_text(size = 14), axis.title.y = element_text(size = 16))  + 
  geom_signif(
    comparisons = list(c("AFLA", "AFRA"),c("AFLA", "SRLA"),c("AFLA", "SRRA")),
    y_position = c(12.5, 13.5, 14.5),
    map_signif_level = FALSE,
    annotations = "***"
  )



####################
### AF-RA peakset ##
####################

subset <- data2[data2$ids=="AFRA",]
boxplot(subset[,-21], las=2)#By library

#Aggregate data across same sample type: AFLA, AFRA, SRLA, SRRA

#AF-LA
id <- match(AFLA, colnames(subset))
set <- subset[,id]
for (i in 1:dim(set)[2]){
  if(i==1){ag <- set[,i]}
  else{ag <- c(ag, set[,i])}
}

c1 <- ag
i1 <- length(c1)
l1 <- rep("AFLA",i1)

#AF-RA
id <- match(AFRA, colnames(subset))
set <- subset[,id]
for (i in 1:dim(set)[2]){
  if(i==1){ag <- set[,i]}
  else{ag <- c(ag, set[,i])}
}

c2 <- ag
i2 <- length(c2)
l2 <- rep("AFRA",i2)

#SR-LA
id <- match(SRLA, colnames(subset))
set <- subset[,id]
for (i in 1:dim(set)[2]){
  if(i==1){ag <- set[,i]}
  else{ag <- c(ag, set[,i])}
}

c3 <- ag 
i3 <- length(c3)
l3 <- rep("SRLA",i3)

#SR-RA
id <- match(SRRA, colnames(subset))
set <- subset[,id]
for (i in 1:dim(set)[2]){
  if(i==1){ag <- set[,i]}
  else{ag <- c(ag, set[,i])}
}

c4 <- ag 
i4 <- length(c4)
l4 <- rep("SRRA",i4)


df <- as.data.frame(cbind(c(c1,c2,c3,c4),c(l1,l2,l3,l4)))
colnames(df) <- c("fpkm", "sample_type")
df$sample_type <- as.factor(df$sample_type)
df$fpkm <- as.numeric(df$fpkm)
boxplot(fpkm ~ sample_type, data = df, main ="AF-RA")


#stats
wilcox.test(df$fpkm[df$sample_type=="AFRA"], df$fpkm[df$sample_type=="AFLA"], alternative="greater") #p-value < 2.2e-16
wilcox.test(df$fpkm[df$sample_type=="AFRA"], df$fpkm[df$sample_type=="SRLA"], alternative="greater") #p-value < 2.2e-16
wilcox.test(df$fpkm[df$sample_type=="AFRA"], df$fpkm[df$sample_type=="SRRA"], alternative="greater") #p-value < 2.2e-16

p_values <- c(2.2e-16, 2.2e-16, 2.2e-16)
p_adjusted <- p.adjust(p_values, method = "bonferroni")

#📊 Figure 2b (top right - AFRA)
ggplot(df, aes(x=sample_type, y=fpkm, fill=sample_type)) + 
  geom_violin(colour = "black", width=0.6) + scale_y_continuous(name = "FPKM") + 
  scale_fill_manual(values=c("#D85B74","#E3D081","#7298BC","#91C7B1")) +
  theme_classic() + geom_boxplot(colour = "black", fill="#F2FDFD", width=0.2,outlier.shape = NA) +
  labs(title="AF-LA enriched peakset") + 
  theme(axis.text.x = element_text(size = 14), axis.title.y = element_text(size = 16))  + 
  geom_signif(
    comparisons = list(c("AFRA", "AFLA"),c("AFRA", "SRLA"),c("AFRA", "SRRA")),
    y_position = c(12.5, 12.5, 13.5),
    map_signif_level = FALSE,
    annotations = "***"
  )



####################
### SR-LA peakset ##
####################

subset <- data2[data2$ids=="SRLA",]
boxplot(subset[,-21], las=2)#By library

#Aggregate data across same sample type: AFLA, AFRA, SRLA, SRRA

#AF-LA
id <- match(AFLA, colnames(subset))
set <- subset[,id]
for (i in 1:dim(set)[2]){
  if(i==1){ag <- set[,i]}
  else{ag <- c(ag, set[,i])}
}

c1 <- ag
i1 <- length(c1)
l1 <- rep("AFLA",i1)

#AF-RA
id <- match(AFRA, colnames(subset))
set <- subset[,id]
for (i in 1:dim(set)[2]){
  if(i==1){ag <- set[,i]}
  else{ag <- c(ag, set[,i])}
}

c2 <- ag
i2 <- length(c2)
l2 <- rep("AFRA",i2)

#SR-LA
id <- match(SRLA, colnames(subset))
set <- subset[,id]
for (i in 1:dim(set)[2]){
  if(i==1){ag <- set[,i]}
  else{ag <- c(ag, set[,i])}
}

c3 <- ag 
i3 <- length(c3)
l3 <- rep("SRLA",i3)

#SR-RA
id <- match(SRRA, colnames(subset))
set <- subset[,id]
for (i in 1:dim(set)[2]){
  if(i==1){ag <- set[,i]}
  else{ag <- c(ag, set[,i])}
}

c4 <- ag 
i4 <- length(c4)
l4 <- rep("SRRA",i4)


df <- as.data.frame(cbind(c(c1,c2,c3,c4),c(l1,l2,l3,l4)))
colnames(df) <- c("fpkm", "sample_type")
df$sample_type <- as.factor(df$sample_type)
df$fpkm <- as.numeric(df$fpkm)
boxplot(fpkm ~ sample_type, data = df, main ="SR-LA")

#stats
wilcox.test(df$fpkm[df$sample_type=="SRLA"], df$fpkm[df$sample_type=="AFRA"], alternative="greater") #p-value < 2.2e-16
wilcox.test(df$fpkm[df$sample_type=="SRLA"], df$fpkm[df$sample_type=="AFLA"], alternative="greater") #p-value < 2.2e-16
wilcox.test(df$fpkm[df$sample_type=="SRLA"], df$fpkm[df$sample_type=="SRRA"], alternative="greater") #p-value < 2.2e-16


#📊 Figure 2b (bottom left - SRLA)
ggplot(df, aes(x=sample_type, y=fpkm, fill=sample_type)) + 
  geom_violin(colour = "black", width=0.6) + scale_y_continuous(name = "FPKM") + 
  scale_fill_manual(values=c("#D85B74","#E3D081","#7298BC","#91C7B1")) +
  theme_classic() + geom_boxplot(colour = "black", fill="#F2FDFD", width=0.2,outlier.shape = NA) +
  labs(title="SR-LA enriched peakset") + 
  theme(axis.text.x = element_text(size = 14), axis.title.y = element_text(size = 16))  + 
  geom_signif(
    comparisons = list(c("SRLA", "AFRA"),c("SRLA", "SRRA"),c("SRLA", "AFLA")),
    y_position = c(12.5, 12.5, 13.5),
    map_signif_level = FALSE,
    annotations = "***"
  )



####################
### SR-RA peakset ##
####################


subset <- data2[data2$ids=="SRRA",]
boxplot(subset[,-21], las=2)#By library

#Aggregate data across same sample type: AFLA, AFRA, SRLA, SRRA

#AF-LA
id <- match(AFLA, colnames(subset))
set <- subset[,id]
for (i in 1:dim(set)[2]){
  if(i==1){ag <- set[,i]}
  else{ag <- c(ag, set[,i])}
}

c1 <- ag
i1 <- length(c1)
l1 <- rep("AFLA",i1)

#AF-RA
id <- match(AFRA, colnames(subset))
set <- subset[,id]
for (i in 1:dim(set)[2]){
  if(i==1){ag <- set[,i]}
  else{ag <- c(ag, set[,i])}
}

c2 <- ag
i2 <- length(c2)
l2 <- rep("AFRA",i2)

#SR-LA
id <- match(SRLA, colnames(subset))
set <- subset[,id]
for (i in 1:dim(set)[2]){
  if(i==1){ag <- set[,i]}
  else{ag <- c(ag, set[,i])}
}

c3 <- ag 
i3 <- length(c3)
l3 <- rep("SRLA",i3)

#SR-RA
id <- match(SRRA, colnames(subset))
set <- subset[,id]
for (i in 1:dim(set)[2]){
  if(i==1){ag <- set[,i]}
  else{ag <- c(ag, set[,i])}
}

c4 <- ag 
i4 <- length(c4)
l4 <- rep("SRRA",i4)


df <- as.data.frame(cbind(c(c1,c2,c3,c4),c(l1,l2,l3,l4)))
colnames(df) <- c("fpkm", "sample_type")
df$sample_type <- as.factor(df$sample_type)
df$fpkm <- as.numeric(df$fpkm)
boxplot(fpkm ~ sample_type, data = df, main ="SR-RA")


#stats
wilcox.test(df$fpkm[df$sample_type=="SRRA"], df$fpkm[df$sample_type=="AFRA"], alternative="greater") #p-value = 3.858e-15
wilcox.test(df$fpkm[df$sample_type=="SRRA"], df$fpkm[df$sample_type=="AFLA"], alternative="greater") #p-value = 2.653e-05
wilcox.test(df$fpkm[df$sample_type=="SRRA"], df$fpkm[df$sample_type=="SRLA"], alternative="greater") #p-value = 8.044e-07

p_values <- c(3.858e-15, 2.653e-05, 8.044e-07)
p_adjusted <- p.adjust(p_values, method = "bonferroni")


#📊 Figure 2b (bottom right - SRRA)
ggplot(df, aes(x=sample_type, y=fpkm, fill=sample_type)) + 
  geom_violin(colour = "black", width=0.6) + scale_y_continuous(name = "FPKM") + 
  scale_fill_manual(values=c("#D85B74","#E3D081","#7298BC","#91C7B1")) +
  theme_classic() + geom_boxplot(colour = "black", fill="#F2FDFD", width=0.2,outlier.shape = NA) +
  labs(title="SR-LA enriched peakset") + 
  theme(axis.text.x = element_text(size = 14), axis.title.y = element_text(size = 16))  + 
  geom_signif(
    comparisons = list(c("SRRA", "SRLA"),c("SRRA", "AFRA"),c("SRRA", "AFLA")),
    y_position = c(11, 12, 13),
    map_signif_level = FALSE,
    annotations = "***"
  )


#################################################################################################################################################################

## Gene Ontology Analysis 

setwd("/Rodriguez_2025_AFepigenome/H3K27ac_enriched_regions_definition/geneToPeaks")


#GREAT region-to-gene association rule was used with default parameters to associate enriched regions to proximal genes (https://great.stanford.edu/great/public/html/)

#4 group peaksets
SRLA_enriched_peaks_genes <- read.delim(file="./geneToPeaks/SR_LA_enriched_nonoverlappingwithallsets_GenrichGreylistedFiltered.txt", header = TRUE, sep = "")
SRLA_enriched_peaks_genes_names <- SRLA_enriched_peaks_genes[,1]
SRRA_enriched_peaks_genes <- read.delim(file="./geneToPeaks/SR_RA_enriched_nonoverlappingwithallsets_GenrichGreylistedFiltered.txt", header = TRUE, sep = "")
SRRA_enriched_peaks_genes_names <- SRRA_enriched_peaks_genes[,1]
AFLA_enriched_peaks_genes <- read.delim(file="./geneToPeaks/AF_LA_enriched_nonoverlappingwithallsets_GenrichGreylistedFiltered.txt", header = TRUE, sep = "")
AFLA_enriched_peaks_genes_names <- AFLA_enriched_peaks_genes[,1]
AFRA_enriched_peaks_genes <- read.delim(file="./geneToPeaks/AF_RA_enriched_nonoverlappingwithallsets_GenrichGreylistedFiltered.txt", header = TRUE, sep = "")
AFRA_enriched_peaks_genes_names <- AFRA_enriched_peaks_genes[,1]



############## Using Genrich peaksets as background
##############


setwd("/Rodriguez_2025_AFepigenome/Genrich_RegRegions_definition/geneToPeaks")


AFLA_background_genes <- read.delim(file="AFLA_Genrich_genes.txt", header = FALSE, sep = "")
AFLA_background_genes_names <- AFLA_background_genes[,1]

AFRA_background_genes <- read.delim(file="AFRA_Genrich_genes.txt", header = FALSE, sep = "")
AFRA_background_genes_names <- AFRA_background_genes[,1]

SRRA_background_genes <- read.delim(file="SRRA_Genrich_genes_b.txt", header = FALSE, sep = "")
SRRA_background_genes_names <- SRRA_background_genes[,1]

SRLA_background_genes <- read.delim(file="SRLA_Genrich_genes_b.txt", header = FALSE, sep = "")
SRLA_background_genes_names <- SRLA_background_genes[,1]


setwd("/Rodriguez_2025_AFepigenome/H3K27ac_enriched_regions_definition/")



## SRLA

#BP
GO_results <- enrichGO(gene = SRLA_enriched_peaks_genes_names, OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", ont = "BP", universe = SRLA_background_genes_names)
GO_results_SRLA_BP <- as.data.frame(GO_results)
write.csv(GO_results_SRLA_BP,file="GO_results_SRLA_BP_background.csv")
fit <- plot(barplot(GO_results, showCategory = 20))
#MF
GO_results <- enrichGO(gene = SRLA_enriched_peaks_genes_names, OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", ont = "MF", universe = SRLA_background_genes_names)
GO_results_SRLA_MF <- as.data.frame(GO_results)
write.csv(GO_results_SRLA_MF,file="GO_results_SRLA_MF_background.csv")
fit <- plot(barplot(GO_results, showCategory = 20))
#CC
GO_results <- enrichGO(gene = SRLA_enriched_peaks_genes_names, OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", ont = "CC", universe = SRLA_background_genes_names)
GO_results_SRLA_CC <- as.data.frame(GO_results)
write.csv(GO_results_SRLA_CC,file="GO_results_SRLA_CC_background.csv")
fit <- plot(barplot(GO_results, showCategory = 20))


## SRRA

#BP
GO_results <- enrichGO(gene = SRRA_enriched_peaks_genes_names, OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", ont = "BP", universe = SRRA_background_genes_names)
GO_results_SRRA_BP <- as.data.frame(GO_results)
write.csv(GO_results_SRRA_BP,file="GO_results_SRRA_BP_background.csv")
fit <- plot(barplot(GO_results, showCategory = 20))
#MF
GO_results <- enrichGO(gene = SRRA_enriched_peaks_genes_names, OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", ont = "MF", universe = SRRA_background_genes_names)
GO_results_SRRA_MF <- as.data.frame(GO_results)
write.csv(GO_results_SRLA_MF,file="GO_results_SRRA_MF_background.csv")
fit <- plot(barplot(GO_results, showCategory = 20))
#CC
GO_results <- enrichGO(gene = SRRA_enriched_peaks_genes_names, OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", ont = "CC", universe = SRRA_background_genes_names)
GO_results_SRRA_CC <- as.data.frame(GO_results)
write.csv(GO_results_SRRA_CC,file="GO_results_SRRA_CC_background.csv")
fit <- plot(barplot(GO_results, showCategory = 9))


## AFLA

#BP
GO_results <- enrichGO(gene = AFLA_enriched_peaks_genes_names, OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", ont = "BP", universe = AFLA_background_genes_names)
GO_results_AFLA_BP <- as.data.frame(GO_results)
write.csv(GO_results_AFLA_BP,file="GO_results_AFLA_BP_background.csv")
fit <- plot(barplot(GO_results, showCategory = 20))
#MF
GO_results <- enrichGO(gene = AFLA_enriched_peaks_genes_names, OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", ont = "MF", universe = AFLA_background_genes_names)
GO_results_AFLA_MF <- as.data.frame(GO_results)
write.csv(GO_results_AFLA_MF,file="GO_results_AFLA_MF_background.csv")
fit <- plot(barplot(GO_results, showCategory = 20))
#CC
GO_results <- enrichGO(gene = AFLA_enriched_peaks_genes_names, OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", ont = "CC", universe = AFLA_background_genes_names)
GO_results_AFLA_CC <- as.data.frame(GO_results)
write.csv(GO_results_AFLA_CC,file="GO_results_AFLA_CC_background.csv")
fit <- plot(barplot(GO_results, showCategory = 20))




## AFRA

#BP
GO_results <- enrichGO(gene = AFRA_enriched_peaks_genes_names, OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", ont = "BP", universe = AFRA_background_genes_names)
GO_results_AFRA_BP <- as.data.frame(GO_results)
write.csv(GO_results_AFRA_BP,file="GO_results_AFRA_BP_background.csv")
fit <- plot(barplot(GO_results, showCategory = 20))
#MF
GO_results <- enrichGO(gene = AFRA_enriched_peaks_genes_names, OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", ont = "MF", universe = AFRA_background_genes_names)
GO_results_AFRA_MF <- as.data.frame(GO_results)
write.csv(GO_results_AFRA_MF,file="GO_results_AFRA_MF_background.csv")
fit <- plot(barplot(GO_results, showCategory = 20))
#CC
GO_results <- enrichGO(gene = AFRA_enriched_peaks_genes_names, OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", ont = "CC", universe = AFRA_background_genes_names)
GO_results_AFRA_CC <- as.data.frame(GO_results)
write.csv(GO_results_AFRA_CC,file="GO_results_AFRA_CC_background.csv")
fit <- plot(barplot(GO_results, showCategory = 20))


#heatmap plots
GO_results_AFLA_BP_background <- read.csv("GO_results_AFLA_BP_background.csv")
GO_results_AFRA_BP_background <- read.csv("GO_results_AFRA_BP_background.csv")
GO_results_SRRA_BP_background <- read.csv("GO_results_SRRA_BP_background.csv")
GO_results_SRLA_BP_background <- read.csv("GO_results_SRLA_BP_background.csv")

GO_results_AFRA_BP_background['group'] <- 'AF-RA'
GO_results_AFLA_BP_background['group'] <- 'AF-LA'
GO_results_SRRA_BP_background['group'] <- 'SR-RA'
GO_results_SRLA_BP_background['group'] <- 'SR-LA'

GO_results_AFLA_BP_background_GOs <- GO_results_AFLA_BP_background[,3]
GO_results_AFLA_BP_background_logpadj <- -log10(GO_results_AFLA_BP_background[,7])   
GO_results_AFLA_BP_background_group <- GO_results_AFLA_BP_background[,11]
GO_results_AFLA_BP_background_matrix <- cbind(GO_results_AFLA_BP_background_GOs,GO_results_AFLA_BP_background_logpadj,GO_results_AFLA_BP_background_group)

GO_results_AFRA_BP_background_GOs <- GO_results_AFRA_BP_background[,3]
GO_results_AFRA_BP_background_logpadj <- -log10(GO_results_AFRA_BP_background[,7])   
GO_results_AFRA_BP_background_group <- GO_results_AFRA_BP_background[,11]
GO_results_AFRA_BP_background_matrix <- cbind(GO_results_AFRA_BP_background_GOs,GO_results_AFRA_BP_background_logpadj,GO_results_AFRA_BP_background_group)

GO_results_SRLA_BP_background_GOs <- GO_results_SRLA_BP_background[,3]
GO_results_SRLA_BP_background_logpadj <- -log10(GO_results_SRLA_BP_background[,7])   
GO_results_SRLA_BP_background_group <- GO_results_SRLA_BP_background[,11]
GO_results_SRLA_BP_background_matrix <- cbind(GO_results_SRLA_BP_background_GOs,GO_results_SRLA_BP_background_logpadj,GO_results_SRLA_BP_background_group)

GO_results_SRRA_BP_background_GOs <- GO_results_SRRA_BP_background[,3]
GO_results_SRRA_BP_background_logpadj <- -log10(GO_results_SRRA_BP_background[,7])  
GO_results_SRRA_BP_background_group <- GO_results_SRRA_BP_background[,11]
GO_results_SRRA_BP_background_matrix <- cbind(GO_results_SRRA_BP_background_GOs,GO_results_SRRA_BP_background_logpadj,GO_results_SRRA_BP_background_group)


GO_BP_matrix <- rbind(GO_results_AFLA_BP_background_matrix,GO_results_AFRA_BP_background_matrix,GO_results_SRLA_BP_background_matrix,GO_results_SRRA_BP_background_matrix)
colnames(GO_BP_matrix) <- c("GO_term","Log_P-adj","Group")
GO_BP_matrix_df <- as.data.frame(GO_BP_matrix)
GO_BP_matrix_df_unmelted <- reshape(GO_BP_matrix_df, idvar = "GO_term", timevar = "Group", direction = "wide")
GO_BP_matrix_df_unmelted[is.na(GO_BP_matrix_df_unmelted)] <- 0
rownames(GO_BP_matrix_df_unmelted) <- GO_BP_matrix_df_unmelted[,1]
GO_BP_matrix_df_unmelted <- GO_BP_matrix_df_unmelted[,-1]
GO_BP_matrix_df_unmelted[] <- lapply(GO_BP_matrix_df_unmelted, function(x) as.numeric(as.character(x)))

GO_BP_matrix_df_unmelted_t <- t(GO_BP_matrix_df_unmelted)
GO_BP_matrix_df_unmelted_scaled <- apply(GO_BP_matrix_df_unmelted_t, 1, scale)
rownames(GO_BP_matrix_df_unmelted_scaled) <- colnames(GO_BP_matrix_df_unmelted_t)



############

#📊 Figure S2d (top left - AFLA)
GO_BP_matrix_df_unmelted
GO_BP_matrix_df_unmelted_scaled_AFLA <- GO_BP_matrix_df_unmelted_scaled[order(GO_BP_matrix_df_unmelted_scaled[,1], decreasing = TRUE), ]
GO_BP_matrix_df_unmelted_scaled_AFLA_top10 <- head(GO_BP_matrix_df_unmelted_scaled_AFLA,10)
rownames(GO_BP_matrix_df_unmelted_scaled_AFLA_top10) <- c("striated muscle tissue dev.","cardiac muscle tissue dev.","striated muscle cell diff.","heart process","muscle cell diff.","muscle tissue dev.","blood circulation","reg. of heart contraction","heart contraction","muscle system process")
GOs_BP_heatmap_AFLAtop10 <- pheatmap(t(GO_BP_matrix_df_unmelted_scaled_AFLA_top10),color = rev(hcl.colors(50,"PuBu")),fontsize = 15, cellwidth = 40, cellheight = 40, angle_col = 45, border_color = "black", cluster_rows = FALSE, treeheight_row = 0, treeheight_col = 0, labels_row = c("AF-LA","AF-RA","SR-LA","SR-RA"))


###

#📊 Figure S2d (top right - AFRA)
GO_BP_matrix_df_unmelted
GO_BP_matrix_df_unmelted_scaled_AFRA <- GO_BP_matrix_df_unmelted_scaled[order(GO_BP_matrix_df_unmelted_scaled[,2], decreasing = TRUE), ]
GO_BP_matrix_df_unmelted_scaled_AFRA_top10 <- head(GO_BP_matrix_df_unmelted_scaled_AFRA,10)
rownames(GO_BP_matrix_df_unmelted_scaled_AFRA_top10) <- c("muscle tissue dev.","cardiac muscle tissue dev.","striated muscle tissue dev.","striated muscle cell diff.","muscle cell diff.","cardiac muscle cell diff.","muscle cell dev.","cardiocyte diff.","ossification","extracellular matrix org.")
GOs_BP_heatmap_AFRAtop10 <- pheatmap(t(GO_BP_matrix_df_unmelted_scaled_AFRA_top10),color = rev(hcl.colors(50,"PuBu")),fontsize = 15, cellwidth = 40, cellheight = 40, angle_col = 45, border_color = "black", cluster_rows = FALSE, treeheight_row = 0, treeheight_col = 0, labels_row = c("AF-LA","AF-RA","SR-LA","SR-RA"))


####

#📊 Figure S2d (bottom right - SRRA)
GO_BP_matrix_df_unmelted
GO_BP_matrix_df_unmelted_scaled_SRRA <- GO_BP_matrix_df_unmelted_scaled[order(GO_BP_matrix_df_unmelted_scaled[,4], decreasing = TRUE), ]
#remove "regulation of cell communication by electrical coupling involved in cardiac conduction" for representation purposes - excessively long and biological already captured by other enrched ontology terms
GO_BP_matrix_df_unmelted_scaled_SRRA <- GO_BP_matrix_df_unmelted_scaled_SRRA[rownames(GO_BP_matrix_df_unmelted_scaled_SRRA) != "regulation of cell communication by electrical coupling involved in cardiac conduction", ]
GO_BP_matrix_df_unmelted_scaled_SRRA_top10 <- head(GO_BP_matrix_df_unmelted_scaled_SRRA,10)
rownames(GO_BP_matrix_df_unmelted_scaled_SRRA_top10) <- c("cell comm. involved in cardiac cond.","reg. of cell comm. by electr. coupling","cell comm. by electr. coupling (cardiac)","cell comm. by electr. coupling","heart process","heart contraction","cardiac conduction","cardiac chamber development","cardiac ventricle development", "non-canonical Wnt sign. pathway")
GOs_BP_heatmap_SRRAtop10 <- pheatmap(t(GO_BP_matrix_df_unmelted_scaled_SRRA_top10),color = rev(hcl.colors(50,"PuBu")),fontsize = 15, cellwidth = 40, cellheight = 40, angle_col = 45, border_color = "black", cluster_rows = FALSE, treeheight_row = 0, treeheight_col = 0, labels_row = c("AF-LA","AF-RA","SR-LA","SR-RA"))


####


#📊 Figure S2d (bottom left - SRLA)
GO_BP_matrix_df_unmelted
GO_BP_matrix_df_unmelted_scaled_SRLA <- GO_BP_matrix_df_unmelted_scaled[order(GO_BP_matrix_df_unmelted_scaled[,3], decreasing = TRUE), ]
GO_BP_matrix_df_unmelted_scaled_SRLA_top10 <- head(GO_BP_matrix_df_unmelted_scaled_SRLA,10)
rownames(GO_BP_matrix_df_unmelted_scaled_SRLA_top10) <- c("regulation of heart contraction","heart process","heart contraction","muscle system process","reg. of muscle system process","reg. of sodium ion TM transport","reg. of blood circulation","muscle contraction","cardiac conduction", "negative reg. of transerase activity")
GOs_BP_heatmap_SRLAtop10 <- pheatmap(t(GO_BP_matrix_df_unmelted_scaled_SRLA_top10),color = rev(hcl.colors(50,"PuBu")),fontsize = 15, cellwidth = 40, cellheight = 40, angle_col = 45, border_color = "black", cluster_rows = FALSE, treeheight_row = 0, treeheight_col = 0, labels_row = c("AF-LA","AF-RA","SR-LA","SR-RA"))




#################################################################################################################################################################


# Association of H3K27ac enrichment with differential gene expression across sample groups


setwd("/Rodriguez_2025_AFepigenome/H3K27ac_enriched_regions_definition/")

SRLA_enriched_peaks_genes <- read.delim(file="./geneToPeaks/SR_LA_enriched_nonoverlappingwithallsets_GenrichGreylistedFiltered.txt", header = TRUE, sep = "")
SRRA_enriched_peaks_genes <- read.delim(file="./geneToPeaks/SR_RA_enriched_nonoverlappingwithallsets_GenrichGreylistedFiltered.txt", header = TRUE, sep = "")
AFLA_enriched_peaks_genes <- read.delim(file="./geneToPeaks/AF_LA_enriched_nonoverlappingwithallsets_GenrichGreylistedFiltered.txt", header = TRUE, sep = "")
AFRA_enriched_peaks_genes <- read.delim(file="./geneToPeaks/AF_RA_enriched_nonoverlappingwithallsets_GenrichGreylistedFiltered.txt", header = TRUE, sep = "")


AF_enriched_peaks_genes <- read.delim(file="./geneToPeaks/AF_enriched_nonoverlappingwithallsets_GenrichGreylistedFiltered.txt", header = TRUE, sep = "")
SR_enriched_peaks_genes <- read.delim(file="./geneToPeaks/SR_enriched_nonoverlappingwithallsets_GenrichGreylistedFiltered.txt", header = TRUE, sep = "")


# Sets of DEGs
setwd("/Rodriguez_2025_AFepigenome/RNAseq_DGE_analysis/")

## NOTE: The input DEG files used below (e.g., files ending in "SRup", "AFup", "AFRAup", "SRRAup", etc.) were generated by subsetting the full DESeq2 output files according to the direction of differential expression.
## These filtered text files were obtained during the analysis performed in RNAseq_DGE_analysis/RNAseq_DE_Analysis.R

AFvsSR_SRup <- read.delim("DEseq2_results_SRup.txt", header = FALSE, sep = "")
AFvsSR_SRup <- AFvsSR_SRup[-1,]
AFvsSR_SRup_genes <- AFvsSR_SRup[,1]
AFvsSR_AFup <- read.delim("DEseq2_results_AFup.txt", header = FALSE, sep = "")
AFvsSR_AFup <- AFvsSR_AFup[-1,]
AFvsSR_AFup_genes <- AFvsSR_AFup[,1]
SRRAvsAFRA_AFRAup <- read.delim("DESeq2_results_SRRAvsAFRA_AFRAup.txt", header = FALSE, sep = "")
SRRAvsAFRA_AFRAup <- SRRAvsAFRA_AFRAup[-1,]
SRRAvsAFRA_AFRAup_genes <- SRRAvsAFRA_AFRAup[,1]
SRRAvsAFRA_SRRAup <- read.delim("DESeq2_results_SRRAvsAFRA_SRRAup.txt", header = FALSE, sep = "")
SRRAvsAFRA_SRRAup <- SRRAvsAFRA_SRRAup[-1,]
SRRAvsAFRA_SRRAup_genes <- SRRAvsAFRA_SRRAup[,1]
SRLAvsAFLA_SRLAup <- read.delim("DESeq2_results_SRLAvsAFLA_SRLAup.txt", header = FALSE, sep = "")
SRLAvsAFLA_SRLAup <- SRLAvsAFLA_SRLAup[-1,]
SRLAvsAFLA_SRLAup_genes <- SRLAvsAFLA_SRLAup[,1]
SRLAvsAFLA_AFLAup <- read.delim("DESeq2_results_SRLAvsAFLA_AFLAup.txt", header = FALSE, sep = "")
SRLAvsAFLA_AFLAup <- SRLAvsAFLA_AFLAup[-1,]
SRLAvsAFLA_AFLAup_genes <- SRLAvsAFLA_AFLAup[,1]


intersectt <- function (x, y) 
{
  if (is.null(x) || is.null(y)) 
    return(NULL)
  u <- as.vector(x)
  v <- as.vector(y)
  c(u[!duplicated(unclass(u)) & (match(u, v, 0L) > 0L)], v[numeric()])
}


setwd("/Rodriguez_2025_AFepigenome/H3K27ac_enriched_regions_definition/")

#Look at intersect of AF-enriched peaks associated genes and AFvsSR_AFup genes
AFvsSR_AFup_genes_AF_enriched_peaks_genes <- intersectt (AFvsSR_AFup_genes,AF_enriched_peaks_genes[,1])  
AFvsSR_AFup_genes_AF_enriched_peaks_genes
AFvsSR_AFup_genes_AF_enriched_peaks_genes_x <- list(
  A = AFvsSR_AFup_genes, 
  B = AF_enriched_peaks_genes[,1])
AFvsSR_AFup_genes_enriched_peaks_genes_VennDiagram <- ggvenn( # 📊 Figure S4e (right)
  AFvsSR_AFup_genes_AF_enriched_peaks_genes_x, 
  fill_color = c("#76e7cd", "#9b7ede" ),
  stroke_size = 0.5, set_name_size = 4
)
pdf("AFvsSR_AFup_genes_enriched_peaks_genes_VennDiagram.pdf")
AFvsSR_AFup_genes_enriched_peaks_genes_VennDiagram
dev.off()
AF_enriched_peaks_genes
AFvsSR_AFup_genes_enriched_peaks_genes_GREAT_details <- AF_enriched_peaks_genes[AF_enriched_peaks_genes$X. %in% AFvsSR_AFup_genes_AF_enriched_peaks_genes, ]
write.csv(AFvsSR_AFup_genes_enriched_peaks_genes_GREAT_details,file="AFvsSR_genes_AF_enriched_peaks.csv")


#Look at intersect of AF-enriched peaks associated genes and AFvsSR_SRup genes
AFvsSR_SRup_genes_AF_enriched_peaks_genes <- intersectt (AFvsSR_SRup_genes,AF_enriched_peaks_genes[,1])  
AFvsSR_SRup_genes_AF_enriched_peaks_genes
AFvsSR_SRup_genes_AF_enriched_peaks_genes_x <- list(
  A = AFvsSR_SRup_genes, 
  B = AF_enriched_peaks_genes[,1])
AFvsSR_SRup_genes_enriched_peaks_genes_VennDiagram <- ggvenn(
  AFvsSR_SRup_genes_AF_enriched_peaks_genes_x, 
  fill_color = c("#76e7cd", "#9b7ede" ),
  stroke_size = 0.5, set_name_size = 4
)
pdf("AFvsSR_SRup_genes_enriched_peaks_genes_VennDiagram.pdf")
AFvsSR_SRup_genes_enriched_peaks_genes_VennDiagram
dev.off()
AF_enriched_peaks_genes
AFvsSR_SRup_genes_AF_enriched_peaks_genes_GREAT_details <- AF_enriched_peaks_genes[AF_enriched_peaks_genes$X. %in% AFvsSR_SRup_genes_AF_enriched_peaks_genes, ]
write.csv(AFvsSR_SRup_genes_AF_enriched_peaks_genes_GREAT_details,file="AFvsSR_SRup_genes_AF_enriched_peaks_genes.csv")


#Look at intersect of SR-enriched peaks associated genes and AFvsSR_AFup genes
AFvsSR_AFup_genes_SR_enriched_peaks_genes <- intersectt (AFvsSR_AFup_genes,SR_enriched_peaks_genes[,1])  
AFvsSR_AFup_genes_SR_enriched_peaks_genes
AFvsSR_AFup_genes_SR_enriched_peaks_genes_x <- list(
  A = AFvsSR_AFup_genes, 
  B = SR_enriched_peaks_genes[,1])
AFvsSR_AFup_genes_SR_enriched_peaks_genes_VennDiagram <- ggvenn(
  AFvsSR_AFup_genes_SR_enriched_peaks_genes_x, 
  fill_color = c("#76e7cd", "#9b7ede" ),
  stroke_size = 0.5, set_name_size = 4
)
pdf("AFvsSR_AFup_genes_SR_enriched_peaks_genes_VennDiagram.pdf")
AFvsSR_AFup_genes_SR_enriched_peaks_genes_VennDiagram
dev.off()
SR_enriched_peaks_genes
AFvsSR_AFup_genes_SR_enriched_peaks_genes_GREAT_details <- SR_enriched_peaks_genes[SR_enriched_peaks_genes$X. %in% AFvsSR_AFup_genes_SR_enriched_peaks_genes, ]
write.csv(AFvsSR_AFup_genes_SR_enriched_peaks_genes_GREAT_details,file="AFvsSR_AFup_genes_SR_enriched_peaks_genes.csv")

#Look at intersect of SR-enriched peaks associated genes and AFvsSR_SRup genes
AFvsSR_SRup_genes_SR_enriched_peaks_genes <- intersectt (AFvsSR_SRup_genes,SR_enriched_peaks_genes[,1])  
AFvsSR_SRup_genes_SR_enriched_peaks_genes
AFvsSR_SRup_genes_SR_enriched_peaks_genes_x <- list(
  A = AFvsSR_SRup_genes, 
  B = SR_enriched_peaks_genes[,1])
AFvsSR_SRup_genes_SR_enriched_peaks_genes_VennDiagram <- ggvenn( # 📊 Figure S4e (left)
  AFvsSR_SRup_genes_SR_enriched_peaks_genes_x, 
  fill_color = c("#76e7cd", "#9b7ede" ),
  stroke_size = 0.5, set_name_size = 4
)
pdf("AFvsSR_SRup_genes_SR_enriched_peaks_genes_VennDiagram.pdf")
AFvsSR_SRup_genes_SR_enriched_peaks_genes_VennDiagram
dev.off()
SR_enriched_peaks_genes
AFvsSR_SRup_genes_SR_enriched_peaks_genes_GREAT_details <- SR_enriched_peaks_genes[SR_enriched_peaks_genes$X. %in% AFvsSR_SRup_genes_SR_enriched_peaks_genes, ]
write.csv(AFvsSR_SRup_genes_SR_enriched_peaks_genes_GREAT_details,file="AFvsSR_SRup_genes_SR_enriched_peaks_genes.csv")


#Look at intersect of AFRA-enriched peaks associated genes and SRRAvsAFRA_SRRAup genes
SRRAvsAFRA_SRRAup_genes_AFRA_enriched_peaks_genes <- intersectt (SRRAvsAFRA_SRRAup_genes,AFRA_enriched_peaks_genes[,1])  
SRRAvsAFRA_SRRAup_genes_AFRA_enriched_peaks_genes
SRRAvsAFRA_SRRAup_genes_AFRA_enriched_peaks_genes_x <- list(
  A = SRRAvsAFRA_SRRAup_genes, 
  B = AFRA_enriched_peaks_genes[,1])
SRRAvsAFRA_SRRAup_genes_AFRA_enriched_peaks_genes_VennDiagram <- ggvenn(
  SRRAvsAFRA_SRRAup_genes_AFRA_enriched_peaks_genes_x, 
  fill_color = c("#76e7cd", "#9b7ede" ),
  stroke_size = 0.5, set_name_size = 4
)
pdf("SRRAvsAFRA_SRRAup_genes_AFRA_enriched_peaks_genes_VennDiagram.pdf")
SRRAvsAFRA_SRRAup_genes_AFRA_enriched_peaks_genes_VennDiagram
dev.off()
AFRA_enriched_peaks_genes
SRRAvsAFRA_SRRAup_genes_AFRA_enriched_peaks_genes_GREAT_details <- AFRA_enriched_peaks_genes[AFRA_enriched_peaks_genes$X. %in% SRRAvsAFRA_SRRAup_genes_AFRA_enriched_peaks_genes, ]
write.csv(SRRAvsAFRA_SRRAup_genes_AFRA_enriched_peaks_genes_GREAT_details,file="SRRAvsAFRA_SRRAup_genes_AFRA_enriched_peaks_genes.csv")

#Look at intersect of AFRA-enriched peaks associated genes and SRRAvsAFRA_AFRAup genes
SRRAvsAFRA_AFRAup_genes_AFRA_enriched_peaks_genes <- intersectt (SRRAvsAFRA_AFRAup_genes,AFRA_enriched_peaks_genes[,1])  
SRRAvsAFRA_AFRAup_genes_AFRA_enriched_peaks_genes
SRRAvsAFRA_AFRAup_genes_AFRA_enriched_peaks_genes_x <- list(
  A = SRRAvsAFRA_AFRAup_genes, 
  B = AFRA_enriched_peaks_genes[,1])
SRRAvsAFRA_AFRAup_genes_AFRA_enriched_peaks_genes_VennDiagram <- ggvenn( # 📊 Figure 3b (bottom right)
  SRRAvsAFRA_AFRAup_genes_AFRA_enriched_peaks_genes_x, 
  fill_color = c("#76e7cd", "#9b7ede" ),
  stroke_size = 0.5, set_name_size = 4
)
pdf("SRRAvsAFRA_AFRAup_genes_AFRA_enriched_peaks_genes_VennDiagram.pdf")
SRRAvsAFRA_AFRAup_genes_AFRA_enriched_peaks_genes_VennDiagram
dev.off()
AFRA_enriched_peaks_genes
SRRAvsAFRA_AFRAup_genes_AFRA_enriched_peaks_genes_GREAT_details <- AFRA_enriched_peaks_genes[AFRA_enriched_peaks_genes$X. %in% SRRAvsAFRA_AFRAup_genes_AFRA_enriched_peaks_genes, ]
write.csv(SRRAvsAFRA_AFRAup_genes_AFRA_enriched_peaks_genes_GREAT_details,file="SRRAvsAFRA_AFRAup_genes_AFRA_enriched_peaks_genes.csv")


#Look at intersect of AFLA-enriched peaks associated genes and SRLAvsAFLA_SRLAup genes
SRLAvsAFLA_SRLAup_genes_AFLA_enriched_peaks_genes <- intersectt (SRLAvsAFLA_SRLAup_genes,AFLA_enriched_peaks_genes[,1])  
SRLAvsAFLA_SRLAup_genes_AFLA_enriched_peaks_genes
SRLAvsAFLA_SRLAup_genes_AFLA_enriched_peaks_genes_x <- list(
  A = SRLAvsAFLA_SRLAup_genes, 
  B = AFLA_enriched_peaks_genes[,1])
SRLAvsAFLA_SRLAup_genes_AFLA_enriched_peaks_genes_VennDiagram <- ggvenn(
  SRLAvsAFLA_SRLAup_genes_AFLA_enriched_peaks_genes_x, 
  fill_color = c("#76e7cd", "#9b7ede" ),
  stroke_size = 0.5, set_name_size = 4
)
pdf("SRLAvsAFLA_SRLAup_genes_AFLA_enriched_peaks_genes_VennDiagram.pdf")
SRLAvsAFLA_SRLAup_genes_AFLA_enriched_peaks_genes_VennDiagram
dev.off()
AFLA_enriched_peaks_genes
SRLAvsAFLA_SRLAup_genes_AFLA_enriched_peaks_genes_GREAT_details <- AFLA_enriched_peaks_genes[AFLA_enriched_peaks_genes$X. %in% SRLAvsAFLA_SRLAup_genes_AFLA_enriched_peaks_genes, ]
write.csv(SRLAvsAFLA_SRLAup_genes_AFLA_enriched_peaks_genes_GREAT_details,file="SRLAvsAFLA_SRLAup_genes_AFLA_enriched_peaks_genes.csv")


#Look at intersect of AFLA-enriched peaks associated genes and SRLAvsAFLA_AFLAup genes
SRLAvsAFLA_AFLAup_genes_AFLA_enriched_peaks_genes <- intersectt (SRLAvsAFLA_AFLAup_genes,AFLA_enriched_peaks_genes[,1])  
SRLAvsAFLA_AFLAup_genes_AFLA_enriched_peaks_genes
SRLAvsAFLA_AFLAup_genes_AFLA_enriched_peaks_genes_x <- list(
  A = SRLAvsAFLA_AFLAup_genes, 
  B = AFLA_enriched_peaks_genes[,1])
SRLAvsAFLA_AFLAup_genes_AFLA_enriched_peaks_genes_VennDiagram <- ggvenn( # 📊 Figure 3b (bottom left)
  SRLAvsAFLA_AFLAup_genes_AFLA_enriched_peaks_genes_x, 
  fill_color = c("#76e7cd", "#9b7ede" ),
  stroke_size = 0.5, set_name_size = 4
)
pdf("SRLAvsAFLA_AFLAup_genes_AFLA_enriched_peaks_genes_VennDiagram.pdf")
SRLAvsAFLA_AFLAup_genes_AFLA_enriched_peaks_genes_VennDiagram
dev.off()
AFLA_enriched_peaks_genes
SRLAvsAFLA_AFLAup_genes_AFLA_enriched_peaks_genes_GREAT_details <- AFLA_enriched_peaks_genes[AFLA_enriched_peaks_genes$X. %in% SRLAvsAFLA_AFLAup_genes_AFLA_enriched_peaks_genes, ]
write.csv(SRLAvsAFLA_AFLAup_genes_AFLA_enriched_peaks_genes_GREAT_details,file="SRLAvsAFLA_AFLAup_genes_AFLA_enriched_peaks_genes.csv")


#Look at intersect of SRRA-enriched peaks associated genes and SRRAvsAFRA_SRRAup genes
SRRAvsAFRA_SRRAup_genes_SRRA_enriched_peaks_genes <- intersectt (SRRAvsAFRA_SRRAup_genes,SRRA_enriched_peaks_genes[,1])  
SRRAvsAFRA_SRRAup_genes_SRRA_enriched_peaks_genes
SRRAvsAFRA_SRRAup_genes_SRRA_enriched_peaks_genes_x <- list(
  A = SRRAvsAFRA_SRRAup_genes, 
  B = SRRA_enriched_peaks_genes[,1])
SRRAvsAFRA_SRRAup_genes_SRRA_enriched_peaks_genes_VennDiagram <- ggvenn( # 📊 Figure 3b (top right)
  SRRAvsAFRA_SRRAup_genes_SRRA_enriched_peaks_genes_x, 
  fill_color = c("#76e7cd", "#9b7ede" ),
  stroke_size = 0.5, set_name_size = 4
)
pdf("SRRAvsAFRA_SRRAup_genes_SRRA_enriched_peaks_genes_VennDiagram.pdf")
SRRAvsAFRA_SRRAup_genes_SRRA_enriched_peaks_genes_VennDiagram
dev.off()
SRRA_enriched_peaks_genes
SRRAvsAFRA_SRRAup_genes_SRRA_enriched_peaks_genes_GREAT_details <- SRRA_enriched_peaks_genes[SRRA_enriched_peaks_genes$X. %in% SRRAvsAFRA_SRRAup_genes_SRRA_enriched_peaks_genes, ]
write.csv(SRRAvsAFRA_SRRAup_genes_SRRA_enriched_peaks_genes_GREAT_details,file="SRRAvsAFRA_SRRAup_genes_SRRA_enriched_peaks_genes_GREAT_details.csv")


#Look at intersect of SRRA-enriched peaks associated genes and SRRAvsAFRA_AFRAup genes
SRRAvsAFRA_AFRAup_genes_SRRA_enriched_peaks_genes <- intersectt (SRRAvsAFRA_AFRAup_genes,SRRA_enriched_peaks_genes[,1])  
SRRAvsAFRA_AFRAup_genes_SRRA_enriched_peaks_genes
SRRAvsAFRA_AFRAup_genes_SRRA_enriched_peaks_genes_x <- list(
  A = SRRAvsAFRA_AFRAup_genes, 
  B = SRRA_enriched_peaks_genes[,1])
SRRAvsAFRA_AFRAup_genes_SRRA_enriched_peaks_genes_VennDiagram <- ggvenn(
  SRRAvsAFRA_AFRAup_genes_SRRA_enriched_peaks_genes_x, 
  fill_color = c("#76e7cd", "#9b7ede" ),
  stroke_size = 0.5, set_name_size = 4
)
pdf("SRRAvsAFRA_AFRAup_genes_SRRA_enriched_peaks_genes_VennDiagram.pdf")
SRRAvsAFRA_AFRAup_genes_SRRA_enriched_peaks_genes_VennDiagram
dev.off()
SRRA_enriched_peaks_genes
SRRAvsAFRA_AFRAup_genes_SRRA_enriched_peaks_genes_GREAT_details <- SRRA_enriched_peaks_genes[SRRA_enriched_peaks_genes$X. %in% SRRAvsAFRA_AFRAup_genes_SRRA_enriched_peaks_genes, ]
write.csv(SRRAvsAFRA_AFRAup_genes_SRRA_enriched_peaks_genes_GREAT_details,file="SRRAvsAFRA_AFRAup_genes_SRRA_enriched_peaks_genes.csv")

#Look at intersect of SRLA-enriched peaks associated genes and SRLAvsAFLA_SRLAup genes
SRLAvsAFLA_SRLAup_genes_SRLA_enriched_peaks_genes <- intersectt (SRLAvsAFLA_SRLAup_genes,SRLA_enriched_peaks_genes[,1])  
SRLAvsAFLA_SRLAup_genes_SRLA_enriched_peaks_genes
SRLAvsAFLA_SRLAup_genes_SRLA_enriched_peaks_genes_x <- list(
  A = SRLAvsAFLA_SRLAup_genes, 
  B = SRLA_enriched_peaks_genes[,1])
SRLAvsAFLA_SRLAup_genes_SRLA_enriched_peaks_genes_VennDiagram <- ggvenn( # 📊 Figure 3b (top left)
  SRLAvsAFLA_SRLAup_genes_SRLA_enriched_peaks_genes_x, 
  fill_color = c("#76e7cd", "#9b7ede" ),
  stroke_size = 0.5, set_name_size = 4
)
pdf("SRLAvsAFLA_SRLAup_genes_SRLA_enriched_peaks_genes_VennDiagram.pdf")
SRLAvsAFLA_SRLAup_genes_SRLA_enriched_peaks_genes_VennDiagram
dev.off()
SRLA_enriched_peaks_genes
SRLAvsAFLA_SRLAup_genes_SRLA_enriched_peaks_genes_GREAT_details <- SRLA_enriched_peaks_genes[SRLA_enriched_peaks_genes$X. %in% SRLAvsAFLA_SRLAup_genes_SRLA_enriched_peaks_genes, ]
write.csv(SRLAvsAFLA_SRLAup_genes_SRLA_enriched_peaks_genes_GREAT_details,file="SRLAvsAFLA_SRLAup_genes_SRLA_enriched_peaks_genes.csv")


#Look at intersect of SRLA-enriched peaks associated genes and SRLAvsAFLA_AFLAup genes
SRLAvsAFLA_AFLAup_genes_SRLA_enriched_peaks_genes <- intersectt (SRLAvsAFLA_AFLAup_genes,SRLA_enriched_peaks_genes[,1])  
SRLAvsAFLA_AFLAup_genes_SRLA_enriched_peaks_genes
SRLAvsAFLA_AFLAup_genes_SRLA_enriched_peaks_genes_x <- list(
  A = SRLAvsAFLA_AFLAup_genes, 
  B = SRLA_enriched_peaks_genes[,1])
SRLAvsAFLA_AFLAup_genes_SRLA_enriched_peaks_genes_VennDiagram <- ggvenn(
  SRLAvsAFLA_AFLAup_genes_SRLA_enriched_peaks_genes_x, 
  fill_color = c("#76e7cd", "#9b7ede" ),
  stroke_size = 0.5, set_name_size = 4
)
pdf("SRLAvsAFLA_AFLAup_genes_SRLA_enriched_peaks_genes_VennDiagram.pdf")
SRLAvsAFLA_AFLAup_genes_SRLA_enriched_peaks_genes_VennDiagram
dev.off()
SRLA_enriched_peaks_genes
SRLAvsAFLA_AFLAup_genes_SRLA_enriched_peaks_genes_GREAT_details <- SRLA_enriched_peaks_genes[SRLA_enriched_peaks_genes$X. %in% SRLAvsAFLA_AFLAup_genes_SRLA_enriched_peaks_genes, ]
write.csv(SRLAvsAFLA_AFLAup_genes_SRLA_enriched_peaks_genes_GREAT_details,file="SRLAvsAFLA_AFLAup_genes_SRLA_enriched_peaks_genes.csv")


###### Fisher's Exact Test


#AFvsSR_AFup_genes_AF_enriched_peaks_genes
a <- 58
b <- 1052-58
c <- 502-58
d <- 18370-a-b-c
AFvsSR_AFup_genes_AF_enriched_peaks_genes_df <- data.frame(
  "Group1" = c(a,b),
  "Group2" = c(c,d),
  row.names = c("Cat1", "Cat2"),
  stringsAsFactors = FALSE
)
AFvsSR_AFup_genes_AF_enriched_peaks_genes_df <- fisher.test(AFvsSR_AFup_genes_AF_enriched_peaks_genes_df)
AFvsSR_AFup_genes_AF_enriched_peaks_genes_df
AFvsSR_AFup_genes_AF_enriched_peaks_genes_df$p.value


#AFvsSR_AFup_genes_SR_enriched_peaks_genes
a <- 23
b <- 1148-23
c <- 502-23
d <- 18370-a-b-c
AFvsSR_AFup_genes_SR_enriched_peaks_genes_df <- data.frame(
  "Group1" = c(a,b),
  "Group2" = c(c,d),
  row.names = c("Cat1", "Cat2"),
  stringsAsFactors = FALSE
)
AFvsSR_AFup_genes_SR_enriched_peaks_genes_df <- fisher.test(AFvsSR_AFup_genes_SR_enriched_peaks_genes_df)
AFvsSR_AFup_genes_SR_enriched_peaks_genes_df
AFvsSR_AFup_genes_SR_enriched_peaks_genes_df$p.value


#AFvsSR_SRup_genes_AF_enriched_peaks_genes
a <- 38
b <- 1052-38
c <- 744-38
d <- 18370-a-b-c
AFvsSR_SRup_genes_AF_enriched_peaks_genes_df <- data.frame(
  "Group1" = c(a,b),
  "Group2" = c(c,d),
  row.names = c("Cat1", "Cat2"),
  stringsAsFactors = FALSE
)
AFvsSR_SRup_genes_AF_enriched_peaks_genes_df <- fisher.test(AFvsSR_SRup_genes_AF_enriched_peaks_genes_df)
AFvsSR_SRup_genes_AF_enriched_peaks_genes_df
AFvsSR_SRup_genes_AF_enriched_peaks_genes_df$p.value


#AFvsSR_SRup_genes_SR_enriched_peaks_genes
a <- 110
b <- 1148-110
c <- 744-110
d <- 18370-a-b-c
AFvsSR_SRup_genes_SR_enriched_peaks_genes_df <- data.frame(
  "Group1" = c(a,b),
  "Group2" = c(c,d),
  row.names = c("Cat1", "Cat2"),
  stringsAsFactors = FALSE
)
AFvsSR_SRup_genes_SR_enriched_peaks_genes_df <- fisher.test(AFvsSR_SRup_genes_SR_enriched_peaks_genes_df)
AFvsSR_SRup_genes_SR_enriched_peaks_genes_df
AFvsSR_SRup_genes_SR_enriched_peaks_genes_df$p.value


#######


#SRRAvsAFRA_SRRAup_genes_AFRA_enriched_peaks_genes
a <- 5
b <- 1100-5
c <- 164-5
d <- 18370-a-b-c
SRRAvsAFRA_SRRAup_genes_AFRA_enriched_peaks_genes_df <- data.frame(
  "Group1" = c(a,b),
  "Group2" = c(c,d),
  row.names = c("Cat1", "Cat2"),
  stringsAsFactors = FALSE
)
SRRAvsAFRA_SRRAup_genes_AFRA_enriched_peaks_genes_df <- fisher.test(SRRAvsAFRA_SRRAup_genes_AFRA_enriched_peaks_genes_df)
SRRAvsAFRA_SRRAup_genes_AFRA_enriched_peaks_genes_df
SRRAvsAFRA_SRRAup_genes_AFRA_enriched_peaks_genes_df$p.value


#SRRAvsAFRA_AFRAup_genes_AFRA_enriched_peaks_genes
a <- 10
b <- 1100-10
c <- 61-10
d <- 18370-a-b-c
SRRAvsAFRA_AFRAup_genes_AFRA_enriched_peaks_genes_df <- data.frame(
  "Group1" = c(a,b),
  "Group2" = c(c,d),
  row.names = c("Cat1", "Cat2"),
  stringsAsFactors = FALSE
)
SRRAvsAFRA_AFRAup_genes_AFRA_enriched_peaks_genes_df <- fisher.test(SRRAvsAFRA_AFRAup_genes_AFRA_enriched_peaks_genes_df)
SRRAvsAFRA_AFRAup_genes_AFRA_enriched_peaks_genes_df
SRRAvsAFRA_AFRAup_genes_AFRA_enriched_peaks_genes_df$p.value


#SRLAvsAFLA_SRLAup_genes_AFLA_enriched_peaks_genes
a <- 9
b <- 1024-9
c <- 118-9
d <- 18370-a-b-c
SRLAvsAFLA_SRLAup_genes_AFLA_enriched_peaks_genes_df <- data.frame(
  "Group1" = c(a,b),
  "Group2" = c(c,d),
  row.names = c("Cat1", "Cat2"),
  stringsAsFactors = FALSE
)
SRLAvsAFLA_SRLAup_genes_AFLA_enriched_peaks_genes_df <- fisher.test(SRLAvsAFLA_SRLAup_genes_AFLA_enriched_peaks_genes_df)
SRLAvsAFLA_SRLAup_genes_AFLA_enriched_peaks_genes_df
SRLAvsAFLA_SRLAup_genes_AFLA_enriched_peaks_genes_df$p.value


#SRLAvsAFLA_AFLAup_genes_AFLA_enriched_peaks_genes
a <- 9
b <- 1024-9
c <- 89-9
d <- 18370-a-b-c
SRLAvsAFLA_AFLAup_genes_AFLA_enriched_peaks_genes_df <- data.frame(
  "Group1" = c(a,b),
  "Group2" = c(c,d),
  row.names = c("Cat1", "Cat2"),
  stringsAsFactors = FALSE
)
SRLAvsAFLA_SRLAup_genes_AFLA_enriched_peaks_genes_df <- fisher.test(SRLAvsAFLA_AFLAup_genes_AFLA_enriched_peaks_genes_df)
SRLAvsAFLA_SRLAup_genes_AFLA_enriched_peaks_genes_df
SRLAvsAFLA_SRLAup_genes_AFLA_enriched_peaks_genes_df$p.value


#SRRAvsAFRA_SRRAup_genes_SRRA_enriched_peaks_genes
a <- 21
b <- 1146-a
c <- 164-a
d <- 18370-a-b-c
SRRAvsAFRA_SRRAup_genes_SRRA_enriched_peaks_genes_df <- data.frame(
  "Group1" = c(a,b),
  "Group2" = c(c,d),
  row.names = c("Cat1", "Cat2"),
  stringsAsFactors = FALSE
)
SRRAvsAFRA_SRRAup_genes_SRRA_enriched_peaks_genes_df <- fisher.test(SRRAvsAFRA_SRRAup_genes_SRRA_enriched_peaks_genes_df)
SRRAvsAFRA_SRRAup_genes_SRRA_enriched_peaks_genes_df
SRRAvsAFRA_SRRAup_genes_SRRA_enriched_peaks_genes_df$p.value



#SRRAvsAFRA_AFRAup_genes_SRRA_enriched_peaks_genes
a <- 5
b <- 1146-a
c <- 61-a
d <- 18370-a-b-c
SRRAvsAFRA_AFRAup_genes_SRRA_enriched_peaks_genes_df <- data.frame(
  "Group1" = c(a,b),
  "Group2" = c(c,d),
  row.names = c("Cat1", "Cat2"),
  stringsAsFactors = FALSE
)
SRRAvsAFRA_AFRAup_genes_SRRA_enriched_peaks_genes_df <- fisher.test(SRRAvsAFRA_AFRAup_genes_SRRA_enriched_peaks_genes_df)
SRRAvsAFRA_AFRAup_genes_SRRA_enriched_peaks_genes_df
SRRAvsAFRA_AFRAup_genes_SRRA_enriched_peaks_genes_df$p.value


#SRLAvsAFLA_SRLAup_genes_SRLA_enriched_peaks_genes
a <- 19
b <- 1069-a
c <- 118-a
d <- 18370-a-b-c
SRLAvsAFLA_SRLAup_genes_SRLA_enriched_peaks_genes_df <- data.frame(
  "Group1" = c(a,b),
  "Group2" = c(c,d),
  row.names = c("Cat1", "Cat2"),
  stringsAsFactors = FALSE
)
SRLAvsAFLA_SRLAup_genes_SRLA_enriched_peaks_genes_df <- fisher.test(SRLAvsAFLA_SRLAup_genes_SRLA_enriched_peaks_genes_df)
SRLAvsAFLA_SRLAup_genes_SRLA_enriched_peaks_genes_df
SRLAvsAFLA_SRLAup_genes_SRLA_enriched_peaks_genes_df$p.value


#SRLAvsAFLA_AFLAup_genes_SRLA_enriched_peaks_genes
a <- 4
b <- 1069-a
c <- 89-a
d <- 18370-a-b-c
SRLAvsAFLA_AFLAup_genes_SRLA_enriched_peaks_genes_df <- data.frame(
  "Group1" = c(a,b),
  "Group2" = c(c,d),
  row.names = c("Cat1", "Cat2"),
  stringsAsFactors = FALSE
)
SRLAvsAFLA_AFLAup_genes_SRLA_enriched_peaks_genes_df <- fisher.test(SRLAvsAFLA_AFLAup_genes_SRLA_enriched_peaks_genes_df)
SRLAvsAFLA_AFLAup_genes_SRLA_enriched_peaks_genes_df
SRLAvsAFLA_AFLAup_genes_SRLA_enriched_peaks_genes_df$p.value


#### Fisher correlation plots ####


setwd("/Rodriguez_2025_AFepigenome/H3K27ac_enriched_regions_definition/")

fihersexacttest_rnaseqchipseqolaps_subsetting_AFvsSR=read.csv("rnaseqchipseqolaps_subsetting_AFvsSR_log.csv", header=TRUE)
fihersexacttest_rnaseqchipseqolaps_subsetting_AFvsSR
fihersexacttest_rnaseqchipseqolaps_subsetting_AFvsSR <- melt(fihersexacttest_rnaseqchipseqolaps_subsetting_AFvsSR)




DEGs_sets_order <- c("LA - AF up","LA - SR up","RA - AF up","RA - SR up")

# 📊 Figure 3a
overlap_heatmap <- ggplot(data = fihersexacttest_rnaseqchipseqolaps_subsetting, aes(x=factor(X,level=DEGs_sets_order), y=factor(variable, level=peaksets_orderforheatmaps),
                                                                                    fill=value)) +
  geom_tile(colour = "black", linewidth = 1) + scale_fill_gradientn(colours=c("#FCFFE1", "#A8DABD", "#58BAB0", "#265D98")) + theme(axis.text=element_text(size=10),axis.text.x=element_text(angle = -90, hjust = 0))
overlap_heatmap +scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0))


# 📊 Figure S4d (right) - as above but with colour scale to match AFvsSR heatmap
overlap_heatmap <- ggplot(data = fihersexacttest_rnaseqchipseqolaps_subsetting, aes(x=factor(X,level=DEGs_sets_order), y=factor(variable, level=peaksets_orderforheatmaps),
                                                                                    fill=value)) +
  geom_tile(colour = "black", linewidth = 1) + scale_fill_gradientn(colours=c("#FCFFE1", "#A8DABD", "#58BAB0", "#265D98"),limits = c(0, 17)) + theme(axis.text=element_text(size=10),axis.text.x=element_text(angle = -90, hjust = 0))
overlap_heatmap +scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0))



DEGs_sets_order_AFvsSR <- c("AFvsSR_AFup","AFvsSR_SRup")

# 📊 Figure S4d (left)
overlap_heatmap <- ggplot(data = fihersexacttest_rnaseqchipseqolaps_subsetting_AFvsSR, aes(x=factor(X,level=DEGs_sets_order_AFvsSR), y=factor(variable, level=peaksets_orderforheatmaps),
                                                                                           fill=value)) +
  geom_tile(colour = "black", linewidth = 1) + scale_fill_gradientn(colours=c("#FCFFE1", "#A8DABD", "#58BAB0", "#265D98"),limits = c(0, 17)) + theme(axis.text=element_text(size=10),axis.text.x=element_text(angle = -90, hjust = 0))

overlap_heatmap +scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) 




#################################################################################################################################################################


# Association of H3K27ac enrichment with differentially methylated regions (DMRs)







