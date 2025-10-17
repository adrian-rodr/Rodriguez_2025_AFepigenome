




####################################################
####################################################
############ Differential Analysis of ##############
############ ChIP-seq data            ##############
####################################################
####################################################


library(TCseq)
library(DiffBind)
library(dplyr)
library(ggplot2)


##########
# Differential enrichment analysis with TCseq
##########



setwd("Rodriguez_2025_AFepigenome/ChIP-seq_basic_processing")

s12b_SR_LA_peaks.narrowPeak <- read.csv('12b_SR_LA_peaks.narrowPeak', sep = '\t', header = FALSE)
s12b_SR_RA_peaks.narrowPeak <- read.csv('12b_SR_RA_peaks.narrowPeak', sep = '\t', header = FALSE)
s18_AF_RA_peaks.narrowPeak<- read.csv('18_AF_RA_peaks.narrowPeak', sep = '\t', header = FALSE)
s20_SR_LA_peaks.narrowPeak <- read.csv('20_SR_LA_peaks.narrowPeak', sep = '\t', header = FALSE)
s20_SR_RA_peaks.narrowPeak <- read.csv('20_SR_RA_peaks.narrowPeak', sep = '\t', header = FALSE)
s2231_AF_LA_peaks.narrowPeak <- read.csv('2231_AF_LA_peaks.narrowPeak', sep = '\t', header = FALSE)
s2231_AF_RA_peaks.narrowPeak <- read.csv('2231_AF_RA_peaks.narrowPeak', sep = '\t', header = FALSE)
s25_AF_RA_peaks.narrowPeak <- read.csv('25_AF_RA_peaks.narrowPeak', sep = '\t', header = FALSE)
s25_AF_LA_peaks.narrowPeak <- read.csv('25_AF_LA_peaks.narrowPeak', sep = '\t', header = FALSE)
s27_SR_LA_peaks.narrowPeak <- read.csv('27_SR_LA_peaks.narrowPeak', sep = '\t', header = FALSE)
s27_SR_RA_peaks.narrowPeak <- read.csv('27_SR_RA_peaks.narrowPeak', sep = '\t', header = FALSE)
s28_AF_RA_peaks.narrowPeak <- read.csv('28_AF_RA_peaks.narrowPeak', sep = '\t', header = FALSE)
s28_AF_LA_peaks.narrowPeak <- read.csv('28_AF_LA_peaks.narrowPeak', sep = '\t', header = FALSE)
s313_AF_LA_peaks.narrowPeak <- read.csv('313_AF_LA_peaks.narrowPeak', sep = '\t', header = FALSE)
s313_AF_RA_peaks.narrowPeak <- read.csv('313_AF_RA_peaks.narrowPeak', sep = '\t', header = FALSE)
s31_AF_LA_peaks.narrowPeak <- read.csv('31_AF_LA_peaks.narrowPeak', sep = '\t', header = FALSE)
s31_AF_RA_peaks.narrowPeak <- read.csv('31_AF_RA_peaks.narrowPeak', sep = '\t', header = FALSE)
s33_SR_RA_peaks.narrowPeak <- read.csv('33_SR_RA_peaks.narrowPeak', sep = '\t', header = FALSE)
s45_SR_LA_peaks.narrowPeak <- read.csv('45_SR_LA_peaks.narrowPeak', sep = '\t', header = FALSE)
s47_AF_LA_peaks.narrowPeak <- read.csv('47_AF_LA_peaks.narrowPeak', sep = '\t', header = FALSE)


narrowPeak_df_Oct24 <- bind_rows(s12b_SR_LA_peaks.narrowPeak,
                                 s12b_SR_RA_peaks.narrowPeak,
                                 s18_AF_RA_peaks.narrowPeak,
                                 s20_SR_LA_peaks.narrowPeak, 
                                 s20_SR_RA_peaks.narrowPeak, 
                                 s2231_AF_LA_peaks.narrowPeak, 
                                 s2231_AF_RA_peaks.narrowPeak,
                                 s25_AF_RA_peaks.narrowPeak,
                                 s25_AF_LA_peaks.narrowPeak,
                                 s27_SR_LA_peaks.narrowPeak,
                                 s27_SR_RA_peaks.narrowPeak, 
                                 s28_AF_RA_peaks.narrowPeak, 
                                 s28_AF_LA_peaks.narrowPeak, 
                                 s313_AF_LA_peaks.narrowPeak, 
                                 s313_AF_RA_peaks.narrowPeak, 
                                 s31_AF_LA_peaks.narrowPeak, 
                                 s31_AF_RA_peaks.narrowPeak, 
                                 s33_SR_RA_peaks.narrowPeak, 
                                 s45_SR_LA_peaks.narrowPeak, 
                                 s47_AF_LA_peaks.narrowPeak, 
)


library("TCseq")
gf_Oct24 <- peakreference (data = narrowPeak_df_Oct24)

save(gf_Oct24 ,file="gf_Oct24.RData")

load("gf_Oct24.RData")


experiment_BAMfiles_H3K27ac_AFvsSR <- read.csv(file='FINAL_DATASET_TCSeq_H3K27ac_OCT24.txt', sep='\t', header = TRUE)
library("TCseq")
tca_AFvsSR_H3K27ac_Oct24 <- TCA(design = experiment_BAMfiles_H3K27ac_AFvsSR, genomicFeature = gf_Oct24)
tca_AFvsSR_H3K27ac_Oct24
tca_AFvsSR_H3K27ac_Oct24 <- countReads(tca_AFvsSR_H3K27ac_Oct24)
tca_AFvsSR_H3K27ac_Oct24
save(tca_AFvsSR_H3K27ac_Oct24, file = "tca_AFvsSR_H3K27ac_Oct24.RData")


# Explore AF vs SR regions (differential analysis)
library("TCseq")
load("tca_AFvsSR_H3K27ac_Oct24.RData")
tca_AFvsSR_H3K27ac_Oct24
tca_AFvsSR_H3K27ac_Oct24
tca_AFvsSR_H3K27ac_Oct24.DB <- DBresult(tca_AFvsSR_H3K27ac_Oct24, group1 = "AF", group2 = "SR") ## Outputs the AFvsSR differentially bound regions (as shown in Supplementary Table 2).
#Volcano plot of AF vs SR regions
library(Repitools)
tca_AFvsSR_H3K27ac_Oct24.DB_GRanges <- unlist(tca_AFvsSR_H3K27ac_Oct24.DB)
tca_AFvsSR_H3K27ac_Oct24.DB_GRanges_df <- annoGR2DF(tca_AFvsSR_H3K27ac_Oct24.DB_GRanges)
library(dplyr)
tca_AFvsSR_H3K27ac_Oct24.DB_GRanges_df = mutate(tca_AFvsSR_H3K27ac_Oct24.DB_GRanges_df, UPORDOWN = ifelse(tca_AFvsSR_H3K27ac_Oct24.DB_GRanges_df$paj>0.05,"non", ifelse(tca_AFvsSR_H3K27ac_Oct24.DB_GRanges_df$logFC < 0,"AF","SR")))
AFvsSR_volcanoplot_tcseq <- ggplot(tca_AFvsSR_H3K27ac_Oct24.DB_GRanges_df, aes(logFC, -log(paj,10)))  +
  geom_point(aes(col = UPORDOWN), size = 1) +
  xlab(expression("log"[2]*"FC")) + 
  ylab(expression("-log"[10]*"P-adj")) +
  scale_color_manual(values = c("#F2977D", "lightgrey", "#3FAAC0")) +
  guides(colour = guide_legend(override.aes = list(size=1.5))) + theme(panel.background = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
AFvsSR_volcanoplot_tcseq


##########
# Differential enrichment analysis with DiffBind
##########


library("DiffBind")
#Constructs a DBA (Differential binding affinity) object from a sample sheet
AF_final_dataset_H3K27ac_Dec23_dba <- dba(sampleSheet="AF_final_dataset_H3K27ac_Dec23.csv")
AF_final_dataset_H3K27ac_Dec23_dba
#
AF_final_dataset_H3K27ac_Nov23_dba <- dba.count(AF_final_dataset_H3K27ac_Nov23_dba, summits=FALSE)
AF_final_dataset_H3K27ac_Nov23_dba
#
save(AF_final_dataset_H3K27ac_Dec23_dba, file="AF_final_dataset_H3K27ac_Dec23_dba.RData")
load("AF_final_dataset_H3K27ac_Dec23_dba.RData")
dba.peakset(AF_final_dataset_H3K27ac_Dec23_dba)
AF_final_dataset_H3K27ac_Dec23_dba
AF_final_dataset_H3K27ac_Dec23_dba <- dba(AF_final_dataset_H3K27ac_Dec23_dba)


#TISSUE
AF_final_dataset_H3K27ac_Dec23_dba <- dba.contrast(AF_final_dataset_H3K27ac_Dec23_dba, categories=DBA_TISSUE, minMembers=2 )
AF_final_dataset_H3K27ac_Dec23_dba <- dba.analyze(AF_final_dataset_H3K27ac_Dec23_dba,bBlacklist=FALSE,bGreylist=FALSE)
AF_final_dataset_H3K27ac_Dec23_dba.DB <- dba.report(AF_final_dataset_H3K27ac_Dec23_dba)
sum(AF_final_dataset_H3K27ac_Dec23_dba.DB$Fold>0) #638
sum(AF_final_dataset_H3K27ac_Dec23_dba.DB$Fold<0) #246
save(AF_final_dataset_H3K27ac_Dec23_dba.DB, file="AF_final_dataset_H3K27ac_Dec23_dba.DB") ## Outputs the AFvsSR differentially bound regions (as shown in 📎 Supplementary Table 2).
write.csv(AF_final_dataset_H3K27ac_Dec23_dba.DB, file="DiffLA")

#CONDITION
AF_final_dataset_H3K27ac_Dec23_dba_condition <- dba.contrast(AF_final_dataset_H3K27ac_Dec23_dba_condition, categories=DBA_CONDITION, minMembers=2 )
AF_final_dataset_H3K27ac_Dec23_dba_condition <- dba.analyze(AF_final_dataset_H3K27ac_Dec23_dba_condition,bBlacklist=FALSE,bGreylist=FALSE)
AF_final_dataset_H3K27ac_Dec23_dba_condition.DB <- dba.report(AF_final_dataset_H3K27ac_Dec23_dba_condition)
sum(AF_final_dataset_H3K27ac_Dec23_dba_condition.DB$Fold>0) #0
sum(AF_final_dataset_H3K27ac_Dec23_dba_condition.DB$Fold<0) #3
save(AF_final_dataset_H3K27ac_Dec23_dba_condition.DB, file="AF_final_dataset_H3K27ac_Dec23_dba_condition.DB") ## Outputs the AFvsSR differentially bound regions (as shown in 📎 Supplementary Table 2).
write.csv(AF_final_dataset_H3K27ac_Dec23_dba_condition.DB, file="AFvsSR_final_dataset_H3K27ac_Dec23_dba_condition.DB.csv")



