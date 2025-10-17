


#################################################
########### TF enrichment analysis ##############
#################################################

library(ggplot2)
library(pheatmap)
library(ReMapEnrich)
library(rtracklayer)

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### Representation of HOMER TF enrichment analysis results  ### 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 


## Generation of plots to represent data obtained from HOMER_enrichment_analysis.sh

## H3K27ac ChIP-seq data

AF_RA_TFBmotifs_knownResults <- read.delim("AFRA_MotifOutput_GenrichBackground/knownResults.txt",header=TRUE) 
AF_LA_TFBmotifs_knownResults <- read.delim("AFLA_MotifOutput_GenrichBackground/knownResults.txt",header=TRUE)
SR_RA_TFBmotifs_knownResults <- read.delim("SRRA_MotifOutput_GenrichBackground/knownResults.txt",header=TRUE)
SR_LA_TFBmotifs_knownResults <- read.delim("SRLA_MotifOutput_GenrichBackground/knownResults.txt",header=TRUE)

AF_RA_TFBmotifs_knownResults['group'] <- 'AF-RA'
AF_LA_TFBmotifs_knownResults['group'] <- 'AF-LA'
SR_RA_TFBmotifs_knownResults['group'] <- 'SR-RA'
SR_LA_TFBmotifs_knownResults['group'] <- 'SR-LA'
SR_TFBmotifs_knownResults['group'] <- 'SR'
AF_TFBmotifs_knownResults['group'] <- 'AF'

#use p values for representation in heatmap 
AF_RA_TFBmotifs_knownResults_motifs <- AF_RA_TFBmotifs_knownResults[,1]
AF_RA_TFBmotifs_knownResults_logpval <- -log10(AF_RA_TFBmotifs_knownResults[,3])
AF_RA_TFBmotifs_knownResults_group <- AF_RA_TFBmotifs_knownResults[,10]
AF_RA_TFBmotifs_matrix <- cbind(AF_RA_TFBmotifs_knownResults_motifs,AF_RA_TFBmotifs_knownResults_group,AF_RA_TFBmotifs_knownResults_logpval)

AF_LA_TFBmotifs_knownResults_motifs <- AF_LA_TFBmotifs_knownResults[,1]
AF_LA_TFBmotifs_knownResults_logpval <- -log10(AF_LA_TFBmotifs_knownResults[,3])
AF_LA_TFBmotifs_knownResults_group <- AF_LA_TFBmotifs_knownResults[,10]
AF_LA_TFBmotifs_matrix <- cbind(AF_LA_TFBmotifs_knownResults_motifs,AF_LA_TFBmotifs_knownResults_group,AF_LA_TFBmotifs_knownResults_logpval)

SR_RA_TFBmotifs_knownResults_motifs <- SR_RA_TFBmotifs_knownResults[,1]
SR_RA_TFBmotifs_knownResults_logpval <- -log10(SR_RA_TFBmotifs_knownResults[,3])
SR_RA_TFBmotifs_knownResults_group <- SR_RA_TFBmotifs_knownResults[,10]
SR_RA_TFBmotifs_matrix <- cbind(SR_RA_TFBmotifs_knownResults_motifs,SR_RA_TFBmotifs_knownResults_group,SR_RA_TFBmotifs_knownResults_logpval)

SR_LA_TFBmotifs_knownResults_motifs <- SR_LA_TFBmotifs_knownResults[,1]
SR_LA_TFBmotifs_knownResults_logpval <- -log10(SR_LA_TFBmotifs_knownResults[,3])
SR_LA_TFBmotifs_knownResults_group <- SR_LA_TFBmotifs_knownResults[,10]
SR_LA_TFBmotifs_matrix <- cbind(SR_LA_TFBmotifs_knownResults_motifs,SR_LA_TFBmotifs_knownResults_group,SR_LA_TFBmotifs_knownResults_logpval)

TFBmotifs_matrix <- rbind(AF_RA_TFBmotifs_matrix,AF_LA_TFBmotifs_matrix,SR_RA_TFBmotifs_matrix,SR_LA_TFBmotifs_matrix)
colnames(TFBmotifs_matrix) <- c("motifs","group","logpval")
TFBmotifs_matrix_df <- as.data.frame(TFBmotifs_matrix)
TFBmotifs_matrix_unmelted <- reshape(TFBmotifs_matrix_df, idvar = "motifs", timevar = "group", direction = "wide")

#AFRA sign motifs  #Sign pval < 10^-39
AFRA_sign_motifs<-TFBmotifs_matrix_unmelted[TFBmotifs_matrix_unmelted$`logpval.AF-RA` > 90,]
#AFLA sign motifs
AFLA_sign_motifs<-TFBmotifs_matrix_unmelted[TFBmotifs_matrix_unmelted$`logpval.AF-LA` > 90,]
#SRRA sign motifs
SRRA_sign_motifs<-TFBmotifs_matrix_unmelted[TFBmotifs_matrix_unmelted$`logpval.SR-RA` > 90,]
#SRLA sign motifs
SRLA_sign_motifs<-TFBmotifs_matrix_unmelted[TFBmotifs_matrix_unmelted$`logpval.SR-LA` > 90,]

combined_unique <- unique(c(AFRA_sign_motifs$motifs, AFLA_sign_motifs$motifs, SRRA_sign_motifs$motifs, SRLA_sign_motifs$motifs))

TFBmotifs_matrix_unmelted_forheatmap <- subset(TFBmotifs_matrix_unmelted, motifs %in% combined_unique)
dim(TFBmotifs_matrix_unmelted_forheatmap)

rownames(TFBmotifs_matrix_unmelted_forheatmap) <- TFBmotifs_matrix_unmelted_forheatmap[,1]
TFBmotifs_matrix_unmelted_forheatmap <- TFBmotifs_matrix_unmelted_forheatmap[,-1]
TFBmotifs_matrix_unmelted_forheatmap[] <- lapply(TFBmotifs_matrix_unmelted_forheatmap, function(x) as.numeric(as.character(x)))

# z-score normalization within each group
TFBmotifs_matrix_unmelted_forheatmap_t <- t(TFBmotifs_matrix_unmelted_forheatmap)
TFBmotifs_matrix_unmelted_forheatmap_scaled <- apply(TFBmotifs_matrix_unmelted_forheatmap_t, 1, scale)
rownames(TFBmotifs_matrix_unmelted_forheatmap_scaled) <- colnames(TFBmotifs_matrix_unmelted_forheatmap_t)
TFBmotifs_matrix_heatmap <- pheatmap(TFBmotifs_matrix_unmelted_forheatmap_scaled,color = hcl.colors(50, "Turku"),fontsize = 7, cellwidth = 10, cellheight = 10)
TFBmotifs_matrix_heatmap <- pheatmap(t(TFBmotifs_matrix_unmelted_forheatmap_scaled),color = rev(hcl.colors(50, "BuPu")),fontsize = 12, cellwidth = 16, cellheight = 16, border_color = "white", angle_col = 45)

# Simplify names fro TFs to improve readability of the plot
TFs_names <- read.csv("TFs_simplified_names_Genrichbackground.csv", header=FALSE)
TFs_names <- TFs_names[,2]
rownames(TFBmotifs_matrix_unmelted_forheatmap_scaled) <- TFs_names
TFBmotifs_matrix_heatmap <- pheatmap(t(TFBmotifs_matrix_unmelted_forheatmap_scaled),color = rev(hcl.colors(50, "BuPu")),fontsize = 12, cellwidth = 16, cellheight = 16, border_color = "white", angle_col = 45)
#change column names 
TFBmotifs_matrix_heatmap


# 📊 Figure S3a
TFBmotifs_matrix_heatmap <- pheatmap(t(TFBmotifs_matrix_unmelted_forheatmap_scaled),color = rev(hcl.colors(50, "BuPu")),fontsize = 12, cellwidth = 16, cellheight = 16, border_color = "white", angle_col = 45, cluster_rows = FALSE, cluster_cols = FALSE)
new_order <- c("logpval.AF-LA", "logpval.AF-RA", "logpval.SR-LA","logpval.SR-RA")
TFBmotifs_matrix_unmelted_forheatmap_scaled <- TFBmotifs_matrix_unmelted_forheatmap_scaled[, new_order]
TFBmotifs_matrix_heatmap <- pheatmap(t(TFBmotifs_matrix_unmelted_forheatmap_scaled),color = rev(hcl.colors(50, "BuPu")),fontsize = 12, cellwidth = 16, cellheight = 16, border_color = "white", angle_col = 45, cluster_rows = FALSE, cluster_cols = TRUE)



## AFvsSR DMRs from EPICv2 array data

AF_Hypermeth_TFBmotifs <- read.delim("DMRs_hypermeth/homerResults.txt",header=TRUE)
AF_Hypometh_TFBmotifs <- read.delim("DMRs_hypometh/homerResults.txt",header=TRUE)
AF_Hypermeth_TFBmotifs['group'] <- 'Hypermethylated_DMRs'
AF_Hypometh_TFBmotifs['group'] <- 'Hypomethylated_DRMs'

AF_Hypometh_TFBmotifs$Best.Match <- factor(AF_Hypometh_TFBmotifs$Best.Match, levels = rev(AF_Hypometh_TFBmotifs$Best.Match))
AF_Hypermeth_TFBmotifs$Best.Match <- factor(AF_Hypermeth_TFBmotifs$Best.Match, levels = rev(AF_Hypermeth_TFBmotifs$Best.Match))

# Extract y limits from global data
color_limits <- range(-(AF_TFBmotifs$Log.P.value ))
size_limits  <- range(AF_TFBmotifs$Perc_of_Targets)

# Define a fixed colour scale
df1 <- AF_TFBmotifs %>% filter(Best.Match %in% AF_Hypermeth_TFBmotifs$Best.Match)
df2 <- AF_TFBmotifs %>% filter(Best.Match %in% AF_Hypometh_TFBmotifs$Best.Match)


# 📊 Figure S7a - HOMER de novo TFBS enrichment analysis over DMRs (AF-hypermethylated and AF-hypomethylated)


ggplot(df1, aes(y = Best.Match, x = 1, size = Perc_of_Targets, fill = -(Log.P.value))) +
  geom_point(shape = 21, color = "black") +
  scale_size(range = c(5, 20), name = "% Regions", limits = size_limits) +
  scale_fill_gradient(low = "lightblue", high = "darkblue", name = "-log P-value", limits = color_limits) +
  theme_minimal() +
  #scale_size_continuous(limits = size_limits) +
  #scale_colour_gradient(limits = color_limits, low = "lightblue", high = "darkblue") +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_blank()
  ) 


ggplot(df2, aes(y = Best.Match, x = 1, size = Perc_of_Targets, fill = -(Log.P.value))) +
  geom_point(shape = 21, color = "black") +
  scale_size(range = c(5, 20), name = "% Regions", limits = size_limits) +
  scale_fill_gradient(low = "lightblue", high = "darkblue", name = "-log P-value", limits = color_limits) +
  theme_minimal() +
  #scale_size_continuous(limits = size_limits) +
  #scale_colour_gradient(limits = color_limits, low = "lightblue", high = "darkblue") +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_blank()
  ) 



### ### ### ### ### ### 
### ReMap Analysis  ### 
### ### ### ### ### ###

library(ReMapEnrich)

#Load the remap catalog, downloaded from: https://remap.univ-amu.fr/download_page
remapCatalog_hg38_cardiacbiotype <- bedToGranges("ReMap2022_hg38_cardiacbiotype.bed") #This is a curated selection of cardiac TF binding sites derived from publicly available ChIP-seq experimental data in relevant cardiac-related samples (including human cardiac tissue or stem-cell derived models)

#Load chromosome sizes from UCSC
ch <- downloadUcscChromSizes("hg38")




## H3K27ac ChIP-seq data

# NOTE: generate bed files for enriched sets in "Rodriguez_2025_AFepigenome/H3K27ac_enriched_regions_definition/" folder for each group with chr prefix in first column to make it compatible with ReMap input
# E.g.: sed 's/^/chr/' AF_LA_enriched_nonoverlappingwithallsets_GenrichGreylistedFiltered.bed > AF_LA_enriched_nonoverlappingwithallsets_GenrichGreylistedFiltered_chr.bed

AFLA_peaks <- bedToGranges("Rodriguez_2025_AFepigenome/H3K27ac_enriched_regions_definition/AF_LA_enriched_nonoverlappingwithallsets_GenrichGreylistedFiltered_chr.bed")
AFRA_peaks <- bedToGranges("Rodriguez_2025_AFepigenome/H3K27ac_enriched_regions_definition/AF_RA_enriched_nonoverlappingwithallsets_GenrichGreylistedFiltered_chr.bed")
SRRA_peaks <- bedToGranges("Rodriguez_2025_AFepigenome/H3K27ac_enriched_regions_definition/SR_RA_enriched_nonoverlappingwithallsets_GenrichGreylistedFiltered_chr.bed")
SRLA_peaks <- bedToGranges("Rodriguez_2025_AFepigenome/H3K27ac_enriched_regions_definition/SR_LA_enriched_nonoverlappingwithallsets_GenrichGreylistedFiltered_chr.bed")


res_AFLA <- enrichment(AFLA_peaks, remapCatalog_hg38_cardiacbiotype, chromSizes = ch, shuffles=10, nCores = 5)
res_AFRA <- enrichment(AFRA_peaks, remapCatalog_hg38_cardiacbiotype, chromSizes = ch, shuffles=10, nCores = 5)
res_SRLA <- enrichment(SRLA_peaks, remapCatalog_hg38_cardiacbiotype, chromSizes = ch, shuffles=10, nCores = 5)
res_SRRA <- enrichment(SRRA_peaks, remapCatalog_hg38_cardiacbiotype, chromSizes = ch, shuffles=10, nCores = 5)
res_AF <- enrichment(AF_peaks, remapCatalog_hg38_cardiacbiotype, chromSizes = ch, shuffles=10, nCores = 5)
res_SR <- enrichment(SR_peaks, remapCatalog_hg38_cardiacbiotype, chromSizes = ch, shuffles=10, nCores = 5)


res_AFLA <- subset(res_AFLA, q.value<0.05)
res_AFRA <- subset(res_AFRA, q.value<0.05)
res_SRLA <- subset(res_SRLA, q.value<0.05)
res_SRRA <- subset(res_SRRA, q.value<0.05)

res_AFLA['group'] <- 'AF-LA'
res_AFRA['group'] <- 'AF-RA'
res_SRLA['group'] <- 'SR-LA'
res_SRRA['group'] <- 'SR-RA'

res_AFLA_names <- res_AFLA[,1]
res_AFLA_qsign <- res_AFLA[,8]
res_AFLA_nbolaps <- res_AFLA[,2]
res_AFLA_group <- res_AFLA[,12]
AFLA_matrix <- cbind(res_AFLA_names,res_AFLA_qsign,res_AFLA_nbolaps,res_AFLA_group)

res_AFRA_names <- res_AFRA[,1]
res_AFRA_qsign <- res_AFRA[,8]
res_AFRA_nbolaps <- res_AFRA[,2]
res_AFRA_group <- res_AFRA[,12]
AFRA_matrix <- cbind(res_AFRA_names,res_AFRA_qsign,res_AFRA_nbolaps,res_AFRA_group)

res_SRLA_names <- res_SRLA[,1]
res_SRLA_qsign <- res_SRLA[,8]
res_SRLA_nbolaps <- res_SRLA[,2]
res_SRLA_group <- res_SRLA[,12]
SRLA_matrix <- cbind(res_SRLA_names,res_SRLA_qsign,res_SRLA_nbolaps,res_SRLA_group)

res_SRRA_names <- res_SRRA[,1]
res_SRRA_qsign <- res_SRRA[,8]
res_SRRA_nbolaps <- res_SRRA[,2]
res_SRRA_group <- res_SRRA[,12]
SRRA_matrix <- cbind(res_SRRA_names,res_SRRA_qsign,res_SRRA_nbolaps,res_SRRA_group)

TFdatasets_matrix <- rbind(AFRA_matrix,AFLA_matrix,SRRA_matrix,SRLA_matrix)
colnames(TFdatasets_matrix) <- c("names","qsign","nbolaps","group")
TFdatasets_matrix_df <- as.data.frame(TFdatasets_matrix)

# Remove rows with these datasets (they have NA values)
#GSE96691.CTCF.cardiomyocyte_CTCF_KO
#GSE60699.HEY2.CM7-1_cardiomyocyte
#GSE97761.TCF7L2.cardiac-ventricle_Wntup

TFdatasets_matrix_df <- subset(TFdatasets_matrix_df, names != "GSE96691.CTCF.cardiomyocyte_CTCF_KO")
TFdatasets_matrix_df <- subset(TFdatasets_matrix_df, names != "GSE60699.HEY2.CM7-1_cardiomyocyte")
TFdatasets_matrix_df <- subset(TFdatasets_matrix_df, names != "GSE97761.TCF7L2.cardiac-ventricle_Wntup")

TFdatasets_matrix_df[, c("qsign", "nbolaps")] <- lapply(TFdatasets_matrix_df[, c("qsign", "nbolaps")], as.numeric)

# Change names to simplify for plots representation of TF ChIP-seq datasets
TFdatasets_matrix_df$names
original_names <- TFdatasets_matrix_df$names[1:33]
keys <- original_names
values <- c("ETS1_CM_d9",
            "CTCF_LV_a",
            "CTCF_heart_a",
            "MED1_CM5",
            "TBX5_CM1",
            "TBX5_CM",
            "ETS1_CM_d14",
            "ETS1_CM_d0",
            "CTCF_LV_b",
            "MED1_CM7",
            "CTCF_LV_c",
            "ETS1_CM_d2",
            "ESRRG_CM",
            "ETS1_CM_d5",
            "TBX5_CM5",
            "CTCF_heart_b",
            "CTCF_cardiac-a",
            "CTCF_cardiac-b",
            "TBX5_7",
            "GATA4_CM5",
            "GATA4_CM7",
            "CTCF_Cardiac_fb",
            "MED1_CM1",
            "MED1_CM",
            "CTCF_RA",
            "CTCF_heart_c",
            "CTCF_heart_d",
            "GATA4_CM1",
            "GATA4_CM",
            "CTCF_heart_e",
            "HEY2_CM7_Dox",
            "HEY1_CM_Dox",
            "ISL1_CM")

# Assign keys to values
named_vector <- setNames(values, keys)

TFdatasets_matrix_df$names <- named_vector[match(TFdatasets_matrix_df$names, names(named_vector))]


# 📊 Figure S3b
ggplot(data = TFdatasets_matrix_df, aes(x = names, y = group, 
                                        color = qsign, size = nbolaps)) + 
  geom_point() +
  scale_color_gradient(low = "red", high = "blue") +
  theme_bw() + 
  ylab("") + 
  xlab("") + 
  ggtitle("ChIP-seq TFs enrichment analysis") + theme(axis.text.x = element_text(angle = 45, hjust = 1))






## AFvsSR DMRs from EPICv2 array data

DMRs_hypermethyalted <- bedToGranges("dmrs_dmrff_betanorm_Chris_DMRs_hypermeth.bed")
DMRs_hypomethylated <- bedToGranges("dmrs_dmrff_betanorm_Chris_DMRs_hypometh.bed")

# Use whole genome as background
# We are analyzing only a small set of regions due to limited available data for cardiac-relevant samples, thuse using collecive peaks (Genrich sets) as background would substantially reduce statistical power. 


res_DMRs_hypermethylated <- enrichment(DMRs_hypermethyalted, remapCatalog_hg38_cardiacbiotype, chromSizes = ch, shuffles=10, nCores = 5)
res_DMRs_hypomethylated <- enrichment(DMRs_hypomethylated, remapCatalog_hg38_cardiacbiotype, chromSizes = ch, shuffles=10, nCores = 5)


res_DMRs_hypermethylated<- na.omit(res_DMRs_hypermethylated)
res_DMRs_hypomethylated<- na.omit(res_DMRs_hypomethylated)


res_DMRs_hypermethylated['group'] <- 'DMRs_hypermeth'
res_DMRs_hypomethylated['group'] <- 'DMRs_hypometh'


res_DMRs_hypermethylated_names <- res_DMRs_hypermethylated[,1]
res_DMRs_hypermethylated_qsign <- res_DMRs_hypermethylated[,8]
res_DMRs_hypermethylated_nbolaps <- res_DMRs_hypermethylated[,2]
res_DMRs_hypermethylated_group <- res_DMRs_hypermethylated[,12]
DMRs_hypermethylated_matrix <- cbind(res_DMRs_hypermethylated_names,res_DMRs_hypermethylated_qsign,res_DMRs_hypermethylated_nbolaps,res_DMRs_hypermethylated_group)

res_DMRs_hypomethylated_names <- res_DMRs_hypomethylated[,1]
res_DMRs_hypomethylated_qsign <- res_DMRs_hypomethylated[,8]
res_DMRs_hypomethylated_nbolaps <- res_DMRs_hypomethylated[,2]
res_DMRs_hypomethylated_group <- res_DMRs_hypomethylated[,12]
DMRs_hypomethylated_matrix <- cbind(res_DMRs_hypomethylated_names,res_DMRs_hypomethylated_qsign,res_DMRs_hypomethylated_nbolaps,res_DMRs_hypomethylated_group)


TFdatasets_matrix <- rbind(DMRs_hypermethylated_matrix,DMRs_hypomethylated_matrix)
colnames(TFdatasets_matrix) <- c("names","qsign","nbolaps","group")
TFdatasets_matrix_df <- as.data.frame(TFdatasets_matrix)

res_DMRs_hypermethylated_sign_BindingSites <- res_DMRs_hypermethylated[res_DMRs_hypermethylated$q.value <0.05,]
res_DMRs_hypomethylated_sign_BindingSites <- res_DMRs_hypomethylated[res_DMRs_hypomethylated$q.value <0.05,]


combined_unique <- unique(c(rownames(res_DMRs_hypermethylated_sign_BindingSites), rownames(res_DMRs_hypomethylated_sign_BindingSites)))


TFdatasets_matrix_df_forheatmap <- subset(TFdatasets_matrix_df, names %in% combined_unique)
dim(TFdatasets_matrix_df_forheatmap)


TFdatasets_matrix_df_forheatmap[, c("qsign", "nbolaps")] <- lapply(TFdatasets_matrix_df_forheatmap[, c("qsign", "nbolaps")], as.numeric)


d<-ggplot(data = TFdatasets_matrix_df_forheatmap, aes(x = names, y = group, 
                                                      color = qsign, size = nbolaps)) + 
  geom_point() +
  scale_color_gradient(low = "red", high = "blue") +
  theme_bw() + 
  ylab("") + 
  xlab("") + 
  ggtitle("ChIP-seq TFs enrichment analysis") + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + coord_fixed(ratio = 2) 
d

# simplify names of datasets to improve readibility of plot
simplifiednames <- c("CTCF - Cardiac FB",
                     "CTCF - Cardiac muscle",
                     "CTCF - RA",
                     "CTCF - Heart_a",
                     "CTCF - Heart_b",
                     "CTCF - Heart_c",
                     "CTCF - LV",
                     "CTCF - Heart_d",
                     "CTCF - Cardiac muscle",
                     "CTCF - LV_b",
                     "CTCF - Heart_e",
                     "CTCF - LV_c",
                     "ESRRG - iPSC-CM",
                     "ETS1 - hESC-CM_d0",
                     "ETS1 - hESC-CM_d14",
                     "ETS1 - hESC-CM_d2",
                     "ETS1 - hESC-CM_d5",
                     "ETS1 - hESC-CM_d9",
                     "MED1 - iPSC-CM",
                     "TBX5 - iPSC-CM",
                     "TBX5 - iPSC-CM_d1",
                     "TBX5 - iPSC-CM_d5",
                     "TBX5 - iPSC-CM_d7"
)

length(simplifiednames)

# 📊 Figure S7d
d<-ggplot(data = TFdatasets_matrix_df_forheatmap_saved2, aes(x = names, y = group, 
                                                             color = qsign, size = nbolaps)) + 
  geom_point() +
  scale_color_gradient(low = "red", high = "blue") +
  theme_bw() + 
  ylab("") + 
  xlab("") + 
  ggtitle("ChIP-seq TFs enrichment analysis") + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + coord_fixed(ratio = 2) +  scale_x_discrete(labels = simplifiednames)
d








### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
###  MEME - FIMO - AF-associated TFs motif occurrences across DMRs  ### 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 

# The FIMO online tool from the MEME suite (v.5.5.7) (match p-value < 0.001) was used to annotate motifs over DMRs and the background set of regions. 
# Binding motifs for TFs of interest were obtained from JASPAR database - download frequency matrix in MEME-compatible format. 
# A hypergeometric test is used to calculate statistical significance of ocurrences across hypermethyalted and hypomethylated DMRs 


AFassocTFsoverHypometh_0001 <- read.table(file = 'AFassocTFsoverHypometh_0001.tsv', sep = '\t', header = TRUE)
AFassocTFsoverHypermeth_0001 <- read.table(file = 'AFassocTFsoverHypermeth_0001.tsv', sep = '\t', header = TRUE)
AFassocTFsoverbackground <- read.table(file = 'AFassocTFsoverbackground.tsv', sep = '\t', header = TRUE)

ETV1_overHypometh <- AFassocTFsoverHypometh_0001[AFassocTFsoverHypometh_0001$motif_alt_id == 'ETV1',]
TBX3_overHypometh <- AFassocTFsoverHypometh_0001[AFassocTFsoverHypometh_0001$motif_alt_id == 'TBX3',]
TBX5_overHypometh <- AFassocTFsoverHypometh_0001[AFassocTFsoverHypometh_0001$motif_alt_id == 'TBX5',]
NKX25_overHypometh <- AFassocTFsoverHypometh_0001[AFassocTFsoverHypometh_0001$motif_alt_id == 'NKX2-5',]
CREB2_overHypometh <- AFassocTFsoverHypometh_0001[AFassocTFsoverHypometh_0001$motif_alt_id == 'ATF4',]
PITX2_overHypometh <- AFassocTFsoverHypometh_0001[AFassocTFsoverHypometh_0001$motif_alt_id == 'PITX2',]
PRRX1_overHypometh <- AFassocTFsoverHypometh_0001[AFassocTFsoverHypometh_0001$motif_alt_id == 'PRRX1',]
HAND2_overHypometh <- AFassocTFsoverHypometh_0001[AFassocTFsoverHypometh_0001$motif_alt_id == 'HAND2',]
ZFHX3_overHypometh <- read.table(file = 'ZFHX3_over_hypometh.tsv', sep = '\t', header = TRUE)

ETV1_overHypometh <- ETV1_overHypometh[ , c("sequence_name", "start", "stop")] #46
TBX3_overHypometh <- TBX3_overHypometh[ , c("sequence_name", "start", "stop")] #41
TBX5_overHypometh <- TBX5_overHypometh[ , c("sequence_name", "start", "stop")] #50
NKX25_overHypometh <- NKX25_overHypometh[ , c("sequence_name", "start", "stop")] #22
CREB2_overHypometh <- CREB2_overHypometh[ , c("sequence_name", "start", "stop")] #10
PITX2_overHypometh <- PITX2_overHypometh[ , c("sequence_name", "start", "stop")] #16
PRRX1_overHypometh <- PRRX1_overHypometh[ , c("sequence_name", "start", "stop")] #19
HAND2_overHypometh <- HAND2_overHypometh[ , c("sequence_name", "start", "stop")] #29
ZFHX3_overHypometh <- ZFHX3_overHypometh[ , c("sequence_name", "start", "stop")] #6

ETV1_overHypermeth <- AFassocTFsoverHypermeth_0001[AFassocTFsoverHypermeth_0001$motif_alt_id == 'ETV1',]
TBX3_overHypermeth <- AFassocTFsoverHypermeth_0001[AFassocTFsoverHypermeth_0001$motif_alt_id == 'TBX3',]
TBX5_overHypermeth <- AFassocTFsoverHypermeth_0001[AFassocTFsoverHypermeth_0001$motif_alt_id == 'TBX5',]
NKX25_overHypermeth <- AFassocTFsoverHypermeth_0001[AFassocTFsoverHypermeth_0001$motif_alt_id == 'NKX2-5',]
CREB2_overHypermeth <- AFassocTFsoverHypermeth_0001[AFassocTFsoverHypermeth_0001$motif_alt_id == 'ATF4',]
PITX2_overHypermeth <- AFassocTFsoverHypermeth_0001[AFassocTFsoverHypermeth_0001$motif_alt_id == 'PITX2',]
PRRX1_overHypermeth <- AFassocTFsoverHypermeth_0001[AFassocTFsoverHypermeth_0001$motif_alt_id == 'PRRX1',]
HAND2_overHypermeth <- AFassocTFsoverHypermeth_0001[AFassocTFsoverHypermeth_0001$motif_alt_id == 'HAND2',]
ZFHX3_overHypermeth <- read.table(file = 'ZFHX3_over_hypermeth.tsv', sep = '\t', header = TRUE)

ETV1_overHypermeth <- ETV1_overHypermeth[ , c("sequence_name", "start", "stop")] #33
TBX3_overHypermeth <- TBX3_overHypermeth[ , c("sequence_name", "start", "stop")] #68
TBX5_overHypermeth <- TBX5_overHypermeth[ , c("sequence_name", "start", "stop")] #82
NKX25_overHypermeth <- NKX25_overHypermeth[ , c("sequence_name", "start", "stop")] #39
CREB2_overHypermeth <- CREB2_overHypermeth[ , c("sequence_name", "start", "stop")] #17
PITX2_overHypermeth <- PITX2_overHypermeth[ , c("sequence_name", "start", "stop")] #20
PRRX1_overHypermeth <- PRRX1_overHypermeth[ , c("sequence_name", "start", "stop")] #6
HAND2_overHypermeth <- HAND2_overHypermeth[ , c("sequence_name", "start", "stop")] #26
ZFHX3_overHypermeth <- ZFHX3_overHypermeth[ , c("sequence_name", "start", "stop")] #17

ETV1_overBackground <- AFassocTFsoverbackground[AFassocTFsoverbackground$motif_alt_id == 'ETV1',]
TBX3_overBackground <- AFassocTFsoverbackground[AFassocTFsoverbackground$motif_alt_id == 'TBX3',]
TBX5_overBackground <- AFassocTFsoverbackground[AFassocTFsoverbackground$motif_alt_id == 'TBX5',]
NKX25_overBackground <- AFassocTFsoverbackground[AFassocTFsoverbackground$motif_alt_id == 'NKX2-5',]
CREB2_overBackground <- AFassocTFsoverbackground[AFassocTFsoverbackground$motif_alt_id == 'ATF4',]
PITX2_overBackground <- AFassocTFsoverbackground[AFassocTFsoverbackground$motif_alt_id == 'PITX2',]
PRRX1_overBackground <- AFassocTFsoverbackground[AFassocTFsoverbackground$motif_alt_id == 'PRRX1',]
HAND2_overBackground <- AFassocTFsoverbackground[AFassocTFsoverbackground$motif_alt_id == 'HAND2',]
ZFHX3_overBackground <- read.table(file = 'ZFHX3_motif_overbackgrou.tsv', sep = '\t', header = TRUE)

ETV1_overBackground <- ETV1_overBackground[ , c("sequence_name", "start", "stop")] #769
TBX3_overBackground <- TBX3_overBackground[ , c("sequence_name", "start", "stop")] #1165
TBX5_overBackground <- TBX5_overBackground[ , c("sequence_name", "start", "stop")] #1338
NKX25_overBackground <- NKX25_overBackground[ , c("sequence_name", "start", "stop")] #765
CREB2_overBackground <- CREB2_overBackground[ , c("sequence_name", "start", "stop")] #582
PITX2_overBackground <- PITX2_overBackground[ , c("sequence_name", "start", "stop")] #427
PRRX1_overBackground <- PRRX1_overBackground[ , c("sequence_name", "start", "stop")] #320
HAND2_overBackground <- HAND2_overBackground[ , c("sequence_name", "start", "stop")] #636
ZFHX3_overBackground <- ZFHX3_overBackground[ , c("sequence_name", "start", "stop")] #477

dmrs_dmrff_betanorm_Chris_DMRs_hypometh_gr <- import("dmrs_dmrff_betanorm_Chris_DMRs_hypometh.bed", format = "BED")
dmrs_dmrff_betanorm_Chris_DMRs_hypermeth_gr <- import("dmrs_dmrff_betanorm_Chris_DMRs_hypermeth.bed", format = "BED")
DMRs_Chris_allcandidateregionsbackground_gr <- import("DMRs_Chris_allcandidateregionsbackground.bed", format = "BED")




hits <- findOverlaps(dmrs_dmrff_betanorm_Chris_DMRs_hypometh_gr, ETV1_overHypometh_gr)
length(unique(queryHits(hits))) #25
hits <- findOverlaps(dmrs_dmrff_betanorm_Chris_DMRs_hypometh_gr, TBX3_overHypometh_gr)
length(unique(queryHits(hits))) #21
hits <- findOverlaps(dmrs_dmrff_betanorm_Chris_DMRs_hypometh_gr, TBX5_overHypometh_gr)
length(unique(queryHits(hits))) #23
hits <- findOverlaps(dmrs_dmrff_betanorm_Chris_DMRs_hypometh_gr, NKX25_overHypometh_gr)
length(unique(queryHits(hits))) #17
hits <- findOverlaps(dmrs_dmrff_betanorm_Chris_DMRs_hypometh_gr, CREB2_overHypometh_gr)
length(unique(queryHits(hits))) #9
hits <- findOverlaps(dmrs_dmrff_betanorm_Chris_DMRs_hypometh_gr, HAND2_overHypometh_gr)
length(unique(queryHits(hits))) #18
hits <- findOverlaps(dmrs_dmrff_betanorm_Chris_DMRs_hypometh_gr, ZFHX3_overHypometh_gr)
length(unique(queryHits(hits))) #6
hits <- findOverlaps(dmrs_dmrff_betanorm_Chris_DMRs_hypometh_gr, PRRX1_overHypometh_gr)
length(unique(queryHits(hits))) #12
hits <- findOverlaps(dmrs_dmrff_betanorm_Chris_DMRs_hypometh_gr, PITX2_overHypometh_gr)
length(unique(queryHits(hits))) #11
hits <- findOverlaps(dmrs_dmrff_betanorm_Chris_DMRs_hypometh_gr, HAND2_overHypometh_gr)
length(unique(queryHits(hits))) #18


hits <- findOverlaps(dmrs_dmrff_betanorm_Chris_DMRs_hypermeth_gr, ETV1_overHypermeth_gr)
length(unique(queryHits(hits))) #23
hits <- findOverlaps(dmrs_dmrff_betanorm_Chris_DMRs_hypermeth_gr, TBX3_overHypermeth_gr)
length(unique(queryHits(hits))) #37
hits <- findOverlaps(dmrs_dmrff_betanorm_Chris_DMRs_hypermeth_gr, TBX5_overHypermeth_gr)
length(unique(queryHits(hits))) #45
hits <- findOverlaps(dmrs_dmrff_betanorm_Chris_DMRs_hypermeth_gr, NKX25_overHypermeth_gr)
length(unique(queryHits(hits))) #22
hits <- findOverlaps(dmrs_dmrff_betanorm_Chris_DMRs_hypermeth_gr, CREB2_overHypermeth_gr)
length(unique(queryHits(hits))) #15
hits <- findOverlaps(dmrs_dmrff_betanorm_Chris_DMRs_hypermeth_gr, HAND2_overHypermeth_gr)
length(unique(queryHits(hits))) #19
hits <- findOverlaps(dmrs_dmrff_betanorm_Chris_DMRs_hypermeth_gr, ZFHX3_overHypermeth_gr)
length(unique(queryHits(hits))) #13
hits <- findOverlaps(dmrs_dmrff_betanorm_Chris_DMRs_hypermeth_gr, PRRX1_overHypermeth_gr)
length(unique(queryHits(hits))) #4
hits <- findOverlaps(dmrs_dmrff_betanorm_Chris_DMRs_hypermeth_gr, PITX2_overHypermeth_gr)
length(unique(queryHits(hits))) #4
hits <- findOverlaps(dmrs_dmrff_betanorm_Chris_DMRs_hypermeth_gr, HAND2_overHypermeth_gr)
length(unique(queryHits(hits))) #19

hits <- findOverlaps(DMRs_Chris_allcandidateregionsbackground_gr, ETV1_overBackground_gr)
length(unique(queryHits(hits))) #562
hits <- findOverlaps(DMRs_Chris_allcandidateregionsbackground_gr, TBX3_overBackground_gr)
length(unique(queryHits(hits))) #729
hits <- findOverlaps(DMRs_Chris_allcandidateregionsbackground_gr, TBX5_overBackground_gr)
length(unique(queryHits(hits))) #790
hits <- findOverlaps(DMRs_Chris_allcandidateregionsbackground_gr, NKX25_overBackground_gr)
length(unique(queryHits(hits))) #566
hits <- findOverlaps(DMRs_Chris_allcandidateregionsbackground_gr, CREB2_overBackground_gr)
length(unique(queryHits(hits))) #471
hits <- findOverlaps(DMRs_Chris_allcandidateregionsbackground_gr, HAND2_overBackground_gr)
length(unique(queryHits(hits))) #482
hits <- findOverlaps(DMRs_Chris_allcandidateregionsbackground_gr, ZFHX3_overBackground_gr)
length(unique(queryHits(hits))) #323
hits <- findOverlaps(DMRs_Chris_allcandidateregionsbackground_gr, PRRX1_overBackground_gr)
length(unique(queryHits(hits))) #200
hits <- findOverlaps(DMRs_Chris_allcandidateregionsbackground_gr, PITX2_overBackground_gr)
length(unique(queryHits(hits))) #339
hits <- findOverlaps(DMRs_Chris_allcandidateregionsbackground_gr, HAND2_overBackground_gr)
length(unique(queryHits(hits))) #482


## Statistical test of FIMO occurrences of AF-associated TFs (after ZOOPS method)

####
## Over hypermeth
####



## ZFHX3
motif <- "ZFHX3"
# Input values
N <- 50610   # Total number of regions
K <- 323 # Regions with motif (total)
n <- 91     # Foreground regions
k <- 13     # Foreground regions with motif
# Hypergeometric test: P(X >= k)
# lower.tail = FALSE → P(X >= k)
p_value <- phyper(k - 1, K, N - K, n, lower.tail = FALSE)
# Print result
cat("P-value for motif enrichment:", p_value, "\n")
# Optional: fold enrichment
foreground_fraction <- k / n
background_fraction <- (K - k) / (N - n)
fold_enrichment <- foreground_fraction / background_fraction
cat("Fold enrichment:", round(fold_enrichment, 2), "\n")
new_row1 <- data.frame(motif = motif, Regions_with_motif = k, perc_regions = (k/n)*100, pval = p_value)
df <- bind_rows(new_row1)


## PRRX1
# Input values
N <- 50610   # Total number of regions
K <- 200 # Regions with motif (total)
n <- 91     # Foreground regions
k <- 4     # Foreground regions with motif
# Hypergeometric test: P(X >= k)
# lower.tail = FALSE → P(X >= k)
p_value <- phyper(k - 1, K, N - K, n, lower.tail = FALSE)
# Print result
cat("P-value for motif enrichment:", p_value, "\n")
# Optional: fold enrichment
foreground_fraction <- k / n
background_fraction <- (K - k) / (N - n)
fold_enrichment <- foreground_fraction / background_fraction
cat("Fold enrichment:", round(fold_enrichment, 2), "\n")
new_row2 <- data.frame(motif, Regions_with_motif = k, perc_regions = (k/n)*100, pval = p_value)
df <- bind_rows(new_row1, new_row2)


## NKX2-5
motif <- "NKX2-5"
# Input values
N <- 50610   # Total number of regions
K <- 566 # Regions with motif (total)
n <- 91     # Foreground regions
k <- 22     # Foreground regions with motif
# Hypergeometric test: P(X >= k)
# lower.tail = FALSE → P(X >= k)
p_value <- phyper(k - 1, K, N - K, n, lower.tail = FALSE)
# Print result
cat("P-value for motif enrichment:", p_value, "\n")
# Optional: fold enrichment
foreground_fraction <- k / n
background_fraction <- (K - k) / (N - n)
fold_enrichment <- foreground_fraction / background_fraction
cat("Fold enrichment:", round(fold_enrichment, 2), "\n")
new_row3 <- data.frame(motif, Regions_with_motif = k, perc_regions = (k/n)*100, pval = p_value)
df <- bind_rows(new_row1, new_row2, new_row3)


## PITX2
motif <- "PITX2"
# Input values
N <- 50610   # Total number of regions
K <- 339 # Regions with motif (total)
n <- 91     # Foreground regions
k <- 16     # Foreground regions with motif
# Hypergeometric test: P(X >= k)
# lower.tail = FALSE → P(X >= k)
p_value <- phyper(k - 1, K, N - K, n, lower.tail = FALSE)
# Print result
cat("P-value for motif enrichment:", p_value, "\n")
# Optional: fold enrichment
foreground_fraction <- k / n
background_fraction <- (K - k) / (N - n)
fold_enrichment <- foreground_fraction / background_fraction
cat("Fold enrichment:", round(fold_enrichment, 2), "\n")
new_row4 <- data.frame(motif, Regions_with_motif = k, perc_regions = (k/n)*100, pval = p_value)
df <- bind_rows(new_row1, new_row2, new_row3, new_row4)


## TBX5
motif <- "TBX5"
# Input values
N <- 50610   # Total number of regions
K <- 790 # Regions with motif (total)
n <- 91     # Foreground regions
k <- 45     # Foreground regions with motif
# Hypergeometric test: P(X >= k)
# lower.tail = FALSE → P(X >= k)
p_value <- phyper(k - 1, K, N - K, n, lower.tail = FALSE)
# Print result
cat("P-value for motif enrichment:", p_value, "\n")
# Optional: fold enrichment
foreground_fraction <- k / n
background_fraction <- (K - k) / (N - n)
fold_enrichment <- foreground_fraction / background_fraction
cat("Fold enrichment:", round(fold_enrichment, 2), "\n")
new_row5 <- data.frame(motif, Regions_with_motif = k, perc_regions = (k/n)*100, pval = p_value)
df <- bind_rows(new_row1, new_row2, new_row3, new_row4, new_row5)


## TBX3
motif <- "TBX3"
# Input values
N <- 50610   # Total number of regions
K <- 729 # Regions with motif (total)
n <- 91     # Foreground regions
k <- 37     # Foreground regions with motif
# Hypergeometric test: P(X >= k)
# lower.tail = FALSE → P(X >= k)
p_value <- phyper(k - 1, K, N - K, n, lower.tail = FALSE)
# Print result
cat("P-value for motif enrichment:", p_value, "\n")
# Optional: fold enrichment
foreground_fraction <- k / n
background_fraction <- (K - k) / (N - n)
fold_enrichment <- foreground_fraction / background_fraction
cat("Fold enrichment:", round(fold_enrichment, 2), "\n")
new_row6 <- data.frame(motif, Regions_with_motif = k, perc_regions = (k/n)*100, pval = p_value)
df <- bind_rows(new_row1, new_row2, new_row3, new_row4, new_row5, new_row6)


## HAND2
motif <- "HAND2"
# Input values
N <- 50610   # Total number of regions
K <- 482 # Regions with motif (total)
n <- 91     # Foreground regions
k <- 19     # Foreground regions with motif
# Hypergeometric test: P(X >= k)
# lower.tail = FALSE → P(X >= k)
p_value <- phyper(k - 1, K, N - K, n, lower.tail = FALSE)
# Print result
cat("P-value for motif enrichment:", p_value, "\n")
# Optional: fold enrichment
foreground_fraction <- k / n
background_fraction <- (K - k) / (N - n)
fold_enrichment <- foreground_fraction / background_fraction
cat("Fold enrichment:", round(fold_enrichment, 2), "\n")
new_row7 <- data.frame(motif, Regions_with_motif = k, perc_regions = (k/n)*100, pval = p_value)
df <- bind_rows(new_row1, new_row2, new_row3, new_row4, new_row5, new_row6, new_row7)


## ETV1
motif <- "ETV1"
# Input values
N <- 50610   # Total number of regions
K <- 562 # Regions with motif (total)
n <- 91     # Foreground regions
k <- 23     # Foreground regions with motif
# Hypergeometric test: P(X >= k)
# lower.tail = FALSE → P(X >= k)
p_value <- phyper(k - 1, K, N - K, n, lower.tail = FALSE)
# Print result
cat("P-value for motif enrichment:", p_value, "\n")
# Optional: fold enrichment
foreground_fraction <- k / n
background_fraction <- (K - k) / (N - n)
fold_enrichment <- foreground_fraction / background_fraction
cat("Fold enrichment:", round(fold_enrichment, 2), "\n")
new_row8 <- data.frame(motif, Regions_with_motif = k, perc_regions = (k/n)*100, pval = p_value)
df <- bind_rows(new_row1, new_row2, new_row3, new_row4, new_row5, new_row6, new_row7, new_row8)


## CREB2
motif <- "CREB2"
# Input values
N <- 50610   # Total number of regions
K <- 471 # Regions with motif (total)
n <- 91     # Foreground regions
k <- 15     # Foreground regions with motif
# Hypergeometric test: P(X >= k)
# lower.tail = FALSE → P(X >= k)
p_value <- phyper(k - 1, K, N - K, n, lower.tail = FALSE)
# Print result
cat("P-value for motif enrichment:", p_value, "\n")
# Optional: fold enrichment
foreground_fraction <- k / n
background_fraction <- (K - k) / (N - n)
fold_enrichment <- foreground_fraction / background_fraction
cat("Fold enrichment:", round(fold_enrichment, 2), "\n")
new_row9 <- data.frame(motif, Regions_with_motif = k, perc_regions = (k/n)*100, pval = p_value)
df <- bind_rows(new_row1, new_row2, new_row3, new_row4, new_row5, new_row6, new_row7, new_row8, new_row9)



####
## Over hypometh
####


## ZFHX3
motif <- "ZFHX3"
# Input values
N <- 50610   # Total number of regions
K <- 323 # Regions with motif (total)
n <- 57     # Foreground regions
k <- 6     # Foreground regions with motif
# Hypergeometric test: P(X >= k)
# lower.tail = FALSE → P(X >= k)
p_value <- phyper(k - 1, K, N - K, n, lower.tail = FALSE)
# Print result
cat("P-value for motif enrichment:", p_value, "\n")
# Optional: fold enrichment
foreground_fraction <- k / n
background_fraction <- (K - k) / (N - n)
fold_enrichment <- foreground_fraction / background_fraction
cat("Fold enrichment:", round(fold_enrichment, 2), "\n")
new_row1 <- data.frame(motif = motif, Regions_with_motif = k, perc_regions = (k/n)*100, pval = p_value)
df2 <- bind_rows(new_row1)


## PRRX1
motif <- "PRRX1"
# Input values
N <- 50610   # Total number of regions
K <- 200 # Regions with motif (total)
n <- 57     # Foreground regions
k <- 12     # Foreground regions with motif
# Hypergeometric test: P(X >= k)
# lower.tail = FALSE → P(X >= k)
p_value <- phyper(k - 1, K, N - K, n, lower.tail = FALSE)
# Print result
cat("P-value for motif enrichment:", p_value, "\n")
# Optional: fold enrichment
foreground_fraction <- k / n
background_fraction <- (K - k) / (N - n)
fold_enrichment <- foreground_fraction / background_fraction
cat("Fold enrichment:", round(fold_enrichment, 2), "\n")
new_row2 <- data.frame(motif, Regions_with_motif = k, perc_regions = (k/n)*100, pval = p_value)
df2 <- bind_rows(new_row1, new_row2)


## NKX2-5
motif <- "NKX2-5"
# Input values
N <- 50610   # Total number of regions
K <- 566 # Regions with motif (total)
n <- 57     # Foreground regions
k <- 17     # Foreground regions with motif
# Hypergeometric test: P(X >= k)
# lower.tail = FALSE → P(X >= k)
p_value <- phyper(k - 1, K, N - K, n, lower.tail = FALSE)
# Print result
cat("P-value for motif enrichment:", p_value, "\n")
# Optional: fold enrichment
foreground_fraction <- k / n
background_fraction <- (K - k) / (N - n)
fold_enrichment <- foreground_fraction / background_fraction
cat("Fold enrichment:", round(fold_enrichment, 2), "\n")
new_row3 <- data.frame(motif, Regions_with_motif = k, perc_regions = (k/n)*100, pval = p_value)
df2 <- bind_rows(new_row1, new_row2, new_row3)


## PITX2
motif <- "PITX2"
# Input values
N <- 50610   # Total number of regions
K <- 339 # Regions with motif (total)
n <- 57     # Foreground regions
k <- 11     # Foreground regions with motif
# Hypergeometric test: P(X >= k)
# lower.tail = FALSE → P(X >= k)
p_value <- phyper(k - 1, K, N - K, n, lower.tail = FALSE)
# Print result
cat("P-value for motif enrichment:", p_value, "\n")
# Optional: fold enrichment
foreground_fraction <- k / n
background_fraction <- (K - k) / (N - n)
fold_enrichment <- foreground_fraction / background_fraction
cat("Fold enrichment:", round(fold_enrichment, 2), "\n")
new_row4 <- data.frame(motif, Regions_with_motif = k, perc_regions = (k/n)*100, pval = p_value)
df2 <- bind_rows(new_row1, new_row2, new_row3, new_row4)


## TBX5
motif <- "TBX5"
# Input values
N <- 50610   # Total number of regions
K <- 790 # Regions with motif (total)
n <- 57     # Foreground regions
k <- 23     # Foreground regions with motif
# Hypergeometric test: P(X >= k)
# lower.tail = FALSE → P(X >= k)
p_value <- phyper(k - 1, K, N - K, n, lower.tail = FALSE)
# Print result
cat("P-value for motif enrichment:", p_value, "\n")
# Optional: fold enrichment
foreground_fraction <- k / n
background_fraction <- (K - k) / (N - n)
fold_enrichment <- foreground_fraction / background_fraction
cat("Fold enrichment:", round(fold_enrichment, 2), "\n")
new_row5 <- data.frame(motif, Regions_with_motif = k, perc_regions = (k/n)*100, pval = p_value)
df2 <- bind_rows(new_row1, new_row2, new_row3, new_row4, new_row5)


## TBX3
motif <- "TBX3"
# Input values
N <- 50610   # Total number of regions
K <- 729 # Regions with motif (total)
n <- 57     # Foreground regions
k <- 21     # Foreground regions with motif
# Hypergeometric test: P(X >= k)
# lower.tail = FALSE → P(X >= k)
p_value <- phyper(k - 1, K, N - K, n, lower.tail = FALSE)
# Print result
cat("P-value for motif enrichment:", p_value, "\n")
# Optional: fold enrichment
foreground_fraction <- k / n
background_fraction <- (K - k) / (N - n)
fold_enrichment <- foreground_fraction / background_fraction
cat("Fold enrichment:", round(fold_enrichment, 2), "\n")
new_row6 <- data.frame(motif, Regions_with_motif = k, perc_regions = (k/n)*100, pval = p_value)
df2 <- bind_rows(new_row1, new_row2, new_row3, new_row4, new_row5, new_row6)


## HAND2
motif <- "HAND2"
# Input values
N <- 50610   # Total number of regions
K <- 482 # Regions with motif (total)
n <- 57     # Foreground regions
k <- 18     # Foreground regions with motif
# Hypergeometric test: P(X >= k)
# lower.tail = FALSE → P(X >= k)
p_value <- phyper(k - 1, K, N - K, n, lower.tail = FALSE)
# Print result
cat("P-value for motif enrichment:", p_value, "\n")
# Optional: fold enrichment
foreground_fraction <- k / n
background_fraction <- (K - k) / (N - n)
fold_enrichment <- foreground_fraction / background_fraction
cat("Fold enrichment:", round(fold_enrichment, 2), "\n")
new_row7 <- data.frame(motif, Regions_with_motif = k, perc_regions = (k/n)*100, pval = p_value)
df2 <- bind_rows(new_row1, new_row2, new_row3, new_row4, new_row5, new_row6, new_row7)


## ETV1
motif <- "ETV1"
# Input values
N <- 50610   # Total number of regions
K <- 562 # Regions with motif (total)
n <- 57     # Foreground regions
k <- 25     # Foreground regions with motif
# Hypergeometric test: P(X >= k)
# lower.tail = FALSE → P(X >= k)
p_value <- phyper(k - 1, K, N - K, n, lower.tail = FALSE)
# Print result
cat("P-value for motif enrichment:", p_value, "\n")
# Optional: fold enrichment
foreground_fraction <- k / n
background_fraction <- (K - k) / (N - n)
fold_enrichment <- foreground_fraction / background_fraction
cat("Fold enrichment:", round(fold_enrichment, 2), "\n")
new_row8 <- data.frame(motif, Regions_with_motif = k, perc_regions = (k/n)*100, pval = p_value)
df2 <- bind_rows(new_row1, new_row2, new_row3, new_row4, new_row5, new_row6, new_row7, new_row8)


## CREB2
motif <- "CREB2"
# Input values
N <- 50610   # Total number of regions
K <- 471 # Regions with motif (total)
n <- 57     # Foreground regions
k <- 9     # Foreground regions with motif
# Hypergeometric test: P(X >= k)
# lower.tail = FALSE → P(X >= k)
p_value <- phyper(k - 1, K, N - K, n, lower.tail = FALSE)
# Print result
cat("P-value for motif enrichment:", p_value, "\n")
# Optional: fold enrichment
foreground_fraction <- k / n
background_fraction <- (K - k) / (N - n)
fold_enrichment <- foreground_fraction / background_fraction
cat("Fold enrichment:", round(fold_enrichment, 2), "\n")
new_row9 <- data.frame(motif, Regions_with_motif = k, perc_regions = (k/n)*100, pval = p_value)
df2 <- bind_rows(new_row1, new_row2, new_row3, new_row4, new_row5, new_row6, new_row7, new_row8, new_row9)


df$status <- "hypermeth"
df2$status <- "hypometh"

df_combined <- rbind(df, df2)
df_binom_combined <- rbind(df_binom, df_binom2)

padj <- p.adjust(df_combined$pval, method="bonferroni", n=dim(df_combined)[1]*2) 

df_hypermeth <- subset(df_combined, status == "hypermeth")
df_hypometh <- subset(df_combined, status == "hypometh")

combined <- cbind(df_hypermeth,df_hypometh)
color_limits <- range(combined$logpadj)
size_limits  <- range(combined$perc_regions)



# 📊 Figure 7c - Hypermethylated DMRs
ggplot(df_hypermeth, aes(y = motif, x = 1, size = perc_regions, fill = logpadj)) +
  geom_point(shape = 21, color = "black") +
  scale_size(range = c(5, 20), name = "% Regions", limits = size_limits) +
  scale_fill_gradient(low = "lightblue", high = "darkblue", name = "-log P-value",  limits = color_limits) +
  theme_minimal() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_blank()
  ) +
  ggtitle("AF-associated TFs enrichment over Hypermeth. DMRs")




# 📊 Figure 7c - hypomethylated DMRs
ggplot(df_hypometh, aes(y = motif, x = 1, size = perc_regions, fill = logpadj)) +
  geom_point(shape = 21, color = "black") +
  scale_size(range = c(5, 20), name = "% Region", limits = size_limits) +
  scale_fill_gradient(low = "lightblue", high = "darkblue", name = "-log P-value", limits = color_limits) +
  theme_minimal() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_blank()
  ) +
  ggtitle("AF-associated TFs enrichment over Hypometh. DMRs")














