


library(reshape2)
library(ggplot2)
library(dplyr)


##########################################
##### Genrich peaksets exploration #######
##########################################
##### Definition of reg. regions #########
##########################################

## Genrich tool (https://github.com/jsh58/Genrich) was used to call peaks for multiple replicates in each sample group collectively.
## https://github.com/jsh58/Genrich


setwd("Rodriguez_2025_AFepigenome/Genrich_RegRegions_definition/")


AF_LA_Genrich_greylistExcluded_g10 <- as.data.frame(read.table("AF_LA_Genrich_greylistExcluded_g10.bed",header = FALSE, sep="\t",stringsAsFactors=FALSE))
SR_LA_Genrich_greylistExcluded_g10 <- as.data.frame(read.table("SR_LA_Genrich_greylistExcluded_g10.bed",header = FALSE, sep="\t",stringsAsFactors=FALSE))
SR_RA_Genrich_greylistExcluded_g10 <- as.data.frame(read.table("SR_RA_Genrich_greylistExcluded_g10.bed",header = FALSE, sep="\t",stringsAsFactors=FALSE))
AF_RA_Genrich_greylistExcluded_g10 <- as.data.frame(read.table("AF_RA_Genrich_greylistExcluded_g10.bed",header = FALSE, sep="\t",stringsAsFactors=FALSE))

AF_RA_Genrich_greylistExcluded_g10_length <- AF_RA_Genrich_greylistExcluded_g10$V3 - AF_RA_Genrich_greylistExcluded_g10$V2
AF_LA_Genrich_greylistExcluded_g10_length <- AF_LA_Genrich_greylistExcluded_g10$V3 - AF_LA_Genrich_greylistExcluded_g10$V2
SR_RA_Genrich_greylistExcluded_g10_length <- SR_RA_Genrich_greylistExcluded_g10$V3 - SR_RA_Genrich_greylistExcluded_g10$V2
SR_LA_Genrich_greylistExcluded_g10_length <- SR_LA_Genrich_greylistExcluded_g10$V3 - SR_LA_Genrich_greylistExcluded_g10$V2


AF_RA_Genrich_greylistExcluded_g10_length <- melt(AF_RA_Genrich_greylistExcluded_g10_length)
AF_RA_Genrich_greylistExcluded_g10_length$class <- "AF_RA"

AF_LA_Genrich_greylistExcluded_g10_length <- melt(AF_LA_Genrich_greylistExcluded_g10_length)
AF_LA_Genrich_greylistExcluded_g10_length$class <- "AF_LA"

SR_RA_Genrich_greylistExcluded_g10_length <- melt(SR_RA_Genrich_greylistExcluded_g10_length)
SR_RA_Genrich_greylistExcluded_g10_length$class <- "SR_RA"

SR_LA_Genrich_greylistExcluded_g10_length <- melt(SR_LA_Genrich_greylistExcluded_g10_length)
SR_LA_Genrich_greylistExcluded_g10_length$class <- "SR_LA"


peaksets_lengths <- rbind(AF_RA_Genrich_greylistExcluded_g10_length,AF_LA_Genrich_greylistExcluded_g10_length,SR_RA_Genrich_greylistExcluded_g10_length,SR_LA_Genrich_greylistExcluded_g10_length)

# 📊 Figure S1e
ggplot(peaksets_lengths, aes(x=class, y=value)) + geom_boxplot(aes(fill=class),outlier.shape=NA) + scale_y_continuous("Peak length (bp)",limits = c(0, 1500)) + scale_x_discrete("Subgroup") + theme_classic() + theme(axis.text.x = element_text(angle=45, hjust = 1)) + scale_fill_manual(values = c("#D85B74","#E3D081","#7298bc","#91C7B1"))                                                                                                                                                                                                                                                                                                                                                                                                                                                                            
ggplot(peaksets_lengths, aes(x=class, y=value)) + geom_boxplot(aes(fill=class)) + scale_y_continuous("Peak length (bp)") + scale_x_discrete("Subgroup") + theme_classic() + theme(axis.text.x = element_text(angle=45, hjust = 1)) + scale_fill_manual(values = c("#D85B74","#E3D081","#7298bc","#91C7B1"))                                                                                                                                                                                                                                                                                                                                                                                                                                                                    

# Calulate average peakset_length
df_avg <- peaksets_lengths %>%
  group_by(class) %>%
  summarise(avg_value = mean(value, na.rm = TRUE))
print(df_avg)


Number_of_peaks<-read.csv("Number_of_peaks_Genrich.csv")

# 📊 Figure S1d
ggplot(Number_of_peaks, aes(x=(factor(Group)), y=Enriched.peaks, fill=X)) + 
  geom_bar(stat="identity", position="stack", color="black") +
  scale_fill_manual(values = c("grey","grey","grey","grey","#D85B74","#E3D081","#7298bc","#91c7b1")) + theme_classic() + scale_y_continuous("Sample group peaks")


## Regulatory regions distance to TSSs

# 📊 Figure S1f
# The TSS distance plots were produced using the Python script control_plots.py.
# To recreate the Python enivronment: 
# module load anaconda3
# conda env create -f bedtools_py.yaml
# conda activate plots_env
# python control_plots.py \
#   -i AF_LA_Genrich_greylistExcluded_g10.bed AF_RA_Genrich_greylistExcluded_g10.bed SR_LA_Genrich_greylistExcluded_g10.bed SR_RA_Genrich_greylistExcluded_g10.bed \
#   -l AFLA-Atria AFRA-Atria SRLA-Atria SRRA-Atria \
#   -t "AF dataset - H3K27ac" \
#   -o out-barplot.svg out-boxplot.svg out-tss.svg  \
#   -tss TSS.biomart.Homo_sapiens.bed \
#   --groupby "Tissues=Atria"

# The above will output figure in .svg format







