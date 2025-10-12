

library(ggsignif)
library(ggplot2)
library(dplyr)

###############################################
###############################################
############ QPCR data plots ##################
###############################################
###############################################


RTqPCR_220824 <- read.table("RTqPCR_220824_dataanalysis.txt", sep="\t",header=TRUE)
RTqPCR_220824_BMP10 <- RTqPCR_220824[which (RTqPCR_220824$Gene == "BMP10"),]
RTqPCR_220824_PITX2c <- RTqPCR_220824[which (RTqPCR_220824$Gene == "PITX2"),]
RTqPCR_220824_MYH7 <- RTqPCR_220824[which (RTqPCR_220824$Gene == "MYH7"),]
RTqPCR_220824_MYH6 <- RTqPCR_220824[which (RTqPCR_220824$Gene == "MYH6"),]
RTqPCR_220824_NPPA <- RTqPCR_220824[which (RTqPCR_220824$Gene == "NPPA"),]
RTqPCR_220824_NPPB <- RTqPCR_220824[which (RTqPCR_220824$Gene == "NPPB"),]
RTqPCR_220824_VASH1 <- RTqPCR_220824[which (RTqPCR_220824$Gene == "VASH1"),]
RTqPCR_220824_HNF4A <- RTqPCR_220824[which (RTqPCR_220824$Gene == "HNF4A"),]
RTqPCR_220824_MT1X <- RTqPCR_220824[which (RTqPCR_220824$Gene == "MT1X"),]
RTqPCR_220824_ANGPTL2 <- RTqPCR_220824[which (RTqPCR_220824$Gene == "ANGPTL2"),]
RTqPCR_220824_STRN <- RTqPCR_220824[which (RTqPCR_220824$Gene == "STRN"),]
RTqPCR_220824_LRRC32 <- RTqPCR_220824[which (RTqPCR_220824$Gene == "LRRC32"),]
RTqPCR_220824_RGS6 <- RTqPCR_220824[which (RTqPCR_220824$Gene == "RGS6"),]


data_analysis_dataframe_withdata <- read.csv("data_analysis_dataframe_RTqpcr291124.csv")
data_analysis_dataframe_withdata <- data_analysis_dataframe_withdata[!is.na(data_analysis_dataframe_withdata$Avg_ct), ]
data_analysis_dataframe_withdata$Anat_side <- NA
data_analysis_dataframe_withdata$Disease_status <- NA
data_analysis_dataframe_withdata$Individual_ID <- NA
data_analysis_dataframe_withdata$Disease_AnatSide <- NA
data_analysis_dataframe_withdata$Anat_side[data_analysis_dataframe_withdata$Sample %in% c("Sample1","Sample3","Sample5", "Sample7", "Sample9")] <- "RA"
data_analysis_dataframe_withdata$Anat_side[data_analysis_dataframe_withdata$Sample %in% c("Sample2","Sample4","Sample6", "Sample8", "Sample10")] <- "LA"
data_analysis_dataframe_withdata$Disease_status[data_analysis_dataframe_withdata$Sample %in% c("Sample1","Sample2","Sample7", "Sample8", "Sample9","Sample10")] <- "AF"
data_analysis_dataframe_withdata$Disease_status[data_analysis_dataframe_withdata$Sample %in% c("Sample3","Sample4","Sample5", "Sample6")] <- "SR"
data_analysis_dataframe_withdata$Individual_ID[data_analysis_dataframe_withdata$Sample %in% c("Sample1")] <- "BCVR-9231"
data_analysis_dataframe_withdata$Individual_ID[data_analysis_dataframe_withdata$Sample %in% c("Sample2")] <- "BCVR-9231"
data_analysis_dataframe_withdata$Individual_ID[data_analysis_dataframe_withdata$Sample %in% c("Sample3")] <- "BCVR-5890"
data_analysis_dataframe_withdata$Individual_ID[data_analysis_dataframe_withdata$Sample %in% c("Sample4")] <- "BCVR-5890"
data_analysis_dataframe_withdata$Individual_ID[data_analysis_dataframe_withdata$Sample %in% c("Sample5")] <- "BCVR-5888"
data_analysis_dataframe_withdata$Individual_ID[data_analysis_dataframe_withdata$Sample %in% c("Sample6")] <- "BCVR-5888"
data_analysis_dataframe_withdata$Individual_ID[data_analysis_dataframe_withdata$Sample %in% c("Sample7")] <- "BCVR-3913"
data_analysis_dataframe_withdata$Individual_ID[data_analysis_dataframe_withdata$Sample %in% c("Sample8")] <- "BCVR-3913"
data_analysis_dataframe_withdata$Individual_ID[data_analysis_dataframe_withdata$Sample %in% c("Sample9")] <- "BCVR-1638"
data_analysis_dataframe_withdata$Individual_ID[data_analysis_dataframe_withdata$Sample %in% c("Sample10")] <- "BCVR-1638"
data_analysis_dataframe_withdata$Disease_AnatSide[data_analysis_dataframe_withdata$Sample %in% c("Sample2","Sample8","Sample10")] <- "AF-LA"
data_analysis_dataframe_withdata$Disease_AnatSide[data_analysis_dataframe_withdata$Sample %in% c("Sample1","Sample7","Sample9")] <- "AF-RA"
data_analysis_dataframe_withdata$Disease_AnatSide[data_analysis_dataframe_withdata$Sample %in% c("Sample4","Sample6")] <- "SR-LA"
data_analysis_dataframe_withdata$Disease_AnatSide[data_analysis_dataframe_withdata$Sample %in% c("Sample3","Sample5")] <- "SR-RA"
data_analysis_dataframe_withdata_ANGPT1 <- data_analysis_dataframe_withdata[which (data_analysis_dataframe_withdata$Marker == "ANGPT1"),]
data_analysis_dataframe_withdata_BLM <- data_analysis_dataframe_withdata[which (data_analysis_dataframe_withdata$Marker == "BLM"),]
data_analysis_dataframe_withdata_SLC43A3 <- data_analysis_dataframe_withdata[which (data_analysis_dataframe_withdata$Marker == "SLC43A3"),]
data_analysis_dataframe_withdata_NR2F1 <- data_analysis_dataframe_withdata[which (data_analysis_dataframe_withdata$Marker == "NR2F1"),]
data_analysis_dataframe_withdata_FOXS1 <- data_analysis_dataframe_withdata[which (data_analysis_dataframe_withdata$Marker == "FOXS1"),]



data_analysis_dataframe_withdata <- read.csv("data_analysis_dataframe_RTqpcr270625.csv")
data_analysis_dataframe_withdata <- data_analysis_dataframe_withdata[!is.na(data_analysis_dataframe_withdata$Avg_Ct), ]
data_analysis_dataframe_withdata$Anat_side <- NA
data_analysis_dataframe_withdata$Disease_status <- NA
data_analysis_dataframe_withdata$Individual_ID <- NA
data_analysis_dataframe_withdata$Disease_AnatSide <- NA
data_analysis_dataframe_withdata$Anat_side[data_analysis_dataframe_withdata$Sample %in% c("Sample1","Sample3","Sample5", "Sample7", "Sample9")] <- "RA"
data_analysis_dataframe_withdata$Anat_side[data_analysis_dataframe_withdata$Sample %in% c("Sample2","Sample4","Sample6", "Sample8", "Sample10")] <- "LA"
data_analysis_dataframe_withdata$Disease_status[data_analysis_dataframe_withdata$Sample %in% c("Sample1","Sample2","Sample7", "Sample8", "Sample9","Sample10")] <- "AF"
data_analysis_dataframe_withdata$Disease_status[data_analysis_dataframe_withdata$Sample %in% c("Sample3","Sample4","Sample5", "Sample6")] <- "SR"
data_analysis_dataframe_withdata$Individual_ID[data_analysis_dataframe_withdata$Sample %in% c("Sample1")] <- "BCVR-9231"
data_analysis_dataframe_withdata$Individual_ID[data_analysis_dataframe_withdata$Sample %in% c("Sample2")] <- "BCVR-9231"
data_analysis_dataframe_withdata$Individual_ID[data_analysis_dataframe_withdata$Sample %in% c("Sample3")] <- "BCVR-5890"
data_analysis_dataframe_withdata$Individual_ID[data_analysis_dataframe_withdata$Sample %in% c("Sample4")] <- "BCVR-5890"
data_analysis_dataframe_withdata$Individual_ID[data_analysis_dataframe_withdata$Sample %in% c("Sample5")] <- "BCVR-5888"
data_analysis_dataframe_withdata$Individual_ID[data_analysis_dataframe_withdata$Sample %in% c("Sample6")] <- "BCVR-5888"
data_analysis_dataframe_withdata$Individual_ID[data_analysis_dataframe_withdata$Sample %in% c("Sample7")] <- "BCVR-3913"
data_analysis_dataframe_withdata$Individual_ID[data_analysis_dataframe_withdata$Sample %in% c("Sample8")] <- "BCVR-3913"
data_analysis_dataframe_withdata$Individual_ID[data_analysis_dataframe_withdata$Sample %in% c("Sample9")] <- "BCVR-1638"
data_analysis_dataframe_withdata$Individual_ID[data_analysis_dataframe_withdata$Sample %in% c("Sample10")] <- "BCVR-1638"
data_analysis_dataframe_withdata$Disease_AnatSide[data_analysis_dataframe_withdata$Sample %in% c("Sample2","Sample8","Sample10")] <- "AF-LA"
data_analysis_dataframe_withdata$Disease_AnatSide[data_analysis_dataframe_withdata$Sample %in% c("Sample1","Sample7","Sample9")] <- "AF-RA"
data_analysis_dataframe_withdata$Disease_AnatSide[data_analysis_dataframe_withdata$Sample %in% c("Sample4","Sample6")] <- "SR-LA"
data_analysis_dataframe_withdata$Disease_AnatSide[data_analysis_dataframe_withdata$Sample %in% c("Sample3","Sample5")] <- "SR-RA"
data_analysis_dataframe_withdata_SCX <- data_analysis_dataframe_withdata[which (data_analysis_dataframe_withdata$Marker == "SCX"),]




################ RT-qPCR ######################



###### ###### ###### ###### ###### ###### ###### ###### ###### 
###### Genes detected as AF-upregulated in main cohort  ###### 
###### ###### ###### ###### ###### ###### ###### ###### ###### 


## Quality control markers

RTqPCR_220824_cardiac_markers <- rbind(RTqPCR_220824_MYH6, RTqPCR_220824_NPPA, RTqPCR_220824_MYH7,RTqPCR_220824_HNF4A)
RTqPCR_220824_cardiac_markers$Individual_ID <- factor(RTqPCR_220824_cardiac_markers$Individual_ID)

# 📊 Figure S5b

RTqPCR_220824_cardiac_markers_PLOT_TBPnorm <- ggplot(RTqPCR_220824_cardiac_markers , aes(y = Relative_expression_TBPnorm, x= factor(Gene,level=c("MYH6","NPPA","MYH7","HNF4A")), fill = Individual_ID)) + 
  geom_boxplot(colour = "black", fill = "#F2FDFD", width=0.4) +
  geom_dotplot(binaxis = "y", stackdir='center', dotsize = 1.5) +
  scale_y_continuous(name = "Relative Expression") +
  scale_x_discrete(name = "Genes") +
  scale_fill_manual(values = c("BCVR-1638"="#FFF2CC","BCVR-3913"="#E4A3A1","BCVR-5888"="#c8ee90","BCVR-5890"="#56BED3","BCVR-9231"="#ffcccb"),
                    breaks = c("BCVR-9231","BCVR-1638","BCVR-3913","BCVR-5890","BCVR-5888") )+#,
  # labels = c("BCVR-9231 (AF)","BCVR-1638 (AF)","BCVR-3913 (AF)","BCVR-5890 (SR)","BCVR-5888 (SR)") )+
  theme_classic() +
  theme(aspect.ratio = 1/3) +
  theme(axis.text.x = element_text(size = 10), axis.text.y =element_text(size = 10) )
RTqPCR_220824_cardiac_markers_PLOT_TBPnorm


RTqPCR_220824_PITX2c_PLOT_TBPnorm <- ggplot(RTqPCR_220824_PITX2c , aes(y = Relative_expression_TBPnorm, x= factor(Anatomical_side,level=c("RA","LA")), fill = Individual_ID)) + 
  geom_boxplot(colour = "black", fill = "#F2FDFD", width=0.4) +
  geom_dotplot(binaxis = "y", stackdir='center', dotsize = 1.8) +                 
  scale_y_continuous(name = "Relative Expression") +
  scale_x_discrete(name = "Anatomical side",expand = c(1, 1)) +
  scale_fill_manual(values = c("BCVR-1638"="#FFF2CC","BCVR-3913"="#E4A3A1","BCVR-5888"="#c8ee90","BCVR-5890"="#56BED3","BCVR-9231"="#ffcccb"),
                    breaks = c("BCVR-9231","BCVR-1638","BCVR-3913","BCVR-5890","BCVR-5888") )+#,
  theme_classic() +
  theme(aspect.ratio = 1/3) +
  theme(axis.text.x = element_text(size = 10), axis.text.y =element_text(size = 10) )
RTqPCR_220824_PITX2c_PLOT_TBPnorm + 
  geom_signif(
    comparisons = list(c("RA", "LA")),
    y_position = c(0.55),
    map_signif_level = FALSE,
    annotations = " "
  )

RTqPCR_220824_PITX2c_RA<- RTqPCR_220824_PITX2c[RTqPCR_220824_PITX2c$Anatomical_side == "RA", ]
RTqPCR_220824_PITX2c_LA<- RTqPCR_220824_PITX2c[RTqPCR_220824_PITX2c$Anatomical_side == "LA", ]

t.test(RTqPCR_220824_PITX2c_LA$Relative_expression_TBPnorm, RTqPCR_220824_PITX2c_RA$Relative_expression_TBPnorm, alternative = "greater", var.equal = TRUE)
#p val = 0.003793


RTqPCR_220824_BMP10_PLOT_TBPnorm <- ggplot(RTqPCR_220824_BMP10 , aes(y = Relative_expression_TBPnorm, x= factor(Anatomical_side,level=c("RA","LA")), fill = Individual_ID)) + 
  geom_boxplot(colour = "black", fill = "#F2FDFD", width=0.4) +
  geom_dotplot(binaxis = "y", stackdir='center', dotsize = 1.8) +                
  scale_y_continuous(name = "Relative Expression") +
  scale_x_discrete(name = "Anatomical side",expand = c(1, 1)) +
  scale_fill_manual(values = c("BCVR-1638"="#FFF2CC","BCVR-3913"="#E4A3A1","BCVR-5888"="#c8ee90","BCVR-5890"="#56BED3","BCVR-9231"="#ffcccb"),
                    breaks = c("BCVR-9231","BCVR-1638","BCVR-3913","BCVR-5890","BCVR-5888") )+#,
  theme_classic() +
  theme(aspect.ratio = 1/3) +
  theme(axis.text.x = element_text(size = 10), axis.text.y =element_text(size = 10) )
RTqPCR_220824_BMP10_PLOT_TBPnorm + 
  geom_signif(
    comparisons = list(c("RA", "LA")),
    y_position = c(410),
    map_signif_level = FALSE,
    annotations = " "
  )

RTqPCR_220824_BMP10_RA<- RTqPCR_220824_BMP10[RTqPCR_220824_BMP10$Anatomical_side == "RA", ]
RTqPCR_220824_BMP10_LA<- RTqPCR_220824_BMP10[RTqPCR_220824_BMP10$Anatomical_side == "LA", ]

t.test(RTqPCR_220824_BMP10_RA$Relative_expression_TBPnorm, RTqPCR_220824_BMP10_LA$Relative_expression_TBPnorm, alternative = "greater", var.equal = TRUE)
#p val = 0.003111



# 📊 Figure 4 barplots

## NPPB

RTqPCR_220824_NPPB_AF <- RTqPCR_220824_NPPB[RTqPCR_220824_NPPB$Disease_status == "AF", ]
RTqPCR_220824_NPPB_SR <- RTqPCR_220824_NPPB[RTqPCR_220824_NPPB$Disease_status == "SR", ]

t_test_result <- t.test(RTqPCR_220824_NPPB_AF$Relative_expression_TBPnorm, RTqPCR_220824_NPPB_SR$Relative_expression_TBPnorm, alternative = "greater")
t_test_result

RTqPCR_220824_NPPB_TBPnorm <- ggplot(RTqPCR_220824_NPPB, 
                                     aes(x = factor(Disease_status, levels = c("AF", "SR")), 
                                         y = Relative_expression_TBPnorm)) +
  stat_summary(aes(fill = Disease_status), fun = mean, geom = "bar", colour = "black", width = 0.4) +
  geom_dotplot(aes(fill = Individual_ID), binaxis = "y", stackdir = "center", dotsize = 1.3) +
    scale_fill_manual(
    name = "Individual_ID",  
    values = c(
      "AF" = "#FDE6E0", 
      "SR" = "#F2FDFD",
      "BCVR-1638" = "#FFF2CC", 
      "BCVR-3913" = "#E4A3A1", 
      "BCVR-5888" = "#c8ee90", 
      "BCVR-5890" = "#56BED3", 
      "BCVR-9231" = "#ffcccb"
    ),
    breaks = c("AF", "SR", "BCVR-9231", "BCVR-1638", "BCVR-3913", "BCVR-5890", "BCVR-5888")
  ) +
  scale_y_continuous(name = "Relative Expression") +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 18), 
    axis.title.y = element_text(size = 18),
    axis.text.y = element_text(size = 16)
  )
RTqPCR_220824_NPPB_TBPnorm + geom_signif(
  comparisons = list(c("AF", "SR")), 
  map_signif_level = FALSE, 
  step_increase = 0.1,     
  textsize = 4  ,
  y_position = 66, 
  annotations = "p = 0.011"
)


## VASH1

RTqPCR_220824_VASH1_AF <- RTqPCR_220824_VASH1[RTqPCR_220824_VASH1$Disease_status == "AF", ]
RTqPCR_220824_VASH1_SR <- RTqPCR_220824_VASH1[RTqPCR_220824_VASH1$Disease_status == "SR", ]

t_test_result <- t.test(RTqPCR_220824_VASH1_AF$Relative_expression_TBPnorm, RTqPCR_220824_VASH1_SR$Relative_expression_TBPnorm, alternative = "greater")
t_test_result
#pval 0.05

RTqPCR_220824_VASH1_TBPnorm <- ggplot(RTqPCR_220824_VASH1, 
                                      aes(x = factor(Disease_status, levels = c("AF", "SR")), 
                                          y = Relative_expression_TBPnorm)) +
  stat_summary(aes(fill = Disease_status), fun = mean, geom = "bar", colour = "black", width = 0.4) +
  geom_dotplot(aes(fill = Individual_ID), binaxis = "y", stackdir = "center", dotsize = 1.3) +
  
  # Custom bar colors for Disease_status
  scale_fill_manual(
    name = "Individual_ID", 
    values = c(
      "AF" = "#FDE6E0", 
      "SR" = "#F2FDFD",
      "BCVR-1638" = "#FFF2CC", 
      "BCVR-3913" = "#E4A3A1", 
      "BCVR-5888" = "#c8ee90", 
      "BCVR-5890" = "#56BED3", 
      "BCVR-9231" = "#ffcccb"
    ),
    breaks = c("AF", "SR", "BCVR-9231", "BCVR-1638", "BCVR-3913", "BCVR-5890", "BCVR-5888")
  ) +
  
  scale_y_continuous(name = "Relative Expression") +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 18), 
    axis.title.y = element_text(size = 18),
    axis.text.y = element_text(size = 16)
  )
RTqPCR_220824_VASH1_TBPnorm + geom_signif(
  comparisons = list(c("AF", "SR")), # Pairwise comparisons
  map_signif_level = FALSE, # Display significance levels (*, **, ***)
  step_increase = 0.1,     # Adjust spacing between brackets
  textsize = 4  ,
  y_position = 2.1, # Adjust text size
  annotations = "p = 0.05"
)



## ANGPTL2

RTqPCR_220824_ANGPTL2_AF <- RTqPCR_220824_VASH1[RTqPCR_220824_ANGPTL2$Disease_status == "AF", ]
RTqPCR_220824_ANGPTL2_SR <- RTqPCR_220824_VASH1[RTqPCR_220824_ANGPTL2$Disease_status == "SR", ]

t_test_result <- t.test(RTqPCR_220824_ANGPTL2_AF$Relative_expression_TBPnorm, RTqPCR_220824_ANGPTL2_SR$Relative_expression_TBPnorm, alternative = "greater")
t_test_result
#pval 0.05059

RTqPCR_220824_ANGPTL2_TBPnorm <- ggplot(RTqPCR_220824_ANGPTL2, 
                                        aes(x = factor(Disease_status, levels = c("AF", "SR")), 
                                            y = Relative_expression_TBPnorm)) +
  stat_summary(aes(fill = Disease_status), fun = mean, geom = "bar", colour = "black", width = 0.4) +
  geom_dotplot(aes(fill = Individual_ID), binaxis = "y", stackdir = "center", dotsize = 1.3) +
    scale_fill_manual(
    name = "Individual_ID",  
    values = c(
      "AF" = "#FDE6E0", 
      "SR" = "#F2FDFD",
      "BCVR-1638" = "#FFF2CC", 
      "BCVR-3913" = "#E4A3A1", 
      "BCVR-5888" = "#c8ee90", 
      "BCVR-5890" = "#56BED3", 
      "BCVR-9231" = "#ffcccb"
    ),
    breaks = c("AF", "SR", "BCVR-9231", "BCVR-1638", "BCVR-3913", "BCVR-5890", "BCVR-5888")
  ) +
  
  scale_y_continuous(name = "Relative Expression") +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 18), 
    axis.title.y = element_text(size = 18),
    axis.text.y = element_text(size = 16)
  ) + geom_signif(
    comparisons = list(c("AF", "SR")),
    map_signif_level = FALSE, 
    step_increase = 0.1,    
    textsize = 4  ,
    y_position = 4.6, 
    annotations = "p = 0.05"
  )


## STRN

RTqPCR_220824_STRN_AF <- RTqPCR_220824_STRN[RTqPCR_220824_STRN$Disease_status == "AF", ]
RTqPCR_220824_STRN_SR <- RTqPCR_220824_STRN[RTqPCR_220824_STRN$Disease_status == "SR", ]

t_test_result <- t.test(RTqPCR_220824_STRN_AF$Relative_expression_GAPDHnorm, RTqPCR_220824_STRN_SR$Relative_expression_GAPDHnorm, alternative = "greater")
t_test_result
#pval 0.008994

RTqPCR_220824_STRN_TBPnorm <- ggplot(RTqPCR_220824_STRN, 
                                     aes(x = factor(Disease_status, levels = c("AF", "SR")), 
                                         y = Relative_expression_TBPnorm)) +
  stat_summary(aes(fill = Disease_status), fun = mean, geom = "bar", colour = "black", width = 0.4) +
  geom_dotplot(aes(fill = Individual_ID), binaxis = "y", stackdir = "center", dotsize = 1.3) +
    scale_fill_manual(
    name = "Individual_ID",  
    values = c(
      "AF" = "#FDE6E0", 
      "SR" = "#F2FDFD",
      "BCVR-1638" = "#FFF2CC", 
      "BCVR-3913" = "#E4A3A1", 
      "BCVR-5888" = "#c8ee90", 
      "BCVR-5890" = "#56BED3", 
      "BCVR-9231" = "#ffcccb"
    ),
    breaks = c("AF", "SR", "BCVR-9231", "BCVR-1638", "BCVR-3913", "BCVR-5890", "BCVR-5888")
  ) +
  
  scale_y_continuous(name = "Relative Expression") +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 18), 
    axis.title.y = element_text(size = 18),
    axis.text.y = element_text(size = 16)
  ) + geom_signif(
    comparisons = list(c("AF", "SR")), 
    map_signif_level = FALSE, 
    step_increase = 0.1,     
    textsize = 4  ,
    y_position = 5.1, 
    annotations = "p = 0.006"
  )



## LRRC32

RTqPCR_220824_LRRC32_AF <- RTqPCR_220824_LRRC32[RTqPCR_220824_LRRC32$Disease_status == "AF", ]
RTqPCR_220824_LRRC32_SR <- RTqPCR_220824_LRRC32[RTqPCR_220824_LRRC32$Disease_status == "SR", ]

t_test_result <- t.test(RTqPCR_220824_LRRC32_AF$Relative_expression_TBPnorm, RTqPCR_220824_LRRC32_SR$Relative_expression_TBPnorm, alternative = "greater")
t_test_result
#pval 0.02609


RTqPCR_220824_LRRC32_TBPnorm <- ggplot(RTqPCR_220824_LRRC32, 
                                       aes(x = factor(Disease_status, levels = c("AF", "SR")), 
                                           y = Relative_expression_TBPnorm)) +
  stat_summary(aes(fill = Disease_status), fun = mean, geom = "bar", colour = "black", width = 0.4) +
  geom_dotplot(aes(fill = Individual_ID), binaxis = "y", stackdir = "center", dotsize = 1.3) +
    scale_fill_manual(
    name = "Individual_ID",  # temporary name for legend (gets overridden below)
    values = c(
      "AF" = "#FDE6E0", 
      "SR" = "#F2FDFD",
      "BCVR-1638" = "#FFF2CC", 
      "BCVR-3913" = "#E4A3A1", 
      "BCVR-5888" = "#c8ee90", 
      "BCVR-5890" = "#56BED3", 
      "BCVR-9231" = "#ffcccb"
    ),
    breaks = c("AF", "SR", "BCVR-9231", "BCVR-1638", "BCVR-3913", "BCVR-5890", "BCVR-5888")
  ) +
  
  scale_y_continuous(name = "Relative Expression") +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 18), 
    axis.title.y = element_text(size = 18),
    axis.text.y = element_text(size = 16)
  )+ geom_signif(
    comparisons = list(c("AF", "SR")), 
    map_signif_level = FALSE, 
    step_increase = 0.1,     
    textsize = 4  ,
    y_position = 1.02, 
    annotations = "p = 0.026"
  )


## SLC43A3

t_test_result <- t.test(data_analysis_dataframe_withdata_SLC43A3_AF$Relative_expression_TBPnorm, data_analysis_dataframe_withdata_SLC43A3_SR$Relative_expression_TBPnorm, alternative = "greater")
t_test_result
#pval  0.4918

RTqPCR_SLC43A3_TBPnorm <- ggplot(data_analysis_dataframe_withdata_SLC43A3, 
                                 aes(x = factor(Disease_status, levels = c("AF", "SR")), 
                                     y = Relative_expression_TBPnorm)) +
  stat_summary(aes(fill = Disease_status), fun = mean, geom = "bar", colour = "black", width = 0.4) +
  geom_dotplot(aes(fill = Individual_ID), binaxis = "y", stackdir = "center", dotsize = 1.3) +
    scale_fill_manual(
    name = "Individual_ID",  
    values = c(
      "AF" = "#FDE6E0", 
      "SR" = "#F2FDFD",
      "BCVR-1638" = "#FFF2CC", 
      "BCVR-3913" = "#E4A3A1", 
      "BCVR-5888" = "#c8ee90", 
      "BCVR-5890" = "#56BED3", 
      "BCVR-9231" = "#ffcccb"
    ),
    breaks = c("AF", "SR", "BCVR-9231", "BCVR-1638", "BCVR-3913", "BCVR-5890", "BCVR-5888")
  ) +
  
  scale_y_continuous(name = "Relative Expression") +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 18), 
    axis.title.y = element_text(size = 18),
    axis.text.y = element_text(size = 16)
  )+ geom_signif(
    comparisons = list(c("AF", "SR")), 
    map_signif_level = FALSE, 
    step_increase = 0.1,    
    textsize = 4  ,
    y_position = 0.63, 
    annotations = "p = 0.49"
  )




## FOXS1

t_test_result <- t.test(data_analysis_dataframe_withdata_FOXS1_AF$Relative_expression_TBPnorm, data_analysis_dataframe_withdata_FOXS1_SR$Relative_expression_TBPnorm, alternative = "greater")
t_test_result
#pval 0.659


RTqPCR_FOXS1_TBPnorm <- ggplot(data_analysis_dataframe_withdata_FOXS1, 
                               aes(x = factor(Disease_status, levels = c("AF", "SR")), 
                                   y = Relative_expression_TBPnorm)) +
  stat_summary(aes(fill = Disease_status), fun = mean, geom = "bar", colour = "black", width = 0.4) +
  geom_dotplot(aes(fill = Individual_ID), binaxis = "y", stackdir = "center", dotsize = 1.3) +
  
  # Custom bar colors for Disease_status
  scale_fill_manual(
    name = "Individual_ID",  # temporary name for legend (gets overridden below)
    values = c(
      "AF" = "#FDE6E0", 
      "SR" = "#F2FDFD",
      "BCVR-1638" = "#FFF2CC", 
      "BCVR-3913" = "#E4A3A1", 
      "BCVR-5888" = "#c8ee90", 
      "BCVR-5890" = "#56BED3", 
      "BCVR-9231" = "#ffcccb"
    ),
    breaks = c("AF", "SR", "BCVR-9231", "BCVR-1638", "BCVR-3913", "BCVR-5890", "BCVR-5888")
  ) +
  
  scale_y_continuous(name = "Relative Expression") +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 18), 
    axis.title.y = element_text(size = 18),
    axis.text.y = element_text(size = 16)
  ) + geom_signif(
    comparisons = list(c("AF", "SR")), # Pairwise comparisons
    map_signif_level = FALSE, # Display significance levels (*, **, ***)
    step_increase = 0.1,     # Adjust spacing between brackets
    textsize = 4  ,
    y_position = 1.21, # Adjust text size
    annotations = "p = 0.66"
  )


## SCX


data_analysis_dataframe_withdata_SCX_AF <- data_analysis_dataframe_withdata_SCX[data_analysis_dataframe_withdata_SCX$Disease_status == "AF", ]
data_analysis_dataframe_withdata_SCX_SR <- data_analysis_dataframe_withdata_SCX[data_analysis_dataframe_withdata_SCX$Disease_status == "SR", ]
t_test_result <- t.test(data_analysis_dataframe_withdata_SCX_AF$Relative_expression_TBPnorm, data_analysis_dataframe_withdata_SCX_SR$Relative_expression_TBPnorm, alternative = "greater")


RTqPCR_SCX_TBPnorm <- ggplot(data_analysis_dataframe_withdata_SCX, 
                             aes(x = factor(Disease_status, levels = c("AF", "SR")), 
                                 y = Relative_expression_TBPnorm)) +
  stat_summary(aes(fill = Disease_status), fun = mean, geom = "bar", colour = "black", width = 0.4) +
  geom_dotplot(aes(fill = Individual_ID), binaxis = "y", stackdir = "center", dotsize = 1.3) +
    scale_fill_manual(
    name = "Individual_ID", 
    values = c(
      "AF" = "#FDE6E0", 
      "SR" = "#F2FDFD",
      "BCVR-1638" = "#FFF2CC", 
      "BCVR-3913" = "#E4A3A1", 
      "BCVR-5888" = "#c8ee90", 
      "BCVR-5890" = "#56BED3", 
      "BCVR-9231" = "#ffcccb"
    ),
    breaks = c("AF", "SR", "BCVR-9231", "BCVR-1638", "BCVR-3913", "BCVR-5890", "BCVR-5888")
  ) +
  
  scale_y_continuous(name = "Relative Expression") +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 19), 
    axis.title.y = element_text(size = 18),
    axis.text.y = element_text(size = 16) 
  )
RTqPCR_SCX_TBPnorm + geom_signif(
  comparisons = list(c("AF", "SR")),
  map_signif_level = FALSE, 
  step_increase = 0.1,    
  textsize = 4  ,
  y_position = 1.22, # Adjust text size
  annotations = "p = 0.038"
  
)




###### Genes detected as SR-upregulated in main cohort


## MT1X

RTqPCR_220824_MT1X_AF <- RTqPCR_220824_MT1X[RTqPCR_220824_MT1X$Disease_status == "AF", ]
RTqPCR_220824_MT1X_SR <- RTqPCR_220824_MT1X[RTqPCR_220824_MT1X$Disease_status == "SR", ]

t_test_result <- t.test(RTqPCR_220824_MT1X_SR$Relative_expression_TBPnorm, RTqPCR_220824_MT1X_AF$Relative_expression_TBPnorm, alternative = "greater")
t_test_result
#pval  0.2


RTqPCR_220824_MT1X_TBPnorm <- ggplot(RTqPCR_220824_MT1X, 
                                     aes(x = factor(Disease_status, levels = c("AF", "SR")), 
                                         y = Relative_expression_TBPnorm)) +
  stat_summary(aes(fill = Disease_status), fun = mean, geom = "bar", colour = "black", width = 0.4) +
  geom_dotplot(aes(fill = Individual_ID), binaxis = "y", stackdir = "center", dotsize = 1.3) +
    scale_fill_manual(
    name = "Individual_ID",  
    values = c(
      "AF" = "#FDE6E0", 
      "SR" = "#F2FDFD",
      "BCVR-1638" = "#FFF2CC", 
      "BCVR-3913" = "#E4A3A1", 
      "BCVR-5888" = "#c8ee90", 
      "BCVR-5890" = "#56BED3", 
      "BCVR-9231" = "#ffcccb"
    ),
    breaks = c("AF", "SR", "BCVR-9231", "BCVR-1638", "BCVR-3913", "BCVR-5890", "BCVR-5888")
  ) +
  
  scale_y_continuous(name = "Relative Expression") +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 21), 
    axis.title.y = element_text(size = 21),
    axis.text.y = element_text(size = 16) 
  )
RTqPCR_220824_MT1X_TBPnorm + geom_signif(
  comparisons = list(c("AF", "SR")), 
  map_signif_level = FALSE, 
  step_increase = 0.1,     
  textsize = 6  ,
  y_position = 1.38, 
  annotations = "p = 0.2"
)


## RGS6

RTqPCR_220824_RGS6_AF <- RTqPCR_220824_RGS6[RTqPCR_220824_RGS6$Disease_status == "AF", ]
RTqPCR_220824_RGS6_SR <- RTqPCR_220824_RGS6[RTqPCR_220824_RGS6$Disease_status == "SR", ]

t_test_result <- t.test(RTqPCR_220824_RGS6_SR$Relative_expression_TBPnorm, RTqPCR_220824_RGS6_AF$Relative_expression_TBPnorm, alternative = "greater")
t_test_result
#pval  0.04699

RTqPCR_220824_RGS6_TBPnorm <- ggplot(RTqPCR_220824_RGS6, 
                                     aes(x = factor(Disease_status, levels = c("AF", "SR")), 
                                         y = Relative_expression_TBPnorm)) +
  stat_summary(aes(fill = Disease_status), fun = mean, geom = "bar", colour = "black", width = 0.4) +
  geom_dotplot(aes(fill = Individual_ID), binaxis = "y", stackdir = "center", dotsize = 1.3) +
    scale_fill_manual(
    name = "Individual_ID",  
    values = c(
      "AF" = "#FDE6E0", 
      "SR" = "#F2FDFD",
      "BCVR-1638" = "#FFF2CC", 
      "BCVR-3913" = "#E4A3A1", 
      "BCVR-5888" = "#c8ee90", 
      "BCVR-5890" = "#56BED3", 
      "BCVR-9231" = "#ffcccb"
    ),
    breaks = c("AF", "SR", "BCVR-9231", "BCVR-1638", "BCVR-3913", "BCVR-5890", "BCVR-5888")
  ) +
  
  scale_y_continuous(name = "Relative Expression") +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 21), 
    axis.title.y = element_text(size = 21),
    axis.text.y = element_text(size = 16)
  )
RTqPCR_220824_RGS6_TBPnorm + geom_signif(
  comparisons = list(c("AF", "SR")), 
  map_signif_level = FALSE, 
  step_increase = 0.1,
  textsize = 6  ,
  y_position = 0.78,
  annotations = "p = 0.047"
)


## ANGPT1

t_test_result <- t.test(data_analysis_dataframe_withdata_ANGPT1_SR$Relative_expression_TBPnorm, data_analysis_dataframe_withdata_ANGPT1_AF$Relative_expression_TBPnorm, alternative = "greater")
t_test_result
#pval  0.7054

RTqPCR_ANGPT1_TBPnorm <- ggplot(data_analysis_dataframe_withdata_ANGPT1, 
                                aes(x = factor(Disease_status, levels = c("AF", "SR")), 
                                    y = Relative_expression_TBPnorm)) +
  stat_summary(aes(fill = Disease_status), fun = mean, geom = "bar", colour = "black", width = 0.4) +
  geom_dotplot(aes(fill = Individual_ID), binaxis = "y", stackdir = "center", dotsize = 1.3) +
    scale_fill_manual(
    name = "Individual_ID",  
    values = c(
      "AF" = "#FDE6E0", 
      "SR" = "#F2FDFD",
      "BCVR-1638" = "#FFF2CC", 
      "BCVR-3913" = "#E4A3A1", 
      "BCVR-5888" = "#c8ee90", 
      "BCVR-5890" = "#56BED3", 
      "BCVR-9231" = "#ffcccb"
    ),
    breaks = c("AF", "SR", "BCVR-9231", "BCVR-1638", "BCVR-3913", "BCVR-5890", "BCVR-5888")
  ) +
  
  scale_y_continuous(name = "Relative Expression") +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 21), 
    axis.title.y = element_text(size = 21),
    axis.text.y = element_text(size = 16) 
  )
RTqPCR_ANGPT1_TBPnorm + geom_signif(
  comparisons = list(c("AF", "SR")), 
  map_signif_level = FALSE, 
  step_increase = 0.1,     
  textsize = 6  ,
  y_position = 1.78, 
  annotations = "p = 0.6425"
)



## BLM 

t_test_result <- t.test(data_analysis_dataframe_withdata_BLM_SR$Relative_expression_GAPDHnorm, data_analysis_dataframe_withdata_BLM_AF$Relative_expression_GAPDHnorm, alternative = "greater")
t_test_result
#pval  0.8303

RTqPCR_BLM_TBPnorm <- ggplot(data_analysis_dataframe_withdata_BLM, 
                             aes(x = factor(Disease_status, levels = c("AF", "SR")), 
                                 y = Relative_expression_TBPnorm)) +
  stat_summary(aes(fill = Disease_status), fun = mean, geom = "bar", colour = "black", width = 0.4) +
  geom_dotplot(aes(fill = Individual_ID), binaxis = "y", stackdir = "center", dotsize = 1.3) +
    scale_fill_manual(
    name = "Individual_ID",  
    values = c(
      "AF" = "#FDE6E0", 
      "SR" = "#F2FDFD",
      "BCVR-1638" = "#FFF2CC", 
      "BCVR-3913" = "#E4A3A1", 
      "BCVR-5888" = "#c8ee90", 
      "BCVR-5890" = "#56BED3", 
      "BCVR-9231" = "#ffcccb"
    ),
    breaks = c("AF", "SR", "BCVR-9231", "BCVR-1638", "BCVR-3913", "BCVR-5890", "BCVR-5888")
  ) +
  
  scale_y_continuous(name = "Relative Expression") +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 21), 
    axis.title.y = element_text(size = 21),
    axis.text.y = element_text(size = 16)
  )
RTqPCR_BLM_TBPnorm + geom_signif(
  comparisons = list(c("AF", "SR")), 
  map_signif_level = FALSE, 
  step_increase = 0.1, 
  textsize = 6  ,
  y_position = 0.13, 
  annotations = "p = 0.8"
)



## NR2F1

t_test_result <- t.test(data_analysis_dataframe_withdata_NR2F1_SR$Relative_expression_TBPnorm, data_analysis_dataframe_withdata_NR2F1_AF$Relative_expression_TBPnorm, alternative = "greater")
t_test_result
#pval  0.2849

RTqPCR_NR2F1_TBPnorm <- ggplot(data_analysis_dataframe_withdata_NR2F1, 
                               aes(x = factor(Disease_status, levels = c("AF", "SR")), 
                                   y = Relative_expression_TBPnorm)) +
  stat_summary(aes(fill = Disease_status), fun = mean, geom = "bar", colour = "black", width = 0.4) +
  geom_dotplot(aes(fill = Individual_ID), binaxis = "y", stackdir = "center", dotsize = 1.3) +
  scale_fill_manual(
    name = "Individual_ID",  
    values = c(
      "AF" = "#FDE6E0", 
      "SR" = "#F2FDFD",
      "BCVR-1638" = "#FFF2CC", 
      "BCVR-3913" = "#E4A3A1", 
      "BCVR-5888" = "#c8ee90", 
      "BCVR-5890" = "#56BED3", 
      "BCVR-9231" = "#ffcccb"
    ),
    breaks = c("AF", "SR", "BCVR-9231", "BCVR-1638", "BCVR-3913", "BCVR-5890", "BCVR-5888")
  ) +
  
  scale_y_continuous(name = "Relative Expression") +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 21), 
    axis.title.y = element_text(size = 21),
    axis.text.y = element_text(size = 16)
  )+ geom_signif(
    comparisons = list(c("AF", "SR")), 
    map_signif_level = FALSE, 
    step_increase = 0.1,     
    textsize = 6  ,
    y_position = 4.05, 
    annotations = "p = 0.28"
  )



###################### ###################### ###################### ###################### ###################### ###################### ###################### ###################### ###################### ###################### ######################

data_analysis_dataframe_withdata <- read.csv("data_analysis_dataframe_ChIPqPCR230625.csv")
data_analysis_dataframe_withdata <- data_analysis_dataframe_withdata[!is.na(data_analysis_dataframe_withdata$Avg_Ct), ]
data_analysis_dataframe_withdata$Anat_side <- NA
data_analysis_dataframe_withdata$Disease_status <- NA
data_analysis_dataframe_withdata$Individual_ID <- NA
data_analysis_dataframe_withdata$Disease_AnatSide <- NA
data_analysis_dataframe_withdata$Anat_side[data_analysis_dataframe_withdata$Sample %in% c("Sample1","Sample3","Sample5", "Sample7", "Sample9")] <- "RA"
data_analysis_dataframe_withdata$Anat_side[data_analysis_dataframe_withdata$Sample %in% c("Sample2","Sample4","Sample6", "Sample8", "Sample10")] <- "LA"
data_analysis_dataframe_withdata$Disease_status[data_analysis_dataframe_withdata$Sample %in% c("Sample1","Sample2","Sample7", "Sample8", "Sample9","Sample10")] <- "AF"
data_analysis_dataframe_withdata$Disease_status[data_analysis_dataframe_withdata$Sample %in% c("Sample3","Sample4","Sample5", "Sample6")] <- "SR"
data_analysis_dataframe_withdata$Individual_ID[data_analysis_dataframe_withdata$Sample %in% c("Sample1")] <- "BCVR-9231"
data_analysis_dataframe_withdata$Individual_ID[data_analysis_dataframe_withdata$Sample %in% c("Sample2")] <- "BCVR-9231"
data_analysis_dataframe_withdata$Individual_ID[data_analysis_dataframe_withdata$Sample %in% c("Sample3")] <- "BCVR-5890"
data_analysis_dataframe_withdata$Individual_ID[data_analysis_dataframe_withdata$Sample %in% c("Sample4")] <- "BCVR-5890"
data_analysis_dataframe_withdata$Individual_ID[data_analysis_dataframe_withdata$Sample %in% c("Sample5")] <- "BCVR-5888"
data_analysis_dataframe_withdata$Individual_ID[data_analysis_dataframe_withdata$Sample %in% c("Sample6")] <- "BCVR-5888"
data_analysis_dataframe_withdata$Individual_ID[data_analysis_dataframe_withdata$Sample %in% c("Sample7")] <- "BCVR-3913"
data_analysis_dataframe_withdata$Individual_ID[data_analysis_dataframe_withdata$Sample %in% c("Sample8")] <- "BCVR-3913"
data_analysis_dataframe_withdata$Individual_ID[data_analysis_dataframe_withdata$Sample %in% c("Sample9")] <- "BCVR-1638"
data_analysis_dataframe_withdata$Individual_ID[data_analysis_dataframe_withdata$Sample %in% c("Sample10")] <- "BCVR-1638"
data_analysis_dataframe_withdata$Disease_AnatSide[data_analysis_dataframe_withdata$Sample %in% c("Sample2","Sample8","Sample10")] <- "AF-LA"
data_analysis_dataframe_withdata$Disease_AnatSide[data_analysis_dataframe_withdata$Sample %in% c("Sample1","Sample7","Sample9")] <- "AF-RA"
data_analysis_dataframe_withdata$Disease_AnatSide[data_analysis_dataframe_withdata$Sample %in% c("Sample4","Sample6")] <- "SR-LA"
data_analysis_dataframe_withdata$Disease_AnatSide[data_analysis_dataframe_withdata$Sample %in% c("Sample3","Sample5")] <- "SR-RA"
data_analysis_dataframe_withdata_SCX_A <- data_analysis_dataframe_withdata[which (data_analysis_dataframe_withdata$Marker == "SCX_A"),]
data_analysis_dataframe_withdata_SCX_B <- data_analysis_dataframe_withdata[which (data_analysis_dataframe_withdata$Marker == "SCX_B"),]
data_analysis_dataframe_withdata_NPPB_AF008 <- data_analysis_dataframe_withdata[which (data_analysis_dataframe_withdata$Marker == "NPPB_AF008"),]
data_analysis_dataframe_withdata_SCX_A <- data_analysis_dataframe_withdata_SCX_A[!is.na(data_analysis_dataframe_withdata_SCX_A$Norm_input), ]
data_analysis_dataframe_withdata_SCX_B <- data_analysis_dataframe_withdata_SCX_B[!is.na(data_analysis_dataframe_withdata_SCX_B$Norm_input), ]
data_analysis_dataframe_withdata_NPPB_AF008 <- data_analysis_dataframe_withdata_NPPB_AF008[!is.na(data_analysis_dataframe_withdata_NPPB_AF008$Norm_input), ]



data_analysis_dataframe_withdata <- read.csv("data_analysis_dataframe_plate6.csv")
data_analysis_dataframe_withdata <- data_analysis_dataframe_withdata[!is.na(data_analysis_dataframe_withdata$Avg_Ct), ]
data_analysis_dataframe_withdata$Anat_side <- NA
data_analysis_dataframe_withdata$Disease_status <- NA
data_analysis_dataframe_withdata$Individual_ID <- NA
data_analysis_dataframe_withdata$Disease_AnatSide <- NA
data_analysis_dataframe_withdata$Anat_side[data_analysis_dataframe_withdata$Sample %in% c("Sample1","Sample3","Sample5", "Sample7", "Sample9")] <- "RA"
data_analysis_dataframe_withdata$Anat_side[data_analysis_dataframe_withdata$Sample %in% c("Sample2","Sample4","Sample6", "Sample8", "Sample10")] <- "LA"
data_analysis_dataframe_withdata$Disease_status[data_analysis_dataframe_withdata$Sample %in% c("Sample1","Sample2","Sample7", "Sample8", "Sample9","Sample10")] <- "AF"
data_analysis_dataframe_withdata$Disease_status[data_analysis_dataframe_withdata$Sample %in% c("Sample3","Sample4","Sample5", "Sample6")] <- "SR"
data_analysis_dataframe_withdata$Individual_ID[data_analysis_dataframe_withdata$Sample %in% c("Sample1")] <- "BCVR-9231"
data_analysis_dataframe_withdata$Individual_ID[data_analysis_dataframe_withdata$Sample %in% c("Sample2")] <- "BCVR-9231"
data_analysis_dataframe_withdata$Individual_ID[data_analysis_dataframe_withdata$Sample %in% c("Sample3")] <- "BCVR-5890"
data_analysis_dataframe_withdata$Individual_ID[data_analysis_dataframe_withdata$Sample %in% c("Sample4")] <- "BCVR-5890"
data_analysis_dataframe_withdata$Individual_ID[data_analysis_dataframe_withdata$Sample %in% c("Sample5")] <- "BCVR-5888"
data_analysis_dataframe_withdata$Individual_ID[data_analysis_dataframe_withdata$Sample %in% c("Sample6")] <- "BCVR-5888"
data_analysis_dataframe_withdata$Individual_ID[data_analysis_dataframe_withdata$Sample %in% c("Sample7")] <- "BCVR-3913"
data_analysis_dataframe_withdata$Individual_ID[data_analysis_dataframe_withdata$Sample %in% c("Sample8")] <- "BCVR-3913"
data_analysis_dataframe_withdata$Individual_ID[data_analysis_dataframe_withdata$Sample %in% c("Sample9")] <- "BCVR-1638"
data_analysis_dataframe_withdata$Individual_ID[data_analysis_dataframe_withdata$Sample %in% c("Sample10")] <- "BCVR-1638"
data_analysis_dataframe_withdata$Disease_AnatSide[data_analysis_dataframe_withdata$Sample %in% c("Sample2","Sample8","Sample10")] <- "AF-LA"
data_analysis_dataframe_withdata$Disease_AnatSide[data_analysis_dataframe_withdata$Sample %in% c("Sample1","Sample7","Sample9")] <- "AF-RA"
data_analysis_dataframe_withdata$Disease_AnatSide[data_analysis_dataframe_withdata$Sample %in% c("Sample4","Sample6")] <- "SR-LA"
data_analysis_dataframe_withdata$Disease_AnatSide[data_analysis_dataframe_withdata$Sample %in% c("Sample3","Sample5")] <- "SR-RA"
data_analysis_dataframe_withdata_NR2F1_B <- data_analysis_dataframe_withdata[which (data_analysis_dataframe_withdata$Marker == "NR2F1_B"),]
data_analysis_dataframe_withdata_STRN <- data_analysis_dataframe_withdata[which (data_analysis_dataframe_withdata$Marker == "STRN"),]
data_analysis_dataframe_withdata_VASH1_B <- data_analysis_dataframe_withdata[which (data_analysis_dataframe_withdata$Marker == "VASH1_B"),]
data_analysis_dataframe_withdata_NR2F1_B <- data_analysis_dataframe_withdata_NR2F1_B[!is.na(data_analysis_dataframe_withdata_NR2F1_B$Norm_input), ]
data_analysis_dataframe_withdata_STRN <- data_analysis_dataframe_withdata_STRN[!is.na(data_analysis_dataframe_withdata_STRN$Norm_input), ]
data_analysis_dataframe_withdata_VASH1_B <- data_analysis_dataframe_withdata_VASH1_B[!is.na(data_analysis_dataframe_withdata_VASH1_B$Norm_input), ]


data_analysis_dataframe_withdata <- read.csv("data_analysis_dataframe_plate2.csv")
data_analysis_dataframe_withdata <- data_analysis_dataframe_withdata[!is.na(data_analysis_dataframe_withdata$Avg_Ct), ]
data_analysis_dataframe_withdata$Anat_side <- NA
data_analysis_dataframe_withdata$Disease_status <- NA
data_analysis_dataframe_withdata$Individual_ID <- NA
data_analysis_dataframe_withdata$Disease_AnatSide <- NA
data_analysis_dataframe_withdata$Anat_side[data_analysis_dataframe_withdata$Sample %in% c("Sample1","Sample3","Sample5", "Sample7", "Sample9")] <- "RA"
data_analysis_dataframe_withdata$Anat_side[data_analysis_dataframe_withdata$Sample %in% c("Sample2","Sample4","Sample6", "Sample8", "Sample10")] <- "LA"
data_analysis_dataframe_withdata$Disease_status[data_analysis_dataframe_withdata$Sample %in% c("Sample1","Sample2","Sample7", "Sample8", "Sample9","Sample10")] <- "AF"
data_analysis_dataframe_withdata$Disease_status[data_analysis_dataframe_withdata$Sample %in% c("Sample3","Sample4","Sample5", "Sample6")] <- "SR"
data_analysis_dataframe_withdata$Individual_ID[data_analysis_dataframe_withdata$Sample %in% c("Sample1")] <- "BCVR-9231"
data_analysis_dataframe_withdata$Individual_ID[data_analysis_dataframe_withdata$Sample %in% c("Sample2")] <- "BCVR-9231"
data_analysis_dataframe_withdata$Individual_ID[data_analysis_dataframe_withdata$Sample %in% c("Sample3")] <- "BCVR-5890"
data_analysis_dataframe_withdata$Individual_ID[data_analysis_dataframe_withdata$Sample %in% c("Sample4")] <- "BCVR-5890"
data_analysis_dataframe_withdata$Individual_ID[data_analysis_dataframe_withdata$Sample %in% c("Sample5")] <- "BCVR-5888"
data_analysis_dataframe_withdata$Individual_ID[data_analysis_dataframe_withdata$Sample %in% c("Sample6")] <- "BCVR-5888"
data_analysis_dataframe_withdata$Individual_ID[data_analysis_dataframe_withdata$Sample %in% c("Sample7")] <- "BCVR-3913"
data_analysis_dataframe_withdata$Individual_ID[data_analysis_dataframe_withdata$Sample %in% c("Sample8")] <- "BCVR-3913"
data_analysis_dataframe_withdata$Individual_ID[data_analysis_dataframe_withdata$Sample %in% c("Sample9")] <- "BCVR-1638"
data_analysis_dataframe_withdata$Individual_ID[data_analysis_dataframe_withdata$Sample %in% c("Sample10")] <- "BCVR-1638"
data_analysis_dataframe_withdata$Disease_AnatSide[data_analysis_dataframe_withdata$Sample %in% c("Sample2","Sample8","Sample10")] <- "AF-LA"
data_analysis_dataframe_withdata$Disease_AnatSide[data_analysis_dataframe_withdata$Sample %in% c("Sample1","Sample7","Sample9")] <- "AF-RA"
data_analysis_dataframe_withdata$Disease_AnatSide[data_analysis_dataframe_withdata$Sample %in% c("Sample4","Sample6")] <- "SR-LA"
data_analysis_dataframe_withdata$Disease_AnatSide[data_analysis_dataframe_withdata$Sample %in% c("Sample3","Sample5")] <- "SR-RA"
data_analysis_dataframe_withdata_LRRC32_AFLA0328 <- data_analysis_dataframe_withdata[which (data_analysis_dataframe_withdata$Marker == "LRRC32_AFLA0328"),]
data_analysis_dataframe_withdata_BLM_SRLA617 <- data_analysis_dataframe_withdata[which (data_analysis_dataframe_withdata$Marker == "BLM_SRLA617"),]
data_analysis_dataframe_withdata_LRRC32_AFLA0328 <- data_analysis_dataframe_withdata_LRRC32_AFLA0328[!is.na(data_analysis_dataframe_withdata_LRRC32_AFLA0328$Norm_input), ]
data_analysis_dataframe_withdata_BLM_SRLA617 <- data_analysis_dataframe_withdata_BLM_SRLA617[!is.na(data_analysis_dataframe_withdata_BLM_SRLA617$Norm_input), ]


data_analysis_dataframe_withdata <- read.csv("data_analysis_dataframe_plate3.csv")
data_analysis_dataframe_withdata <- data_analysis_dataframe_withdata[!is.na(data_analysis_dataframe_withdata$Avg_Ct), ]
data_analysis_dataframe_withdata$Anat_side <- NA
data_analysis_dataframe_withdata$Disease_status <- NA
data_analysis_dataframe_withdata$Individual_ID <- NA
data_analysis_dataframe_withdata$Disease_AnatSide <- NA
data_analysis_dataframe_withdata$Anat_side[data_analysis_dataframe_withdata$Sample %in% c("Sample1","Sample3","Sample5", "Sample7", "Sample9")] <- "RA"
data_analysis_dataframe_withdata$Anat_side[data_analysis_dataframe_withdata$Sample %in% c("Sample2","Sample4","Sample6", "Sample8", "Sample10")] <- "LA"
data_analysis_dataframe_withdata$Disease_status[data_analysis_dataframe_withdata$Sample %in% c("Sample1","Sample2","Sample7", "Sample8", "Sample9","Sample10")] <- "AF"
data_analysis_dataframe_withdata$Disease_status[data_analysis_dataframe_withdata$Sample %in% c("Sample3","Sample4","Sample5", "Sample6")] <- "SR"
data_analysis_dataframe_withdata$Individual_ID[data_analysis_dataframe_withdata$Sample %in% c("Sample1")] <- "BCVR-9231"
data_analysis_dataframe_withdata$Individual_ID[data_analysis_dataframe_withdata$Sample %in% c("Sample2")] <- "BCVR-9231"
data_analysis_dataframe_withdata$Individual_ID[data_analysis_dataframe_withdata$Sample %in% c("Sample3")] <- "BCVR-5890"
data_analysis_dataframe_withdata$Individual_ID[data_analysis_dataframe_withdata$Sample %in% c("Sample4")] <- "BCVR-5890"
data_analysis_dataframe_withdata$Individual_ID[data_analysis_dataframe_withdata$Sample %in% c("Sample5")] <- "BCVR-5888"
data_analysis_dataframe_withdata$Individual_ID[data_analysis_dataframe_withdata$Sample %in% c("Sample6")] <- "BCVR-5888"
data_analysis_dataframe_withdata$Individual_ID[data_analysis_dataframe_withdata$Sample %in% c("Sample7")] <- "BCVR-3913"
data_analysis_dataframe_withdata$Individual_ID[data_analysis_dataframe_withdata$Sample %in% c("Sample8")] <- "BCVR-3913"
data_analysis_dataframe_withdata$Individual_ID[data_analysis_dataframe_withdata$Sample %in% c("Sample9")] <- "BCVR-1638"
data_analysis_dataframe_withdata$Individual_ID[data_analysis_dataframe_withdata$Sample %in% c("Sample10")] <- "BCVR-1638"
data_analysis_dataframe_withdata$Disease_AnatSide[data_analysis_dataframe_withdata$Sample %in% c("Sample2","Sample8","Sample10")] <- "AF-LA"
data_analysis_dataframe_withdata$Disease_AnatSide[data_analysis_dataframe_withdata$Sample %in% c("Sample1","Sample7","Sample9")] <- "AF-RA"
data_analysis_dataframe_withdata$Disease_AnatSide[data_analysis_dataframe_withdata$Sample %in% c("Sample4","Sample6")] <- "SR-LA"
data_analysis_dataframe_withdata$Disease_AnatSide[data_analysis_dataframe_withdata$Sample %in% c("Sample3","Sample5")] <- "SR-RA"
data_analysis_dataframe_withdata_ANGPT1_SRRA306 <- data_analysis_dataframe_withdata[which (data_analysis_dataframe_withdata$Marker == "ANGPT1_SRRA306"),]
data_analysis_dataframe_withdata_ANGPTL2_AFRA1622 <- data_analysis_dataframe_withdata[which (data_analysis_dataframe_withdata$Marker == "ANGPTL2_AFRA1622"),]
data_analysis_dataframe_withdata_ANGPT1_SRRA306 <- data_analysis_dataframe_withdata_ANGPT1_SRRA306[!is.na(data_analysis_dataframe_withdata_ANGPT1_SRRA306$Norm_input), ]
data_analysis_dataframe_withdata_ANGPTL2_AFRA1622 <- data_analysis_dataframe_withdata_ANGPTL2_AFRA1622[!is.na(data_analysis_dataframe_withdata_ANGPTL2_AFRA1622$Norm_input), ]



data_analysis_dataframe_withdata <- read.csv("data_analysis_dataframe_plate4_b.csv")
data_analysis_dataframe_withdata <- data_analysis_dataframe_withdata[!is.na(data_analysis_dataframe_withdata$Avg_Ct), ]
data_analysis_dataframe_withdata$Anat_side <- NA
data_analysis_dataframe_withdata$Disease_status <- NA
data_analysis_dataframe_withdata$Individual_ID <- NA
data_analysis_dataframe_withdata$Disease_AnatSide <- NA
data_analysis_dataframe_withdata$Anat_side[data_analysis_dataframe_withdata$Sample %in% c("Sample1","Sample3","Sample5", "Sample7", "Sample9")] <- "RA"
data_analysis_dataframe_withdata$Anat_side[data_analysis_dataframe_withdata$Sample %in% c("Sample2","Sample4","Sample6", "Sample8", "Sample10")] <- "LA"
data_analysis_dataframe_withdata$Disease_status[data_analysis_dataframe_withdata$Sample %in% c("Sample1","Sample2","Sample7", "Sample8", "Sample9","Sample10")] <- "AF"
data_analysis_dataframe_withdata$Disease_status[data_analysis_dataframe_withdata$Sample %in% c("Sample3","Sample4","Sample5", "Sample6")] <- "SR"
data_analysis_dataframe_withdata$Individual_ID[data_analysis_dataframe_withdata$Sample %in% c("Sample1")] <- "BCVR-9231"
data_analysis_dataframe_withdata$Individual_ID[data_analysis_dataframe_withdata$Sample %in% c("Sample2")] <- "BCVR-9231"
data_analysis_dataframe_withdata$Individual_ID[data_analysis_dataframe_withdata$Sample %in% c("Sample3")] <- "BCVR-5890"
data_analysis_dataframe_withdata$Individual_ID[data_analysis_dataframe_withdata$Sample %in% c("Sample4")] <- "BCVR-5890"
data_analysis_dataframe_withdata$Individual_ID[data_analysis_dataframe_withdata$Sample %in% c("Sample5")] <- "BCVR-5888"
data_analysis_dataframe_withdata$Individual_ID[data_analysis_dataframe_withdata$Sample %in% c("Sample6")] <- "BCVR-5888"
data_analysis_dataframe_withdata$Individual_ID[data_analysis_dataframe_withdata$Sample %in% c("Sample7")] <- "BCVR-3913"
data_analysis_dataframe_withdata$Individual_ID[data_analysis_dataframe_withdata$Sample %in% c("Sample8")] <- "BCVR-3913"
data_analysis_dataframe_withdata$Individual_ID[data_analysis_dataframe_withdata$Sample %in% c("Sample9")] <- "BCVR-1638"
data_analysis_dataframe_withdata$Individual_ID[data_analysis_dataframe_withdata$Sample %in% c("Sample10")] <- "BCVR-1638"
data_analysis_dataframe_withdata$Disease_AnatSide[data_analysis_dataframe_withdata$Sample %in% c("Sample2","Sample8","Sample10")] <- "AF-LA"
data_analysis_dataframe_withdata$Disease_AnatSide[data_analysis_dataframe_withdata$Sample %in% c("Sample1","Sample7","Sample9")] <- "AF-RA"
data_analysis_dataframe_withdata$Disease_AnatSide[data_analysis_dataframe_withdata$Sample %in% c("Sample4","Sample6")] <- "SR-LA"
data_analysis_dataframe_withdata$Disease_AnatSide[data_analysis_dataframe_withdata$Sample %in% c("Sample3","Sample5")] <- "SR-RA"
data_analysis_dataframe_withdata_MT1X_11 <- data_analysis_dataframe_withdata[which (data_analysis_dataframe_withdata$Marker == "MT1X_11"),]
data_analysis_dataframe_withdata_RGS6_443 <- data_analysis_dataframe_withdata[which (data_analysis_dataframe_withdata$Marker == "RGS6"),]
data_analysis_dataframe_withdata_MT1X_11 <- data_analysis_dataframe_withdata_MT1X_11[!is.na(data_analysis_dataframe_withdata_MT1X_11$Norm_input), ]
data_analysis_dataframe_withdata_RGS6_443 <- data_analysis_dataframe_withdata_RGS6_443[!is.na(data_analysis_dataframe_withdata_RGS6_443$Norm_input), ]



data_analysis_dataframe_withdata <- read.csv("data_analysis_dataframe_plate5.csv")
data_analysis_dataframe_withdata <- data_analysis_dataframe_withdata[!is.na(data_analysis_dataframe_withdata$Avg_Ct), ]
data_analysis_dataframe_withdata$Anat_side <- NA
data_analysis_dataframe_withdata$Disease_status <- NA
data_analysis_dataframe_withdata$Individual_ID <- NA
data_analysis_dataframe_withdata$Disease_AnatSide <- NA
data_analysis_dataframe_withdata$Anat_side[data_analysis_dataframe_withdata$Sample %in% c("Sample1","Sample3","Sample5", "Sample7", "Sample9")] <- "RA"
data_analysis_dataframe_withdata$Anat_side[data_analysis_dataframe_withdata$Sample %in% c("Sample2","Sample4","Sample6", "Sample8", "Sample10")] <- "LA"
data_analysis_dataframe_withdata$Disease_status[data_analysis_dataframe_withdata$Sample %in% c("Sample1","Sample2","Sample7", "Sample8", "Sample9","Sample10")] <- "AF"
data_analysis_dataframe_withdata$Disease_status[data_analysis_dataframe_withdata$Sample %in% c("Sample3","Sample4","Sample5", "Sample6")] <- "SR"
data_analysis_dataframe_withdata$Individual_ID[data_analysis_dataframe_withdata$Sample %in% c("Sample1")] <- "BCVR-9231"
data_analysis_dataframe_withdata$Individual_ID[data_analysis_dataframe_withdata$Sample %in% c("Sample2")] <- "BCVR-9231"
data_analysis_dataframe_withdata$Individual_ID[data_analysis_dataframe_withdata$Sample %in% c("Sample3")] <- "BCVR-5890"
data_analysis_dataframe_withdata$Individual_ID[data_analysis_dataframe_withdata$Sample %in% c("Sample4")] <- "BCVR-5890"
data_analysis_dataframe_withdata$Individual_ID[data_analysis_dataframe_withdata$Sample %in% c("Sample5")] <- "BCVR-5888"
data_analysis_dataframe_withdata$Individual_ID[data_analysis_dataframe_withdata$Sample %in% c("Sample6")] <- "BCVR-5888"
data_analysis_dataframe_withdata$Individual_ID[data_analysis_dataframe_withdata$Sample %in% c("Sample7")] <- "BCVR-3913"
data_analysis_dataframe_withdata$Individual_ID[data_analysis_dataframe_withdata$Sample %in% c("Sample8")] <- "BCVR-3913"
data_analysis_dataframe_withdata$Individual_ID[data_analysis_dataframe_withdata$Sample %in% c("Sample9")] <- "BCVR-1638"
data_analysis_dataframe_withdata$Individual_ID[data_analysis_dataframe_withdata$Sample %in% c("Sample10")] <- "BCVR-1638"
data_analysis_dataframe_withdata$Disease_AnatSide[data_analysis_dataframe_withdata$Sample %in% c("Sample2","Sample8","Sample10")] <- "AF-LA"
data_analysis_dataframe_withdata$Disease_AnatSide[data_analysis_dataframe_withdata$Sample %in% c("Sample1","Sample7","Sample9")] <- "AF-RA"
data_analysis_dataframe_withdata$Disease_AnatSide[data_analysis_dataframe_withdata$Sample %in% c("Sample4","Sample6")] <- "SR-LA"
data_analysis_dataframe_withdata$Disease_AnatSide[data_analysis_dataframe_withdata$Sample %in% c("Sample3","Sample5")] <- "SR-RA"
data_analysis_dataframe_withdata_SLC43A3 <- data_analysis_dataframe_withdata[which (data_analysis_dataframe_withdata$Marker == "SLC43A3"),]
data_analysis_dataframe_withdata_FOXS1 <- data_analysis_dataframe_withdata[which (data_analysis_dataframe_withdata$Marker == "FOXS1"),]
data_analysis_dataframe_withdata_FGF1 <- data_analysis_dataframe_withdata[which (data_analysis_dataframe_withdata$Marker == "FGF1"),]
data_analysis_dataframe_withdata_SLC43A3 <- data_analysis_dataframe_withdata_SLC43A3[!is.na(data_analysis_dataframe_withdata_SLC43A3$Norm_input), ]
data_analysis_dataframe_withdata_FOXS1 <- data_analysis_dataframe_withdata_FOXS1[!is.na(data_analysis_dataframe_withdata_FOXS1$Norm_input), ]




data_analysis_dataframe_withdata <- read.csv("data_analysis_dataframe.csv")
data_analysis_dataframe_withdata <- data_analysis_dataframe_withdata[!is.na(data_analysis_dataframe_withdata$Avg_Ct), ]
data_analysis_dataframe_withdata$Anat_side <- NA
data_analysis_dataframe_withdata$Disease_status <- NA
data_analysis_dataframe_withdata$Individual_ID <- NA
data_analysis_dataframe_withdata$Disease_AnatSide <- NA
data_analysis_dataframe_withdata$Anat_side[data_analysis_dataframe_withdata$Sample %in% c("Sample1","Sample3","Sample5", "Sample7", "Sample9")] <- "RA"
data_analysis_dataframe_withdata$Anat_side[data_analysis_dataframe_withdata$Sample %in% c("Sample2","Sample4","Sample6", "Sample8", "Sample10")] <- "LA"
data_analysis_dataframe_withdata$Disease_status[data_analysis_dataframe_withdata$Sample %in% c("Sample1","Sample2","Sample7", "Sample8", "Sample9","Sample10")] <- "AF"
data_analysis_dataframe_withdata$Disease_status[data_analysis_dataframe_withdata$Sample %in% c("Sample3","Sample4","Sample5", "Sample6")] <- "SR"
data_analysis_dataframe_withdata$Individual_ID[data_analysis_dataframe_withdata$Sample %in% c("Sample1")] <- "BCVR-9231"
data_analysis_dataframe_withdata$Individual_ID[data_analysis_dataframe_withdata$Sample %in% c("Sample2")] <- "BCVR-9231"
data_analysis_dataframe_withdata$Individual_ID[data_analysis_dataframe_withdata$Sample %in% c("Sample3")] <- "BCVR-5890"
data_analysis_dataframe_withdata$Individual_ID[data_analysis_dataframe_withdata$Sample %in% c("Sample4")] <- "BCVR-5890"
data_analysis_dataframe_withdata$Individual_ID[data_analysis_dataframe_withdata$Sample %in% c("Sample5")] <- "BCVR-5888"
data_analysis_dataframe_withdata$Individual_ID[data_analysis_dataframe_withdata$Sample %in% c("Sample6")] <- "BCVR-5888"
data_analysis_dataframe_withdata$Individual_ID[data_analysis_dataframe_withdata$Sample %in% c("Sample7")] <- "BCVR-3913"
data_analysis_dataframe_withdata$Individual_ID[data_analysis_dataframe_withdata$Sample %in% c("Sample8")] <- "BCVR-3913"
data_analysis_dataframe_withdata$Individual_ID[data_analysis_dataframe_withdata$Sample %in% c("Sample9")] <- "BCVR-1638"
data_analysis_dataframe_withdata$Individual_ID[data_analysis_dataframe_withdata$Sample %in% c("Sample10")] <- "BCVR-1638"
data_analysis_dataframe_withdata$Disease_AnatSide[data_analysis_dataframe_withdata$Sample %in% c("Sample2","Sample8","Sample10")] <- "AF-LA"
data_analysis_dataframe_withdata$Disease_AnatSide[data_analysis_dataframe_withdata$Sample %in% c("Sample1","Sample7","Sample9")] <- "AF-RA"
data_analysis_dataframe_withdata$Disease_AnatSide[data_analysis_dataframe_withdata$Sample %in% c("Sample4","Sample6")] <- "SR-LA"
data_analysis_dataframe_withdata$Disease_AnatSide[data_analysis_dataframe_withdata$Sample %in% c("Sample3","Sample5")] <- "SR-RA"
data_analysis_dataframe_withdata_NPPB_genebody <- data_analysis_dataframe_withdata[which (data_analysis_dataframe_withdata$Marker == "NPPB_genebody"),]
data_analysis_dataframe_withdata_NPPB_AFLA0013 <- data_analysis_dataframe_withdata[which (data_analysis_dataframe_withdata$Marker == "NPPB_AFLA0013"),]
data_analysis_dataframe_withdata_NPPB_genebody <- data_analysis_dataframe_withdata_NPPB_genebody[!is.na(data_analysis_dataframe_withdata_NPPB_genebody$Norm_input), ]
data_analysis_dataframe_withdata_NPPB_AFLA0013 <- data_analysis_dataframe_withdata_NPPB_AFLA0013[!is.na(data_analysis_dataframe_withdata_NPPB_AFLA0013$Norm_input), ]


################ ChIP-qPCR ######################


# 📊 Figures S5c and S5d


#### Regions detected as AF-enriched in main cohort 

### VASH1 


t_test_result <- t.test(data_analysis_dataframe_withdata_VASH1_B_AF$RelEnrich_Norm_GAPDH_Input, data_analysis_dataframe_withdata_VASH1_B_SR$RelEnrich_Norm_GAPDH_Input, alternative = "greater")
t_test_result
#pval  0.2829


ChIPqPCR_VASH1_B_GAPDHInputnorm <- ggplot(data_analysis_dataframe_withdata_VASH1_B,    # 📊 Figures 5a
                                          aes(x = factor(Disease_status, levels = c("AF", "SR")), 
                                              y = RelEnrich_Norm_GAPDH_Input)) +
  stat_summary(aes(fill = Disease_status), fun = mean, geom = "bar", colour = "black", width = 0.4) +
  geom_dotplot(aes(fill = Individual_ID), binaxis = "y", stackdir = "center", dotsize = 1.3) +
    scale_fill_manual(
    name = "Individual_ID",  
    values = c(
      "AF" = "#FDE6E0", 
      "SR" = "#F2FDFD",
      "BCVR-1638" = "#FFF2CC", 
      "BCVR-3913" = "#E4A3A1", 
      "BCVR-5888" = "#c8ee90", 
      "BCVR-5890" = "#56BED3", 
      "BCVR-9231" = "#ffcccb"
    ),
    breaks = c("AF", "SR", "BCVR-9231", "BCVR-1638", "BCVR-3913", "BCVR-5890", "BCVR-5888")
  ) +
  
  scale_y_continuous(name = "Normalized ChIP Enrichment") +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 21), 
    axis.title.y = element_text(size = 21)
  )+ geom_signif(
    comparisons = list(c("AF", "SR")), 
    map_signif_level = FALSE,
    step_increase = 0.1,    
    textsize = 4  ,
    y_position = 4.3, # Adjust text size
    annotations = "p = 0.2829"
  )


### ANGPTL2


t_test_result <- t.test(data_analysis_dataframe_withdata_ANGPTL2_AFRA1622_AF$RelEnrich_Norm_GAPDH_Input, data_analysis_dataframe_withdata_ANGPTL2_AFRA1622_SR$RelEnrich_Norm_GAPDH_Input, alternative = "greater")
t_test_result
#pval  0.1459


ChIPqPCR_ANGPTL2_AFRA1622_GAPDHInputnorm <- ggplot(data_analysis_dataframe_withdata_ANGPTL2_AFRA1622, 
                                                   aes(x = factor(Disease_status, levels = c("AF", "SR")), 
                                                       y = RelEnrich_Norm_GAPDH_Input)) +
  stat_summary(aes(fill = Disease_status), fun = mean, geom = "bar", colour = "black", width = 0.4) +
  geom_dotplot(aes(fill = Individual_ID), binaxis = "y", stackdir = "center", dotsize = 1.3) +
    scale_fill_manual(
    name = "Individual_ID",  
    values = c(
      "AF" = "#FDE6E0", 
      "SR" = "#F2FDFD",
      "BCVR-1638" = "#FFF2CC", 
      "BCVR-3913" = "#E4A3A1", 
      "BCVR-5888" = "#c8ee90", 
      "BCVR-5890" = "#56BED3", 
      "BCVR-9231" = "#ffcccb"
    ),
    breaks = c("AF", "SR", "BCVR-9231", "BCVR-1638", "BCVR-3913", "BCVR-5890", "BCVR-5888")
  ) +
  
  scale_y_continuous(name = "Normalized ChIP Enrichment") +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 14), 
    axis.title.y = element_text(size = 14),
    axis.text.y = element_text(size = 14) 
  ) 


ChIPqPCR_ANGPTL2_AFRA1622_GAPDHInputnorm + geom_signif(
  comparisons = list(c("AF", "SR")),
  map_signif_level = FALSE, 
  step_increase = 0.1,     
  textsize = 4  ,
  y_position = 4.8, 
  annotations = "p = 0.1459")



### SLC43A3

t_test_result <- t.test(data_analysis_dataframe_withdata_SLC43A3_AF$RelEnrich_Norm_GAPDH_Input, data_analysis_dataframe_withdata_SLC43A3_SR$RelEnrich_Norm_GAPDH_Input, alternative = "greater")
t_test_result
#pval  0.2201

ChIPqPCR_SLC43A3_GAPDHInputnorm <- ggplot(data_analysis_dataframe_withdata_SLC43A3, 
                                          aes(x = factor(Disease_status, levels = c("AF", "SR")), 
                                              y = RelEnrich_Norm_GAPDH_Input)) +
  stat_summary(aes(fill = Disease_status), fun = mean, geom = "bar", colour = "black", width = 0.4) +
  geom_dotplot(aes(fill = Individual_ID), binaxis = "y", stackdir = "center", dotsize = 1.3) +
    scale_fill_manual(
    name = "Individual_ID",  
    values = c(
      "AF" = "#FDE6E0", 
      "SR" = "#F2FDFD",
      "BCVR-1638" = "#FFF2CC", 
      "BCVR-3913" = "#E4A3A1", 
      "BCVR-5888" = "#c8ee90", 
      "BCVR-5890" = "#56BED3", 
      "BCVR-9231" = "#ffcccb"
    ),
    breaks = c("AF", "SR", "BCVR-9231", "BCVR-1638", "BCVR-3913", "BCVR-5890", "BCVR-5888")
  ) +
  
  scale_y_continuous(name = "Normalized ChIP Enrichment") +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 14), 
    axis.title.y = element_text(size = 14),
    axis.text.y = element_text(size = 14) 
  ) 


ChIPqPCR_SLC43A3_GAPDHInputnorm + geom_signif(
  comparisons = list(c("AF", "SR")), 
  map_signif_level = FALSE, 
  step_increase = 0.1,     
  textsize = 4  ,
  y_position = 1, # Adjust text size
  annotations = "p = 0.2201")



### FOXS1


t_test_result <- t.test(data_analysis_dataframe_withdata_FOXS1_AF$RelEnrich_Norm_GAPDH_Input, data_analysis_dataframe_withdata_FOXS1_SR$RelEnrich_Norm_GAPDH_Input, alternative = "greater")
t_test_result
#pval  0.2444

ChIPqPCR_FOXS1_GAPDHInputnorm <- ggplot(data_analysis_dataframe_withdata_FOXS1, 
                                        aes(x = factor(Disease_status, levels = c("AF", "SR")), 
                                            y = RelEnrich_Norm_GAPDH_Input)) +
  stat_summary(aes(fill = Disease_status), fun = mean, geom = "bar", colour = "black", width = 0.4) +
  geom_dotplot(aes(fill = Individual_ID), binaxis = "y", stackdir = "center", dotsize = 1.6) +
    scale_fill_manual(
    name = "Individual_ID",  
    values = c(
      "AF" = "#FDE6E0", 
      "SR" = "#F2FDFD",
      "BCVR-1638" = "#FFF2CC", 
      "BCVR-3913" = "#E4A3A1", 
      "BCVR-5888" = "#c8ee90", 
      "BCVR-5890" = "#56BED3", 
      "BCVR-9231" = "#ffcccb"
    ),
    breaks = c("AF", "SR", "BCVR-9231", "BCVR-1638", "BCVR-3913", "BCVR-5890", "BCVR-5888")
  ) +
  
  scale_y_continuous(name = "Normalized ChIP Enrichment") +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 14), 
    axis.title.y = element_text(size = 14),
    axis.text.y = element_text(size = 14) 
  ) 


ChIPqPCR_FOXS1_GAPDHInputnorm + geom_signif(
  comparisons = list(c("AF", "SR")),
  map_signif_level = FALSE, 
  step_increase = 0.1,    
  textsize = 4  ,
  y_position = 0.79, 
  annotations = "p = 0.2444")



### STRN

t_test_result <- t.test(data_analysis_dataframe_withdata_STRN_AF$RelEnrich_Norm_GAPDH_Input, data_analysis_dataframe_withdata_STRN_SR$RelEnrich_Norm_GAPDH_Input, alternative = "greater")
t_test_result
#pval  0.2727

ChIPqPCR_STRN_GAPDHInputnorm <- ggplot(data_analysis_dataframe_withdata_STRN, 
                                       aes(x = factor(Disease_status, levels = c("AF", "SR")), 
                                           y = RelEnrich_Norm_GAPDH_Input)) +
  stat_summary(aes(fill = Disease_status), fun = mean, geom = "bar", colour = "black", width = 0.4) +
  geom_dotplot(aes(fill = Individual_ID), binaxis = "y", stackdir = "center", dotsize = 1.3) +
    scale_fill_manual(
    name = "Individual_ID",  
    values = c(
      "AF" = "#FDE6E0", 
      "SR" = "#F2FDFD",
      "BCVR-1638" = "#FFF2CC", 
      "BCVR-3913" = "#E4A3A1", 
      "BCVR-5888" = "#c8ee90", 
      "BCVR-5890" = "#56BED3", 
      "BCVR-9231" = "#ffcccb"
    ),
    breaks = c("AF", "SR", "BCVR-9231", "BCVR-1638", "BCVR-3913", "BCVR-5890", "BCVR-5888")
  ) +
  
  scale_y_continuous(name = "Normalized ChIP Enrichment") +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 14), 
    axis.title.y = element_text(size = 14),
    axis.text.y = element_text(size = 14) 
  ) 
ChIPqPCR_STRN_GAPDHInputnorm + geom_signif(
  comparisons = list(c("AF", "SR")), 
  map_signif_level = FALSE, 
  step_increase = 0.1,    
  textsize = 4  ,
  y_position = 1.5,
  annotations = "p = 0.2727")



### LRRC32

t_test_result <- t.test(data_analysis_dataframe_withdata_LRRC32_AFLA0328_AF$RelEnrich_Norm_GAPDH_Input, data_analysis_dataframe_withdata_LRRC32_AFLA0328_SR$RelEnrich_Norm_GAPDH_Input, alternative = "greater")
t_test_result
#pval  0.1041

ChIPqPCR_LRRC32_GAPDHInputnorm <- ggplot(data_analysis_dataframe_withdata_LRRC32_AFLA0328, 
                                         aes(x = factor(Disease_status, levels = c("AF", "SR")), 
                                             y = RelEnrich_Norm_GAPDH_Input)) +
  stat_summary(aes(fill = Disease_status), fun = mean, geom = "bar", colour = "black", width = 0.4) +
  geom_dotplot(aes(fill = Individual_ID), binaxis = "y", stackdir = "center", dotsize = 1.6) +
    scale_fill_manual(
    name = "Individual_ID",  
    values = c(
      "AF" = "#FDE6E0", 
      "SR" = "#F2FDFD",
      "BCVR-1638" = "#FFF2CC", 
      "BCVR-3913" = "#E4A3A1", 
      "BCVR-5888" = "#c8ee90", 
      "BCVR-5890" = "#56BED3", 
      "BCVR-9231" = "#ffcccb"
    ),
    breaks = c("AF", "SR", "BCVR-9231", "BCVR-1638", "BCVR-3913", "BCVR-5890", "BCVR-5888")
  ) +
  
  scale_y_continuous(name = "Normalized ChIP Enrichment") +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 14), 
    axis.title.y = element_text(size = 14),
    axis.text.y = element_text(size = 14) 
  ) 


ChIPqPCR_LRRC32_GAPDHInputnorm + geom_signif(
  comparisons = list(c("AF", "SR")), 
  map_signif_level = FALSE, 
  step_increase = 0.1,     
  textsize = 4  ,
  y_position = 9, 
  annotations = "p = 0.1041")


### SCX 


data_analysis_dataframe_withdata_SCX_B_AF <- data_analysis_dataframe_withdata_SCX_B[data_analysis_dataframe_withdata_SCX_B$Disease_status == "AF", ]
data_analysis_dataframe_withdata_SCX_B_SR <- data_analysis_dataframe_withdata_SCX_B[data_analysis_dataframe_withdata_SCX_B$Disease_status == "SR", ]
t_test_result <- t.test(data_analysis_dataframe_withdata_SCX_B_AF$RelEnrich_Norm_GAPDH_Input, data_analysis_dataframe_withdata_SCX_B_SR$RelEnrich_Norm_GAPDH_Input, alternative = "greater")


ChIPqPCR_SCX_B_GAPDHInputnorm <- ggplot(data_analysis_dataframe_withdata_SCX_B, 
                                        aes(x = factor(Disease_status, levels = c("AF", "SR")), 
                                            y = RelEnrich_Norm_GAPDH_Input)) +
  stat_summary(aes(fill = Disease_status), fun = mean, geom = "bar", colour = "black", width = 0.4) +
  geom_dotplot(aes(fill = Individual_ID), binaxis = "y", stackdir = "center", dotsize = 1.3) +
    scale_fill_manual(
    name = "Individual_ID",  
    values = c(
      "AF" = "#FDE6E0", 
      "SR" = "#F2FDFD",
      "BCVR-1638" = "#FFF2CC", 
      "BCVR-3913" = "#E4A3A1", 
      "BCVR-5888" = "#c8ee90", 
      "BCVR-5890" = "#56BED3", 
      "BCVR-9231" = "#ffcccb"
    ),
    breaks = c("AF", "SR", "BCVR-9231", "BCVR-1638", "BCVR-3913", "BCVR-5890", "BCVR-5888")
  ) +
  
  scale_y_continuous(name = "Normalized ChIP Enrichment") +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 14), 
    axis.title.y = element_text(size = 14),
    axis.text.y = element_text(size = 14) 
  ) 

ChIPqPCR_SCX_B_GAPDHInputnorm + geom_signif(
  comparisons = list(c("AF", "SR")), 
  map_signif_level = FALSE, 
  step_increase = 0.1,    
  textsize = 4  ,
  y_position = 13,
  annotations = "p = 0.1484")


### NPPB - Region AFLA008


data_analysis_dataframe_withdata_NPPB_AF008_AF <- data_analysis_dataframe_withdata_NPPB_AF008[data_analysis_dataframe_withdata_NPPB_AF008$Disease_status == "AF", ]
data_analysis_dataframe_withdata_NPPB_AF008_SR <- data_analysis_dataframe_withdata_NPPB_AF008[data_analysis_dataframe_withdata_NPPB_AF008$Disease_status == "SR", ]
t_test_result <- t.test(data_analysis_dataframe_withdata_NPPB_AF008_AF$RelEnrich_Norm_GAPDH_Input, data_analysis_dataframe_withdata_NPPB_AF008_SR$RelEnrich_Norm_GAPDH_Input, alternative = "greater")


ChIPqPCR_NPPB_AF008_GAPDHInputnorm <- ggplot(data_analysis_dataframe_withdata_NPPB_AF008,    # 📊 Figure 5c 
                                             aes(x = factor(Disease_status, levels = c("AF", "SR")), 
                                                 y = RelEnrich_Norm_GAPDH_Input)) +
  stat_summary(aes(fill = Disease_status), fun = mean, geom = "bar", colour = "black", width = 0.4) +
  geom_dotplot(aes(fill = Individual_ID), binaxis = "y", stackdir = "center", dotsize = 1.3) +
    scale_fill_manual(
    name = "Individual_ID",  
    values = c(
      "AF" = "#FDE6E0", 
      "SR" = "#F2FDFD",
      "BCVR-1638" = "#FFF2CC", 
      "BCVR-3913" = "#E4A3A1", 
      "BCVR-5888" = "#c8ee90", 
      "BCVR-5890" = "#56BED3", 
      "BCVR-9231" = "#ffcccb"
    ),
    breaks = c("AF", "SR", "BCVR-9231", "BCVR-1638", "BCVR-3913", "BCVR-5890", "BCVR-5888")
  ) +
  
  scale_y_continuous(name = "Normalized ChIP Enrichment") +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 14), 
    axis.title.y = element_text(size = 14),
    axis.text.y = element_text(size = 14) 
  ) 


ChIPqPCR_NPPB_AF008_GAPDHInputnorm + geom_signif(
  comparisons = list(c("AF", "SR")), 
  map_signif_level = FALSE, 
  step_increase = 0.1,     
  textsize = 4  ,
  y_position = 2.5, 
  annotations = "p = 0.1333")


### NPPB - Region AFLA0013


#Run ANOVA
data_analysis_dataframe_withdata_NPPB_AFLA0013

res_aov <- aov(RelEnrich_Norm_GAPDH_Input.1 ~ Disease_AnatSide,
               data = data_analysis_dataframe_withdata_NPPB_AFLA0013
)
summary(res_aov)

data_analysis_dataframe_withdata_NPPB_AFLA0013_PLOT_GAPDHInputnorm <- ggplot(data_analysis_dataframe_withdata_NPPB_AFLA0013,    # 📊 Figures 5c
                                                                             aes(x = factor(Disease_AnatSide,level=c("AF-LA","AF-RA","SR-LA","SR-RA")), 
                                                                                 y = RelEnrich_Norm_GAPDH_Input.1)) +
  stat_summary(aes(fill = Disease_AnatSide), fun = mean, geom = "bar", colour = "black", width = 0.4) +
  geom_dotplot(aes(fill = Individual_ID), binaxis = "y", stackdir = "center", dotsize = 1.8) +
  scale_fill_manual(
    name = "Individual_ID",  
    values = c(
      "AF-LA" = "#F9D4DA", 
      "AF-RA" = "#F7F1D0",
      "SR-LA" = "#D5E4F0",
      "SR-RA" = "#E3F3EC",
      "BCVR-1638" = "#FFF2CC", 
      "BCVR-3913" = "#E4A3A1", 
      "BCVR-5888" = "#c8ee90", 
      "BCVR-5890" = "#56BED3", 
      "BCVR-9231" = "#ffcccb"
    ),
    breaks = c("AF", "SR", "BCVR-9231", "BCVR-1638", "BCVR-3913", "BCVR-5890", "BCVR-5888")
  ) +
  
  scale_y_continuous(name = "Normalized ChIP Enrichment") +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 14), 
    axis.title.y = element_text(size = 14),
    axis.text.y = element_text(size = 14) 
  ) 
ChIPqPCR_SCX_B_GAPDHInputnorm + geom_signif(
  comparisons = list(c("AF", "SR")), 
  map_signif_level = FALSE, 
  step_increase = 0.1,     
  textsize = 4  ,
  y_position = 13, 
  annotations = "p = 0.1484")



#### Regions detected as SR-enriched in main cohort 


### MT1X

data_analysis_dataframe_withdata_MT1X_AF <- data_analysis_dataframe_withdata_NR2F1_B[data_analysis_dataframe_withdata_MT1X$Disease_status == "AF", ]
data_analysis_dataframe_withdata_MT1X_SR <- data_analysis_dataframe_withdata_NR2F1_B[data_analysis_dataframe_withdata_MT1X$Disease_status == "SR", ]


t_test_result <- t.test(data_analysis_dataframe_withdata_MT1X_11_AF$RelEnrich_Norm_GAPDH_Input, data_analysis_dataframe_withdata_MT1X_11_SR$RelEnrich_Norm_GAPDH_Input, alternative = "less")
t_test_result
#pval  0.02019

ChIPqPCR_MT1X_GAPDHInputnorm <- ggplot(data_analysis_dataframe_withdata_MT1X_11,       # 📊 Figures 5b
                                       aes(x = factor(Disease_status, levels = c("AF", "SR")), 
                                           y = RelEnrich_Norm_GAPDH_Input)) +
  stat_summary(aes(fill = Disease_status), fun = mean, geom = "bar", colour = "black", width = 0.4) +
  geom_dotplot(aes(fill = Individual_ID), binaxis = "y", stackdir = "center", dotsize = 1.3) +
    scale_fill_manual(
    name = "Individual_ID",  
    values = c(
      "AF" = "#FDE6E0", 
      "SR" = "#F2FDFD",
      "BCVR-1638" = "#FFF2CC", 
      "BCVR-3913" = "#E4A3A1", 
      "BCVR-5888" = "#c8ee90", 
      "BCVR-5890" = "#56BED3", 
      "BCVR-9231" = "#ffcccb"
    ),
    breaks = c("AF", "SR", "BCVR-9231", "BCVR-1638", "BCVR-3913", "BCVR-5890", "BCVR-5888")
  ) +
  
  scale_y_continuous(name = "Normalized ChIP Enrichment") +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 18), 
    axis.title.y = element_text(size = 18)
  ) + geom_signif(
    comparisons = list(c("AF", "SR")), 
    map_signif_level = FALSE, 
    step_increase = 0.1,     
    textsize = 4  ,
    y_position = 0.24, # Adjust text size
    annotations = "p = 0.02019"
  )



### NR2F1


data_analysis_dataframe_withdata_NR2F1_B_AF <- data_analysis_dataframe_withdata_NR2F1_B[data_analysis_dataframe_withdata_NR2F1_B$Disease_status == "AF", ]
data_analysis_dataframe_withdata_NR2F1_B_SR <- data_analysis_dataframe_withdata_NR2F1_B[data_analysis_dataframe_withdata_NR2F1_B$Disease_status == "SR", ]

t_test_result <- t.test(data_analysis_dataframe_withdata_NR2F1_B_AF$RelEnrich_Norm_GAPDH_Input, data_analysis_dataframe_withdata_NR2F1_B_SR$RelEnrich_Norm_GAPDH_Input, alternative = "less")
t_test_result
#pval  0.7479

ChIPqPCR_NR2F1_B_GAPDHInputnorm <- ggplot(data_analysis_dataframe_withdata_NR2F1_B, 
                                          aes(x = factor(Disease_status, levels = c("AF", "SR")), 
                                              y = RelEnrich_Norm_GAPDH_Input)) +
  stat_summary(aes(fill = Disease_status), fun = mean, geom = "bar", colour = "black", width = 0.4) +
  geom_dotplot(aes(fill = Individual_ID), binaxis = "y", stackdir = "center", dotsize = 1.3) +
    scale_fill_manual(
    name = "Individual_ID", 
    values = c(
      "AF" = "#FDE6E0", 
      "SR" = "#F2FDFD",
      "BCVR-1638" = "#FFF2CC", 
      "BCVR-3913" = "#E4A3A1", 
      "BCVR-5888" = "#c8ee90", 
      "BCVR-5890" = "#56BED3", 
      "BCVR-9231" = "#ffcccb"
    ),
    breaks = c("AF", "SR", "BCVR-9231", "BCVR-1638", "BCVR-3913", "BCVR-5890", "BCVR-5888")
  ) +
  
  scale_y_continuous(name = "Normalized ChIP Enrichment") +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 14), 
    axis.title.y = element_text(size = 14),
    axis.text.y = element_text(size = 14) 
  ) 


ChIPqPCR_NR2F1_B_GAPDHInputnorm + geom_signif(
  comparisons = list(c("AF", "SR")), 
  map_signif_level = FALSE, 
  step_increase = 0.1,    
  textsize = 4  ,
  y_position = 1.37, 
  annotations = "p = 0.7479")


### BLM

data_analysis_dataframe_withdata_BLM_SRLA617_AF <- data_analysis_dataframe_withdata_BLM_SRLA617[data_analysis_dataframe_withdata_BLM_SRLA617$Disease_status == "AF", ]
data_analysis_dataframe_withdata_BLM_SRLA617_SR <- data_analysis_dataframe_withdata_BLM_SRLA617[data_analysis_dataframe_withdata_BLM_SRLA617$Disease_status == "SR", ]

t_test_result <- t.test(data_analysis_dataframe_withdata_BLM_SRLA617_AF$RelEnrich_Norm_GAPDH_Input, data_analysis_dataframe_withdata_BLM_SRLA617_SR$RelEnrich_Norm_GAPDH_Input, alternative = "less")
t_test_result
#pval  0.1189

ChIPqPCR_BLM_SRLA617_GAPDHInputnorm <- ggplot(data_analysis_dataframe_withdata_BLM_SRLA617, 
                                              aes(x = factor(Disease_status, levels = c("AF", "SR")), 
                                                  y = RelEnrich_Norm_GAPDH_Input)) +
  stat_summary(aes(fill = Disease_status), fun = mean, geom = "bar", colour = "black", width = 0.4) +
  geom_dotplot(aes(fill = Individual_ID), binaxis = "y", stackdir = "center", dotsize = 1.3) +
    scale_fill_manual(
    name = "Individual_ID", 
    values = c(
      "AF" = "#FDE6E0", 
      "SR" = "#F2FDFD",
      "BCVR-1638" = "#FFF2CC", 
      "BCVR-3913" = "#E4A3A1", 
      "BCVR-5888" = "#c8ee90", 
      "BCVR-5890" = "#56BED3", 
      "BCVR-9231" = "#ffcccb"
    ),
    breaks = c("AF", "SR", "BCVR-9231", "BCVR-1638", "BCVR-3913", "BCVR-5890", "BCVR-5888")
  ) +
  
  scale_y_continuous(name = "Normalized ChIP Enrichment") +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 14), 
    axis.title.y = element_text(size = 14),
    axis.text.y = element_text(size = 14) 
  ) 
ChIPqPCR_BLM_SRLA617_GAPDHInputnorm + geom_signif(
  comparisons = list(c("AF", "SR")), 
  map_signif_level = FALSE, 
  step_increase = 0.1,    
  textsize = 4  ,
  y_position = 0.27, 
  annotations = "p = 0.1189")



### RGS6

data_analysis_dataframe_withdata_RGS6_443_AF <- data_analysis_dataframe_withdata_RGS6_443[data_analysis_dataframe_withdata_RGS6_443$Disease_status == "AF", ]
data_analysis_dataframe_withdata_RGS6_443_SR <- data_analysis_dataframe_withdata_RGS6_443[data_analysis_dataframe_withdata_RGS6_443$Disease_status == "SR", ]

t_test_result <- t.test(data_analysis_dataframe_withdata_RGS6_443_AF$RelEnrich_Norm_GAPDH_Input, data_analysis_dataframe_withdata_RGS6_443_SR$RelEnrich_Norm_GAPDH_Input, alternative = "less")
t_test_result
#pval  0.9441
ChIPqPCR_RGS6_443_GAPDHInputnorm <- ggplot(data_analysis_dataframe_withdata_RGS6_443, 
                                           aes(x = factor(Disease_status, levels = c("AF", "SR")), 
                                               y = RelEnrich_Norm_GAPDH_Input)) +
  stat_summary(aes(fill = Disease_status), fun = mean, geom = "bar", colour = "black", width = 0.4) +
  geom_dotplot(aes(fill = Individual_ID), binaxis = "y", stackdir = "center", dotsize = 1.3) +
    scale_fill_manual(
    name = "Individual_ID", 
    values = c(
      "AF" = "#FDE6E0", 
      "SR" = "#F2FDFD",
      "BCVR-1638" = "#FFF2CC", 
      "BCVR-3913" = "#E4A3A1", 
      "BCVR-5888" = "#c8ee90", 
      "BCVR-5890" = "#56BED3", 
      "BCVR-9231" = "#ffcccb"
    ),
    breaks = c("AF", "SR", "BCVR-9231", "BCVR-1638", "BCVR-3913", "BCVR-5890", "BCVR-5888")
  ) +
  
  scale_y_continuous(name = "Normalized ChIP Enrichment") +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 14), 
    axis.title.y = element_text(size = 14),
    axis.text.y = element_text(size = 14) 
  ) 
ChIPqPCR_RGS6_443_GAPDHInputnorm + geom_signif(
  comparisons = list(c("AF", "SR")), 
  map_signif_level = FALSE, 
  step_increase = 0.1,    
  textsize = 4  ,
  y_position = 19,
  annotations = "p = 0.9441")



### ANGPT1

data_analysis_dataframe_withdata_ANGPT1_SRRA306_AF <- data_analysis_dataframe_withdata_ANGPT1_SRRA306[data_analysis_dataframe_withdata_ANGPT1_SRRA306$Disease_status == "AF", ]
data_analysis_dataframe_withdata_ANGPT1_SRRA306_SR <- data_analysis_dataframe_withdata_ANGPT1_SRRA306[data_analysis_dataframe_withdata_ANGPT1_SRRA306$Disease_status == "SR", ]

t_test_result <- t.test(data_analysis_dataframe_withdata_ANGPT1_SRRA306_AF$RelEnrich_Norm_GAPDH_Input, data_analysis_dataframe_withdata_ANGPT1_SRRA306_SR$RelEnrich_Norm_GAPDH_Input, alternative = "less")
t_test_result
#pval  0.5824

ChIPqPCR_ANGPT1_SRRA306_GAPDHInputnorm <- ggplot(data_analysis_dataframe_withdata_ANGPT1_SRRA306, 
                                                 aes(x = factor(Disease_status, levels = c("AF", "SR")), 
                                                     y = RelEnrich_Norm_GAPDH_Input)) +
  stat_summary(aes(fill = Disease_status), fun = mean, geom = "bar", colour = "black", width = 0.4) +
  geom_dotplot(aes(fill = Individual_ID), binaxis = "y", stackdir = "center", dotsize = 1.3) +
    scale_fill_manual(
    name = "Individual_ID",  # temporary name for legend (gets overridden below)
    values = c(
      "AF" = "#FDE6E0", 
      "SR" = "#F2FDFD",
      "BCVR-1638" = "#FFF2CC", 
      "BCVR-3913" = "#E4A3A1", 
      "BCVR-5888" = "#c8ee90", 
      "BCVR-5890" = "#56BED3", 
      "BCVR-9231" = "#ffcccb"
    ),
    breaks = c("AF", "SR", "BCVR-9231", "BCVR-1638", "BCVR-3913", "BCVR-5890", "BCVR-5888")
  ) +
  
  scale_y_continuous(name = "Normalized ChIP Enrichment") +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 14), 
    axis.title.y = element_text(size = 14),
    axis.text.y = element_text(size = 14) 
  ) 
ChIPqPCR_ANGPT1_SRRA306_GAPDHInputnorm + geom_signif(
  comparisons = list(c("AF", "SR")), 
  map_signif_level = FALSE, 
  step_increase = 0.1,    
  textsize = 4  ,
  y_position = 0.235, 
  annotations = "p = 0.5824")







