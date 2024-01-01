library(ggplot2)
library(gridExtra)

#specify path for where DESeq2 outputs files are stored for each treatment
#DESEq2 outputs correspond to differential expression between the sexes within a mating treatment
path_diffExp_trt <- "C:/Users/mishr/OneDrive/Desktop/ProjectII/Annotated_R_scripts_and_data/Diff.Expression.Between.Sexes.Within.Mating.Treatment/"

#specify path for where the SBGE data from Osada et al. is stored
path_for_SBGE_data <- "C:/Users/mishr/OneDrive/Desktop/ProjectII/Annotated_R_scripts_and_data/Sex.Biased.Gene.Expression.Osada_et_al_2017/"

###########################################################################
## Plots of log2FC in Trt1 vs Trt2 ########################################
###########################################################################
#define function to make scatterplots of log2FC male/female for treatment_1 vs treatment_2
SexBiasPlotsBetweenTreatments <- function(trt1, trt2, tissue){
  #load DESeq2 output files for each treatment
  #make a dataframe of log2FC for each treatment
  trt1.file <- read.table(paste0(path_diffExp_trt, trt1, ".", tissue, ".tsv"), header = T, sep = '\t')
  trt2.file <- read.table(paste0(path_diffExp_trt, trt2, ".", tissue, ".tsv"), header = T, sep = '\t')
  df <- cbind.data.frame(Sex.Bias.Trt1 = trt1.file$log2FoldChange, Sex.Bias.Trt2 = trt2.file$log2FoldChange)
  
  #make plot with LOESS fit
  p <- ggplot(df, aes(x = Sex.Bias.Trt1, y = Sex.Bias.Trt2)) + 
    theme_minimal() +
    geom_point(size = 0.5, colour = "grey78") + 
    geom_smooth(method = "loess", lwd = 1, se = T, level = 0.95, span = 0.5, colour = "blue", fill = "green") +
    geom_abline(slope = 1, intercept = 0, colour = "brown", lwd = 1) +
    geom_hline(yintercept = 0) +
    geom_vline(xintercept = 0) +
    labs(x = paste0("Sex Bias in ", trt1), y = paste0("Sex Bias in ", trt2)) +
    theme(legend.title = element_blank(),
          plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(size = 12, hjust = 0.5), 
          axis.text.y = element_text(size = 12, hjust = 1),
          axis.title.x = element_text(margin = margin(0,10,0,0), size = 14),
          axis.title.y = element_text(margin = margin(0,10,0,0), size = 14)) 
  
  return(p) 
}

#use above function to make the scatterplots, arrange into a grid
b1 <- SexBiasPlotsBetweenTreatments("MCabs", "MCsim", "Body")
b2 <- SexBiasPlotsBetweenTreatments("MCabs", "MCcom", "Body")
b3 <- SexBiasPlotsBetweenTreatments("MCsim", "MCcom", "Body")

h1 <- SexBiasPlotsBetweenTreatments("MCabs", "MCsim", "Head")
h2 <- SexBiasPlotsBetweenTreatments("MCabs", "MCcom", "Head")
h3 <- SexBiasPlotsBetweenTreatments("MCsim", "MCcom", "Head")

grid.arrange(b1, b2, b3, h1, h2, h3, nrow = 2)



###########################################################################
## Plots of log2FC from Osada data vs each of our mating treatment ########
###########################################################################

SexBiasPlotsVsOsada <- function(treatment, tissue){
  #load DESeq2 output file for mating treatment and sex bias data from Osada et al
  sexbias <- read.table(paste0(path_for_SBGE_data, "Sex.Biased.Genes.", tissue, ".tsv"), header = T, sep = '\t')
  trt.file <- read.table(paste0(path_diffExp_trt, treatment, ".", tissue, ".tsv"), header = T, sep = '\t')
  
  #find common genes and subset above dataframes
  gene.list <- intersect(sexbias$geneID, trt.file$FlyBaseID)
  sexbias <- sexbias[sexbias$geneID %in% gene.list,]
  trt.file <- trt.file[trt.file$FlyBaseID %in% gene.list,]
  
  #make a dataframe of log2FC 
  df <- cbind.data.frame(Sex.Bias.Osada = sexbias$log2FoldChange, Sex.Bias.Trt = trt.file$log2FoldChange)
  
  #make plot with LOESS fit
  p <- ggplot(df, aes(x = Sex.Bias.Osada, y = Sex.Bias.Trt)) + 
    theme_minimal() +
    geom_point(size = 0.5, colour = "grey78") + 
    geom_smooth(method = "loess", lwd = 1, se = T, level = 0.95, span = 0.5, colour = "blue", fill = "green") +
    geom_abline(slope = 1, intercept = 0, colour = "brown", lwd = 1) +
    geom_hline(yintercept = 0) +
    geom_vline(xintercept = 0) +
    labs(x = paste0("Sex Bias in Osada et al.'s data"), y = paste0("Sex Bias in ", treatment)) +
    theme(legend.title = element_blank(),
          plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(size = 12, hjust = 0.5), 
          axis.text.y = element_text(size = 12, hjust = 1),
          axis.title.x = element_text(margin = margin(0,10,0,0), size = 14),
          axis.title.y = element_text(margin = margin(0,10,0,0), size = 14)) 
  
  return(p)  
}

body1 <- SexBiasPlotsVsOsada("MCabs", "Body")
body2 <- SexBiasPlotsVsOsada("MCsim", "Body")
body3 <- SexBiasPlotsVsOsada("MCcom", "Body")

head1 <- SexBiasPlotsVsOsada("MCabs", "Head")
head2 <- SexBiasPlotsVsOsada("MCsim", "Head")
head3 <- SexBiasPlotsVsOsada("MCcom", "Head")

grid.arrange(body1, body2, body3, head1, head2, head3, nrow = 2)


###########################################################################
## Correlation of log2FC in Trt1 vs Trt2 ##################################
###########################################################################
SexBiasCorrelationBetweenTreatments <- function(trt1, trt2, tissue){
  #load DESeq2 output files for each treatment
  #make a dataframe of log2FC for each treatment
  trt1.file <- read.table(paste0(path_diffExp_trt, trt1, ".", tissue, ".tsv"), header = T, sep = '\t')
  trt2.file <- read.table(paste0(path_diffExp_trt, trt2, ".", tissue, ".tsv"), header = T, sep = '\t')
  df <- cbind.data.frame(Sex.Bias.Trt1 = trt1.file$log2FoldChange, Sex.Bias.Trt2 = trt2.file$log2FoldChange)
  
  #estimate correlation between sex bias
  cor.log2FC <- cor.test(df$Sex.Bias.Trt1, df$Sex.Bias.Trt2)
  
  return(cor.log2FC) 
}

SexBiasCorrelationBetweenTreatments("MCabs", "MCsim", "Body")
SexBiasCorrelationBetweenTreatments("MCabs", "MCcom", "Body")
SexBiasCorrelationBetweenTreatments("MCsim", "MCcom", "Body")

SexBiasCorrelationBetweenTreatments("MCabs", "MCsim", "Head")
SexBiasCorrelationBetweenTreatments("MCabs", "MCcom", "Head")
SexBiasCorrelationBetweenTreatments("MCsim", "MCcom", "Head")


###########################################################################
## Correlation of log2FC from Osada data vs each of our mating treatment ##
###########################################################################
SexBiasCorrelationVsOsada <- function(treatment, tissue){
  #load DESeq2 output file for mating treatment and sex bias data from Osada et al
  sexbias <- read.table(paste0(path_for_SBGE_data, "Sex.Biased.Genes.", tissue, ".tsv"), header = T, sep = '\t')
  trt.file <- read.table(paste0(path_diffExp_trt, treatment, ".", tissue, ".tsv"), header = T, sep = '\t')
  
  #find common genes and subset above dataframes
  gene.list <- intersect(sexbias$geneID, trt.file$FlyBaseID)
  sexbias <- sexbias[sexbias$geneID %in% gene.list,]
  trt.file <- trt.file[trt.file$FlyBaseID %in% gene.list,]
  
  #make a dataframe of log2FC 
  df <- cbind.data.frame(Sex.Bias.Osada = sexbias$log2FoldChange, Sex.Bias.Trt = trt.file$log2FoldChange)
  
  #estimate correlation between sex bias
  cor.log2FC <- cor.test(df$Sex.Bias.Osada, df$Sex.Bias.Trt)
  
  return(cor.log2FC)
}

SexBiasCorrelationVsOsada("MCabs", "Body")
SexBiasCorrelationVsOsada("MCsim", "Body")
SexBiasCorrelationVsOsada("MCcom", "Body")

SexBiasCorrelationVsOsada("MCabs", "Head")
SexBiasCorrelationVsOsada("MCsim", "Head")
SexBiasCorrelationVsOsada("MCcom", "Head")
