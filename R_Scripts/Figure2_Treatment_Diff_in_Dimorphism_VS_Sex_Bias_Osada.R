library(ggplot2)
library(gridExtra)

#set path to load gonad specificity data
#this data may be used to make versions of figures without genes with highly gonad-specific expression
path_1 <- 'C:/Users/mishr/OneDrive/Desktop/ProjectII/Annotated_R_scripts_and_data/Gonad.Specific.Expression/'
GSI.females <- read.table(paste0(path_1, "Gonad.Specificity.Index.Females.tsv"), header = T)
GSI.males <- read.table(paste0(path_1, "Gonad.Specificity.Index.Males.tsv"), header = T)

gonad.specific.genes <- union(GSI.females[which(GSI.females$Gonad.Specificity.Index >=0.95),]$geneID, 
                              GSI.males[which(GSI.males$Gonad.Specificity.Index >=0.95),]$geneID)

#specify path for where DESeq2 outputs files are stored for each treatment
#DESEq2 outputs correspond to differential expression between the sexes within a mating treatment
path_diffExp_trt <- "C:/Users/mishr/OneDrive/Desktop/ProjectII/Annotated_R_scripts_and_data/Diff.Expression.Between.Sexes.Within.Mating.Treatment/"

#specify path for where the SBGE data from Osada et al. is stored
path_for_SBGE_data <- "C:/Users/mishr/OneDrive/Desktop/ProjectII/Annotated_R_scripts_and_data/Sex.Biased.Gene.Expression.Osada_et_al_2017/"


##########################################################################################################
### function to make Fig.2 : scatterplots of Treatment Difference in Dimorphism vs               #########
###  Sex Bias from Osada et al with a LOESS regression applied for fitted line                   #########
##########################################################################################################

#comment out relevant lines to include gonad-specific genes

#NOTE 1: in this version, the Y-axis is (sign of sex bias from Osada data)*(Treatment Difference in Dimorphism), in
#order to make plot interpretation easier
#NOTE 2: the y-axis is also limited from -1 to 1, to make visualisation easier

LoessRegressionPlots_MainVersion <- function(treatment1, treatment2, tissue){
  
  #load SBGE data obtained from external dataset (Osada et al 2017) 
  sexbias <- read.table(paste0(path_for_SBGE_data, "Sex.Biased.Genes.", tissue, ".tsv"), header = T, sep = '\t')
  
  #load DESeq2 outputs for differential expression between the sexes for a given mating treatment
  trt1.file <- read.table(paste0(path_diffExp_trt, treatment1, ".", tissue, ".tsv"), header = T, sep = '\t')
  trt2.file <- read.table(paste0(path_diffExp_trt, treatment2, ".", tissue, ".tsv"), header = T, sep = '\t')
  
  #only retain genes that are common to all of the above three files
  gene.list <- intersect(sexbias$geneID, intersect(trt1.file$FlyBaseID, trt2.file$FlyBaseID))
  sexbias <- sexbias[sexbias$geneID %in% gene.list,]
  trt1.file <- trt1.file[trt1.file$FlyBaseID %in% gene.list,]
  trt2.file <- trt2.file[trt2.file$FlyBaseID %in% gene.list,]
  
  #for loop to assign a sign (-1 or 1) to sex bias values from data in Osada et al 2017
  sex.bias.sign <- c()
  sign <- NULL
  for (i in 1:nrow(sexbias)){
    if(sexbias$log2FoldChange[i] > 0){
      sign = 1
    } else {
      sign = -1
    }
    sex.bias.sign <- c(sex.bias.sign, sign)
  }
  
  
  #make dataframe of sex bias and difference in sex effect between treatment pair (transformed with sign of sex bias)
  trt.diff.in.dimorphism <- sex.bias.sign*(trt1.file$log2FoldChange - trt2.file$log2FoldChange)
  df <- cbind.data.frame(geneID = gene.list, Trt.Diff.in.Dimorphism = trt.diff.in.dimorphism, Sex.Bias.Osada = sexbias$log2FoldChange)
  
  #COMMENT THIS OUT TO RUN WITH GONAD-SPECIFIC GENES
  #df <- df[!(df$geneID %in% gonad.specific.genes),]
  
  #make plot
  p <- ggplot(df, aes(x = Sex.Bias.Osada, y = Trt.Diff.in.Dimorphism)) + 
    theme_classic() +
    geom_point(size = 0.5, colour = "grey78") + 
    geom_smooth(method = "loess", lwd = 0.7, se = T, level = 0.95, span = 0.5, colour = "blue", fill = "green") +
    geom_smooth(method = "loess", lwd = 0.7, se = T, level = 0.95, span = 0.5, colour = "blue", fill = "green") +
    geom_hline(yintercept = 0, colour = "black", lwd = 0.5, linetype = "solid") +
    #geom_vline(xintercept = -5, colour = "black", lwd = 0.5, linetype = "dashed") +
    geom_vline(xintercept = -0.5, colour = "black", lwd = 0.5, linetype = "dashed") +
    geom_vline(xintercept = 0.5, colour = "black", lwd = 0.5, linetype = "dashed") +
    #geom_vline(xintercept = 5, colour = "black", lwd = 0.5, linetype = "dashed") +
    labs(x = "Sex Bias", y = "Difference in Sex Effect", title = paste0(treatment1, " - ", treatment2, ", ", tissue)) +
    theme(legend.title = element_blank(),
          plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(size = 12, hjust = 0.5), 
          axis.text.y = element_text(size = 12, hjust = 1),
          axis.title.x = element_text(margin = margin(0,10,0,0), size = 14),
          axis.title.y = element_text(margin = margin(0,10,0,0), size = 14)) + coord_cartesian(ylim = c(-1,1))
  
  return(p)  
}

b1 <- LoessRegressionPlots_MainVersion("MCabs", "MCsim", "Body")
b2 <- LoessRegressionPlots_MainVersion("MCabs", "MCcom", "Body")
b3 <- LoessRegressionPlots_MainVersion("MCsim", "MCcom", "Body")

h1 <- LoessRegressionPlots_MainVersion("MCabs", "MCsim", "Head")
h2 <- LoessRegressionPlots_MainVersion("MCabs", "MCcom", "Head")
h3 <- LoessRegressionPlots_MainVersion("MCsim", "MCcom", "Head")

grid.arrange(b1, b2, b3, nrow = 3)
grid.arrange(b1, b2, b3, h1, h2, h3, nrow = 2)



##########################################################################################################
### function to make an alternate version Fig.2 : without transformation applied to Y-axis       #########
##########################################################################################################

#comment out relevant lines to include gonad-specific genes

#NOTE 1: in this alternate version, the Y-axis is (Treatment Difference in Dimorphism), without transforming it
#around the axis by multiplying by sign of sex bias from Osada et al [this version was NOT used in final draft]
#NOTE 2: the y-axis is again limited from -1 to 1, to make visualisation easier

LoessRegressionPlot <- function(treatment1, treatment2, tissue){
  #load SBGE data obtained from external dataset (Osada et al 2017) 
  sexbias <- read.table(paste0(path_for_SBGE_data, "Sex.Biased.Genes.", tissue, ".tsv"), header = T, sep = '\t')
  
  #load DESeq2 outputs for differential expression between the sexes for a given mating treatment
  trt1.file <- read.table(paste0(path_diffExp_trt, treatment1, ".", tissue, ".tsv"), header = T, sep = '\t')
  trt2.file <- read.table(paste0(path_diffExp_trt, treatment2, ".", tissue, ".tsv"), header = T, sep = '\t')
  
  #only retain genes that are common to all of the above three files
  gene.list <- intersect(sexbias$geneID, intersect(trt1.file$FlyBaseID, trt2.file$FlyBaseID))
  sexbias <- sexbias[sexbias$geneID %in% gene.list,]
  trt1.file <- trt1.file[trt1.file$FlyBaseID %in% gene.list,]
  trt2.file <- trt2.file[trt2.file$FlyBaseID %in% gene.list,]
  
  #make dataframe of sex bias and difference in sex effect between treatment pair
  trt.diff.in.dimorphism <- trt1.file$log2FoldChange - trt2.file$log2FoldChange
  df <- cbind.data.frame(geneID = gene.list, Trt.Diff.in.Dimorphism = trt.diff.in.dimorphism, Sex.Bias.Osada = sexbias$log2FoldChange)
  
  #COMMENT THIS OUT TO RUN WITH GONAD-SPECIFIC GENES
  #df <- df[!(df$geneID %in% gonad.specific.genes),]
  
  #df <- df[!(df$geneID %in% gonad.specific.genes),]
  
  #make plot
  p <- ggplot(df, aes(x = Sex.Bias.Osada, y = Trt.Diff.in.Dimorphism)) + 
    theme_classic() +
    geom_point(size = 0.5, colour = "grey78") + 
    geom_smooth(method = "loess", lwd = 0.7, se = T, level = 0.95, span = 0.5, colour = "blue", fill = "green") +
    geom_hline(yintercept = 0, colour = "black", lwd = 0.5, linetype = "solid") +
    geom_vline(xintercept = -5, colour = "black", lwd = 0.5, linetype = "dashed") +
    geom_vline(xintercept = -0.5, colour = "black", lwd = 0.5, linetype = "dashed") +
    geom_vline(xintercept = 0.5, colour = "black", lwd = 0.5, linetype = "dashed") +
    geom_vline(xintercept = 5, colour = "black", lwd = 0.5, linetype = "dashed") +
    labs(x = "Sex Bias", y = "Difference in Sex Effect", title = paste0(treatment1, " - ", treatment2, ", ", tissue)) +
    theme(legend.title = element_blank(),
          plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(size = 12, hjust = 0.5), 
          axis.text.y = element_text(size = 12, hjust = 1),
          axis.title.x = element_text(margin = margin(0,10,0,0), size = 14),
          axis.title.y = element_text(margin = margin(0,10,0,0), size = 14)) + coord_cartesian(ylim = c(-1,1)) 
  #annotate(geom="text", x=-1.5, y=0.9, label="FB", color="black", fontface="bold") +
  #annotate(geom="text", x=1.5, y=0.9, label="MB", color="black", fontface="bold") 
  #annotate(geom="text", x=-0.4, y=1, label="UB", color="black")
  
  return(p)  
}

b1 <- LoessRegressionPlot("MCabs", "MCsim", "Body")
b2 <- LoessRegressionPlot("MCabs", "MCcom", "Body")
b3 <- LoessRegressionPlot("MCsim", "MCcom", "Body")

h1 <- LoessRegressionPlot("MCabs", "MCsim", "Head")
h2 <- LoessRegressionPlot("MCabs", "MCcom", "Head")
h3 <- LoessRegressionPlot("MCsim", "MCcom", "Head")

grid.arrange(b1, b2, b3, h1, h2, h3, nrow = 2)


