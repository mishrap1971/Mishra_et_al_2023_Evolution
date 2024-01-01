#####################################################
###  for sex-biased gene expression #################
#####################################################

#specify path for where DESeq2 outputs files are stored for each treatment
#DESEq2 outputs correspond to differential expression between the sexes within a mating treatment
path_diffExp_trt <- "C:/Users/mishr/OneDrive/Desktop/ProjectII/Annotated_R_scripts_and_data/Diff.Expression.Between.Sexes.Within.Mating.Treatment/"

#define function to find % of genes with significant sex-biased gene expression (p-adj < 0.05) in 
#all treatments within a tissue (body/head)
#function also performs Fisher's exact test of proportion of sig. sex-biased genes in the 3 treatments
#NOTE: we exclude "unbiased" genes from analysis (i.e., |log2FC| < 0.5)

Proportion.of.genes.with.sig.SexBias <- function(tissue){
  MCabs.df <- read.table(paste0(path_diffExp_trt, 'MCabs.', tissue, '.tsv'), sep = '\t', header = T)
  MCsim.df <- read.table(paste0(path_diffExp_trt, 'MCsim.', tissue, '.tsv'), sep = '\t', header = T)
  MCcom.df <- read.table(paste0(path_diffExp_trt, 'MCcom.', tissue, '.tsv'), sep = '\t', header = T)
  
  sig.genes.MCabs <- MCabs.df[which(MCabs.df$padj < 0.05 & abs(MCabs.df$log2FoldChange) > 0.5),]
  sig.genes.MCsim <- MCsim.df[which(MCsim.df$padj < 0.05 & abs(MCsim.df$log2FoldChange) > 0.5),]
  sig.genes.MCcom <- MCcom.df[which(MCcom.df$padj < 0.05 & abs(MCcom.df$log2FoldChange) > 0.5),]
  
  nonsig.genes.MCabs <- MCabs.df[which(MCabs.df$padj >= 0.05 & abs(MCabs.df$log2FoldChange) > 0.5),]
  nonsig.genes.MCsim <- MCsim.df[which(MCsim.df$padj >= 0.05 & abs(MCsim.df$log2FoldChange) > 0.5),]
  nonsig.genes.MCcom <- MCcom.df[which(MCcom.df$padj >= 0.05 & abs(MCcom.df$log2FoldChange) > 0.5),]
  
  #fisher's exact test
  sig.genes <- c(nrow(sig.genes.MCabs), nrow(sig.genes.MCsim), nrow(sig.genes.MCcom))
  nonsig.genes <- c(nrow(nonsig.genes.MCabs), nrow(nonsig.genes.MCsim), nrow(nonsig.genes.MCcom))
  df1 <- rbind.data.frame(sig.genes, nonsig.genes)
  colnames(df1) <- c("MCabs", "MCsim", "MCcom") 
  row.names(df1) <- c("sig.genes", "nonsig.genes")
  test <- fisher.test(df1, workspace=2e9)
  
  return(df1)
  
  #uncomment the below line if fisher's test p-values are required
  #return(test$p.value)
  
}

Proportion.of.genes.with.sig.SexBias("Head")
Proportion.of.genes.with.sig.SexBias("Body")


################################################
###  for sex-specific splicing #################
################################################
#Set global variables
FDRThreshold = 0.1
mappedReadsThreshold = 50

#set path to folder where JunctionSeq outputs are stored (for Diff. Splicing between Treatments within one sex) 
path_for_JS_data <- "C:/Users/mishr/OneDrive/Desktop/ProjectII/Annotated_R_scripts_and_data/JunctionSeq_Differential_Splicing_Outputs/Between.Sex.Within.A.Treatment/"

#define function to find % of genes with significant sex-specifc splicing (p-adj < 0.1) in 
#all treatments within a tissue (body/head)
#function also performs Fisher's exact test of proportion of sig. sex-specifically spliced genes in the 3 treatments
Proportion.of.genes.with.sig.SexSpecificSplicing <- function(tissue){
  #load JunctionSeq outputs for each treatment
  MCabs.JS.df <- read.table(paste0(path_for_JS_data, 'MCabs.', tissue, '/allGenes.results.txt'), sep = "\t", header = TRUE)
  MCsim.JS.df <- read.table(paste0(path_for_JS_data, 'MCsim.', tissue, '/allGenes.results.txt'), sep = "\t", header = TRUE)
  MCcom.JS.df <- read.table(paste0(path_for_JS_data, 'MCcom.', tissue, '/allGenes.results.txt'), sep = "\t", header = TRUE)
  
  #remove NAs from dataframes
  MCabs.JS.df <- MCabs.JS.df[!(is.na(MCabs.JS.df[,23])) | ! (is.na(MCabs.JS.df[,24])),]
  MCsim.JS.df <- MCsim.JS.df[!(is.na(MCsim.JS.df[,23])) | ! (is.na(MCsim.JS.df[,24])),]
  MCcom.JS.df <- MCcom.JS.df[!(is.na(MCcom.JS.df[,23])) | ! (is.na(MCcom.JS.df[,24])),]
  
  #remove all genes with fewer than 50 reads mapping to exons in males and females
  MCabs.JS.df <- MCabs.JS.df[MCabs.JS.df[,23] > mappedReadsThreshold & MCabs.JS.df[,24] > mappedReadsThreshold,]
  MCsim.JS.df <- MCsim.JS.df[MCsim.JS.df[,23] > mappedReadsThreshold & MCsim.JS.df[,24] > mappedReadsThreshold,]
  MCcom.JS.df <- MCcom.JS.df[MCcom.JS.df[,23] > mappedReadsThreshold & MCcom.JS.df[,24] > mappedReadsThreshold,]
  
  #obtain list of common genes and retain them in the 3 treatments
  gene.list <- intersect(MCabs.JS.df$geneID, intersect(MCsim.JS.df$geneID, MCcom.JS.df$geneID))
  MCabs.JS.df <- MCabs.JS.df[MCabs.JS.df$geneID %in% gene.list,]
  MCsim.JS.df <- MCsim.JS.df[MCsim.JS.df$geneID %in% gene.list,]
  MCcom.JS.df <- MCcom.JS.df[MCcom.JS.df$geneID %in% gene.list,]
  
  #genes with significant differential splicing in each treatment
  #MCabs
  sig.genes.MCabs <- unique(MCabs.JS.df$geneID[MCabs.JS.df$geneWisePadj < FDRThreshold])
  nonsig.genes.MCabs <- unique(MCabs.JS.df$geneID[MCabs.JS.df$geneWisePadj >= FDRThreshold])
  
  #MCsim
  sig.genes.MCsim <- unique(MCsim.JS.df$geneID[MCsim.JS.df$geneWisePadj < FDRThreshold])
  nonsig.genes.MCsim <- unique(MCsim.JS.df$geneID[MCsim.JS.df$geneWisePadj >= FDRThreshold])
  
  #MCcom
  sig.genes.MCcom <- unique(MCcom.JS.df$geneID[MCcom.JS.df$geneWisePadj < FDRThreshold])
  nonsig.genes.MCcom <- unique(MCcom.JS.df$geneID[MCcom.JS.df$geneWisePadj >= FDRThreshold])
  
  #fisher's exact test
  sig.genes <- c(length(sig.genes.MCabs), length(sig.genes.MCsim), length(sig.genes.MCcom))
  nonsig.genes <- c(length(nonsig.genes.MCabs), length(nonsig.genes.MCsim), length(nonsig.genes.MCcom))
  df2 <- rbind.data.frame(sig.genes, nonsig.genes)
  colnames(df2) <- c("MCabs", "MCsim", "MCcom") 
  row.names(df2) <- c("sig.genes", "nonsig.genes")
  test <- fisher.test(df2, workspace=2e9)
  
  #return (test$p.value)
  return(df2)
}

Proportion.of.genes.with.sig.SexSpecificSplicing("Head")
Proportion.of.genes.with.sig.SexSpecificSplicing("Body")