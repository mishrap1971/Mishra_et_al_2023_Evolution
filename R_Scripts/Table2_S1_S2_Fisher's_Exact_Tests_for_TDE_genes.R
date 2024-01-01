#This script takes the Treatment Effect dataframes to perform Fisher's exacts tests, outputs from
#which were used to make Tables 2, S1 and S2

#TDE genes = Treatment Differentially Expressed Genes

######################################################################################################
### Table S1: Proportion of TDE genes for each sex in each pairwise comparison of treatments #########
######################################################################################################

#first define path for DESeq2 outputs for differential expression between treatments
#we use the raw DESeq2 output and not the "Treatment Effect Dataframes" because the latter 
#excludes genes for which we lacked SBGE (sex-biased gene expression) data
path_1 <- 'C:/Users/mishr/OneDrive/Desktop/ProjectII/Annotated_R_scripts_and_data/Diff.Expression.Between.Mating.Treatments/'

#define function to yield %TDE genes in each pairwise comparison of treatments
#e.g., sex = "Female", tissue = "Male", treatment.comparison = "MCsim.vs.MCcom"
Percentage.of.TDE.genes <- function(sex, tissue, treatment.comparison){
  #load treatment effect dataframe
  trtEffect.df <- na.omit(read.table(paste0(path_1, sex, ".", tissue, ".", treatment.comparison ,".tsv"), header = T))
  
  #total number of genes
  total.genes <- nrow(trtEffect.df)
  
  #number of TDE genes for the given pairwise treatment comparison
  number.of.TDE.genes <- length(which(trtEffect.df$padj < 0.1))
  
  #percentage of TDE genes out of total
  percentage.TDE.genes <- 100*(number.of.TDE.genes/total.genes)
  
  result <- c(paste0("Total number of genes = ", total.genes),
              paste0("Number of TDE genes = ", number.of.TDE.genes),
              paste0("% of TDE genes = ", percentage.TDE.genes))
  
  return(result)
}


######################################################################################################
### Table 2 & S2: Fisher's exact tests of proportion of TDE genes in each sex bias category ##########
######################################################################################################

#first define path for where Treatment Effect dataframes are stored
path_2 <- 'C:/Users/mishr/OneDrive/Desktop/ProjectII/Annotated_R_scripts_and_data/Treatment.Effect.Dataframes/'

#then set path and load gonad specificity data
#this data will be used to make Table S2, which is a version of Table 2 without genes with 
#highly gonad-specific expression
path_3 <- 'C:/Users/mishr/OneDrive/Desktop/ProjectII/Annotated_R_scripts_and_data/Gonad.Specific.Expression/'
GSI.females <- read.table(paste0(path_3, "Gonad.Specificity.Index.Females.tsv"), header = T)
GSI.males <- read.table(paste0(path_3, "Gonad.Specificity.Index.Males.tsv"), header = T)

gonad.specific.genes <- union(GSI.females[which(GSI.females$Gonad.Specificity.Index >=0.95),]$geneID, 
                              GSI.males[which(GSI.males$Gonad.Specificity.Index >=0.95),]$geneID)


#Define function to run all the Fisher's exact tests for a given sample (e.g., male head)
#there are two kinds of Fisher's exact tests being performed: 1) to see if proportion of TDE genes in
#unbiased genes is sig. different from sex-biased genes, 2) to see if proportion of TDE genes in
#male-biased genes is sig. different from female-biased genes. 

#NOTE - treatment.comparison argument should be in the format "Trt1.vs.Trt2" (e.g., "MCabs.vs.MCcom)
FishersTests.of.Proportion.of.TDE.genes.in.SB.categories <- function(sex, tissue, treatment.comparison){
  #load treatment effect dataframe
  sample <- na.omit(read.table(paste0(path_2, sex, ".", tissue, ".log2FCMatingTreatment.Dataframe.tsv"), header = T, sep = '\t'))
  
  #COMMENT THIS OUT TO RUN WITH GONAD-SPECIFIC GENES (like it is in Figure 1; to make figure S1 uncomment it)
  #sample <- sample[!sample$geneID %in% gonad.specific.genes,]
  
  #reduce the above dataframe into displaying only sex bias category and adjusted p-value associated with 
  #log2FC treatment1/treatment2
  sample <- sample[,c(paste0("padj.", treatment.comparison), "sex.bias.category")]
  
  #make separate dataframes for TDE genes and non-TDE genes
  #TDE genes = "Treatment Differentially Expressed" genes, i.e., where log2 FC (trt1/trt2) has p-adjusted < 0.1
  TDE.df <- sample[which(sample[,1] < 0.1),]
  nonTDE.df <- sample[which(sample[,1] >= 0.1),]
  
  #divide TDE/non-TDE genes in unbiased (UB), sex-biased (SB), female-biased (FB) and male-biased (MB) genes
  unbiased.TDE <- TDE.df[which(TDE.df$sex.bias.category == "Unbiased"),]
  unbiased.nonTDE <- nonTDE.df[which(nonTDE.df$sex.bias.category == "Unbiased"),]
  sexbiased.TDE <- TDE.df[which(TDE.df$sex.bias.category == "Female-biased" | 
                                  TDE.df$sex.bias.category == "Male-biased"),]
  sexbiased.nonTDE <- nonTDE.df[which(nonTDE.df$sex.bias.category == "Female-biased" | 
                                        nonTDE.df$sex.bias.category == "Male-biased"),]
  fbiased.TDE <- TDE.df[which(TDE.df$sex.bias.category == "Female-biased"),]
  fbiased.nonTDE <- nonTDE.df[which(nonTDE.df$sex.bias.category == "Female-biased"),]
  mbiased.TDE <- TDE.df[which(TDE.df$sex.bias.category == "Male-biased"),]
  mbiased.nonTDE <- nonTDE.df[which(nonTDE.df$sex.bias.category == "Male-biased"),]
  
  #make a dataframe condensing the above information
  fraction.of.TDE.genes <- cbind.data.frame("fraction of TDE in UB genes" = nrow(unbiased.TDE)/(nrow(unbiased.TDE) + nrow(unbiased.nonTDE)),
                                            "No. of TDE in UB genes" = nrow(unbiased.TDE), 
                                            "No. of nonTDE in UB genes" = nrow(unbiased.nonTDE),
                                            "fraction of TDE in SB genes" = nrow(sexbiased.TDE)/(nrow(sexbiased.TDE) + nrow(sexbiased.nonTDE)),
                                            "No. of TDE in SB genes" = nrow(sexbiased.TDE), 
                                            "No. of nonTDE in SB genes" = nrow(sexbiased.nonTDE),
                                            "fraction of TDE in FB genes" = nrow(fbiased.TDE)/(nrow(fbiased.TDE) + nrow(fbiased.nonTDE)),
                                            "No. of TDE in FB genes" = nrow(fbiased.TDE), 
                                            "No. of nonTDE in FB genes" = nrow(fbiased.nonTDE),
                                            "fraction of TDE in MB genes" = nrow(mbiased.TDE)/(nrow(mbiased.TDE) + nrow(mbiased.nonTDE)),
                                            "No. of TDE in MB genes" = nrow(mbiased.TDE), 
                                            "No. of nonTDE in MB genes" = nrow(mbiased.nonTDE))
  
  
  ##now perform the Fisher's exact tests
  #2x2 Fisher's exact test between UB and SB genes
  TDE.df1 <- c(nrow(unbiased.TDE), nrow(sexbiased.TDE))
  nonTDE.df1 <- c(nrow(unbiased.nonTDE), nrow(sexbiased.nonTDE))
  df1 <- rbind.data.frame(TDE.df1, nonTDE.df1)
  colnames(df1) <- c("UB", "SB") 
  row.names(df1) <- c("TDE.df", "nonTDE.df")
  UB.vs.SB.pVal <- fisher.test(df1)$p.value
  
  #2x2 Fisher's exact test between MB and FB genes
  TDE.df2 <- c(nrow(fbiased.TDE), nrow(mbiased.TDE))
  nonTDE.df2 <- c(nrow(fbiased.nonTDE), nrow(mbiased.nonTDE))
  df2 <- rbind.data.frame(TDE.df2, nonTDE.df2)
  colnames(df2) <- c("FB", "MB") 
  row.names(df2) <- c("TDE.df", "nonTDE.df")
  MB.vs.FB.pVal <- fisher.test(df2)$p.value
  
  pValues <- cbind.data.frame("UB.vs.SB.pVal" = UB.vs.SB.pVal, "MB.vs.FB.pVal" = MB.vs.FB.pVal)
  results <- cbind.data.frame(fraction.of.TDE.genes, pValues)
  
  return(results)
}
