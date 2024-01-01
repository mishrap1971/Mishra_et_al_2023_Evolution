#This script performs much of the analysis involving TDS (Treatment Differentially Spliced) Genes, 
#from Table S3 to S6

#Set global variables
FDRThreshold = 0.1
mappedReadsThreshold = 50

#set path to folder where JunctionSeq outputs are stored (for Diff. Splicing between Treatments within one sex) 
path_for_JS_data <- "C:/Users/mishr/OneDrive/Desktop/ProjectII/Annotated_R_scripts_and_data/JunctionSeq_Differential_Splicing_Outputs/Within.Sex.Between.Treatments/"


######################################################################################################
### Table S3: Proportion of TDS genes for each sex in each pairwise comparison of treatments #########
######################################################################################################

#define function to yield %TDS genes in each pairwise comparison of treatments
#e.g., sex = "Female", tissue = "Male", treatment.comparison = "MCsim.vs.MCcom"
Percentage.of.TDS.genes <- function(sex, tissue, treatment.comparison){
  #load JunctionSeq output
  JS.output <- read.table(paste0(path_for_JS_data, treatment.comparison, ".", sex, ".", tissue, "/allGenes.results.txt"), sep = "\t", header = TRUE)
  
  #remove rows for genes where expression data is "NA" for either treatment or 
  #expression is less than the coverage cut-off (mappedReadsThreshold) 
  JS.output <- JS.output[!(is.na(JS.output[,23]) | is.na(JS.output[,24])),]
  JS.output <- JS.output[JS.output[,23] > mappedReadsThreshold & JS.output[,24] > mappedReadsThreshold,]
  
  #total number of genes
  total.genes <- length(unique(JS.output$geneID))
  
  #unique list of geneIDs with significant differential splicing (padj < 0.1)
  sig.genes.JS.output <- unique(JS.output$geneID[JS.output$geneWisePadj < FDRThreshold])
  number.of.TDS.genes <- length(sig.genes.JS.output)
  
  #percentage of TDS genes out of total
  percentage.TDS.genes <- 100*(length(sig.genes.JS.output)/length(unique(JS.output$geneID)))
  
  result <- c(paste0("Total number of genes = ", total.genes),
              paste0("Number of TDS genes = ", number.of.TDS.genes),
              paste0("% of TDS genes = ", percentage.TDS.genes))
  
  return(result)
}

Percentage.of.TDS.genes("Female", "Body", "MCabs.vs.MCsim")
Percentage.of.TDS.genes("Female", "Body", "MCabs.vs.MCcom")
Percentage.of.TDS.genes("Female", "Body", "MCsim.vs.MCcom")

Percentage.of.TDS.genes("Male", "Body", "MCabs.vs.MCsim")
Percentage.of.TDS.genes("Male", "Body", "MCabs.vs.MCcom")
Percentage.of.TDS.genes("Male", "Body", "MCsim.vs.MCcom")

Percentage.of.TDS.genes("Female", "Head", "MCabs.vs.MCsim")
Percentage.of.TDS.genes("Female", "Head", "MCabs.vs.MCcom")
Percentage.of.TDS.genes("Female", "Head", "MCsim.vs.MCcom")

Percentage.of.TDS.genes("Male", "Head", "MCabs.vs.MCsim")
Percentage.of.TDS.genes("Male", "Head", "MCabs.vs.MCcom")
Percentage.of.TDS.genes("Male", "Head", "MCsim.vs.MCcom")




############################################################################################################
## make dataframes that combine Differential Splicing between mating treatments (from JunctionSeq)        ##
## with sex bias/SBGE (Osada) data and Treatment Effect Dataframes  [for Tables S4-S5]                    ##
############################################################################################################
#set path and load SBGE data
path_for_SBGE_data <- "C:/Users/mishr/OneDrive/Desktop/ProjectII/Annotated_R_scripts_and_data/Sex.Biased.Gene.Expression.Osada_et_al_2017/"
sexbias.body <- read.table(paste0(path_for_SBGE_data, "Sex.Biased.Genes.Body.tsv"), header = T)
sexbias.head <- read.table(paste0(path_for_SBGE_data, "Sex.Biased.Genes.Head.tsv"), header = T)

#set path for saving output dataframes
path_save <- "C:/Users/mishr/OneDrive/Desktop/ProjectII/Annotated_R_scripts_and_data/Diff.Splicing.Dataframes/"

#define function to make dataframes
Dataframes_DiffSplicing_between_Trts <- function(sex, tissue){
  #load sex bias data
  sexbias.df <- read.table(paste0(path_for_SBGE_data, "Sex.Biased.Genes.", tissue, ".tsv"), header = T)
  
  #load JunctionSeq outputs for each of the pairwise treatment comparisons
  #remove rows for genes where expression data is "NA" for either treatment or 
  #expression is less than the coverage cut-off (mappedReadsThreshold) 
  MCabs.vs.MCsim <- read.table(paste0(path_for_JS_data, "MCabs.vs.MCsim.", sex, ".", tissue, "/allGenes.results.txt"), sep = "\t", header = TRUE)
  MCabs.vs.MCsim <- MCabs.vs.MCsim[!(is.na(MCabs.vs.MCsim[,23]) | is.na(MCabs.vs.MCsim[,24])),]
  MCabs.vs.MCsim <- MCabs.vs.MCsim[MCabs.vs.MCsim[,23] > mappedReadsThreshold & MCabs.vs.MCsim[,24] > mappedReadsThreshold,]
  
  MCabs.vs.MCcom <- read.table(paste0(path_for_JS_data, "MCabs.vs.MCcom.", sex, ".", tissue, "/allGenes.results.txt"), sep = "\t", header = TRUE)
  MCabs.vs.MCcom <- MCabs.vs.MCcom[!(is.na(MCabs.vs.MCcom[,23]) | is.na(MCabs.vs.MCcom[,24])),]
  MCabs.vs.MCcom <- MCabs.vs.MCcom[MCabs.vs.MCcom[,23] > mappedReadsThreshold & MCabs.vs.MCcom[,24] > mappedReadsThreshold,]
  
  MCsim.vs.MCcom <- read.table(paste0(path_for_JS_data, "MCsim.vs.MCcom.", sex, ".", tissue, "/allGenes.results.txt"), sep = "\t", header = TRUE)
  MCsim.vs.MCcom <- MCsim.vs.MCcom[!(is.na(MCsim.vs.MCcom[,23]) | is.na(MCsim.vs.MCcom[,24])),]
  MCsim.vs.MCcom <- MCsim.vs.MCcom[MCsim.vs.MCcom[,23] > mappedReadsThreshold & MCsim.vs.MCcom[,24] > mappedReadsThreshold,]
  
  
  #Obtain a list of genes with significant and non-significant differential splicing for each 
  #treatment comparison, then make simplified dataframes
  sig.MCabs.vs.MCsim <- unique(MCabs.vs.MCsim$geneID[MCabs.vs.MCsim$geneWisePadj < FDRThreshold])
  nonsig.MCabs.vs.MCsim <- unique(MCabs.vs.MCsim$geneID[MCabs.vs.MCsim$geneWisePadj > FDRThreshold])
  sig.MCabs.vs.MCsim.df <- cbind.data.frame(geneID = sig.MCabs.vs.MCsim, Diff.Splicing = rep("Significant", length(sig.MCabs.vs.MCsim)))
  nonsig.MCabs.vs.MCsim.df <- cbind.data.frame(geneID=nonsig.MCabs.vs.MCsim, Diff.Splicing=rep("Not Significant", length(nonsig.MCabs.vs.MCsim)))
  MCabs.vs.MCsim.df <- rbind.data.frame(sig.MCabs.vs.MCsim.df, nonsig.MCabs.vs.MCsim.df)
  
  sig.MCabs.vs.MCcom <- unique(MCabs.vs.MCcom$geneID[MCabs.vs.MCcom$geneWisePadj < FDRThreshold])
  nonsig.MCabs.vs.MCcom <- unique(MCabs.vs.MCcom$geneID[MCabs.vs.MCcom$geneWisePadj > FDRThreshold])
  sig.MCabs.vs.MCcom.df <- cbind.data.frame(geneID = sig.MCabs.vs.MCcom, Diff.Splicing = rep("Significant", length(sig.MCabs.vs.MCcom)))
  nonsig.MCabs.vs.MCcom.df <- cbind.data.frame(geneID=nonsig.MCabs.vs.MCcom, Diff.Splicing=rep("Not Significant", length(nonsig.MCabs.vs.MCcom)))
  MCabs.vs.MCcom.df <- rbind.data.frame(sig.MCabs.vs.MCcom.df, nonsig.MCabs.vs.MCcom.df)
  
  sig.MCsim.vs.MCcom <- unique(MCsim.vs.MCcom$geneID[MCsim.vs.MCcom$geneWisePadj < FDRThreshold])
  nonsig.MCsim.vs.MCcom <- unique(MCsim.vs.MCcom$geneID[MCsim.vs.MCcom$geneWisePadj > FDRThreshold])
  sig.MCsim.vs.MCcom.df <- cbind.data.frame(geneID = sig.MCsim.vs.MCcom, Diff.Splicing = rep("Significant", length(sig.MCsim.vs.MCcom)))
  nonsig.MCsim.vs.MCcom.df <- cbind.data.frame(geneID=nonsig.MCsim.vs.MCcom, Diff.Splicing=rep("Not Significant", length(nonsig.MCsim.vs.MCcom)))
  MCsim.vs.MCcom.df <- rbind.data.frame(sig.MCsim.vs.MCcom.df, nonsig.MCsim.vs.MCcom.df)
  
  #get list of genes tested in all 3 differential splicing comparisons + SBGE analysis; retain those genes
  gene.list <- Reduce(intersect, list(sexbias.df$geneID, MCabs.vs.MCsim.df$geneID, MCabs.vs.MCcom.df$geneID, MCsim.vs.MCcom.df$geneID))  
  sexbias.df <- sexbias.df[sexbias.df$geneID %in% gene.list,]
  MCabs.vs.MCsim.df <- MCabs.vs.MCsim.df[MCabs.vs.MCsim.df$geneID %in% gene.list,]
  MCabs.vs.MCcom.df <- MCabs.vs.MCcom.df[MCabs.vs.MCcom.df$geneID %in% gene.list,]
  MCsim.vs.MCcom.df <- MCsim.vs.MCcom.df[MCsim.vs.MCcom.df$geneID %in% gene.list,]
  
  #sort above dataframes by geneID
  sexbias.df <- sexbias.df[order(sexbias.df$geneID),]
  MCabs.vs.MCsim.df <- MCabs.vs.MCsim.df[order(MCabs.vs.MCsim.df$geneID),]
  MCabs.vs.MCcom.df <- MCabs.vs.MCcom.df[order(MCabs.vs.MCcom.df$geneID),]
  MCsim.vs.MCcom.df <- MCsim.vs.MCcom.df[order(MCsim.vs.MCcom.df$geneID),]
  
  
  #make final dataframe
  df <- cbind.data.frame(sexbias.df, 
                         Diff.Splicing.MCabs.vs.MCsim = MCabs.vs.MCsim.df$Diff.Splicing,
                         Diff.Splicing.MCabs.vs.MCcom = MCabs.vs.MCcom.df$Diff.Splicing, 
                         Diff.Splicing.MCsim.vs.MCcom = MCsim.vs.MCcom.df$Diff.Splicing)
  
  #save to file
  write.table(df, paste0(path_save, sex, ".", tissue, '.DiffSplicing.SBGE.dataframe.tsv'), quote = F, sep = '\t', col.names = T, row.names = F)
  
  return(df)
}

Dataframes_DiffSplicing_between_Trts("Female", "Body")
Dataframes_DiffSplicing_between_Trts("Male", "Body")
Dataframes_DiffSplicing_between_Trts("Female", "Head")
Dataframes_DiffSplicing_between_Trts("Male", "Head")


######################################################################################################
### Tables S4 & S5: Fisher's exact tests of proportion of TDS genes in each sex bias category ########
######################################################################################################

#first define path for where Differential Splicing dataframes are stored
path_dataframes <- "C:/Users/mishr/OneDrive/Desktop/ProjectII/Annotated_R_scripts_and_data/Diff.Splicing.Dataframes/"

#then set path and load gonad specificity data
#this data will be used to make Table S5, which is a version of Table S4 without genes with 
#highly gonad-specific expression
path_3 <- 'C:/Users/mishr/OneDrive/Desktop/ProjectII/Annotated_R_scripts_and_data/Gonad.Specific.Expression/'
GSI.females <- read.table(paste0(path_3, "Gonad.Specificity.Index.Females.tsv"), header = T)
GSI.males <- read.table(paste0(path_3, "Gonad.Specificity.Index.Males.tsv"), header = T)

gonad.specific.genes <- union(GSI.females[which(GSI.females$Gonad.Specificity.Index >=0.95),]$geneID, 
                              GSI.males[which(GSI.males$Gonad.Specificity.Index >=0.95),]$geneID)


#Define function to run all the Fisher's exact tests for a given sample (e.g., male head)
#there are two kinds of Fisher's exact tests being performed: 1) to see if proportion of TDS genes in
#unbiased genes is sig. different from sex-biased genes, 2) to see if proportion of TDS genes in
#male-biased genes is sig. different from female-biased genes. 

#NOTE - treatment.comparison argument should be in the format "Trt1.vs.Trt2" (e.g., "MCabs.vs.MCcom)
FishersTests.of.Proportion.of.TDS.genes.in.SB.categories <- function(sex, tissue, treatment.comparison){
  #load treatment effect dataframe
  sample <- na.omit(read.table(paste0(path_dataframes, sex, ".", tissue, ".DiffSplicing.SBGE.dataframe.tsv"), header = T, sep = '\t'))
  
  #COMMENT THIS OUT TO RUN WITH GONAD-SPECIFIC GENES (for Table S5)
  #sample <- sample[!sample$geneID fraction ofinfraction of gonad.specific.genes,]
  
  #reduce the above dataframe into displaying only sex bias category and adjusted p-value associated with 
  #log2FC treatment1/treatment2
  sample <- sample[,c(paste0("Diff.Splicing.", treatment.comparison), "Sex.Bias.Category")]
  
  #make separate dataframes for TDS genes and non-TDS genes
  #TDS genes = "Treatment Differentially Spliced" genes, i.e., where p-adjusted < 0.1
  TDS.df <- sample[which(sample[,1] == "Significant"),]
  nonTDS.df <- sample[which(sample[,1] == "Not Significant"),]
  
  #divide DE/non-DE genes in unbiased (UB), sex-biased (SB), female-biased (FB) and male-biased (MB) genes
  unbiased.TDS <- TDS.df[which(TDS.df$Sex.Bias.Category == "Unbiased"),]
  unbiased.nonTDS <- nonTDS.df[which(nonTDS.df$Sex.Bias.Category == "Unbiased"),]
  sexbiased.TDS <- TDS.df[which(TDS.df$Sex.Bias.Category == "Female-biased" | 
                                  TDS.df$Sex.Bias.Category == "Male-biased"),]
  sexbiased.nonTDS <- nonTDS.df[which(nonTDS.df$Sex.Bias.Category == "Female-biased" | 
                                        nonTDS.df$Sex.Bias.Category == "Male-biased"),]
  fbiased.TDS <- TDS.df[which(TDS.df$Sex.Bias.Category == "Female-biased"),]
  fbiased.nonTDS <- nonTDS.df[which(nonTDS.df$Sex.Bias.Category == "Female-biased"),]
  mbiased.TDS <- TDS.df[which(TDS.df$Sex.Bias.Category == "Male-biased"),]
  mbiased.nonTDS <- nonTDS.df[which(nonTDS.df$Sex.Bias.Category == "Male-biased"),]
  
  #make a dataframe condensing the above information
  fraction.of.TDS.genes <- cbind.data.frame("fraction of TDS in UB genes" = nrow(unbiased.TDS)/(nrow(unbiased.TDS) + nrow(unbiased.nonTDS)),
                                            "No. of TDS in UB genes" = nrow(unbiased.TDS), 
                                            "No. of nonTDS in UB genes" = nrow(unbiased.nonTDS),
                                            "fraction of TDS in SB genes" = nrow(sexbiased.TDS)/(nrow(sexbiased.TDS) + nrow(sexbiased.nonTDS)),
                                            "No. of TDS in SB genes" = nrow(sexbiased.TDS), 
                                            "No. of nonTDS in SB genes" = nrow(sexbiased.nonTDS),
                                            "fraction of TDS in FB genes" = nrow(fbiased.TDS)/(nrow(fbiased.TDS) + nrow(fbiased.nonTDS)),
                                            "No. of TDS in FB genes" = nrow(fbiased.TDS), 
                                            "No. of nonTDS in FB genes" = nrow(fbiased.nonTDS),
                                            "fraction of TDS in MB genes" = nrow(mbiased.TDS)/(nrow(mbiased.TDS) + nrow(mbiased.nonTDS)),
                                            "No. of TDS in MB genes" = nrow(mbiased.TDS), 
                                            "No. of nonTDS in MB genes" = nrow(mbiased.nonTDS))
  
  
  ##now perform the Fisher's exact tests
  #2x2 Fisher's exact test between UB and SB genes
  TDS.df1 <- c(nrow(unbiased.TDS), nrow(sexbiased.TDS))
  nonTDS.df1 <- c(nrow(unbiased.nonTDS), nrow(sexbiased.nonTDS))
  df1 <- rbind.data.frame(TDS.df1, nonTDS.df1)
  colnames(df1) <- c("UB", "SB") 
  row.names(df1) <- c("TDS.df", "nonTDS.df")
  UB.vs.SB.pVal <- fisher.test(df1)$p.value
  
  #2x2 Fisher's exact test between MB and FB genes
  TDS.df2 <- c(nrow(fbiased.TDS), nrow(mbiased.TDS))
  nonTDS.df2 <- c(nrow(fbiased.nonTDS), nrow(mbiased.nonTDS))
  df2 <- rbind.data.frame(TDS.df2, nonTDS.df2)
  colnames(df2) <- c("FB", "MB") 
  row.names(df2) <- c("TDS.df", "nonTDS.df")
  MB.vs.FB.pVal <- fisher.test(df2)$p.value
  
  pValues <- cbind.data.frame("UB.vs.SB.pVal" = UB.vs.SB.pVal, "MB.vs.FB.pVal" = MB.vs.FB.pVal)
  results <- cbind.data.frame(fraction.of.TDS.genes, pValues)
  
  return(results)
}

FishersTests.of.Proportion.of.TDS.genes.in.SB.categories("Female", "Body", "MCabs.vs.MCsim")
FishersTests.of.Proportion.of.TDS.genes.in.SB.categories("Female", "Body", "MCabs.vs.MCcom")
FishersTests.of.Proportion.of.TDS.genes.in.SB.categories("Female", "Body", "MCabs.vs.MCcom")

FishersTests.of.Proportion.of.TDS.genes.in.SB.categories("Female", "Head", "MCabs.vs.MCsim")
FishersTests.of.Proportion.of.TDS.genes.in.SB.categories("Female", "Head", "MCabs.vs.MCcom")
FishersTests.of.Proportion.of.TDS.genes.in.SB.categories("Female", "Head", "MCabs.vs.MCcom")

FishersTests.of.Proportion.of.TDS.genes.in.SB.categories("Male", "Body", "MCabs.vs.MCsim")
FishersTests.of.Proportion.of.TDS.genes.in.SB.categories("Male", "Body", "MCabs.vs.MCcom")
FishersTests.of.Proportion.of.TDS.genes.in.SB.categories("Male", "Body", "MCabs.vs.MCcom")

FishersTests.of.Proportion.of.TDS.genes.in.SB.categories("Male", "Head", "MCabs.vs.MCsim")
FishersTests.of.Proportion.of.TDS.genes.in.SB.categories("Male", "Head", "MCabs.vs.MCcom")
FishersTests.of.Proportion.of.TDS.genes.in.SB.categories("Male", "Head", "MCabs.vs.MCcom")


######################################################################################################
### Tables S6: Fisher's exact tests of overlap between of TDE and TDS genes in each pairwise  ########
###         comparison of mating treatments                                                   ########
######################################################################################################

#define path for where Treatment Effect dataframes are stored
path_dataframes <- 'C:/Users/mishr/OneDrive/Desktop/ProjectII/Annotated_R_scripts_and_data/Treatment.Effect.Dataframes/'

#function to perform Fisher's exact tests
Fishers.Tests.Overlap.in.TDE.and.TDS.genes <- function(sex, tissue, treatment.comparison){
  #load JunctionSeq output
  JS.output <- read.table(paste0(path_for_JS_data, treatment.comparison, ".", sex, ".", tissue, "/allGenes.results.txt"), sep = "\t", header = TRUE)
  
  #load Treatment Effect dataframe
  trtEffect.df <- na.omit(read.table(paste0(path_dataframes, sex, ".", tissue, ".log2FCMatingTreatment.Dataframe.tsv"), header = T, sep = '\t'))
  trtEffect.df <- trtEffect.df[,c("geneID", paste0("padj.", treatment.comparison))]
  
  #remove rows for genes where expression data is "NA" for either treatment or 
  #expression is less than the coverage cut-off (mappedReadsThreshold) 
  JS.output <- JS.output[!(is.na(JS.output[,23]) | is.na(JS.output[,24])),]
  JS.output <- JS.output[JS.output[,23] > mappedReadsThreshold & JS.output[,24] > mappedReadsThreshold,]
  
  #make a list of genes tested in both analyses (Diff. Expression & Diff. Splicing) and retain
  #only those genes within both dataframes
  gene.list <- intersect(unique(trtEffect.df$geneID), unique(JS.output$geneID))
  trtEffect.df <- trtEffect.df[trtEffect.df$geneID %in% gene.list,]
  JS.output <- JS.output[JS.output$geneID %in% gene.list,]
  
  #percent significantly differentially spliced genes; sig. adjusted p-value < 0.1
  sig.genes.diff.splicing <- unique(JS.output$geneID[JS.output$geneWisePadj < FDRThreshold])
  nonsig.genes.diff.splicing <- unique(JS.output$geneID[JS.output$geneWisePadj >= FDRThreshold])
  
  #percent significantly differential expressed genes; sig. adjusted p-value < 0.1
  sig.genes.diff.exp <- trtEffect.df$geneID[trtEffect.df[,2] >= 0.1]
  nonsig.genes.diff.exp <- trtEffect.df$geneID[trtEffect.df[,2] >= 0.1]
  
  sig.in.both <- intersect(sig.genes.diff.splicing, sig.genes.diff.exp)
  sig.in.splicing.not.expression <- intersect(sig.genes.diff.splicing, nonsig.genes.diff.exp)
  sig.in.expression.not.splicing <- intersect(sig.genes.diff.exp, nonsig.genes.diff.exp)
  sig.in.neither <- intersect(nonsig.genes.diff.splicing, nonsig.genes.diff.exp)
  
  (length(sig.in.both)/length(gene.list))*100
  (length(sig.in.splicing.not.expression)/length(gene.list))*100
  
  #fisher's exact test
  sig.genes.splicing <- c(length(sig.in.both), length(sig.in.splicing.not.expression))
  nonsig.genes.splicing <- c(length(sig.in.expression.not.splicing), length(sig.in.neither))
  df <- rbind.data.frame(sig.genes.splicing, nonsig.genes.splicing)
  colnames(df) <- c("Sig.expression", "NonSig.expression") 
  row.names(df) <- c("Sig.splicing", "NonSig.splicing")
  test <- fisher.test(df)
  
  results <- c(paste0("Total number of genes = ", length(gene.list)),
               paste0("Number of TDS genes = ", length(sig.genes.diff.splicing)),
               paste0("Number of TDE genes = ", length(sig.genes.diff.exp)),
               paste0("Number of non-TDE genes = ", length(nonsig.genes.diff.exp)),
               paste0("Number of genes that are TDE and TDS = ", length(sig.in.both)),
               paste0("Number of genes that are non-TDE and TDS = ", length(sig.in.splicing.not.expression)),
               paste0("Fisher's Test p-value = ", test$p.value))
  
  return(results)
  
}

