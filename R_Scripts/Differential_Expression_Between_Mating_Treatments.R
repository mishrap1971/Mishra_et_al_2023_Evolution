library(DESeq2)

#set path for folder where HTSeq-count outputs are stored
path <- "C:/Users/mishr/OneDrive/Desktop/ProjectII/Annotated_R_scripts_and_data/HTSeq.count.files"
setwd(path)

###########################################################################################################
#function to perform differential expression for between two mating treatments for a certain sex (M/F) and 
#tissue type (head/whole body)
###########################################################################################################

#specify path to save outputs into 
path_save_1 <- "C:/Users/mishr/OneDrive/Desktop/ProjectII/Annotated_R_scripts_and_data/Diff.Expression.Between.Mating.Treatments/"

DiffExp_Between_MatingTrts <- function(sex_tissue, treatment1, treatment2, output_name){
  #load htseq-count files for treatment1 and treatment2 for the given combination of sex and tissue type
  count.files.treatment1 <- list.files(pattern = paste0(treatment1, ".*", sex_tissue))
  count.files.treatment2 <- list.files(pattern = paste0(treatment2, ".*", sex_tissue))
  count.files <- c(count.files.treatment1, count.files.treatment2)
  
  #assign sample name by trimming the count file's extension (here it's ".tsv")
  samplenames <- gsub('.{4}$', '', count.files)
  
  #set "factors" for the differential expression model
  sampleCondition.matingtreatment <- factor(c(rep(treatment1, 4), rep(treatment2, 4)))
  sampleCondition.population <- factor(c(paste0(treatment1, "_1"), paste0(treatment1, "_2"),
                                         paste0(treatment1, "_3"), paste0(treatment1, "_4"),
                                         paste0(treatment2, "_1"), paste0(treatment2, "_2"),
                                         paste0(treatment2, "_3"), paste0(treatment2, "_4")))
  
  #construct condition table
  sampleTable <- data.frame(sampleName = samplenames, fileName = count.files, 
                            matingtreatment = sampleCondition.matingtreatment,
                            population = sampleCondition.population)
  
  #build DESeqDataSet from htseq count files
  #omit genes where mean read counts over all the 8 samples is <50
  ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable, design= ~matingtreatment)
  keep.genes <- rowMeans(counts(ddsHTSeq)) >= 50
  ddsHTSeq <- ddsHTSeq[keep.genes,]
  
  #run DESeq2
  dds.diffExp <- DESeq(ddsHTSeq)
  
  #name the DESeq2 results
  #set the "contrast" such that the log2 FoldChange in the output is treatment1/treatment2
  results.diffExp = results(dds.diffExp, contrast=c("matingtreatment", treatment1, treatment2)) 
  
  #extract the FlyBaseIDs and set it as the first column for ease
  results.diffExp$FlyBaseID <- rownames(results.diffExp)
  results.diffExp <- results.diffExp[,c(7,1:6)]
  
  #save output to file
  write.table(results.diffExp, paste0(path_save_1, output_name, ".tsv"), sep = '\t', quote = F, row.names = F, col.names = T)
  
  return(results.diffExp)
}

#run above function for all samples, except those that involve MCsim_Female_Body
DiffExp_Between_MatingTrts("Female_head", "MCabs", "MCsim", "Female.Head.MCabs.vs.MCsim")
DiffExp_Between_MatingTrts("Male_body", "MCabs", "MCsim", "Male.Body.MCabs.vs.MCsim")
DiffExp_Between_MatingTrts("Male_head", "MCabs", "MCsim", "Male.Head.MCabs.vs.MCsim")

DiffExp_Between_MatingTrts("Female_head", "MCsim", "MCcom", "Female.Head.MCsim.vs.MCcom")
DiffExp_Between_MatingTrts("Male_body", "MCsim", "MCcom", "Male.Body.MCsim.vs.MCcom")
DiffExp_Between_MatingTrts("Male_head", "MCsim", "MCcom", "Male.Head.MCsim.vs.MCcom")

DiffExp_Between_MatingTrts("Female_body", "MCabs", "MCcom", "Female.Body.MCabs.vs.MCcom")
DiffExp_Between_MatingTrts("Female_head", "MCabs", "MCcom", "Female.Head.MCabs.vs.MCcom")
DiffExp_Between_MatingTrts("Male_body", "MCabs", "MCcom", "Male.Body.MCabs.vs.MCcom")
DiffExp_Between_MatingTrts("Male_head", "MCabs", "MCcom", "Male.Head.MCabs.vs.MCcom")


#################################################################################################################
#A replicate population from MCsim_Female_Body was contaminated with male tissue, and                ##
#hence had to be removed. Differential Expression involving this particular sample need an                     ##
#edited version of the previous function.                                                                      ##
#################################################################################################################

#######################################################
#   For MCabs vs MCsim, Female Body   ####
#######################################################
count.files.MCabs <- list.files(pattern = paste0("MCabs", ".*", "Female_body"))
count.files.MCsim <- list.files(pattern = paste0("MCsim", ".*", "Female_body"))
count.files <- c(count.files.MCabs, count.files.MCsim)

#assign sample name by trimming the count file's extension (here it's ".tsv")
samplenames <- gsub('.{4}$', '', count.files)

#set "factors" for the differential expression model
sampleCondition.matingtreatment <- factor(c(rep("MCabs", 4), rep("MCsim", 3)))
sampleCondition.population <- factor(c("MCabs_1", "MCabs_2", "MCabs_3", "MCabs_4",
                                       "MCsim_2", "MCsim_3", "MCsim_4"))

#construct condition table
sampleTable <- data.frame(sampleName = samplenames, fileName = count.files, 
                          matingtreatment = sampleCondition.matingtreatment,
                          population = sampleCondition.population)

#build DESeqDataSet from htseq count files
#omit genes where mean read counts over all the 8 samples is <50
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable, design= ~matingtreatment)
keep.genes <- rowMeans(counts(ddsHTSeq)) >= 50
ddsHTSeq <- ddsHTSeq[keep.genes,]

#run DESeq2
dds.diffExp <- DESeq(ddsHTSeq)

#name the DESeq2 results
#set the "contrast" such that the log2 FoldChange in the output is MCabs/MCsim
results.diffExp = results(dds.diffExp, contrast=c("matingtreatment", "MCabs", "MCsim")) 

#extract the FlyBaseIDs and set it as the first column for ease
results.diffExp$FlyBaseID <- rownames(results.diffExp)
results.diffExp <- results.diffExp[,c(7,1:6)]

#save output to file
write.table(results.diffExp, paste0(path_save_1, "Female.Body.MCabs.vs.MCsim.tsv"), 
            sep = '\t', quote = F, row.names = F, col.names = T)

################################################################
#   For MCsim vs MCcom , Female Body   ####
################################################################
count.files.MCsim <- list.files(pattern = paste0("MCsim", ".*", "Female_body"))
count.files.MCcom <- list.files(pattern = paste0("MCcom", ".*", "Female_body"))
count.files <- c(count.files.MCsim, count.files.MCcom)

#assign sample name by trimming the count file's extension (here it's ".tsv")
samplenames <- gsub('.{4}$', '', count.files)

#set "factors" for the differential expression model
sampleCondition.matingtreatment <- factor(c(rep("MCsim", 3), rep("MCcom", 4)))
sampleCondition.population <- factor(c("MCsim_2", "MCsim_3", "MCsim_4",
                                       "MCcom_1", "MCcom_2", "MCcom_3", 
                                       "MCcom_4"))

#construct condition table
sampleTable <- data.frame(sampleName = samplenames, fileName = count.files, 
                          matingtreatment = sampleCondition.matingtreatment,
                          population = sampleCondition.population)

#build DESeqDataSet from htseq count files
#omit genes where mean read counts over all the 8 samples is <50
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable, design= ~matingtreatment)
keep.genes <- rowMeans(counts(ddsHTSeq)) >= 50
ddsHTSeq <- ddsHTSeq[keep.genes,]

#run DESeq2
dds.diffExp <- DESeq(ddsHTSeq)

#name the DESeq2 results
#set the "contrast" such that the log2 FoldChange in the output is MCabs/MCsim
results.diffExp = results(dds.diffExp, contrast=c("matingtreatment", "MCsim", "MCcom")) 

#extract the FlyBaseIDs and set it as the first column for ease
results.diffExp$FlyBaseID <- rownames(results.diffExp)
results.diffExp <- results.diffExp[,c(7,1:6)]

#save output to file
write.table(results.diffExp, paste0(path_save_1, "Female.Body.MCsim.vs.MCcom.tsv"), 
            sep = '\t', quote = F, row.names = F, col.names = T)

################################################################################################################
################################################################################################################
################################################################################################################

#Just to simplify downstream codes: need to make dataframes for each tissue-sex combination (e.g., female body), 
#containing the following - treatment effect (i.e., log2FC trt1/trt2) for all 3 pairwise comparison of 
#treatments, plus sex-bias category and log2FC male/female for each gene (from external data, Osada et al 2017) 

#before defining function to make such dataframes, define path for: a) where the external 
#Sex-Biased Gene Expression data is stored, b) where DESeq2 outputs are stored, and
#c) define path to save output dataframe
path_for_SBGE_data <- "C:/Users/mishr/OneDrive/Desktop/ProjectII/Annotated_R_scripts_and_data/Sex.Biased.Gene.Expression.Osada_et_al_2017/"
path_save_2 <- "C:/Users/mishr/OneDrive/Desktop/ProjectII/Annotated_R_scripts_and_data/Treatment.Effect.Dataframes/"
path_for_DESeq2_data <- "C:/Users/mishr/OneDrive/Desktop/ProjectII/Annotated_R_scripts_and_data/Diff.Expression.Between.Mating.Treatments/"

TreatmentEffectDataframes <- function(sex, tissue){
  #load SBGE data for tissue (body/head)
  sbge.df <- read.table(paste0(path_for_SBGE_data, 'Sex.Biased.Genes.', tissue, '.tsv'), header = T)
  
  #load DESeq2 outputs for each pairwise comparison of mating treatments
  MCabs.vs.MCsim <- read.table(paste0(path_for_DESeq2_data, sex, ".", tissue, ".MCabs.vs.MCsim.tsv"), header = T)
  MCabs.vs.MCcom <- read.table(paste0(path_for_DESeq2_data, sex, ".", tissue, ".MCabs.vs.MCcom.tsv"), header = T)
  MCsim.vs.MCcom <- read.table(paste0(path_for_DESeq2_data, sex, ".", tissue, ".MCsim.vs.MCcom.tsv"), header = T)
  
  #find genes that have data in all of the above four dataframes
  common.genes <- intersect(MCabs.vs.MCsim$FlyBaseID, intersect(MCabs.vs.MCcom$FlyBaseID, MCsim.vs.MCcom$FlyBaseID))
  common.genes <- intersect(sbge.df$geneID, common.genes)
  
  #subset the 4 dataframes to have common genes from above
  MCabs.vs.MCsim <- MCabs.vs.MCsim[MCabs.vs.MCsim$FlyBaseID %in% common.genes,]
  MCabs.vs.MCcom <- MCabs.vs.MCcom[MCabs.vs.MCcom$FlyBaseID %in% common.genes,]
  MCsim.vs.MCcom <- MCsim.vs.MCcom[MCsim.vs.MCcom$FlyBaseID %in% common.genes,]
  sbge.df <- sbge.df[sbge.df$geneID %in% common.genes,]
  
  #construct dataframe with the following: 1) FlyBase geneID, 2) log2FC (MCabs/MCsim), 3) p-value(MCabs/MCsim),
  #4) adjusted p-value(MCabs/MCsim), 5) log2FC (MCabs/MCcom), 6) p-value(MCabs/MCcom), 
  #7) adjusted p-value(MCabs/MCcom), 8) log2FC (MCsim/MCcom), 9) p-value(MCsim/MCcom),
  #10) adjusted p-value(MCsim/MCcom), 11) log2FC (male/female) from Osada et al. data,
  #12) sex bias category from Osada et al. data
  df <- cbind.data.frame(geneID = common.genes, log2FC.MCabs.vs.MCsim = MCabs.vs.MCsim$log2FoldChange, 
                         pVal.MCabs.vs.MCsim = MCabs.vs.MCsim$pvalue, padj.MCabs.vs.MCsim = MCabs.vs.MCsim$padj, 
                         log2FC.MCabs.vs.MCcom = MCabs.vs.MCcom$log2FoldChange, pVal.MCabs.vs.MCcom = MCabs.vs.MCcom$pvalue, 
                         padj.MCabs.vs.MCcom = MCabs.vs.MCcom$padj, log2FC.MCsim.vs.MCcom = MCsim.vs.MCcom$log2FoldChange,
                         pval.MCsim.vs.MCcom = MCsim.vs.MCcom$pvalue, padj.MCsim.vs.MCcom = MCsim.vs.MCcom$padj, 
                         log2FC.sex.Osada = sbge.df$log2FoldChange,
                         sex.bias.category = sbge.df$Sex.Bias.Category)
  
  #save above dataframe to file
  write.table(df, paste0(path_save_2, sex, ".", tissue, ".log2FCMatingTreatment.Dataframe.tsv"), 
              sep = '\t', quote = F, row.names = F, col.names = T)
  
  return(df)
}

TreatmentEffectDataframes("Female", "Body")
TreatmentEffectDataframes("Female", "Head")
TreatmentEffectDataframes("Male", "Body")
TreatmentEffectDataframes("Male", "Head")


