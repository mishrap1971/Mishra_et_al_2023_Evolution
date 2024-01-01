library(DESeq2)

#set path for folder where HTSeq-count outputs are stored
path <- "C:/Users/mishr/OneDrive/Desktop/ProjectII/Annotated_R_scripts_and_data/HTSeq.count.files"
setwd(path)

########################################################################################################
#function to perform differential expression between the two sexes for a certain mating treatment in a 
#certain tissue type (head/whole body)
########################################################################################################

#specify path to save outputs into 
path_save <- "C:/Users/mishr/OneDrive/Desktop/ProjectII/Annotated_R_scripts_and_data/Diff.Expression.Between.Sexes.Within.Mating.Treatment/"

DiffExp_Between_Sexes <- function(treatment, tissue){
  #load htseq-count files for the two sexes for a given combination of mating-treatment and tissue type
  count.files.males <- list.files(pattern = paste0(treatment, ".*", "Male_" ,tissue))
  count.files.females <- list.files(pattern = paste0(treatment, ".*", "Female_" ,tissue))
  count.files <- c(count.files.males, count.files.females)
  
  #assign sample name by trimming the count file's extension (".tsv")
  samplenames <- gsub('.{4}$', '', count.files)
  
  #set "factors" for the differential expression model
  sampleCondition.Sex <- factor(c(rep("Male", 4), rep("Female", 4)))
  sampleCondition.Population <- factor(rep(c("MCabs_1", "MCabs_2", "MCabs_3", "MCabs_4"), 2))
  
  #construct condition table
  sampleTable <- data.frame(sampleName = samplenames, fileName = count.files, 
                            sex = sampleCondition.Sex,
                            population = sampleCondition.Population)
  
  #build DESeqDataSet from htseq count files
  #omit genes where mean read counts over all the 8 samples is <50
  ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable, design= ~sex + population)
  keep.genes <- rowMeans(counts(ddsHTSeq)) >= 50
  ddsHTSeq <- ddsHTSeq[keep.genes,]
  
  #run DESeq2
  dds.diffExp <- DESeq(ddsHTSeq)
  
  #name the DESeq2 results
  #set the "contrast" such that the log2 FoldChange in the output is Males/Females
  results.diffExp = results(dds.diffExp, contrast=c("sex", "Male", "Female")) 
  
  #extract the FlyBaseIDs and set it as the first column for ease
  results.diffExp$FlyBaseID <- rownames(results.diffExp)
  results.diffExp <- results.diffExp[,c(7,1:6)]
  
  #save output to file
  write.table(results.diffExp, paste0(path_save, treatment, "_", tissue, "_Diff_Exp_Males_vs_Females.tsv"), 
              sep = '\t', quote = F, row.names = F, col.names = T)
  
  return(results.diffExp)
}

#run above function for all samples, except those that involve MCsim_Female_Body
DiffExp_Between_Sexes("MCabs", "body")
DiffExp_Between_Sexes("MCcom", "body")

DiffExp_Between_Sexes("MCabs", "head")
DiffExp_Between_Sexes("MCsim", "head")
DiffExp_Between_Sexes("MCcom", "head")

#################################################################################################################
#A replicate population from MCsim_Female_Body was contaminated with male tissue, and                          ##
#hence had to be removed. Differential Expression involving this particular sample need an                     ##
#edited version of the previous function.                                                                      ##
#################################################################################################################

#load htseq-count files for the two sexes for a given combination of mating-treatment and tissue type
count.files.males <- list.files(pattern = paste0("MCsim", ".*", "Male_" , "body"))
count.files.females <- list.files(pattern = paste0("MCsim", ".*", "Female_" , "body"))
count.files <- c(count.files.males, count.files.females)

#assign sample name by trimming the count file's extension (".tsv")
samplenames <- gsub('.{4}$', '', count.files)

#set "factors" for the differential expression model
sampleCondition.Sex <- factor(c(rep("Male", 4), rep("Female", 3)))
sampleCondition.Population <- factor(c("MCsim_1", "MCsim_2", 
                                       "MCsim_3", "MCsim_4",
                                       "MCsim_2", 
                                       "MCsim_3", "MCsim_4"))

#construct condition table
sampleTable <- data.frame(sampleName = samplenames, fileName = count.files, 
                          sex = sampleCondition.Sex,
                          population = sampleCondition.Population)

#build DESeqDataSet from htseq count files
#omit genes where mean read counts over all the 8 samples is <50
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable, design= ~ sex + population)
keep.genes <- rowMeans(counts(ddsHTSeq)) >= 50
ddsHTSeq <- ddsHTSeq[keep.genes,]

#run DESeq2
dds.diffExp <- DESeq(ddsHTSeq)

#name the DESeq2 results
#set the "contrast" such that the log2 FoldChange in the output is Males/Females
results.diffExp = results(dds.diffExp, contrast=c("sex", "Male", "Female")) 

#extract the FlyBaseIDs and set it as the first column for ease
results.diffExp$FlyBaseID <- rownames(results.diffExp)
results.diffExp <- results.diffExp[,c(7,1:6)]

#save output to file
write.table(results.diffExp, paste0(path_save, "MCsim", "_", "body", "_Diff_Exp_Males_vs_Females.tsv"), 
            sep = '\t', quote = F, row.names = F, col.names = T)
