#This R script includes all the components needed to make Figures 1, S1 & S2
#The first section includes the code for how to make dotplots of Treatment Effect vs Sex Bias Category (Fig 1/S1)
#The second section includes code to perform one-sample permutation tests for whether the Treatment Effect 
#is significantly different from 0 (p-values from the second section were used to add the asterisks seen 
#atop the dotplots in Fig.1 and Fig. S1)
#The third section includes code to make Fig. S2: Treatment Effect vs Sex Bias (the latter as 
#a continuous variable), with a LOESS fitted line

library(plyr)
library(ggplot2)
library(gridExtra)
library(RColorBrewer)
library(ggpubr)
library(boot)


##############################################################################################
### SECTION 1 ################################################################################
##############################################################################################

#set path to load gonad specificity data
#this data will be used to make versions of figures/tables without genes with highly gonad-specific expression
path_1 <- 'C:/Users/mishr/OneDrive/Desktop/ProjectII/Annotated_R_scripts_and_data/Gonad.Specific.Expression/'
GSI.females <- read.table(paste0(path_1, "Gonad.Specificity.Index.Females.tsv"), header = T)
GSI.males <- read.table(paste0(path_1, "Gonad.Specificity.Index.Males.tsv"), header = T)

gonad.specific.genes <- union(GSI.females[which(GSI.females$Gonad.Specificity.Index >=0.95),]$geneID, 
                              GSI.males[which(GSI.males$Gonad.Specificity.Index >=0.95),]$geneID)

##########################################################################################################
#function to make dotplots of Treatment Effect vs Sex Bias Category for each pairwise comparion of   #####
#          mating treatments, for a given sex/tissue combination                                     #####
##########################################################################################################

#before defining function, set path for dataframes with condensed DESeq2 outputs (between mating trts) and
#sex bias data from Osada et al.'s data
path_2 <- 'C:/Users/mishr/OneDrive/Desktop/ProjectII/Annotated_R_scripts_and_data/Treatment.Effect.Dataframes/'

#function to yield dotplots
TreatmentEffect.Vs.SexBiasCategory.DotPlots <- function(sex, tissue){
  
  #load treatment effect dataframes
  sample <- na.omit(read.table(paste0(path_2, sex, ".", tissue, ".log2FCMatingTreatment.Dataframe.tsv"), header = T, sep = '\t'))
  
  #COMMENT THIS OUT TO RUN WITH GONAD-SPECIFIC GENES (like it is in Figure 1; to make figure S1 uncomment it)
  #sample <- sample[!sample$geneID %in% gonad.specific.genes,]
  
  #make smaller dataframes for each pairwise comparison of treatments, just including 
  #TreatmentEffect (log2FC Trt1/Trt2) and Sex Bias Category
  
  #dataframe MCabs vs MCsim
  unbiased.MCabs.Vs.MCsim <- cbind.data.frame(SexBias = rep("Unbiased", length(which(sample$sex.bias.category == "Unbiased"))),
                                    TreatmentEffect = sample[which(sample$sex.bias.category == "Unbiased"),]$log2FC.MCabs.vs.MCsim)
  fbiased.MCabs.Vs.MCsim <- cbind.data.frame(SexBias=rep("Female-biased", length(which(sample$sex.bias.category == "Female-biased"))),
                                   TreatmentEffect = sample[which(sample$sex.bias.category == "Female-biased"),]$log2FC.MCabs.vs.MCsim)
  mbiased.MCabs.Vs.MCsim <- cbind.data.frame(SexBias=rep("Male-biased", length(which(sample$sex.bias.category == "Male-biased"))),
                                   TreatmentEffect = sample[which(sample$sex.bias.category == "Male-biased"),]$log2FC.MCabs.vs.MCsim)
  df.MCabs.Vs.MCsim <- rbind.data.frame(unbiased.MCabs.Vs.MCsim, fbiased.MCabs.Vs.MCsim, mbiased.MCabs.Vs.MCsim)
  df.MCabs.Vs.MCsim$SexBias <- factor(df.MCabs.Vs.MCsim$SexBias, levels = c("Female-biased", "Unbiased", "Male-biased"))
  
  #dataframe for MCabs vs MCcom
  unbiased.MCabs.Vs.MCcom <- cbind.data.frame(SexBias = rep("Unbiased", length(which(sample$sex.bias.category == "Unbiased"))),
                                              TreatmentEffect = sample[which(sample$sex.bias.category == "Unbiased"),]$log2FC.MCabs.vs.MCcom)
  fbiased.MCabs.Vs.MCcom <- cbind.data.frame(SexBias=rep("Female-biased", length(which(sample$sex.bias.category == "Female-biased"))),
                                             TreatmentEffect = sample[which(sample$sex.bias.category == "Female-biased"),]$log2FC.MCabs.vs.MCcom)
  mbiased.MCabs.Vs.MCcom <- cbind.data.frame(SexBias=rep("Male-biased", length(which(sample$sex.bias.category == "Male-biased"))),
                                             TreatmentEffect = sample[which(sample$sex.bias.category == "Male-biased"),]$log2FC.MCabs.vs.MCcom)
  df.MCabs.Vs.MCcom <- rbind.data.frame(unbiased.MCabs.Vs.MCcom, fbiased.MCabs.Vs.MCcom, mbiased.MCabs.Vs.MCcom)
  df.MCabs.Vs.MCcom$SexBias <- factor(df.MCabs.Vs.MCcom$SexBias, levels = c("Female-biased", "Unbiased", "Male-biased"))
  
  #dataframe for MCsim vs MCcom
  unbiased.MCsim.Vs.MCcom <- cbind.data.frame(SexBias = rep("Unbiased", length(which(sample$sex.bias.category == "Unbiased"))),
                                              TreatmentEffect = sample[which(sample$sex.bias.category == "Unbiased"),]$log2FC.MCsim.vs.MCcom)
  fbiased.MCsim.Vs.MCcom <- cbind.data.frame(SexBias=rep("Female-biased", length(which(sample$sex.bias.category == "Female-biased"))),
                                             TreatmentEffect = sample[which(sample$sex.bias.category == "Female-biased"),]$log2FC.MCsim.vs.MCcom)
  mbiased.MCsim.Vs.MCcom <- cbind.data.frame(SexBias=rep("Male-biased", length(which(sample$sex.bias.category == "Male-biased"))),
                                             TreatmentEffect = sample[which(sample$sex.bias.category == "Male-biased"),]$log2FC.MCsim.vs.MCcom)
  df.MCsim.Vs.MCcom <- rbind.data.frame(unbiased.MCsim.Vs.MCcom, fbiased.MCsim.Vs.MCcom, mbiased.MCsim.Vs.MCcom)
  df.MCsim.Vs.MCcom$SexBias <- factor(df.MCsim.Vs.MCcom$SexBias, levels = c("Female-biased", "Unbiased", "Male-biased"))
  
  ##########################################################################
  #         summarise above dataframes with plyr
  #########################################################################
  df.MCabs.Vs.MCsim.summary <- ddply(df.MCabs.Vs.MCsim, c("SexBias"), summarise, N = length(TreatmentEffect), mean = mean(TreatmentEffect),
                                     sd = sd(TreatmentEffect), se = sd/sqrt(N))
  df.MCabs.Vs.MCcom.summary <- ddply(df.MCabs.Vs.MCcom, c("SexBias"), summarise, N = length(TreatmentEffect), mean = mean(TreatmentEffect),
                                     sd = sd(TreatmentEffect), se = sd/sqrt(N))
  df.MCsim.Vs.MCcom.summary <- ddply(df.MCsim.Vs.MCcom, c("SexBias"), summarise, N = length(TreatmentEffect), mean = mean(TreatmentEffect),
                                     sd = sd(TreatmentEffect), se = sd/sqrt(N))
  
  #make line plots for each df
  p1 <- ggplot(df.MCabs.Vs.MCsim.summary, aes(x = SexBias, y = mean, group = SexBias)) + 
    geom_point(size = 3)+
    geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width=.4,
                  position = position_dodge(0.05)) +
    theme_classic() +
    labs(x = "", y = "") +
    theme(legend.title = element_blank(),
          legend.text = element_text(color = "black", hjust = 0, size = 11.8),
          axis.text.x = element_text(size = 13, hjust = 0.5), 
          axis.text.y = element_text(size = 12, hjust = 1),
          axis.title.y = element_text(margin = margin(0,10,0,0), size = 14)) +
    scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 10)) +
    geom_abline(intercept = 0, slope = 0,  linewidth = 1, linetype="dashed", color = "grey48")
  p1
  
  p2 <- ggplot(df.MCabs.Vs.MCcom.summary, aes(x = SexBias, y = mean, group = SexBias)) + 
    geom_point(size = 3)+
    geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width=.4,
                  position = position_dodge(0.05)) +
    theme_classic() +
    labs(x = "", y = "") +
    theme(legend.title = element_blank(),
          legend.text = element_text(color = "black", hjust = 0, size = 11.8),
          axis.text.x = element_text(size = 13, hjust = 0.5), 
          axis.text.y = element_text(size = 12, hjust = 1),
          axis.title.y = element_text(margin = margin(0,10,0,0), size = 14)) +
    scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 10)) +
    geom_abline(intercept = 0, slope = 0,  linewidth = 1, linetype="dashed", color = "grey48") 
  p2
  
  p3 <- ggplot(df.MCsim.Vs.MCcom.summary, aes(x = SexBias, y = mean, group = SexBias)) + 
    geom_point(size = 3)+
    geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width=.4,
                  position = position_dodge(0.05)) +
    theme_classic() +
    labs(x = "", y = "") +
    theme(legend.title = element_blank(),
          legend.text = element_text(color = "black", hjust = 0, size = 11.8),
          axis.text.x = element_text(size = 13, hjust = 0.5), 
          axis.text.y = element_text(size = 12, hjust = 1),
          axis.title.y = element_text(margin = margin(0,10,0,0), size = 14)) +
    scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 10)) +
    geom_abline(intercept = 0, slope = 0,  linewidth = 1, linetype="dashed", color = "grey48") 
  p3
  
  #width x height = 1100 x 350
  grid.plot <- grid.arrange(p1, p2, p3, nrow=1)
  
  return(grid.plot)
}

#make dotplots; aggregate and edit in Powerpoint to make final figure
#NOTE - final figure has asterisks to denote which comparisons has a TreatmentEffect significantly different 
#from zero, with the help of one-sample permutation tests (refer to code in the section succeeding this)
DotPlot.female.body <- TreatmentEffect.Vs.SexBiasCategory.DotPlots("Female", "Body")
DotPlot.female.head <- TreatmentEffect.Vs.SexBiasCategory.DotPlots("Female", "Head")
DotPlot.male.body <- TreatmentEffect.Vs.SexBiasCategory.DotPlots("Male", "Body")
DotPlot.male.head <- TreatmentEffect.Vs.SexBiasCategory.DotPlots("Male", "Head")


##############################################################################################
### SECTION 2 ################################################################################
##############################################################################################
#set path to load gonad specificity data
path_1 <- 'C:/Users/mishr/OneDrive/Desktop/ProjectII/Annotated_R_scripts_and_data/Gonad.Specific.Expression/'
GSI.females <- read.table(paste0(path_1, "Gonad.Specificity.Index.Females.tsv"), header = T)
GSI.males <- read.table(paste0(path_1, "Gonad.Specificity.Index.Males.tsv"), header = T)

gonad.specific.genes <- union(GSI.females[which(GSI.females$Gonad.Specificity.Index >=0.95),]$geneID, 
                              GSI.males[which(GSI.males$Gonad.Specificity.Index >=0.95),]$geneID)


#########################################################################################
## functions to perform one-sample permutation tests of treatment effect for a ##########
## given sample (e.g., male head)                                             ###########
#########################################################################################
#before defining function, set path for folder with treatment dataframes 
path_2 <- 'C:/Users/mishr/OneDrive/Desktop/ProjectII/Annotated_R_scripts_and_data/Treatment.Effect.Dataframes/'

#define function to obtain mean treatment effect and bootstrap 95% Conf. Intervals from 
#an input list of trtEffects (function returns a one-row dataframe)
Bootstrap.CIs <- function(list.of.Trt.Effects){
  mean.TreatmentEffect <- mean(list.of.Trt.Effects)
  meanfun <- function(data, i){
    d <- data[i]
    return(mean(d))   
  }
  
  #obtain bootstrap confidence intervals
  bootstrap.mean <- boot(list.of.Trt.Effects, statistic=meanfun, R=10000)
  bootstrap <- boot.ci(bootstrap.mean, conf=0.95, type="bca")
  bootstrap.confInts <- bootstrap$bca[4:5]
  
  #make a dataframe with mean Treatment Effect and its CIs
  result <- cbind.data.frame(mean.TreatmentEffect = mean.TreatmentEffect, Bootstrap.Lower.Conf.Int = bootstrap.confInts[1],
                             Bootstrap.Upper.Conf.Int = bootstrap.confInts[2])
  
  return(result)
}

#define function to perform one-sample permutation tests, to test whether the Treatment Effect (for a given
#sex bias category in a pairwise treatment comparison) is significantly greater than zero. Function yields
#a p-value.
PermutationTest.of.TreatmentEffect <- function(list.of.Trt.Effects){
  mean.TreatmentEffect <- mean(list.of.Trt.Effects)
  
  #make empty vector to store permuted values of treatment effects
  permuted.TrtEffects <- c()
  
  #make empty vector to store test statistic for each iteration
  test.statistics <- c()
  
  #for loop to perform 10000 permutations
  for (i in 1:10000){
    permuted.TrtEffects <- sample(list.of.Trt.Effects, length(list.of.Trt.Effects), replace = F)
    permuted.TrtEffects <- c(-1*permuted.TrtEffects[1:(length(list.of.Trt.Effects)/2)],
                             permuted.TrtEffects[(length(list.of.Trt.Effects)/2 + 1):(length(list.of.Trt.Effects))])
    mean.TrtEffect.permuted <- mean(permuted.TrtEffects)
    test.statistics <- c(test.statistics, mean.TrtEffect.permuted)
  }
  
  pvalue <- sum(abs(test.statistics) > abs(mean.TreatmentEffect))/10000
  
  return(pvalue)
}

#Use the previous two functions on each sex bias category and pairwise treatment comparison within a given 
#sample (e.g., female head)
#NOTE - takes a few minutes to run because bootstrapping takes time
Perform.Bootstrapping.and.PermTest <- function(sex, tissue){
  #load treatment effect dataframes
  sample <- na.omit(read.table(paste0(path_2, sex, ".", tissue, ".log2FCMatingTreatment.Dataframe.tsv"), header = T, sep = '\t'))
  
  #COMMENT THIS OUT TO RUN FUNCTIONS WITH GONAD-SPECIFIC GENES
  #sample <- sample[!sample$geneID %in% gonad.specific.genes,]
  
  #make separate dataframes for each sex bias category
  unbiased <- sample[which(sample$sex.bias.category == "Unbiased"),]
  fbiased <- sample[which(sample$sex.bias.category == "Female-biased"),]
  mbiased <- sample[which(sample$sex.bias.category == "Male-biased"),]
  
  #run the bootstraping function on each pairwise treatment comparison for all the sex bias categories, the
  #combine (using rbind) into a simple dataframe
  unbiased.MCabs.vs.MCsim.boot <- Bootstrap.CIs(unbiased$log2FC.MCabs.vs.MCsim)
  fbiased.MCabs.vs.MCsim.boot <- Bootstrap.CIs(fbiased$log2FC.MCabs.vs.MCsim)
  mbiased.MCabs.vs.MCsim.boot <- Bootstrap.CIs(mbiased$log2FC.MCabs.vs.MCsim)
  
  unbiased.MCabs.vs.MCcom.boot <- Bootstrap.CIs(unbiased$log2FC.MCabs.vs.MCcom)
  fbiased.MCabs.vs.MCcom.boot <- Bootstrap.CIs(fbiased$log2FC.MCabs.vs.MCcom)
  mbiased.MCabs.vs.MCcom.boot <- Bootstrap.CIs(mbiased$log2FC.MCabs.vs.MCcom)
  
  unbiased.MCsim.vs.MCcom.boot <- Bootstrap.CIs(unbiased$log2FC.MCsim.vs.MCcom)
  fbiased.MCsim.vs.MCcom.boot <- Bootstrap.CIs(fbiased$log2FC.MCsim.vs.MCcom)
  mbiased.MCsim.vs.MCcom.boot <- Bootstrap.CIs(mbiased$log2FC.MCsim.vs.MCcom)
  
  df.boot <- rbind.data.frame(unbiased.MCabs.vs.MCsim.boot, fbiased.MCabs.vs.MCsim.boot,
                              mbiased.MCabs.vs.MCsim.boot, unbiased.MCabs.vs.MCcom.boot,
                              fbiased.MCabs.vs.MCcom.boot, mbiased.MCabs.vs.MCcom.boot,
                              unbiased.MCsim.vs.MCcom.boot, fbiased.MCsim.vs.MCcom.boot,
                              mbiased.MCsim.vs.MCcom.boot)
  
  #run the permutation test function on each pairwise treatment comparison for all the 
  #sex bias categories, then combined the p-values into one vector
  unbiased.MCabs.vs.MCsim.PermTest <- PermutationTest.of.TreatmentEffect(unbiased$log2FC.MCabs.vs.MCsim)
  fbiased.MCabs.vs.MCsim.PermTest <- PermutationTest.of.TreatmentEffect(fbiased$log2FC.MCabs.vs.MCsim)
  mbiased.MCabs.vs.MCsim.PermTest <- PermutationTest.of.TreatmentEffect(mbiased$log2FC.MCabs.vs.MCsim)
  
  unbiased.MCabs.vs.MCcom.PermTest <- PermutationTest.of.TreatmentEffect(unbiased$log2FC.MCabs.vs.MCcom)
  fbiased.MCabs.vs.MCcom.PermTest <- PermutationTest.of.TreatmentEffect(fbiased$log2FC.MCabs.vs.MCcom)
  mbiased.MCabs.vs.MCcom.PermTest <- PermutationTest.of.TreatmentEffect(mbiased$log2FC.MCabs.vs.MCcom)
  
  unbiased.MCsim.vs.MCcom.PermTest <- PermutationTest.of.TreatmentEffect(unbiased$log2FC.MCsim.vs.MCcom)
  fbiased.MCsim.vs.MCcom.PermTest <- PermutationTest.of.TreatmentEffect(fbiased$log2FC.MCsim.vs.MCcom)
  mbiased.MCsim.vs.MCcom.PermTest <- PermutationTest.of.TreatmentEffect(mbiased$log2FC.MCsim.vs.MCcom)
  
  all.pvalues.PermTest <- c(unbiased.MCabs.vs.MCsim.PermTest, fbiased.MCabs.vs.MCsim.PermTest,
                            mbiased.MCabs.vs.MCsim.PermTest, unbiased.MCabs.vs.MCcom.PermTest,
                            fbiased.MCabs.vs.MCcom.PermTest, mbiased.MCabs.vs.MCcom.PermTest,
                            unbiased.MCsim.vs.MCcom.PermTest, fbiased.MCsim.vs.MCcom.PermTest,
                            mbiased.MCsim.vs.MCcom.PermTest)
  
  #make a dataframe with the outputs from above
  sample_name <- paste0(Sex, "_", Tissue)
  sex_bias_categories <- rep(c("Unbiased", "Female-biased", "Male-biased"), 3)
  treatment_comparison <- c(rep("MCabs_vs_MCsim", 3), rep("MCabs_vs_MCcom", 3), rep("MCsim_vs_MCcom", 3))
  
  df <- cbind.data.frame(SampleName = rep(sample_name, 9), Sex.Bias.Category = sex_bias_categories,
                         Comparison = treatment_comparison, df.boot, 
                         pVal.OneSamplePermutationTest = all.pvalues.PermTest)
  
  return(df)
}

Perform.Bootstrapping.and.PermTest("Female", "Body")
Perform.Bootstrapping.and.PermTest("Male", "Body")
Perform.Bootstrapping.and.PermTest("Female", "Head")
Perform.Bootstrapping.and.PermTest("Male", "Head")


##############################################################################################
### SECTION 3 ################################################################################
##############################################################################################

#set path to load gonad specificity data
path_1 <- 'C:/Users/mishr/OneDrive/Desktop/ProjectII/Annotated_R_scripts_and_data/Gonad.Specific.Expression/'
GSI.females <- read.table(paste0(path_1, "Gonad.Specificity.Index.Females.tsv"), header = T)
GSI.males <- read.table(paste0(path_1, "Gonad.Specificity.Index.Males.tsv"), header = T)

gonad.specific.genes <- union(GSI.females[which(GSI.females$Gonad.Specificity.Index >=0.95),]$geneID, 
                              GSI.males[which(GSI.males$Gonad.Specificity.Index >=0.95),]$geneID)


##########################################################################################################
#function to make dotplots of Treatment Effect vs Sex Bias Category for each pairwise comparion of   #####
#          mating treatments, for a given sex/tissue combination                                     #####
##########################################################################################################

#before defining function, set path for dataframes with condensed DESeq2 outputs (between mating trts) and
#sex bias data from Osada et al.'s data
path_2 <- 'C:/Users/mishr/OneDrive/Desktop/ProjectII/Annotated_R_scripts_and_data/Treatment.Effect.Dataframes/'

#function to yield dotplots
#comment out the relevant lines to run without gonad-specific genes
#change to xlim = c(-2,2) for head samples, as well as coordinates of labels 'FB' and 'MB'
TreatmentEffect.Vs.SexBias.LOESS.fit <- function(sex, tissue){
  
  #load treatment effect dataframes
  sample <- na.omit(read.table(paste0(path_2, sex, ".", tissue, ".log2FCMatingTreatment.Dataframe.tsv"), header = T, sep = '\t'))
  
  #COMMENT THIS OUT TO RUN WITH GONAD-SPECIFIC GENES (like it is in Figure 1; to make figure S1 uncomment it)
  #sample <- sample[!sample$geneID %in% gonad.specific.genes,]
  
  #make scatterplots
  #specify "loess" as smoothing method with the window for local regressions ("span") set to 0.5 & 
  #conf. interval set to 0.95
  p1 <- ggplot(sample, aes(x = log2FC.sex.Osada, y = log2FC.MCabs.vs.MCsim)) + 
    theme_classic() +
    geom_point(size = 0.5, colour = "grey78") + 
    geom_smooth(method = "loess", lwd = 0.7, se = T, level = 0.95, span = 0.5, colour = "blue", fill = "green") +
    geom_hline(yintercept = 0, colour = "black", lwd = 0.5, linetype = "dashed") +
    geom_vline(xintercept = -0.5, colour = "black", lwd = 0.5, linetype = "dashed") +
    geom_vline(xintercept = 0.5, colour = "black", lwd = 0.5, linetype = "dashed") +
    annotate(geom="text", x=-4.5, y=0.9, label="FB", color="black", fontface="bold") +
    annotate(geom="text", x=4.5, y=0.9, label="MB", color="black", fontface="bold") +
    labs(x = "", y = "Treatment Effect", title = "") +
    theme(legend.title = element_blank(),
          plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(size = 12, hjust = 0.5), 
          axis.text.y = element_text(size = 12, hjust = 1),
          axis.title.x = element_text(margin = margin(0,10,0,0), size = 14),
          axis.title.y = element_text(margin = margin(0,10,0,0), size = 14)) + coord_cartesian(ylim = c(-1,1))
  #coord_cartesian(xlim = c(-2,2), ylim = c(-1,1))
  
  p2 <- ggplot(sample, aes(x = log2FC.sex.Osada, y = log2FC.MCabs.vs.MCcom)) + 
    theme_classic() +
    geom_point(size = 0.5, colour = "grey78") + 
    geom_smooth(method = "loess", lwd = 0.7, se = T, level = 0.95, span = 0.5, colour = "blue", fill = "green") +
    geom_hline(yintercept = 0, colour = "black", lwd = 0.5, linetype = "dashed") +
    geom_vline(xintercept = -0.5, colour = "black", lwd = 0.5, linetype = "dashed") +
    geom_vline(xintercept = 0.5, colour = "black", lwd = 0.5, linetype = "dashed") +
    annotate(geom="text", x=-4.5, y=0.9, label="FB", color="black", fontface="bold") +
    annotate(geom="text", x=4.5, y=0.9, label="MB", color="black", fontface="bold") +
    labs(x = "Sex Bias", y = "", title = "") +
    theme(legend.title = element_blank(),
          plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(size = 12, hjust = 0.5), 
          axis.text.y = element_text(size = 12, hjust = 1),
          axis.title.x = element_text(margin = margin(0,10,0,0), size = 14),
          axis.title.y = element_text(margin = margin(0,10,0,0), size = 14))  + coord_cartesian(ylim = c(-1,1))
  #coord_cartesian(xlim = c(-2,2), ylim = c(-1,1))
  
  p3 <- ggplot(sample, aes(x = log2FC.sex.Osada, y = log2FC.MCsim.vs.MCcom)) + 
    theme_classic() +
    geom_point(size = 0.5, colour = "grey78") + 
    geom_smooth(method = "loess", lwd = 0.7, se = T, level = 0.95, span = 0.5, colour = "blue", fill = "green") +
    geom_hline(yintercept = 0, colour = "black", lwd = 0.5, linetype = "dashed") +
    geom_vline(xintercept = -0.5, colour = "black", lwd = 0.5, linetype = "dashed") +
    geom_vline(xintercept = 0.5, colour = "black", lwd = 0.5, linetype = "dashed") +
    annotate(geom="text", x=-4.5, y=0.9, label="FB", color="black", fontface="bold") +
    annotate(geom="text", x=4.5, y=0.9, label="MB", color="black", fontface="bold") +
    labs(x = "", y = "", title = "") +
    theme(legend.title = element_blank(),
          plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(size = 12, hjust = 0.5), 
          axis.text.y = element_text(size = 12, hjust = 1),
          axis.title.x = element_text(margin = margin(0,10,0,0), size = 14),
          axis.title.y = element_text(margin = margin(0,10,0,0), size = 14))  + coord_cartesian(ylim = c(-1,1))
  #coord_cartesian(xlim = c(-2,2), ylim = c(-1,1))
  
  p <- grid.arrange(p1, p2, p3, nrow = 1)
  
  return(p)  
  
}

p1 <- TreatmentEffect.Vs.SexBias.LOESS.fit("Female", "Body")  
p2 <- TreatmentEffect.Vs.SexBias.LOESS.fit("Male", "Body")  
p3 <- TreatmentEffect.Vs.SexBias.LOESS.fit("Female", "Head")  
p4 <- TreatmentEffect.Vs.SexBias.LOESS.fit("Male", "Head")  

#save as TIFF file
grid.arrange(p1, p2, p3, p4, nrow = 4)

