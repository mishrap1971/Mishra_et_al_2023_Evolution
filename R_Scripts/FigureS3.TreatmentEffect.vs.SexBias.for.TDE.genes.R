#This script was used to make the scatterplots in Fig S3
#Note: Only TDE genes are depicted, and the SBGE data (external, from Osada et al 2017) is slightly curtailed
#to aid visualisation

library(ggplot2)
library(gridExtra)

#specify path for: a) where the external Sex-Biased Gene Expression data is stored, and
#b) where Treatment Effect dataframes are stored
path_for_SBGE_data <- "C:/Users/mishr/OneDrive/Desktop/ProjectII/Annotated_R_scripts_and_data/Sex.Biased.Gene.Expression.Osada_et_al_2017/"
path_trtEffect_data <- 'C:/Users/mishr/OneDrive/Desktop/ProjectII/Annotated_R_scripts_and_data/Treatment.Effect.Dataframes/'

#load sex bias data for head and whole body, estimated from data in Osada et al 2017
sexbias.body <- read.table(paste0(path_for_SBGE_data, "Sex.Biased.Genes.Body.tsv"), header = T)
sexbias.head <- read.table(paste0(path_for_SBGE_data, "Sex.Biased.Genes.Head.tsv"), header = T)

#curtail the sex bias values in each tissue, because we don't want final plots to have sex bias huge outliers
#for body: remove genes whose sex bias falls outside the 95% Conf Int (lower CI = -2.840484, upper CI = 11.723497)
#for head: only retain genes where -2<log2FC<2
sexbias.body <- sexbias.body[which(sexbias.body$log2FoldChange > -2.840484 & sexbias.body$log2FoldChange < 11.723497),]
sexbias.head <- sexbias.head[which(sexbias.head$log2FoldChange > -2 & sexbias.head$log2FoldChange < 2),]

#load treatment effect dataframes for all samples
female.body <- read.table(paste0(path_trtEffect_data, "Female.Body.log2FCMatingTreatment.Dataframe.tsv"), header = T, sep ='\t')
female.head <- read.table(paste0(path_trtEffect_data, "Female.Head.log2FCMatingTreatment.Dataframe.tsv"), header = T, sep ='\t')
male.body <- read.table(paste0(path_trtEffect_data, "Male.Body.log2FCMatingTreatment.Dataframe.tsv"), header = T, sep ='\t')
male.head <- read.table(paste0(path_trtEffect_data, "Male.Head.log2FCMatingTreatment.Dataframe.tsv"), header = T, sep ='\t')

#retain genes for which we have sex bias data as well as data in Treatment Effect dataframes 
female.body <- female.body[female.body$geneID %in% sexbias.body$geneID,]
male.body <- male.body[male.body$geneID %in% sexbias.body$geneID,]
female.head <- female.head[female.head$geneID %in% sexbias.head$geneID,]
male.head <- male.head[male.head$geneID %in% sexbias.head$geneID,]

#############################################################################################
## make scatterplots of Treatment Effect vs Sex Bias (log2FC M/F) for TDE genes #############
#############################################################################################

#################################################
############# female body  ######################
#################################################
sexbias.female.body <- sexbias.body[sexbias.body$geneID %in% female.body$geneID,]
female.body$log2FC.Sex.Osada <- sexbias.female.body$log2FoldChange
sig.MCabs.vs.MCsim <- female.body[which(female.body$padj.1 < 0.1),][c(1,2,11)]
sig.MCabs.vs.MCcom <- female.body[which(female.body$padj.2 < 0.1),][c(1,5,11)]
sig.MCsim.vs.MCcom <- female.body[which(female.body$padj.3 < 0.1),][c(1,8,11)]

p1 <- ggplot(sig.MCabs.vs.MCsim, aes(log2FC.sex.Osada, log2FC.MCabs.vs.MCsim)) +
  geom_point(size = 1, show.legend = FALSE) +
  labs(x = '', 
       y = 'Treatment Effect (Trt 1 relative to Trt 2)',  
       tag = "A") +
  theme_classic() +
  theme(axis.text.x = element_text(margin = margin(5,0,0,0), color = "black", size = 11),
        axis.text.y = element_text(margin = margin(0,5,0,0), color = "black", size = 11), 
        axis.title.x = element_text(margin = margin(10,0,0,0), size = 14),
        axis.title.y = element_text(margin = margin(0,10,0,0), size = 14), 
        plot.tag = element_text(face = "bold", size = 14),
        plot.margin = margin(6,6,0,6)) + 
  geom_hline(yintercept = 0, color = "black", lwd = 1) +
  geom_vline(xintercept = 0, color = "black", lwd = 1) + ylim(-3,4)
p1

p2 <- ggplot(sig.MCabs.vs.MCcom, aes(log2FC.sex.Osada, log2FC.MCabs.vs.MCcom)) +
  geom_point(size = 1, show.legend = FALSE) +
  labs(x = 'Sex Bias', 
       y = '',  
       tag = "B") +
  theme_classic() +
  theme(axis.text.x = element_text(margin = margin(5,0,0,0), color = "black", size = 11),
        axis.text.y = element_text(margin = margin(0,5,0,0), color = "black", size = 11), 
        axis.title.x = element_text(margin = margin(10,0,0,0), size = 14),
        axis.title.y = element_text(margin = margin(0,10,0,0), size = 14), 
        plot.tag = element_text(face = "bold", size = 14),
        plot.margin = margin(6,6,0,6)) + 
  geom_hline(yintercept = 0, color = "black", lwd = 1) +
  geom_vline(xintercept = 0, color = "black", lwd = 1) + ylim(-3,4)
p2

p3 <- ggplot(sig.MCsim.vs.MCcom, aes(log2FC.sex.Osada, log2FC.MCsim.vs.MCcom)) +
  geom_point(size = 1, show.legend = FALSE) +
  labs(x = '', 
       y = '',  
       tag = "C") +
  theme_classic() +
  theme(axis.text.x = element_text(margin = margin(5,0,0,0), color = "black", size = 11),
        axis.text.y = element_text(margin = margin(0,5,0,0), color = "black", size = 11), 
        axis.title.x = element_text(margin = margin(10,0,0,0), size = 14),
        axis.title.y = element_text(margin = margin(0,10,0,0), size = 14), 
        plot.tag = element_text(face = "bold", size = 14),
        plot.margin = margin(6,6,0,6)) + 
  geom_hline(yintercept = 0, color = "black", lwd = 1) +
  geom_vline(xintercept = 0, color = "black", lwd = 1) + ylim(-3,4)
p3

#################################################
############# male body  ######################
#################################################
sexbias.male.body <- sexbias.body[sexbias.body$geneID %in% male.body$geneID,]
male.body$log2FC.Sex.Osada <- sexbias.male.body$log2FoldChange
sig.MCabs.vs.MCsim <- male.body[which(male.body$padj.1 < 0.1),][c(1,2,11)]
sig.MCabs.vs.MCcom <- male.body[which(male.body$padj.2 < 0.1),][c(1,5,11)]
sig.MCsim.vs.MCcom <- male.body[which(male.body$padj.3 < 0.1),][c(1,8,11)]

p4 <- ggplot(sig.MCabs.vs.MCsim, aes(log2FC.sex.Osada, log2FC.MCabs.vs.MCsim)) +
  geom_point(size = 1, show.legend = FALSE) +
  labs(x = '', 
       y = '',  
       tag = "D") +
  theme_classic() +
  theme(axis.text.x = element_text(margin = margin(5,0,0,0), color = "black", size = 11),
        axis.text.y = element_text(margin = margin(0,5,0,0), color = "black", size = 11), 
        axis.title.x = element_text(margin = margin(10,0,0,0), size = 14),
        axis.title.y = element_text(margin = margin(0,10,0,0), size = 14), 
        plot.tag = element_text(face = "bold", size = 14),
        plot.margin = margin(6,6,0,6)) + 
  geom_hline(yintercept = 0, color = "black", lwd = 1) +
  geom_vline(xintercept = 0, color = "black", lwd = 1) + ylim(-2,3)
p4

p5 <- ggplot(sig.MCabs.vs.MCcom, aes(log2FC.sex.Osada, log2FC.MCabs.vs.MCcom)) +
  geom_point(size = 1, show.legend = FALSE) +
  labs(x = '', 
       y = '',  
       tag = "E") +
  theme_classic() +
  theme(axis.text.x = element_text(margin = margin(5,0,0,0), color = "black", size = 11),
        axis.text.y = element_text(margin = margin(0,5,0,0), color = "black", size = 11), 
        axis.title.x = element_text(margin = margin(10,0,0,0), size = 14),
        axis.title.y = element_text(margin = margin(0,10,0,0), size = 14), 
        plot.tag = element_text(face = "bold", size = 14),
        plot.margin = margin(6,6,0,6)) + 
  geom_hline(yintercept = 0, color = "black", lwd = 1) +
  geom_vline(xintercept = 0, color = "black", lwd = 1) + ylim(-2,3)
p5

p6 <- ggplot(sig.MCsim.vs.MCcom, aes(log2FC.sex.Osada, log2FC.MCsim.vs.MCcom)) +
  geom_point(size = 1, show.legend = FALSE) +
  labs(x = '', 
       y = '',  
       tag = "F") +
  theme_classic() +
  theme(axis.text.x = element_text(margin = margin(5,0,0,0), color = "black", size = 11),
        axis.text.y = element_text(margin = margin(0,5,0,0), color = "black", size = 11), 
        axis.title.x = element_text(margin = margin(10,0,0,0), size = 14),
        axis.title.y = element_text(margin = margin(0,10,0,0), size = 14), 
        plot.tag = element_text(face = "bold", size = 14),
        plot.margin = margin(6,6,0,6)) + 
  geom_hline(yintercept = 0, color = "black", lwd = 1) +
  geom_vline(xintercept = 0, color = "black", lwd = 1) + ylim(-2,3)
p6

#################################################
############# female head  ######################
#################################################
sexbias.female.head <- sexbias.head[sexbias.head$geneID %in% female.head$geneID,]
female.head$log2FC.Sex.Osada <- sexbias.female.head$log2FoldChange
sig.MCabs.vs.MCsim <- female.head[which(female.head$padj.1 < 0.1),][c(1,2,11)]
sig.MCabs.vs.MCcom <- female.head[which(female.head$padj.2 < 0.1),][c(1,5,11)]
sig.MCsim.vs.MCcom <- female.head[which(female.head$padj.3 < 0.1),][c(1,8,11)]

p7 <- ggplot(sig.MCabs.vs.MCsim, aes(log2FC.sex.Osada, log2FC.MCabs.vs.MCsim)) +
  geom_point(size = 1, show.legend = FALSE) +
  labs(x = '', 
       y = '',  
       tag = "G") +
  theme_classic() +
  theme(axis.text.x = element_text(margin = margin(5,0,0,0), color = "black", size = 11),
        axis.text.y = element_text(margin = margin(0,5,0,0), color = "black", size = 11), 
        axis.title.x = element_text(margin = margin(10,0,0,0), size = 14),
        axis.title.y = element_text(margin = margin(0,10,0,0), size = 14), 
        plot.tag = element_text(face = "bold", size = 14),
        plot.margin = margin(6,6,0,6)) + 
  geom_hline(yintercept = 0, color = "black", lwd = 1) +
  geom_vline(xintercept = 0, color = "black", lwd = 1) + ylim(-2,4)
p7

p8 <- ggplot(sig.MCabs.vs.MCcom, aes(log2FC.sex.Osada, log2FC.MCabs.vs.MCcom)) +
  geom_point(size = 1, show.legend = FALSE) +
  labs(x = '', 
       y = '',  
       tag = "H") +
  theme_classic() +
  theme(axis.text.x = element_text(margin = margin(5,0,0,0), color = "black", size = 11),
        axis.text.y = element_text(margin = margin(0,5,0,0), color = "black", size = 11), 
        axis.title.x = element_text(margin = margin(10,0,0,0), size = 14),
        axis.title.y = element_text(margin = margin(0,10,0,0), size = 14), 
        plot.tag = element_text(face = "bold", size = 14),
        plot.margin = margin(6,6,0,6)) + 
  geom_hline(yintercept = 0, color = "black", lwd = 1) +
  geom_vline(xintercept = 0, color = "black", lwd = 1) + ylim(-2,4)
p8

p9 <- ggplot(sig.MCsim.vs.MCcom, aes(log2FC.sex.Osada, log2FC.MCsim.vs.MCcom)) +
  geom_point(size = 1, show.legend = FALSE) +
  labs(x = '', 
       y = '',  
       tag = "I") +
  theme_classic() +
  theme(axis.text.x = element_text(margin = margin(5,0,0,0), color = "black", size = 11),
        axis.text.y = element_text(margin = margin(0,5,0,0), color = "black", size = 11), 
        axis.title.x = element_text(margin = margin(10,0,0,0), size = 14),
        axis.title.y = element_text(margin = margin(0,10,0,0), size = 14), 
        plot.tag = element_text(face = "bold", size = 14),
        plot.margin = margin(6,6,0,6)) + 
  geom_hline(yintercept = 0, color = "black", lwd = 1) +
  geom_vline(xintercept = 0, color = "black", lwd = 1) + ylim(-2,4)
p9


#################################################
############# male head  ######################
#################################################

sexbias.male.head <- sexbias.head[sexbias.head$geneID %in% male.head$geneID,]
male.head$log2FC.Sex.Osada <- sexbias.male.head$log2FoldChange
sig.MCabs.vs.MCsim <- male.head[which(male.head$padj.1 < 0.1),][c(1,2,11)]
sig.MCabs.vs.MCcom <- male.head[which(male.head$padj.2 < 0.1),][c(1,5,11)]
sig.MCsim.vs.MCcom <- male.head[which(male.head$padj.3 < 0.1),][c(1,8,11)]

p10 <- ggplot(sig.MCabs.vs.MCsim, aes(log2FC.sex.Osada, log2FC.MCabs.vs.MCsim)) +
  geom_point(size = 1, show.legend = FALSE) +
  labs(x = '', 
       y = '',  
       tag = "J") +
  theme_classic() +
  theme(axis.text.x = element_text(margin = margin(5,0,0,0), color = "black", size = 11),
        axis.text.y = element_text(margin = margin(0,5,0,0), color = "black", size = 11), 
        axis.title.x = element_text(margin = margin(10,0,0,0), size = 14),
        axis.title.y = element_text(margin = margin(0,10,0,0), size = 14), 
        plot.tag = element_text(face = "bold", size = 14),
        plot.margin = margin(6,6,0,6)) + 
  geom_hline(yintercept = 0, color = "black", lwd = 1) +
  geom_vline(xintercept = 0, color = "black", lwd = 1) + ylim(-2,3)
p10

p11 <- ggplot(sig.MCabs.vs.MCcom, aes(log2FC.sex.Osada, log2FC.MCabs.vs.MCcom)) +
  geom_point(size = 1, show.legend = FALSE) +
  labs(x = '', 
       y = '',  
       tag = "K") +
  theme_classic() +
  theme(axis.text.x = element_text(margin = margin(5,0,0,0), color = "black", size = 11),
        axis.text.y = element_text(margin = margin(0,5,0,0), color = "black", size = 11), 
        axis.title.x = element_text(margin = margin(10,0,0,0), size = 14),
        axis.title.y = element_text(margin = margin(0,10,0,0), size = 14), 
        plot.tag = element_text(face = "bold", size = 14),
        plot.margin = margin(6,6,0,6)) + 
  geom_hline(yintercept = 0, color = "black", lwd = 1) +
  geom_vline(xintercept = 0, color = "black", lwd = 1) + ylim(-2,3)
p11

p12 <- ggplot(sig.MCsim.vs.MCcom, aes(log2FC.sex.Osada, log2FC.MCsim.vs.MCcom)) +
  geom_point(size = 1, show.legend = FALSE) +
  labs(x = '', 
       y = '',  
       tag = "L") +
  theme_classic() +
  theme(axis.text.x = element_text(margin = margin(5,0,0,0), color = "black", size = 11),
        axis.text.y = element_text(margin = margin(0,5,0,0), color = "black", size = 11), 
        axis.title.x = element_text(margin = margin(10,0,0,0), size = 14),
        axis.title.y = element_text(margin = margin(0,10,0,0), size = 14), 
        plot.tag = element_text(face = "bold", size = 14),
        plot.margin = margin(6,6,0,6)) + 
  geom_hline(yintercept = 0, color = "black", lwd = 1) +
  geom_vline(xintercept = 0, color = "black", lwd = 1) + ylim(-2,3)
p12

#save with dimensions 980 x 680
grid.arrange(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, nrow=4)





