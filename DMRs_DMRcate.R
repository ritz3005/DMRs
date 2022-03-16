#####Part 1#####################################################################
options(stringsAsFactors = FALSE)
rm(list=ls())

#packages
library(openxlsx)
library(dplyr)
library(tibble)
library(ggplot2)
#library(tidyverse)
library(parallel)
library(readxl)

library(openxlsx)
library(affycoretools)
library(BiocStyle)
library(xtable)
library(minfi)
library(IlluminaHumanMethylationEPICmanifest)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
data(SNPs.147CommonSingle)
library(sva)
library(limma)
library(Glimma)
library(RColorBrewer)
library(DMRcate)
library(Homo.sapiens)
library(missMethyl)
library(DT)
library(Gviz)
#install.packages("Gviz")

#source("https://bioconductor.org/biocLite.R")
#biocLite("AnnotationHub")
#biocLite("shiny")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#install.packages("BiocManager")

#BiocManager::install("AnnotationHub")
#BiocManager::install("Glimma")
#BiocManager::install("DMRcate")
#BiocManager::install("shiny")

##paths###

idatdir <-  "/mnt/nas/global/ame/users/Ritika/Projects/FELICITy_Project/EPIC_DATA_norm/EPIC data/infos/ScanData/"
#list.files(idatdir)
phenodir <-  "/mnt/nas/global/ame/users/Ritika/Projects/FELICITy_Project/FELICITy_data_anal_Jim/Updated/"

rdir <- "/mnt/nas/global/ame/users/Ritika/Projects/FELICITy_Project/EPIC_DATA_norm/Results_models_FELICITy/"
list.files(phenodir)

basedir <- "/mnt/nas/global/ame/users/Ritika/Projects/FELICITy_Project/EPIC_DATA_norm/EPIC data/"


samps <- read.xlsx(paste0(phenodir, "EPIC_plate_FELICITy_21090919.xlsx"))
samps$Patient.Code <- fixIt(samps$Patient.Code)
#samplesheet
targets <- read.metharray.sheet(idatdir, "M00936_Pl1_3_n114_Felicity.csv") 


targets$Sample_Name <- fixIt(targets$Sample_Name)
samps <- samps[match(targets$Sample_Name, samps$Patient.Code),]

##Phenotype data
phenos <- read.xlsx(paste0(phenodir, "Felicity_phenotypes_114.xlsx"))

phenos <- phenos[match(targets$Sample_Name, phenos$Code),]
phenos$Group <- factor(phenos$Group)

samps <- cbind(samps[,-c(11,13)], phenos[,-1])
samps$Stress.group <- factor(samps$Stress.group)
samps$Illumina.Plate <- factor(samps$Illumina.Plate)
samps$AChE.BChE <- samps$fet_AChE/samps$fet_BChE
samps$Stress.group <- relevel(samps$Stress.group, "2")
samps$Gender <- factor(samps$Gender)
samps$Smoking <- factor(samps$Smoking)

##loading normalised data
load(paste0(idatdir, "methdataall.Rdata"))
dat 

annoEPIC = getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
head(annoEPIC)

#matching cpg site ronames of mvals and annoepic
annoEPICSub <- annoEPIC[match(rownames(getM(eset.nosex)),annoEPIC$Name),                    
                        c(1:4,12:19,22:ncol(annoEPIC))]                       


################Stress:- adjusted##########################################################################################################

##basic variable checks##############
sum(is.na(samps$stress))
#[1] 7


sum(is.na(samps$Autoimmunediseases))#0
sum(is.na(samps$Gestational_Diabetes)) #0

table(samps$Smoking)
table(samps$Autoimmunediseases)
table(samps$Gestational_Diabetes)

dim(model.matrix(~stress+ Gender + Smoking + Autoimmunediseases +  Gestational_age_birth + Illumina.Plate, samps))
#[1] 107  8

###designing the model matrix using limma - important for DMP and DMR analysis#########
design <- model.matrix(~ stress  + Gender + Smoking + Autoimmunediseases + Gestational_Diabetes + Gestational_age_birth +Illumina.Plate, samps)

sv <- sva(getM(eset.nosex)[,!is.na(samps$stress)], design,vfilter = 5e4) #3

design <- cbind(design, sv$sv)


###DMP analysis using M values#####ignore###
fit <- lmFit(getM(eset.nosex)[,!is.na(samps$stress)], design)
fit2 <- eBayes(fit)
dim(fit2)
colnames(fit2)
topTable(fit2, 2)
DMPs <- topTable(fit2, coef=2,confint = TRUE, genelist=annoEPICSub, number= 808554) 

bonfc <- 0.05/nrow(fit2)
bonfc
#[1] 6.183879e-08
sum(DMPs$P.Value < bonfc) #1
sum(DMPs$adj.P.Val < 0.05) # 4
DMPs10 <- topTable(fit2, 2,confint = TRUE, genelist = annoEPICSub)

save(DMPs, samps,  file = paste0(rdir,"DMPS_stress_adjustedmod.Rdata"))

##getting beta values##################################
bval <- getBeta(eset.nosex)
bval[1:5, 1:5]

##DNA Methylated Regions for stress ##Mvalues################################################################
library(DMRcate)

myAnnotation <- cpg.annotate(object = getM(eset.nosex)[,!is.na(samps$Cortisol)], datatype = "array", 
                             what = "M", 
                             analysis.type = "differential", 
                             design = design, 
                             contrasts = FALSE, 
                             coef = 2, 
                             arraytype = "EPIC")


#Your contrast returned 4 individually significant probes; a small but real effect. 
#Consider manually setting the value of pcutoff to return more DMRs, but be warned that doing this increases 
#the risk of Type I errors.


str(myAnnotation)
#List of 7
#$ ID    : chr [1:808554] "cg14817997" "cg26928153" "cg16269199" "cg13869341" ...
#$ stat  : num [1:808554] 0.456 -0.75 0.437 -0.461 0.651 ...
#$ CHR   : chr [1:808554] "chr1" "chr1" "chr1" "chr1" ...
#$ pos   : int [1:808554] 10525 10848 10850 15865 18827 29407 29425 68849 68889 69591 ...
#$ betafc: num [1:808554] 0.00021 -0.000104 0.000222 -0.000228 0.000316 ...
#$ indfdr: num [1:808554] 0.944 0.898 0.947 0.944 0.914 ...
#$ is.sig: logi [1:808554] FALSE FALSE FALSE FALSE FALSE FALSE ...
#- attr(*, "row.names")= int [1:808554] 1 2 3 4 5 6 7 8 9 10 ...
#- attr(*, "class")= chr "annot"

# DMR analysis
DMRs <- dmrcate(myAnnotation, lambda=1000, C=2)
results.ranges <- extractRanges(DMRs)

results.ranges
#GRanges object with 1 range and 6 metadata columns:
#  seqnames            ranges strand |   no.cpgs              minfdr
#<Rle>         <IRanges>  <Rle> | <integer>           <numeric>
#  [1]     chr7 41476044-41476263      * |         6 8.23896505824553e-16
#Stouffer            maxbetafc           meanbetafc
#<numeric>            <numeric>            <numeric>
#  [1] 0.790128683778966 0.000184393685962849 0.000142019580021224
#overlapping.promoters
#<character>
#  [1]             AR4D-001
#-------
#  seqinfo: 1 sequence from an unspecified genome; no seqlengths


# visualization
dmr.table <- data.frame(results.ranges)

#seqnames    start      end width strand no.cpgs       minfdr  Stouffer
#1     chr7 41476044 41476263   220      *       6 8.238965e-16 0.7901287
#maxbetafc   meanbetafc overlapping.promoters
#1 0.0001843937 0.0001420196             AR4D-001



write.table(dmr.table, "/mnt/nas/global/ame/users/Ritika/Projects/DMR1_.csv", 
            sep=",", row.names=FALSE)

# draw the plot for the top DMR
basedir <- #specify base dir to save images
setwd(basedir)

phen.col <- c(rep("orange", 38), rep("blue", 38))

pdf(paste0(basedir,"dmrc.pdf"), onefile=T)   
par(mfrow=c(1,1))
plot <- DMR.plot(ranges = results.ranges, dmr = 1, CpGs = getBeta(eset.nosex), phen.col=phen.col, 
                 what = "Beta", arraytype = "EPIC", genome = "hg19")
print(plot)
dev.off()



#####extracting the cpgs present in the DMR#####################

RR <- as.data.frame(results.ranges)
RR$DMRID <-rownames(RR)
row.names(RR) = NULL
RR$DMRNO <- rownames (RR)
row.names(RR) = RR$DMRID
RR$DMRID = NULL
RR <-RR[order(RR$minfdr), , drop = FALSE]
RR

#Now you need to pull the CpG info from the dmr output that was used to make the results ranges
cgID <- as.data.frame(DMRs$input)

#Look at your RR file, choose a DMR that looks good, and take note of the DMRNO; Run
DMRNUM <- readline(prompt = "What is your DMR Number:") ##fill 1

#Enter the number into the console and hit enter, then run these lines 
#and it should spit out a table listing the probes as well as other useful info

assign(paste0("DMR_",DMRNUM), subset(subset(RR,DMRNO==DMRNUM)))

assign(paste0("DMR_",DMRNUM,"_probelist"), 
       subset(cgID, cgID$CHR==assign(paste0("DMR_",DMRNUM), 
                                     subset(subset(RR,DMRNO==DMRNUM)))$seqnames & cgID$pos>assign(paste0("DMR_",DMRNUM), 
                                                                                                  subset(subset(RR,DMRNO==DMRNUM)))$start-1 & cgID$pos<assign(paste0("DMR_",DMRNUM),
                                                                                                                                                              subset(subset(RR,DMRNO==DMRNUM)))$end+1))

DMR_1_probelist
ID  weights  CHR      pos        betafc    indfdr is.sig
#341630 cg01065373 3.584847 chr7 41476044 1.843937e-04 0.7150366  FALSE

raw          fdr  sig
#341630 2.021358e-21 8.238965e-16 TRUE


DMR_cpgs <- merge(DMR_1_probelist, DMPs, by.x = "ID", by.y = "Name") 

write.table(DMR_cpgs, "/Codes/codes_2020_Jan/Model_codes/DMRcpgs.csv",  sep=",", row.names=FALSE)

###############################################################


