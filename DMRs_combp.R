path <- "/mnt/nas/global/ame/users/Ritika/Projects/FELICITy_Project/EPIC_DATA_norm/Codes/DMRs_combp/"


#########################################################
#Recreating the results using the R version of comb-p
source(paste0(path,"source_dmrs.R"))
library(parallel)

#preparing the data	
#A dataframe with colname name "chr","start", "end","p" and "probe", 
#indicating chromosome (1,2,3,...,X,Y), chromosome start and end position, P value and probe names

cmbp_df   <- dplyr::rename(dmrff_df,chr=CHR,p=p.value,probe=IlmnID) %>% mutate(start=MAPINFO, end=MAPINFO+1)

#the same grid: 
comb_p_res <- lapply(c(500,5000), function(mg) {
  
  unlink("resu_combp.csv")
  
  is_success <- tryCatch({combp_local(cmbp_df,dist.cutoff=mg,bin.size=310,seed=1e-3,region_plot=F,mht_plot=F,nCores=25,verbose=TRUE)},error=function(e) return(c(NA,NA,NA)))
  
  if (identical(is_success,c(NA,NA,NA))) return(NULL)
  
  #reading in the results
  combp_res <- tryCatch({read.csv("resu_combp.csv",sep=",",header=T,stringsAsFactors=F)},error=function(e) return(NULL))
  
  if (is.null(combp_res)) return(NULL)
  
  return(combp_res) })

save(comb_p_res,file=combp_outfile)

return(T)

})
##################Paths#########################################################################################
###

path <- "/mnt/nas/global/ame/users/Ritika/Projects/FELICITy_Project/EPIC_DATA_norm/Codes/DMRs_combp/"

setwd("/mnt/nas/global/ame/users/Ritika/Projects/FELICITy_Project/EPIC_DATA_norm/Codes/DMRs_combp/")

#############Variable of interest:- PDQ############################################
#Recreating the results using the R version of comb-p
source(paste0(path,"source_dmrs.R"))
library(parallel)

rdir <-  "/mnt/nas/global/ame/users/Ritika/Projects/FELICITy_Project/EPIC_DATA_norm/EPIC data/"

load(paste0(rdir, "DMPS_pdq.Rdata"))

rownames(DMPs1) <- DMPs1$Name

##uploading the annotation file of EPIC### 
ppath <- "/mnt/nas/global/ame/users/Ritika/Projects/KORA_Age_Project/KORA_AGE_data/Full_data/Data/"
anno <- read.csv(paste0(ppath,"MethylationEPIC_v-1-0_B4.csv"),skip=7,fill=T , stringsAsFactors=F) 

##selecting the variables
anno <- anno[anno$IlmnID %in% rownames(DMPs1), c("IlmnID","CHR","MAPINFO", "UCSC_RefGene_Name")]
dim(anno)
rownames(anno) <- anno$IlmnID

head(anno)

#Matching the columsn and rows
anno <- anno[rownames(DMPs1),]
colnames(anno)

##Merging Annotation file and DMPs 
DMR_dff <- cbind(anno,DMPs1) 

DMR_dff  <- (DMR_dff [,-c(5,17,18)])

###using Comb commands to check for DMRs#############################################

cmbp_df   <- dplyr::rename(DMR_dff, chr=CHR, p=BC_pvals, probe=IlmnID) %>% mutate(start=MAPINFO, end=MAPINFO+1)


c <- combp_local(cmbp_df,dist.cutoff=500,bin.size=310,seed=0.05,region_plot=F,mht_plot=F,nCores=25,verbose=TRUE)

#the same grid: 
comb_p_res <- lapply(c(500,1000), function(mg) {
  
  unlink("resu_combp.csv")
  
  is_success <- tryCatch({combp_local(cmbp_df,dist.cutoff=mg,bin.size=310,seed=1e-3,region_plot=F,mht_plot=F,nCores=25,verbose=TRUE)},error=function(e) return(c(NA,NA,NA)))
  
  if (identical(is_success,c(NA,NA,NA))) return(NULL)
  
  #reading in the results
  combp_res <- tryCatch({read.csv("resu_combp.csv",sep=",",header=T,stringsAsFactors=F)},error=function(e) return(NULL))
  
  if (is.null(combp_res)) return(NULL)
  
  return(combp_res) })

save(comb_p_res,file=combp_outfile) ##make your own combfile