
#######################################################Load Packages########################################################

library(rrBLUP)
library(data.table)
library(tidyverse)

#########################################################Input data#########################################################

args = commandArgs(trailingOnly = TRUE)
MARKER <- args[1] #genotype data set "/data/users/mjahani/JOON_PAV/pav_gwas/SAM.5maf.noheader.sed12.rrblup.csv"
PHENOTYPE <- args[2] #phenotypes"/data/users/mjahani/JOON_PAV/pav_gwas/phenotype.csv"
SAMPLES <-  args[3] #list of samples same order as rrblup input
SNPs <-  args[4] #SNP IDs same order as rrblup input
SAVE_DIR <- args[5] #directory to save result "/data/users/mjahani/JOON_PAV/pav_gwas/result"

##########################################################Read data#########################################################

Markers_impute <- fread(MARKER,
                        header = F)#read genotype file
phenotype <- read.csv(PHENOTYPE,
                      header = T,
                      row.names = 1) 

samples_list <- fread(SAMPLES,
          header = F) %>%
  rename(SAMPLE = V1)

SNP_list <- fread(SNPs,
                      header = F) %>%
  rename(SNP_ID = V1)
##################################################Fit the prediction Model##################################################
phenotype %>% 
  colnames() -> TRAIT_ID

for (trait in 1:4) {
  
GP_MODEL <- mixed.solve(phenotype[,trait], #phenotyoic data
                        Z=Markers_impute, #marker data
                        K=NULL, 
                        method="REML",
                        SE=FALSE, 
                        return.Hinv=FALSE)#calculation of the model
BLUPS <- as.matrix(GP_MODEL$u) #extract matix of BLUPS
BLUE <- as.vector(GP_MODEL$beta) #extract the BLUE value
PREDICTED_PHENOTYPE <- (as.matrix(Markers_impute) %*% (BLUPS)) + BLUE #calculate breeding values

BLUPS %>% 
  as.data.frame() %>% 
  rename(!!TRAIT_ID[trait] := V1) %>%
  cbind(SNP_list,.) -> SNP_list

PREDICTED_PHENOTYPE %>%
  as.data.frame() %>%
  rename(!!TRAIT_ID[trait] := V1) %>%
  cbind(samples_list,.) -> samples_list
  
rm(GP_MODEL,
   BLUPS,
   BLUE,
   PREDICTED_PHENOTYPE)

}

fwrite(samples_list,
       paste0(SAVE_DIR,"/Predicted_phenotype"),
       sep = "\t",
       col.names = T)

fwrite(SNP_list,
       paste0(SAVE_DIR,"/BLUPS"),
       sep = "\t",
       col.names = T)


