library(rrBLUP)
library(data.table)
library(dplyr)
library(tidyr)
library(foreach)
library(doParallel)

args = commandArgs(trailingOnly = TRUE)

MARKER <- args[1] #"/DATA/home/mjahani/LINKADE_DRAG/new_method/GP/SAM_introgression_donor_2nd_GERMPLASM_maf0.03_rrBLUP.in"
PHENOTYPE <- args[2] #"/DATA/home/mjahani/LINKADE_DRAG/new_method/GP/phenotype_common_georgia_corrected.csv"
IDs <- args[3]
FRQs <- args[4]
DONOR <- args[5]
SAVE_DIR <- args[6]


# MARKER <- "/DATA/home/mjahani/LINKADE_DRAG/new_method/GP/SAM_introgression_donor_2nd_GERMPLASM_maf0.03_rrBLUP.in"
# PHENOTYPE <- "/DATA/home/mjahani/LINKADE_DRAG/new_method/GP/phenotype_common_georgia_corrected.csv"
# IDs <- "/DATA/home/mjahani/LINKADE_DRAG/new_method/GP/SAM_introgression_donor_2nd_GERMPLASM_maf0.03_SNP_list"
# FRQs <- "/DATA/home/mjahani/LINKADE_DRAG/new_method/GP/FREQ_VS_BLUP/frquency/SAM_introgression_donor_2nd_GERMPLASM_maf0.03_introgression.frq"
# DONOR <- "2nd_GERMPLASM"
# SAVE_DIR <- "/DATA/home/mjahani/LINKADE_DRAG/new_method/GP/FREQ_VS_BLUP"

##########################################################Read data#########################################################

Markers_impute <- fread(MARKER,
                        header = F)#read genotype file

phenotype <- read.csv(PHENOTYPE,
                      header = T,
                      row.names = 1) 

VARIAN_IDS <- fread(IDs,
                              header = F) %>% 
  rename(VARIANT_ID = V1)

VARIANT_FREQUENCY <- read.csv(FRQs,
                            header = F,
                            sep="\t") %>% 
  rename(POSITION = V1,
         FRQ = V2) %>%
  separate(POSITION, into = c("CHR","START_END"), remove = F, sep = ":") %>%
  separate(START_END, into = c("START","END"), remove = F, sep = "-") %>%
  mutate(CHR = as.numeric(gsub("Ha412HOChr","",CHR))) %>%
  mutate(VARIANT_ID = paste0(CHR,":",START)) %>%
  select(VARIANT_ID,
         FRQ)
##################################################Fit the prediction Model##################################################

phenotype %>% 
  colnames() -> TRAIT_ID


registerDoParallel(cores=60)

  BLUPS <- foreach(trait = 1:length(TRAIT_ID), .combine='rbind') %dopar% {

                as.data.frame(
                  as.matrix(
                    mixed.solve(phenotype[,trait], #phenotyoic data
                                Z=Markers_impute, #marker data
                                K=NULL, 
                                method="REML",
                                SE=FALSE, 
                                return.Hinv=FALSE)$u)) %>%
                  mutate(TRAIT = TRAIT_ID[trait]) %>%
                  rename(BLUP = V1) %>% 
                  cbind(.,VARIAN_IDS) %>%
                  full_join(.,VARIANT_FREQUENCY) %>%
                  select(TRAIT,
                         BLUP,
                         FRQ)
}
 
 
  
  BLUPS %>%
  mutate(TRAIT = gsub("plant_biomass","biomass",TRAIT)) %>%
    mutate(TRAIT = gsub("plant_height","height",TRAIT)) %>%
    mutate(TRAIT = gsub("seed_lxw","seed_size",TRAIT)) %>%
    separate(TRAIT, into = c("TRAIT","LOCATION"), sep = "__", remove = T) %>%
    mutate(LOCATION = gsub("iowa","IA",LOCATION)) %>%
    mutate(LOCATION = gsub("georgia","GA",LOCATION)) %>%
    mutate(LOCATION = gsub("UBC","BC",LOCATION)) %>%
    fwrite(paste0(SAVE_DIR,
                  "/",
                  DONOR,
                  "_ALL_INTROGRESSION_EFFECTS_FRQUENCY"),
           col.names = T,
           sep = ",",
           quote = F)
   
  
  RESULT <- matrix(nrow = length(TRAIT_ID),
                   ncol = 4)
  for (trait in 1:length(TRAIT_ID)) {

      FRQU <- pull(filter(BLUPS,TRAIT == TRAIT_ID[trait]),FRQ)
      BLUPs <- pull(filter(BLUPS,TRAIT == TRAIT_ID[trait]),BLUP)
      fit=lm(FRQU ~ BLUPs)
      RESULT[trait,1] <- TRAIT_ID[trait] #trait name
      RESULT[trait,2] <- summary(fit)$r.squared #R-squared
      RESULT[trait,3] <- summary(fit)$coefficients[2,1] #Beta1 (slope)
      RESULT[trait,4] <- summary(fit)$coefficients[2,4] #Pvalue of Beta 1
      rm(FRQU,
         BLUPs)
  }
  
  RESULT %>%
    as.data.frame() %>%
    rename(TRAIT = V1,
           R_SQUARED = V2,
           SLOPE = V3,
           P_SLOPE = V4) %>%
    mutate(TRAIT = gsub("plant_biomass","biomass",TRAIT)) %>%
    mutate(TRAIT = gsub("plant_height","height",TRAIT)) %>%
    mutate(TRAIT = gsub("seed_lxw","seed_size",TRAIT)) %>%
    separate(TRAIT, into = c("TRAIT","LOCATION"), sep = "__", remove = T) %>%
    mutate(LOCATION = gsub("iowa","IA",LOCATION)) %>%
    mutate(LOCATION = gsub("georgia","GA",LOCATION)) %>%
    mutate(LOCATION = gsub("UBC","BC",LOCATION)) %>%
    mutate(SIGNIFICANT = ifelse(as.numeric(P_SLOPE) <= 0.05,"YES","NO")) %>%
    fwrite(paste0(SAVE_DIR,
                  "/",
                  DONOR,
                  "_linear_model_result"),
           col.names = T,
           sep = ",",
           quote = F)
    
    
  