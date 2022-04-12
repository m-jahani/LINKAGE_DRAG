

#######################################################Load Packages########################################################

library(rrBLUP)
library(data.table)
library(dplyr)
library(foreach)
library(doParallel)

#########################################################Input data#########################################################

args = commandArgs(trailingOnly = TRUE)
MARKER <- args[1] #genotype data set "/data/users/mjahani/JOON_PAV/pav_gwas/SAM.5maf.noheader.sed12.rrblup.csv"
PHENOTYPE <- args[2] #phenotypes"/data/users/mjahani/JOON_PAV/pav_gwas/phenotype.csv"
SAVE_DIR <- args[3] #directory to save result "/data/users/mjahani/JOON_PAV/pav_gwas/result"
PERMUTATION <- 100
registerDoParallel(cores=10000)
##########################################################Read data#########################################################

Markers_impute <- fread(MARKER,
                        header = F)#read genotype file
phenotype <- read.csv(PHENOTYPE,
                      header = T,
                      row.names = 1) 

##################################################Fit the prediction Model##################################################

phenotype %>% 
  colnames() -> TRAIT_ID

final <- data.frame(a=seq(1,as.numeric(PERMUTATION),1),b="PERM") %>% mutate(data_type=paste0(b,a)) %>% select(data_type) %>% rbind(.,"observation")
P_VALUES <- data.frame(P=c("LARGER","SMALLER"))
for (trait in 1:length(TRAIT_ID)) {
  trait_result <- data.frame(data_type = "observation")
  #observation
  as.data.frame(
    as.matrix(
      mixed.solve(phenotype[,trait], #phenotyoic data
                  Z=Markers_impute, #marker data
                  K=NULL, 
                  method="REML",
                  SE=FALSE, 
                  return.Hinv=FALSE)$u)) %>% #calculation of BLUP
    summarise(!!TRAIT_ID[trait] := mean(V1)) %>%
    cbind(trait_result,.) -> trait_result
  
 # for (PERM in 1:10000) {

     PERMUT_result <- foreach(PERM= 1:as.numeric(PERMUTATION), .combine='rbind') %dopar% {
  as.data.frame(
    as.matrix(
      mixed.solve(phenotype[,trait], #phenotyoic data
                  Z=Markers_impute[sample(nrow(Markers_impute)),] , #shuufled marker data
                  K=NULL, 
                  method="REML",
                  SE=FALSE, 
                  return.Hinv=FALSE)$u)) %>%
    summarise(!!TRAIT_ID[trait] := mean(V1)) %>%
      mutate(data_type = paste0("PERM",PERM)) %>%
      select(data_type,everything(.)) 
  }
  
     trait_result %>%
       rbind(.,PERMUT_result) %>%
       full_join(.,final) -> final
       

       filter(PERMUT_result,TRAIT_ID[trait] < as.numeric(select(trait_result,TRAIT_ID[trait]))) %>%
       summarise(P = n()/as.numeric(PERMUTATION)) %>% 
         mutate(P_type="LARGER") %>%
         mutate(TRAIT=TRAIT_ID[trait]) %>%
         mutate(Significant=ifelse(P<0.05,"Significant","non_Significant")) %>%
         select(TRAIT,P_type,P,Significant) %>%
         rbind(.,
               filter(PERMUT_result,TRAIT_ID[trait] > as.numeric(select(trait_result,TRAIT_ID[trait]))) %>%
                 summarise(P = n()/as.numeric(PERMUTATION)) %>% 
                 mutate(P_type="SMALLER") %>%
                 mutate(TRAIT=TRAIT_ID[trait]) %>%
                 select(TRAIT,P_type,TRAIT,P) %>%
                 mutate(Significant=ifelse(P<0.05,"Significant","non_Significant")) %>%
                 select(TRAIT,P_type,P,Significant)) %>%
         fwrite(paste0(SAVE_DIR,"/P_VALUES"),
                sep = "\t",
                col.names = T,
                quote = F,
                append = T)
         

     
     rm(trait_result,
        PERMUT_result)
  
  
}

fwrite(final,
       paste0(SAVE_DIR,"/average_BLUP_PERMUTAION_OBSERVATION"),
       sep = "\t",
       col.names = T)


