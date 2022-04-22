library(dplyr)
library(tidyr)
library(data.table)
library(foreach)
library(doParallel)

GWAS_DIR="/DATA/home/mjahani/LINKADE_DRAG/new_method/GWAS/result"
SAVE_DIR="/DATA/home/mjahani/LINKADE_DRAG/new_method/GWAS/result"

list.files(path=GWAS_DIR) %>% 
  as.data.frame() %>% filter(grepl(".ps",.)) %>% 
  mutate(list=gsub("SAM_introgression_donor_","",.))  %>% 
  separate(list, into=c("DONOR_TRAIT","LOCATION"),sep="__",remove=T) %>% 
  mutate(LOCATION=gsub(".ps","",LOCATION)) %>% 
  mutate(DONOR_TRAIT=gsub("2nd_GERMPLASM_maf0.03_","Secondary_Germplasm__",DONOR_TRAIT)) %>% 
  mutate(DONOR_TRAIT=gsub("ANNUUS_maf0.03_","Wild_Annuus__",DONOR_TRAIT)) %>% 
  separate(DONOR_TRAIT, into=c("DONOR","TRAIT"),sep="__",remove=T) %>% 
  mutate(TRAIT = gsub("plant_biomass","biomass",TRAIT)) %>%
  mutate(TRAIT = gsub("plant_height","height",TRAIT)) %>%
  mutate(TRAIT = gsub("seed_lxw","seed_size",TRAIT)) %>%
  mutate(LOCATION = gsub("iowa","IA",LOCATION)) %>%
  mutate(LOCATION = gsub("georgia","GA",LOCATION)) %>%
  mutate(LOCATION = gsub("UBC","BC",LOCATION)) -> DATA_LIST

  registerDoParallel(cores=60)
  result <- foreach(i = 1:nrow(DATA_LIST), .combine='rbind') %dopar% {
  fread(DATA_LIST[i,1],header=F) %>% 
    select(ID = V1,
           Beta = V2,
           P = V4) %>%
    mutate(SIGN = ifelse(Beta>=0,"POSITIVE","NEGATIVE")) %>%
    mutate(threshold = ifelse(DATA_LIST[i,2] == "Secondary_Germplasm", 2.880814,4.792041)) %>%
    separate(ID, into = c("chr","START_END"), sep = ":", remove = T) %>%
    separate(START_END, into = c("start","end"), sep = "-", remove = T) %>%
    mutate(chr = gsub("Ha412HOChr","Chr",chr)) %>%
    filter((log10(P)*-1) >= threshold) %>%
    mutate(TRAIT = DATA_LIST[i,3]) %>%
    mutate(DONOR = DATA_LIST[i,2]) %>%
    select(chr,
           start,
           end,
           SIGN,
           TRAIT,
           DONOR) 
}

  


  result %>%
    filter(DONOR=="Wild_Annuus") %>%
    filter(SIGN=="POSITIVE") %>%
    select(-DONOR,
           -SIGN) %>%
    distinct(chr,
             start,
             end,
             TRAIT) %>%
    group_by(chr,
             start,
             end) %>%
    mutate(id = row_number()) %>%
    ungroup() %>%
    mutate(id=gsub(1,"one",id)) %>%
    mutate(id=gsub(2,"two",id)) %>%
    spread(id,
           TRAIT) %>%
    mutate(group=paste0(one," & ",two)) %>%
    select(-one,
           -two) %>%
    mutate(group=gsub(" & NA","",group)) %>%
fwrite(paste0(SAVE_DIR,"/Wild_Annuus_POSITIVE_SHINY.csv"),
           sep = ",",
           col.names = T,
           quote = F)

  result %>%
    filter(DONOR=="Wild_Annuus") %>%
    filter(SIGN=="NEGATIVE") %>%
    select(-DONOR,
           -SIGN) %>%
    distinct(chr,
             start,
             end,
             TRAIT) %>%
    group_by(chr,
             start,
             end) %>%
    mutate(id = row_number()) %>%
    ungroup() %>%
    mutate(id=gsub(1,"one",id)) %>%
    mutate(id=gsub(2,"two",id)) %>%
    spread(id,
           TRAIT) %>%
    mutate(group=paste0(one," & ",two)) %>%
    select(-one,
           -two) %>%
    mutate(group=gsub(" & NA","",group)) %>%
    fwrite(paste0(SAVE_DIR,"/Wild_Annuus_NEGATIVE_SHINY.csv"),
           sep = ",",
           col.names = T,
           quote = F)

  result %>%
    filter(DONOR=="Secondary_Germplasm") %>%
    filter(SIGN=="POSITIVE") %>%
    select(-DONOR,
           -SIGN) %>%
    distinct(chr,
             start,
             end,
             TRAIT) %>%
    rename(group = TRAIT) %>%
    fwrite(paste0(SAVE_DIR,"/Secondary_Germplasm_POSITIVE_SHINY.csv"),
           sep = ",",
           col.names = T,
           quote = F)


  
  
  