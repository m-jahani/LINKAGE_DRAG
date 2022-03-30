library(data.table)
library(tidyverse)

PATH="/DATA/home/mjahani/LINKADE_DRAG/GWAS/result/"

list.files(path = PATH) %>%
  as.data.frame() %>% 
  filter(grepl("*.ps$",.)) -> file_list

fread(paste0(PATH, as.character(a[i,1]))) %>%
  select(SNP_ID = V1) -> RESULT_Beta
for (i in 1:nrow(file_list)) {
  fread(paste0(PATH, as.character(a[i,1]))) %>%
    select(SNP_ID = V1,
           !!gsub(".ps","",gsub("Annuus.ann_fullsam.tranche90_snps_bi_AN50_beagle_AF99_","",as.character(a[i,1]))) := V2) %>%
    full_join(.,RESULT_Beta) -> RESULT_Beta
}

fwrite(RESULT_Beta,paste0(PATH,"RESULT_Beta"),
       col.names = T,
       sep = "\t",
       quote = F)

rm(RESULT_Beta)

fread(paste0(PATH, as.character(a[i,1]))) %>%
  select(SNP_ID = V1) -> RESULT_P
for (i in 1:nrow(file_list)) {
  fread(paste0(PATH, as.character(a[i,1]))) %>%
    select(SNP_ID = V1,
           !!gsub(".ps","",gsub("Annuus.ann_fullsam.tranche90_snps_bi_AN50_beagle_AF99_","",as.character(a[i,1]))) := V4) %>%
    full_join(.,RESULT_P) -> RESULT_P
}

fwrite(RESULT_P,paste0(PATH,"RESULT_P"),
       col.names = T,
       sep = "\t",
       quote = F)

rm(RESULT_P)