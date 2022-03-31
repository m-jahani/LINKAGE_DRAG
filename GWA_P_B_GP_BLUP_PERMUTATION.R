library(data.table)
library(dplyr)
library(tidyr)
library(doParallel)
library(foreach)

CPU=50
args = commandArgs(trailingOnly = TRUE)

DATA <- args[1] #"/DATA/home/mjahani/LINKADE_DRAG/introg_effect/BLUPS_intg_frq"
SAVE_DIR <- args[2] #"/DATA/home/mjahani/LINKADE_DRAG/introg_effect"

#data fro effect of each SNP inGWAS or GP
fread(DATA,
      header = T) -> intg_effects


int_frq <- c("annuus.frq0.03","annuus.frq0.04","annuus.frq0.05","other.frq0.03","other.frq0.04","other.frq0.05")

#average effect for target regions
BLUP_average_traits_frq <-  NULL
for (i in 1:length(int_frq)) {
  intg_effects %>%
    filter(!is.na(.data[[int_frq[i]]])) %>%
    select(starts_with("BLUP_")) %>%
    rename_with(~ tolower(gsub("BLUP_", "",.x, fixed = TRUE))) %>% 
    summarise_if(is.numeric, mean, na.rm = TRUE) %>%
    mutate(INTREGRESSION = as.character(int_frq[i]),
           DATA = "REAL",
           TYPE = "MEAN",
           SOURCE = "GP_BLUP") %>% 
    rbind(.,BLUP_average_traits_frq) -> BLUP_average_traits_frq
}

BLUP_count_traits_frq <-  NULL
for (i in 1:length(int_frq)) {
  intg_effects %>%
    filter(!is.na(.data[[int_frq[i]]])) %>%
    select(starts_with("BLUP_")) %>%
    rename_with(~ tolower(gsub("BLUP_", "",.x, fixed = TRUE))) %>% 
    gather(trait,value)  %>%
    group_by(trait) %>%
    summarise(count=n(), summ=sum(value>0)) %>%
    ungroup() %>%
    mutate(positive_pro = summ/count) %>%
    select(-summ,
           -count) %>% 
    spread(trait,positive_pro) %>%
    mutate(INTREGRESSION = as.character(int_frq[i]),
           DATA = "REAL",
           TYPE = "POS_RATIO",
           SOURCE = "GP_BLUP") %>% 
    rbind(.,BLUP_count_traits_frq) -> BLUP_count_traits_frq
}


BETA_count_traits_frq <-  NULL
for (i in 1:length(int_frq)) {
  intg_effects %>%
    filter(!is.na(.data[[int_frq[i]]])) %>%
    select(starts_with("Beta_")) %>%
    rename_with(~ tolower(gsub("Beta_", "",.x, fixed = TRUE))) %>% 
    gather(trait,value)  %>%
    group_by(trait) %>%
    summarise(count=n(), summ=sum(value>0)) %>%
    ungroup() %>%
    mutate(positive_pro = summ/count) %>%
    select(-summ,
           -count) %>% 
    spread(trait,positive_pro) %>%
    mutate(INTREGRESSION = as.character(int_frq[i]),
           DATA = "REAL",
           TYPE = "POS_RATIO",
           SOURCE = "GWA_Beta") %>% 
    rbind(.,BETA_count_traits_frq) -> BETA_count_traits_frq
}


P_count_traits_frq <-  NULL
for (i in 1:length(int_frq)) {
  intg_effects %>%
    filter(!is.na(.data[[int_frq[i]]])) %>%
    select(starts_with("P_")) %>%
    rename_with(~ tolower(gsub("P_", "",.x, fixed = TRUE))) %>% 
    gather(trait,value)  %>%
    group_by(trait) %>%
    summarise(count=n(), summ=sum((log10(value)*-1)>4)) %>%
    ungroup() %>%
    mutate(positive_pro = summ/count) %>%
    select(-summ,
           -count) %>% 
    spread(trait,positive_pro) %>%
    mutate(INTREGRESSION = as.character(int_frq[i]),
           DATA = "REAL",
           TYPE = "OUTLIER_P",
           SOURCE = "GWA_P") %>% 
    rbind(.,P_count_traits_frq) -> P_count_traits_frq
}


BLUP_average_traits_frq %>%
  fwrite(paste0(SAVE_DIR,
                "/",
                "BLUP_average_permutation_and_realdata_result"),
         col.names = T,
         append = T,
         sep = "\t",
         quote = F,
         na = "NA")

BLUP_count_traits_frq %>%
  fwrite(paste0(SAVE_DIR,
                "/",
                "BLUP_positive_ratio_permutation_and_realdata_result"),
         col.names = T,
         append = T,
         sep = "\t",
         quote = F,
         na = "NA")

BETA_count_traits_frq %>%
  fwrite(paste0(SAVE_DIR,
                "/",
                "BETA_positive_ratio_permutation_and_realdata_result"),
         col.names = T,
         append = T,
         sep = "\t",
         quote = F,
         na = "NA")

P_count_traits_frq %>%
  fwrite(paste0(SAVE_DIR,
                "/",
                "P_oulier_permutation_and_realdata_result"),
         col.names = T,
         append = T,
         sep = "\t",
         quote = F,
         na = "NA")

rm(BLUP_average_traits_frq,
   BLUP_count_traits_frq,
   BETA_count_traits_frq,
   P_count_traits_frq,
   int_frq)

#PERMUTATION##############################################################
frq_pair <- list(c("annuus.frq0.03","other.frq0.03"),
                 c("annuus.frq0.04","other.frq0.04"),
                 c("annuus.frq0.05","other.frq0.05"),
                 c("other.frq0.03","annuus.frq0.03"),
                 c("other.frq0.04","annuus.frq0.04"),
                 c("other.frq0.05","annuus.frq0.05")) 

for (i in 1:length(frq_pair)) { #loop over frequencies
#Average effect for backgroung region n=10000
  #bed file for targeted regions
  intg_effects %>%
    filter(!is.na(.data[[frq_pair[[i]][[1]]]])) %>%
    distinct(.data[[frq_pair[[i]][[1]]]]) %>%
    separate(.data[[frq_pair[[i]][[1]]]], into= c("chr","start_end"), sep=":", remove = T) %>% 
    separate(start_end, into= c("start","end"), sep="-", remove = T) %>% 
    select(chr,
           start,
           end) %>%
    fwrite(paste0(SAVE_DIR,
                  "/",
                  frq_pair[[i]][[1]],
                  ".bed"),
           sep = "\t",
           col.names = F,
           quote = F)
  #bed file for the other gene pool to exclude
  intg_effects %>%
    filter(!is.na(.data[[frq_pair[[i]][[2]]]])) %>%
    distinct(.data[[frq_pair[[i]][[2]]]]) %>%
    separate(.data[[frq_pair[[i]][[2]]]], into= c("chr","start_end"), sep=":", remove = T) %>% 
    separate(start_end, into= c("start","end"), sep="-", remove = T) %>% 
    select(chr,
           start,
           end) %>%
    fwrite(paste0(SAVE_DIR,
                  "/",
                  frq_pair[[i]][[2]],
                  ".bed"),
           sep = "\t",
           col.names = F,
           quote = F)
  #permutation for picking background regions and calculate average effect
  registerDoParallel(cores=CPU)
  #for (PER in 1:10000) { 
  foreach(PER=1:10000, .combine='rbind', .errorhandling='stop') %dopar% {#permutation start
    #using bedtools to pick random regions (same size as target regions) and exclude second genepool regions
    system(paste0("bedtools shuffle -i ",
                  SAVE_DIR,
                  "/",
                  frq_pair[[i]][[1]],
                  ".bed -g ",
                  SAVE_DIR,
                  "/",
                  "HA412v2_chromosome.bed -excl ",
                  SAVE_DIR,
                  "/",
                  frq_pair[[i]][[2]],
                  ".bed -noOverlapping > ",
                  SAVE_DIR,
                  "/",
                  frq_pair[[i]][[1]],
                  "_",
                  frq_pair[[i]][[2]],
                  "_",
                  PER,
                  "perm.bed")) 
    
    system(paste0("bedtools intersect -wa -wb -a ",
                  SAVE_DIR,
                  "/",
                  frq_pair[[i]][[1]],
                  "_",
                  frq_pair[[i]][[2]],
                  "_",
                  PER,
                  "perm.bed -b ",
                  SAVE_DIR,
                  "/ALL_SNPS.bed | ", 
                  paste("awk '{print $4","\t","$5","\t","$1",":","$2","-","$3}' > ",sep = '"'),
                  SAVE_DIR,
                  "/",
                  frq_pair[[i]][[1]],
                  "_",
                  frq_pair[[i]][[2]],
                  "_",
                  PER,
                  "perm_SNP_OVERLAP.bed"))
             
    
    #read bedtools result for random regions
    assign(paste0("SNP_GENE",PER),
           fread(paste0(SAVE_DIR,
                        "/",
                        frq_pair[[i]][[1]],
                        "_",
                        frq_pair[[i]][[2]],
                        "_",
                        PER,
                        "perm_SNP_OVERLAP.bed"),
                 header = F) %>%
             rename(CHR = V1,
                    POS = V2,
                    Gene_ID = V3))
    
    system(paste0("rm ",
                  SAVE_DIR,
                  "/",
                  frq_pair[[i]][[1]],
                  "_",
                  frq_pair[[i]][[2]],
                  "_",
                  PER,
                  "perm.bed ",
                  SAVE_DIR,
                  "/",
                  frq_pair[[i]][[1]],
                  "_",
                  frq_pair[[i]][[2]],
                  "_",
                  PER,
                  "perm_SNP_OVERLAP.bed"))

    #calculate average of each trait effect for random regions
    intg_effects %>%
      filter(SNP_ID %in% pull(select(mutate(get(paste0("SNP_GENE",PER)),SNP_ID=paste0(CHR,":",POS)),SNP_ID),SNP_ID)) %>%
      select(starts_with("BLUP_")) %>%
      rename_with(~ tolower(gsub("BLUP_", "",.x, fixed = TRUE))) %>% 
      summarise_if(is.numeric, mean, na.rm = TRUE) %>%
      mutate(INTREGRESSION = as.character(frq_pair[[i]][[1]]),
             DATA = paste0("PERM",PER),
             TYPE = "MEAN",
             SOURCE = "GP_BLUP") %>%
      fwrite(
        paste0(SAVE_DIR,
               "/",
               "BLUP_average_permutation_and_realdata_result"),
        col.names = F,
        append = T,
        sep = "\t",
        quote = F,
        na = "NA")
    
    intg_effects %>%
      filter(SNP_ID %in% pull(select(mutate(get(paste0("SNP_GENE",PER)),SNP_ID=paste0(CHR,":",POS)),SNP_ID),SNP_ID)) %>%
      select(starts_with("BLUP_")) %>%
      rename_with(~ tolower(gsub("BLUP_", "",.x, fixed = TRUE))) %>% 
      gather(trait,value)  %>%
      group_by(trait) %>%
      summarise(count=n(), summ=sum(value>0)) %>%
      ungroup() %>%
      mutate(positive_pro = summ/count) %>%
      select(-summ,
             -count) %>% 
      spread(trait,positive_pro) %>%
      mutate(INTREGRESSION = as.character(frq_pair[[i]][[1]]),
             DATA =  paste0("PERM",PER),
             TYPE = "POS_RATIO",
             SOURCE = "GP_BLUP") %>% 
      fwrite(
        paste0(SAVE_DIR,
               "/",
               "BLUP_positive_ratio_permutation_and_realdata_result"),
        col.names = F,
        append = T,
        sep = "\t",
        quote = F,
        na = "NA")

    intg_effects %>%
      filter(SNP_ID %in% pull(select(mutate(get(paste0("SNP_GENE",PER)),SNP_ID=paste0(CHR,":",POS)),SNP_ID),SNP_ID)) %>%
      select(starts_with("Beta_")) %>%
      rename_with(~ tolower(gsub("Beta_", "",.x, fixed = TRUE))) %>% 
      gather(trait,value)  %>%
      group_by(trait) %>%
      summarise(count=n(), summ=sum(value>0)) %>%
      ungroup() %>%
      mutate(positive_pro = summ/count) %>%
      select(-summ,
             -count) %>% 
      spread(trait,positive_pro) %>%
      mutate(INTREGRESSION = as.character(frq_pair[[i]][[1]]),
             DATA = paste0("PERM",PER),
             TYPE = "POS_RATIO",
             SOURCE = "GWA_Beta") %>% 
      fwrite(
        paste0(SAVE_DIR,
               "/",
               "BETA_positive_ratio_permutation_and_realdata_result"),
        col.names = F,
        append = T,
        sep = "\t",
        quote = F,
        na = "NA")
  

    intg_effects %>%
      filter(SNP_ID %in% pull(select(mutate(get(paste0("SNP_GENE",PER)),SNP_ID=paste0(CHR,":",POS)),SNP_ID),SNP_ID)) %>%
      select(starts_with("P_")) %>%
      rename_with(~ tolower(gsub("P_", "",.x, fixed = TRUE))) %>% 
      gather(trait,value)  %>%
      group_by(trait) %>%
      summarise(count=n(), summ=sum((log10(value)*-1)>4)) %>%
      ungroup() %>%
      mutate(positive_pro = summ/count) %>%
      select(-summ,
             -count) %>% 
      spread(trait,positive_pro) %>%
      mutate(INTREGRESSION = as.character(frq_pair[[i]][[1]]),
             DATA = paste0("PERM",PER),
             TYPE = "OUTLIER_P",
             SOURCE = "GWA_P") %>% 
      fwrite(
        paste0(SAVE_DIR,
               "/",
               "P_oulier_permutation_and_realdata_result"),
        col.names = F,
        append = T,
        sep = "\t",
        quote = F,
        na = "NA")
    
    rm(list = c(paste0("SNP_GENE",PER)))
  }#permutation end
  
  #delete bed files 
  system(paste0("rm ",
                SAVE_DIR,
                "/",
                frq_pair[[i]][[1]],
                ".bed ",
                SAVE_DIR,
                "/",
                frq_pair[[i]][[2]],
                ".bed "))
}#frequencies loop end
