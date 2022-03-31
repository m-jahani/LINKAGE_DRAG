library(data.table)
library(dplyr)
library(tidyr)
library(doParallel)
library(foreach)

CPU=6
#DATA <- "/DATA/home/mjahani/LINKADE_DRAG/introg_effect/BLUPS_intg_frq"
#SAVE_DIR <- "/DATA/home/mjahani/LINKADE_DRAG/introg_effect"

args = commandArgs(trailingOnly = TRUE)

DATA <- args[1] #"/DATA/home/mjahani/LINKADE_DRAG/introg_effect/BLUPS_intg_frq"
SAVE_DIR <- args[2] #"/DATA/home/mjahani/LINKADE_DRAG/introg_effect"

#data fro effect of each SNP inGWAS or GP
fread(DATA,
      header = T) -> intg_effects
#Extract SNP list
intg_effects %>% 
  distinct(SNP_ID) %>%
  separate(SNP_ID, 
           into = c("CHR","POS"),
           sep = ":",
           remove = T)-> SNPS

int_frq <- c("annuus.frq0.01","annuus.frq0.02","annuus.frq0.03","annuus.frq0.04","annuus.frq0.05","other.frq0.01","other.frq0.02","other.frq0.03","other.frq0.04","other.frq0.05")

#average effect for target regions
average_traits_frq <-  NULL
for (i in 1:length(int_frq)) {
intg_effects %>%
  filter(!is.na(.data[[int_frq[i]]])) %>%
  summarise_if(is.numeric, mean, na.rm = TRUE) %>%
  mutate(INTREGRESSION = as.character(int_frq[i]),
         DATA = "REAL",
         TYPE = "MEAN") %>% 
  rbind(.,average_traits_frq) -> average_traits_frq
}

count_traits_frq <-  NULL
for (i in 1:length(int_frq)) {
  intg_effects %>%
    filter(!is.na(.data[[int_frq[i]]])) %>%
    select(-contains(".frq")) %>%
    gather(trait,value,2:47)  %>%
    group_by(trait) %>%
    summarise(count=n(), summ=sum(value>0)) %>%
    ungroup() %>%
    mutate(positive_pro = summ/count) %>%
    select(-summ,
           -count) %>% 
    spread(trait,positive_pro) %>%
    mutate(INTREGRESSION = as.character(int_frq[i]),
           DATA = "REAL",
           TYPE = "COUNT") %>% 
    rbind(.,count_traits_frq) -> count_traits_frq
}


average_traits_frq %>%
  fwrite(paste0(SAVE_DIR,
                "/",
                "average_permutation_and_realdata_result"),
         col.names = T,
         append = T,
         sep = "\t",
         quote = F,
         na = "NA")

count_traits_frq %>%
  fwrite(paste0(SAVE_DIR,
                "/",
                "count_permutation_and_realdata_result"),
         col.names = T,
         append = T,
         sep = "\t",
         quote = F,
         na = "NA")

rm(average_traits_frq,
   count_traits_frq)

frq_pair <- list(c("annuus.frq0.03","other.frq0.03"),
               c("annuus.frq0.04","other.frq0.04"),
               c("annuus.frq0.05","other.frq0.05"),
               c("other.frq0.03","annuus.frq0.03"),
               c("other.frq0.04","annuus.frq0.04"),
               c("other.frq0.05","annuus.frq0.05")) 

#Average effect for backgroung region n=10000
for (i in 1:length(int_frq)) {
  #bed file for targeted regions
  intg_effects %>%
    filter(!is.na(.data[[frq_pair[[i]][[1]]]])) %>%
    distinct(.data[[frq_pair[[i]][[1]]]]) %>%
    separate(.data[[frq_pair[[i]][[1]]]], into= c("chr","start_end"), sep=":", remove = F) %>% 
    separate(start_end, into= c("start","end"), sep="-", remove = T) %>% 
    select(chr,
           start,
           end,
           frq_pair[[i]][[1]]) %>%
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
    separate(.data[[frq_pair[[i]][[2]]]], into= c("chr","start_end"), sep=":", remove = F) %>% 
    separate(start_end, into= c("start","end"), sep="-", remove = T) %>% 
    select(chr,
           start,
           end,
           frq_pair[[i]][[2]]) %>%
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

   #read bedtools result for random regions
   assign(paste0("bed",PER),
          fread(paste0(SAVE_DIR,
               "/",
               frq_pair[[i]][[1]],
               "_",
               frq_pair[[i]][[2]],
               "_",
               PER,
               "perm.bed"),
        header = F)) 
  
   system(paste0("rm ",
                 SAVE_DIR,
                 "/",
                 frq_pair[[i]][[1]],
                 "_",
                 frq_pair[[i]][[2]],
                 "_",
                 PER,
                 "perm.bed"))
   
  #find overlapping SNPs with the random regions
   assign(
     paste0("SNP_GENE",PER),
     foreach(i=1:nrow(get(paste0("bed",PER))), .combine='rbind', .errorhandling='stop') %dopar% { 
      mutate(select(filter(SNPS,
                           CHR == as.character(get(paste0("bed",PER))[i,1]),
                           (POS >= as.numeric(get(paste0("bed",PER))[i,2]) &
                              POS <= as.numeric(get(paste0("bed",PER))[i,3]))
      ),
      CHR,POS),
      Gene_ID=as.character(get(paste0("bed",PER))[i,4]))})
   
   rm(list = c(paste0("bed",PER)))
   
  #calculate average of each trait effect for random regions
 intg_effects %>%
   filter(SNP_ID %in% pull(select(mutate(get(paste0("SNP_GENE",PER)),SNP_ID=paste0(CHR,":",POS)),SNP_ID),SNP_ID)) %>%
   summarise_if(is.numeric, mean, na.rm = TRUE) %>%
   mutate(INTREGRESSION = as.character(frq_pair[[i]][[1]]),
          DATA = paste("PERM",PER)) %>%
   fwrite(
     paste0(SAVE_DIR,
            "/",
            "average_permutation_and_realdata_result"),
          col.names = F,
          append = T,
          sep = "\t",
          quote = F,
          na = "NA")
 
 intg_effects %>%
   filter(SNP_ID %in% pull(select(mutate(get(paste0("SNP_GENE",PER)),SNP_ID=paste0(CHR,":",POS)),SNP_ID),SNP_ID)) %>%
     select(-contains(".frq")) %>%
     gather(trait,value,2:47)  %>%
     group_by(trait) %>%
     summarise(count=n(), summ=sum(value>0)) %>%
     ungroup() %>%
     mutate(positive_pro = summ/count) %>%
     select(-summ,
            -count) %>% 
     spread(trait,positive_pro) %>%
     mutate(INTREGRESSION = as.character(frq_pair[[i]][[1]]),
            DATA = paste("PERM",PER),
            TYPE = "COUNT") %>%
   fwrite(
     paste0(SAVE_DIR,
            "/",
            "count_permutation_and_realdata_result"),
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
}

rm(intg_effects,
   SNPS,
   frq_pair)
