library(data.table)
library(tidyverse)
library(foreach)
library(doParallel)


args = commandArgs(trailingOnly = TRUE)
TPED <- args[1]
TFAM <- args[2]
SAVE_DIR <- args[3]


fread(TPED) -> data
fread(TFAM) %>% 
  select(V1) -> IDs

registerDoParallel(cores=2)
seq(5,ncol(data),2) -> SAMPLES
SAMPLES ->  SAMPLES1

for (i in SAMPLES) {
  SAMPLES1[SAMPLES1 != i] -> SAMPLES1
  for (j in SAMPLES1) {
    result <- select(data,1:4)
    for (r in 1:100) {
      data %>%
        select(1:4,
               hap1 =i,
               hap2 = i+1) %>%
        mutate(randomizer = sample(2, size = nrow(data), replace = TRUE)) %>%
        mutate(gamet1 = ifelse(randomizer == 1,hap1,hap2)) %>%
        select(1:4,
               gamet1) %>%
        full_join(.,
                  select(mutate(mutate(select(data,
                                              1:4,
                                              hap1 = j,
                                              hap2 = j+1),
                                       randomizer = sample(2, size = nrow(data), replace = TRUE)),
                                gamet2 = ifelse(randomizer == 1,hap1,hap2)),
                         1:4,
                         gamet2)) %>%
        select(1:4,
               !!paste0(IDs[which(SAMPLES==i),1],
                        "_",
                        IDs[which(SAMPLES==j),1],
                        "_hap1_r",r)  := gamet1,
               !!paste0(IDs[which(SAMPLES==i),1],
                        "_",
                        IDs[which(SAMPLES==j),1],
                        "_hap2_r",r) := gamet2) %>%
        full_join(.,result) -> result
    }
    
    fwrite(result,
           paste0(SAVE_DIR,
                  "/",
                  IDs[which(SAMPLES==i),1],
                  "_",
                  IDs[which(SAMPLES==j),1],
                  ".tped"),
           sep = "\t",
           col.names = F)
    
    result %>%
      select(1:4) %>%
      fwrite(paste0(SAVE_DIR,
                    "/",
                    IDs[which(SAMPLES==i),1],
                    "_",
                    IDs[which(SAMPLES==j),1],
                    ".map"),
             sep = "\t",
             col.names = F)
    
    data.frame(matrix(ncol = 2,nrow = 100)) %>% 
      mutate(X1=paste0(IDs[which(SAMPLES==i),1],
                      "_",
                      IDs[which(SAMPLES==j),1],
                      "_R_",
                      row_number())) %>%
      mutate(X2 = X1,
             X3 = 0,
             X4 = 0,
             X5 = 0,
             X6 = -9) %>%
      fwrite(paste0(SAVE_DIR,
                    "/",
                    IDs[which(SAMPLES==i),1],
                    "_",
                    IDs[which(SAMPLES==j),1],
                    ".tfam"),
             sep = "\t",
             col.names = F)
  }
  
}


