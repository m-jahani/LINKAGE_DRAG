####################################################### Load Packages########################################################

library(rrBLUP)
library(data.table)
library(dplyr)
library(foreach)
library(doParallel)

######################################################### Input data#########################################################

args <- commandArgs(trailingOnly = TRUE)
MARKER <- args[1] # genotype data set "/data/users/mjahani/JOON_PAV/pav_gwas/SAM.5maf.noheader.sed12.rrblup.csv"
PHENOTYPE <- args[2] # phenotypes"/data/users/mjahani/JOON_PAV/pav_gwas/phenotype.csv"
SAVE_DIR <- args[3] # directory to save result "/data/users/mjahani/JOON_PAV/pav_gwas/result"
trait <- as.numeric(args[4]) # array in slurm
DONOR <- args[5]
PERMUTATION <- 10000
registerDoParallel(cores = 45)
########################################################## Read data#########################################################

Markers_impute <- fread(MARKER,
  header = F
) # read genotype file
phenotype <- read.csv(PHENOTYPE,
  header = T,
  row.names = 1
)

################################################## Fit the prediction Model##################################################

phenotype %>%
  colnames() -> TRAIT_ID

final <- data.frame(a = seq(1, as.numeric(PERMUTATION), 1), b = "PERM") %>%
  mutate(data_type = paste0(b, a)) %>%
  select(data_type) %>%
  rbind(., "observation")
P_VALUES <- data.frame(P = c("LARGER", "SMALLER"))
# for (trait in 1:length(TRAIT_ID)) {
# trait_result <- data.frame(data_type = "observation")
# observation
as.data.frame(
  as.matrix(
    mixed.solve(phenotype[, trait], # phenotyoic data
      Z = Markers_impute, # marker data
      K = NULL,
      method = "REML",
      SE = FALSE,
      return.Hinv = FALSE
    )$u
  )
) -> BLUPs # calculation of BLUPs
BLUPs %>%
  summarise(!!TRAIT_ID[trait] := mean(V1)) %>%
  mutate(data_type = "observation") %>%
  select(data_type, everything()) %>%
  fwrite(paste0(SAVE_DIR, "/Average_PERMUT_result_", DONOR, "_", TRAIT_ID[trait]),
    sep = "\t",
    col.names = F,
    quote = F,
    append = T
  )

BLUPs %>%
  mutate(sign = ifelse(V1 > 0, "pos", "not_pos")) %>%
  group_by(sign) %>%
  tally() %>%
  ungroup() %>%
  mutate(!!TRAIT_ID[trait] := n / sum(n)) %>%
  filter(sign == "pos") %>%
  mutate(data_type = "observation") %>%
  select(data_type, TRAIT_ID[trait]) %>%
  fwrite(paste0(SAVE_DIR, "/postive_prop_PERMUT_result_", DONOR, "_", TRAIT_ID[trait]),
    sep = "\t",
    col.names = F,
    quote = F,
    append = T
  )

rm(BLUPs)
# for (PERM in 1:10000) {

foreach(PERM = 1:as.numeric(PERMUTATION), .combine = "rbind") %dopar% {
  assign(paste0("BLUP", PERM), as.data.frame(
    as.matrix(
      mixed.solve(phenotype[, trait], # phenotypic data
        Z = Markers_impute[sample(nrow(Markers_impute)), ], # shuffled marker data
        K = NULL,
        method = "REML",
        SE = FALSE,
        return.Hinv = FALSE
      )$u
    )
  ))


  get(paste0("BLUP", PERM)) %>%
    summarise(!!TRAIT_ID[trait] := mean(V1)) %>%
    mutate(data_type = paste0("PERM", PERM)) %>%
    select(data_type, everything(.)) %>%
    fwrite(paste0(SAVE_DIR, "/Average_PERMUT_result_", DONOR, "_", TRAIT_ID[trait]),
      sep = "\t",
      col.names = F,
      quote = F,
      append = T
    )
  get(paste0("BLUP", PERM)) %>%
    mutate(sign = ifelse(V1 > 0, "pos", "not_pos")) %>%
    group_by(sign) %>%
    tally() %>%
    ungroup() %>%
    mutate(!!TRAIT_ID[trait] := n / sum(n)) %>%
    filter(sign == "pos") %>%
    mutate(data_type = paste0("PERM", PERM)) %>%
    select(data_type, TRAIT_ID[trait]) %>%
    fwrite(paste0(SAVE_DIR, "/postive_prop_PERMUT_result_", DONOR, "_", TRAIT_ID[trait]),
      sep = "\t",
      col.names = F,
      quote = F,
      append = T
    )

  rm(list = paste0("BLUP", PERM))
}


fread(paste0(SAVE_DIR, "/Average_PERMUT_result_", DONOR, "_", TRAIT_ID[trait]), header = F) %>%
  filter(V1 != "observation") %>%
  mutate(V2 = as.numeric(V2)) -> AVERAGE_PERMUT_result

fread(paste0(SAVE_DIR, "/Average_PERMUT_result_", DONOR, "_", TRAIT_ID[trait]), header = F) %>%
  filter(V1 == "observation") %>%
  mutate(V2 = as.numeric(V2)) -> AVERAGE_trait_result

fread(paste0(SAVE_DIR, "/postive_prop_PERMUT_result_", DONOR, "_", TRAIT_ID[trait]), header = F) %>%
  filter(V1 != "observation") %>%
  mutate(V2 = as.numeric(V2)) -> postive_prop_PERMUT_result

fread(paste0(SAVE_DIR, "/postive_prop_PERMUT_result_", DONOR, "_", TRAIT_ID[trait]), header = F) %>%
  filter(V1 == "observation") %>%
  mutate(V2 = as.numeric(V2)) -> postive_prop_trait_result


filter(AVERAGE_PERMUT_result, V2 < as.numeric(select(AVERAGE_trait_result, V2))) %>%
  summarise(P = n() / as.numeric(nrow(AVERAGE_PERMUT_result))) %>%
  mutate(P_type = "LARGER") %>%
  mutate(TRAIT = TRAIT_ID[trait]) %>%
  mutate(Significant = ifelse(P < 0.05, "Significant", "non_Significant")) %>%
  select(TRAIT, P_type, P, Significant) %>%
  rbind(
    .,
    filter(AVERAGE_PERMUT_result, V2 > as.numeric(select(AVERAGE_trait_result, V2))) %>%
      summarise(P = n() / as.numeric(nrow(AVERAGE_PERMUT_result))) %>%
      mutate(P_type = "SMALLER") %>%
      mutate(TRAIT = TRAIT_ID[trait]) %>%
      select(TRAIT, P_type, TRAIT, P) %>%
      mutate(Significant = ifelse(P < 0.05, "Significant", "non_Significant")) %>%
      select(TRAIT, P_type, P, Significant)
  ) %>%
  fwrite(paste0(SAVE_DIR, "/P_VALUES_AVERAGE_", DONOR),
    sep = "\t",
    col.names = F,
    quote = F,
    append = T
  )

filter(postive_prop_PERMUT_result, V2 < as.numeric(select(postive_prop_trait_result, V2))) %>%
  summarise(P = n() / as.numeric(nrow(postive_prop_PERMUT_result))) %>%
  mutate(P_type = "LARGER") %>%
  mutate(TRAIT = TRAIT_ID[trait]) %>%
  mutate(Significant = ifelse(P < 0.05, "Significant", "non_Significant")) %>%
  select(TRAIT, P_type, P, Significant) %>%
  rbind(
    .,
    filter(postive_prop_PERMUT_result, V2 > as.numeric(select(postive_prop_trait_result, V2))) %>%
      summarise(P = n() / as.numeric(nrow(postive_prop_PERMUT_result))) %>%
      mutate(P_type = "SMALLER") %>%
      mutate(TRAIT = TRAIT_ID[trait]) %>%
      select(TRAIT, P_type, TRAIT, P) %>%
      mutate(Significant = ifelse(P < 0.05, "Significant", "non_Significant")) %>%
      select(TRAIT, P_type, P, Significant)
  ) %>%
  fwrite(paste0(SAVE_DIR, "/P_VALUES_POST_PROP_", DONOR),
    sep = "\t",
    col.names = F,
    quote = F,
    append = T
  )



rm(
  AVERAGE_PERMUT_result,
  AVERAGE_trait_result,
  postive_prop_PERMUT_result,
  postive_prop_trait_result
)