library(rrBLUP)
library(data.table)
library(dplyr)
library(tidyr)
args <- commandArgs(trailingOnly = TRUE)

MARKER <- args[1] # "/DATA/home/mjahani/LINKADE_DRAG/new_method/GP/SAM_introgression_donor_2nd_GERMPLASM_maf0.03_rrBLUP.in"
MARKER_HOMOZ <- args[2] # "/DATA/home/mjahani/LINKADE_DRAG/new_method/GP/highly_homozygous/SAM_introgression_donor_2nd_GERMPLASM_maf0.03_234homozygous_rrBLUP.in"
PHENOTYPE <- args[3] # "/DATA/home/mjahani/LINKADE_DRAG/new_method/GP/phenotype_common_georgia_corrected.csv"
PHENOTYPE_HOMOZ <- args[4] # "/DATA/home/mjahani/LINKADE_DRAG/new_method/GP/highly_homozygous/phenotype_common_georgia_corrected_highly_homozygous.csv"
DONOR <- args[5] # "2nd_GERMPLASM"
SAVE_DIR <- args[6] # "/DATA/home/mjahani/LINKADE_DRAG/new_method/GP/highly_homozygous/result"

########################################################## Read data#########################################################

Markers_impute <- fread(MARKER,
  header = F
) # read genotype file
phenotype <- read.csv(PHENOTYPE,
  header = T,
  row.names = 1
)

Markers_homoz_impute <- fread(MARKER_HOMOZ,
  header = F
) # read genotype file
phenotype_homoz <- read.csv(PHENOTYPE_HOMOZ,
  header = T,
  row.names = 1
)

################################################## Fit the prediction Model##################################################

phenotype %>%
  colnames() -> TRAIT_ID

BLUPS <- NULL
result <- NULL

for (trait in 1:length(TRAIT_ID)) {
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
  ) %>%
    mutate(SMPLES = "full_SAM") %>%
    mutate(TRAIT = TRAIT_ID[trait]) %>%
    rename(BLUP = V1) %>%
    rbind(., BLUPS) -> BLUPS # calculation of BLUPs


  as.data.frame(
    as.matrix(
      mixed.solve(phenotype_homoz[, trait], # phenotyoic data
        Z = Markers_homoz_impute, # marker data
        K = NULL,
        method = "REML",
        SE = FALSE,
        return.Hinv = FALSE
      )$u
    )
  ) %>%
    mutate(SMPLES = "highly_homozygous234") %>%
    mutate(TRAIT = TRAIT_ID[trait]) %>%
    rename(BLUP = V1) %>%
    rbind(., BLUPS) -> BLUPS # calculation of BLUPs

  BLUPS %>%
    filter(SMPLES == "full_SAM", TRAIT == TRAIT_ID[trait]) %>%
    mutate(sign = ifelse(BLUP > 0, "pos", "not_pos")) %>%
    group_by(sign) %>%
    tally() %>%
    ungroup() %>%
    filter(sign == "pos") %>%
    pull(n) -> full_SAM_pos_frq

  BLUPS %>%
    filter(SMPLES == "highly_homozygous234", TRAIT == TRAIT_ID[trait]) %>%
    mutate(sign = ifelse(BLUP > 0, "pos", "not_pos")) %>%
    group_by(sign) %>%
    tally() %>%
    ungroup() %>%
    filter(sign == "pos") %>%
    pull(n) -> highly_homozygous234_pos_frq


  BLUPS %>%
    filter(TRAIT == TRAIT_ID[trait]) %>%
    group_by(
      TRAIT,
      SMPLES
    ) %>%
    summarize(AVERAGE = mean(BLUP)) %>%
    ungroup() %>%
    spread(SMPLES, AVERAGE) %>%
    mutate(P = (t.test(
      pull(filter(BLUPS, SMPLES == "full_SAM", TRAIT == TRAIT_ID[trait]), BLUP),
      pull(filter(BLUPS, SMPLES == "highly_homozygous234", TRAIT == TRAIT_ID[trait]), BLUP)
    )$p.value)) %>%
    mutate(STAT = "Average_effect") %>%
    select(
      TRAIT,
      STAT,
      full_SAM,
      highly_homozygous234,
      P
    ) %>%
    rbind(., result) -> result

  BLUPS %>%
    filter(TRAIT == TRAIT_ID[trait]) %>%
    group_by(
      TRAIT,
      SMPLES
    ) %>%
    summarize(POSTIVE_COUNT = sum(BLUP > 0) / n()) %>%
    ungroup() %>%
    spread(SMPLES, POSTIVE_COUNT) %>%
    mutate(P = (prop.test(
      x = c(as.numeric(full_SAM_pos_frq), as.numeric(highly_homozygous234_pos_frq)),
      n = c(as.numeric(ncol(Markers_impute)), as.numeric(ncol(Markers_homoz_impute))),
      alternative = "two.sided"
    )$p.value)) %>%
    mutate(STAT = "POSTIVE_FRQ") %>%
    select(
      TRAIT,
      STAT,
      full_SAM,
      highly_homozygous234,
      P
    ) %>%
    rbind(., result) -> result

  rm(
    full_SAM_pos_frq,
    highly_homozygous234_pos_frq
  )
}


fwrite(BLUPS,
  paste0(
    SAVE_DIR,
    "/",
    DONOR,
    "_ALL_INTROGRESSION_EFFECTS"
  ),
  col.names = T,
  sep = ",",
  quote = F
)


data.frame(
  TRAIT = c(
    "height", "head_weight", # Strong heterosis
    "stem_weight", "stem_diameter", "biomass", # Moderate heterosis
    "seed_weight", "seed_size", "leaf_weight", "leaf_area", "head_diameter", # Weak heterosis
    "stigma_antho", "oil", "leaf_sla", "dtf", "disk_antho", "branching"
  ), # No heterosis
  Heterosis = c(
    "Strong_heterosis", "Strong_heterosis",
    "Moderate_heterosis", "Moderate_heterosis", "Moderate_heterosis",
    "Weak_heterosis", "Weak_heterosis", "Weak_heterosis", "Weak_heterosis", "Weak_heterosis",
    "No_heterosis", "No_heterosis", "No_heterosis", "No_heterosis", "No_heterosis", "No_heterosis"
  )
) -> TRAIT_HET


result %>%
  mutate(TRAIT = gsub("plant_biomass", "biomass", TRAIT)) %>%
  mutate(TRAIT = gsub("plant_height", "height", TRAIT)) %>%
  mutate(TRAIT = gsub("seed_lxw", "seed_size", TRAIT)) %>%
  separate(TRAIT, into = c("TRAIT", "LOCATION"), sep = "__", remove = T) %>%
  mutate(LOCATION = gsub("iowa", "IA", LOCATION)) %>%
  mutate(LOCATION = gsub("georgia", "GA", LOCATION)) %>%
  mutate(LOCATION = gsub("UBC", "BC", LOCATION)) %>%
  mutate(FULL_SAM_VS_HOMOZYGOUS = ifelse((as.numeric(full_SAM) > as.numeric(highly_homozygous234)) & P < 0.05,
    "GRATER_SIGINIFICANT",
    ifelse((as.numeric(full_SAM) < as.numeric(highly_homozygous234)) & P < 0.05,
      "SMALLER_SIGINIFICANT",
      "NOT_SIGNIFICANT"
    )
  )) %>%
  full_join(., TRAIT_HET) %>%
  select(
    TRAIT,
    LOCATION,
    Heterosis,
    STAT,
    full_SAM,
    highly_homozygous234,
    P,
    FULL_SAM_VS_HOMOZYGOUS
  ) %>%
  fwrite(paste0(
    SAVE_DIR,
    "/",
    DONOR,
    "_het_VS_hom_introgression_effect_postive_prop"
  ),
  col.names = T,
  sep = ",",
  quote = F
  )