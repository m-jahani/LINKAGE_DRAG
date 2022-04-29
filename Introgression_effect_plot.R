library(data.table)
library(tidyverse)
library(plyr)

PATH <- "/Users/mojtabajahani/Documents/Projects/Linkage_drag/New_method/GP/introgression_impact_new"
Pvalue_annuus <- "/Users/mojtabajahani/Documents/Projects/Linkage_drag/New_method/GP/introgression_impact_new/P_VALUES_AVERAGE_WILD_ANNUUS"
Pvalue_2nd <- "/Users/mojtabajahani/Documents/Projects/Linkage_drag/New_method/GP/introgression_impact_new/P_VALUES_AVERAGE_Secondary_Germplasm"

list.files(path = PATH) %>%
  as.data.frame() %>% 
  mutate(FILE=as.character(.)) %>% 
  select(FILE) %>% 
  filter(!grepl("P_VALUE",FILE)) %>% 
  separate(FILE,into = c("START","LOCATION"),sep = "__",remove = F) %>% 
  mutate(LOCATION = gsub("_z","",LOCATION)) %>% 
  separate(START,into = c("DATA","DONOR_TRAIT"),sep = "_PERMUT_result_",remove = T) %>% 
  filter(DATA == "Average") %>% 
  select(-DATA) %>%
  mutate(DONOR_TRAIT=gsub("WILD_ANNUUS","WILD_ANNUUS_",DONOR_TRAIT)) %>%
  mutate(DONOR_TRAIT=gsub("SECONDARY_GERMPLASM","SECONDARY_GERMPLASM_",DONOR_TRAIT)) %>% 
  separate(DONOR_TRAIT, into = c("DONOR","TRAIT"), sep = "__", remove = T) %>% 
  mutate(TRAIT= gsub("plant_biomass","biomass",TRAIT)) %>% 
  mutate(TRAIT= gsub("plant_height","height",TRAIT)) %>%
  mutate(TRAIT= gsub("seed_lxw","seed_size",TRAIT)) %>%
  mutate(LOCATION = gsub("iowa","Iowa",LOCATION)) %>%
  mutate(LOCATION= gsub("georgia","Georgia",LOCATION)) %>%
  mutate(LOCATION= gsub("UBC","BC",LOCATION)) -> DATA_MAP

fread(Pvalue_annuus) %>% 
  separate(V1, into = c("TRAIT","LOCATION"), sep = "__", remove = T) %>% 
  mutate(LOCATION = gsub("_z","",LOCATION)) %>% 
  mutate(DONOR = "WILD_ANNUUS") %>% 
  select(DONOR,
         TRAIT,
         LOCATION,
         TEST = V4,
         DIRECTION = V2) %>% 
  filter(TEST == "Significant") -> WILD_ANNUUS_TEST
  
fread(Pvalue_2nd) %>% 
  separate(V1, into = c("TRAIT","LOCATION"), sep = "__", remove = T) %>%
  mutate(LOCATION = gsub("_z","",LOCATION)) %>% 
  mutate(DONOR = "SECONDARY_GERMPLASM") %>% 
  select(DONOR,
         TRAIT,
         LOCATION,
         TEST = V4,
         DIRECTION = V2) %>% 
  filter(TEST == "Significant") %>% 
  rbind(.,WILD_ANNUUS_TEST) %>% 
  mutate(TRAIT= gsub("plant_biomass","biomass",TRAIT)) %>%
  mutate(TRAIT= gsub("plant_height","height",TRAIT)) %>%
  mutate(TRAIT= gsub("seed_lxw","seed_size",TRAIT)) %>%
  mutate(LOCATION = gsub("iowa","Iowa",LOCATION)) %>%
  mutate(LOCATION= gsub("georgia","Georgia",LOCATION)) %>%
  mutate(LOCATION= gsub("UBC","BC",LOCATION)) -> test
 
 rm(WILD_ANNUUS_TEST)
# 
# Result <- matrix(nrow=nrow(DATA_MAP),ncol = 5)
# 
# for (i in 1:nrow(DATA_MAP)) {
#   Result[i,1] <- DATA_MAP[i,2]
#   Result[i,2] <- DATA_MAP[i,3]
#   Result[i,3] <- DATA_MAP[i,4]
#   
#   fread(paste0(PATH,
#                "/",
#                as.character(DATA_MAP[i,1]))) -> tmp_data
# 
#   tmp_data %>%
#     filter(V1 == "observation") %>%
#     pull(V2) -> Result[i,4]
#   
#   tmp_data %>% 
#     filter(V1 != "observation") %>% 
#     summarize(Average = mean(V2)) %>%
#     pull(Average) -> Result[i,5]
#   
#   rm(tmp_data)
# }
#  
# Result %>% 
#   as.data.frame() %>% 
#   select(DONOR = V1,
#          TRAIT = V2,
#          LOCATION = V3,
#          Observation = V4,
#          Permutation = V5) %>% 
#   left_join(.,test)  %>% 
#   mutate(TEST = ifelse(is.na(TEST),"Non_Significant",TEST)) %>%
#   mutate(DIRECTION = ifelse(is.na(DIRECTION),"EQUAL",DIRECTION)) %>%
#   gather(TYPE,VALUE,Observation,Permutation) %>%
#   mutate(VALUE=as.numeric(VALUE)) %>%
#   mutate(TRAIT_LOC = paste0(TRAIT,"_",LOCATION)) -> FINAL_DATA

# FINAL_DATA %>% 
#   mutate(COLOR=ifelse(TYPE == "Permutation","gray",
#                 ifelse(TEST == "Non_Significant","yellow",
#                        ifelse(DIRECTION == "LARGER","blue",
#                               ifelse(DIRECTION == "SMALLER", "red","yellow"))))) %>%
#   ggplot(.,aes(x = TRAIT_LOC,
#                y = VALUE,
#                color = COLOR,
#                shape = TYPE)) +
#   geom_point() +
#   theme_cowplot() +
#   theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1)) +
#   scale_shape_manual(values = c(20,1),
#                      name = "Data type",
#                      labels = c("Observation","Permutation (Average)")) +
#   scale_color_manual(values = c("black","blue","gray","red"),
#                      name = "TEST",
#                      labels = c("Non_significant","Significantly Larger","Permutation", "Significantly Smaller")) +
#   ylab("Standard Average effect size of introgressions") +
#   xlab("Trait-Location") 


DATA_MAP %>% 
  mutate(TRAIT_LOC = paste0(TRAIT,"_",LOCATION)) -> DATA_PLOT


all_point_data <- NULL
for (i in 1:nrow(DATA_PLOT)) {
  fread(paste0(PATH,"/",as.character(DATA_PLOT[i,1]))) %>% 
    mutate(TYPE = ifelse(V1 == "observation","Observation","Permutation")) %>% 
    mutate(DONOR = DATA_PLOT[i,2],
           TRAIT = DATA_PLOT[i,3],
           LOCATION = DATA_PLOT[i,4],
           TRAIT_LOC = DATA_PLOT[i,5]) %>%
    select(DONOR,
           TRAIT,
           LOCATION,
           TRAIT_LOC,
           TYPE,
           VALUE = V2) %>% 
    rbind(.,all_point_data) -> all_point_data
}
  
all_point_data %>% 
  left_join(.,test) %>% 
  mutate(TEST = ifelse(is.na(TEST),"Non_Significant",TEST)) %>%
  mutate(DIRECTION = ifelse(is.na(DIRECTION),"EQUAL",DIRECTION)) %>% 
  mutate(COLOR=ifelse(TYPE == "Permutation","gray",
                      ifelse(TEST == "Non_Significant","black",
                             ifelse(DIRECTION == "LARGER","blue",
                                    ifelse(DIRECTION == "SMALLER", "red","yellow")))))  -> all_point_data

ggplot(all_point_data) +
  geom_jitter(data = filter(all_point_data,TYPE != "Observation"),
              aes(x = TRAIT_LOC, y = VALUE, color = COLOR ),
              shape = 1,
              size = 2,
              position = position_jitter(width = 0.21)) +

  geom_point(data = filter(all_point_data,TYPE == "Observation"),
             aes(x = TRAIT_LOC, y = VALUE, color = COLOR, shape = TYPE ),
                          shape = 20,
                          size = 5) +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        plot.title = element_text(hjust = 0.5),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.position='right',
        text = element_text(size = 20)
  ) +
  scale_color_manual(values = c("yellow","blue","gray","red"),
                     name = "",
                     labels = c("No significant difference Observation","Significantly Larger Observation","Null Distribution", "Significantly Smaller Observation")) +
  ylab("Standard Average effect size of introgressions") +
  xlab("Trait-Location") +
  facet_wrap(~DONOR,nrow=2,scales = "free_y")






