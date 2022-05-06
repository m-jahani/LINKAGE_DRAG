library(data.table)
library(tidyverse)


PATH <- "/Users/mojtabajahani/Documents/Projects/Linkage_drag/New_method/GP/introgression_impact_new"
#Pvalue_annuus <- "/Users/mojtabajahani/Documents/Projects/Linkage_drag/New_method/GP/introgression_impact_new/P_VALUES_AVERAGE_WILD_ANNUUS"
#Pvalue_2nd <- "/Users/mojtabajahani/Documents/Projects/Linkage_drag/New_method/GP/introgression_impact_new/P_VALUES_AVERAGE_Secondary_Germplasm"



list.files(path = PATH) %>%
  as.data.frame() %>% 
  mutate(FILE=as.character(.)) %>% 
  select(FILE) %>% 
  filter(grepl("Average_PERMUT_result",FILE)) %>% 
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
  mutate(LOCATION = gsub("iowa","IA",LOCATION)) %>%
  mutate(LOCATION= gsub("georgia","GA",LOCATION)) %>%
  mutate(LOCATION= gsub("UBC","BC",LOCATION)) %>%
  mutate(TRAIT_LOC = paste0(TRAIT,"_",LOCATION)) -> DATA_MAP



all_point_data <- NULL
for (i in 1:nrow(DATA_MAP)) {
  fread(paste0(PATH,"/",as.character(DATA_MAP[i,1]))) %>% 
    filter(V1 == "observation") %>% 
    mutate(DONOR = DATA_MAP[i,2],
           TRAIT = DATA_MAP[i,3],
           LOCATION = DATA_MAP[i,4],
           TRAIT_LOC = DATA_MAP[i,5]) %>%
    select(DONOR,
           TRAIT,
           LOCATION,
           TRAIT_LOC,
           VALUE = V2) %>% 
    rbind(.,all_point_data) -> all_point_data
}

min(all_point_data$VALUE)

all_point_data %>% 
  ggplot(.,aes(x=TRAIT_LOC,y=VALUE,color=DONOR)) +
  geom_point(size=3) +
  scale_color_manual(values = c("red","blue"),
                     name = "Donor",
                     labels = c("Secondary gemplasm","Wild annuus")
                    ) +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        plot.title = element_text(hjust = 0.5),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.position='right',
        text = element_text(size = 15),
  ) +
  annotate("rect", xmin=0, xmax=3.5, ymin=-Inf, ymax=+Inf, alpha=0.08, fill="yellow") +
  annotate("rect", xmin=6.5, xmax=9.5, ymin=-Inf, ymax=+Inf, alpha=0.08, fill="yellow") +
  annotate("rect", xmin=12.5, xmax=15.5, ymin=-Inf, ymax=+Inf, alpha=0.08, fill="yellow") +
  annotate("rect", xmin=18.5, xmax=21.5, ymin=-Inf, ymax=+Inf, alpha=0.08, fill="yellow") +
  annotate("rect", xmin=23.5, xmax=25.5, ymin=-Inf, ymax=+Inf, alpha=0.08, fill="yellow") +
  annotate("rect", xmin=28.5, xmax=31.5, ymin=-Inf, ymax=+Inf, alpha=0.08, fill="yellow") +
  annotate("rect", xmin=34.5, xmax=37.5, ymin=-Inf, ymax=+Inf, alpha=0.08, fill="yellow") +
  annotate("rect", xmin=40.5, xmax=43.5, ymin=-Inf, ymax=+Inf, alpha=0.08, fill="yellow") +
  ylab("Standard Average effect size of introgressions") 

  

