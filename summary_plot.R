library(data.table)
library(tidyverse)
library(cowplot)
library(PNWColors)
library(ggbeeswarm)



fread("/Users/mojtabajahani/Documents/Projects/Linkage_drag/permutation_introg_effect/BLUP_average_permutation_and_realdata_result") %>% 
  gather(TRAIT,VALUE,1:46) %>%
  mutate(TRAIT=gsub("plant_biomass","biomass",TRAIT)) %>% 
  mutate(TRAIT=gsub("plant_height","height",TRAIT)) %>%
  mutate(TRAIT=gsub("seed_lxw","seed_size",TRAIT)) %>% 
  rename(INTROGRESSIONS = INTREGRESSION) %>%
  mutate(INTROGRESSIONS=gsub("other.","Secondary_Germplasm_",INTROGRESSIONS)) %>%
  mutate(INTROGRESSIONS=gsub("frq","Frequency>",INTROGRESSIONS)) %>%
  mutate(INTROGRESSIONS=gsub("annuus.","Annuus_",INTROGRESSIONS)) -> DATA 
 
DATA %>% distinct(INTROGRESSIONS,TRAIT,TYPE,SOURCE) -> ING_TRAIT


P_VALUE <- NULL
for (i in 1:nrow(ING_TRAIT)) {
  
  as.numeric(pull(filter(DATA,
                         INTROGRESSIONS == as.character(ING_TRAIT[i,1]),
                         TRAIT ==  as.character(ING_TRAIT[i,4]),
                         DATA == "REAL"),
                  VALUE)) -> REAL_VALUE
  
  DATA %>%  
    filter(INTROGRESSIONS == as.character(ING_TRAIT[i,1]),
           TRAIT ==  as.character(ING_TRAIT[i,4]),
           !DATA == "REAL") %>% 
    filter(VALUE >  REAL_VALUE ) %>% 
    summarise(P=n()/10000) %>%
    mutate(INTROGRESSIONS = as.character(ING_TRAIT[i,1]),
           TRAIT = as.character(ING_TRAIT[i,4]),
           TYPE = as.character(ING_TRAIT[i,2]),
           SOURCE = as.character(ING_TRAIT[i,3]),
           TEST=ifelse(P<=0.05,"Significantly_Larger",
                       ifelse(P>=0.95,"Significantly_Smaller","No_Significant_Difference"))) %>%
    separate(TRAIT, into = c("TRAITS","LOCATION"), sep = "__", remove = T) %>%
    rbind(.,P_VALUE) -> P_VALUE
    
  rm(REAL_VALUE)
}

DATA %>% separate(TRAIT, into = c("TRAITS","LOCATION"), sep = "__", remove = T) -> DATA
  
DATA %>% filter(DATA == "REAL") %>% full_join(.,P_VALUE) -> P_VALUE

# DATA %>% 
#   mutate(DATA = gsub("PERM\\d+","PERM",DATA)) %>%  
#   group_by(INTROGRESSIONS,DATA,TYPE,SOURCE,TRAIT) %>% 
#   summarise(MEAN=mean(VALUE)) %>%
#   ungroup() %>% 
#   full_join(.,P_VALUE) %>% 
#   separate(TRAIT, into=c("TRAITS","LOCATION"), sep = "__",remove = T) %>% 
#   mutate(TEST=ifelse(P<=0.05,"Significantly_Larger",ifelse(P>=0.95,"Significantly_Smaller","No_Significant_Difference"))) -> RESULT
# 
# rm(ING_TRAIT,
#    P_VALUE,
#    i)


# colours <- pnw_palette("Bay",5)
# 
# 
# 
# order <- c("branching","disk_antho","dtf","leaf_sla", "oil", "stigma_antho",
#            "head_diameter","leaf_area","leaf_weight", "seed_size","seed_weight",
#            "biomass","stem_diameter","stem_weight","head_weight","height")
# 
# 
# RESULT %>% 
#   #filter(!TRAITS %in% c("leaf_sla","leaf_area")) %>% 
#   mutate(TRAITS = gsub("_"," ",TRAITS)) %>%
#   mutate(TRAITS = gsub("leaf sla","SLA",TRAITS)) %>%
#   mutate(TRAITS = gsub("dtf","DTF",TRAITS)) %>%
#   mutate(LOCATION = gsub("ubc","UBC",LOCATION)) %>%
#   mutate(LOCATION = gsub("georgia","Georgia",LOCATION)) %>%
#   mutate(LOCATION = gsub("iowa","Iowa",LOCATION)) %>%
#   filter(INTROGRESSIONS == "other.frq0.03") -> plot_DATA #,DATA=="REAL"
# 
# 
# plot_DATA%>%
#   ggplot(.,aes(x=LOCATION,y=MEAN,color=DATA,shape=TEST)) + #fct_relevel(TRAITS,order)
#     geom_point(size=3) +
#     facet_wrap(~TRAITS) +
#     theme_cowplot() +
#     theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
#     xlab("Trait") +
#     ylab("BLUP") +
#     scale_color_manual(values=pnw_palette("Bay",5)[c(1,3)],
#                        name="Data",
#                        labels=c("Observrd BLUPs effect","Average of 10000 Permutation")) +
#     annotate("rect",xmin=0.5,xmax=6.5,ymin=min(plot_DATA$MEAN),ymax=max(plot_DATA$MEAN),fill=pnw_palette("Winter",4)[4],alpha=0.3) +
#     annotate("rect",xmin=6.5,xmax=11.5,ymin=min(plot_DATA$MEAN),ymax=max(plot_DATA$MEAN),fill=pnw_palette("Winter",4)[3],alpha=0.3) +
#     annotate("rect",xmin=11.5,xmax=14.5,ymin=min(plot_DATA$MEAN),ymax=max(plot_DATA$MEAN),fill=pnw_palette("Winter",4)[2],alpha=0.3) +
#     annotate("rect",xmin=14.5,xmax=16.5,ymin=min(plot_DATA$MEAN),ymax=max(plot_DATA$MEAN),fill=pnw_palette("Winter",4)[1],alpha=0.3) 


DATA %>% distinct(SOURCE,INTROGRESSIONS) -> INTRO

for (i in 1:nrow(INTRO)) {
  DATA %>% 
    filter(INTROGRESSIONS == INTRO[i,1],
           DATA != "REAL") -> PERMUTATION_POINTS
  
  P_VALUE %>% 
    filter(INTROGRESSIONS == INTRO[i,1],
           DATA == "REAL") -> ACTUAL_POINTS
  
  PERMUTATION_POINTS %>%
    ggplot(.,aes(x=LOCATION,y=VALUE)) +
    geom_violin() +
    geom_point(data = ACTUAL_POINTS, size=2, aes(y = VALUE,color = TEST)) +
    facet_wrap(~TRAITS) + 
    theme_cowplot() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    xlab("Trait") +
    ylab(INTRO[i,2]) +
    ggtitle(INTRO[i,1]) +
    scale_color_manual(values = pnw_palette("Bay",8)[c(4,6,8)],
                       name="Test",
                       labels=c("No Significant difference","Significantly larger","Significantly smaller"))
    ggsave(paste0("/Users/mojtabajahani/Documents/Projects/Linkage_drag/permutation_introg_effect/",INTRO[i,2],"_",INTRO[i,1],".pdf"))
  
  rm(PERMUTATION_POINTS,ACTUAL_POINTS)
  
}

P_VALUE %>% 
  separate(INTROGRESSIONS,into = c("donor","frequency"),sep = "_Frequency", remove = T) %>%
  #filter(frequency == ">0.05") %>%
  select(donor,
         frequency,
         TRAITS,
         LOCATION,
         TEST) %>%
  ggplot(.,aes(x=frequency,y="TEST",shape=TEST,color=LOCATION)) +
  geom_point()+
  #geom_quasirandom() +
  geom_beeswarm(size=5,dodge.width=0.5) + 
  facet_grid(TRAITS~donor) +
  theme_cowplot()




pl_data %>%
  ggplot(.,aes(x=LOCATION,y=VALUE)) + #fct_relevel(TRAITS,order)
  geom_violin() +
  geom_point(data = ACTUAL_POINTS, size=2, aes(y=MEAN,color=TEST)) +
  facet_wrap(~TRAITS) + 
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  xlab("Trait") +
  ylab("BLUP") +
  ggtitle("other.frq0.03")+
  scale_color_manual(values=pnw_palette("Bay",8)[c(4,6,8)],
                     name="Test",
                     labels=c("No Significant difference","Significantly larger","Significantly smaller"))



  filter(pl_data,INTROGRESSIONS=="other.frq0.03",DATA=="REAL")

  min(plot_DATA$MEAN)

beta_GP_simulation %>% 

  ggplot(.,aes(x=fct_relevel(trait,order),y=beta1,color=location)) +
  geom_point(size=3) +
  geom_hline(yintercept=0,linetype="dotted") +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  xlab("Trait") +
  ylab("Standard Beta Coefficient (slope)") +
  scale_color_manual(values=pnw_palette("Bay",5)[c(1,3,5)],
                     name="Location",
                     labels=c("Georgia","Iowa","BC")) +
  annotate("rect",xmin=0.5,xmax=6.5,ymin=-1,ymax=1.5,fill=pnw_palette("Winter",4)[4],alpha=0.3) +
  annotate("rect",xmin=6.5,xmax=11.5,ymin=-1,ymax=1.5,fill=pnw_palette("Winter",4)[3],alpha=0.3) +
  annotate("rect",xmin=11.5,xmax=14.5,ymin=-1,ymax=1.5,fill=pnw_palette("Winter",4)[2],alpha=0.3) +
  annotate("rect",xmin=14.5,xmax=16.5,ymin=-1,ymax=1.5,fill=pnw_palette("Winter",4)[1],alpha=0.3) +
  geom_point(size=3) +
  coord_cartesian(ylim=c(-0.025,0.025)) 

  
