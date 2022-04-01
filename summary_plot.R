library(data.table)
library(tidyverse)
library(cowplot)

fread("/Users/mojtabajahani/Documents/Projects/Linkage_drag/permutation_introg_effect/BLUP_average_permutation_and_realdata_result") %>% 
  gather(TRAIT,VALUE,1:46) %>%
  mutate(TRAIT=gsub("plant_biomass","biomass",TRAIT)) %>% 
  mutate(TRAIT=gsub("plant_height","height",TRAIT)) %>%
  mutate(TRAIT=gsub("seed_lxw","seed_size",TRAIT))-> DATA 
 
DATA %>% distinct(INTREGRESSION,TRAIT,TYPE,SOURCE) -> ING_TRAIT


P_VALUE <- NULL
for (i in 1:nrow(ING_TRAIT)) {
  
  as.numeric(pull(filter(DATA,
                         INTREGRESSION == as.character(ING_TRAIT[i,1]),
                         TRAIT ==  as.character(ING_TRAIT[i,4]),
                         DATA == "REAL"),
                  VALUE)) -> REAL_VALUE
  
  DATA %>%  
    filter(INTREGRESSION == as.character(ING_TRAIT[i,1]),
           TRAIT ==  as.character(ING_TRAIT[i,4]),
           !DATA == "REAL") %>% 
    filter(VALUE >  REAL_VALUE ) %>% 
    summarise(P=n()/10000) %>%
    mutate(INTREGRESSION = as.character(ING_TRAIT[i,1]),
           TRAIT = as.character(ING_TRAIT[i,4]),
           TYPE = as.character(ING_TRAIT[i,2]),
           SOURCE = as.character(ING_TRAIT[i,3])) %>%
    rbind(.,P_VALUE) -> P_VALUE
    
  rm(REAL_VALUE)
}


DATA %>% 
  mutate(DATA = gsub("PERM\\d+","PERM",DATA)) %>%  
  group_by(INTREGRESSION,DATA,TYPE,SOURCE,TRAIT) %>% 
  summarise(MEAN=mean(VALUE)) %>%
  ungroup() %>% 
  full_join(.,P_VALUE) %>% 
  separate(TRAIT, into=c("TRAITS","LOCATION"), sep = "__",remove = T) %>% 
  mutate(TEST=ifelse(P<=0.05,"Significantly_Larger",ifelse(P>=0.95,"Significantly_Smaller","No_Significant_Difference"))) -> RESULT

rm(DATA,
   ING_TRAIT,
   P_VALUE,
   i)



library(cowplot)
library(PNWColors)

colours <- pnw_palette("Bay",5)



order <- c("branching","disk_antho","dtf","leaf_sla", "oil", "stigma_antho",
           "head_diameter","leaf_area","leaf_weight", "seed_size","seed_weight",
           "biomass","stem_diameter","stem_weight","head_weight","height")
RESULT %>% distinct(TRAITS,LOCATION) %>% mutate(b=1) %>% spread(LOCATION,b)-> a


RESULT %>%
  mutate(TRAITS = gsub("_"," ",TRAITS)) %>%
  mutate(TRAITS = gsub("leaf sla","SLA",TRAITS)) %>%
  mutate(TRAITS = gsub("dtf","DTF",TRAITS)) %>%
  filter(INTREGRESSION == "other.frq0.03",DATA=="REAL") -> plot_DATA
plot_DATA%>%
  ggplot(.,aes(x=fct_relevel(TRAITS,order),y=MEAN,color=LOCATION,shape=TEST)) +
  #geom_boxplot(outlier.shape = NA) +
  #geom_jitter(width = 0.3 ,  alpha= 0.7, color= "#E0E0E0",fill=NA, size = 0.9 )
    geom_point(size=3) +
    geom_hline(yintercept=0,linetype="dotted") +
    theme_cowplot() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    xlab("Trait") +
    ylab("BLUP") +
    scale_color_manual(values=pnw_palette("Bay",5)[c(1,3,5)],
                       name="Location",
                       labels=c("Georgia","Iowa","BC")) +
    annotate("rect",xmin=0.5,xmax=6.5,ymin=min(plot_DATA$MEAN),ymax=max(plot_DATA$MEAN),fill=pnw_palette("Winter",4)[4],alpha=0.3) +
    annotate("rect",xmin=6.5,xmax=11.5,ymin=min(plot_DATA$MEAN),ymax=max(plot_DATA$MEAN),fill=pnw_palette("Winter",4)[3],alpha=0.3) +
    annotate("rect",xmin=11.5,xmax=14.5,ymin=min(plot_DATA$MEAN),ymax=max(plot_DATA$MEAN),fill=pnw_palette("Winter",4)[2],alpha=0.3) +
    annotate("rect",xmin=14.5,xmax=16.5,ymin=min(plot_DATA$MEAN),ymax=max(plot_DATA$MEAN),fill=pnw_palette("Winter",4)[1],alpha=0.3) 


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

  
