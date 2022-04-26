library(data.table)
library(tidyverse)
library(cowplot)
library(PNWColors)

#calculates Zscore of blups and draw plots
args = commandArgs(trailingOnly = TRUE)
INTRO_FRQ_ANNUUS <- args[1] #"/Users/mojtabajahani/Documents/Projects/Linkage_drag/New_method/GP/FRQ_VS_BLUP_NEW/WILD_ANNUUS_ALL_INTROGRESSION_EFFECTS_FRQUENCY_ID"
INTRO_FRQ_Secondry <- args[2] #"/Users/mojtabajahani/Documents/Projects/Linkage_drag/New_method/GP/FRQ_VS_BLUP_NEW/Secondary_Germplasm_ALL_INTROGRESSION_EFFECTS_FRQUENCY_ID"
SAVE_DIR <- args[3] #"/Users/mojtabajahani/Documents/Projects/Linkage_drag/New_method/GP/FRQ_VS_BLUP_NEW/"

args = commandArgs(trailingOnly = TRUE)
INTRO_FRQ_ANNUUS <- "/Users/mojtabajahani/Documents/Projects/Linkage_drag/New_method/GP/FRQ_VS_BLUP_NEW/WILD_ANNUUS_ALL_INTROGRESSION_EFFECTS_FRQUENCY_ID"
INTRO_FRQ_Secondry <- "/Users/mojtabajahani/Documents/Projects/Linkage_drag/New_method/GP/FRQ_VS_BLUP_NEW/Secondary_Germplasm_ALL_INTROGRESSION_EFFECTS_FRQUENCY_ID"
SAVE_DIR <- "/Users/mojtabajahani/Documents/Projects/Linkage_drag/New_method/GP/FRQ_VS_BLUP_NEW"


fread(INTRO_FRQ_ANNUUS) %>% 
  mutate(DONOR="Wild_Annuus") %>%
  select(DONOR,
         TRAIT,
         LOCATION,
         BLUP,
         FRQ) -> INTRO_FRQ_ANNUUS_DATA

fread(INTRO_FRQ_Secondry) %>% 
  mutate(DONOR="Secondary_Germplasm") %>%
  select(DONOR,
         TRAIT,
         LOCATION,
         BLUP,
         FRQ) -> INTRO_FRQ_Secondry_DATA

rbind(INTRO_FRQ_ANNUUS_DATA,
      INTRO_FRQ_Secondry_DATA) ->  INTRO_FRQ_DATA

rm(INTRO_FRQ_ANNUUS_DATA,
      INTRO_FRQ_Secondry_DATA)

INTRO_FRQ_DATA %>%
  distinct(DONOR,
           TRAIT,
           LOCATION) -> DONOR_TRAIT_LOCA

SCALED_BLUP <- NULL
for (i in 1:nrow(DONOR_TRAIT_LOCA)) {
  INTRO_FRQ_DATA %>%
    filter(DONOR == as.character(DONOR_TRAIT_LOCA[i,1]),
           TRAIT == as.character(DONOR_TRAIT_LOCA[i,2]),
           LOCATION == as.character(DONOR_TRAIT_LOCA[i,3])) %>% 
    mutate_at(vars(BLUP), list(ZBLUP = ~as.vector(scale(.,center = T,scale = T)))) %>%
    rbind(.,SCALED_BLUP) -> SCALED_BLUP
}
  
mutate(trait_loc = paste0(TRAIT," (",LOCATION,")"))

SCALED_BLUP %>%
  filter(DONOR == "Wild_Annuus") %>%
  ggplot(., aes(FRQ, ZBLUP,color = LOCATION)) +
  geom_smooth(method = "lm",
              se = F) +
  xlab("Frequency of introgression") +
  ylab("Effect size of introgression (Z score)") +
  geom_hline(yintercept=0,linetype="dotted") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        plot.title = element_text(hjust = 0.5),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.position='right'
  ) +
  scale_color_manual(values=pnw_palette("Bay",5)[c(1,3,5)],
                     name="Location",
                     labels=c("Georgia","Iowa","BC"))+
  facet_wrap(~TRAIT) +
  ggtitle("Introgressions from Wild annuus")

ggsave(paste0(SAVE_DIR,"/Wild_Annuus_SLOPE_LINE_PLOT.PDF"))


SCALED_BLUP %>%
  filter(DONOR == "Secondary_Germplasm") %>%
  ggplot(., aes(FRQ, ZBLUP,color = LOCATION)) +
  geom_smooth(method = "lm",
              se = F) +
  xlab("Frequency of introgression") +
  ylab("Effect size of introgression (Z score)") +
  geom_hline(yintercept=0,linetype="dotted") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        plot.title = element_text(hjust = 0.5),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.position='right'
  ) +
  scale_color_manual(values=pnw_palette("Bay",5)[c(1,3,5)],
                     name="Location",
                     labels=c("Georgia","Iowa","BC"))+
  facet_wrap(~TRAIT) +
  ggtitle("Introgressions from Secondary germplasm")

ggsave(paste0(SAVE_DIR,"/Secondary_Gemplasm_SLOPE_LINE_PLOT.PDF"))


RESULT <- matrix(nrow = nrow(DONOR_TRAIT_LOCA),
                 ncol = 6)

for (trait in 1:nrow(DONOR_TRAIT_LOCA)) {
  
  FRQU <- pull(filter(SCALED_BLUP,DONOR == as.character(DONOR_TRAIT_LOCA[trait,1]),
                      TRAIT == as.character(DONOR_TRAIT_LOCA[trait,2]),
                      LOCATION == as.character(DONOR_TRAIT_LOCA[trait,3])),FRQ)
  
  ZBLUPs <- pull(filter(SCALED_BLUP,DONOR == as.character(DONOR_TRAIT_LOCA[trait,1]),
                        TRAIT == as.character(DONOR_TRAIT_LOCA[trait,2]),
                        LOCATION == as.character(DONOR_TRAIT_LOCA[trait,3])),ZBLUP)
  fit=lm(FRQU ~ ZBLUPs)
  RESULT[trait,1] <-   as.character(DONOR_TRAIT_LOCA[trait,1]) #DONOT ID
  RESULT[trait,2] <-   as.character(DONOR_TRAIT_LOCA[trait,2]) #trait id
  RESULT[trait,3] <-   as.character(DONOR_TRAIT_LOCA[trait,3]) #location id
  RESULT[trait,4] <- summary(fit)$r.squared #R-squared
  RESULT[trait,5] <- summary(fit)$coefficients[2,1] #Beta1 (slope)
  RESULT[trait,6] <- summary(fit)$coefficients[2,4] #Pvalue of Beta 1
  rm(FRQU,
     ZBLUPs)
}

RESULT %>%
  as.data.frame() %>%
  mutate(R_SQUARED = as.numeric(V4),
         SLOPE = as.numeric(V5),
         P_SLOPE = as.numeric(V6)) %>%
  mutate(TEST=ifelse(P_SLOPE <=  0.05,"Significant", "Non_Significant")) %>%
  select(DONOR = V1,
         TRAIT = V2,
         LOCATION = V3,
         R_SQUARED,
         SLOPE,
         P_SLOPE,
         TEST) -> RESULT

RESULT %>% 
  ggplot(.,aes(x=TRAIT,
               y=SLOPE,
               color=LOCATION,
               shape=DONOR)) +
  geom_point(size=5) +
  geom_hline(yintercept=0,linetype="dotted") +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        text = element_text(size=20)) +
  xlab("Trait") +
  ylab("Standard Beta Coefficient (slope)") +
  scale_color_manual(values=pnw_palette("Bay",5)[c(1,3,5)],
                     name="Location",
                     labels=c("Georgia","Iowa","BC")) +
  scale_shape_manual(values = c(8, 19),
                     name = "Donor",
                     labels = c("Secondary Germplasm","Wild Annuus")) +
  coord_cartesian(ylim=c(
    (max(abs(RESULT$SLOPE))+0.01)*-1,
    max(abs(RESULT$SLOPE))+0.01)) 
  
    
ggsave(paste0(SAVE_DIR,"/SLOPE_POINT_PLOT.PDF"), 
       width = 20,
       height = 7)



SCALED_BLUP %>% 
  filter(TRAIT %in% c("branching","oil","dtf"),
         DONOR=="Wild_Annuus") %>% arrange(desc(FRQ)) %>%
  transform(TRAIT=factor(TRAIT,levels=c("branching","oil","dtf"))) %>%
  ggplot(.,aes(FRQ, 
               ZBLUP,
               color = LOCATION)) +
  geom_smooth(method = "lm",
              se = F) +
  xlab("Frequancy of introgression") +
  ylab("Effect size of introgression (Z score)") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        plot.title = element_text(hjust = 0.5),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.position='none',
        text = element_text(size=20)
  ) +
    xlim(0, 0.5)+
  scale_color_manual(values=pnw_palette("Bay",5)[c(1,3,5)],
                     name="Location",
                     labels=c("Georgia","Iowa","BC")) +
  facet_wrap(~TRAIT)
  
  
  ggsave(paste0(SAVE_DIR,"/Wild_Annuus_SLOPE_POINT_PLOT.PDF"), 
         width = 10,
         height = 3)
    
  
  SCALED_BLUP %>% 
    filter(TRAIT %in% c("branching","oil","dtf"),
           DONOR=="Secondary_Germplasm") %>% arrange(desc(FRQ)) %>%
  transform(TRAIT=factor(TRAIT,levels=c("branching","oil","dtf"))) %>%
    ggplot(.,aes(FRQ, 
                 ZBLUP,
                 color = LOCATION)) +
    geom_smooth(method = "lm",
                se = F) +
    xlab("Frequancy of introgression") +
    ylab("Effect size of introgression (Z score)") +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA, size=0.5),
          plot.title = element_text(hjust = 0.5),
          strip.background = element_blank(),
          strip.placement = "outside",
          legend.position='none',
          text = element_text(size=20)
    ) +
    xlim(0, 0.5)+
    scale_color_manual(values=pnw_palette("Bay",5)[c(1,3,5)],
                       name="Location",
                       labels=c("Georgia","Iowa","BC")) +
    facet_wrap(~TRAIT)
  
  
  ggsave(paste0(SAVE_DIR,"/Secondary_Germplasm_SLOPE_POINT_PLOT.PDF"), 
         width = 10,
         height = 3)
  
  