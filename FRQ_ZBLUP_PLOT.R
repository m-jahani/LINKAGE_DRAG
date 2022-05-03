library(data.table)
library(tidyverse)
library(cowplot)
library(PNWColors)

#calculates Zscore of blups and draw plots
args = commandArgs(trailingOnly = TRUE)
INTRO_GRQ <- args[1] #/Users/mojtabajahani/Documents/Projects/Linkage_drag/New_method/GP/FRQ_VS_BLUP/2nd_GERMPLASM_ALL_INTROGRESSION_EFFECTS_FRQUENCY"
DONOR <- args[2] #"Secondary Germplasm"
SAVE_DIR <- args[3] 

# INTRO_GRQ <- "/Users/mojtabajahani/Documents/Projects/Linkage_drag/New_method/GP/FRQ_VS_BLUP/ANNUUS_ALL_INTROGRESSION_EFFECTS_FRQUENCY"
# DONOR <- "WILD ANNUUS"
# SAVE_DIR <- "/Users/mojtabajahani/Documents/Projects/Linkage_drag/New_method/GP/FRQ_VS_BLUP"


fread(INTRO_GRQ) %>%
  mutate(trait_loc = paste0(TRAIT," (",LOCATION,")")) -> DATA

DATA %>%
  distinct(trait_loc) %>%
  pull(trait_loc) -> TRAIT_LOCA

SCALED_BLUP <- NULL
for (i in 1:length(TRAIT_LOCA)) {
  DATA %>%
    filter(trait_loc==TRAIT_LOCA[i]) %>%
    mutate_at(vars(BLUP), list(ZBLUP = ~as.vector(scale(.,center = T,scale = T)))) %>% 
    rbind(.,SCALED_BLUP) -> SCALED_BLUP
}
  


SCALED_BLUP %>%
  ggplot(., aes(x=FRQ, y=ZBLUP)) +
  geom_smooth(method = "lm") +
  xlab("Frequancy of introgression") +
  ylab("Z score Effect size of introgression") +
  geom_hline(yintercept=0,linetype="dotted") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        plot.title = element_text(hjust = 0.5),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.position='none'
  )+
  facet_wrap(~trait_loc)+
  ggtitle(DONOR)

ggsave(paste0(SAVE_DIR,"/",DONOR,"_SLOPE_LINE_PLOT.PDF"))


RESULT <- matrix(nrow = length(TRAIT_LOCA),
                 ncol = 6)
for (trait in 1:length(TRAIT_LOCA)) {
  
  FRQU <- pull(filter(SCALED_BLUP,trait_loc == TRAIT_LOCA[trait]),FRQ)
  ZBLUPs <- pull(filter(SCALED_BLUP,trait_loc == TRAIT_LOCA[trait]),ZBLUP)
  fit=lm(ZBLUPs ~ FRQU)
  RESULT[trait,1] <-  TRAIT_LOCA[trait] #trait location ID
  RESULT[trait,2] <-  pull(distinct(filter(SCALED_BLUP,trait_loc == TRAIT_LOCA[trait]),TRAIT),TRAIT) #trait id
  RESULT[trait,3] <-  pull(distinct(filter(SCALED_BLUP,trait_loc == TRAIT_LOCA[trait]),LOCATION),LOCATION) #location id
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
  select(TRAIT_LOC = V1,
         TRAIT = V2,
         LOCATION = V3,
         R_SQUARED,
         SLOPE,
         P_SLOPE,
         TEST) -> RESULT

RESULT %>%
  ggplot(.,aes(x=TRAIT,y=SLOPE,color=LOCATION,shape=TEST)) +
  geom_point(size=3) +
  geom_hline(yintercept=0,linetype="dotted") +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  xlab("Trait") +
  ylab("Standard Beta Coefficient (slope)") +
  ggtitle(DONOR)+
  scale_color_manual(values=pnw_palette("Bay",5)[c(1,3,5)],
                     name="Location",
                     labels=c("Georgia","Iowa","BC")) +
  scale_shape_manual(values=c(8, 19))+
  coord_cartesian(ylim=c(
    (max(abs(RESULT$SLOPE))+0.01)*-1,
    max(abs(RESULT$SLOPE))+0.01)) 
    
ggsave(paste0(SAVE_DIR,"/",DONOR,"_SLOPE_POINT_PLOT.PDF"))
    
