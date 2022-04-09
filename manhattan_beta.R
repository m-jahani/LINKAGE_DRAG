library(dplyr)
library(ggplot2)
library(tidyr)
library(data.table)

args = commandArgs(trailingOnly = TRUE)

pheno1  <- args[1]
pheno2  <- args[2]
map  <- args[3]
save_dir <- args[4]


# pheno1  <- "/Users/mojtabajahani/Downloads/SAM_introgression_donor_ANNUUS_maf0.03_biomass__UBC.ps"
# pheno2  <- "/Users/mojtabajahani/Downloads/SAM_introgression_donor_2nd_GERMPLASM_maf0.03_biomass__UBC.ps"
# map  <- "/Users/mojtabajahani/Downloads/SAM_introgression_donor_ANNUUS.map"
# save_dir <- "/Users/mojtabajahani/Downloads"

ben_treshold_0.1_annuua <- 4.792041
ben_treshold_0.1_2ndgermplasm <- 2.880814

gsub(".ps","",gsub("/DATA/home/mjahani/LINKADE_DRAG/new_method/GWAS/result/SAM_introgression_donor_ANNUUS_maf0.03_","",pheno1)) -> TRAIT

pheno1_name <- "ANNUUS"
pheno2_name <- "2nd_GERMPLASM"

fread(as.character(pheno1),
      header = F) %>% 
  select(SNP = V1 ,
         Beta = V2,
         P = V4 ) %>% 
  separate(SNP,
           into = c("CHR","START_END") ,
           sep = ":" ,
           remove = F) %>%
    separate(START_END,
             into = c("START","END") ,
             sep = "-" ,
             remove = T) %>%
  mutate(BP = as.numeric(START)) %>%
  mutate(CHR = gsub("Ha412HOChr","CH",CHR)) %>%
  mutate(PL=ifelse(Beta>=0,(log10(P)*-1),log10(P))) %>%
  mutate(Donor = pheno1_name) %>% 
    select(Donor,
           CHR,
           BP,
           PL) %>%
      rbind(.,
            fread(as.character(pheno2),
                  header = F) %>% 
              select(SNP = V1 ,
                     Beta = V2,
                     P = V4 ) %>% 
              separate(SNP,
                       into = c("CHR","START_END") ,
                       sep = ":" ,
                       remove = F) %>%
              separate(START_END,
                       into = c("START","END") ,
                       sep = "-" ,
                       remove = T) %>%
              mutate(BP = as.numeric(START)) %>%
              mutate(CHR = gsub("Ha412HOChr","CH",CHR)) %>%
              mutate(PL=ifelse(Beta>=0,(log10(P)*-1),log10(P))) %>%
              mutate(Donor = pheno2_name) %>% 
              select(Donor,
                     CHR,
                     BP,
                     PL) )-> gwasResults
 
fread(as.character(map),
      header = F) %>% 
  select(SNP = V2) %>%
  separate(SNP,
           into = c("CHR","START_END") ,
           sep = ":" ,
           remove = F) %>%
  separate(START_END,
           into = c("BP","END") ,
           sep = "-" ,
           remove = T) %>%
  mutate(CHR = gsub("Ha412HOChr","CH",CHR)) %>%
  mutate(BP = as.numeric(BP)) %>%
  select(SNP,
         CHR,
         BP) %>%
  
  # Compute chromosome size
  group_by(CHR) %>%
  summarise(chr_len=max(BP)) %>%
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(gwasResults, ., by=c("CHR"="CHR")) %>%
  
  # Add a cumulative position of each SNP
  arrange(CHR, BP) %>%
  mutate( BPcum=BP+tot) -> don
rm(gwasResults)

axisdf = don %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

don %>% 
  group_by(CHR) %>% 
  summarize(start=min(BPcum),end=max(BPcum)) %>%
  ungroup() %>%
  mutate(row_num = row_number()) %>%
  filter(row_num %% 2 == 1) %>%
  select(start,end) -> CHR_Boundary

max(don$PL)-> MAXLP
min(don$PL)-> MINLP
max(don$BPcum)-> MAXPOS
min(don$BPcum)-> MINPOS
  
  
ggplot(don, aes(x=BPcum, y=PL, color=Donor)) +
  
  # Show all points
  geom_point(size=1) + #alpha=0.8, 
  scale_color_manual(values = c("#FF8000", "#3399FF")) + 
  #custom X axis:
  scale_x_continuous(label = axisdf$CHR, breaks= axisdf$center ) +
  scale_y_continuous(breaks = seq(floor(min(don$PL)), ceiling(max(don$PL)), by = 2)) +
  geom_hline(yintercept=as.numeric(ben_treshold_0.1_annuua), linetype="dashed",color = "red", size=0.2) +
  geom_hline(yintercept=as.numeric(ben_treshold_0.1_2ndgermplasm), linetype="dashed",color = "green", size=0.2) +
  geom_hline(yintercept=as.numeric(ben_treshold_0.1_annuua)*-1, linetype="dashed",color = "red", size=0.2) +
  geom_hline(yintercept=as.numeric(ben_treshold_0.1_2ndgermplasm)*-1, linetype="dashed",color = "green", size=0.2) +
  geom_hline(yintercept=0, linetype="solid",color = "black", size=0.5) +
  annotate("rect", xmin=as.numeric(CHR_Boundary[1,1]), xmax=as.numeric(CHR_Boundary[1,2]), ymin=-Inf, ymax=+Inf, alpha=0.09, fill="gray") +
  annotate("rect", xmin=as.numeric(CHR_Boundary[2,1]), xmax=as.numeric(CHR_Boundary[2,2]), ymin=-Inf, ymax=+Inf, alpha=0.09, fill="gray") +
  annotate("rect", xmin=as.numeric(CHR_Boundary[3,1]), xmax=as.numeric(CHR_Boundary[3,2]), ymin=-Inf, ymax=+Inf, alpha=0.09, fill="gray") +
  annotate("rect", xmin=as.numeric(CHR_Boundary[4,1]), xmax=as.numeric(CHR_Boundary[4,2]), ymin=-Inf, ymax=+Inf, alpha=0.09, fill="gray") +
  annotate("rect", xmin=as.numeric(CHR_Boundary[5,1]), xmax=as.numeric(CHR_Boundary[5,2]), ymin=-Inf, ymax=+Inf, alpha=0.09, fill="gray") +
  annotate("rect", xmin=as.numeric(CHR_Boundary[6,1]), xmax=as.numeric(CHR_Boundary[6,2]), ymin=-Inf, ymax=+Inf, alpha=0.09, fill="gray") +
  annotate("rect", xmin=as.numeric(CHR_Boundary[7,1]), xmax=as.numeric(CHR_Boundary[7,2]), ymin=-Inf, ymax=+Inf, alpha=0.09, fill="gray") +
  annotate("rect", xmin=as.numeric(CHR_Boundary[8,1]), xmax=as.numeric(CHR_Boundary[8,2]), ymin=-Inf, ymax=+Inf, alpha=0.09, fill="gray") +
  annotate("rect", xmin=as.numeric(CHR_Boundary[9,1]), xmax=as.numeric(CHR_Boundary[9,2]), ymin=-Inf, ymax=+Inf, alpha=0.09, fill="gray") +
  annotate("text", x = MAXPOS+300000000, y = MINLP-1.75, label = "ANNUUS",size=5) +
  annotate("text", x = MAXPOS+300000000, y = MINLP-2.50, label = "2nd_GERMPLASM",size=5) +

  # Custom the theme:
  theme_bw() +
  theme(
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    strip.background = element_blank(),
    text=element_text(size=18)
  )+
  ggtitle(as.character(TRAIT))+   
  labs(
    y = "log10 P-Value * Beta sign",
    x = "Position") -> P

if (nrow(filter(don,Donor == "ANNUUS", abs(PL) <= ben_treshold_0.1_annuua ))>0) {
  P + geom_segment(data = filter(don,Donor == "ANNUUS", abs(PL) <= ben_treshold_0.1_annuua ),aes(x = as.numeric(BPcum), y = min(don$PL)-1.5, xend = as.numeric(BPcum)+1000, yend = min(don$PL)-2),col="black") -> P }
if (nrow(filter(don,Donor == "ANNUUS", abs(PL) > ben_treshold_0.1_annuua , PL > 0))) {
P + geom_segment(data = filter(don,Donor == "ANNUUS", abs(PL) > ben_treshold_0.1_annuua , PL > 0),aes(x = as.numeric(BPcum), y = min(don$PL)-1.5, xend = as.numeric(BPcum)+1000, yend = min(don$PL)-2),col="green")  -> P}
  if (nrow(filter(don,Donor == "ANNUUS", abs(PL) > ben_treshold_0.1_annuua , PL < 0)) >0) {
  P + geom_segment(data = filter(don,Donor == "ANNUUS", abs(PL) > ben_treshold_0.1_annuua , PL < 0),aes(x = as.numeric(BPcum), y = min(don$PL)-1.5, xend = as.numeric(BPcum)+1000, yend = min(don$PL)-2),col="red")  -> P}
    if (nrow(filter(don,Donor == "2nd_GERMPLASM", abs(PL) <= ben_treshold_0.1_2ndgermplasm))) {
  P + geom_segment(data = filter(don,Donor == "2nd_GERMPLASM", abs(PL) <= ben_treshold_0.1_2ndgermplasm),aes(x = as.numeric(BPcum), y = min(don$PL)-2.25, xend = as.numeric(BPcum)+1000, yend = min(don$PL)-2.75),col="black")  -> P}
      if (nrow(filter(don,Donor == "2nd_GERMPLASM", abs(PL) > ben_treshold_0.1_2ndgermplasm , PL > 0))) {
  P + geom_segment(data = filter(don,Donor == "2nd_GERMPLASM", abs(PL) > ben_treshold_0.1_2ndgermplasm , PL > 0), aes(x = as.numeric(BPcum), y = min(don$PL)-2.25, xend = as.numeric(BPcum)+1000, yend = min(don$PL)-2.75),col="green")  -> P}
        if (nrow(filter(don,Donor == "2nd_GERMPLASM", abs(PL) > ben_treshold_0.1_2ndgermplasm , PL < 0))) {
  P + geom_segment(data = filter(don,Donor == "2nd_GERMPLASM", abs(PL) > ben_treshold_0.1_2ndgermplasm , PL < 0), aes(x = as.numeric(BPcum), y = min(don$PL)-2.25, xend = as.numeric(BPcum)+1000, yend = min(don$PL)-2.75),col="red")  -> P}
          if (nrow(filter(don,Donor == "2nd_GERMPLASM", abs(PL) > ben_treshold_0.1_2ndgermplasm , PL < 0))) {
  P + geom_segment(data = filter(don,Donor == "2nd_GERMPLASM", abs(PL) > ben_treshold_0.1_2ndgermplasm , PL < 0), aes(x = as.numeric(BPcum), y = min(don$PL)-2.25, xend = as.numeric(BPcum)+1000, yend = min(don$PL)-2.75),col="red")  -> P}

  ggsave(paste0(save_dir,"/",as.character(TRAIT),".pdf"),width=22,height=10)




