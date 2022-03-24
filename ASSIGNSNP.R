#To assign snps to intervals
#runs as: Rscript GENEFINDER.R SNP_list file.bed output
#SNP_list: CHR"\t"POS #without header
#file.bed: CHR"\t"START_POS\"t"END_POS"\t"interval_ID #without header
#output: name and directory to save result in #CHR"\t”POS”\t”interval_ID format


library(foreach)
library(doParallel)
library(tidyverse)
library(data.table)

args = commandArgs(trailingOnly = TRUE)


fread(args[1],
      sep = "\t",
      header = F) %>%
  rename(CHR=V1,
         POS=V2) -> SNPS

fread(args[2],
      sep = "\t",
      header = F) -> BED

registerDoParallel(cores=48)
SNP_GENE <- foreach(i=1:nrow(BED), .combine='rbind', .errorhandling='stop') %dopar% {
  mutate(select(filter(SNPS,
                       CHR == as.character(BED[i,1]),
                       (POS >= as.numeric(BED[i,2]) &
                          POS <= as.numeric(BED[i,3]))
  ),
  CHR,POS),
  Gene_ID=as.character(BED[i,4]))}

SNPS %>%
  left_join(.,SNP_GENE) %>%
  fwrite(args[3],
         col.names = F,
         sep = "\t",
         na = "NA",
         quote = F)