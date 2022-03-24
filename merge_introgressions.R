library(data.table)
library(tidyverse)
library(IRanges)

fread("/Users/mojtabajahani/Documents/Projects/Linkage_drag/pcadmix.xxx.regions.txt") %>% 
  separate(hap, into=c("sample","haplotype"),sep = "_") %>% 
  mutate(LOOP = paste(CHR,sample,anc,sep = "_")) -> DATA

#loop vector over chromosome sample and ancestry
DATA %>% distinct(LOOP) %>% pull(LOOP) -> CHR_sample_anc

#merge introgression of haplotypes
hap_merged_data <- NULL
for (i in 1:length(CHR_sample_anc)) {
  DATA %>% 
    filter(LOOP == CHR_sample_anc[i]) %>%
    select(start,end)  -> coordinates
  IRanges::reduce(
    IRanges(pull(coordinates,start),pull(coordinates,end)) 
  ) %>%
    as.data.frame() %>%
    mutate(CHR_sample_anc = CHR_sample_anc[i]) %>%
    separate(CHR_sample_anc , 
             into =c("CHR","sample","anc"), 
             sep ="_" , 
             remove = T) %>%
    select(CHR,
           start,
           end,
           sample,
           anc) %>%
    rbind(.,hap_merged_data) -> hap_merged_data
  
}

rm(DATA,
   CHR_sample_anc,
   coordinates,
   i)


hap_merged_data %>%
  mutate(LOOP = paste(CHR,anc,sep = "_")) -> DATA2

#loop vector over chromosome and ancestry
DATA2 %>% distinct(LOOP) %>% pull(LOOP) -> CHR_anc

#merge introgression of samples to common regions for all samples
sample_merged_data <- NULL
for (i in 1:length(CHR_anc)) {
  DATA2 %>% 
    filter(LOOP == CHR_anc[i]) %>%
    select(start,end)  -> coordinates
  
  IRanges::reduce(
    IRanges(pull(coordinates,start),pull(coordinates,end)) 
  ) %>%
    as.data.frame() %>%
    mutate(CHR_anc = CHR_anc[i]) %>%
    separate(CHR_anc , 
             into =c("CHR","anc"), 
             sep ="_" , 
             remove = T) %>%
    select(CHR,
           start,
           end,
           anc) -> merged_data
  
  findOverlaps(
    IRanges(pull(coordinates,start), pull(coordinates,end)) ,
    IRanges(pull(merged_data,start), pull(merged_data,end))
  ) %>% as.data.frame() -> overlap
  
  DATA2 %>% 
    filter(LOOP == CHR_anc[i])  %>%
    mutate(queryHits = row_number()) %>% 
    left_join(.,overlap) %>% 
    left_join(.,
              select(
                mutate(merged_data,subjectHits = row_number()),
                merged_start = start,
                merged_end = end,
                subjectHits
              )
    ) %>%
    select(sample,
           CHR,
           start,
           end,
           anc,
           merged_start,
           merged_end) %>%
    rbind(.,sample_merged_data) -> sample_merged_data
}

rm(CHR_anc,
   coordinates,
   DATA2,
   hap_merged_data,
   merged_data,
   overlap,
   i)
#find introgression that are much larger than those in same group
sample_merged_data %>%
  mutate(merge_length=(merged_end-merged_start)) %>%
  mutate(intg_length=(end-start)) %>%
  left_join(.,
            summarise(
              group_by(
                mutate(
                  sample_merged_data,intg_length=(end-start)),
                merged_start,merged_end),
              mean_intg_length=mean(intg_length))) %>%
  mutate(merge_intgmean_diff=(intg_length-mean_intg_length)) %>%
  mutate(prop_residual=merge_intgmean_diff/mean_intg_length) %>%
  mutate(proper_merge=ifelse(prop_residual > 2,"LARGE","GOOD")) %>% #twice bigger than the mean length in group is considered as large
  arrange(CHR,anc,merged_start) -> DATA_TMP1

rm(sample_merged_data)
#vector of groups with LARGE introgression
DATA_TMP1 %>% filter(proper_merge=="LARGE") %>%
  mutate(FILT = paste(anc,CHR,merged_start,merged_end,sep = "_")) %>%
  distinct(FILT) %>%
  pull(FILT) -> TO_FILER
#groups without LARGE introgression (the first part of merged introgression)
DATA_TMP1 %>%
  mutate(FILT=paste(anc,CHR,merged_start,merged_end,sep = "_")) %>%
  filter(!FILT %in% TO_FILER) %>%
  select(-FILT) -> DONE_part1

#groups with LARGE introgression 
DATA_TMP1 %>%
  mutate(FILT = paste(anc,CHR,merged_start,merged_end,sep = "_")) %>%
  filter(FILT %in% TO_FILER) %>%
  select(-FILT) -> GOOD_BAD

#good size introgression in the groups with LARGE introgression
GOOD_BAD %>%
  filter(proper_merge=="GOOD") %>%
  select(CHR,start,end,sample,anc) %>%
  mutate(LOOP = paste(CHR,anc,sep = "_")) -> DATA3

DATA3 %>% distinct(LOOP) %>% pull(LOOP) -> CHR_anc

rm(DATA_TMP1,
   TO_FILER)

DONE_part2 <- NULL
for (i in 1:length(CHR_anc)) {
  DATA3 %>% 
    filter(LOOP == CHR_anc[i]) %>%
    select(start,end)  -> coordinates
  
  IRanges::reduce(
    IRanges(pull(coordinates,start),pull(coordinates,end)) 
  ) %>%
    as.data.frame() %>%
    mutate(CHR_anc = CHR_anc[i]) %>%
    separate(CHR_anc , 
             into =c("CHR","anc"), 
             sep ="_" , 
             remove = T) %>%
    select(CHR,
           start,
           end,
           anc) -> merged_data
  
  findOverlaps(
    IRanges(pull(coordinates,start), pull(coordinates,end)) ,
    IRanges(pull(merged_data,start), pull(merged_data,end))
  ) %>% as.data.frame() -> overlap
  
  DATA3 %>% 
    filter(LOOP == CHR_anc[i])  %>%
    mutate(queryHits = row_number()) %>% 
    left_join(.,overlap) %>% 
    left_join(.,
              select(
                mutate(merged_data,subjectHits = row_number()),
                merged_start = start,
                merged_end = end,
                subjectHits
              )
    ) %>%
    select(sample,
           CHR,
           start,
           end,
           anc,
           merged_start,
           merged_end) %>%
    rbind(.,DONE_part2) -> DONE_part2
}

rm(CHR_anc,
   coordinates,
   DATA3,
   merged_data,
   overlap,
   i)

#extract the merged introgression of the groups with large introgression (after removing large introgression) - (#merged introgression of goods after new merging without LARGE)
DONE_part2 %>% 
  mutate(FILT2 = paste(anc,CHR,merged_start,merged_end,sep = "_")) %>% 
  distinct(FILT2,.keep_all = T) %>% 
  select(-FILT2) %>%
  mutate(LOOP = paste(CHR,anc,sep = "_")) -> merged_good_intg

#large introgression
GOOD_BAD %>%
  filter(proper_merge == "LARGE") %>% 
  select(CHR,start,end,sample,anc) %>%
  mutate(LOOP = paste(CHR,anc,sep = "_")) -> DATA4


rm(GOOD_BAD)

DATA4 %>% distinct(LOOP) %>% pull(LOOP) -> CHR_anc
#overlap between merged introgression and large introgression
DONE_part3 <- NULL
for (i in 1:length(CHR_anc)) {
  DATA4 %>% 
    filter(LOOP == CHR_anc[i]) %>% 
    select(start,end)  -> coordinates
  
  merged_good_intg %>%
    filter(LOOP == CHR_anc[i]) %>%
    select(start=merged_start,end=merged_end) -> merged_good_intg_vec
  
  findOverlaps(
    IRanges(pull(coordinates,start), pull(coordinates,end)) ,
    IRanges(pull(merged_good_intg_vec,start), pull(merged_good_intg_vec,end))
  ) %>% as.data.frame() -> overlap
  
  
  DATA4 %>% 
    filter(LOOP == CHR_anc[i]) %>% 
    mutate(queryHits = row_number()) %>% 
    select(sample,
           CHR,
           start,
           end,
           anc,
           queryHits) %>% 
    left_join(.,overlap) %>% 
    left_join(.,
              select(
                mutate(merged_good_intg_vec,subjectHits = row_number()),
                merged_start = start,
                merged_end = end,
                subjectHits
              )
    ) %>%
    select(sample,
           CHR,
           start,
           end,
           anc,
           merged_start,
           merged_end) %>%
    mutate(merged_start=ifelse(is.na(merged_start),start,merged_start)) %>%
    mutate(merged_end=ifelse(is.na(merged_end),end,merged_end)) %>%
    filter(!is.na(merged_start)) -> DATA5 
  
  #break down the large introgression to two part: 1. same size of the rest in the group 2. the rest of introgression length
  DATA5 %>% 
    select(sample,
           CHR,
           anc,
           new_start = merged_start,
           new_end = merged_end) -> DATA6
  
  
  DATA5 %>% mutate(loop=row_number()) -> DATA5
  DATA7 <- NULL
  for (i in 1:nrow(DATA5)) {
    
    DATA5 %>%
      filter(loop == i) %>%
      mutate(new_start =
               ifelse(start >= merged_start && (end > merged_end && start <= merged_end ) ,#first condition
                      merged_end,#if first condition correct
                      ifelse((start < merged_start && merged_start < end) && (end <= merged_end ),#second condition
                             start,#if second condition correct
                             ifelse(start >= merged_start && end <= merged_end,#third condition
                                    "REMOVE",#if third condition correct
                                    ifelse(end <= merged_start && end < merged_end,#fourth condition
                                           "REMOVE",
                                           ifelse((start < merged_start && end > merged_start) && (start < merged_end && end > merged_end) ,#fifth condition
                                                  start,
                                                  "FAIL"))#if fifth condition is wrong some thing has a problem
                             )))) %>%
      mutate(new_end =
               ifelse(
                 start >= merged_start && (end > merged_end && start <= merged_end ) ,#first condition
                 end,#if first condition correct
                 ifelse((start < merged_start && merged_start < end) && (end <= merged_end ),#second condition
                        merged_start,#if second condition correct
                        ifelse(start >= merged_start && end <= merged_end,#third condition
                               "REMOVE",#if third condition correct
                               ifelse(end <= merged_start && end < merged_end,#fourth condition
                                      "REMOVE",
                                      ifelse((start < merged_start && end > merged_start) && (start < merged_end && end > merged_end) ,#fifth condition
                                             merged_start,
                                             "FAIL"))#if fifth condition is wrong some thing has a problem
                        ))))  %>%
      filter(new_start !="REMOVE") %>%
      rbind(.,DATA7) -> DATA7
    
    DATA5 %>% 
      filter(loop == i) %>%
      mutate(new_start =
               ifelse(
                 start >= merged_start && (end > merged_end && start <= merged_end ) ,#first condition
                 merged_end,#if first condition correct
                 ifelse((start < merged_start && merged_start < end) && (end <= merged_end ),#second condition
                        start,#if second condition correct
                        ifelse(start >= merged_start && end <= merged_end,#third condition
                               "REMOVE",#if third condition correct
                               ifelse(end <= merged_start && end < merged_end,#fourth condition
                                      "REMOVE",
                                      ifelse((start < merged_start && end > merged_start) && (start < merged_end && end > merged_end) ,#fifth condition
                                             merged_end,
                                             "FAIL"))#if fifth condition is wrong some thing has a problem
                        )))) %>%
      mutate(new_end =
               ifelse(
                 start >= merged_start && (end > merged_end && start <= merged_end ) ,#first condition
                 end,#if first condition correct
                 ifelse((start < merged_start && merged_start < end) && (end <= merged_end ),#second condition
                        merged_start,#if second condition correct
                        ifelse(start >= merged_start && end <= merged_end,#third condition
                               "REMOVE",#if third condition correct
                               ifelse(end <= merged_start && end < merged_end,#fourth condition
                                      "REMOVE",
                                      ifelse((start < merged_start && end > merged_start) && (start < merged_end && end > merged_end) ,#fifth condition
                                             end,
                                             "FAIL"))#if fifth condition is wrong some thing has a problem
                        ))))  %>%
      filter(new_start !="REMOVE") %>%
      rbind(.,DATA7) -> DATA7
  }
  
  DATA7 %>%
    distinct(sample,
             CHR,
             start,
             end,
             anc,
             merged_start,
             merged_end,
             new_start,
             new_end) %>% 
    select(sample,
           CHR,
           anc,
           new_start,
           new_end) %>% 
    rbind(DATA6) %>% 
    rbind(.,
          DONE_part3) -> DONE_part3
  
  DONE_part3 %>% filter(is.na(new_start))
  DONE_part3 %>% filter(new_start=="FAIL")
  
  
}


rm(CHR_anc,
   coordinates,
   DATA4,
   DATA5,
   DATA6,
   DATA7,
   merged_good_intg,
   merged_good_intg_vec,
   overlap,
   i)


rbind(
  rbind(
    select(DONE_part1,
           sample,
           CHR,
           anc,
           new_start = merged_start,
           new_end = merged_end),
    select(DONE_part2,
           sample,
           CHR,
           anc,
           new_start = merged_start,
           new_end = merged_end)),
  DONE_part3) -> INTREGRESSION_SAM 

rm(DONE_part1,
   DONE_part2,
   DONE_part3)

data.frame(anc=c("1","2"),
           donor=c("ANN","OTHER"))-> anc_donor

  INTREGRESSION_SAM %>% 
    group_by(anc,CHR,new_start,new_end) %>%
    tally() %>%
    ungroup() %>%
    filter(n > 15) %>% #5 precent frequency
    select(-n) %>% 
    left_join(.,anc_donor)  %>% 
    mutate(ID = paste0(CHR,":",new_start,"-",new_end,"_",donor)) %>%
    select(CHR,
           start = new_start,
           end = new_end,
           ID)  -> INTREGRESSION_SAM_maf0.5
  
  fwrite(INTREGRESSION_SAM_maf0.5,"/Users/mojtabajahani/Documents/Projects/Linkage_drag/INTREGRESSION_SAM_maf0.5.csv",
         sep = ",",
         quote = F,
         col.names = T)
  
  rm(anc_donor,
     INTREGRESSION_SAM)