library(data.table)
library(dplyr)


args = commandArgs(trailingOnly = TRUE)
INTFRQ <- args[1] #introg.frequency
RESULT <- args[2] #GWA,GP, or other results
SAVE <- args[3] #directory and name to save the merged result

fread(INTFRQ,
      header=T) -> INTFRQ_FILE

  fread(RESULT,
        header=T) -> RESULT_FILE
  
  full_join(INTFRQ_FILE,RESULT_FILE) %>%
    
    fwrite(args[3],
           col.names = T,
           quote = F,
           sep = "\t",
           na = "NA")
    