library(tidyverse)
library(plotly)
library(caret)
library(stats)
library(dplyr)



zeo_ref <- read.table("zeolites_DB.tsv", header = T, sep = "\t", stringsAsFactors = F ) 

numeric_cols  <- zeo_ref %>%
                 select_if(is.numeric) %>%
                 colnames()

for (col in numeric_cols){
  zeo_ref[col] <- 0
  
}

zeo_ref_uniq <- do.call(cbind,   lapply(zeo_ref,   unique)) %>% 
                data.frame()

write.table(zeo_ref_uniq, "zeolite_template.tsv", sep = "\t", quote = F, row.names = FALSE)


