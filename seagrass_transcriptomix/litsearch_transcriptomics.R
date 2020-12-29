library(dplyr)
library(ggplot2)
library(ggraph)
library(igraph)
library(readr)
library(litsearchr)

naive_results <- import_results(file="DATA/scopus_seagrass.ris")

nrow(naive_results)

colnames(naive_results)

naive_results$year <- as.numeric(naive_results$year)

naive_results <- naive_results[order(naive_results$year),]


naive_results[,c("title", "year")]


keywords <- extract_terms(keywords=naive_results[, "keywords"], method="tagged", min_n=1)

title_terms <- extract_terms(text=naive_results[, "title"], method="fakerake", min_freq=3, min_n=2)

terms <- unique(c(keywords, title_terms))

docs <- paste(naive_results[, "title"], naive_results[, "abstract"])



