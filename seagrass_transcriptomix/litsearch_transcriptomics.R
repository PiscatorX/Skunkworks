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

dfm <- create_dfm(elements=docs, features=terms)

g <- create_network(dfm, min_studies=3)

ggraph(g, layout="stress") +
      coord_fixed() +
      expand_limits(x=c(-3, 3)) +
      geom_edge_link(aes(alpha=weight)) +
      geom_node_point(shape="circle filled", fill="white") +
      geom_node_text(aes(label=name), hjust="outward", check_overlap=TRUE) +
      guides(edge_alpha=FALSE)


strengths <- strength(g)

term_strengths <- data.frame(term=names(strengths), strength=strengths, row.names=NULL) %>%
                   mutate(rank=rank(strength, ties.method="min")) %>%
                   arrange(strength) 
             

term_strengths

cutoff_fig <- ggplot(term_strengths, aes(x=rank, y=strength, label=term)) +
              geom_line() +
              geom_point() +
              geom_text(data=filter(term_strengths, rank>5), hjust="right", nudge_y=20, check_overlap=TRUE)

cutoff_fig

cutoff_cum <- find_cutoff(g, method="cumulative", percent=0.8)

cutoff_cum

cutoff_fig +
     geom_hline(yintercept=cutoff_cum, linetype="dashed")


term_strengths %>%
               arrange(desc(strength)) %>%
               select(-rank) %>%
               filter(strength > cutoff_cum)
