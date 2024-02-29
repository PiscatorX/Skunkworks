#1. Load packages
library(phyloseq)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(vegan)

#setwd("C:/Emma/Masters 2023/Statistics/Chapter 1/Dummy data")

#2. Import dummy data and edit accordingly
    #2.1 Import otu table
      asv_sheet <- read.csv("C:/Emma/Masters 2023/Statistics/Chapter 1/Accumulation curves/asv_sheet.csv")
      colnames(asv_sheet)[1] <- "Row_ID"
      asv_sheet <- column_to_rownames(asv_sheet, "Row_ID")

    #2.2 Import taxa table
      tax_sheet <- read.csv("C:/Emma/Masters 2023/Statistics/Chapter 1/Dummy data/tax_sheet.csv", row.names=1)
      tax_sheet <- subset(tax_sheet, select = "superkingdom", "phylum", "class", "oder","family", "genus", "species")

    colnames(tax_sheet) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")

    #2.3 Import sample metadata table
      sample_sheet <- read.csv("C:/Emma/Masters 2023/Statistics/Chapter 1/Dummy data/sample_sheet.csv")
      sample_sheet <- subset(sample_sheet, select = "Site", "Method")

class(tax_sheet)
class(asv_sheet)

#3. Make phyloseq object
OTU = otu_table(asv_sheet, taxa_are_rows = TRUE)
TAX = tax_table(tax_sheet)

view(OTU)
view(TAX)

physeq <- phyloseq(otu_table(asv_sheet, taxa_are_rows = TRUE),
                   tax_table(tax_sheet))

physeq

#4. Make some graphs 
#4.1 Plot bar
plot_bar(physeq, fill = "Class")

#4.2 Accumulation curves (rarefaction curve)

# Rarefaction analysis using vegan
rarefaction_results <- rarefy(physeq, sample = 100)

# Plot rarefaction curves
plot(rarefaction_results, col = "blue", lty = 1, xlab = "Sequencing Depth", ylab = "OTUs")

#4.3 Accumulation curves (rarefaction curve)
acc_curve <- specaccum(physeq)
plot(acc_curve, xlab = "Number of Samples", ylab = "OTUs")