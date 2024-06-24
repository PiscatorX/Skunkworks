#/usr/bin/env Rscript
library(stringr)
library(dada2)

fnFs <- sort(Sys.glob("eDNA-0046_seqdata_demultiplexed/seqdata_demultiplexed_12S*/*_R1.fastq.gz"))
fnRs <- sort(Sys.glob("eDNA-0046_seqdata_demultiplexed/seqdata_demultiplexed_12S*/*_R2.fastq.gz"))

lapply(list(fnFs, fnRs), length)

#This removes the unknown files
fnFs <- fnFs[stringr::str_detect(fnFs, "unknown_R1.fastq.gz", negate = TRUE)]
fnRs <- fnRs[stringr::str_detect(fnRs, "unknown_R2.fastq.gz", negate = TRUE)]

lapply(list(fnFs, fnRs), length)

#remove non-zero size files
non_zero_fnFs <- sapply(fnFs, function(file) file.info(file)$size > 1)
fnFs <- fnFs[non_zero_fnFs]
non_zero_fnRs <- sapply(fnRs, function(file) file.info(file)$size > 1)
fnRs <- fnRs[non_zero_fnRs]

lapply(list(fnFs, fnRs), length)

# Extract sample names
# We use the "_R" to extract sample names
sample_names_F <- sapply(strsplit(basename(fnFs), "_R"), `[`, 1)
sample_names_R <- sapply(strsplit(basename(fnRs), "_R"), `[`, 1)
(sample_names <- intersect(sample_names_F,  sample_names_R))

length(sample_names)

fnFs <- fnFs[sample_names_F  %in%  sample_names]
fnRs <- fnRs[sample_names_R  %in%  sample_names]

lapply(list(fnFs, fnRs), length)

# Place filtered files in filtered/ subdirectory
path = "MiSeq_output/"
filtFs <- file.path(path, "filtered", paste0(sample_names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample_names, "_R_filt.fastq.gz"))
names(filtFs) <- sample_names
names(filtRs) <- sample_names

lapply(list(fnFs, filtFs, fnRs, filtRs), length)

data.frame(fnFs, filtFs, fnRs, filtRs)

#Relaxed filtering criteria for truncLen and maxEE
#This writes out the filenames created above to disk
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(0,0),
              maxN=0, maxEE=c(2,5), truncQ=2, rm.phix=TRUE,
              compress=TRUE, multithread=FALSE) # On Windows set multithread=FALSE

str(out)


filt_trim <- data.frame(out) |>
             dplyr::arrange(desc(reads.in)) |>
             tibble::rownames_to_column(var = "read_ID") |>
             tidyr::separate(read_ID, sep = "_R", into = c("sample_name","ext")) |>
             dplyr::select(-ext) |>
             dplyr::arrange(desc(reads.out)) 

filt_trim


lapply(list(fnFs, filtFs, fnRs, filtRs), length)

#some reads failed so we have to remove them from our list
filtFs <- filtFs[file.exists(filtFs)]
filtRs <- filtRs[file.exists(filtRs)]

#check our counts
#lots of reads lost
lapply(list(fnFs, filtFs, fnRs, filtRs), length)

errF <- learnErrors(filtFs, multithread = TRUE)
errR <- learnErrors(filtRs, multithread = TRUE)

dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)

#Merge reads
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)

#Generate sequence table
seqtab <- makeSequenceTable(mergers)

#Remove chimeras
seqtab_nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)

#check Reads failing
sum(seqtab_nochim)/sum(seqtab)


#Track reads through the pipeline
getN <- function(x) sum(getUniques(x))
dada_filtering_cols <- lapply(list(sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab_nochim)), data.frame) 
dada_filtering <- do.call(cbind, dada_filtering_cols) 
                
colnames(dada_filtering) <- c("denoisedF", "denoisedR", "merged", "nonchim")

(filtering_steps <- dada_filtering |>
                   data.frame() |>
                   tibble::rownames_to_column(var = "sample_name"))


(dada_filtering <- merge(filt_trim, filtering_steps, all.x = TRUE) |>
                   dplyr::arrange(desc(nonchim)))

write.table(dada_filtering, "dada_filtering.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

write.table(seqtab_nochim, "12S_otu_table.tsv", sep = "\t")

#Assign taxonomy
#https://www.reference-midori.info/download.php
#MIDORI2_UNIQ_NUC_GB259_srRNA_DADA2.fasta.gz
#only subsample of the database is used for speed
taxa <- assignTaxonomy(seqtab_nochim, "MIDORI2_subsample.fasta.gz", multithread=TRUE)

write.table(taxa, "12S_taxa_table.tsv", sep = "\t")
