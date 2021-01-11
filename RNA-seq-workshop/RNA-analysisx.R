library(GenomicAlignments)
library(BiocParallel)
library(GenomicFeatures)
library(Rsamtools)
library(airway)
library(dplyr)


register(SerialParam())

#This is toy data for test practise
airway_dir  <- system.file("extdata", package = "airway", mustWork = T)

list.files(airway_dir)

csvfile <-  file.path(airway_dir, "sample_table.csv")

(sample_table <-  read.csv(csvfile, row.names = 1))

nrow(sample_table)

sample_table %>% group_by(dex) %>% summarise(Count= n())

(filenames <- file.path(airway_dir, paste0(sample_table$Run,"_subset.bam")))

file.exists(filenames)

(bam_files <- BamFileList(filenames,  yieldSize = 2000000))

seqinfo(bam_files[1])

gtf <- file.path(airway_dir,  "Homo_sapiens.GRCh37.75_subset.gtf")
(gtf_txdb <- makeTxDbFromGFF(gtf, format = "gtf",  circ_seqs = character()))

(ebg <- exonsBy(gtf_txdb, by="gene"))


se <- summarizeOverlaps(features = ebg,
                        reads = bam_files,
                        mode = "Union",
                        singleEnd = F,
                        ignore.strand = T,
                        fragments = T)
                        


se

dim(se)

assay(se)

assayNames(se)

colSums(assay(se))

rowRanges(se)

str(metadata(rowRanges(se)))

colData(se)

(colData(se) <- DataFrame(sample_table))

se$dex


se$cell

#first level of a factor be the reference level 

se$dex <- relevel(se$dex, "untrt")



