# pre-processing data for shiny app

library(tximport)
library(DESeq2)
library(RColorBrewer)
library(dplyr)
library(tidyr)
library(ggplot2)
library(RColorBrewer)
library(pheatmap)

# set wd to /Analysis

meta_all <- read.csv("meta_all.csv", row.names = 1) # make sure this is up to date

# reading in kallisto output files 

files <- file.path("..", "tak3.1_kallisto", meta_all$study_accession, 
                   meta_all$run_accession, 
                   paste(meta_all$run_accession,".tsv", sep="")) # files to be pooled this can be modified for individual studies

tx2gene <- read.table("../tak3.1_kallisto/PRJDB4420/DRR050343/DRR050343.tsv", header=T)
tx2gene <- tibble(tx2gene$target_id, gsub(".1", "", tx2gene$target_id, fixed=T))

colnames(tx2gene) <- c("TXNAME", "GENEID")

names(files) <- meta_all$run_accession

txi.all <- tximport(files, type = "kallisto", tx2gene = tx2gene) # pooling all samples together

saveRDS(txi.all, file="txi.all.RData")
txi_all <- readRDS("txi.all.RData")

# next we import the txi object into DESeq2
de_all <- DESeqDataSetFromTximport(txi.all, meta_all, ~ condition)

saveRDS(de_all, file="de.all.RData")
de_all <- readRDS("de.all.RData")

vsd <- vst(de_all) # variance stabalised gene count

saveRDS(vsd, file="vst.all.RData")
vsd <- readRDS("vst.all.RData")


