# this script takes the name of a study and computes runs through a "standard" DESeq2 analysis workflow

library(tximport)
library(DESeq2)
library(RColorBrewer)
library(dplyr)
library(tidyr)
library(ggplot2)
library(RColorBrewer)
library(pheatmap)
library(GGally)
library(gridExtra)
library(gdata)

study <- c("MT_CAM_NOP") # study input name, change this depending on what you are interested in
  
# read in meta-data
meta <- read.table(paste("../Metadata/",study,".txt", sep=""), header=T, sep = "\t")
# format/add additional columns 
rownames(meta) <- meta$run_accession
meta <- separate(meta, condition, c("genotype", "strain", "tissue", "age"), sep = "_", remove = F)
meta$genotype.tissue <- interaction(meta$strain, meta$tissue, sep=".") 
meta$condition.study <- interaction(meta$condition, meta$study_accession,  sep=".")
meta$condition.rep <- interaction(meta$condition, meta$rep,  sep=".")

# create file paths
files <- file.path("..", "tak3.1_kallisto", study, 
                     meta$run_accession, 
                     "abundance.tsv") # files to be pooled this can be modified for individual studies

# building transcript-gene list from one example file
tx2gene <- read.table("../tak3.1_kallisto/PRJDB4420/DRR050343/DRR050343.tsv", header=T)
tx2gene <- tibble(tx2gene$target_id, gsub(".1", "", tx2gene$target_id, fixed=T))

colnames(tx2gene) <- c("TXNAME", "GENEID")

names(files) <- meta$run_accession

txi <- tximport(files, type = "kallisto", tx2gene = tx2gene) # pooling all samples together

dir.create(study)

write.csv(txi$abundance, paste(study,"/",study,"_abundance.csv", sep=""))
write.csv(txi$abundance, paste(study,"/",study,"_counts.csv", sep=""))

# import into DESeq2

de <- DESeqDataSetFromTximport(txi, meta, ~ condition)

length(counts(de))
keep <- rowSums(counts(de)) >= 10 # filtering out genes with < 5 counts in total
de <- de[keep,] 
length(counts(de))
keep2 <- apply(counts(de), 1, function(row) all(row !=0 )) # filtering any gene that has no expression in any condition
de <- de[keep2,] 
length(counts(de))

vsd <- vst(de)
rld <- rlog(de)

# overview plots, PCA and heatmaps

# PCA plot
pdf(paste(study,"/",study,"_PCA.pdf", sep=""))
plotPCA(vsd, intgroup=c("condition"))
dev.off()

# heatmap of sample distances
sample_dist <- dist(t(assay(rld)))
dist_matrix <- as.matrix(sample_dist)
rownames(dist_matrix) <- vsd$condition.rep
colnames(dist_matrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

pdf(paste(study,"/",study,"_dist.pdf", sep=""))
pheatmap(dist_matrix,
         clustering_distance_rows=sample_dist,
         clustering_distance_cols=sample_dist,
         col=colors)
dev.off()

# correlation matrix

pdf(paste(study,"/",study,"_corr.pdf", sep=""), width=1.5*length(vsd$run_accession), height = 1.5*length(vsd$run_accession))
ggpairs(as.data.frame(assay(vsd))) + theme_light()
dev.off()

# differential gene expression

de <- DESeq(de) 


# gene dispersion
pdf(paste(study,"/",study,"_disp.pdf", sep=""))
plotDispEsts(de)
dev.off()

# heatmap of variable genes 

topVarGenes <- head( order( rowVars( assay(rld) ), decreasing=TRUE ), 50 )

pdf(paste(study,"/",study,"_hm_top50_var.pdf", sep=""))

annot <-  as.data.frame(colData(de)[,c("condition", "rep")])
annot$rep <- as.factor(annot$rep)
pheatmap(assay(vsd)[ topVarGenes, ], cluster_rows=T, show_rownames=T, annotation_col = annot)
dev.off()


# Differential gene results/plotting

conditions_u <- unique(de$condition)
conditions_u <- as.character(conditions_u)

conditions_name <- combn(as.character(unique(de$condition)),2, FUN=paste, collapse='.')

conditions_int <- combn(as.character(unique(de$condition)),2)

p <- list()

for (x in 1:length(conditions_int[1,])) {

  res <- paste("res",x,sep="")
    
  assign(res, results(de, contrast = c("condition",conditions_int[1,x],conditions_int[2,x])))
  
  #plotMA(get(paste("res",x, sep = "")), alpha=0.01, main=conditions_name[x])
  
  # MAplot
  p[[x]] <- ggplot(as.data.frame(get(paste("res",x, sep = ""))), 
         aes( baseMean,log2FoldChange, color=-log10(padj)>2)) + 
    geom_point(size=0.5) + scale_x_log10() + scale_color_manual(values=c("grey50", "red3")) + 
    theme_light() + ggtitle(conditions_name[x])
  
  # Vulcan plot
  p[[x+length(conditions_name)]] <- ggplot(as.data.frame(get(paste("res",x, sep = ""))), 
                                           aes(log2FoldChange, -log10(padj), color=-log10(padj)>2)) + 
    geom_point(size=0.5) + scale_color_manual(values=c("grey50", "red3")) +
    annotate("text", x=0,y=-10, label = paste("DE genes at FDR = 0.01: ",length(which(get(paste("res",x, sep = ""))$padj<0.01)))) +
    theme_light() + ggtitle(conditions_name[x]) 
  
  # export each result file as seperate csv
  
  write.csv(get(paste("res",x, sep = "")), paste(study,"/",conditions_name[x],".csv",sep=""))
  
  # plot heatmap of top50 sig genes for each comparison
  
  top50 <- head(arrange(tibble::rownames_to_column(as.data.frame(get(paste("res",x, sep = ""))),var="Phytozome.ID"), padj),50) 
  top50_genes <- top50$Phytozome.ID
  
  pdf(paste(study,"/",conditions_name[x],"_top50_sig.pdf", sep="")) 
  pheatmap(assay(vsd)[top50_genes,], cluster_rows=T, show_rownames=T,
           cluster_cols=T,annotation_col = annot)
  dev.off()
  
  }

# generate a plotting matrix for MA and Vulcan plots
m <- matrix(0,length(unique(de$condition)),length(unique(de$condition)))

upperTriangle(m,diag = F,byrow=T) <- c(1:length(conditions_name))
  
lowerTriangle(m,diag = F,byrow=T) <- c((length(conditions_name) + 1:length(conditions_name)))

diag(m) <- NA

# plotting command for forementioned matrix of plots
pdf(paste(study,"/",study,"_MA_VA.pdf", sep=""), 
    width=8*length(unique(de$condition)), height = 5*length(unique(de$condition)))
grid.arrange(grobs= p, layout_matrix = m)
dev.off()

dev.off()

