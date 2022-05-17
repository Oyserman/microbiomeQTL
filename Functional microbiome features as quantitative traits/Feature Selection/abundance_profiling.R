library(ggplot2)
library(dplyr)
library(reshape2)
library(gplots)

norm_log_counts <- read.delim("final_assembly/feature_selection/results/norm_log_common_contig_conunts.tsv")

rhizosphere_contig <- read.delim("final_assembly/feature_selection/results/union_Rhizosphere_contig.txt", header = F)
rhizosphere_contig <- rhizosphere_contig$V1

rhizosphere_count <- norm_log_counts[norm_log_counts$Taxa.and.Samples %in% rhizosphere_contig,]

taxon <- read.delim("final_assembly/kraken/2020_april/contig_10k_taxonomy.tsv", header = F)
rhizosphere_taxon <- taxon[taxon$V1 %in% rhizosphere_contig,]

#df_rhizosphere_count <- melt(rhizosphere_count, value.name = "count", variable.name = "sample")
#df_rhizosphere_count <- t(df_rhizosphere_count)

rownames(rhizosphere_count) <- rhizosphere_count$Taxa.and.Samples
rhizosphere_count <- rhizosphere_count[,-1]
mat <- as.matrix(sapply(rhizosphere_count, as.numeric))
rownames(mat) <- rownames(rhizosphere_count)

hist(rhizosphere_count$X249_P_L001)

hist(rhizosphere_count$BULK1.2_L001)

pdf(file = "heatmap_count_with_colorkey.pdf", width = 12, height = 12)
#heatmap(mat, scale = "none",cexRow = 0.2, cexCol = 0.5)
heatmap.2(mat,scale = "none", cexRow = 0.2, cexCol = 0.5,
          trace = "none", density.info = "none")
dev.off()
