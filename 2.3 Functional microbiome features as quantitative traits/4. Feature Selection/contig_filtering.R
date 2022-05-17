#!/usr/bin/env Rscript

# program setup -----------------------------------------------------------
suppressPackageStartupMessages({
  library(tidyverse)
  library(ggpubr)
  library(RColorBrewer)
  library(cowplot)
  library(metagenomeSeq)
  library(ggplot2)
  library(dplyr)
  library(reshape2)
})

theme_set(theme_pubr(border = T, legend = "right"))
set.seed(56753)

# setwd("/mnt/mfs/bioinfdata/ngs2/ME2/raaijmakers_group/RIL2019_analysis")

# prepare the combined contig depth table--------------- 
## 'paste' command has been used on sorted tables in bash, remove the extra contig_id column here
#combined_contig_depth <- read.table("final_assembly/depth/combined_contig.tsv", header = T, sep = "\t")
#even_indexes<-seq(2,228,2)
#combined_contig_depth <- combined_contig_depth[, c(1, even_indexes)]
#write.table(combined_contig_depth, file = "final_assembly/depth/clean_combined_contig.tsv",sep = "\t", row.names = F, col.names = T, quote = F)

#combine_df <- read.delim("final_assembly/depth/clean_combined_contig.tsv")
#taxonomy <- read.delim("final_assembly/kraken/2020_april/contig_10k_taxonomy.tsv", header = F)
#colnames(taxonomy) <- c("contig_id", "taxon", "length")
#overview <- full_join(taxonomy,combine_df, by="contig_id")
#rownames(overview) <- overview$contig_id
#sum(is.na(overview))

#write.table(overview, file = "final_assembly/feature_selection/raw_contig_taxon_depth_overview", sep = "\t", col.names = T, row.names = F, quote = F)

###------
overview <- read.delim("final_assembly/feature_selection/raw_contig_taxon_depth_overview")
rownames(overview) <- overview$contig_id

# extract depth data
combine_contig_counts <- overview[, -c(1:3)]

# extract phenotype data
sample_id <- colnames(combine_contig_counts)
accession_metadata <- data.frame(group=c(rep("RIL", 96), rep("BULK",7), rep("modern", 6), rep("wild",5)), row.names=sample_id)

# extract feature data
contig_taxonomy <- overview[, c(1:3)]

# remove contigs that are absent in all samples (zero)
combine_contigs_to_keep <- rowSums(combine_contig_counts) > 0
sum(combine_contigs_to_keep)

## create a MRexperiment object
## AnnotatedDataFrames needed by metagenomeSeq MR experiments
phenotype_data <- AnnotatedDataFrame(accession_metadata)
feature_data <- AnnotatedDataFrame(contig_taxonomy)
# row.names(feature_data@data) <- feature_data@data$contig_id

## create the MR exp with remaining contigs
combine_contig_mr_exp <- newMRexperiment(
  combine_contig_counts[combine_contigs_to_keep, ],
  phenoData = phenotype_data,
  featureData = feature_data
)

head(libSize(combine_contig_mr_exp))
expSummary(combine_contig_mr_exp)

# metagenomeSeq cumulative sum scaling (CSS)
css_quantile <- cumNormStat(combine_contig_mr_exp)
contig_css <- cumNorm(combine_contig_mr_exp, css_quantile)
#contig_css@phenoData@data
# zero-inflated gaussian mixture model
accession_group <- pData(contig_css)$group
normalisation_factor <- normFactors(contig_css)
normalisation_factor <- log2(normalisation_factor / median(normalisation_factor) + 1)
model_mtx <- model.matrix(~ accession_group + normalisation_factor)
model_control <- zigControl(maxit = 20, verbose = T)
model_fit <- fitZig(contig_css, model_mtx, control = model_control, useCSSoffset = T)
effective_samples <- calculateEffectiveSamples(model_fit)
effective_samples_minimum <- ave(effective_samples)[[1]]

# statistics of effective sample size
#plot(effective_samples)
#plot(density(effective_samples))
#max(effective_samples)
#min(effective_samples)
#mean(effective_samples)

# filter rare contigs based on the effective sample size
rare_contigs <- unname(which(rowSums(MRcounts(contig_css)) < effective_samples_minimum))
contig_css_common <- contig_css[-rare_contigs, ]

normalised_counts_no_log <- MRcounts(contig_css, norm = T)[-rare_contigs,]
dim(normalised_counts_no_log)
#plot(normalised_counts[1,])
#export normalised common contig counts
# exportMat(normalised_counts, file = file.path("final_assembly/feature_selection/results/", "norm_common_contig_conunts.ts"))



#export normalised common contig counts with log2 transformation
# exportMat(contig_css_common_normalised, file = file.path("final_assembly/feature_selection/results/", "norm_log_common_contig_conunts.tsv"))

# ZIG on rare-filtered ASVs to find rhizosphere ASVs
fz_model <- model.matrix(~ accession_group)
#colnames(fz_model) <- c("BULK", "modern", "RIL","wild")
colnames(fz_model) <- levels(accession_group)
fz <- fitZig(contig_css_common, fz_model, control = model_control)

contrast_matrix_modern <- makeContrasts(BULK - modern, levels = fz@fit$design)
contrast_matrix_wild <- makeContrasts(BULK - wild, levels = fz@fit$design)
contrast_matrix_parent <- makeContrasts(BULK - (modern + wild), levels = fz@fit$design)
contrast_matrix_all <- makeContrasts(BULK - (modern + RIL + wild ), levels = fz@fit$design)
contrast_matrix_ril <- makeContrasts(BULK - RIL, levels = fz@fit$design)

fz_contrast_modern <- contrasts.fit(fz@fit, contrast_matrix_modern)
fz_contrast_wild <- contrasts.fit(fz@fit, contrast_matrix_wild)
fz_contrast_parent <- contrasts.fit(fz@fit, contrast_matrix_parent)
fz_contrast_all <- contrasts.fit(fz@fit, contrast_matrix_all)
fz_contrast_ril <- contrasts.fit(fz@fit, contrast_matrix_ril)

fz_bayes_modern <- eBayes(fz_contrast_modern)
fz_bayes_wild <- eBayes(fz_contrast_wild)
fz_bayes_parent <- eBayes(fz_contrast_parent)
fz_bayes_all <- eBayes(fz_contrast_all)
fz_bayes_ril <- eBayes(fz_contrast_ril)

head(fz_bayes_ril$t)

fz_padj_modern <- p.adjust(fz_bayes_modern$p.value, method = "hochberg")
fz_padj_wild <- p.adjust(fz_bayes_wild$p.value, method = "hochberg")
fz_padj_parent <- p.adjust(fz_bayes_parent$p.value, method = "hochberg")
fz_padj_all <- p.adjust(fz_bayes_all$p.value, method = "hochberg")
fz_padj_ril <- p.adjust(fz_bayes_ril$p.value, method = "hochberg")

# select significantly rhizosphere-enriched contigs in each comparison
# (t < 0: bulk -  rhizosphere < 0 )
modern_sig_enriched <- fz_padj_modern <= 0.01 & fz_bayes_modern$t < 0
sum(modern_sig_enriched)

wild_sig_enriched <- fz_padj_wild <= 0.01 & fz_bayes_wild$t < 0
sum(wild_sig_enriched)

parent_sig_enriched <- fz_padj_parent <= 0.01 & fz_bayes_parent$t < 0
sum(parent_sig_enriched)

all_sig_enriched <- fz_padj_all <= 0.01 & fz_bayes_all$t < 0
sum(all_sig_enriched)

ril_sig_enriched <- fz_padj_ril <= 0.01 & fz_bayes_ril$t < 0
sum(ril_sig_enriched)

sum(modern_sig_enriched | wild_sig_enriched | parent_sig_enriched | all_sig_enriched | ril_sig_enriched)

### rhizoshere contigs: enriched in the rhizosphere in either comparison 
# (BULK vs. modern | BULK vs. wild | BULK vs. parent | BULK vs. RIL | BULK vs. all)
contig_css_filtered_rhi <- contig_css_common[modern_sig_enriched | wild_sig_enriched | parent_sig_enriched | all_sig_enriched | ril_sig_enriched,]
write.table(rownames(contig_css_filtered_rhi), file = "final_assembly/feature_selection/20200511/rhi_contig_id.txt", sep = "\n", row.names = F, col.names = F, quote = F)

### calculate log 2 fold change
# 1.from un-normalised data
raw_ril_counts <- contig_css_filtered_rhi@assayData$counts[,accession_metadata$group == "RIL"]
raw_bulk_counts <- contig_css_filtered_rhi@assayData$counts[,accession_metadata$group == "BULK"]
0 %in% apply(raw_bulk_counts,1,mean)
0 %in% apply(raw_ril_counts,1,mean)
# no zero is present in the values
df_raw_avg_counts <- data.frame(apply(raw_bulk_counts,1,mean) , apply(raw_ril_counts,1,mean))
colnames(df_raw_avg_counts) <- c("bulk", "ril")
raw_RIL_BULK_log2FC <- log2(df_raw_avg_counts$ril/df_raw_avg_counts$bulk )
df_raw_avg_counts$log2FC <- raw_RIL_BULK_log2FC
sum(raw_RIL_BULK_log2FC > 2)

plot(df_raw_avg_counts$bulk,df_raw_avg_counts$ril, pch=24, cex=0.5)
points(df_raw_avg_counts$bulk[df_raw_avg_counts$log2FC>2],df_raw_avg_counts$ril[df_raw_avg_counts$log2FC>2], pch=24, cex=0.5, col="pink")

#contig_css_filtered_rhi_4_folds <- contig_css_filtered_rhi[raw_RIL_BULK_log2FC > 2, ]
#rhi_4_folds_contig_id <- rownames(contig_css_filtered_rhi_4_folds)

### calculate log 2 fold change
# 2.from normalised data without log2 transformation
contig_css_filtered_rhi_norm_no_log <- MRcounts(contig_css_filtered_rhi, norm = TRUE, log = FALSE) %>%
  na.omit()

norm_no_log_ril_counts <- contig_css_filtered_rhi_norm_no_log[,accession_metadata$group == "RIL"]
norm_no_log_bulk_counts <- contig_css_filtered_rhi_norm_no_log[,accession_metadata$group == "BULK"]
0 %in% apply(norm_no_log_ril_counts,1,mean)
0 %in% apply(norm_no_log_bulk_counts,1,mean)
# no zero is present in the values
df_norm_no_log_avg_counts <- data.frame(apply(norm_no_log_bulk_counts,1,mean) , apply(norm_no_log_ril_counts,1,mean))
colnames(df_norm_no_log_avg_counts) <- c("bulk", "ril")
norm_no_log_RIL_BULK_log2FC <- log2(df_norm_no_log_avg_counts$ril/df_norm_no_log_avg_counts$bulk )
df_norm_no_log_avg_counts$log2FC <- norm_no_log_RIL_BULK_log2FC
sum(norm_no_log_RIL_BULK_log2FC > 2)


#### plot ---------------------------
plot(df_norm_no_log_avg_counts$bulk,df_norm_no_log_avg_counts$ril, cex=0.8, xlab = "BULK", ylab = "RIL")
plot(df_norm_no_log_avg_counts$bulk,df_norm_no_log_avg_counts$ril, cex=0.8,xlab = "BULK", ylab = "RIL", ylim = c(0,10))

points(df_norm_no_log_avg_counts$bulk[df_norm_no_log_avg_counts$log2FC>2],df_norm_no_log_avg_counts$ril[df_norm_no_log_avg_counts$log2FC>2], cex=0.8, col="green")

hist(df_norm_no_log_avg_counts$ril, breaks = 1000, xlim = c(0,10))

select_points <- (df_norm_no_log_avg_counts$bulk>1|df_norm_no_log_avg_counts$ril>1) & (df_norm_no_log_avg_counts$log2FC>2)
sum(df_norm_no_log_avg_counts$bulk>1|df_norm_no_log_avg_counts$ril>1)
sum(df_norm_no_log_avg_counts$log2FC>2)
sum(select_points)

points(df_norm_no_log_avg_counts$bulk[select_points], df_norm_no_log_avg_counts$ril[select_points],cex=0.8, col="blue")

selected_rhi_contig_id <- rownames(df_norm_no_log_avg_counts)[select_points]
write.table(selected_rhi_contig_id, file = "final_assembly/feature_selection/20200511/filtered_rhi_contig_id", sep = "\n", row.names = F, col.names = F, quote = F)


#### check contig length, taxonomy ----------------------------------
#sum(rownames(contig_css_filtered_rhi)[select_points] == rownames(df_norm_no_log_avg_counts)[select_points])
contig_css_filtered_rhi_filter <- contig_css_filtered_rhi[select_points]
taxon_filtered_rhi <- as.character(contig_css_filtered_rhi_filter@featureData@data$taxon)
len_filtered_rhi <-contig_css_filtered_rhi_filter@featureData@data$length
summary(len_filtered_rhi)
#hist(len_filtered_rhi,breaks = 1249)#, xlim = c(0,50000))
plot(density(len_filtered_rhi))
feature_filtered_rhi <- cbind(selected_rhi_contig_id,taxon_filtered_rhi,len_filtered_rhi)
write.table(feature_filtered_rhi, file = "final_assembly/feature_selection/20200511/filtered_rhi_feature", sep = "\t", row.names = F, col.names = F, quote = F)

contig_css_filtered_rhi_filter_count <- MRcounts(contig_css_filtered_rhi_filter, norm = TRUE, log = TRUE) %>%
  na.omit()
dim(contig_css_filtered_rhi_filter_count)
write.table(contig_css_filtered_rhi_filter_count, file = "final_assembly/feature_selection/20200511/filtered_rhi_count", sep = "\t", row.names = T, col.names = T, quote = F)

#norm_counts_long <- melt(contig_css_filtered_rhi_norm_no_log)
#head(norm_counts_long)
#colnames(norm_counts_long) <- c("contig_id","sample_id","normalised_counts" )
#plot(norm_counts_long$normalised_counts)

### calculate log 2 fold change
# 3.from normalised data with log2 transformation
contig_css_filtered_rhi_norm_log <- MRcounts(contig_css_filtered_rhi, norm = TRUE, log = TRUE) %>%
  na.omit()

norm_log_ril_counts <- contig_css_filtered_rhi_norm_log[,accession_metadata$group == "RIL"]
norm_log_bulk_counts <- contig_css_filtered_rhi_norm_log[,accession_metadata$group == "BULK"]
0 %in% apply(norm_log_ril_counts,1,mean)
0 %in% apply(norm_log_bulk_counts,1,mean)
# no zero is present in the values
df_norm_log_avg_counts <- data.frame(apply(norm_log_bulk_counts,1,mean) , apply(norm_log_ril_counts,1,mean))
colnames(df_norm_log_avg_counts) <- c("bulk", "ril")
norm_log_RIL_BULK_log2FC <- log2(df_norm_log_avg_counts$ril/df_norm_log_avg_counts$bulk )
df_norm_log_avg_counts$log2FC <- norm_log_RIL_BULK_log2FC
sum(norm_log_RIL_BULK_log2FC > 2)

plot(df_norm_log_avg_counts$bulk,df_norm_log_avg_counts$ril, pch=3, cex=0.5)
points(df_norm_log_avg_counts$bulk[df_norm_log_avg_counts$log2FC>2],df_norm_log_avg_counts$ril[df_norm_log_avg_counts$log2FC>2], pch=3, cex=0.5, col="pink")

### visualization
#df_raw_avg_counts$type <- "raw"
#df_norm_no_log_avg_counts$type <- "norm_no_log"
#df_norm_log_avg_counts$type <- "norm_log"


#### plot 
par(mfcol=c(1,1))
raw_bulk <- hist(log2(df_raw_avg_counts$bulk))
raw_ril <- hist(log2(df_raw_avg_counts$ril))
plot(raw_ril, col=rgb(0,0,1,1/4), xlim = c(-10,10), xlab = "log2(raw counts)", main="Histogram of log2(raw counts)")
plot(raw_bulk, col=rgb(1,0,0,1/4), add= T)
legend(-10, 600, legend=c("RIL", "BULK"),
       fill=c(rgb(0,0,1,1/4), rgb(1,0,0,1/4)), cex = 0.8)

norm_no_log_bulk <- hist(log2(df_norm_no_log_avg_counts$bulk))
norm_no_log_ril <- hist(log2(df_norm_no_log_avg_counts$ril))
plot(norm_no_log_ril, col=rgb(0,0,1,1/4), xlim = c(-10,10), xlab = "log2(normalised counts)", main="Histogram of log2(normalised counts)")
plot(norm_no_log_bulk, col=rgb(1,0,0,1/4), add= T)
legend(-10, 600, legend=c("RIL", "BULK"),
       fill=c(rgb(0,0,1,1/4), rgb(1,0,0,1/4)), cex = 0.8)

norm_log_bulk <- hist(log2(df_norm_log_avg_counts$bulk))
norm_log_ril <- hist(log2(df_norm_log_avg_counts$ril))
plot(norm_log_ril, col=rgb(0,0,1,1/4), xlim = c(-10,10), xlab = "log2(normalised counts to log2 transform scale)", main="Histogram of log2(normalised counts to log2 transform scale)")
plot(norm_log_bulk, col=rgb(1,0,0,1/4), add= T)
legend(-10, 600, legend=c("RIL", "BULK"),
       fill=c(rgb(0,0,1,1/4), rgb(1,0,0,1/4)), cex = 0.8)

### make dataframe and export
###

modern_contig_id <- rownames(contig_css_filtered_modern)
wild_contig_id <- rownames(contig_css_filtered_wild)
parent_contig_id <- rownames(contig_css_filtered_parent)
ril_contig_id <- rownames(contig_css_filtered_ril)
all_contig_id <- rownames(contig_css_filtered_all)
#length(union(modern_contig_id, wild_contig_id) %>% union(parent_contig_id) %>% union(ril_contig_id) %>% union(all_contig_id))


# final normalised and filtered contig counts
contig_css_filtered_modern <- contig_css_common[modern_sig_enriched, ]
contig_css_filtered_wild <- contig_css_common[wild_sig_enriched, ]
contig_css_filtered_parent <- contig_css_common[parent_sig_enriched, ]
contig_css_filtered_all <- contig_css_common[all_sig_enriched, ]
contig_css_filtered_ril <- contig_css_common[ril_sig_enriched, ]

contig_css_common_normalised <- MRcounts(contig_css_common, norm = TRUE, log = TRUE) %>%
  na.omit()
contig_css_filtered_normalised_modern <- MRcounts(contig_css_filtered_modern, norm = TRUE, log = TRUE) %>%
  na.omit()
contig_css_filtered_normalised_wild <- MRcounts(contig_css_filtered_wild, norm = TRUE, log = TRUE) %>%
  na.omit()
contig_css_filtered_normalised_parent <- MRcounts(contig_css_filtered_parent, norm = TRUE, log = TRUE) %>%
  na.omit()
contig_css_filtered_normalised_all <- MRcounts(contig_css_filtered_all, norm = TRUE, log = TRUE) %>%
  na.omit()
contig_css_filtered_normalised_ril <- MRcounts(contig_css_filtered_ril, norm = TRUE, log = TRUE) %>%
  na.omit()
contig_css_filtered_normalised_union <- MRcounts(contig_css_filtered_union, norm = TRUE, log = TRUE) %>%
  na.omit()
##! contig_css_filtered_normalised_rhi_4_folds <- MRcounts(contig_css_filtered_rhi_4_folds, norm = TRUE, log = TRUE) %>%
  ##na.omit()
contig_css_filtered_normalised_rhi <- MRcounts(contig_css_filtered_rhi, norm = TRUE, log = TRUE) %>%
  na.omit()

## calculate fold change on normalised counts
norm_ril_counts <- contig_css_filtered_normalised_rhi[,accession_metadata$group == "RIL"]
norm_bulk_counts <- contig_css_filtered_normalised_rhi[,accession_metadata$group == "BULK"]
0 %in% apply(norm_bulk_counts,1,mean)
0 %in% apply(norm_ril_counts,1,mean)
# no zero is present in the values
norm_RIL_BULK_fold_change <- log2( apply(norm_ril_counts,1,mean) / apply(norm_bulk_counts,1,mean) )
sum(norm_RIL_BULK_fold_change > 2)

df_norm_abudance_profile <- data.frame(apply(norm_bulk_counts,1,mean), apply(norm_ril_counts,1,mean))
colnames(df_norm_abudance_profile) <- c("bulk","ril")
plot(df_norm_abudance_profile$bulk, df_norm_abudance_profile$ril)
points(df_norm_abudance_profile$bulk[df_norm_abudance_profile$bulk>1 |df_norm_abudance_profile$ril>1], df_norm_abudance_profile$ril[df_norm_abudance_profile$bulk>1 |df_norm_abudance_profile$ril>1 ], col="green")
points(df_norm_abudance_profile$bulk[norm_RIL_BULK_fold_change>2], df_norm_abudance_profile$ril[norm_RIL_BULK_fold_change>2], col="blue")
points(df_norm_abudance_profile$bulk[select_fc_meaningful], df_norm_abudance_profile$ril[select_fc_meaningful], col="red")

select_fc_meaningful <- norm_RIL_BULK_fold_change > 2 & (df_norm_abudance_profile$bulk>1 |df_norm_abudance_profile$ril>1)
sum(select_fc_meaningful)
sum(norm_RIL_BULK_fold_change > 2 & (df_norm_abudance_profile$bulk>1 |df_norm_abudance_profile$ril>1))

hist(log(df_norm_abudance_profile$bulk))
hist(log(df_norm_abudance_profile$ril))


dim(contig_css_filtered_normalised_modern)
dim(contig_css_filtered_normalised_wild)
dim(contig_css_filtered_normalised_parent)
dim(contig_css_filtered_normalised_all)
dim(contig_css_filtered_normalised_ril)
dim(contig_css_filtered_normalised_union)




#------------
# plot heatmap
pdf(file = "heatmap_sig_rhi_count.pdf", width = 12, height = 12)
mat_contig_css_filtered_rhi_norm_log <- as.matrix(contig_css_filtered_rhi_norm_log)
heatmap(mat_contig_css_filtered_rhi_norm_log, scale = "row",cexRow = 0.2, cexCol = 0.5)
dev.off()
#heatmap.2(mat_contig_css_filtered_rhi_norm_log,scale = "none", cexRow = 0.2, cexCol = 0.5,
          #trace = "none", density.info = "none")
#plot(contig_css_filtered_normalised_ril_e[1,])
#sum(rownames(contig_css_filtered_normalised1) %in% rownames(contig_css_filtered_normalised2))

#rhizoshere_contig_numerical_code1 <- rownames(contig_css_filtered_normalised1)
#rhizoshere_contig_numerical_code2 <- rownames(contig_css_filtered_normalised2)
#rhizoshere_contig_numerical_code3 <- rownames(contig_css_filtered_normalised3)
#rhizoshere_contig_numerical_code4 <- rownames(contig_css_filtered_normalised4)

#rhizosphere_contig_BvsM <- overview[rownames(overview) %in% rhizoshere_contig_numerical_code1, "contig_id"]
#rhizosphere_contig_BvsW <- overview[rownames(overview) %in% rhizoshere_contig_numerical_code2, "contig_id"]
#rhizosphere_contig_BvsMW <- overview[rownames(overview) %in% rhizoshere_contig_numerical_code3, "contig_id"]
#rhizosphere_contig_BvsMWR <- overview[rownames(overview) %in% rhizoshere_contig_numerical_code4, "contig_id"]

#write.table(rhizosphere_contig_BvsM, file = "final_assembly/feature_selection/Rhizosphere_contig_BvsM.txt", sep = "\n", row.names = F, col.names = F, quote = F)
#write.table(rhizosphere_contig_BvsW, file = "final_assembly/feature_selection/Rhizosphere_contig_BvsW.txt", sep = "\n", row.names = F, col.names = F, quote = F)
#write.table(rhizosphere_contig_BvsMW, file = "final_assembly/feature_selection/Rhizosphere_contig_BvsMW.txt", sep = "\n", row.names = F, col.names = F, quote = F)
#write.table(rhizosphere_contig_BvsMWR, file = "final_assembly/feature_selection/Rhizosphere_contig_BvsMWR.txt", sep = "\n", row.names = F, col.names = F, quote = F)

