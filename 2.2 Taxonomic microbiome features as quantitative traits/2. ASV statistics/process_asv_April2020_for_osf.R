#!/usr/bin/env Rscript

# program setup -----------------------------------------------------------
suppressPackageStartupMessages({
  library(tidyverse)
  library(ggpubr)
  library(RColorBrewer)
  library(cowplot)
  library(metagenomeSeq)
  library(magrittr)
  library(dplyr)
  library(phyloseq)
  library(gridExtra)
})

theme_set(theme_pubr(border = T, legend = "right"))
set.seed(56753)


# data loading ------------------------------------------------------------
# All ASV from RIL experiment, but also from the 8 genotypes, with Chloroplast and Mitochondria still present
asv_output <- read.delim("result_April2020/asv_16s.csv")
# asv_18s_output <- read.delim("asv_18s.csv")

# make the taxonomy table by taking the first columns, and rename Rhizobium to an easier name
asv_taxonomy <- asv_output[, 1:8] %>%
  mutate(genus = sub("Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium", "Rhizobium", genus, fixed = T))
asv_counts <- as.matrix(asv_output[, -c(1:8)])
accession_metadata <- data.frame(accession = colnames(asv_counts)) %>%
  mutate(
    accession = as.character(accession),
    group = as.factor(case_when(
      grepl("_P$", accession, ignore.case = T) ~ "Pink",
      grepl("_Y$", accession, ignore.case = T) ~ "Yellow",
      grepl("_BULK", accession, ignore.case = T) ~ "Bulk",
      grepl("_MM[\\d]$", accession, ignore.case = T, perl = T) ~ "Modern",
      grepl("_P[\\d]$", accession, ignore.case = T, perl = T) ~ "Wild",
      TRUE ~ "Other"
    )),
    sample_id = as.numeric(sub("^[X]?([\\d]+).*", "\\1", accession, perl = T)),
    ril_id = as.numeric(ifelse(group %in% c("Pink", "Yellow"), sub(".*_([\\d]+).*", "\\1", accession, perl = T), NA))
  )
row.names(accession_metadata) <- accession_metadata$accession

####################################################
# Seperate the 8 genotypes from the RIL experiment #
####################################################

# ASV pre-filtering -------------------------------------------------------
# remove other genotypes; Y:253,226,259,275; P:217,237,249 (too long to germinate)
# Low reads X94_267_P   X190_254_Y   X124_245_Y X52_238_P   X148_256_Y
# sort(colSums(asv_css_filtered_normalised), decreasing = TRUE)

asv_output[, -c(1:8)] %>% colSums() %>% sort(decreasing = TRUE) %>% plot()
mean_seq_depth <- asv_output[, -c(1:8)] %>% colSums() %>% mean()
sd_seq_depth <-asv_output[, -c(1:8)] %>% colSums() %>% sd()
abline(h= mean_seq_depth-(2*sd_seq_depth))

Samples_to_remove <- names(which(asv_output[, -c(1:8)] %>% colSums() < mean_seq_depth-(2*sd_seq_depth)))

cat("The following samples were removed from the analysis because of low sequencing depth:\n", Samples_to_remove,"\n")

Grew_long_Y <-c(253, 226, 259, 275)
Grew_long_P <-c(217, 237, 249)
# Based on 2 sd cut-off: # 124_245_Y 148_256_Y 190_254_Y  52_238_P  94_267_P

Low_reads_Y <-c(245, 256, 254)
Low_reads_P <-c(238, 267)

genotypes_to_keep_RIL <- !(
  accession_metadata$group == "Other" |
  (accession_metadata$group == "Yellow" & accession_metadata$ril_id %in% c(Grew_long_Y,Low_reads_Y)) |
  (accession_metadata$group == "Pink" & accession_metadata$ril_id %in% c(Grew_long_P, Low_reads_P))
)

# remove 8 genotypes samples
asv_counts_ril <- asv_counts[,genotypes_to_keep_RIL]
accession_metadata_ril <- accession_metadata[genotypes_to_keep_RIL, ] %>%
  mutate(group = droplevels(group))

# remove ASVs that are chloroplasts, mitochondria, or not bacteria
asvs_to_keep <- !(
  is.na(asv_taxonomy$kingdom) |
  (asv_taxonomy$kingdom != "Bacteria") |
  (!is.na(asv_taxonomy$order) & asv_taxonomy$order == "Chloroplast") |
  (!is.na(asv_taxonomy$family) & asv_taxonomy$family == "Mitochondria")
)

# remove zero-count ASVs
asvs_to_keep <- asvs_to_keep & (rowSums(asv_counts_ril) > 0)

# filter ASV table
asv_taxonomy_f <- asv_taxonomy[asvs_to_keep, ]
# asv_taxonomy_f$asv <- 1:nrow(asv_taxonomy_f)
# saveRDS(asv_taxonomy_f, file = "result_April2020/asv_taxonomy.rds", compress = "xz")

########################################################################################################################
#

# ASV CSS normalisation ---------------------------------------------------
# ADFs needed by metagenomeSeq MR experiments
phenotype_data <- AnnotatedDataFrame(data.frame(group = accession_metadata_ril$group, row.names = accession_metadata_ril$accession))
feature_data <- AnnotatedDataFrame(asv_taxonomy_f)
# feature_data@data$asv <- 1:nrow(feature_data)
row.names(feature_data@data) <- feature_data@data$asv
# create the MR exp with remaining ASVs
asv_mr_exp <- newMRexperiment(
  asv_counts_ril[asvs_to_keep, ],
  phenoData = phenotype_data,
  featureData = feature_data
)

# metagenomeSeq cumulative sum scaling (CSS)
css_quantile <- cumNormStat(asv_mr_exp)
asv_css <- cumNorm(asv_mr_exp, css_quantile)

# zero-inflated gaussian mixture model
accession_group <- pData(asv_css)$group
normalisation_factor <- normFactors(asv_css)
normalisation_factor <- log2(normalisation_factor / median(normalisation_factor) + 1)
model_mtx <- model.matrix(~ accession_group + normalisation_factor)
model_control <- zigControl(maxit = 20, verbose = T)
model_fit <- fitZig(asv_css, model_mtx, control = model_control, useCSSoffset = T)
effective_samples <- calculateEffectiveSamples(model_fit)
effective_samples_minimum <- ave(effective_samples)[[1]]
# filter rare ASVs based on the effective sample size
rare_asvs <- unname(which(rowSums(MRcounts(asv_css)) < effective_samples_minimum))
normalised_counts <- as.data.frame(MRcounts(asv_css, norm = T))[-rare_asvs,]

# ZIG on rare-filtered ASVs to find rhizosphere ASVs
asv_css_common <- asv_css[-rare_asvs, ]
fz_model <- model.matrix(~ accession_group)
colnames(fz_model) <- levels(accession_group)
fz <- fitZig(asv_css_common, fz_model, control = model_control)

contrast_matrix <- makeContrasts(Bulk - (Yellow + Pink + Modern + Wild), levels = fz@fit$design)
fz_contrast <- contrasts.fit(fz@fit, contrast_matrix)
fz_bayes <- eBayes(fz_contrast)
fz_padj <- p.adjust(fz_bayes$p.value, method = "hochberg")

# Get negative t scores, and significant adjusted p-values
rhizosphere_enriched_pos <- intersect(which(fz_bayes$t<0), which(fz_padj <= 0.01))

# Rhizosphere enriched
asv_filtered_rhizosphere <- asv_css_common[rhizosphere_enriched_pos,]

# Calculate log(abundance, 2) fold change
RIL_BULK_fold_change1 <- log((apply(asv_filtered_rhizosphere@assayData$counts[,accession_metadata_ril$group%in% c("Pink","Yellow") ], 1, mean)+1) / (apply(asv_filtered_rhizosphere@assayData$counts[,accession_metadata_ril$group=="Bulk"],1,mean)+1),2)
plot(sort(RIL_BULK_fold_change1))
abline(h=1)

# Identify which have a log(abundance, 2) fold change >2
# which(RIL_BULK_fold_change>2) %in%
length(which(RIL_BULK_fold_change1>1))

# Log2 fold enriched rhizosphere taxa
asv_enriched_rhizosphere <- asv_filtered_rhizosphere[which(RIL_BULK_fold_change1>1) ,]

#  asv_css_common_normalised <- MRcounts(asv_css_common, norm = TRUE, log = TRUE) %>%  na.omit()

#
asv_css_enriched_rhizosphere <- MRcounts(asv_enriched_rhizosphere, norm = TRUE, log = TRUE) %>%
  na.omit()


rhizoshere_asv_numerical_code <- rownames(asv_css_enriched_rhizosphere)

write.table(paste("asv",rhizoshere_asv_numerical_code,sep=""), file = "result_April2020/Rhizosphere_ASV.txt", sep = "\n", row.names = F, col.names = T, quote = F)

cat(sprintf(
  "%.1f%% ASVs filtered",
  100 - nrow(asv_css_enriched_rhizosphere) / nrow(asv_output) * 100
))



########################################################################
#
########################################################################

# make cut-offs for: stochastic, core, and variable -----------------------
core_cutoff <- 0.9
rare_cutoff <- 0.5

css_statistics_rhiz <- as.data.frame(t(apply(asv_css_enriched_rhizosphere[, accession_metadata_ril$group != "Bulk"], 1, function(x){c(
  "zeroes" = sum(x == 0),
  "nz" = sum(x > 0),
  "mean" = mean(x, na.rm = T),
  "nz_mean" = mean(x[x > 0]),
  "sum" = sum(x, na.rm = T),
  "n" = length(x),
  "min" = min(x, na.rm = T),
  "max" = max(x, na.rm = T),
  "norm_sd" = sd(x, na.rm = T) / max(x, na.rm = T),
  "norm_sd_nz" = sd(x[x > 0]) / max(x, na.rm = T)
)})), row.names = row.names(asv_css_enriched_rhizosphere)) %>%
  rownames_to_column("asv") %>%
  mutate_all(~replace(., is.na(.), 0)) %>%
  mutate(
    asv = as.numeric(asv),
    ecology = as.factor(case_when(
      nz / n > core_cutoff ~ "Core",
      nz / n < rare_cutoff ~ "Stochastic",
      TRUE ~ "Flexible"
    ))
  ) #%>% left_join(asv_taxonomy, by = "asv")


css_statistics_rhiz$asv[which(css_statistics_rhiz$ecology %in% c("Core", "Flexible"))]

cf_asv_pos_in_counts_table <- rownames(asv_mr_exp@assayData$counts) %in% css_statistics_rhiz$asv[which(css_statistics_rhiz$ecology %in% c("Core", "Flexible"))]
cf_ratio_of_raw_counts <- colSums(asv_mr_exp@assayData$counts[cf_asv_pos_in_counts_table,]) / colSums(asv_mr_exp@assayData$counts)
plot(cf_ratio_of_raw_counts)
plot(cf_ratio_of_raw_counts ~ accession_group)

ave(cf_ratio_of_raw_counts[accession_group=="Bulk"])[1]
ave(cf_ratio_of_raw_counts[accession_group%in%c("Pink","Yellow","Wild","Pink")])[1]

# Figure 2 b
# visualise groups
ggplot(css_statistics_rhiz, aes(x = norm_sd_nz, y = nz / n, size = (sum/max(sum))^2, fill = ecology)) +
  geom_point(alpha = 0.25, shape = 21, color = "black") +
  geom_hline(yintercept = core_cutoff, linetype = "dotted") +
  annotate("text", x = 0, y = core_cutoff, label = sprintf("%.0f%%", core_cutoff*100), hjust = -0.2, vjust = -0.4) +
  geom_hline(yintercept = rare_cutoff, linetype = "dotted") +
  annotate("text", x = 0, y = rare_cutoff, label = sprintf("%.0f%%", rare_cutoff*100), hjust = -0.2, vjust = -0.4) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0, 0, 0.05), labels = scales::percent_format(accuracy = 1)) +
  labs(
    title = "ASV separation based on presence/absence",
    x = "Normalised standard deviation",
    y = "% non-zero counts",
    size = "norm.sum^2"
  )


dev.off()
###
p5 <- plotOrd(asv_css_enriched_rhizosphere[,which(accession_group!="Bulk")], usePCA = T, useDist = T, bg = accession_group[which(accession_group!="Bulk")], pch = 21,  main = "All Rhizosphere ASV")
legend("topright", levels(accession_group)[c(2,4,3,5)], text.col = c(2,4,3,5), cex=0.75)

# Exporting of tables -----------------------------------------------------

# Print out asv counts (after removing mitochondria etc.) and taxonomy table for RIL experiment
# The core rhizosphere
# write.table(paste("asv",css_statistics_rhiz$asv[which(css_statistics_rhiz$ecology=="Core")],sep=""), file = "result_April2020/Rhizosphere_Core_ASV.txt", sep = "\n", row.names = F, col.names = T, quote = F)

# Supplemental table # 5085 ASV
# Raw counts for all ASV except chloroplast and mitochondria and 0's
as.data.frame(asv_mr_exp@assayData$counts) %>%
  rowid_to_column("aasv") %>%
  write.table(., file = "result_April2020/supplemental_tables/TS1_asv_counts.csv", sep = "\t", row.names = F, col.names = substring(colnames(.), 2), quote = F)

# Supplemental table # Taxonomy for all asv
write.table(asv_mr_exp@featureData@data, file = "result_April2020/supplemental_tables/TS2_asv_taxonomy.csv", sep = "\t", row.names = F, col.names = T, quote = F)

# CSS Normalized, filtered for effective sample size
as.data.frame(asv_css_common_normalised) %>%
  rownames_to_column("aasv") %>%
  write.table(., file = "result_April2020/supplemental_tables/TS3_asv_normalized_filtered.csv", sep = "\t", row.names = F, col.names = substring(colnames(.), 2), quote = F)

# CSS Normalized, filtered for effective sample size, and Rhizosphere ASV (benjamini hockberg 0.01), log2 fold change >2
as.data.frame(asv_css_enriched_rhizosphere) %>%
  rownames_to_column("aasv") %>%
  write.table(., file = "result_April2020/supplemental_tables/TS4_asv_rhizosphere.csv", sep = "\t", row.names = F, col.names = substring(colnames(.), 2), quote = F)


# Statistics for rhizosphere ASV
write.table(css_statistics_rhiz, file = "result_April2020/supplemental_tables/TS5_rhizosphere_statistics.csv", sep = "\t", row.names = F, col.names = T, quote = F)

# phyloseq ----------------------------------------------------------------

# First explore the taxonomy of the full dataset
# 'css_common' is all remaining asv after filtering for 'rare' asv based on the 'effective_samples_minimum'
# First make a sample data frame. Order it and create a sample order variable to make things pretty
ps_sampledata <- data.frame(asv_css_common@phenoData@data)
order_samples <- c(which(ps_sampledata$group=="Bulk"),which(ps_sampledata$group=="Modern"),which(ps_sampledata$group=="Wild"),which(ps_sampledata$group=="Pink"),which(ps_sampledata$group=="Yellow"))
new_ps_sampledata <- data.frame(ps_sampledata[order_samples,])
rownames(new_ps_sampledata) <-rownames(ps_sampledata)[order_samples]
colnames(new_ps_sampledata) <- "Treatment"
new_ps_sampledata$sample_order <- 1:213

# Next create the ASV table for Phyloseq
ASV_matrix_for_phyloseq <- asv_css_common@assayData$counts[,as.numeric(order_samples)]
ASV_matrix_for_phyloseq <- MRcounts(asv_css_common, norm = TRUE, log = FALSE) %>%
  na.omit()
# And finally the phyloseq table
taxonomy_matrix_for_phyloseq <- as.matrix(asv_css_common@featureData@data)[,-1]

# a subset for the bulk ASVs
Bulk_ASV_matrix_for_phyloseq <- ASV_matrix_for_phyloseq[-rhizosphere_enriched_pos,]
Bulk_taxonomy_matrix_for_phyloseq <- as.matrix(asv_css_common@featureData@data)[-rhizosphere_enriched_pos,-1]

# a subset for the Rhizosphere ASVs
Rhizo_ASV_matrix_for_phyloseq <- ASV_matrix_for_phyloseq[rhizosphere_enriched_pos,]
Rhizo_taxonomy_matrix_for_phyloseq <- as.matrix(asv_css_common@featureData@data)[rhizosphere_enriched_pos,-1]

# A subset for just the RIL
P_Y_columns <- which(ps_sampledata$group %in% c("Pink","Yellow"))

# Now we create the first figures,
# The top 50 Phylum, Class Family and Genus
# a general overview of all 1700 ASV

ps_all <- phyloseq(
  otu_table(ASV_matrix_for_phyloseq, taxa_are_rows = T),
  sample_data(new_ps_sampledata),
  tax_table(taxonomy_matrix_for_phyloseq)
)

ps_all <- transform_sample_counts(ps_all, function(x) x / sum(x) )

re <- rownames(ps_all@otu_table) %in% as.character(css_statistics_rhiz$asv[css_statistics_rhiz$ecology %in% c("Core","Flexible")])
ps_rhizo <- subset_taxa(ps_all, re)

##################
# Core, Flexible #
##################

unique(ps_rhizo@tax_table[,3])

ps_Alphaproteobacteria_r = subset_taxa(ps_rhizo, class == "Alphaproteobacteria")
ps_Actinobacteria_r = subset_taxa(ps_rhizo, class == "Actinobacteria")
ps_Gammaproteobacteria_r = subset_taxa(ps_rhizo, class == "Gammaproteobacteria")
ps_Bacteroidia_r = subset_taxa(ps_rhizo, class == "Bacteroidia")
ps_Saccharimonadia_r = subset_taxa(ps_rhizo, class == "Saccharimonadia")
ps_Verrucomicrobiae_r = subset_taxa(ps_rhizo, class == "Verrucomicrobiae")
ps_Bacilli_r = subset_taxa(ps_rhizo, class == "Bacilli")
ps_Acidimicrobiia_r = subset_taxa(ps_rhizo, class == "Acidimicrobiia")
ps_SVBA_r = subset_taxa(ps_rhizo, class %in% c("Acidimicrobiia","Bacilli","Verrucomicrobiae","Saccharimonadia"))

bp1_c <- plot_bar(ps_Gammaproteobacteria_r, fill = "genus", facet_grid=~genus) + geom_bar(stat = "identity") + theme(axis.text.x=element_blank(), legend.position = "none", axis.ticks.x=element_blank(), axis.title.x=element_blank()) + guides(colour = FALSE) + ggtitle("Gammaproteobacteria")
bp2_c <- plot_bar(ps_Alphaproteobacteria_r, fill = "genus", facet_grid=~genus) + geom_bar(stat= "identity") + theme(axis.text.x=element_blank(), legend.position = "none", axis.ticks.x=element_blank(), axis.title.x=element_blank()) + guides(colour = FALSE) + ggtitle("Alphaproteobacteria")
bp3_c <- plot_bar(ps_Actinobacteria_r, fill = "genus", facet_grid=~genus) + geom_bar(stat= "identity") + theme(axis.text.x=element_blank(), legend.position = "none", axis.ticks.x=element_blank(), axis.title.x=element_blank()) + guides(colour = FALSE) + ggtitle("Actinobacteria")
bp4_c <- plot_bar(ps_Bacteroidia_r, fill = "genus", facet_grid=~genus) + geom_bar(stat= "identity") + theme(axis.text.x=element_blank(), legend.position = "none", axis.ticks.x=element_blank(), axis.title.x=element_blank()) + guides(colour = FALSE) + ggtitle("Bacteroidia")
# bp5_c <- plot_bar(ps_Saccharimonadia_r, fill = "order", facet_grid=~order) + geom_bar(stat= "identity") + theme(axis.text.x=element_blank(), legend.position = "none", axis.ticks.x=element_blank()) + guides(colour = FALSE) + ggtitle("Rhizosphere \nSaccharimonadia")
# bp6_c <- plot_bar(ps_Verrucomicrobiae_r, fill = "order", facet_grid=~order) + geom_bar(stat= "identity") + theme(axis.text.x=element_blank(), legend.position = "none", axis.ticks.x=element_blank()) + guides(colour = FALSE) + ggtitle("Rhizosphere \nVerrucomicrobiae")
# bp7_c <- plot_bar(ps_Bacilli_r, fill = "order", facet_grid=~order) + geom_bar(stat= "identity") + theme(axis.text.x=element_blank(), legend.position = "none", axis.ticks.x=element_blank()) + guides(colour = FALSE) + ggtitle("Rhizosphere \nBacilli")
# bp8_c <- plot_bar(ps_Acidimicrobiia_r, fill = "order", facet_grid=~order) + geom_bar(stat= "identity") + theme(axis.text.x=element_blank(), legend.position = "none", axis.ticks.x=element_blank()) + guides(colour = FALSE) + ggtitle("Rhizosphere \nAcidimicrobiia")
bp9_c <- plot_bar(ps_SVBA_r, fill = "genus", facet_grid=~order) + geom_bar(stat= "identity") + theme(axis.text.x=element_blank(), legend.position = "none", axis.ticks.x=element_blank(), axis.title.x=element_blank()) + guides(colour = FALSE) + ggtitle("Saccharimonadia, Verrucomicrobiae, Bacilli, Acidimicrobiia")


# A nice figure that is not in the manuscript showing all the different families and genera abundances
grid.arrange(bp1_c, bp2_c, bp3_c, bp4_c, bp9_c, nrow = 5)
# grid.arrange(bp5_c, bp6_c, bp7_c, bp8_c, nrow = 2)

ps_all@sam_data$color <- string.to.colors(ps_all@sam_data$Treatment)

# Figure 2 A (An ugly version, the colors were adjusted in Adobe illustrator for the final manuscript)
ps_all.ord <- ordinate(ps_all, "PCoA", "bray")
p1 = plot_ordination(ps_all, ps_all.ord, type="samples", color="Treatment", title="PCoA, 1703 ASV") + scale_colour_manual(values = c("red", "blue", "pink", "green", "yellow"))
print(p1)

# Same as above, but only using the rhizosphere asv (66 taxa).
ps_rhizo.ord <- ordinate(ps_rhizo, "PCoA", "bray")
p2 = plot_ordination(ps_rhizo, ps_rhizo.ord, type="samples", color="Treatment", title="Accession") + scale_colour_manual(values = c("red", "blue", "pink", "green", "yellow"))
print(p2)
