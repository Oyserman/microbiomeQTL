#!/usr/bin/env Rscript

# program setup -----------------------------------------------------------
  library(tidyverse)
  library(ggpubr)
  library(RColorBrewer)
  library(qtl2)
  library(cowplot)
  library(yaml)
  library(multcompView)
  library(magrittr)

argv <- commandArgs(trailingOnly = T)
theme_set(theme_pubr(border = T, legend = "right"))
max_treads <- 4
setwd("/data/archive/ngs2/ME2/raaijmakers_group/RIL2019_analysis/16S_processing/for_OSF/")

# load R/qtl2 data --------------------------------------------------------

rqtl_phenofile_r <- "./replicated_flexiblecore.csv"

# Load the yaml files for the microQTL
rqtl_rnaseq_map <- "./rnaseq_map.yaml"

# to be deleted?
# rqtl_rnaseq_map_pdw <- "./rnaseq_map_pdw.yaml"


taxonomy_file <- "./asv_taxonomy.rds"
asv_taxonomy <- readRDS(taxonomy_file)

# Treat 0's as NA if 75% or more of samples are non-zero

############ Classic Phenotypes First #############

zeros <- TRUE

# Curate the metadata

metadata_file_r <- "./RIL_metadata_replicates_5_19_2019.txt"

metadata_r <- read.table(file=metadata_file_r, row.names = 1, header = TRUE)[,-3]
metadata_r$Color <- as.character(metadata_r$Color)

# Make the variables numeric
metadata_r$Color[which(metadata_r$Color=="P")] <- 1
metadata_r$Color[which(metadata_r$Color=="Y")] <- 2
metadata_r$Color <- as.numeric(metadata_r$Color)

# Create a separate variable for harvest day
metadata_r_Harvest_Day <- metadata_r$Harvest_Day
names(metadata_r_Harvest_Day) <- rownames(metadata_r)

# Convert to CSV
pdw_phenotypes <- read.table("RIL_metadata_replicates_5_19_2019.txt")
pdw_phenotypes$V2 <- c("Color", metadata_r$Color)
pdw_phenotypes <- pdw_phenotypes[,-4]
# write.table(pdw_phenotypes, "RIL_metadata_replicates_5_19_2019.csv", sep=",", col.names=FALSE, row.names=FALSE, quote=FALSE)


######### Defining outputs names
#########
output_basename_r <- tools::file_path_sans_ext(rqtl_phenofile_r)
output_lodscores_r <- sprintf("%s_lod.csv", output_basename_r)
output_plodscores_r <- sprintf("%s_plod.csv", output_basename_r)
output_peaks_r <- sprintf("%s_peaks.csv", output_basename_r)
output_heritability_r <- sprintf("%s_heritability.csv", output_basename_r)
output_heatmap_r <-sprintf("%s_lod.pdf", output_basename_r)
#########

# Run rQTL2
# pmap_1 <- read_csv("Tomaat_Micro_biome_RNA_Seq_map/pmap.csv")

# microbiome as phenotypes
in_cross_r <- read_cross2(rqtl_rnaseq_map, quiet = F)
in_cross_p <- read_cross2(rqtl_rnaseq_map, quiet = F)
in_cross_y <- read_cross2(rqtl_rnaseq_map, quiet = F)

in_cross_p$pheno <- in_cross_r$pheno[1:89,]
in_cross_y$pheno <- in_cross_r$pheno[90:176,]


# read the yaml files
rqtl_yaml_r <- read_yaml(rqtl_rnaseq_map)

# cat(in_cross)
# check_c <- check_cross2(in_cross_r)

# A variable that may be played with: determining when 0's are informative. Here we choose to treat all 0's as uninformative.
# For metaenomes, this is less important. This a key, practical, discussion point
# minimum_samples_r <- floor(dim(in_cross_r$pheno)[1] / 10 * 2.5)
minimum_samples_r <- dim(in_cross_r$pheno)[1]
minimum_samples_p <- dim(in_cross_p$pheno)[1]
minimum_samples_y <- dim(in_cross_y$pheno)[1]

# minimum_samples_n <- dim(in_cross_r_new$pheno)[1]

# Function to replaces 0's with NA
replace_zero_with_NA <- function(phenotypes, ms) {
missing_data_as_NA <- colSums(phenotypes == 0) <= ms
subset_in_cross_pheno <- phenotypes[,missing_data_as_NA]
subset_in_cross_pheno[subset_in_cross_pheno == 0] <- NA
phenotypes[,missing_data_as_NA] <- subset_in_cross_pheno
return(phenotypes)
}

# replace with zeros with NA, thereby excluding an accession from the analysis
in_cross_r$pheno <- replace_zero_with_NA(in_cross_r$pheno, minimum_samples_r)
in_cross_p$pheno <- replace_zero_with_NA(in_cross_p$pheno, minimum_samples_p)
in_cross_y$pheno <- replace_zero_with_NA(in_cross_y$pheno, minimum_samples_y)


# run R/qtl2 --------------------------------------------------------------
# pr <- calc_genoprob(in_cross, in_cross$gmap, error_prob = 1e-4, cores = max_treads)

# Insert pseudomarkers for the gmap used in each analysis
map_r <- insert_pseudomarkers(in_cross_r$gmap, step = 1, cores = 1)
map_p <- insert_pseudomarkers(in_cross_p$gmap, step = 1, cores = 1)
map_y <- insert_pseudomarkers(in_cross_y$gmap, step = 1, cores = 1)

# make metadata -----------------------------------------------------------

p_m_n <- names(unlist(map_r, use.names=TRUE, recursive = TRUE))

p_m_n_chr <- NULL
p_m_n_marker_id <- NULL
for (i in 1:length(p_m_n)) {
  p_m_n_chr[as.numeric(i)] <- strsplit(p_m_n[as.numeric(i)], "[.]")[[1]][1]
 # p_m_n_marker_id[i] <- strsplit(p_m_n[i], "[.]")[[1]][2]
  p_m_n_marker_id[i] <- gsub("^.*\\.","", p_m_n[as.numeric(i)])
}

new_gmap <- as.data.frame(cbind(p_m_n_marker_id, p_m_n_chr, unlist(map_r, use.names=TRUE, recursive = TRUE)))
colnames(new_gmap) <- c("marker_id","chr","pos")
new_gmap$marker_id <- as.factor(new_gmap$marker_id)
new_gmap$chr <- as.integer(new_gmap$chr)
new_gmap$pos <- as.numeric(new_gmap$pos)

chr_metadata <- new_gmap %>%
  arrange(chr, pos) %>%
  rowid_to_column("index") %>%
  group_by(chr) %>%
  summarise(num_markers = n(), start = min(index), end = max(index), width = end - start) %>%
  remove_missing() %>%
  arrange(chr) %>%
  mutate(chr = factor(chr, ordered = T))

############

# without pseudomarkers (Pseudomarkers are used in final analysis!)
# pr_pdw_r <- calc_genoprob(in_cross_pdw_r, error_prob = 1e-4, cores = max_treads)
# pr_n <- calc_genoprob(in_cross_r_new, map_n, error_prob = 1e-4, cores = max_treads)

# With Pseudomarkers as done in the final analysis
pr_r <- calc_genoprob(in_cross_r, map_r, error_prob = 1e-4, cores = max_treads)
pr_p <- calc_genoprob(in_cross_p, map_p, error_prob = 1e-4, cores = max_treads)
pr_y <- calc_genoprob(in_cross_y, map_y, error_prob = 1e-4, cores = max_treads)

# grid <- calc_grid(in_cross$gmap, step=1)
# pr_grid <- probs_to_grid(pr, grid)
# kinship_grid <- calc_kinship(pr_grid)
# apr <- genoprob_to_alleleprob(pr)

cat("Calculating a kinship matrix...\n")
# kinship_n <- calc_kinship(pr_n)

kinship_r <- calc_kinship(pr_r)

# herit_n <- est_herit(in_cross_r_new$pheno, kinship_n, addcovar=metadata_r)

###### Estimating the heritability of 16S abundances using Replicate, Harvest Day, PDW, Leaves, Rhizosphere, and Soil_Before
herit_r <- est_herit(in_cross_r$pheno, kinship_r, addcovar=metadata_r)
herit_p <- est_herit(in_cross_p$pheno, kinship_r, addcovar=metadata_r[,-1])
herit_y <- est_herit(in_cross_y$pheno, kinship_r, addcovar=metadata_r[,-1])

# It is possible to filter based on heritability, however here we do not.
herit_thres <- 0.0
# l_herit_n <- length(which(herit_n>=herit_thres))
l_herit_r <- length(which(herit_r>=herit_thres))
l_herit_p <- length(which(herit_p>=herit_thres))
l_herit_y <- length(which(herit_y>=herit_thres))

# out_n <- scan1(pr_n, in_cross_r_new$pheno, addcovar = metadata_r)
out_r <- scan1(pr_r, in_cross_r$pheno[,which(herit_r>=herit_thres)], addcovar = metadata_r)
out_p <- scan1(pr_p, in_cross_p$pheno[,which(herit_p>=herit_thres)], addcovar = metadata_r[,-1])
out_y <- scan1(pr_y, in_cross_y$pheno[,which(herit_y>=herit_thres)], addcovar = metadata_r[,-1])


# operm_n <- scan1perm(pr_n, in_cross_r_new$pheno, addcovar = metadata_r, n_perm=10)
operm_r <- scan1perm(pr_r, in_cross_r$pheno[,which(herit_r>=herit_thres)], addcovar = metadata_r, n_perm = 1000, cores = 12)
operm_p <- scan1perm(pr_p, in_cross_p$pheno[,which(herit_p>=herit_thres)], addcovar = metadata_r[,-1], n_perm = 1000, cores = 12)
operm_y <- scan1perm(pr_y, in_cross_y$pheno[,which(herit_y>=herit_thres)], addcovar = metadata_r[,-1], n_perm = 1000, cores = 12)

# qlod_sum_n <- quantile(operm_n, c(0.8,0.85,0.9,0.95))
qlod_sum_r <- quantile(operm_r, c(0.8,0.85,0.9,0.95))
qlod_sum_p <- quantile(operm_p, c(0.8,0.85,0.9,0.95))
qlod_sum_y <- quantile(operm_y, c(0.8,0.85,0.9,0.95))
#
# WRITE OUT permutation scores
# write.table(t(operm), file = output_plodscores, sep = "\t", quote = F, row.names = T, col.names = F)

##########


# IMPORTANT, using the permutation threshold from the combined global analysis

# lod_inclusion_threshold_n <- as.numeric(qlod_sum_n[1])
lod_inclusion_threshold_r <- as.numeric(qlod_sum_r[1])
lod_inclusion_threshold_p <- as.numeric(qlod_sum_p[1])
lod_inclusion_threshold_y <- as.numeric(qlod_sum_y[1])

cat("Finding peaks...\n")

# found_peaks_n <- find_peaks(out_n, map = map_n, threshold = lod_inclusion_threshold_n) %>

# Both replicates combined
found_peaks_r <- find_peaks(out_r, map = map_r, threshold = lod_inclusion_threshold_r, prob = 0.95) %>%
  mutate(asv = as.numeric(substring(lodcolumn, 4))) %>%
  left_join(asv_taxonomy, by = "asv")

# Only one set of replicates included in the shotgun metagenomics
found_peaks_p <- find_peaks(out_r, map = map_r, threshold = lod_inclusion_threshold_r, prob = 0.95) %>%
  mutate(asv = as.numeric(substring(lodcolumn, 4))) %>%
  left_join(asv_taxonomy, by = "asv")


Add_herit_and_effects <- function(found_peaks_X, pr_X, in_cross_X, metadata_X, map_X, herit_X) {
  for (i in 1:nrow(found_peaks_X)){
    contig_id <- found_peaks_X$lodcolumn[i]
    chr <- as.character(found_peaks_X$chr[i])
    pos <- found_peaks_X$pos[i]
    coefs <- scan1coef(pr_X[, chr], in_cross_X$pheno[,contig_id], addcovar = metadata_X)
    coef <- coefs[which(rownames(coefs) == find_marker(map_X, chr = chr, pos = pos)),][1]
    found_peaks_X$effect[i] <-  coef / mean(in_cross_X$pheno[,contig_id],na.rm=TRUE)
    found_peaks_X$heritability[i] <- herit_X[contig_id]
    #print(coef[1])
    if(coef[1]>0){
      print(paste(contig_id,"AA, modern"))
      found_peaks_X$allele[i] <- "modern"
    }
    else{
      print(paste(contig_id, "BB, wild"))
      found_peaks_X$allele[i] <- "wild"
    }
  }
 return(found_peaks_X)
}

found_peaks_r <- Add_herit_and_effects(found_peaks_r, pr_r, in_cross_r, metadata_r, map_r, herit_r)
found_peaks_p <- Add_herit_and_effects(found_peaks_p, pr_p, in_cross_p, metadata_r, map_p, herit_p)

# TO DO
# cat("Writing QTL peaks to", output_peaks, "...\n")
write.table(found_peaks_r, file = "r_peaks.txt", sep = "\t", quote = F, row.names = F, col.names = T)
write.table(found_peaks_p, file = "p_peaks.txt", sep = "\t", quote = F, row.names = F, col.names = T)

# plotting example ----------------------------------------------------------------

ymx2 <- max(in_cross_r$pheno,na.rm=TRUE)
column_id <- 2
plot(out_r, map_r, lodcolumn=column_id, col=rgb(0,0,1,0.5), ylim = c(0,ymx+0.2), main = colnames(out_r)[column_id])
abline(h= c(lod_inclusion_threshold_r), col=c(rgb(0,0,1,0.5)))
