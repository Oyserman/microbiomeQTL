#!/usr/bin/env Rscript

# VERSION FOR RUNNING ON THE SERVER

# setwd("D:/Documents/WUR/MSc Thesis MBF/Analyses/QTL_ECvdW")
setwd("/lustre/BIF/nobackup/wijk059/working_folder_elmar/R_QTL_on_server")
.libPaths(c("./R-libraries", .libPaths()))

if (!file.exists("./QTL_analysis_output")){dir.create(file.path("./QTL_analysis_output"))}
# program setup -----------------------------------------------------------
# install.packages("tidyverse", repos = "http://cran.us.r-project.org")
# install.packages("ggpubr", repos = "http://cran.us.r-project.org")
# install.packages("RColorBrewer", repos = "http://cran.us.r-project.org")
# install.packages("qtl2", repos = "http://cran.us.r-project.org")
# install.packages("cowplot", repos = "http://cran.us.r-project.org")
# install.packages("yaml", repos = "http://cran.us.r-project.org")
# install.packages("multcompView", repos = "http://cran.us.r-project.org")

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggpubr)
  library(RColorBrewer)
  library(qtl2)
  library(cowplot)
  library(yaml)
  library(multcompView)
})

argv <- commandArgs(trailingOnly = T)
theme_set(theme_pubr(border = T, legend = "right"))
max_treads <- 24

# Usage:
# Rscript asv_qtl_April2020 pink_flexiblecore.csv

# load R/qtl2 data --------------------------------------------------------

temp_snp_count_file_name <- "./all_selected_snps_MM_P_bed_ECvdW_0based_bed_210617_extracted_snp_allele_counts.tsv"

temp_snp_count_file = read.csv(file=temp_snp_count_file_name, header = TRUE, sep = "\t")
temp_snp_count_file$id = paste(temp_snp_count_file$contig, "pos", temp_snp_count_file$position, temp_snp_count_file$allele, sep="_")

# create correct naming format and filter out files columns that are not needed
# change names: so that the .1 of the metadata matches the p... that has no .1 but does have .2
temp_sample_names = colnames(temp_snp_count_file)

sample_names_to_have_point_1_included = list()

for (i in 1:length(temp_sample_names)){
  sample_name = temp_sample_names[i]
  if (startsWith(sample_name, "X")){
    sample_name = str_split(sample_name, "_")[[1]][1]
    sample_name = str_split(sample_name, "X")[[1]][2]
    sample_name = paste("p", sample_name, sep = "")
    if (endsWith(sample_name, ".2")){
      sample_names_to_have_point_1_included = c(sample_names_to_have_point_1_included, str_split(sample_name, "\\.")[[1]][1])
    }
    temp_sample_names[i] = sample_name
  }
}

temp_sample_names[which(temp_sample_names %in% sample_names_to_have_point_1_included)] = paste(temp_sample_names[which(temp_sample_names %in% sample_names_to_have_point_1_included)], ".1", sep = "")

colnames(temp_snp_count_file) = temp_sample_names

snp_count_file = subset(temp_snp_count_file, select = c(which(temp_sample_names == "id"),which(startsWith(temp_sample_names[4:length(temp_sample_names)], "p"))+3))

write.table(t(snp_count_file), file="ril_snp_allele_count.csv", sep= ",", quote = F, row.names = T, col.names = F)

rqtl_phenofile = "ril_snp_allele_count.csv"

# Load the yaml files for the microQTL and for the PDW
rqtl_rnaseq_map <- "rnaseq_map.yaml"



# rqtl_rnaseq_map_pdw <- "./rnaseq_map_pdw.yaml"


# taxonomy_file <- "result_April2020/asv_taxonomy.rds"
# asv_taxonomy <- readRDS(taxonomy_file)

# Treat 0's as NA if 75% or more of samples are non-zero

############ Classic Phenotypes First #############

zeros <- TRUE

# Curate the metadata

metadata_file <- "RIL_metadata_replicates_5_19_2019.txt"

metadata <- read.table(file=metadata_file, row.names = 1, header = TRUE)[,-3]

metadata = metadata[which(metadata$Color == "P"),] # only keep P as we do not use Y
metadata = subset(metadata, select = -(Color))


######### Defining outputs names
#########
output_basename <- tools::file_path_sans_ext(rqtl_phenofile)
output_lodscores <- sprintf("%s_lod.csv", output_basename)
output_plodscores <- sprintf("%s_plod.csv", output_basename)
output_peaks <- sprintf("%s_peaks.csv", output_basename)
output_heritability <- sprintf("%s_heritability.csv", output_basename)
output_heatmap <-sprintf("%s_lod.pdf", output_basename)
#########

# Run rQTL2
# pmap_1 <- read_csv("Tomaat_Micro_biome_RNA_Seq_map/pmap.csv")

# Now the micro phenotypes
in_cross <- read_cross2(rqtl_rnaseq_map, quiet = F)

# read the yaml files
rqtl_yaml <- read_yaml(rqtl_rnaseq_map)


# run R/qtl2 --------------------------------------------------------------
# pr <- calc_genoprob(in_cross, in_cross$gmap, error_prob = 1e-4, cores = max_treads)

# Insert pseudomarkers for the gmap used in each analysis
#####
#
map <- insert_pseudomarkers(in_cross$gmap, step = 1, cores = 1)

# make metadata -----------------------------------------------------------
# SKIP
#####

# p_m_n <- names(unlist(map, use.names=TRUE, recursive = TRUE))
#
# p_m_n_chr <- NULL
# p_m_n_marker_id <- NULL
# for (i in 1:length(p_m_n)) {
#   p_m_n_chr[as.numeric(i)] <- strsplit(p_m_n[as.numeric(i)], "[.]")[[1]][1]
#  # p_m_n_marker_id[i] <- strsplit(p_m_n[i], "[.]")[[1]][2]
#   p_m_n_marker_id[i] <- gsub("^.*\\.","", p_m_n[as.numeric(i)])
# }
#
# new_gmap <- as.data.frame(cbind(p_m_n_marker_id, p_m_n_chr, unlist(map_n, use.names=TRUE, recursive = TRUE)))
# colnames(new_gmap) <- c("marker_id","chr","pos")
# new_gmap$marker_id <- as.factor(new_gmap$marker_id)
# new_gmap$chr <- as.integer(new_gmap$chr)
# new_gmap$pos <- as.numeric(new_gmap$pos)
#
# chr_metadata <- new_gmap %>%
#   arrange(chr, pos) %>%
#   rowid_to_column("index") %>%
#   group_by(chr) %>%
#   summarise(num_markers = n(), start = min(index), end = max(index), width = end - start) %>%
#   remove_missing() %>%
#   arrange(chr) %>%
#   mutate(chr = factor(chr, ordered = T))

############

# ECvdW: use this:
# without pseudomarkers (Pseudomarkers are used in final analysis)
# pr_pdw_r <- calc_genoprob(in_cross_pdw_r, error_prob = 1e-4, cores = max_treads)
# pr_n <- calc_genoprob(in_cross_r_new, map_n, error_prob = 1e-4, cores = max_treads)

pr <- calc_genoprob(in_cross, map, error_prob = 1e-4, cores = max_treads) # pseudomarkers
# pr <- calc_genoprob(in_cross, error_prob = 1e-4, cores = max_treads) # no pseudomarkers


# grid <- calc_grid(in_cross$gmap, step=1)
# pr_grid <- probs_to_grid(pr, grid)
# kinship_grid <- calc_kinship(pr_grid)
# apr <- genoprob_to_alleleprob(pr)

cat("Calculating a kinship matrix...\n")
# kinship_n <- calc_kinship(pr_n)

kinship <- calc_kinship(pr)

# herit_n <- est_herit(in_cross_r_new$pheno, kinship_n, addcovar=metadata_r)

###### Estimating the heritability of 16S abundances using Replicate, Harvest Day, PDW, Leaves, Rhizosphere, and Soil_Before
###### FOr the pink and yellow, as above but without the replicate
herit <- est_herit(in_cross$pheno, kinship, addcovar=metadata)
length(herit)
plot(herit)

# It is possible to filter based on heritability, however here we do not.
herit_thres <- 0.0
# l_herit_n <- length(which(herit_n>=herit_thres))
l_herit <- length(which(herit>=herit_thres))


# TODO
# Write out a single table with all heritability scores
herit_summary = summary(herit)
write.table(herit, file = "QTL_analysis_output/Heritability_210618.txt", sep = "\t", quote = F, row.names = T, col.names = T)

# out <- scan1(pr, in_cross$pheno[,which(herit>herit_thres)], kinship = kinship, model = "normal")


#####
# out_n <- scan1(pr_n, in_cross_r_new$pheno, addcovar = metadata_r)
out <- scan1(pr, in_cross$pheno[,which(herit>=herit_thres)], addcovar = metadata, cores = max_treads)

# operm_n <- scan1perm(pr_n, in_cross_r_new$pheno, addcovar = metadata_r, n_perm=10)
operm <- scan1perm(pr, in_cross$pheno[,which(herit>=herit_thres)], addcovar = metadata, n_perm = 1000, cores = max_treads) # change later to n_perm = 1000

# qlod_sum_n <- quantile(operm_n, c(0.8,0.85,0.9,0.95))
qlod_sum <- quantile(operm, c(0.8,0.85,0.9,0.95))

save.image(file = "global_environment_QTL_210618.RData")

lod_scores <- as.data.frame(out) %>%
  rownames_to_column("marker") %>%
  rowid_to_column("index") %>%
  separate(marker, into = c("chr", "pos"), remove = F, sep = "([\\.])") %>%
  pivot_longer(-c(1:4), names_to = "phenotype", values_to = "lod_score", names_ptypes = list(phenotype = factor()))
# separate(phenotype, into = c("contig", "temp_pos"), remove = F, sep = "(_pos_)") %>%
# separate(temp_pos, into = c("m_pos", "allele"), remove = T, sep = "(_)")

lod_scores[!(startsWith(lod_scores$chr, "c")),c("chr", "pos")] <- NA
dim(lod_scores)
dim(lod_scores[which(!(is.na(lod_scores$pos))),])
lod_scores[which(!(is.na(lod_scores$pos))),]

# lod_scores[!(lod_scores$chr %in% 1:12), c("pos", "chr")] <- NA
#
# test_out = as.data.frame(out[1:500,])
# dim(test_out)
# test_out = rownames_to_column(test_out, "marker")
# test_out = rowid_to_column(test_out, "index")
# test_out = separate(test_out, marker, into = c("chr", "pos"), remove = F, sep = "([\\.])")
# test_out = pivot_longer(test_out, -c(1:4), names_to = "phenotype", values_to = "lod_score", names_ptypes = list(phenotype = factor()))
# test_out = separate(test_out, phenotype, into = c("contig", "m_pos"), remove = F, sep = "(_pos_)")
# test_out = separate(test_out, m_pos, into = c("m_pos", "allele"), remove = T, sep = "(_)")
# dim(test_out)
# test_out[!(startsWith(test_out$chr, "c")),c("chr", "pos")] <- NA
# View(test_out[1:100,])
# test_out[which(!(is.na(test_out$chr))),]


output_lodscores = "QTL_analysis_output/output_lodscores_210618.tsv"
cat("Writing LOD scores to", output_lodscores, "...\n")

#########
# plot(out_y,map_y)

write.table(lod_scores, file = output_lodscores, sep = "\t", quote = F, row.names = F, col.names = T)

# infer LOD peak threshold using permutation testing
# out_perm <- scan1perm(pr, in_cross$pheno, n_perm = 1000, cores = max_treads)
# perm_thresh <- as.data.frame(summary(out_perm, alpha = 0.05)) %>%
#   gather("phenotype", "threshold")


# IMPORTANT, only using the permutation threshold from the combined global analysis

# lod_inclusion_threshold_n <- as.numeric(qlod_sum_n[1])
lod_inclusion_threshold <- as.numeric(qlod_sum[4]) # [1] is the 80% quantile, c(0.8,0.85,0.9,0.95) is the list, so [4] is the 95%

cat("Finding peaks...\n")


# found_peaks_n <- find_peaks(out_n, map = map_n, threshold = lod_inclusion_threshold_n) %>%
#   mutate(asv = as.numeric(substring(lodcolumn, 4))) %>%
#   left_join(asv_taxonomy, by = "asv")

# found_peaks <- find_peaks(out, map = map, threshold = lod_inclusion_threshold, prob = 0.95) %>%
#   mutate(asv = as.numeric(substring(lodcolumn, 4))) %>%
#   left_join(asv_taxonomy, by = "asv")
lod_inclusion_threshold
found_peaks <- find_peaks(out, map = map, threshold = lod_inclusion_threshold, prob = 0.95)
dim(found_peaks)

for (i in 1:nrow(found_peaks)){
  # i=1 ###TEST
  snp_allele_pos <- found_peaks$lodcolumn[i]
  chr <- as.character(found_peaks$chr[i])
  pos <- found_peaks$pos[i]
  coefs <- scan1coef(pr[,chr], in_cross$pheno[,snp_allele_pos], addcovar = metadata)
  coef <- coefs[which(rownames(coefs) == find_marker(map, chr = chr, pos = pos)),][1]
  found_peaks$effect[i] <-  coef / mean(in_cross$pheno[,snp_allele_pos],na.rm=TRUE)
  found_peaks$heritability[i] <- herit[snp_allele_pos]
  
  # print(coef[1])
  
  
  if(coef[1]>0){
    # print(paste(contig_id,"AA", "modern"))
    found_peaks$allele[i] <- "modern"
  }
  else{
    # print(paste(contig_id, "BB", "wild"))
    found_peaks$allele[i] <- "wild"
  }
}
dim(found_peaks)


write.table(found_peaks, file = "QTL_analysis_output/found_peaks_210618.txt", sep = "\t", quote = F, row.names = F, col.names = T)

# from here on, continue the rest of the analysis
