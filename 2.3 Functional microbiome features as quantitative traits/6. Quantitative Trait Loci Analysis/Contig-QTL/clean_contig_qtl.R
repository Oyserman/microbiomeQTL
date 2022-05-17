#!/usr/bin/env Rscript

# program setup -----------------------------------------------------------
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
max_treads <- 4

# Usage:
#

# load R/qtl2 data --------------------------------------------------------

setwd("data/final_assembly/qtl/metagenome_qtl/taxo_level/20200601/")
rqtl_phenofile_p <- "phenofile_filtered_rhi_norm_log_count_p.csv"

pfile <- read.csv(rqtl_phenofile_p)

rqtl_control_pink <- "rnaseq_map_control.yaml"
contig_taxonomy <- read.delim("~/data/final_assembly/kraken/2020_may/contig_10kb_taxonomy.tsv",stringsAsFactors=FALSE)

#old_map_geno <- read.csv("../20200518/sl23_geno.csv")
rnaseq_map_geno <- read.csv("Tomaat_Micro_biome_RNA_Seq_map/geno.csv")
############ Import metadata #############
metadata_p <- read.table("pRIL_metadata_5_19_2019_formated.txt", sep = "\t",row.names = 1, header = T)

############ Check 0 in the counts #############
#filtered_rhi_count <- read.csv("phenofile_filtered_rhi_norm_log_count_p.csv")
#sum(colSums(filtered_rhi_count==0, na.rm = T) >=1)
#plot(colSums(filtered_rhi_count==0, na.rm = T)/nrow(filtered_rhi_count)*100, ylim = c(0, 0.5))


######### Read QTL data from a set of files listed in the control.yaml #############
in_cross_p <- read_cross2(rqtl_control_pink, quiet = F)
check_c <- check_cross2(in_cross_p)

library(qtl)
plotMap(in_cross_p$gmap)
plot.cross(in_cross_p)


#rqtl_yaml_p <- read_yaml(rqtl_control_pink)

#minimum_samples_p <- dim(in_cross_p$pheno)[1]
#in_cross_p$pheno <- replace_zero_with_NA(in_cross_p$pheno, minimum_samples_p)


######### Run R/qtl2 #############
# pr <- calc_genoprob(in_cross_p, in_cross$gmap, error_prob = 1e-4, cores = max_treads)

#### If pseudomarkers are used:
map_p <- insert_pseudomarkers(in_cross_p$gmap, step = 1, cores = 1)

#### make metadata of the new gmap with pseudomarkers used -----------------------------------------------

#p_m_n <- names(unlist(map_p, use.names=TRUE, recursive = TRUE))

# p_m_n_chr <- NULL
# p_m_n_marker_id <- NULL
# for (i in 1:length(p_m_n)) {
#   p_m_n_chr[as.numeric(i)] <- strsplit(p_m_n[as.numeric(i)], "[.]")[[1]][1]
#   # p_m_n_marker_id[i] <- strsplit(p_m_n[i], "[.]")[[1]][2]
#   p_m_n_marker_id[i] <- gsub("^.*\\.","", p_m_n[as.numeric(i)])
# }
# 
# new_gmap <- as.data.frame(cbind(p_m_n_marker_id, p_m_n_chr, unlist(map_p, use.names=TRUE, recursive = TRUE)))
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
#### -----------------------------------------------

set.seed(20200518)
#### Calculate conditional genotype probabilities-----------------------------------------------
pr_p <- calc_genoprob(in_cross_p, map_p, error_prob = 1e-4, cores = max_treads)
#pr_p_no_pseudo <- calc_genoprob(in_cross_p, error_prob = 1e-4, cores = max_treads)

# grid <- calc_grid(in_cross$gmap, step=1)
# pr_grid <- probs_to_grid(pr, grid)
# kinship_grid <- calc_kinship(pr_grid)

# apr <- genoprob_to_alleleprob(pr)

# cat("Calculating a kinship matrix...\n")

#### Calculate a kinship matrix for heritability estimation-----------------------------------------------
kinship_p <- calc_kinship(pr_p)
#kinship_p_no_pseudo <- calc_kinship(pr_p_no_pseudo)

herit_p <- est_herit(in_cross_p$pheno, kinship_p, addcovar=metadata_p)
#herit_p_no_psudo <- est_herit(in_cross_p$pheno, kinship_p_no_pseudo, addcovar=metadata_p)


plot(herit_p)
library(ggplot2)
rhiz_contig <- colnames(in_cross_p$pheno)
rhiz_taxonomy <- contig_taxonomy[contig_taxonomy$contig_id %in% rhiz_contig,]
df_heritability <- data.frame(herit_p)
df_heritability$contig_id <- rownames(df_heritability)
df_heri_taxo <- left_join(df_heritability, rhiz_taxonomy, by = "contig_id")
ggplot(df_heri_taxo, aes(x=contig_id, y=herit_p, color=class)) + 
  geom_point() +
  labs(x = "Contig",y = "Heritability",colour = "Class",
    title = "Heritability of rhizosphere enriched contigs")+
  geom_abline(intercept = 0.2, slope = 0, color="black", linetype="dashed", size=0.8)+
  theme_light()+ 
  theme(
    plot.title = element_text(face = "bold", size = 12),
    axis.ticks.y = element_line(colour = "grey70", size = 0.2),
    axis.ticks.x = element_blank())
#plot(herit_p_no_psudo)
#l_herit_p_no_pseudo <- length(which(herit_p_no_psudo>=herit_thres))

#abline(h=0.2, lty=2, col="red")
herit_thres <- 0.2
l_herit_p <- length(which(herit_p>=herit_thres))
)
# TO DO #
# Write out a single table with all heritability scores
#df_herit_p<- data.frame(herit_p)
#df_herit_p$contig_id <- rownames(df_herit_p)
#write.table(df_herit_p, file = "results/Heritability_rhi_contig.txt", sep = "\t", quote = F, row.names = F, col.names = T)

#### Perform a genome scan-----------------------------------------------
# out <- scan1(pr, in_cross$pheno[,which(herit>herit_thres)], kinship = kinship, model = "normal")
out_p <- scan1(pr_p, in_cross_p$pheno, addcovar = metadata_p)
out_p_heri <- scan1(pr_p, in_cross_p$pheno[,which(herit_p>=herit_thres)], addcovar = metadata_p)
#out_p_no_pseudo <- scan1(pr_p_no_pseudo, in_cross_p$pheno[,which(herit_p_no_psudo>=herit_thres)], addcovar = metadata_p)

#plot
#par(mar=c(5.1, 4.1, 1.1, 1.1))
#ymx <- maxlod(out_p) # overall maximum LOD score
#plot(out_p, map_p, lodcolumn=1, col="slateblue", ylim=c(0, ymx*1.02))
#plot(out_p, map_p, lodcolumn=2, col="violetred", add=TRUE)
#legend("topleft", lwd=2, col=c("slateblue", "violetred"), colnames(out_p), bg="gray90")

#### Run permutation test-----------------------------------------------
operm_p <- scan1perm(pr_p, in_cross_p$pheno, addcovar = metadata_p, n_perm=1000,cores = 8)
operm_p_heri <- scan1perm(pr_p, in_cross_p$pheno[,which(herit_p>=herit_thres)], addcovar = metadata_p, n_perm=1000,cores = 20)

#operm_p_pseudo <- scan1perm(pr_p_no_pseudo, in_cross_p$pheno[,which(herit_p_no_psudo>=herit_thres)], addcovar = metadata_p, n_perm=10)

qlod_sum_p <- quantile(operm_p, c(0.8,0.85,0.9,0.95))
qlod_sum_p_heri <- quantile(operm_p_heri, c(0.8,0.85,0.9,0.95))

#qlod_sum_p_pseudo <- quantile(operm_p_pseudo, c(0.8,0.85,0.9,0.95))


# TO DO #
#! WRITE OUT pernutation scores
#! write.table(t(operm), file = output_plodscores, sep = "\t", quote = F, row.names = T, col.names = F)

# lod_scores_p_test <- as.data.frame(out_p) %>%
#   rownames_to_column("marker") %>%
#   rowid_to_column("index") %>%
#   separate(marker, into = c("pos", "chr"), remove = F, sep = "([\\.\\-])") %>%
#   pivot_longer(-c(1:4), names_to = "phenotype", values_to = "lod_score", names_ptypes = list(phenotype = factor()))
# 
# lod_scores_p_test[!(lod_scores_p_test$chr %in% 1:12), c("pos", "chr")] <- NA

# cat("Writing LOD scores to", output_lodscores, "...\n")


plot(out_p,map_p)
#plot(out_p_no_pseudo,in_cross_p$gmap)


# write.table(lod_scores, file = output_lodscores, sep = "\t", quote = F, row.names = F, col.names = T)


# IMPORTANT, only using the permutation threshold from the combined global analysis

## use 0.95 quantile threshold
lod_inclusion_threshold_p <- as.numeric(qlod_sum_p[4])
lod_inclusion_threshold_p_heri <- as.numeric(qlod_sum_p_heri[4])


#### Find lod peaks-----------------------------------------------
found_peaks_p <- find_peaks(out_p, map = map_p, threshold = lod_inclusion_threshold_p,prob = 0.95) %>%
  mutate(contig_id = lodcolumn) %>%
  left_join(contig_taxonomy, by = "contig_id")

found_peaks_p_heri <- find_peaks(out_p_heri, map = map_p, threshold = lod_inclusion_threshold_p_heri,prob = 0.95) %>%
  mutate(contig_id = lodcolumn) %>%
  left_join(contig_taxonomy, by = "contig_id")

found_peaks_p$replicate <- rep("pink",717)

# found_peaks_p_no_pseudo <- find_peaks(out_p_no_pseudo, map =in_cross_p$gmap, threshold = lod_inclusion_threshold_p) %>%
#   mutate(contig_id = lodcolumn) %>%
#   left_join(contig_taxonomy, by = "contig_id")

write.table(found_peaks_p,"pertumation1000/found_peak_no_heritability_filter.tsv", quote = F, col.names = T, row.names = F, sep = "\t")
#write.table(found_peaks_p,"results/found_peak_CI.tsv", quote = F, col.names = T, row.names = F, sep = "\t")

#found_peaks_p_old_map <- read.delim("../20200518/results/lod_peaks.tsv")

length(unique(found_peaks_p$lodcolumn))
length(unique(found_peaks_p_heri$lodcolumn))

#length(unique(found_peaks_p_old_map$lodcolumn))
#sum(unique(found_peaks_p$lodcolumn) %in% unique(found_peaks_p_old_map$lodcolumn))

#plot(found_peaks_p_old_map$chr)
table(as.numeric(found_peaks_p$chr))


# TO DO
# cat("Writing QTL peaks to", output_peaks, "...\n")
# write.table(found_peaks, file = output_peaks, sep = "\t", quote = F, row.names = F, col.names = T)

# plotting ----------------------------------------------------------------

uniqe_phyla_with_peaks <- as.character(unique(found_peaks_p$phylum))
# [1] "Actinobacteria"              "Proteobacteria"              "Candidatus_Saccharibacteria"
# [4] "unassigned" 
uniqe_classes_with_peaks <- as.character(unique(found_peaks_p$class))
# [1] "Actinobacteria"      "Gammaproteobacteria" "Alphaproteobacteria"

uniqe_orders_with_peaks <- as.character(unique(found_peaks_p$order))
# [1] "Streptomycetales" "Cellvibrionales"  "Sphingomonadales" "Rhizobiales"      "Xanthomonadales" 
# [6] "unassigned"       "Pseudomonadales"  "Caulobacterales" 

# Within the Actinobacteria
which(found_peaks_p$class==uniqe_classes_with_peaks[1])
# Within the Actinobacteria,309 Streptomyces are represented
table(as.character(found_peaks_p$genus[which(found_peaks_p$class==uniqe_classes_with_peaks[1])]))
# Streptomyces 
# 649
length(unique(found_peaks_p$contig_id[found_peaks_p$genus=="Streptomyces"]))

# Within the Alphaproteobacteria
which(found_peaks_p$class==uniqe_classes_with_peaks[2])
# Within the Gammaproteobacteria,309 Streptomyces are represented
table(as.character(found_peaks_p$genus[which(found_peaks_p$class==uniqe_classes_with_peaks[2])]))
# Altererythrobacter    Novosphingobium       Sphingomonas       Sphingopyxis 
# 1                  1                  1                 81 

# Within the Gammaproteobacteria
which(found_peaks_p$class==uniqe_classes_with_peaks[3])

# Within the Gammaproteobacteria,309 Streptomyces are represented
table(as.character(found_peaks_p$genus[which(found_peaks_p$class==uniqe_classes_with_peaks[3])]))
# Cellvibrio 
# 151 
###### 20200610 END #######

write.table(unique(found_peaks_p$contig_id),"results/peak_contig_id.txt", quote = F, col.names = F, row.names = F, sep = "\n")
write.table(unique(found_peaks_p$contig_id[found_peaks_p$genus=="Streptomyces"]),"results/strep_peak_contig_id.txt", quote = F, col.names = F, row.names = F, sep = "\n")


##### plot

# plot LOD curves for all unique phenotypes (contigs)
unique_contig_peak <- unique(found_peaks_p$lodcolumn)
table(found_peaks_p$genus)

pdf("all_peaks.pdf", width = 10.2, height = 4.3) 
par(mar=c(4,4,4,1))
plot(out_p, map_p, lodcolumn=unique_contig_peak[1], col=rgb(0,0,0,0.01), ylim = c(0,ymx+0.2),main="Overview of metageome QTL peaks")
for (i in unique_contig_peak[2:476]){
  plot(out_p, map_p, lodcolumn=i, col=rgb(0,0,0,0.01), add=TRUE)}
abline(h= lod_inclusion_threshold_p , col=c(rgb(0,0,0,0.5)), lty = 2)
dev.off()

# plot LOD curves for streptomtces contigs
unique_strep_contig_peak <- unique(found_peaks_p[found_peaks_p$genus=="Streptomyces","lodcolumn"])
pdf("streptomyces_peaks.pdf", width = 10.2, height = 4.3) 
plot(out_p, map_p, lodcolumn=unique_strep_contig_peak[1], col=rgb(0,0,1,0.005), ylim = c(0,ymx+0.2), main="Streptomyces")
for (i in unique_strep_contig_peak[2:300]){
  plot(out_p, map_p, lodcolumn=i, col=rgb(0,0,1,0.005), add=TRUE)}
abline(h= lod_inclusion_threshold_p , col=c(rgb(0,0,1,0.5)), lty = 2)
dev.off()

# plot LOD curves for cellvibrio contigs
pdf("Cellvibrio_peaks.pdf", width = 10.2, height = 4.3) 
unique_cellv_contig_peak <- unique(found_peaks_p[found_peaks_p$genus=="Cellvibrio","lodcolumn"])
plot(out_p, map_p, lodcolumn=unique_cellv_contig_peak[1], col=rgb(1,0,0,0.005), ylim = c(0,ymx+0.2), main="Cellvibrio")
for (i in unique_cellv_contig_peak[2:105]){
  plot(out_p, map_p, lodcolumn=i, col=rgb(1,0,0,0.005), add=TRUE)}
abline(h= lod_inclusion_threshold_p , col=c(rgb(1,0,0,0.5)), lty = 2)
dev.off()

pdf("colored_all_peaks.pdf", width = 10.2, height = 4.3) 
unique_genus <- unique(found_peaks_p$genus)
cols<-rainbow(length(unique_genus),alpha = 0.8)
names(cols) <- unique_genus
cols[1:3] <- adjustcolor(cols[1:3], alpha.f = 0.1)
par(mar=c(5, 4, 4, 6)+0.1)
plot(out_p, map_p, lodcolumn=unique_contig_peak[1],ylim = c(0,maxlod(out_p)+0.2), col=cols[1],main="Overview of metageome QTL peaks")
for (i in 2:length(unique_contig_peak)){
  contig_id <- unique_contig_peak[i]
  genus <- found_peaks_p$genus[found_peaks_p$lodcolumn==contig_id]
  color <- cols[genus]
  plot(out_p, map_p, lodcolumn=contig_id, col=color, add=TRUE)}
abline(h= lod_inclusion_threshold_p , col=c(rgb(0,0,0,0.5)), lty = 2)
op <- par(cex=0.6)
legend("right", inset = c(-0.15,0), legend = unique_genus, xpd = TRUE, 
       horiz = FALSE, col = rainbow(12), lty = 1, bty = "n")
par(op)
dev.off()

for (i in 1:nrow(found_peaks_p)){
  contig_id <- found_peaks_p$lodcolumn[i]
  chr <- as.character(found_peaks_p$chr[i])
  pos <- found_peaks_p$pos[i]
  coefs <- scan1coef(pr_p[, chr], in_cross_p$pheno[,contig_id], addcovar = metadata_p)
  coef <- coefs[which(rownames(coefs) == find_marker(map_p, chr = chr, pos = pos)),][1]
  found_peaks_p$effect[i] <-  coef/ mean(in_cross_p$pheno[,contig_id],na.rm=TRUE)
  found_peaks_p$heritability[i] <- herit_p[contig_id]
  #print(coef[1])
  #if(coef[1]>0) print("Non-negative number") else print("Negative number")
  if(coef[1]>0){
    print(paste(contig_id,"AA, modern"))
    found_peaks_p$allele[i] <- "modern"
  }
  else{
    print(paste(contig_id, "BB, wild"))
    found_peaks_p$allele[i] <- "wild"
  }
}

test_coefs_1 <- scan1coef(pr_p[, chr], in_cross_p$pheno[,contig_id], addcovar = metadata_p)
test_coefs_2 <-scan1coef(pr_p[, chr], in_cross_p$pheno[,contig_id])

test_1 <- test_coefs_1[which(rownames(coefs) == find_marker(map_p, chr = chr, pos = pos)),][1]
test_2 <- test_coefs_2[which(rownames(coefs) == find_marker(map_p, chr = chr, pos = pos)),][1]

test_1/ mean(in_cross_p$pheno[,contig_id],na.rm=TRUE)
test_2/ mean(in_cross_p$pheno[,contig_id],na.rm=TRUE)

#---for report
found_peaks_contig <- found_peaks_p
found_peaks_contig$effect <- round(found_peaks_contig$effect,2 )
found_peaks_contig$heritability <- round(found_peaks_contig$heritability,2 )
found_peaks_contig$lod <- round(found_peaks_contig$lod,2 )

write.table(found_peaks_contig, "peak_contig_for_report.tsv", quote=F, sep = "\t",row.names=F)


#------summarize the peaks
found_peaks_p$chr_pos <- paste(found_peaks_p$chr,found_peaks_p$pos, sep = "_")
genus_list <- as.character(unique(found_peaks_p$genus))
df_chr_pos <- data.frame(table(as.character(found_peaks_p$chr_pos[found_peaks_p$genus==genus_list[1]])))
df_chr_pos$genus <- genus_list[1]
colnames(df_chr_pos)[1] <- "chr_pos"
for (i in 2:11){
  genus_id <- genus_list[i]
  chr_pos <- as.character(found_peaks_p$chr_pos[found_peaks_p$genus==genus_id])
  new_chr_pos <- data.frame(table(chr_pos))
  new_chr_pos$genus <- genus_id
  df_chr_pos <- rbind(df_chr_pos,new_chr_pos)
}
df_chr_pos$chr_pos <- as.character(df_chr_pos$chr_pos)
df_chr_pos <- separate(df_chr_pos,col=chr_pos,into=c("chr","pos"),sep = "_")
#df_chr_pos$chr_pos <- paste(df_chr_pos$chr,df_chr_pos$pos,sep = "_")
df_chr_pos$chr <- as.numeric(df_chr_pos$chr )

df_chr_pos_sorted <- df_chr_pos[
  with(df_chr_pos, order(chr, pos)),
]
colnames(df_chr_pos_sorted) <- c("Chromosome", "Position (Mbp)", "# QTLs", "Linked bacetial genus")
write.table(df_chr_pos_sorted, "Desktop/summary_chr_pos_taxon.tsv", row.names = F, sep = "\t")
read.delim("Desktop/summary_chr_pos_taxon.tsv")

df_chr_pos_sorted_lod <- df_chr_pos_sorted
df_chr_pos_sorted_lod$chr_pos <- paste(df_chr_pos_sorted_lod[,1],df_chr_pos_sorted_lod[,2], sep = "_")
peak_lod_taxon$chr_pos <- paste(peak_lod_taxon$chr,peak_lod_taxon$pos, sep = "_")
