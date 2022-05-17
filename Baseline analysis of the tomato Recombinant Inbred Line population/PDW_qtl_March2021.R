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

# LOAD THE IMAGE

argv <- commandArgs(trailingOnly = T)
theme_set(theme_pubr(border = T, legend = "right"))
max_treads <- 4
setwd("~/PDW_QTL/")
# Load the Khan and current PDW datasets
K_2012 <- read.csv("kazmi2012_raw_data_reformatted_forQTL.txt", header = TRUE)
K_2012_abr <- K_2012[,which(colnames(K_2012) %in% c("id", "dry_shoot","dry_root"))]
K_2012_abr$id <- paste("k",K_2012_abr$id, sep="")
# write.csv(K_2012_abr,"~/PDW_QTL/khan_pdw_2012.txt", header=TRUE, row.names = FALSE)


# load R/qtl2 data --------------------------------------------------------

rqtl_phenofile <- "~/PDW_QTL/khan_pdw_2012.txt"
rqtl_phenofile <- "~/PDW_QTL/RIL_metadata_replicates_5_19_2019.csv"

# Load the yaml files for the microQTL and for the PDW
rqtl_rnaseq_map <- "~/PDW_QTL/rnaseq_map.yaml"
rqtl_rnaseq_map <- "~/PDW_QTL/rnaseq_PDW_map.yaml"

############ Classic Phenotypes First #############

######### Defining outputs names
#########
output_basename_r <- tools::file_path_sans_ext(rqtl_phenofile)
output_lodscores_r <- sprintf("%s_lod.csv", rqtl_phenofile)
output_plodscores_r <- sprintf("%s_plod.csv", rqtl_phenofile)
output_peaks_r <- sprintf("%s_peaks.csv", rqtl_phenofile)
output_heritability_r <- sprintf("%s_heritability.csv", rqtl_phenofile)
output_heatmap_r <-sprintf("%s_lod.pdf", rqtl_phenofile)

#########

in_cross_pdw <- read_cross2(rqtl_rnaseq_map, quiet = F)
color_covar <- in_cross_pdw$pheno[,1]
in_cross_pdw$pheno <- in_cross_pdw$pheno[,-1]

# Insert pseudomarkers for the gmap used in each analysis
map <- insert_pseudomarkers(in_cross_pdw$gmap, step = 1, cores = 1)

# make metadata -----------------------------------------------------------

p_m_n <- names(unlist(map, use.names=TRUE, recursive = TRUE))

p_m_n_chr <- NULL
p_m_n_marker_id <- NULL
for (i in 1:length(p_m_n)) {
  p_m_n_chr[as.numeric(i)] <- strsplit(p_m_n[as.numeric(i)], "[.]")[[1]][1]
 # p_m_n_marker_id[i] <- strsplit(p_m_n[i], "[.]")[[1]][2]
  p_m_n_marker_id[i] <- gsub("^.*\\.","", p_m_n[as.numeric(i)])
}

new_gmap <- as.data.frame(cbind(p_m_n_marker_id, p_m_n_chr, unlist(map, use.names=TRUE, recursive = TRUE)))
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

pr <- calc_genoprob(in_cross_pdw, map, error_prob = 1e-4, cores = max_treads)

cat("Calculating a kinship matrix...\n")
# kinship_n <- calc_kinship(pr_n)

kinship <- calc_kinship(pr)

###### FOr the pink and yellow, as above but without the replicate
herit <- est_herit(in_cross_pdw$pheno, kinship, addcovar = color_covar)

# out_n <- scan1(pr_n, in_cross_r_new$pheno, addcovar = metadata_r)
out <- scan1(pr, in_cross_pdw$pheno, addcovar = color_covar)

# plot(out[,4])
# operm_n <- scan1perm(pr_n, in_cross_r_new$pheno, addcovar = metadata_r, n_perm=10)
operm <- scan1perm(pr, in_cross_pdw$pheno, n_perm = 1000, cores = 6, addcovar = color_covar)

# qlod_sum_n <- quantile(operm_n, c(0.8,0.85,0.9,0.95))
qlod_sum <- quantile(operm, c(0.8,0.85,0.9,0.95))

##########
df_out <- as.data.frame(out)

out_lod_scores <- cbind(new_gmap,df_out)

lod_inclusion_threshold_pdw <- as.numeric(qlod_sum[4])

cat("Finding peaks...\n")

# found_peaks <- find_peaks(out, map, threshold = lod_inclusion_threshold) %>%
found_peaks_pdw <- find_peaks(out, map = map, threshold = as.numeric(qlod_sum[4]), prob = 0.95)

calc_effects <- function(found_peaks, pr, covar, incross, map) {

  for (i in 1:nrow(found_peaks)){
  lodcolum_id <- found_peaks$lodcolumn[i]
  chr <- as.character(found_peaks$chr[i])
  pos <- found_peaks$pos[i]
  coefs <- scan1coef(pr[, chr], incross$pheno[,lodcolum_id], addcovar = color_covar)
  coef <- coefs[which(rownames(coefs) == find_marker(map, chr = chr, pos = pos)),][1]
  found_peaks$coef[i] <-  coef
  found_peaks$effect[i] <-  coef / mean(incross$pheno[,lodcolum_id],na.rm=TRUE) * 100
  found_peaks$heritability[i] <- herit[lodcolum_id]
  #print(coef[1])
  #if(coef[1]>0) print(“Non-negative number”) else print(“Negative number”)
  if (coef[1]>0) {
    print(paste(lodcolum_id,"AA, modern",sep=" "))
    found_peaks$allele[i] <- "modern"
  }
  else {
    print(paste(lodcolum_id, "BB, wild"))
    found_peaks$allele[i] <- "wild"
  }
}
  return(found_peaks)
}

found_peaks_pdw_effects <- calc_effects(found_peaks_pdw, pr, covar, in_cross_pdw, map) 

plot(out[,5] , type = "l")


###########################
# PLANT DRY WEIGHT FIGURE #
###########################
close.screen(all = TRUE)

dev.off()

m <- rbind(c(0.2,0.79,0.5,1),
           c(0.2,0.39,0.01,0.5),
           c(0.4,0.59,0.01,0.5),
           c(0.6,0.79,0.01,0.5))

split.screen(m)
ymx <- max(found_peaks_pdw_effects$lod)
ymx2 <- max(in_cross_pdw$pheno[,5],na.rm=TRUE)
screen(1)
par(mar = c(3, 2, 2, 0), cex.axis = 0.75)
plot(out, map, lodcolumn=5, col=rgb(0,0,1,0.5), ylim = c(0,ymx+0.2), ylab="", xlab="")
abline(h= summary(operm)[4], col=c(rgb(0,0,1,0.5)), lwd = 2)
mtext(side = 2, "LOD",line = 1.25, cex = 0.75)
mtext(side = 1, "Chromosome",line = 1.25, cex = 0.75)
mtext(side = 3, "Plant Dry Weight",line = 0.5, cex = 0.75)

# p PDW chr 2 pos 42
gg <- maxmarg(pr, map, chr = found_peaks_pdw_effects$chr[3], pos = found_peaks_pdw_effects$pos[3], return_char=TRUE)
screen(2)
par(mar = c(2, 2, 1, 0))
plot_pxg(gg, in_cross_pdw$pheno[,5], SEmult=2, sort=FALSE,ylim = c(0,ymx2+0.2), col=rgb(1,0,0,alpha=0.2),pch = 19, seg_lwd=2, seg_width=0.4, vlines=NA, cex=0.75, ylab="")
mtext(side = 3, paste("chr", found_peaks_pdw_effects$chr[3],"pos", found_peaks_pdw_effects$pos[3], sep = " "), line = 1, cex = 0.5)
mtext(side = 3, paste("Effect ", abs(round(found_peaks_pdw_effects$effect[3],4)), "%"), line = 0, cex = 0.5)
mtext(side = 2, "Shoot Dry Weight (g)",line = 1.25, cex = 0.75)
mtext(side = 1, "Allele",line = 1.25, cex = 0.75)

gg <- maxmarg(pr, map, chr = found_peaks_pdw_effects$chr[4], pos = found_peaks_pdw_effects$pos[4], return_char=TRUE)
screen(3)
par(mar = c(2, 2, 1, 0))
plot_pxg(gg, in_cross_pdw$pheno[,5], SEmult = 2, sort=FALSE, ylim = c(0,ymx2+0.2), col=rgb(1,0,0,alpha=0.2), pch = 19, seg_lwd=2, seg_width=0.4, vlines=NA, cex=0.75, ylab="")
mtext(side = 3, paste("chr", found_peaks_pdw_effects$chr[4],"pos", found_peaks_pdw_effects$pos[4], sep = " "), line = 1, cex = 0.5)
mtext(side = 3, paste("Effect ", abs(round(found_peaks_pdw_effects$effect[4],4)), "%"), line = 0, cex = 0.5)
mtext(side = 2, "Shoot Dry Weight (g)",line = 1.25, cex = 0.75)
mtext(side = 1, "Allele",line = 1.25, cex = 0.75)

#######################################
#### Is this effect additive? PDW  ####
#######################################

gg1 <- maxmarg(pr, map, chr = found_peaks_pdw_effects$chr[3], pos = found_peaks_pdw_effects$pos[3], return_char=TRUE) %>% as.data.frame()
gg2 <- maxmarg(pr, map, chr = found_peaks_pdw_effects$chr[4], pos = found_peaks_pdw_effects$pos[4], return_char=TRUE) %>% as.data.frame()

gg12 <-cbind(gg1, gg2)
colnames(gg12) <- c("chr2", "chr9")
gg12_meta <-  merge(gg12, in_cross_pdw$pheno, by=0, all=TRUE)

gg12_meta$merged <- paste(gg12_meta$chr2, gg12_meta$chr9, sep="_")

table(gg12_meta$merged)

# chr 2 wild, chr 9 modern
# BB_AA x2, BB_BB x1, AA_AA x1, AA_BB x 0

gg12_meta$additive_alleles <- "zero"
gg12_meta$additive_alleles[which(gg12_meta$merged == names(table(gg12_meta$merged)[1]))] <- "one"
gg12_meta$additive_alleles[which(gg12_meta$merged == names(table(gg12_meta$merged)[2]))] <- "zero"
gg12_meta$additive_alleles[which(gg12_meta$merged == names(table(gg12_meta$merged)[3]))] <- "two"
gg12_meta$additive_alleles[which(gg12_meta$merged == names(table(gg12_meta$merged)[4]))] <- "one"
gg12_meta$additive_alleles <- factor(gg12_meta$additive_alleles, levels = c("zero","one","two"))

res.aov <- aov(PDW ~ additive_alleles, data = gg12_meta)
summary(res.aov)
TUKEY <- TukeyHSD(res.aov)

# plot(TUKEY , las=1 , col="brown")


generate_label_df <- function(TUKEY, variable){
  
  # Extract labels and factor levels from Tukey post-hoc
  Tukey.levels <- TUKEY[[variable]][,4]
  Tukey.labels <- data.frame(multcompLetters(Tukey.levels)['Letters'])
  
  #I need to put the labels in the same order as in the boxplot :
  Tukey.labels$treatment=rownames(Tukey.labels)
  Tukey.labels=Tukey.labels[order(Tukey.labels$treatment) , ]
  return(Tukey.labels)
}

# Apply the function on my dataset
LABELS <- generate_label_df(TUKEY, "additive_alleles")

my_colors <- c(
  rgb(143,199,74,maxColorValue = 255),
  rgb(242,104,34,maxColorValue = 255),
  rgb(111,145,202,maxColorValue = 255)
)

# Draw the basic boxplot
labels_n <- c(paste("n=", length(which(gg12_meta$additive_alleles=="zero"))), 
              paste("n=", length(which(gg12_meta$additive_alleles=="one"))),
              paste("n=", length(which(gg12_meta$additive_alleles=="two"))))

screen(4)
par(mar = c(2, 2, 1, 0))
a <- boxplot(gg12_meta$PDW ~ gg12_meta$additive_alleles, 
             ylab = "Plant Dry Weight (g)",  xlab = "Number of associated alleles", 
             ylim = c(0,1), cex.axis = 0.5, main = "Additivity in Plant Dry Weight",
             col=my_colors, cex.main = 0.5, pch = 1, cex = 0.5)

text(x=c(1,2,3), y = rep(0.8, 3), labels = labels_n, cex = 0.5)

# I want to write the letter over each box. Over is how high I want to write it.
over <- 0.1

#Add the labels
text(c(1:length(unique(gg12_meta$additive_alleles))), over , LABELS[,1], col=my_colors[as.numeric(LABELS[,1])], cex = 0.5)

##############
# Rhizosphere #
##############

close.screen(all = TRUE)

dev.off()

m <- rbind(c(0.2,0.79,0.5,1),
           c(0.2,0.39,0.01,0.5),
           c(0.4,0.59,0.01,0.5),
           c(0.6,0.79,0.01,0.5))

split.screen(m)
ymx <- max(found_peaks_pdw_effects$lod)
ymx2 <- max(in_cross_pdw$pheno[,3],na.rm=TRUE)
screen(1)
par(mar = c(3, 2, 2, 0), cex.axis = 0.75)
plot(out, map, lodcolumn=3, col=rgb(0,0,1,0.5), ylim = c(0,ymx+0.2), ylab="", xlab="")
abline(h= summary(operm)[4], col=c(rgb(0,0,1,0.5)), lwd = 2)
mtext(side = 2, "LOD",line = 1.25, cex = 0.75)
mtext(side = 1, "Chromosome",line = 1.25, cex = 0.75)
mtext(side = 3, "Rhizosphere Mass",line = 0.5, cex = 0.75)

# p PDW chr 2 pos 42
gg <- maxmarg(pr, map, chr = found_peaks_pdw_effects$chr[1], pos = found_peaks_pdw_effects$pos[1], return_char=TRUE)
screen(2)
par(mar = c(2, 2, 1, 0))
plot_pxg(gg, in_cross_pdw$pheno[,3], SEmult=2, sort=FALSE,ylim = c(0,ymx2+0.2), col=rgb(1,0,0,alpha=0.2),pch = 19, seg_lwd=2, seg_width=0.4, vlines=NA, cex=0.75, ylab="")
mtext(side = 3, paste("chr", found_peaks_pdw_effects$chr[1],"pos", found_peaks_pdw_effects$pos[1], sep = " "), line = 1, cex = 0.5)
mtext(side = 3, paste("Effect ", abs(round(found_peaks_pdw_effects$effect[1],4)), "%"), line = 0, cex = 0.5)
mtext(side = 2, "Rhizosphere Mass (g)",line = 1.25, cex = 0.75)
mtext(side = 1, "Allele",line = 1.25, cex = 0.75)

gg <- maxmarg(pr, map, chr = found_peaks_pdw_effects$chr[2], pos = found_peaks_pdw_effects$pos[2], return_char=TRUE)
screen(3)
par(mar = c(2, 2, 1, 0))
plot_pxg(gg, in_cross_pdw$pheno[,3], SEmult = 2, sort=FALSE, ylim = c(0,ymx2+0.2), col=rgb(1,0,0,alpha=0.2), pch = 19, seg_lwd=2, seg_width=0.4, vlines=NA, cex=0.75, ylab="")
mtext(side = 3, paste("chr", found_peaks_pdw_effects$chr[2],"pos", found_peaks_pdw_effects$pos[2], sep = " "), line = 1, cex = 0.5)
mtext(side = 3, paste("Effect ", abs(round(found_peaks_pdw_effects$effect[2],4)), "%"), line = 0, cex = 0.5)
mtext(side = 2, "Rhizosphere Mass (g)",line = 1.25, cex = 0.75)
mtext(side = 1, "Allele",line = 1.25, cex = 0.75)

#######################################
#### Is this effect additive? PDW  ####
#######################################

gg1 <- maxmarg(pr, map, chr = found_peaks_pdw_effects$chr[1], pos = found_peaks_pdw_effects$pos[1], return_char=TRUE) %>% as.data.frame()
gg2 <- maxmarg(pr, map, chr = found_peaks_pdw_effects$chr[2], pos = found_peaks_pdw_effects$pos[2], return_char=TRUE) %>% as.data.frame()

gg12 <-cbind(gg1, gg2)
colnames(gg12) <- c("chr5", "chr9")
gg12_meta <-  merge(gg12, in_cross_pdw$pheno, by=0, all=TRUE)

gg12_meta$merged <- paste(gg12_meta$chr5, gg12_meta$chr9, sep="_")

table(gg12_meta$merged)

# chr 2 wild, chr 9 modern
# BB_AA x2, BB_BB x1, AA_AA x1, AA_BB x 0

gg12_meta$additive_alleles <- "zero"
gg12_meta$additive_alleles[which(gg12_meta$merged == names(table(gg12_meta$merged)[1]))] <- "one"
gg12_meta$additive_alleles[which(gg12_meta$merged == names(table(gg12_meta$merged)[2]))] <- "zero"
gg12_meta$additive_alleles[which(gg12_meta$merged == names(table(gg12_meta$merged)[3]))] <- "two"
gg12_meta$additive_alleles[which(gg12_meta$merged == names(table(gg12_meta$merged)[4]))] <- "one"
gg12_meta$additive_alleles <- factor(gg12_meta$additive_alleles, levels = c("zero","one","two"))

res.aov <- aov(Rhizosphere ~ additive_alleles, data = gg12_meta)
summary(res.aov)
TUKEY <- TukeyHSD(res.aov)

# plot(TUKEY , las=1 , col="brown")

# Apply the function on my dataset
LABELS <- generate_label_df(TUKEY, "additive_alleles")[c(3,1,2),]

my_colors <- c(
  rgb(143,199,74,maxColorValue = 255),
  rgb(242,104,34,maxColorValue = 255),
  rgb(242,104,34,maxColorValue = 255)
)

# Draw the basic boxplot
labels_n <- c(paste("n=", length(which(gg12_meta$additive_alleles=="zero"))), 
              paste("n=", length(which(gg12_meta$additive_alleles=="one"))),
              paste("n=", length(which(gg12_meta$additive_alleles=="two"))))

screen(4)
par(mar = c(2, 2, 1, 0))
a <- boxplot(gg12_meta$Rhizosphere ~ gg12_meta$additive_alleles, 
             ylab = "Rhizosphere Mass (g)",  xlab = "Number of associated alleles", 
             ylim = c(0,5), cex.axis = 0.5, main = "Additivity in Rhizosphere Mass",
             col=my_colors, cex.main = 0.5, pch = 1, cex = 0.5)

text(x=c(1,2,3), y = rep(0.8, 3), labels = labels_n, cex = 0.5)

# I want to write the letter over each box. Over is how high I want to write it.
over <- 0.1

#Add the labels
text(c(1:length(unique(gg12_meta$additive_alleles))), over , LABELS[,1], col=my_colors[as.numeric(LABELS[,1])], cex = 0.5)

# Save the image
setwd("/home/nioo/beno/PDW_QTL")
getwd()
save.image(file = "PDW_qtl_March2021.RData")
load("PDW_qtl_March2021.RData")

