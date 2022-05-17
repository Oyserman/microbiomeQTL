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

setwd("data/final_assembly/qtl/metagenome_qtl/bin_level/20200620/")
rqtl_phenofile_p <- "../phenofile_all_good_bins.csv"

pfile <- read.csv(rqtl_phenofile_p)
pfile$id <- paste("p",pfile$id, sep = "")
pfile$id[70] <- "p276.1"
pfile$id[96] <- "p308.1"


#phenofile$id[!phenofile$id %in% rnaseq_map_geno$id] 
write.table(pfile, "phenofile_bin_depth.csv", sep = ",", quote = F, col.names=T, row.names=F)
phenofile_p <- read.csv("phenofile_bin_depth.csv")


rqtl_control_pink <- "rnaseq_map_control.yaml"

rnaseq_map_geno <- read.csv("Tomaat_Micro_biome_RNA_Seq_map/geno.csv")
bin_taxo <- read.delim("~/data/final_assembly/binning/new/new_checkm_lineage_wf.tsv")
good_bin_taxo <- bin_taxo[bin_taxo$Completeness>90&bin_taxo$Contamination<5,] 
write.table(good_bin_taxo,"checkm_good_bin.tsv",sep = "\t", quote = F, col.names=T, row.names=F)
ggplot(bin_taxo, aes(x=Completeness, y=Contamination)) +
  geom_point(size=1.5, shape=2,color="cornsilk4") +
  geom_hline(yintercept=5, linetype="dashed", 
             color = "steelblue", size=0.8) +
  geom_vline(xintercept = 90, linetype="dashed", 
             color = "steelblue", size=0.8)+
  geom_hline(yintercept=30, linetype="dashed", 
             color = "orange2", size=0.8) +
  geom_vline(xintercept = 70, linetype="dashed", 
             color = "orange2", size=0.8)

############ Import metadata #############
metadata_p <- read.table("pRIL_metadata_5_19_2019_formated.txt", sep = "\t",row.names = 1, header = T)

############ Check 0 in the counts #############
sum(colSums(phenofile_p==0, na.rm = T) >=1)
#plot(colSums(phenofile_p==0, na.rm = T)/nrow(phenofile_p)*100, ylim = c(0, 0.5))

######### Read QTL data from a set of files listed in the control.yaml #############
in_cross_p <- read_cross2(rqtl_control_pink, quiet = F)
check_c <- check_cross2(in_cross_p)

######### Run R/qtl2 #############
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

#plot(herit_p_no_psudo)
#l_herit_p_no_pseudo <- length(which(herit_p_no_psudo>=herit_thres))

#abline(h=0.2, lty=2, col="red")
herit_thres <- 0.2
l_herit_p <- length(which(herit_p>=herit_thres))
# TO DO #
# Write out a single table with all heritability scores
#df_herit_p<- data.frame(herit_p)
#df_herit_p$contig_id <- rownames(df_herit_p)
#write.table(df_herit_p, file = "results/Heritability_rhi_contig.txt", sep = "\t", quote = F, row.names = F, col.names = T)

#### Perform a genome scan-----------------------------------------------
# out <- scan1(pr, in_cross$pheno[,which(herit>herit_thres)], kinship = kinship, model = "normal")
out_p <- scan1(pr_p, in_cross_p$pheno, addcovar = metadata_p)
#out_p_heri <- scan1(pr_p, in_cross_p$pheno[,which(herit_p>=herit_thres)], addcovar = metadata_p)
#out_p_no_pseudo <- scan1(pr_p_no_pseudo, in_cross_p$pheno[,which(herit_p_no_psudo>=herit_thres)], addcovar = metadata_p)

#plot
par(mar=c(5.1, 4.1, 1.1, 1.1))
ymx <- maxlod(out_p) # overall maximum LOD score
#plot(out_p, map_p, lodcolumn=1, col="slateblue", ylim=c(0, ymx*1.02))
#plot(out_p, map_p, lodcolumn=2, col="violetred", add=TRUE)
#legend("topleft", lwd=2, col=c("slateblue", "violetred"), colnames(out_p), bg="gray90")

#### Run permutation test-----------------------------------------------
operm_p <- scan1perm(pr_p, in_cross_p$pheno, addcovar = metadata_p, n_perm=1000,cores = 8)
#operm_p_heri <- scan1perm(pr_p, in_cross_p$pheno[,which(herit_p>=herit_thres)], addcovar = metadata_p, n_perm=1000,cores = 20)

#operm_p_pseudo <- scan1perm(pr_p_no_pseudo, in_cross_p$pheno[,which(herit_p_no_psudo>=herit_thres)], addcovar = metadata_p, n_perm=10)

qlod_sum_p <- quantile(operm_p, c(0.8,0.85,0.9,0.95))
#qlod_sum_p_heri <- quantile(operm_p_heri, c(0.8,0.85,0.9,0.95))

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

#plot(out_p_no_pseudo,in_cross_p$gmap)


# write.table(lod_scores, file = output_lodscores, sep = "\t", quote = F, row.names = F, col.names = T)


# IMPORTANT, only using the permutation threshold from the combined global analysis

## use 0.95 quantile threshold
lod_inclusion_threshold_p <- as.numeric(qlod_sum_p[4])
#lod_inclusion_threshold_p_heri <- as.numeric(qlod_sum_p_heri[4])


#### Find lod peaks-----------------------------------------------
found_peaks_p <- find_peaks(out_p, map = map_p, threshold = lod_inclusion_threshold_p,prob = 0.95)

# found_peaks_p_heri <- find_peaks(out_p_heri, map = map_p, threshold = lod_inclusion_threshold_p_heri,prob = 0.95) %>%
#   mutate(contig_id = lodcolumn) %>%
#   left_join(contig_taxonomy, by = "contig_id")

found_peaks_p$lodcolumn <- found_peaks_p$lodcolumn

# found_peaks_p_no_pseudo <- find_peaks(out_p_no_pseudo, map =in_cross_p$gmap, threshold = lod_inclusion_threshold_p) %>%
#   mutate(contig_id = lodcolumn) %>%
#   left_join(contig_taxonomy, by = "contig_id")

write.table(found_peaks_p,"found_peaks_for_bin_permu1000.tsv", quote = F, col.names = T, row.names = F, sep = "\t")
#write.table(found_peaks_p,"results/found_peak_CI.tsv", quote = F, col.names = T, row.names = F, sep = "\t")

pdf("bin_peaks.pdf", width = 10.2, height = 4.5) 
plot(out_p,map_p,lodcolumn = "bin.254.tsv", ylim = c(0,ymx+0.2), main="bin 254")
abline(h= lod_inclusion_threshold_p , col=c(rgb(0,0,0,0.5)), lty = 2)

plot(out_p,map_p,lodcolumn = "bin.368.tsv", ylim = c(0,ymx+0.2),main="bin 368")
abline(h= lod_inclusion_threshold_p , col=c(rgb(0,0,0,0.5)), lty = 2)

plot(out_p,map_p,lodcolumn = "bin.492.tsv", ylim = c(0,ymx+0.2),main="bin 492")
abline(h= lod_inclusion_threshold_p , col=c(rgb(0,0,0,0.5)), lty = 2)

plot(out_p,map_p,lodcolumn = "bin.516.tsv", ylim = c(0,ymx+0.2), main="bin 516")
abline(h= lod_inclusion_threshold_p , col=c(rgb(0,0,0,0.5)), lty = 2)

plot(out_p,map_p,lodcolumn = "bin.72.tsv", ylim = c(0,ymx+0.2), main="bin 72")
abline(h= lod_inclusion_threshold_p , col=c(rgb(0,0,0,0.5)), lty = 2)

dev.off()

par(mfrow=c(3,2),mar=c(5.4, 4.1, 2, 2))
pdf("~/data/final_assembly/qtl/metagenome_qtl/bin_level/20200620/bin_peaks_in_A_plot2*3.pdf", width = 22.5, height = 12) 
plot(out_p,map_p,lodcolumn = "bin.254.tsv", ylim = c(0,ymx+0.2), main="bin 254")
abline(h= lod_inclusion_threshold_p , col=c(rgb(0,0,0,0.5)), lty = 2)

plot(out_p,map_p,lodcolumn = "bin.368.tsv", ylim = c(0,ymx+0.2),main="bin 368")
abline(h= lod_inclusion_threshold_p , col=c(rgb(0,0,0,0.5)), lty = 2)

plot(out_p,map_p,lodcolumn = "bin.492.tsv", ylim = c(0,ymx+0.2),main="bin 492")
abline(h= lod_inclusion_threshold_p , col=c(rgb(0,0,0,0.5)), lty = 2)

plot(out_p,map_p,lodcolumn = "bin.516.tsv", ylim = c(0,ymx+0.2), main="bin 516")
abline(h= lod_inclusion_threshold_p , col=c(rgb(0,0,0,0.5)), lty = 2)

plot(out_p,map_p,lodcolumn = "bin.72.tsv", ylim = c(0,ymx+0.2), main="bin 72")
abline(h= lod_inclusion_threshold_p , col=c(rgb(0,0,0,0.5)), lty = 2)

dev.off()


bin_taxo$Bin.Id <- paste(bin_taxo$Bin.Id,"tsv", sep = ".")

bin_taxo[bin_taxo$Bin.Id %in% found_peaks_p$lodcolumn,]

setwd("~/data/final_assembly/depth/new_bin_depth/")
read_num <- read.delim("total_read_num.txt", header = F, stringsAsFactors = F)

for (bin in unique(found_peaks_p$lodcolumn)){
  print(bin)
  bin_depth_file <- read.delim(bin, header = F,stringsAsFactors = F)
  bin_depth_input <- cbind(read_num$V1, read_num$V2, bin_depth_file$V2)
  colnames(bin_depth_input) <- c("id","read_num","depth")
  bin_depth_input <- data.frame(bin_depth_input, stringsAsFactors = F)
  bin_depth_input$read_num <- sapply(bin_depth_input$read_num, as.numeric)
  bin_depth_input$depth <- sapply(bin_depth_input$depth, as.numeric)
  bin_depth_input$norm_depth <- bin_depth_input$depth/bin_depth_input$read_num*1000000
  bin_depth_input$group <- c(rep("RIL",96),rep("BULK",7),rep("MM",6),rep("P",5))
  for (i in 1:114){
  bin_depth_input$id[i] <- str_split(bin_depth_input$id,"_")[[i]][1]
  }
  write.table(bin_depth_input, paste("norm",bin, sep = "."), sep = "\t", quote = F, col.names=T, row.names=F)
}

temp <- list.files(path="./",pattern="norm*")
pdf("~/data/final_assembly/qtl/metagenome_qtl/bin_level/20200620/peak_bin_depth_line.pdf",width = 10.2, height = 5.5)

bin254_depth <- read.delim(temp[1],stringsAsFactors = F)
bin368_depth <- read.delim(temp[2],stringsAsFactors = F)
bin492_depth <- read.delim(temp[3],stringsAsFactors = F)
bin516_depth <- read.delim(temp[4],stringsAsFactors = F)
bin72_depth <- read.delim(temp[5],stringsAsFactors = F)

all_bin_depth <- data.frame(bin254_depth$id,bin254_depth$group,bin72_depth$norm_depth,bin254_depth$norm_depth,bin368_depth$norm_depth,bin492_depth$norm_depth,bin516_depth$norm_depth)
colnames(all_bin_depth) <- c("id","group","bin72","bin254","bin368","bin492","bin516")
all_bin_depth <- melt(all_bin_depth)
colnames(all_bin_depth) <- c("id","group","bin","norm_depth")
g <- ggplot(all_bin_depth,aes(id, norm_depth,group=bin,col=bin))+ 
  geom_line()+
  geom_point(aes(shape=group))+
  ggtitle("Relative abundance of genomic bins") +
  ylab("Relative abundance") +
  xlab("Sample") +
  theme(plot.title = element_text(hjust = 0.5), axis.text = element_text(angle = 90,size = 4))
print(g)

dev.off()

g <- ggplot(bin492_depth,aes(id, norm_depth,fill=group))+ 
  geom_bar(stat = "identity") +
  ggtitle("Relative abundance of genomic bins") +
  ylab("Relative abundance") +
  xlab("Sample") +
  theme(plot.title = element_text(hjust = 0.5), axis.text = element_text(angle = 90,size = 4))
print(g)

par(mar=c(5, 4, 4, 8)+0.1)
cols<-rainbow(length(unique(found_peaks_p$lodcolumn)),alpha = 0.8)
plot(out_p, map_p, lodcolumn=unique(found_peaks_p$lodcolumn)[1],lty = 2, ylim = c(0,maxlod(out_p)+0.2))
for (i in 2:length(unique(found_peaks_p$lodcolumn))){
  plot(out_p, map_p, lodcolumn=unique(found_peaks_p$lodcolumn[i]),lty = 2, add=TRUE)}
abline(h= lod_inclusion_threshold_p , col=c(rgb(0,0,0,0.5)), lty = 2)
op <- par(cex=0.6)
legend("right", inset = c(-0.17,0), legend = unique(found_peaks_p$lodcolumn), xpd = TRUE, 
       horiz = FALSE, col = rainbow(12), lty = 1, bty = "n")


for (i in 1:nrow(found_peaks_p)){
  bin_id <- found_peaks_p$lodcolumn[i]
  chr <- as.character(found_peaks_p$chr[i])
  pos <- found_peaks_p$pos[i]
  coefs <- scan1coef(pr_p[, chr], in_cross_p$pheno[,bin_id], addcovar = metadata_p)
  col <- c("slateblue", "violetred", "green3")
  plot(coefs, map_p[chr], columns=c(1,2,7), col=col)
  
  coef <- coefs[which(rownames(coefs) == find_marker(map_p, chr = chr, pos = pos)),][1]
  found_peaks_p$effect[i] <-  coef/ mean(in_cross_p$pheno[,bin_id],na.rm=TRUE)
  found_peaks_p$heritability[i] <- herit_p[bin_id]
  #print(coef[1])
  #if(coef[1]>0) print("Non-negative number") else print("Negative number")
  if(coef[1]>0){
    print(paste(bin_id,"AA, modern"))
    found_peaks_p$allele[i] <- "modern"
  }
  else{
    print(paste(bin_id, "BB, wild"))
    found_peaks_p$allele[i] <- "wild"
  }
}
write.table(found_peaks_p,"found_peaks_for_bin_permu1000_full.tsv", quote = F, col.names = T, row.names = F, sep = "\t")
g <- maxmarg(pr_p, map_p, chr=6, pos=37.961594, return_char=TRUE)
par(mar=c(4.1, 4.1, 0.6, 0.6))
plot_pxg(g, in_cross_p$pheno[,"bin.72.tsv"], ylab="bin.72")

plot_peaks(found_peaks_p, map_p, lod_labels = T)

plot(out_p, map_p$"6", lodcolumn=1,lty = 1, ylim = c(0,maxlod(out_p)+0.2))
for (i in 2:ncol(out_p)){
  plot(out_p, map_p$"6", lodcolumn=i,lty = 2, add=TRUE)}
abline(h= lod_inclusion_threshold_p , col=c(rgb(0,0,0,0.5)), lty = 2)




#---for report
found_peaks_bin <- found_peaks_p
found_peaks_bin$effect <- round(found_peaks_bin$effect,2 )
found_peaks_bin$heritability <- round(found_peaks_bin$heritability,2 )
found_peaks_bin$lod <- round(found_peaks_bin$lod,2 )

write.table(found_peaks_bin, "peak_bin_for_report.tsv", quote=F, sep = "\t",row.names=F)
