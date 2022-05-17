### this version fixes 0-based indexing issues from inStrain -
# bed file for use by samtools depth will be made 1-based, but also a 0-based refence will be made for other use
# GFF file references in this file will be changed to be 0-based
file_info = "210708"

setwd("D:/Documents/WUR/MSc Thesis MBF/Analyses/R_analysis")
# install.packages("ggplot2")
# install.packages("magrittr")
# if (!requireNamespace("BiocManager", quietly = TRUE)){install.packages("BiocManager")}
# BiocManager::install("goseq")

# if (!requireNamespace("BiocManager", quietly=TRUE))
#   install.packages("BiocManager")
# BiocManager::install("topGO")


library(ggplot2)
library(magrittr)
# library(goseq)
library(topGO)
library(GO.db)

# browseVignettes("goseq")

# create output directories
if (!file.exists("./bed_files")){dir.create(file.path("./bed_files"))}
if (!file.exists("./analysis_output")){dir.create(file.path("./analysis_output"))}
if (!file.exists("./snp_data")){dir.create(file.path("./snp_data"))}
if (!file.exists("./figures")){dir.create(file.path("./figures"))}
if (!file.exists("./GO_analysis")){dir.create(file.path("./GO_analysis"))}

# input files -------------------------------
taxonomy = read.table(file = "contig_10kb_taxonomy.tsv", header = TRUE)

IS_compare_file = read.table("IS.COMPARE_comparisonsTable.tsv", header=TRUE) # IS compare output
is_MM_SNVs = read.table("filtered_rhi_contig_id.bin_merged_MM.IS_SNVs.tsv", header=TRUE, sep="\t") # IS profile output
is_P_SNVs = read.table("filtered_rhi_contig_id.bin_merged_P.IS_SNVs.tsv", header=TRUE, sep = "\t")

gff_file = read.table(file="./Prokka_output/filtered_rhi_contig_id_210111_no_fasta.gff", sep = "\t", header=F, quote="\"")
colnames(gff_file) = c("contig_id_short", "source", "feature", "start", "end", "score", "strand", "phase", "attributes")
# convert GFF to 0-based start and end positions to be compatible with inStrain
gff_file$start = gff_file$start - 1
gff_file$end = gff_file$end - 1


# prokka_fna = read.table(file= "./Prokka_output/filtered_rhi_contig_id_210111_fna_210113.fna", sep = "\t", header=F)

prokka_gene_id_list = read.table(file="./Prokka_output/filtered_rhi_contig_id_210111_fna_210113_gene_list.tsv", sep = "\t", header = F, quote="\"")
colnames(prokka_gene_id_list) = c("contig_gene", "prokka_gene_id", "annotation")

eggnog_output_file = read.table(file="./eggNOG_output/rhizosphere_enriched_contig_gene_predictions_prokka_210111_fna_210113_eggNOG_annotations_210206.tsv", sep = "\t", quote="\"")
colnames(eggnog_output_file) = c("Query",	"Seed_ortholog",	"e-value",	"score",	"best_tax_lvl",	"Preferred_name",	"GO_terms",	"EC_number",	"KEGG_KO",	"KEGG_pathway",	"KEGG_module",	"KEGG_reaction",	"KEGG_rclass",	"BRITE",	"KEGG_TC",	"CAZy",	"BiGG_reaction",	"annot_lvl",	"matching_OGs",	"Best_OG",	"COG_cat",	"description")

# contig files for specific analyses:
contigs_streptomyces = read.table(file="./contig_selections/Streptomyces.list", header = F, col.names = "contig")
contigs_cellvibrio = read.table(file="./contig_selections/Cellvibrio.list", header = F, col.names = "contig")

# inStrain compare output analysis ------------------

pop_snps = merge(IS_compare_file[which(IS_compare_file$population_SNPs >0),], taxonomy, by.x = "scaffold", by.y="contig_id")

# inStrain profile output --------------------

is_P_SNVs$snv_pos = paste(is_P_SNVs$scaffold,is_P_SNVs$position, sep="_pos_")
is_MM_SNVs$snv_pos = paste(is_MM_SNVs$scaffold,is_MM_SNVs$position, sep="_pos_")

### Get SNPs that are in MM or P, but not in the other
# the sequencing depth in these loci then needs to be extracted, before the rest of this script can be continued! Use my python script for that
is_MM_SNVs_unique = subset(is_MM_SNVs[which(!(is_MM_SNVs$snv_pos %in% is_P_SNVs$snv_pos)),], select = c("scaffold", "position"))
is_P_SNVs_unique = subset(is_P_SNVs[which(!(is_P_SNVs$snv_pos %in% is_MM_SNVs$snv_pos)),], select = c("scaffold", "position"))

# 1-based BED files for working with samtools
is_MM_SNVs_unique_1based = is_MM_SNVs_unique
is_P_SNVs_unique_1based = is_P_SNVs_unique
is_MM_SNVs_unique_1based$position = is_MM_SNVs_unique_1based$position + 1
is_P_SNVs_unique_1based$position = is_P_SNVs_unique_1based$position + 1

write.table(is_MM_SNVs_unique_1based, file = "./bed_files/is_minmapq_0_MM_SNVs_unique_bed_1based_210519.tsv", quote=FALSE, sep="\t", row.names = F, col.names=F)
write.table(is_P_SNVs_unique_1based, file = "./bed_files/is_minmapq_0_P_SNVs_unique_bed_1based_210519.tsv", quote=FALSE, sep="\t", row.names = F, col.names=F)


# load in the coverage files -----------------
# # old 0-based coverage files
# is_MM_SNVs_unique_P_coverage = read.table("./count_unique_SNPs/MM_unique_snps.merged_P_filtered_rhizosphere_minmapq_0_210114.tsv", header = F)
# is_P_SNVs_unique_MM_coverage = read.table("./count_unique_SNPs/P_unique_snps.merged_MM_filtered_rhizosphere_minmapq_0_210114.tsv", header = F)

# correct coverage files: 1-based indexing of samtools
is_MM_SNVs_unique_P_coverage = read.table("./count_unique_SNPs/MM_unique_snps.merged_P_filtered_rhizosphere_minmapq_0_1based_210519.tsv", header = F)
is_P_SNVs_unique_MM_coverage = read.table("./count_unique_SNPs/P_unique_snps.merged_MM_filtered_rhizosphere_minmapq_0_1based_210519.tsv", header = F)

colnames(is_MM_SNVs_unique_P_coverage) = c("scaffold", "position", "position_coverage")
colnames(is_P_SNVs_unique_MM_coverage) = c("scaffold", "position", "position_coverage")

# make 0-based again
is_MM_SNVs_unique_P_coverage$position = is_MM_SNVs_unique_P_coverage$position - 1
is_P_SNVs_unique_MM_coverage$position = is_P_SNVs_unique_MM_coverage$position - 1

### fuse with original datasets -------------------
# 1. merge MM_unique_P_cov with MM_unique to get info about reference base
# add SNV_pos to merge the sets
is_MM_SNVs_unique_P_coverage$snv_pos = paste(is_MM_SNVs_unique_P_coverage$scaffold,is_MM_SNVs_unique_P_coverage$position, sep="_pos_")
is_P_SNVs_unique_MM_coverage$snv_pos = paste(is_P_SNVs_unique_MM_coverage$scaffold,is_P_SNVs_unique_MM_coverage$position, sep="_pos_")
is_MM_SNVs_unique$snv_pos = paste(is_MM_SNVs_unique$scaffold,is_MM_SNVs_unique$position, sep="_pos_")
is_P_SNVs_unique$snv_pos = paste(is_P_SNVs_unique$scaffold,is_P_SNVs_unique$position, sep="_pos_")
full_is_MM_SNVs_unique_P_coverage = merge(subset(is_MM_SNVs, select = -c(position_coverage)),
                                          subset(is_MM_SNVs_unique_P_coverage, select=c(snv_pos,position_coverage)), by="snv_pos")
full_is_P_SNVs_unique_MM_coverage = merge(subset(is_P_SNVs, select = -c(position_coverage)),
                                          subset(is_P_SNVs_unique_MM_coverage, select=c(snv_pos,position_coverage)), by="snv_pos")

# full_is_MM_SNVs_unique_P_coverage
full_is_MM_SNVs_unique_P_coverage$con_base = full_is_MM_SNVs_unique_P_coverage$ref_base
full_is_MM_SNVs_unique_P_coverage$var_base = full_is_MM_SNVs_unique_P_coverage$ref_base
full_is_MM_SNVs_unique_P_coverage$allele_count = 1
full_is_MM_SNVs_unique_P_coverage$class = "SNV"
for (base in c("A","C","T","G")){ # reset coverage
  full_is_MM_SNVs_unique_P_coverage[[base]] = 0
}

for (row in 1:nrow(full_is_MM_SNVs_unique_P_coverage)){
  full_is_MM_SNVs_unique_P_coverage[[full_is_MM_SNVs_unique_P_coverage$ref_base[row]]][row] = full_is_MM_SNVs_unique_P_coverage$position_coverage[row] # add correct coverage to allele of reference
}

full_is_MM_SNVs_unique_P_coverage[full_is_MM_SNVs_unique_P_coverage$mutation_type == "S" | full_is_MM_SNVs_unique_P_coverage$mutation_type == "N", ]$mutation_type = "S"
full_is_MM_SNVs_unique_P_coverage[full_is_MM_SNVs_unique_P_coverage$mutation_type == "S", ]$mutation = "S:"
full_is_MM_SNVs_unique_P_coverage[full_is_MM_SNVs_unique_P_coverage$mutation_type == "S", ]$cryptic = "False" # all S are never cryptic? Doesn't really matter as I don't use this

full_is_MM_SNVs_unique_P_coverage$ref_freq = 0 # reset frequencies
full_is_MM_SNVs_unique_P_coverage$con_freq = 0
full_is_MM_SNVs_unique_P_coverage$var_freq = 0

full_is_MM_SNVs_unique_P_coverage[full_is_MM_SNVs_unique_P_coverage$position_coverage > 0, ]$ref_freq = 1 # if coverage is present, put frequency to 1
full_is_MM_SNVs_unique_P_coverage[full_is_MM_SNVs_unique_P_coverage$position_coverage > 0, ]$con_freq = 1
full_is_MM_SNVs_unique_P_coverage[full_is_MM_SNVs_unique_P_coverage$position_coverage > 0, ]$var_freq = 1


# full_is_P_SNVs_unique_MM_coverage
full_is_P_SNVs_unique_MM_coverage$con_base = full_is_P_SNVs_unique_MM_coverage$ref_base
full_is_P_SNVs_unique_MM_coverage$var_base = full_is_P_SNVs_unique_MM_coverage$ref_base
full_is_P_SNVs_unique_MM_coverage$allele_count = 1
full_is_P_SNVs_unique_MM_coverage$class = "SNV"
for (base in c("A","C","T","G")){ # reset coverage
  full_is_P_SNVs_unique_MM_coverage[[base]] = 0
}

for (row in 1:nrow(full_is_P_SNVs_unique_MM_coverage)){
  full_is_P_SNVs_unique_MM_coverage[[full_is_P_SNVs_unique_MM_coverage$ref_base[row]]][row] = full_is_P_SNVs_unique_MM_coverage$position_coverage[row] # add correct coverage to allele of reference
}

full_is_P_SNVs_unique_MM_coverage[full_is_P_SNVs_unique_MM_coverage$mutation_type == "S" | full_is_P_SNVs_unique_MM_coverage$mutation_type == "N", ]$mutation_type = "S"
full_is_P_SNVs_unique_MM_coverage[full_is_P_SNVs_unique_MM_coverage$mutation_type == "S", ]$mutation = "S:"
full_is_P_SNVs_unique_MM_coverage[full_is_P_SNVs_unique_MM_coverage$mutation_type == "S", ]$cryptic = "False" # all S are never cryptic

full_is_P_SNVs_unique_MM_coverage$ref_freq = 0 # reset frequencies
full_is_P_SNVs_unique_MM_coverage$con_freq = 0
full_is_P_SNVs_unique_MM_coverage$var_freq = 0

full_is_P_SNVs_unique_MM_coverage[full_is_P_SNVs_unique_MM_coverage$position_coverage > 0, ]$ref_freq = 1 # if coverage is present, put frequency to 1
full_is_P_SNVs_unique_MM_coverage[full_is_P_SNVs_unique_MM_coverage$position_coverage > 0, ]$con_freq = 1
full_is_P_SNVs_unique_MM_coverage[full_is_P_SNVs_unique_MM_coverage$position_coverage > 0, ]$var_freq = 1


# 2. merge with original datasets ----
is_MM_SNVs_final = rbind(is_MM_SNVs, full_is_P_SNVs_unique_MM_coverage)
is_P_SNVs_final = rbind(is_P_SNVs, full_is_MM_SNVs_unique_P_coverage)

# 3. merge final datasets
merged_df = merge(is_P_SNVs_final, is_MM_SNVs_final, by = "snv_pos")
colnames(merged_df) # x is P, y is MM

# 4. fix annotation discrepancies

# bugchecking

# View(merged_df[which(merged_df$mutation_type.x == "" | merged_df$mutation_type.y == ""),])

# View(merged_df[which((merged_df$gene.y == "" & merged_df$mutation_type.y == "") | (merged_df$gene.x == "" & merged_df$mutation_type.x == "")),])

dim(merged_df[which((merged_df$gene.y == "" & merged_df$mutation_type.y == "") | (merged_df$gene.x == "" & merged_df$mutation_type.x == "")),])
dim(is_MM_SNVs[which((is_MM_SNVs$gene == "" & is_MM_SNVs$mutation_type== "")),])

merged_df$mutation_type.x[which(merged_df$gene.x == "" & merged_df$mutation_type.x == "")] = "I"
merged_df$mutation_type.y[which(merged_df$gene.y == "" & merged_df$mutation_type.y == "")] = "I"

dim(merged_df[which((merged_df$gene.y == "" & merged_df$mutation_type.y == "") | (merged_df$gene.x == "" & merged_df$mutation_type.x == "")),])
dim(merged_df[which((merged_df$mutation_type.x == "" | merged_df$mutation_type.y == "")),])

# statistics
dim(is_MM_SNVs)[1]
dim(is_MM_SNVs_unique)[1]
dim(is_MM_SNVs_final)[1]

dim(is_P_SNVs)[1]
dim(is_P_SNVs_unique)[1]
dim(is_P_SNVs_final)[1]

dim(merged_df)[1]
dim(merged_df)[1] - dim(is_MM_SNVs)[1] # in P, not in MM
dim(merged_df)[1] - dim(is_P_SNVs)[1] # in MM, not in P

# plot(merged_df$A.x[1:1000], merged_df$A.y[1:1000])
# plot(merged_df$allele_count.x, merged_df$allele_count.y)


# differences in allele numbers : supplementary information?
length(which(abs(merged_df$allele_count.x - merged_df$allele_count.y) > 0)) # unequal allele numbers
length(which(merged_df$allele_count.x == 1 | merged_df$allele_count.y == 1)) # 1 allele
length(which(merged_df$allele_count.x == 1 & merged_df$allele_count.y == 1)) # both 1 allele
sum(pop_snps$population_SNPs) #  number of population SNPs
length(which(merged_df$allele_count.x == 1 & merged_df$allele_count.y > 1)) # P 1 allele, MM > 1
length(which(merged_df$allele_count.x > 1 & merged_df$allele_count.y == 1)) # MM 1 allele, P > 1
length(which(merged_df$allele_count.x == 2 & merged_df$allele_count.y == 2)) # both 2 allele
length(which(merged_df$allele_count.x > 2 & merged_df$allele_count.y > 2)) # both multiallelic
length(which(merged_df$allele_count.x > 2 & merged_df$allele_count.y == 2)) # MM 2 alleles, P > 3
length(which(merged_df$allele_count.x == 2 & merged_df$allele_count.y > 2)) # P 2 alleles, MM > 3

# plot P vs MM allele count
# hist(merged_df$allele_count.x - merged_df$allele_count.y, freq=TRUE)
# plot(merged_df$ref_freq.x, merged_df$con_freq.x)
# plot(merged_df$ref_freq.x, merged_df$ref_freq.y, xlim=c(0,1), ylim=c(0,1))





# SNP selection -----------------------------------------------------------

# calculate differences in ref_freq, because ref is always the same
merged_df$ref_freq_diff = merged_df$ref_freq.y - merged_df$ref_freq.x

# 95% CI
ref_freq_hi_2sd = mean(merged_df$ref_freq_diff) + 1.96 * sd(merged_df$ref_freq_diff)
ref_freq_lo_2sd = mean(merged_df$ref_freq_diff) - 1.96 * sd(merged_df$ref_freq_diff)

# plot(density(merged_df$ref_freq_diff), xlab = "Reference allele frequency difference (MM - P)", main = "Reference frequency difference distribution")
abline(v = ref_freq_hi_2sd,col="red")
abline(v = ref_freq_lo_2sd,col="red")
text(c(ref_freq_lo_2sd-0.2, ref_freq_hi_2sd+0.2), y = c(0.3,0.3), labels = as.character(c(format(round(ref_freq_lo_2sd, 4), nsmall = 4),format(round(ref_freq_hi_2sd, 4), nsmall = 4))))

# find which snv pos have ref_freq > ref_freq_hi_2sd or ref_freq < ref_freq_lo_2sd
selected_df = merged_df[((merged_df$ref_freq_diff >= ref_freq_hi_2sd) | (merged_df$ref_freq_diff <= ref_freq_lo_2sd)),]
dim(selected_df)

# check if there are SNVs which are in 1 but not in the other
# raw_pop_snp_diff = merged_df[which(merged_df$scaffold.x != merged_df$scaffold.y),]
# raw_pop_snp_diff
# 
# pop_snp_diff = selected_df[which(selected_df$scaffold.x != selected_df$scaffold.y),]
# pop_snp_diff # there are none



# add taxonomy to selected SNVs
selected_df_taxonomy = merge(selected_df, taxonomy, by.x="scaffold.x", by.y="contig_id") # use this for further analysis
merged_df_taxonomy = merge(merged_df, taxonomy, by.x="scaffold.x", by.y="contig_id") # this can be used for comparison later on: taxonomic distribution


# ----------------------------------------------------------------------------------------------------------------------------
# for each SNP, perform Fisher's exact test

# NEW FISHER'S EXACT TEST WITH 2X4 CONTINGENCY TABLE
fisher_p = rep(0, length(selected_df_taxonomy$scaffold.x))
for (row in 1:nrow(selected_df_taxonomy)){
  allele_counts_x = c(selected_df_taxonomy$A.x[row], selected_df_taxonomy$C.x[row], selected_df_taxonomy$G.x[row], selected_df_taxonomy$T.x[row])
  allele_counts_y = c(selected_df_taxonomy$A.y[row], selected_df_taxonomy$C.y[row], selected_df_taxonomy$G.y[row], selected_df_taxonomy$T.y[row])
  fisher_test_result = fisher.test(rbind(allele_counts_x,allele_counts_y))
  fisher_p[row] = fisher_test_result$p.value
}

selected_df_taxonomy = cbind(selected_df_taxonomy, fisher_p)



# multiple testing correction: Benjamini Hochberg
selected_df_taxonomy = selected_df_taxonomy[order(selected_df_taxonomy$fisher_p),]

alpha = 0.01
threshold_list = rep(NA, dim(selected_df_taxonomy)[1])
for (k in 1:nrow(selected_df_taxonomy)){
  threshold_list[k] = k / nrow(selected_df_taxonomy) * alpha
}
selected_df_taxonomy = cbind(selected_df_taxonomy, threshold_list)

accept_list = selected_df_taxonomy$fisher_p <= selected_df_taxonomy$threshold_list

fisher_selected_df = selected_df_taxonomy[1:max(which(accept_list)),]

# statistics
dim(merged_df)[1] # total number of SNPs in raw data
dim(selected_df_taxonomy)[1] # SNPs after taking the 95% CI
dim(fisher_selected_df)[1] # SNPs after fisher's exact test and Benjamini-Hochberg multiple testing correction


##### Writing BED file --------------------
# This BED file will be used to extract the SNP allele coverage from the RILs using a python script
# The QTL analysis is performed in another R script
streptomyces_bed = subset(fisher_selected_df[which(fisher_selected_df$genus == "Streptomyces"),], select = c("scaffold.x", "position.x"))
dim(streptomyces_bed)
colnames(streptomyces_bed) = c("contig_id", "position")

# write.table(streptomyces_bed, file = "./bed_files/streptomyces_selected_snps_bed_210114.tsv", sep = "\t", col.names = F, row.names = F, quote = F)

full_bed = subset(fisher_selected_df, select = c("scaffold.x", "position.x"))
dim(full_bed)
colnames(full_bed) = c("contig_id", "position")
full_bed = full_bed[order(full_bed$contig_id, full_bed$position),]

# write.table(full_bed, file = "./bed_files/all_selected_snps_MM_P_bed_ECvdW_210114.tsv", sep = "\t", col.names = F, row.names = F, quote = F)
write.table(full_bed, file = "./bed_files/all_selected_snps_MM_P_bed_ECvdW_0based_bed_210617.tsv", sep = "\t", col.names = F, row.names = F, quote = F)

# Write BED file with ALL SNPs ------------------------------------

cellvibrio_bed = subset(merged_df[which(merged_df$scaffold.x %in% contigs_cellvibrio$contig),], select = c("scaffold.x", "position.x"))
streptomyces_bed = subset(merged_df[which(merged_df$scaffold.x %in% contigs_streptomyces$contig),], select = c("scaffold.x", "position.x"))
cellvibrio_streptomyces_bed = rbind(cellvibrio_bed, streptomyces_bed)

cellvibrio_bed = full_bed[order(cellvibrio_bed$scaffold.x, cellvibrio_bed$position.x),]
streptomyces_bed = full_bed[order(streptomyces_bed$scaffold.x, streptomyces_bed$position.x),]
cellvibrio_streptomyces_bed = full_bed[order(cellvibrio_streptomyces_bed$scaffold.x, cellvibrio_streptomyces_bed$position.x),]

dim(cellvibrio_bed)
length(which(merged_df_taxonomy$genus == "Cellvibrio"))

dim(streptomyces_bed)
length(which(merged_df_taxonomy$genus == "Streptomyces"))


dim(cellvibrio_streptomyces_bed)
dim(cellvibrio_bed)[1] + dim(streptomyces_bed)[1]

write.table(cellvibrio_bed, file = "./bed_files/cellvibrio_all_SNPs_MM_P_ECvdW_0based_bed_210618.tsv", sep = "\t", col.names = F, row.names = F, quote = F)
write.table(streptomyces_bed, file = "./bed_files/streptomyces_all_snps_MM_P_ECvdW_0based_bed_210618.tsv", sep = "\t", col.names = F, row.names = F, quote = F)
write.table(cellvibrio_streptomyces_bed, file = "./bed_files/cellvibrio_streptomyces_all_snps_MM_P_ECvdW_0based_bed_210618.tsv", sep = "\t", col.names = F, row.names = F, quote = F)



# compare with old bed file -----------------------------------------
full_bed_old = read.table("./bed_files/all_selected_snps_MM_P_bed_ECvdW_0based_bed_210224.tsv", sep = "\t")
colnames(full_bed_old) = c("contig_id", "position")
full_bed_old = full_bed_old[order(full_bed_old$contig_id, full_bed_old$position),]
write.table(full_bed_old, file = "./bed_files/all_selected_snps_MM_P_bed_ECvdW_0based_bed_210224.tsv", sep = "\t", col.names = F, row.names = F, quote = F)

full_bed_old$snv_pos = paste(full_bed_old$contig_id, full_bed_old$position, sep = "_")
full_bed$snv_pos = paste(full_bed$contig_id, full_bed$position, sep = "_")


bed_union = union(full_bed$snv_pos, full_bed_old$snv_pos) # both sets combined
bed_snv_old = full_bed_old$snv_pos %in% bed_union # in old
bed_snv_new = full_bed$snv_pos %in% bed_union # in new

length(bed_snv_old)
length(bed_snv_new)

bed_snv_old_unique = full_bed_old$snv_pos[!(full_bed_old$snv_pos %in% full_bed$snv_pos)]
length(bed_snv_old_unique)

bed_snv_new_unique = full_bed$snv_pos[!(full_bed$snv_pos %in% full_bed_old$snv_pos)]
length(bed_snv_new_unique)



# Add annotation metadata ---------------------------------------
# to SNPs that are selected by Fisher's exact test
# fix the indexing issues: GFF is 1-based, inStrain is 0-based
# make GFF file 0-based (inclusive) and adapt the addition of the annotation and the pos in gene
### Analysis of genes

short_contig_id_list = rep(dim(fisher_selected_df)[1])
for (row in 1:nrow(fisher_selected_df)){
  short_contig_id_list[row] = unlist(strsplit(fisher_selected_df$scaffold.x[row], "_length_"))[1]
}

f_df_genes_annotation = cbind(fisher_selected_df, short_contig_id_list) # create new dataframe to protect the old one

f_df_genes_annotation$feature = rep("", nrow(f_df_genes_annotation))
f_df_genes_annotation$strand = rep("", nrow(f_df_genes_annotation))
f_df_genes_annotation$start = rep("", nrow(f_df_genes_annotation))
f_df_genes_annotation$end = rep("", nrow(f_df_genes_annotation))
f_df_genes_annotation$pos_in_gene = rep("", nrow(f_df_genes_annotation))
f_df_genes_annotation$prokka_gene_id = rep("", nrow(f_df_genes_annotation))
f_df_genes_annotation$annotation = rep("", nrow(f_df_genes_annotation))
f_df_genes_annotation$source = rep("", nrow(f_df_genes_annotation))
f_df_genes_annotation$attributes = rep("", nrow(f_df_genes_annotation))


extra_row_counter = nrow(f_df_genes_annotation)
max_rows = nrow(f_df_genes_annotation)
for (row in 1:max_rows){
  temp_contigs_list = gff_file[which(gff_file$contig_id_short == f_df_genes_annotation$short_contig_id_list[row]),]
  
  if (nrow(temp_contigs_list) > 0 &   # check if there are contigs in GFF file
      any(temp_contigs_list$start <= f_df_genes_annotation[row,]$position.x &  # check if SNP pos is within a region; GFF file is 0-based
          temp_contigs_list$end >= f_df_genes_annotation[row,]$position.x)){
    
    temp_pos_contigs_list = temp_contigs_list[which(temp_contigs_list$start <= f_df_genes_annotation[row,]$position.x &
                                                      temp_contigs_list$end >= f_df_genes_annotation[row,]$position.x),]
    if (nrow(temp_pos_contigs_list)>1){print(paste("ROW:", row, "HAS", nrow(temp_pos_contigs_list), "ANNOTATIONS"))}
    for (ann_i in 1:nrow(temp_pos_contigs_list)){
      if (ann_i > 1){
        extra_row_counter = extra_row_counter + 1
        print(extra_row_counter)
        f_df_genes_annotation[extra_row_counter,] = f_df_genes_annotation[row,]
        f_df_genes_annotation$attributes[extra_row_counter] = temp_pos_contigs_list$attributes[ann_i]
        
        temp_annotation_list = unlist(strsplit(f_df_genes_annotation$attributes[extra_row_counter], ";"))
        f_df_genes_annotation$annotation[extra_row_counter] = unlist(strsplit(temp_annotation_list[length(temp_annotation_list)], "product="))[2]
        f_df_genes_annotation$prokka_gene_id[extra_row_counter] = unlist(strsplit(temp_annotation_list[1], "ID="))[2]
        
        f_df_genes_annotation$feature[extra_row_counter] =    temp_pos_contigs_list$feature[ann_i]
        f_df_genes_annotation$strand[extra_row_counter] =     temp_pos_contigs_list$strand[ann_i]
        f_df_genes_annotation$start[extra_row_counter] =      as.numeric(temp_pos_contigs_list$start[ann_i])
        f_df_genes_annotation$end[extra_row_counter] =        as.numeric(temp_pos_contigs_list$end[ann_i])
        f_df_genes_annotation$source[extra_row_counter] =     temp_pos_contigs_list$source[ann_i]
        
        # incorporate SNV pos in gene
        # positive strand is fine, negative should be reversed!
        if (f_df_genes_annotation$strand[extra_row_counter] == "+"){
          f_df_genes_annotation$pos_in_gene[extra_row_counter] = -as.numeric(f_df_genes_annotation$start[extra_row_counter]) + f_df_genes_annotation$position.x[extra_row_counter] # 0-based: pos_in_gene = pos - start
        } else if (f_df_genes_annotation$strand[extra_row_counter] == "-"){
          f_df_genes_annotation$pos_in_gene[extra_row_counter] = as.numeric(f_df_genes_annotation$end[extra_row_counter]) - f_df_genes_annotation$position.x[extra_row_counter] # 0-based: pos_in_gene = end - pos
        }
        
      } else{
        f_df_genes_annotation$attributes[row] = temp_pos_contigs_list$attributes[ann_i]
        
        temp_annotation_list = unlist(strsplit(f_df_genes_annotation$attributes[row], ";"))
        f_df_genes_annotation$annotation[row] = unlist(strsplit(temp_annotation_list[length(temp_annotation_list)], "product="))[2]
        f_df_genes_annotation$prokka_gene_id[row] = unlist(strsplit(temp_annotation_list[1], "ID="))[2]
        
        f_df_genes_annotation$feature[row] =    temp_pos_contigs_list$feature[ann_i]
        f_df_genes_annotation$strand[row] =     temp_pos_contigs_list$strand[ann_i]
        f_df_genes_annotation$start[row] =      as.numeric(temp_pos_contigs_list$start[ann_i])
        f_df_genes_annotation$end[row] =        as.numeric(temp_pos_contigs_list$end[ann_i])
        f_df_genes_annotation$source[row] =     temp_pos_contigs_list$source[ann_i]
        
        if (f_df_genes_annotation$strand[row] == "+"){
          f_df_genes_annotation$pos_in_gene[row] = -as.numeric(f_df_genes_annotation$start[row]) + f_df_genes_annotation$position.x[row] # the start nt also counts, so add 1
        } else if (f_df_genes_annotation$strand[row] == "-"){
          f_df_genes_annotation$pos_in_gene[row] = as.numeric(f_df_genes_annotation$end[row]) - f_df_genes_annotation$position.x[row] # end does not count, but still add 1 as you include the new start position
        }
      }
    }
  }  
}
f_df_genes_annotation$mutation.x[which(f_df_genes_annotation$mutation.x == "S:")] = paste("S", f_df_genes_annotation$pos_in_gene[which(f_df_genes_annotation$mutation.x == "S:")], sep = ":")
f_df_genes_annotation$mutation.y[which(f_df_genes_annotation$mutation.y == "S:")] = paste("S", f_df_genes_annotation$pos_in_gene[which(f_df_genes_annotation$mutation.y == "S:")], sep = ":")



# Statistics of selected SNVs ------------------------------------------------------------------------

nrow(fisher_selected_df) # total number of selected SNPs
nrow(f_df_genes_annotation) # total number of SNPs, and SNPs being duplicated if they have multiple annotations
nrow(f_df_genes_annotation) - nrow(fisher_selected_df) # number of duplicated SNPs due to multiple annotations
length(unique(f_df_genes_annotation$snv_pos)) # this is fine

length(unique(f_df_genes_annotation$scaffold.x)) # number of contigs in the selected SNP list
length(which(f_df_genes_annotation$attributes != "")) # number of annotations: can be multiple annotations per SNP
length(which(f_df_genes_annotation$attributes == "")) # number of SNPs without annotations

# Prokka
length(unique(gff_file$attributes)) # number of annotations in GFF file
(length(unique(f_df_genes_annotation$attributes))-1) / length(which(f_df_genes_annotation$attributes != "")) * 100 # % of unique annotations in SNP list
length(unique(f_df_genes_annotation$prokka_gene_id)) - 1# number of genes with unique identifier
length(unique(f_df_genes_annotation$annotation)) - 1 # number of unique annotations

# Prokka FNA
length(unique(f_df_genes_annotation$gene.x)) # Prodigal: number of predicted genes: expected to be less than Prokka
length(unique(f_df_genes_annotation$gene.y)) # Prodigal: number of predicted genes: expected to be less than Prokka
length(which(f_df_genes_annotation$gene.x == "")) # Number of SNPs for P without Prodigal genes
length(which(f_df_genes_annotation$gene.y == "")) # Number of SNPs for MM without Prodigal genes
# MM has more annotations: weird that there is an inconsistentency between MM and P

# S/N/I statistics
length(which(fisher_selected_df$mutation_type.x == "S" & fisher_selected_df$mutation_type.y == "S"))
length(which(fisher_selected_df$mutation_type.x == "N" | fisher_selected_df$mutation_type.y == "N"))
length(which(fisher_selected_df$mutation_type.x == "I" & fisher_selected_df$mutation_type.y == "I")) # I in both datasets: used for report
length(which(fisher_selected_df$mutation_type.x == "I" & fisher_selected_df$mutation_type.y == "S"))
length(which(fisher_selected_df$mutation_type.x == "S" & fisher_selected_df$mutation_type.y == "I"))
length(which(fisher_selected_df$mutation_type.x == "M" | fisher_selected_df$mutation_type.y == "M"))
# View(fisher_selected_df[which(fisher_selected_df$mutation_type.x == "M" | fisher_selected_df$mutation_type.y == "M"),])

# check intergenic regions
length(which(f_df_genes_annotation$mutation_type.x == "I"))
length(which(f_df_genes_annotation$mutation_type.y == "I"))
length(which(f_df_genes_annotation$mutation_type.x == "I" & f_df_genes_annotation$mutation_type.y == "S"))
length(which(f_df_genes_annotation$mutation_type.x == "I" & f_df_genes_annotation$mutation_type.y == "N"))
length(which(f_df_genes_annotation$mutation_type.y == "I" & f_df_genes_annotation$mutation_type.x == "S"))
length(which(f_df_genes_annotation$mutation_type.y == "I" & f_df_genes_annotation$mutation_type.x == "N"))


length(which(f_df_genes_annotation$mutation_type.x == "I" & f_df_genes_annotation$mutation_type.y == "I")) # I in both datasets: used for report
length(which(f_df_genes_annotation$mutation_type.x == "I" | f_df_genes_annotation$mutation_type.y == "I")) # I in either dataset
length(which(f_df_genes_annotation$mutation_type.x != "I" & f_df_genes_annotation$mutation_type.y == "I")) # I in MM not in P
length(which(f_df_genes_annotation$mutation_type.y != "I" & f_df_genes_annotation$mutation_type.x == "I")) # I in P not in MM

# View(f_df_genes_annotation[which(f_df_genes_annotation$mutation_type.y == "" | f_df_genes_annotation$mutation_type.x == ""),])


# compare Intergenic I and annotations

length(which(f_df_genes_annotation$mutation_type.x == "I" & f_df_genes_annotation$mutation_type.y == "I" & f_df_genes_annotation$attributes == "")) 
length(which(f_df_genes_annotation$mutation_type.x == "I" & f_df_genes_annotation$mutation_type.y == "I" & f_df_genes_annotation$attributes != ""))

# View(f_df_genes_annotation[which((f_df_genes_annotation$mutation_type.x == "I" | f_df_genes_annotation$mutation_type.y == "I") & f_df_genes_annotation$attributes != ""),])

View(f_df_genes_annotation[which(f_df_genes_annotation$position.x == f_df_genes_annotation$start),])

# check overlap of prokka GFF annotations and the inStrain N/S that should be in genes.
dim(f_df_genes_annotation[which((f_df_genes_annotation$mutation_type.x != "I" & f_df_genes_annotation$mutation_type.y != "I") & f_df_genes_annotation$attributes == ""),])

length(which((f_df_genes_annotation$mutation_type.x == "I" | f_df_genes_annotation$mutation_type.y == "I") & f_df_genes_annotation$attributes != ""))


# Writing Fisher selected DF with annotations -------------

write.table(f_df_genes_annotation, file = "./snp_data/SNV_data_MM_P_selected_Fisher_annotations_ECvdW_210607.tsv", sep = "\t", col.names = T, row.names = F, quote = F)





##### Analysis of selected SNPs ----------------------------------------------------------------------
# get only the nonsynonymous SNPs:
# nonsynonymous_snps_df = data.frame(f_df_genes_annotation[which(f_df_genes_annotation$mutation_type.x == "N" | f_df_genes_annotation$mutation_type.y == "N"),])
# 
# length(unique(nonsynonymous_snps_df$attributes))-1 # number of genes with N mutations in either set
# sort(table(nonsynonymous_snps_df$attributes), decreasing = T)[2:11]
# par(mfrow=c(1,1))
# table(nonsynonymous_snps_df$attributes) %>% sort() %>% plot()




# Analysis of annotations --------------------------------
length(unique(f_df_genes_annotation$attributes)) -1 # number of annotations
length(unique(f_df_genes_annotation$prokka_gene_id)) - 1# number of annotation - contig couples -- approximate number of annotated genes
dim(f_df_genes_annotation)

for (level in c("prokka_gene_id","scaffold.x","annotation")){
  unique_genes = data.frame(unique(f_df_genes_annotation[[level]]))
  dim(unique_genes)
  colnames(unique_genes) = level
  for (row in 1:nrow(unique_genes)){
    unique_genes$SNP_count[row] = length(which(unique_genes[[level]][row] == f_df_genes_annotation[[level]]))
    
    unique_genes$N_either_count[row] = length(which(f_df_genes_annotation$mutation_type.x[which(unique_genes[[level]][row] == f_df_genes_annotation[[level]])] == "N" | 
                                                      f_df_genes_annotation$mutation_type.y[which(unique_genes[[level]][row] == f_df_genes_annotation[[level]])] == "N"))
    unique_genes$S_both_count[row] =   length(which(f_df_genes_annotation$mutation_type.x[which(unique_genes[[level]][row] == f_df_genes_annotation[[level]])] == "S" & 
                                                      f_df_genes_annotation$mutation_type.y[which(unique_genes[[level]][row] == f_df_genes_annotation[[level]])] == "S"))
    unique_genes$I_both_count[row] =   length(which(f_df_genes_annotation$mutation_type.x[which(unique_genes[[level]][row] == f_df_genes_annotation[[level]])] == "I" & 
                                                      f_df_genes_annotation$mutation_type.y[which(unique_genes[[level]][row] == f_df_genes_annotation[[level]])] == "I"))
    unique_genes$discrepancy[row] = unique_genes$SNP_count[row] - sum(unique_genes$N_either_count[row], unique_genes$S_both_count[row], unique_genes$I_both_count[row])
    
    unique_genes$N_either_perc[row] = unique_genes$N_either_count[row] / unique_genes$SNP_count[row] * 100
    unique_genes$S_both_perc[row] = unique_genes$S_both_perc[row] / unique_genes$SNP_count[row] * 100
    unique_genes$N_S_ratio_total[row] = unique_genes$N_either_perc[row] / unique_genes$S_both_perc[row]
    unique_genes$I_both_perc[row] = unique_genes$I_both_count[row] / unique_genes$SNP_count[row] * 100
    
    unique_genes$P_N_perc[row] = length(which(f_df_genes_annotation$mutation_type.x[which(unique_genes[[level]][row] == f_df_genes_annotation[[level]])] == "N")) / unique_genes$SNP_count[row] * 100
    unique_genes$P_S_perc[row] = length(which(f_df_genes_annotation$mutation_type.x[which(unique_genes[[level]][row] == f_df_genes_annotation[[level]])] == "S")) / unique_genes$SNP_count[row] * 100
    unique_genes$P_N_S_ratio[row] = unique_genes$P_N_perc[row] / unique_genes$P_S_perc[row]
    unique_genes$MM_N_perc[row] = length(which(f_df_genes_annotation$mutation_type.y[which(unique_genes[[level]][row] == f_df_genes_annotation[[level]])] == "N")) / unique_genes$SNP_count[row] * 100
    unique_genes$MM_S_perc[row] = length(which(f_df_genes_annotation$mutation_type.y[which(unique_genes[[level]][row] == f_df_genes_annotation[[level]])] == "S")) / unique_genes$SNP_count[row] * 100
    unique_genes$MM_N_S_ratio[row] = unique_genes$MM_N_perc[row] / unique_genes$MM_S_perc[row]
    
    contig_list = unique(f_df_genes_annotation$scaffold.x[which(f_df_genes_annotation[[level]] == unique_genes[[level]][row])])
    unique_genes$contig_count[row] = length(contig_list)
    
    unique_genes$contig_list[row] = paste(contig_list, collapse = ";")
    
    # total gene count is the number of attribute_pos per annotation: contig / attribute / attribute_pos
    gene_list = unique(f_df_genes_annotation$prokka_gene_id[which(f_df_genes_annotation[[level]] == unique_genes[[level]][row])])
    unique_genes$total_gene_count[row] = length(gene_list)
    unique_genes$prokka_gene_id_list[row] = paste(gene_list, collapse = ";")
    
    if (any(level == c("annotation", "prokka_gene_id"))){ # position in gene list
      unique_genes$pos_in_gene[row] = paste(f_df_genes_annotation$pos_in_gene[which(unique_genes[[level]][row] == f_df_genes_annotation[[level]])], collapse = ";")
    }
    
    
    for (tax in c("class","order","family","genus","species")){
      tax_list = f_df_genes_annotation[[tax]][which(f_df_genes_annotation[[level]] == unique_genes[[level]][row])]
      unique_genes[[paste(tax, "_count", sep="")]][row] = length(unique(tax_list))
      unique_genes[[paste(tax, "_list", sep="")]][row] = paste(unique(tax_list), collapse = ";")
    }
  }
  unique_genes = unique_genes[order(unique_genes$N_either_count, decreasing=TRUE),]
  write.csv(unique_genes, file = paste("./analysis_output/filtered_SNP_gene_analysis_210607_",level, ".csv", sep=""), row.names = F)
}






# Contig analysis -------------------------
unique_contigs = data.frame(unique(f_df_genes_annotation$scaffold.x))
colnames(unique_contigs) = c("contig_id")
dim(unique_contigs)[1] # total number of contigs containing the selected SNPs

unique_contigs = merge(unique_contigs, taxonomy)

for (row in 1:nrow(unique_contigs)){
  # SNPs per contig
  unique_contigs$SNP_count[row] = length(which(unique_contigs$contig_id[row] == f_df_genes_annotation$scaffold.x))
  # contig size
  unique_contigs$contig_length[row] = as.numeric(strsplit(strsplit(unique_contigs$contig_id[row], "_length_")[[1]][2], "_cov_")[[1]][1])
}
par(mfrow=c(1,1))
hist(unique_contigs$SNP_count, breaks=seq(min(unique_contigs$SNP_count), max(unique_contigs$SNP_count),1), xlab = "SNP count", main="SNP frequency per contig")

# plot contig length vs SNP counts
par(mfrow=c(1,2))
plot(unique_contigs$contig_length, unique_contigs$SNP_count, xlab="contig size", ylab="SNP count")
plot(log(unique_contigs$contig_length), log(unique_contigs$SNP_count), xlab="log(contig size)", ylab="log(SNP count)")

mean(unique_contigs$SNP_count)




# QTL loci analysis -------------------------------------------------------------------

# load in QTL data: linked to bed file of all selected SNPs
found_peaks_old = read.table("QTL_output/found_peaks_210121.txt", header=TRUE, sep = "\t")
found_peaks = read.table("QTL_output/found_peaks_210618.txt", header=TRUE, sep = "\t")
# parse out SNV_pos column to extract annotations from the 'f_df_genes_annotation' table

length(which(found_peaks$allele == "modern"))
length(which(found_peaks$allele == "wild"))

extra_row_counter = nrow(found_peaks)
for (row in 1:nrow(found_peaks)){
  found_peaks$snv_pos[row] = paste(unlist(strsplit(found_peaks$lodcolumn[row], "_"))[1:8], collapse = "_")
  found_peaks$short_contig_id[row] = short_contig_id = unlist(strsplit(found_peaks$snv_pos[row], "_length_"))[1]
  found_peaks$base[row] = paste(unlist(strsplit(found_peaks$lodcolumn[row], "_"))[9], collapse = "_")
  found_peaks$ref_base[row] = f_df_genes_annotation$ref_base.x[which(found_peaks$snv_pos[row] == f_df_genes_annotation$snv_pos)][1]
  found_peaks$pos_in_gene[row] = "" # initialise for correct order in table
  
  # annotations
  
  temp_peak_annotation_df = f_df_genes_annotation[which(found_peaks$snv_pos[row] == f_df_genes_annotation$snv_pos),]
  
  if (nrow(temp_peak_annotation_df)>1){print(paste("ROW:", row, "HAS", nrow(temp_peak_annotation_df), "ANNOTATIONS"))}
  for (k in 1:nrow(temp_peak_annotation_df)){
    if (k > 1){
      extra_row_counter = extra_row_counter + 1
      print(extra_row_counter)
      
      found_peaks[extra_row_counter,] = found_peaks[row,]
      
      found_peaks$prokka_gene_id[extra_row_counter] = temp_peak_annotation_df$prokka_gene_id[k]
      found_peaks$annotation[extra_row_counter] = temp_peak_annotation_df$annotation[k]
      found_peaks$start[extra_row_counter] =      temp_peak_annotation_df$start[k]
      found_peaks$end[extra_row_counter] =        temp_peak_annotation_df$end[k]
      found_peaks$source[extra_row_counter] =     temp_peak_annotation_df$source[k]
      found_peaks$attributes[extra_row_counter] = temp_peak_annotation_df$attributes[k]
      found_peaks$pos_in_gene[extra_row_counter] =    temp_peak_annotation_df$pos_in_gene[k]
      
    } else if (k == 1){
      found_peaks$prokka_gene_id[row] = temp_peak_annotation_df$prokka_gene_id[k]
      found_peaks$annotation[row] =     temp_peak_annotation_df$annotation[k]
      found_peaks$start[row] =          temp_peak_annotation_df$start[k]
      found_peaks$end[row] =            temp_peak_annotation_df$end[k]
      found_peaks$source[row] =         temp_peak_annotation_df$source[k]
      found_peaks$attributes[row] =     temp_peak_annotation_df$attributes[k]
      found_peaks$pos_in_gene[row] =    temp_peak_annotation_df$pos_in_gene[k]
    }
  }
  
  found_peaks$upstream_gene_id[row] = ""
  found_peaks$upstream_distance[row] = ""
  found_peaks$upstream_gene_annotation[row] = ""
  found_peaks$downstream_gene_id[row] = ""
  found_peaks$downstream_distance[row] = ""
  found_peaks$downstream_gene_annotation[row] = ""
  
  # for intergenic SNPs:find downstream regions of SNP: both + and - in GFF file directly adjacent
  if (found_peaks$prokka_gene_id[row] == ""){
    pos = as.numeric(unlist(strsplit(found_peaks$snv_pos[row], "_pos_"))[2])
    gff_gene_list = gff_file[which(gff_file$contig_id_short == found_peaks$short_contig_id[row]),]
    
    # Downstream and Upstream genes
    if (dim(gff_gene_list)[1] > 0 & (any(gff_gene_list$start > pos) | any(gff_gene_list$end < pos))){ # up-and downstream genes are outside the GFF start-end region, as region is inclusive
      
      gff_neg_strand_genes = gff_gene_list[which(gff_gene_list$strand == "-"),]
      gff_pos_strand_genes = gff_gene_list[which(gff_gene_list$strand == "+"),]
      
      # for downstream gene: find which genes have start > pos and pick the first
      # INCORPORATE + sign
      if (any(gff_pos_strand_genes$start > pos) & nrow(gff_pos_strand_genes)>0){
        
        # add the annotation
        downstream_gene_df = gff_pos_strand_genes[which(gff_pos_strand_genes$end == min(gff_pos_strand_genes$end[which(gff_pos_strand_genes$start > pos)])),]
        found_peaks$downstream_gene_id[row] = unlist(strsplit(unlist(strsplit(downstream_gene_df$attributes, ";"))[1], "ID="))[2]
        found_peaks$downstream_gene_annotation[row] = unlist(strsplit(unlist(strsplit(downstream_gene_df$attributes, ";"))[length(unlist(strsplit(downstream_gene_df$attributes, ";")))], "product="))[2]
        found_peaks$downstream_distance[row] = downstream_gene_df$start - pos
      }
      
      # for upstream gene: find which genes have end < pos and pick the last
      # INCORPORATE - SIGN
      if (any(gff_neg_strand_genes$end < pos) & nrow(gff_neg_strand_genes)>0){
        
        # add the annotation
        upstream_gene_df = gff_neg_strand_genes[which(gff_neg_strand_genes$end == max(gff_neg_strand_genes$end[which(gff_neg_strand_genes$end < pos)])),]
        found_peaks$upstream_gene_id[row] = unlist(strsplit(unlist(strsplit(upstream_gene_df$attributes, ";"))[1], "ID="))[2]
        found_peaks$upstream_gene_annotation[row] = unlist(strsplit(unlist(strsplit(upstream_gene_df$attributes, ";"))[length(unlist(strsplit(upstream_gene_df$attributes, ";")))], "product="))[2]
        found_peaks$upstream_distance[row] = pos - upstream_gene_df$end
      }
    }
  }
}
# mutation type
for (row in 1:nrow(found_peaks)){
  mut_P = f_df_genes_annotation$mutation_type.x[which(f_df_genes_annotation$snv_pos == found_peaks$snv_pos[row])] 
  mut_MM = f_df_genes_annotation$mutation_type.y[which(f_df_genes_annotation$snv_pos == found_peaks$snv_pos[row])]
  
  if (length(mut_P) > 1){
    print(row)
    print(mut_P)
    mut_P = mut_P[1] # some have 2 annotations, but these are always the same
    mut_MM = mut_MM[1]
  }
  
  if (mut_P == "N" | mut_MM == "N"){
    found_peaks$mutation_type[row] = "N"
  } else if (mut_P == "S" & mut_MM == "S"){
    found_peaks$mutation_type[row] = "S"
  } else if (mut_P == "I" & mut_MM == "I"){
    found_peaks$mutation_type[row] = "I"
  } else if (mut_P == "M" & mut_MM == "M"){
    found_peaks$mutation_type[row] = "M"
  } else {
    found_peaks$mutation_type[row] = "D"
  }
}
# taxonomy
for (row in 1:nrow(found_peaks)){
  found_peaks$kingdom[row] = f_df_genes_annotation$kingdom[which(found_peaks$snv_pos[row] == f_df_genes_annotation$snv_pos)][1]
  found_peaks$phylum[row] = f_df_genes_annotation$phylum[which(found_peaks$snv_pos[row] == f_df_genes_annotation$snv_pos)][1]
  found_peaks$class[row] = f_df_genes_annotation$class[which(found_peaks$snv_pos[row] == f_df_genes_annotation$snv_pos)][1]
  found_peaks$order[row] = f_df_genes_annotation$order[which(found_peaks$snv_pos[row] == f_df_genes_annotation$snv_pos)][1]
  found_peaks$family[row] = f_df_genes_annotation$family[which(found_peaks$snv_pos[row] == f_df_genes_annotation$snv_pos)][1]
  found_peaks$genus[row] = f_df_genes_annotation$genus[which(found_peaks$snv_pos[row] == f_df_genes_annotation$snv_pos)][1]
  found_peaks$species[row] = f_df_genes_annotation$species[which(found_peaks$snv_pos[row] == f_df_genes_annotation$snv_pos)][1]
}

dim(found_peaks)

found_peaks = found_peaks[order(found_peaks$lod, decreasing = T),]

# write found peaks file with new annotations
write.table(found_peaks, file = paste("./snp_data/QTL_found_peaks_annotations_", file_info , ".tsv", sep=""), quote = F, row.names = F, sep = "\t")

# statistics
length(which(found_peaks$annotation == "Vitamin B12 transporter BtuB"))

length(unique(found_peaks$lodcolumn[which(found_peaks$allele == "modern")]))
length(unique(found_peaks$lodcolumn[which(found_peaks$allele == "wild")]))
length(which(found_peaks$allele == "modern"))
length(which(found_peaks$allele == "wild"))
length(unique(found_peaks$prokka_gene_id[which(found_peaks$allele == "modern")])) -1
length(unique(found_peaks$prokka_gene_id[which(found_peaks$allele == "wild")])) -1
length(unique(found_peaks$annotation[which(found_peaks$allele == "modern")])) -1
length(unique(found_peaks$annotation[which(found_peaks$allele == "wild")])) -1

qtl_snv_pos_list = unique(found_peaks$snv_pos)
length(qtl_snv_pos_list) # total number of unique SNPs
qtl_N = length(which(fisher_selected_df$mutation_type.x[which(fisher_selected_df$snv_pos %in% qtl_snv_pos_list)] == "N" | 
                       fisher_selected_df$mutation_type.y[which(fisher_selected_df$snv_pos %in% qtl_snv_pos_list)] == "N"))
qtl_S = length(which(fisher_selected_df$mutation_type.x[which(fisher_selected_df$snv_pos %in% qtl_snv_pos_list)] == "S" & 
                       fisher_selected_df$mutation_type.y[which(fisher_selected_df$snv_pos %in% qtl_snv_pos_list)] == "S"))
qtl_I = length(which(fisher_selected_df$mutation_type.x[which(fisher_selected_df$snv_pos %in% qtl_snv_pos_list)] == "I" & 
                       fisher_selected_df$mutation_type.y[which(fisher_selected_df$snv_pos %in% qtl_snv_pos_list)] == "I"))
qtl_M = length(which(fisher_selected_df$mutation_type.x[which(fisher_selected_df$snv_pos %in% qtl_snv_pos_list)] == "M" & 
                       fisher_selected_df$mutation_type.y[which(fisher_selected_df$snv_pos %in% qtl_snv_pos_list)] == "M"))
qtl_D = length(qtl_snv_pos_list) - sum(qtl_N, qtl_S, qtl_I, qtl_M)


# check if bacterial loci can have multiple plant loci
nrow(found_peaks) # QTL peaks
length(unique(found_peaks$lodindex)) # unique bacterial alleles
nrow(found_peaks) - length(unique(found_peaks$lodindex)) # difference: maximum number of bacterial alleles with at least 2 plant loci


# Analysis of QTL-selected SNPs -------------------------------------------

# TODO add a table for downstream / upstream genes

for (level in c("prokka_gene_id", "short_contig_id","annotation")){
  unique_genes_qtl = data.frame(unique(found_peaks[[level]]))
  dim(unique_genes_qtl)
  colnames(unique_genes_qtl) = level
  for (row in 1:nrow(unique_genes_qtl)){
    if (unique_genes_qtl[[level]][row] == ""){unique_genes_qtl[[level]][row] = "None"} # prevent errors if annotation is absent, which occurs once per level
    
    if (level == "prokka_gene_id" | level == "short_contig_id") {
      unique_genes_qtl$annotation[row] = paste(unique(found_peaks$annotation[which(unique_genes_qtl[[level]][row] == found_peaks[[level]])]), collapse = ";")
    }
    
    unique_genes_qtl$QTL_peak_count[row] = length(which(unique_genes_qtl[[level]][row] == found_peaks[[level]]))
    unique_genes_qtl$SNP_count[row] = 0
    unique_genes_qtl$QTL_wild[row] = length(which(found_peaks$allele[which(unique_genes_qtl[[level]][row] == found_peaks[[level]])] == "wild"))
    unique_genes_qtl$QTL_modern[row] = length(which(found_peaks$allele[which(unique_genes_qtl[[level]][row] == found_peaks[[level]])] == "modern"))
    
    # QTL statistics
    qtl_snv_pos_list = unique(found_peaks$snv_pos[which(unique_genes_qtl[[level]][row] == found_peaks[[level]])])
    
    unique_genes_qtl$QTL_N_either_count[row]  = length(which(fisher_selected_df$mutation_type.x[which(fisher_selected_df$snv_pos %in% qtl_snv_pos_list)] == "N" | 
                                                               fisher_selected_df$mutation_type.y[which(fisher_selected_df$snv_pos %in% qtl_snv_pos_list)] == "N"))
    unique_genes_qtl$QTL_S_both_count[row]    = length(which(fisher_selected_df$mutation_type.x[which(fisher_selected_df$snv_pos %in% qtl_snv_pos_list)] == "S" & 
                                                               fisher_selected_df$mutation_type.y[which(fisher_selected_df$snv_pos %in% qtl_snv_pos_list)] == "S"))
    unique_genes_qtl$QTL_I_both_count[row]    = length(which(fisher_selected_df$mutation_type.x[which(fisher_selected_df$snv_pos %in% qtl_snv_pos_list)] == "I" & 
                                                               fisher_selected_df$mutation_type.y[which(fisher_selected_df$snv_pos %in% qtl_snv_pos_list)] == "I"))
    unique_genes_qtl$QTL_discrepancy[row]     = unique_genes_qtl$QTL_peak_count[row] - sum(unique_genes_qtl$QTL_N_either_count[row], unique_genes_qtl$QTL_S_both_count[row], unique_genes_qtl$QTL_I_both_count[row])
    
    unique_genes_qtl$QTL_N_either_perc[row] = unique_genes_qtl$QTL_N_either_count[row] / unique_genes_qtl$QTL_peak_count[row] * 100
    unique_genes_qtl$QTL_S_both_perc[row] = unique_genes_qtl$QTL_S_both_count[row] / unique_genes_qtl$QTL_peak_count[row] * 100
    unique_genes_qtl$QTL_N_S_ratio_total[row] = unique_genes_qtl$QTL_N_either_perc[row] / unique_genes_qtl$QTL_S_both_perc[row]
    unique_genes_qtl$QTL_I_both_perc[row] = unique_genes_qtl$QTL_I_both_count[row] / unique_genes_qtl$QTL_peak_count[row] * 100
    
    # Selected SNVs - NOT QTL SPECIFIC, just for easy statistics
    unique_genes_qtl$SNP_count[row] = length(which(unique_genes_qtl[[level]][row] == f_df_genes_annotation[[level]]))
    
    unique_genes_qtl$N_either_count[row] = length(which(f_df_genes_annotation$mutation_type.x[which(unique_genes_qtl[[level]][row] == f_df_genes_annotation[[level]])] == "N" | 
                                                          f_df_genes_annotation$mutation_type.y[which(unique_genes_qtl[[level]][row] == f_df_genes_annotation[[level]])] == "N"))
    unique_genes_qtl$S_both_count[row] =   length(which(f_df_genes_annotation$mutation_type.x[which(unique_genes_qtl[[level]][row] == f_df_genes_annotation[[level]])] == "S" & 
                                                          f_df_genes_annotation$mutation_type.y[which(unique_genes_qtl[[level]][row] == f_df_genes_annotation[[level]])] == "S"))
    unique_genes_qtl$I_both_count[row] =   length(which(f_df_genes_annotation$mutation_type.x[which(unique_genes_qtl[[level]][row] == f_df_genes_annotation[[level]])] == "I" & 
                                                          f_df_genes_annotation$mutation_type.y[which(unique_genes_qtl[[level]][row] == f_df_genes_annotation[[level]])] == "I"))
    unique_genes_qtl$discrepancy[row] = unique_genes_qtl$SNP_count[row] - sum(unique_genes_qtl$N_either_count[row], unique_genes_qtl$S_both_count[row], unique_genes_qtl$I_both_count[row])
    
    unique_genes_qtl$N_either_perc[row] = unique_genes_qtl$N_either_count[row] / unique_genes_qtl$SNP_count[row] * 100
    unique_genes_qtl$S_both_perc[row] = unique_genes_qtl$S_both_count[row] / unique_genes_qtl$SNP_count[row] * 100
    unique_genes_qtl$N_S_ratio_total[row] = unique_genes_qtl$N_either_perc[row] / unique_genes_qtl$S_both_perc[row]
    unique_genes_qtl$I_both_perc[row] = unique_genes_qtl$I_both_count[row] / unique_genes_qtl$SNP_count[row] * 100
    
    unique_genes_qtl$P_N_perc[row] = length(which(f_df_genes_annotation$mutation_type.x[which(unique_genes_qtl[[level]][row] == f_df_genes_annotation[[level]])] == "N")) / unique_genes_qtl$SNP_count[row] * 100
    unique_genes_qtl$P_S_perc[row] = length(which(f_df_genes_annotation$mutation_type.x[which(unique_genes_qtl[[level]][row] == f_df_genes_annotation[[level]])] == "S")) / unique_genes_qtl$SNP_count[row] * 100
    unique_genes_qtl$P_N_S_ratio[row] = unique_genes_qtl$P_N_perc[row] / unique_genes_qtl$P_S_perc[row]
    unique_genes_qtl$MM_N_perc[row] = length(which(f_df_genes_annotation$mutation_type.y[which(unique_genes_qtl[[level]][row] == f_df_genes_annotation[[level]])] == "N")) / unique_genes_qtl$SNP_count[row] * 100
    unique_genes_qtl$MM_S_perc[row] = length(which(f_df_genes_annotation$mutation_type.y[which(unique_genes_qtl[[level]][row] == f_df_genes_annotation[[level]])] == "S")) / unique_genes_qtl$SNP_count[row] * 100
    unique_genes_qtl$MM_N_S_ratio[row] = unique_genes_qtl$MM_N_perc[row] / unique_genes_qtl$MM_S_perc[row]
    
    # QTL-unique parameters
    plant_chr_list = found_peaks$chr[which(found_peaks[[level]] == unique_genes_qtl[[level]][row])]
    plant_pos_list = found_peaks$pos[which(found_peaks[[level]] == unique_genes_qtl[[level]][row])]
    plant_chr_pos_list = paste(plant_chr_list, "_", plant_pos_list, collapse = ";", sep = "")
    unique_genes_qtl$plant_QTL_loci[row] = plant_chr_pos_list
    
    contig_list = unique(found_peaks$short_contig_id[which(found_peaks[[level]] == unique_genes_qtl[[level]][row])])
    unique_genes_qtl$contig_count[row] = length(contig_list)
    
    unique_genes_qtl$contig_list[row] = paste(contig_list, collapse = ";")
    
    # total gene count is the number of attribute_pos per "level"
    gene_list = unique(found_peaks$prokka_gene_id[which(found_peaks[[level]] == unique_genes_qtl[[level]][row])])
    unique_genes_qtl$total_gene_count[row] = length(gene_list)
    unique_genes_qtl$prokka_gene_id_list[row] = paste(gene_list, collapse = ";")
    
    # add the position in the gene from f_df_genes_annotation and found_peaks
    unique_genes_qtl$pos_in_gene_qtl[row] = paste(found_peaks$pos_in_gene[which(found_peaks[[level]] == unique_genes_qtl[[level]][row])], collapse = ";")
    unique_genes_qtl$pos_in_gene_sel[row] = paste(f_df_genes_annotation$pos_in_gene[which(unique_genes_qtl[[level]][row] == f_df_genes_annotation[[level]])], collapse = ";")
    
    for (tax in c("class","order","family","genus","species")){
      tax_list = found_peaks[[tax]][which(found_peaks[[level]] == unique_genes_qtl[[level]][row])]
      unique_genes_qtl[[paste(tax, "_count", sep="")]][row] = length(unique(tax_list))
      unique_genes_qtl[[paste(tax, "_list", sep="")]][row] = paste(unique(tax_list), collapse = ";")
    }
  }
  unique_genes_qtl = unique_genes_qtl[order(unique_genes_qtl$QTL_peak_count, decreasing=TRUE),]
  write.csv(unique_genes_qtl, file = paste("./analysis_output/qtl_found_peaks_SNP_gene_analysis_",level, "_",file_info, ".csv", sep=""), row.names = F)
}


# --------------------------------------------------------------------------------------
##### Taxonomy analysis
# check snps in QTL / selected SNVs / full dataset per taxonomy

for (tax in c("class","order","family","genus","species")){
  
  unique_tax = data.frame(unique(fisher_selected_df[[tax]]))
  colnames(unique_tax) = tax
  for (row in 1:nrow(unique_tax)){
    unique_tax$QTL_peak_count[row] = length(which(unique_tax[[tax]][row] == found_peaks[[tax]]))
    unique_tax$SNP_count[row] = length(which(unique_tax[[tax]][row] == fisher_selected_df[[tax]]))
    unique_tax$unfiltered_SNP_count[row] = length(which(unique_tax[[tax]][row] == merged_df_taxonomy[[tax]]))
    unique_tax$total_length[row] = 0
    unique_tax$sel_SNP_per_kbp[row] = 0
    unique_tax$raw_SNP_per_kbp[row] = 0
    
    # QTL statistics
    qtl_snv_pos_list = unique(found_peaks$snv_pos[which(unique_tax[[tax]][row] == found_peaks[[tax]])])
    
    unique_tax$QTL_N_either_count[row]  = length(which(fisher_selected_df$mutation_type.x[which(fisher_selected_df$snv_pos %in% qtl_snv_pos_list)] == "N" | 
                                                         fisher_selected_df$mutation_type.y[which(fisher_selected_df$snv_pos %in% qtl_snv_pos_list)] == "N"))
    unique_tax$QTL_S_both_count[row]    = length(which(fisher_selected_df$mutation_type.x[which(fisher_selected_df$snv_pos %in% qtl_snv_pos_list)] == "S" & 
                                                         fisher_selected_df$mutation_type.y[which(fisher_selected_df$snv_pos %in% qtl_snv_pos_list)] == "S"))
    unique_tax$QTL_I_both_count[row]    = length(which(fisher_selected_df$mutation_type.x[which(fisher_selected_df$snv_pos %in% qtl_snv_pos_list)] == "I" & 
                                                         fisher_selected_df$mutation_type.y[which(fisher_selected_df$snv_pos %in% qtl_snv_pos_list)] == "I"))
    unique_tax$QTL_discrepancy[row]     = unique_tax$QTL_peak_count[row] - sum(unique_tax$QTL_N_either_count[row], unique_tax$QTL_S_both_count[row], unique_tax$QTL_I_both_count[row])
    
    unique_tax$QTL_N_either_perc[row] = unique_tax$QTL_N_either_count[row] / unique_tax$QTL_peak_count[row] * 100
    unique_tax$QTL_S_both_perc[row] = unique_tax$QTL_S_both_count[row] / unique_tax$QTL_peak_count[row] * 100
    unique_tax$QTL_N_S_ratio_total[row] = unique_tax$QTL_N_either_perc[row] / unique_tax$QTL_S_both_perc[row]
    unique_tax$QTL_I_both_perc[row] = unique_tax$QTL_I_both_count[row] / unique_tax$QTL_peak_count[row] * 100
    
    unique_tax$QTL_wild[row] = length(which(found_peaks$allele[which(unique_tax[[tax]][row] == found_peaks[[tax]])] == "wild"))
    unique_tax$QTL_modern[row] = length(which(found_peaks$allele[which(unique_tax[[tax]][row] == found_peaks[[tax]])] == "modern"))
    
    # Selected SNVs - NOT QTL SPECIFIC, just for easy statistics
    unique_tax$N_either_count[row] = length(which(fisher_selected_df$mutation_type.x[which(unique_tax[[tax]][row] == fisher_selected_df[[tax]])] == "N" | 
                                                    fisher_selected_df$mutation_type.y[which(unique_tax[[tax]][row] == fisher_selected_df[[tax]])] == "N"))
    unique_tax$S_both_count[row] =   length(which(fisher_selected_df$mutation_type.x[which(unique_tax[[tax]][row] == fisher_selected_df[[tax]])] == "S" & 
                                                    fisher_selected_df$mutation_type.y[which(unique_tax[[tax]][row] == fisher_selected_df[[tax]])] == "S"))
    unique_tax$I_both_count[row] =   length(which(fisher_selected_df$mutation_type.x[which(unique_tax[[tax]][row] == fisher_selected_df[[tax]])] == "I" & 
                                                    fisher_selected_df$mutation_type.y[which(unique_tax[[tax]][row] == fisher_selected_df[[tax]])] == "I"))
    unique_tax$discrepancy[row] = unique_tax$SNP_count[row] - sum(unique_tax$N_either_count[row], unique_tax$S_both_count[row], unique_tax$I_both_count[row])
    
    unique_tax$N_either_perc[row] = unique_tax$N_either_count[row] / unique_tax$SNP_count[row] * 100
    unique_tax$S_both_perc[row] = unique_tax$S_both_count[row] / unique_tax$SNP_count[row] * 100
    unique_tax$N_S_ratio_total[row] = unique_tax$N_either_perc[row] / unique_tax$S_both_perc[row]
    unique_tax$I_both_perc[row] = unique_tax$I_both_count[row] / unique_tax$SNP_count[row] * 100
    
    unique_tax$P_N_perc[row] = length(which(fisher_selected_df$mutation_type.x[which(unique_tax[[tax]][row] == fisher_selected_df[[tax]])] == "N")) / unique_tax$SNP_count[row] * 100
    unique_tax$P_S_perc[row] = length(which(fisher_selected_df$mutation_type.x[which(unique_tax[[tax]][row] == fisher_selected_df[[tax]])] == "S")) / unique_tax$SNP_count[row] * 100
    unique_tax$P_N_S_ratio[row] = unique_tax$P_N_perc[row] / unique_tax$P_S_perc[row]
    unique_tax$MM_N_perc[row] = length(which(fisher_selected_df$mutation_type.y[which(unique_tax[[tax]][row] == fisher_selected_df[[tax]])] == "N")) / unique_tax$SNP_count[row] * 100
    unique_tax$MM_S_perc[row] = length(which(fisher_selected_df$mutation_type.y[which(unique_tax[[tax]][row] == fisher_selected_df[[tax]])] == "S")) / unique_tax$SNP_count[row] * 100
    unique_tax$MM_N_S_ratio[row] = unique_tax$MM_N_perc[row] / unique_tax$MM_S_perc[row]
    
    # contig analysis
    tax_contigs = fisher_selected_df$scaffold.x[which(unique_tax[[tax]][row] == fisher_selected_df[[tax]])]
    for (row2 in 1:length(tax_contigs)){
      unique_tax$total_length[row] = unique_tax$total_length[row] + as.numeric(strsplit(strsplit(tax_contigs[row2], "_length_")[[1]][2], "_cov_")[[1]][1])
    }
    unique_tax$sel_SNP_per_kbp[row] = unique_tax$SNP_count[row] / unique_tax$total_length[row] * 1000
    unique_tax$raw_SNP_per_kbp[row] = unique_tax$unfiltered_SNP_count[row] / unique_tax$total_length[row] * 1000
    
    tax_contig_qtl = found_peaks$short_contig_id[which(unique_tax[[tax]][row] == found_peaks[[tax]])]
    unique_tax$contig_count[row] = length(unique(tax_contig_qtl))
    # unique_tax$contig_list[row] = paste(tax_contigs, collapse = ";") # takes up too much space in the csv file
    
    # QTL-unique parameters
    plant_chr_list = found_peaks$chr[which(found_peaks[[tax]] == unique_tax[[tax]][row])]
    plant_pos_list = found_peaks$pos[which(found_peaks[[tax]] == unique_tax[[tax]][row])]
    plant_chr_pos_list = paste(plant_chr_list, "_", plant_pos_list, collapse = ";", sep = "")
    unique_tax$plant_QTL_loci[row] = plant_chr_pos_list
    
    # total gene count is the number of attribute_pos per "level"
    gene_list = unique(found_peaks$prokka_gene_id[which(found_peaks[[tax]] == unique_tax[[tax]][row])])
    unique_genes_qtl$total_gene_count[row] = length(gene_list)
    unique_genes_qtl$prokka_gene_id_list[row] = paste(gene_list, collapse = ";")
  }
  
  for (row in 1:nrow(unique_tax)){
    unique_tax$QTL_percentage[row] = unique_tax$QTL_peak_count[row] / sum(unique_tax$QTL_peak_count) * 100
    unique_tax$SNP_percentage[row] = unique_tax$SNP_count[row] / sum(unique_tax$SNP_count) * 100
    unique_tax$unfiltered_SNP_percentage[row] = unique_tax$unfiltered_SNP_count[row] / dim(merged_df)[1] * 100
  }
  unique_tax = unique_tax[order(unique_tax$QTL_peak_count, decreasing=TRUE),]
  write.csv(unique_tax, file = paste("./analysis_output/filtered_SNP_summary_taxonomy_",tax, "_", file_info, ".csv", sep=""))
} 




#####
##### Analysis of taxonomy
# data is created above
# TODO
# unique_genus = read.table("filtered_SNP_summary_taxonomy_genus.csv", header=T, sep=",")
# 
# # ggplot(data.frame(unique_genus$P_dN_dS_ratio))
# 
# hist(unique_genus[which(unique_genus$SNP_count>10),]$SNP_per_kbp, xlab = "SNPs per kbp", main = "Genus")
# 
# boxplot(data.frame(unique_genus[which(unique_genus$SNP_count>10),]$P_dN_dS_ratio, unique_genus[which(unique_genus$SNP_count>10),]$MM_dN_dS_ratio), names = c("P", "MM"), ylab = "N / S")
# 





# GO term enrichment analysis --------------------------------------------------------------------------------------------------------------------------

# Preparation of GO annotations

# check eggNOG file
length(which(eggnog_output_file$GO_terms != ""))
nrow(eggnog_output_file)

# make list of "gene universe" and write to map file
rhizosphere_gene_list = data.frame(prokka_gene_id_list)
rhizosphere_gene_list$GO_terms = rep("", nrow(rhizosphere_gene_list))

# eggNOG file is not entirely complete for some reason, so also add the genes from the rhizosphere that are not yet in the eggNOG file

# create map file: ---------------
# add annotations for enrichment analysis and add genes not present in eggnog file
for (i in 1:nrow(prokka_gene_id_list)){
  if (!any(which(eggnog_output_file$Query == prokka_gene_id_list$contig_gene[i]))){
    rhizosphere_gene_list$GO_terms[i] = ""
  } else {
    rhizosphere_gene_list$GO_terms[i] = eggnog_output_file$GO_terms[which(eggnog_output_file$Query == prokka_gene_id_list$contig_gene[i])]
  }
}

# add GO terms to genes with the same annotations
# Find overlap in GO terms for terms with the same annotation that are already present
for (i in 1:length(unique(prokka_gene_id_list$annotation))){
  annotation_name = unique(prokka_gene_id_list$annotation)[i]
  go_set = ""
  go_list = rhizosphere_gene_list$GO_terms[which(rhizosphere_gene_list$annotation == annotation_name)] # not parsed yet
  go_list = go_list[go_list != ""] # filter elements without GO terms
  
  if (length(go_list) == 0){next} # if there are no GO terms, there is nothing to be done
  for (j in 1:length(go_list)){
    go_terms = unlist(strsplit(go_list[j], ",")) # parse go terms
    if (j == 1){go_set = go_terms} # initialise list of GO terms
    go_set = intersect(go_terms, go_set) # intersect of GO terms is the full list
    
  }
  if (length(go_set) == 0){go_set = ""} # in the case none of the GO terms overlap and none are left
  rhizosphere_gene_list$GO_terms[which(rhizosphere_gene_list$annotation == annotation_name & rhizosphere_gene_list$GO_terms == "")] = paste(go_set,collapse = ",")
}

# add GO mapping format: prokka_id
for (i in 1:nrow(prokka_gene_id_list)){
  go_list = paste(unlist(strsplit(rhizosphere_gene_list$GO_terms[i], ",")), collapse = ", ")
  rhizosphere_gene_list$topgo_map[i] = paste(rhizosphere_gene_list$prokka_gene_id[i], "\t", go_list, sep = "", collapse = "")
  rhizosphere_gene_list$contig_id[i] = paste(unlist(strsplit(rhizosphere_gene_list$contig_gene[i], "_"))[1:6], collapse = "_")
}
rhizosphere_gene_list = merge(rhizosphere_gene_list, taxonomy, by = "contig_id")


write.table(file = "./GO_analysis/rhizosphere_enriched_go_universe.map", subset(rhizosphere_gene_list, select = topgo_map), col.names = F, row.names = F, quote = F)




# GO enrichment ---------------------------

geneID2GO = readMappings(file = "./GO_analysis/rhizosphere_enriched_go_universe.map")
GO2geneID = inverseList(geneID2GO)
geneNames = names(geneID2GO)

length(geneNames)
length(GO2geneID)

# perform test on all rhizosphere QTL peaks and full rhizosphere gene universe
interesting_gene_list = unique(found_peaks$prokka_gene_id)

geneList = factor(as.integer(geneNames %in% interesting_gene_list))
names(geneList) = geneNames

# testing
interesting_gene_list[which(!(interesting_gene_list %in% geneNames[geneNames %in% interesting_gene_list]))]
length(which(geneList != 0))
length(unique(found_peaks$prokka_gene_id)) - 1

# create GOdata object
GOdata = new("topGOdata", ontology = c("BP"), allGenes = geneList, annot = annFUN.gene2GO, gene2GO = geneID2GO)

graph(GOdata) # number of nodes is number of GO terms tested
length(genes(GOdata))
numGenes(GOdata)
length(sigGenes(GOdata))
# GOdataGraph = graph(GOdata)

resultFisher  = runTest(GOdata, algorithm = "classic", statistic = "fisher")
# resultKS      = runTest(GOdata, algorithm = "classic", statistic = "ks")
# resultKS.elim = runTest(GOdata, algorithm = "elim", statistic = "ks")

pvalFisher    = score(resultFisher)
# pvalKS        = score(resultKS)
# pvalKS.elim   = score(resultKS.elim)

# hist(pvalFisher, 20, xlab = "p-value")
# hist(pvalKS, 50, xlab = "p-value")
# hist(pvalKS.elim, 50, xlab = "p-value")

# sign_pvalFisher = pvalFisher[pvalFisher < 0.05]
# View(sign_pvalFisher[order(sign_pvalFisher)])

# remove GO terms that are not in interesting_gene_list / or not with their children: topGO does not remove GO terms not in my dataset, so BH will be too stringent
pvalFisher = pvalFisher[order(pvalFisher)]
pval_list = data.frame(names(pvalFisher))
colnames(pval_list) = "go_term"
for (i in 1:nrow(pval_list)){ # intitialise dataframe
  pval_list$pvalFisher[i] = pvalFisher[i]
}

pval_list_new = data.frame(pval_list[0:0,]) # initialise smaller dataframe with only GO terms from interesting genes
for (i in 1:nrow(pval_list)){
  any_go_has_genes = F
  GO_term = pval_list$go_term[i]
  GO_term_genes = GO2geneID[GO_term]
  if (any(interesting_gene_list %in% unlist(unname(GO_term_genes)))){
    any_go_has_genes = T
  } else {
    GO_terms_offspring = unlist(unname(as.list(GOBPOFFSPRING[GO_term])))
    for (j in 1:length(GO_terms_offspring)){
      GO_term_child = GO_terms_offspring[j]
      GO_term_child_genes = GO2geneID[GO_term_child]
      if (any(interesting_gene_list %in% unlist(unname(GO_term_child_genes)))){
        any_go_has_genes = T
        break
      }
    }
  }
  if (any_go_has_genes){
    pval_list_new = rbind(pval_list_new, pval_list[i,])
  }
}

# Benjamini-Hochberg multiple testing correction
for (i in 1:nrow(pval_list_new)){
  pval_list_new$theshold[i] = i / nrow(pval_list_new) * 0.05
}

if (length(which(pval_list_new$pvalFisher <= pval_list_new$theshold)) > 0){
  topNodes_in_graph = max(which(pval_list_new$pvalFisher <= pval_list_new$theshold))
} else {
  topNodes_in_graph = 0
}


# create dataframe
# allRes <- GenTable(GOdata, classicFisher = resultFisher,
#                    classicKS = resultKS, elimKS = resultKS.elim,
#                    orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = topNodes_in_graph) 
allRes <- GenTable(GOdata, classicFisher = resultFisher,
                   orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = topNodes_in_graph)
allRes = cbind(allRes, pval_list_new$theshold[1:topNodes_in_graph])
colnames(allRes) = c(colnames(allRes)[1:(ncol(allRes)-1)], "BH_threshold")

# add genes to GO terms
for (row in 1:nrow(allRes)){
  GO_term = allRes$GO.ID[row]
  GO_term_genes = GO2geneID[GO_term]
  go_interesting_gene_list = c()
  
  # add direct genes with GO term
  if (any(interesting_gene_list %in% unlist(unname(GO_term_genes)))){
    go_interesting_gene_list = c(go_interesting_gene_list, interesting_gene_list[which(interesting_gene_list %in% unlist(unname(GO_term_genes)))])
  }
  
  # add genes that have offspring of GO term
  GO_terms_offspring = unlist(unname(as.list(GOBPOFFSPRING[GO_term])))
  for (j in 1:length(GO_terms_offspring)){
    GO_term_child = GO_terms_offspring[j]
    GO_term_child_genes = GO2geneID[GO_term_child]
    if (any(interesting_gene_list %in% unlist(unname(GO_term_child_genes)))){
      go_interesting_gene_list = c(go_interesting_gene_list, interesting_gene_list[which(interesting_gene_list %in% unlist(unname(GO_term_child_genes)))])
    }
  }
  if (length(go_interesting_gene_list > 0)){
    allRes$gene_list[row] = paste(go_interesting_gene_list, collapse = ";")
    # add annotations
    GO_genes = unique(interesting_gene_list[which(interesting_gene_list %in% go_interesting_gene_list)])
    allRes$annotations[row] = paste(unique(found_peaks$annotation[which(found_peaks$prokka_gene_id %in% GO_genes)]), collapse = ";")
    allRes$genus[row] = paste(unique(found_peaks$genus[which(found_peaks$prokka_gene_id %in% GO_genes)]), collapse = ";")
  } else {
    allRes$gene_list[row] = ""
    allRes$annotations[row] = ""
    allRes$genus[row] = ""
  }
}


# showSigOfNodes(GOdata, score(resultFisher), firstSigNodes = length(pvalFisher[pvalFisher < 0.05]), useInfo = 'all')

printGraph(GOdata, resultFisher, firstSigNodes = topNodes_in_graph, fn.prefix = paste("./GO_analysis/GO_BP_rhizosphere_QTL_peaks_", file_info, sep=""), useInfo = "all", pdfSW = TRUE)

write.table(allRes, file = paste("./GO_analysis/GO_BP_rhizosphere_QTL_peaks_", file_info, ".tsv", sep=""), quote = F, sep = "\t", row.names = F, col.names = T)


# analysis of wild and modern QTLs
# divide interesting gene list in two parts: QTL associated to wild or modern plants
for (allele in c("modern", "wild")){
  
  interesting_gene_list = unique(found_peaks$prokka_gene_id[which(found_peaks$allele == allele)])
  geneList = factor(as.integer(geneNames %in% interesting_gene_list))
  names(geneList) = geneNames
  
  interesting_gene_list[which(!(interesting_gene_list %in% geneNames[geneNames %in% interesting_gene_list]))]
  length(which(geneList != 0))
  length(unique(found_peaks$prokka_gene_id)) - 1
  
  GOdata = new("topGOdata", ontology = c("BP"), allGenes = geneList, annot = annFUN.gene2GO, gene2GO = geneID2GO)
  
  numGenes(GOdata)
  length(sigGenes(GOdata))
  # GOdataGraph = graph(GOdata)
  
  resultFisher  = runTest(GOdata, algorithm = "classic", statistic = "fisher")
  # resultKS      = runTest(GOdata, algorithm = "classic", statistic = "ks")
  # resultKS.elim = runTest(GOdata, algorithm = "elim", statistic = "ks")
  
  pvalFisher    = score(resultFisher)
  # pvalKS        = score(resultKS)
  # pvalKS.elim   = score(resultKS.elim)
  
  # remove GO terms that are not in interesting_gene_list / or not with their children: topGO does not remove GO terms not in my dataset, so BH will be too stringent
  pvalFisher = pvalFisher[order(pvalFisher)]
  pval_list = data.frame(names(pvalFisher))
  colnames(pval_list) = "go_term"
  for (i in 1:nrow(pval_list)){ # intitialise dataframe
    pval_list$pvalFisher[i] = pvalFisher[i]
  }
  
  pval_list_new = data.frame(pval_list[0:0,]) # initialise smaller dataframe with only GO terms from interesting genes
  for (i in 1:nrow(pval_list)){
    any_go_has_genes = F
    GO_term = pval_list$go_term[i]
    GO_term_genes = GO2geneID[GO_term]
    if (any(interesting_gene_list %in% unlist(unname(GO_term_genes)))){
      any_go_has_genes = T
    } else {
      GO_terms_offspring = unlist(unname(as.list(GOBPOFFSPRING[GO_term])))
      for (j in 1:length(GO_terms_offspring)){
        GO_term_child = GO_terms_offspring[j]
        GO_term_child_genes = GO2geneID[GO_term_child]
        if (any(interesting_gene_list %in% unlist(unname(GO_term_child_genes)))){
          any_go_has_genes = T
          break
        }
      }
    }
    if (any_go_has_genes){
      pval_list_new = rbind(pval_list_new, pval_list[i,])
    }
  }
  
  # Benjamini-Hochberg multiple testing correction
  for (i in 1:nrow(pval_list_new)){
    pval_list_new$theshold[i] = i / nrow(pval_list_new) * 0.05
  }
  
  if (length(which(pval_list_new$pvalFisher <= pval_list_new$theshold)) > 0){
    topNodes_in_graph = max(which(pval_list_new$pvalFisher <= pval_list_new$theshold))
  } else {
    topNodes_in_graph = 0
  }
  if (topNodes_in_graph == 0){next}
  
  # allRes <- GenTable(GOdata, classicFisher = resultFisher,
  #                    classicKS = resultKS, elimKS = resultKS.elim,
  #                    orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = topNodes_in_graph)
  allRes <- GenTable(GOdata, classicFisher = resultFisher,
                     orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = topNodes_in_graph)
  allRes = cbind(allRes, pval_list_new$theshold[1:topNodes_in_graph])
  colnames(allRes) = c(colnames(allRes)[1:(ncol(allRes)-1)], "BH_threshold")
  
  # add genes to GO terms
  for (row in 1:nrow(allRes)){
    GO_term = allRes$GO.ID[row]
    GO_term_genes = GO2geneID[GO_term]
    go_interesting_gene_list = c()
    
    # add direct genes with GO term
    if (any(interesting_gene_list %in% unlist(unname(GO_term_genes)))){
      go_interesting_gene_list = c(go_interesting_gene_list, interesting_gene_list[which(interesting_gene_list %in% unlist(unname(GO_term_genes)))])
    }
    
    # add genes that have offspring of GO term
    GO_terms_offspring = unlist(unname(as.list(GOBPOFFSPRING[GO_term])))
    for (j in 1:length(GO_terms_offspring)){
      GO_term_child = GO_terms_offspring[j]
      GO_term_child_genes = GO2geneID[GO_term_child]
      if (any(interesting_gene_list %in% unlist(unname(GO_term_child_genes)))){
        go_interesting_gene_list = c(go_interesting_gene_list, interesting_gene_list[which(interesting_gene_list %in% unlist(unname(GO_term_child_genes)))])
      }
    }
    if (length(go_interesting_gene_list > 0)){
      allRes$gene_list[row] = paste(go_interesting_gene_list, collapse = ";")
      # add annotations
      GO_genes = unique(interesting_gene_list[which(interesting_gene_list %in% go_interesting_gene_list)])
      allRes$annotations[row] = paste(unique(found_peaks$annotation[which(found_peaks$prokka_gene_id %in% GO_genes)]), collapse = ";")
      allRes$genus[row] = paste(unique(found_peaks$genus[which(found_peaks$prokka_gene_id %in% GO_genes)]), collapse = ";")
    } else {
      allRes$gene_list[row] = ""
      allRes$annotations[row] = ""
      allRes$genus[row] = ""
    }
  }
  
  # showSigOfNodes(GOdata, score(resultFisher), firstSigNodes = length(pvalFisher[pvalFisher < 0.05]), useInfo = 'all')
  
  printGraph(GOdata, resultFisher, firstSigNodes = topNodes_in_graph, fn.prefix = paste("./GO_analysis/GO_BP_rhizosphere_QTL_peaks_", allele, "_", file_info, sep = ""), useInfo = "all", pdfSW = TRUE)
  
  write.table(allRes, file = paste("./GO_analysis/GO_BP_rhizosphere_QTL_peaks_", allele, "_", file_info, ".tsv", sep = "")
              , quote = F, sep = "\t", row.names = F, col.names = T)
  
}



# GO enrichment for taxonomic level ----------------------------------------------------

level = "family"

if (!file.exists(paste("./GO_analysis/", level, sep = ""))){dir.create(file.path(paste("./GO_analysis/", level, sep = "")))}

qtl_unique_level = unique(found_peaks[[level]][order(found_peaks[[level]])])

# loop over all genera for GO enrichment analysis
for (i in 1:length(qtl_unique_level)){
  tax = qtl_unique_level[i]
  # tax = "Streptomycetaceae"
  
  # create new universe with all GENES of taxonomy, so do not take from full universe!
  tax_rhizosphere_gene_list = rhizosphere_gene_list$prokka_gene_id[which(rhizosphere_gene_list[[level]] == tax)]
  geneID2GO_tax_in_universe = geneID2GO[which(names(geneID2GO) %in% tax_rhizosphere_gene_list)]
  
  geneNames_tax = names(geneID2GO_tax_in_universe)
  length(geneNames_tax)
  
  for (allele in c("modern", "wild")){
    interesting_gene_list = unique(found_peaks$prokka_gene_id[which((found_peaks[[level]] == tax) & found_peaks$allele == allele)])
    
    geneList = factor(as.integer(geneNames_tax %in% interesting_gene_list))
    names(geneList) = geneNames_tax
    
    interesting_gene_list[which(!(interesting_gene_list %in% geneNames_tax[geneNames_tax %in% interesting_gene_list]))]
    length(which(geneList != 0))
    length(tax_rhizosphere_gene_list)
    
    if (all(lengths(geneID2GO_tax_in_universe[which(geneList != 0)]) == 0)){next} # skip if there are no GO terms at all
    
    GOdata_tax = new("topGOdata", ontology = c("BP"), allGenes = geneList, annot = annFUN.gene2GO, gene2GO = geneID2GO_tax_in_universe)
    
    numGenes(GOdata_tax)
    length(sigGenes(GOdata_tax))
    
    if (length(sigGenes(GOdata_tax)) == 0){next} # skip to next iteration if there are no enriched GO terms to prevent bugs
    
    resultFisher  = runTest(GOdata_tax, algorithm = "classic", statistic = "fisher")
    # resultKS      = runTest(GOdata_tax, algorithm = "classic", statistic = "ks")
    # resultKS.elim = runTest(GOdata_tax, algorithm = "elim", statistic = "ks")
    
    pvalFisher    = score(resultFisher)
    # pvalKS        = score(resultKS)
    # pvalKS.elim   = score(resultKS.elim)
    
    # remove GO terms that are not in interesting_gene_list / or not with their children: topGO does not remove GO terms not in my dataset, so BH will be too stringent
    pvalFisher = pvalFisher[order(pvalFisher)]
    pval_list = data.frame(names(pvalFisher))
    colnames(pval_list) = "go_term"
    for (i in 1:nrow(pval_list)){ # intitialise dataframe
      pval_list$pvalFisher[i] = pvalFisher[i]
    }
    
    pval_list_new = data.frame(pval_list[0:0,]) # initialise smaller dataframe with only GO terms from interesting genes
    for (i in 1:nrow(pval_list)){
      any_go_has_genes = F
      GO_term = pval_list$go_term[i]
      GO_term_genes = GO2geneID[GO_term]
      if (any(interesting_gene_list %in% unlist(unname(GO_term_genes)))){
        any_go_has_genes = T
      } else {
        GO_terms_offspring = unlist(unname(as.list(GOBPOFFSPRING[GO_term])))
        for (j in 1:length(GO_terms_offspring)){
          GO_term_child = GO_terms_offspring[j]
          GO_term_child_genes = GO2geneID[GO_term_child]
          if (any(interesting_gene_list %in% unlist(unname(GO_term_child_genes)))){
            any_go_has_genes = T
            break
          }
        }
      }
      if (any_go_has_genes){
        pval_list_new = rbind(pval_list_new, pval_list[i,])
      }
    }
    
    # Benjamini-Hochberg multiple testing correction
    for (i in 1:nrow(pval_list_new)){
      pval_list_new$theshold[i] = i / nrow(pval_list_new) * 0.05
    }
    
    if (length(which(pval_list_new$pvalFisher <= pval_list_new$theshold)) > 0){
      topNodes_in_graph = max(which(pval_list_new$pvalFisher <= pval_list_new$theshold))
    } else {
      topNodes_in_graph = 0
    }
    
    if (topNodes_in_graph == 0){next} # if there are no significant GO terms, skip to next iteration
    
    # allRes <- GenTable(GOdata_tax, classicFisher = resultFisher,
    #                    classicKS = resultKS, elimKS = resultKS.elim,
    #                    orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = topNodes_in_graph)
    allRes <- GenTable(GOdata_tax, classicFisher = resultFisher,
                       orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = topNodes_in_graph)
    allRes = cbind(allRes, pval_list_new$theshold[1:topNodes_in_graph])
    colnames(allRes) = c(colnames(allRes)[1:(ncol(allRes)-1)], "BH_threshold")
    
    # add genes to GO terms
    for (row in 1:nrow(allRes)){
      GO_term = allRes$GO.ID[row]
      GO_term_genes = GO2geneID[GO_term]
      go_interesting_gene_list = c()
      
      # add direct genes with GO term
      if (any(interesting_gene_list %in% unlist(unname(GO_term_genes)))){
        go_interesting_gene_list = c(go_interesting_gene_list, interesting_gene_list[which(interesting_gene_list %in% unlist(unname(GO_term_genes)))])
      }
      
      # add genes that have offspring of GO term
      GO_terms_offspring = unlist(unname(as.list(GOBPOFFSPRING[GO_term])))
      for (j in 1:length(GO_terms_offspring)){
        GO_term_child = GO_terms_offspring[j]
        GO_term_child_genes = GO2geneID[GO_term_child]
        if (any(interesting_gene_list %in% unlist(unname(GO_term_child_genes)))){
          go_interesting_gene_list = c(go_interesting_gene_list, interesting_gene_list[which(interesting_gene_list %in% unlist(unname(GO_term_child_genes)))])
        }
      }
      if (length(go_interesting_gene_list > 0)){
        allRes$gene_list[row] = paste(go_interesting_gene_list, collapse = ";")
        # add annotations
        GO_genes = unique(interesting_gene_list[which(interesting_gene_list %in% go_interesting_gene_list)])
        allRes$annotations[row] = paste(unique(found_peaks$annotation[which(found_peaks$prokka_gene_id %in% GO_genes)]), collapse = ";")
        allRes$genus[row] = paste(unique(found_peaks$genus[which(found_peaks$prokka_gene_id %in% GO_genes)]), collapse = ";")
      } else {
        allRes$gene_list[row] = ""
        allRes$annotations[row] = ""
        allRes$genus[row] = ""
      }
    }
    
    # showSigOfNodes(GOdata_tax, score(resultFisher), firstSigNodes = length(pvalFisher[pvalFisher < 0.05]), useInfo = 'all')
    
    printGraph(GOdata_tax, resultFisher, firstSigNodes = topNodes_in_graph, fn.prefix = paste("./GO_analysis/", level, "/", "GO_BP_rhizosphere_QTL_peaks_", level, "_", tax, "_", allele, "_", file_info, sep = ""), useInfo = "all", pdfSW = TRUE)
    
    write.table(allRes, file = paste("./GO_analysis/", level, "/", "GO_BP_rhizosphere_QTL_peaks_", level, "_", tax, "_", allele, "_", file_info, ".tsv", sep = "")
                , quote = F, sep = "\t", row.names = F, col.names = T)
    
  }
}



# final analyses ----------------

# find genes associated with GO term
interesting_gene_list = unique(found_peaks$prokka_gene_id[which(found_peaks$allele == "wild")])

GO_term = "GO:0015893"
GO_term_genes = GO2geneID[GO_term]
go_interesting_gene_list = c()

if (any(interesting_gene_list %in% unlist(unname(GO_term_genes)))){
  go_interesting_gene_list = c(go_interesting_gene_list, interesting_gene_list[which(interesting_gene_list %in% unlist(unname(GO_term_genes)))])
}

GO_terms_offspring = unlist(unname(as.list(GOBPOFFSPRING[GO_term])))
for (j in 1:length(GO_terms_offspring)){
  GO_term_child = GO_terms_offspring[j]
  GO_term_child_genes = GO2geneID[GO_term_child]
  if (any(interesting_gene_list %in% unlist(unname(GO_term_child_genes)))){
    go_interesting_gene_list = c(go_interesting_gene_list, interesting_gene_list[which(interesting_gene_list %in% unlist(unname(GO_term_child_genes)))])
  }
}

GO_0015893_peaks = unique(interesting_gene_list[which(interesting_gene_list %in% go_interesting_gene_list)])

length(GO_0015893_peaks)

unique(found_peaks$annotation[which(found_peaks$prokka_gene_id %in% GO_0015893_peaks)])



