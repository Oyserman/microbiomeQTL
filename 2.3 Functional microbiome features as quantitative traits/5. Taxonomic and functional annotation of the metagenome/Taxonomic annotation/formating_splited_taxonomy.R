
splited_taxonomy <- read.csv("~/data/final_assembly/kraken/2020_may/10kb_taxonomy_splited.csv", header = F, row.names = F)
splited_taxonomy_t <- t(splited_taxonomy)
df_splited_taxonomy <- data.frame(splited_taxonomy_t)[-1,]
colnames(df_splited_taxonomy) <- c("contig_id", "kingdom", "phylum", 
                                   "class", "order", "family", "genus", "species")

write.table(df_splited_taxonomy, "~/data/final_assembly/kraken/2020_may/contig_10kb_taxonomy.tsv", sep = "\t", row.names = F, quote = F)

contig_taxonomy <- read.delim("~/data/final_assembly/kraken/2020_may/contig_10kb_taxonomy.tsv")
