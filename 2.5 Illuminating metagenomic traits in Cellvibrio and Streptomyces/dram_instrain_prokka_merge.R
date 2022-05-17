inStrain
dram.instrain
rilORF=read.delim(file = '../../../../Downloads/filtered_rhi_contig_id_210111_no_fasta.gff',
                  header = F)
rilORF$ORFid=str_split(rilORF$V9, ';', simplify = T)[,1] %>% 
  str_remove(., 'ID=')

rilORF$gene= str_match(rilORF$V9, "gene=(.*?);")[,2]
rilORF$COG= str_match(rilORF$V9, "COG:(.*?);")[,2]
rilORF$product= str_match(rilORF$V9, "product=(.*)")[,2]
rilORF$UniProt= str_match(rilORF$V9, "UniProtKB:(.*?);")[,2]


head(rilORF)
rilORF = rilORF %>% mutate(contigID=V1,
                           start=V4, end = V5, orientation=V7) %>% 
  select(ORFid, contigID, start, end, orientation , gene, product, COG, UniProt)

write_tsv(rilORF, '../../../../Desktop/RIL_ORF_all.tsv')  

dram.instrain
head(dramUpdate)
head(rilORF)
rilORF$node_pos=paste(rilORF$contigID, rilORF$start, rilORF$end, sep='_')
dramUpdate=dram %>%
  mutate(start=start_position,
         end=end_position,
         kegg=kegg_hit,
         orientation=if_else(strandedness == 1, '+', '-'),
         pfam=pfam_hits,
         cazy=cazy_hits, contigID = node) %>%
  select(c(node_pos, contigID, start, end, orientation, kegg, pfam, cazy))

test_merge=full_join(rilORF, dramUpdate)

library(RDPutils)
n <- c(1:length(test_merge$node_pos))
orf.names <- make_otu_names(n, otu_format="R")
orf.names = str_replace(orf.names, 'OTU', 'B2R')
head(orf.names)

test_merge$newORF=orf.names
tail(test_merge)
head(test_merge)

taxRIL=read.delim('../contig_10kb_taxonomy.tsv')
taxRIL$contigID=str_split(taxRIL$contig_id, '_length',simplify = T)[,1]
head(test_merge_tax)

test_merge_tax=left_join(test_merge, taxRIL[-1])

test_merge_tax = test_merge_tax %>%
  mutate(MAG=if_else(contigID  %in% dramUpdate$contigID, 'MAG72', ''))

write.table(test_merge_tax,'../microQTL_contig_ORF_table_with_taxonomy.txt', na="", quote = F, sep="\t")
