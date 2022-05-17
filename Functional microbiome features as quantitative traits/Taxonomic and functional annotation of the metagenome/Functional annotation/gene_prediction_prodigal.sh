# activate conda environment for gene prediction by prodigal
conda activate xinya

# predict genes from the first assembly (co-assembly of the parental samples)
pigz -c -d /mnt/nfs/bioinfdata/ngs2/ME2/raaijmakers_group/RIL2019_analysis/coassembly_spades/contigs_1k.fna.gz | \
prodigal -q -f gff -p meta \
-d ~/data/contig_1k_prodigal/contigs_1k.fna \
-a ~/data/contig_1k_prodigal/contigs_1k.faa \
-o ~/data/contig_1k_prodigal/contigs_1k.gff

# create symbotic link of the output
ln -s ~/data/contig_1k_prodigal/ /mnt/nfs/bioinfdata/ngs2/ME2/raaijmakers_group/RIL2019_analysis/prodigal/a1_contigs_1k/

# predict genes from the third assembly (MAPQ<20, unmapped/ambiguously mapped reads etc.)
prodigal -q -f gff -p meta \
-d /mnt/nfs/bioinfdata/ngs2/ME2/raaijmakers_group/RIL2019_analysis/progigal/filtered_a3_contigs_1k.fna
-a /mnt/nfs/bioinfdata/ngs2/ME2/raaijmakers_group/RIL2019_analysis/progigal/filtered_a3_contigs_1k.faa
-o /mnt/nfs/bioinfdata/ngs2/ME2/raaijmakers_group/RIL2019_analysis/progigal/filtered_a3_contigs_1k.gff
-i /mnt/nfs/bioinfdata/ngs2/ME2/raaijmakers_group/RIL2019_analysis/assembly_evaluation/test/a3_rm_cov100.fa
