# activate conda environment for binning
conda activate metabat2

# generate depth matrix that is needed for binning by metabat2
jgi_summarize_bam_contig_depths --outputDepth \
/mnt/nfs/bioinfdata/ngs2/ME2/raaijmakers_group/RIL2019_analysis/final_assembly/binning/depth_matrix.txt \
/mnt/nfs/bioinfdata/ngs2/ME2/raaijmakers_group/RIL2019_analysis/final_assembly/mapping/*.bam

# binning
metabat2 -i /mnt/nfs/bioinfdata/ngs2/ME2/raaijmakers_group/thomg/micro_qtl/spades_coassembly/samtools/contigs_1k.fna.gz \
-a /mnt/nfs/bioinfdata/ngs2/ME2/raaijmakers_group/RIL2019_analysis/final_assembly/binning/depth_matrix.txt \
-o /mnt/nfs/bioinfdata/ngs2/ME2/raaijmakers_group/RIL2019_analysis/final_assembly/binning/new/new_bins/bin -t 20
