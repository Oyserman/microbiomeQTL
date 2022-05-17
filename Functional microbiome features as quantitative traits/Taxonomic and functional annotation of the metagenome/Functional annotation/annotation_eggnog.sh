# activate the conda environment
conda activate eggnog

# annotate sub input faa file
emapper.py  \
-i /mnt/nfs/bioinfdata/home/NIOO/xinyap/data/contig_1k_prodigal/contigs_1k.0.faa \
--output /mnt/nfs/bioinfdata/ngs2/ME2/raaijmakers_group/RIL2019_analysis/final_assembly/eggnog/tmp/contigs_1k.0 \
-m diamond --cpu 20
