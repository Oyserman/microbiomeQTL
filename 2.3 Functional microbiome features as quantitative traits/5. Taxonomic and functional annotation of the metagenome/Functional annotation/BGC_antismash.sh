cd /mnt/nfs/bioinfdata/ngs2/ME2/raaijmakers_group/RIL2019_analysis/final_assembly

# run antiSMASH to predict BGCs in the final assembly
conda activate antismash
antismash -c 12 \
--clusterhmmer --cb-knownclusters --cb-subclusters --cb-general \
--output-dir ./antismash_0519 \
--genefinding-tool prodigal-m \
merged_a13_rm100_10k.fa

# run BiG-SCAPE to create BGC seq-similarity network
conda activate bigscape
python bigscape.py -i ./antismash_0519/ \
-o ./big_scape_0522/ \
--pfam_dir ~/db/pfam/ \
-c 8 --exclude_gbk_str merged \
--cutoffs 0.75 0.80 0.85 -v \
> ./big_scape_0522/run.log

# visualise the network in Cytoscape
# input file: big_scape_mix/network_files/2020-06-09_17-26-58_hybrids_glocal/mix/mix_c0.80.network
