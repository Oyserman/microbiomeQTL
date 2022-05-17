conda activate checkm_env

checkm lineage_wf --tab_table \
--file /mnt/nfs/bioinfdata/ngs2/ME2/raaijmakers_group/RIL2019_analysis/final_assembly/binning/new/new_checkm_lineage_wf.tsv \
--threads 8 --extension "fa" \
/mnt/nfs/bioinfdata/ngs2/ME2/raaijmakers_group/RIL2019_analysis/final_assembly/binning/new/new_bins/ \
/mnt/nfs/bioinfdata/ngs2/ME2/raaijmakers_group/RIL2019_analysis/final_assembly/binning/new/new_checkm_output/
