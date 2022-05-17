# create and activate kofamscan conda environment
# download ko list and hmm profiles
# set config.yml

# run kofamscan
## for A1
exec_annotation -f mapper \
-o /mnt/nfs/bioinfdata/ngs2/ME2/raaijmakers_group/RIL2019_analysis/final_assembly/kofam/test/a1_ko.txt \
/mnt/nfs/bioinfdata/ngs2/ME2/raaijmakers_group/RIL2019_analysis/final_assembly/prodigal/a1_contigs_1k/contigs_1k.faa

## for A3
exec_annotation -f mapper \
-o /mnt/nfs/bioinfdata/ngs2/ME2/raaijmakers_group/RIL2019_analysis/final_assembly/kofam/test/filtered_a3_ko.txt
/mnt/nfs/bioinfdata/ngs2/ME2/raaijmakers_group/RIL2019_analysis/final_assembly/prodigal/filtered_a3_contigs_1k.faa

# only keep annoated ORFs
cd /mnt/nfs/bioinfdata/ngs2/ME2/raaijmakers_group/RIL2019_analysis/final_assembly/kofam/test/ 
grep "K" a1_ko.txt > ko_annotated.txt
grep "K" filtered_a3_ko.txt >> ko_annotated.txt
