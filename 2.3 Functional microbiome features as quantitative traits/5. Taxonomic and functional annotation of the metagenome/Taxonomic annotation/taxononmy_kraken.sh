# activate conda environment for taxonomy classification
conda activate kraken2

# build kraken2 database
$DBNAME=/mnt/nfs/bioinfdata/ngs2/ME2/raaijmakers_group/RIL2019_analysis/final_assembly/kraken/kraken2-microbial-fatfree_202003
kraken2-build --download-library archaea --db $DBNAME
kraken2-build --download-library viral --db $DBNAME
kraken2-build --download-library bacteria --db $DBNAME
kraken2-build --download-library fungi --db $DBNAME
kraken2-build --download-library protozoa --db $DBNAME
kraken2-build --download-library UniVec_Core --db $DBNAME

# run kraken2 without --use-names
kraken2 --threads 20 \
--db /mnt/nfs/bioinfdata/ngs2/ME2/raaijmakers_group/RIL2019_analysis/final_assembly/kraken/kraken2-microbial-fatfree_202003 \
--report /mnt/nfs/bioinfdata/ngs2/ME2/raaijmakers_group/RIL2019_analysis/final_assembly/kraken/2020_may/10kb_without_use_names_report.tsv \
"/mnt/nfs/bioinfdata/ngs2/ME2/raaijmakers_group/RIL2019_analysis/final_assembly/merged_a13_rm100_10k.fa" \
> /mnt/nfs/bioinfdata/ngs2/ME2/raaijmakers_group/RIL2019_analysis/final_assembly/kraken/2020_may/10kb_without_use_names.kraken

##-----translate kraken2 output-----##
# activate the environment
conda deactivate
conda activate taxonkit

# translate taxonomy id in kraken2 standard output
python /mnt/nfs/bioinfdata/ngs2/ME2/raaijmakers_group/RIL2019_analysis/final_assembly/kraken/2020_may/translateKraken2.py \
--krakenout /mnt/nfs/bioinfdata/ngs2/ME2/raaijmakers_group/RIL2019_analysis/final_assembly/kraken/2020_may/10kb_without_use_names.kraken \
--translatedout /mnt/nfs/bioinfdata/ngs2/ME2/raaijmakers_group/RIL2019_analysis/final_assembly/kraken/2020_may/10kb_translated.kraken \
--taxdatadir /mnt/nfs/bioinfdata/ngs2/ME2/raaijmakers_group/RIL2019_analysis/final_assembly/kraken/2020_may/taxdatadir/

# split to 1 taxonomy 1 column
python /mnt/nfs/bioinfdata/ngs2/ME2/raaijmakers_group/RIL2019_analysis/final_assembly/for_ben/code/split_taxonomy_columns.py \
/mnt/nfs/bioinfdata/ngs2/ME2/raaijmakers_group/RIL2019_analysis/final_assembly/kraken/2020_may/10kb_translated.kraken

# convert csv file to tsv in R
/mnt/nfs/bioinfdata/ngs2/ME2/raaijmakers_group/RIL2019_analysis/final_assembly/for_ben/code/formating_splited_taxonomy.R
