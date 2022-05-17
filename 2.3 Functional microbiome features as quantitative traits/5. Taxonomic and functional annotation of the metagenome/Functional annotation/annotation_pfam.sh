# activate the conda environment to run hmmer
conda activate hmmer

# download and prepare pfam database
wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam32.0/Pfam-A.hmm.gz
gunzip Pfam-A.hmm.gz
hmmpress Pfam-A.hmm

# annotate predicted genes against pfam database (search protein domains)
# hmmsearch --tblout out.txt -E 1e-5 --cpu 2 Pfam-A.hmm input_proteins.faa
hmmsearch --tblout"/mnt/nfs/bioinfdata/ngs2/ME2/raaijmakers_group/RIL2019_analysis/final_assembly/pfam/tmp/a3_contigs_1k.8.txt" -E 1e-5 --cpu 4 ~/pfam/Pfam-A.hmm \
"/mnt/nfs/bioinfdata/ngs2/ME2/raaijmakers_group/RIL2019_analysis/final_assembly/prodigal/filtered_a3_contigs_1k.8.faa"

