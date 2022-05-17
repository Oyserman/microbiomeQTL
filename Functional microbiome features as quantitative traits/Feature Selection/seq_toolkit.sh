# remove sequences with length < 10 kbp
# conda activate xinya #
seqtk seq -L 10000 original.fasta

# split a large fasta file to multiple (e.g. 6) small sub-files
# conda activate pyfasta #
pyfasta split -n 6 original.fasta
