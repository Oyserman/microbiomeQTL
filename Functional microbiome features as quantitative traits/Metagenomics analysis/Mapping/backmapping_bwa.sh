# index the reference assembly
bwa index /mnt/nfs/bioinfdata/ngs2/ME2/raaijmakers_group/RIL2019_analysis/final_assembly/merged_a13_rm100_10k.fa
# backmap raw read pairs to the assembly
for f in /mnt/nfs/bioinfdata/ngs2/ME2/raaijmakers_group/RIL2019_raw/RIL2019_raw_Shotgun/renamed/*R1.fq.gz; \
do f1=${f%_R1.fq.gz}; id=${f1##*/}; \
bwa mem -t 20 \
/mnt/nfs/bioinfdata/ngs2/ME2/raaijmakers_group/RIL2019_analysis/assembly_evaluation/test/merged_a3_rm_cov100.fa \
$f ${f1}_R2.fq.gz | \
samtools view -b | \
samtools sort -t 20 -O bam \
-o ~/data/unmapped_ril_assembly/mapping/${id}.bam ; done

# we can use A1, A1+A2 or A1+A3 as reference respectively for benchmarking
bwa index <ref>
for f in /mnt/nfs/bioinfdata/ngs2/ME2/raaijmakers_group/RIL2019_raw/RIL2019_raw_Shotgun/renamed/*R1.fq.gz; \
do f1=${f%_R1.fq.gz}; id=${f1##*/}; \
bwa mem -t 20 <indexed_ref> $f ${f1}_R2.fq.gz | \
samtools view -b | \
samtools sort -t 20 -O bam \
-o ${id}_sorted.bam ; done
