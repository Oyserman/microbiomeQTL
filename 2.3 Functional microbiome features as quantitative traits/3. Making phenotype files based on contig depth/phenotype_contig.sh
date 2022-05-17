# index bam files
for f in /mnt/nfs/bioinfdata/ngs2/ME2/raaijmakers_group/RIL2019_analysis/final_assembly/mapping/*.bam
do samtools index -b -@ 20 $f
done

# calculate genome coverage for all positions in BEDGRAPH format
for f in /mnt/nfs/bioinfdata/ngs2/ME2/raaijmakers_group/RIL2019_analysis/final_assembly/mapping/*.bam
do bedtools genomecov -bga -ibam $f \
> /mnt/nfs/bioinfdata/ngs2/ME2/raaijmakers_group/RIL2019_analysis/final_assembly/depth/$f.bed
done

# use python script to get contig coverage statistics
# format: contig_id	length	total_depth	avg_depth(=total_depth/length)	percentage_coverage
cd /mnt/nfs/bioinfdata/ngs2/ME2/raaijmakers_group/RIL2019_analysis/final_assembly/depth/
for f in *.bed
do id={f%.bed}
python ../for_ben/codes/contig_statistics.py $f > ${id}_contig.tsv
done
