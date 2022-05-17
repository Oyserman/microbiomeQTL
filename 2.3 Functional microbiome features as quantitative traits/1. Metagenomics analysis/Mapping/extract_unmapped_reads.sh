for i in $(seq 205 1 308); \
do samtools bam2fq \
-1 /mnt/nfs/bioinfdata/home/NIOO/xinyap/unmapped_paired_RIL_reads/r1/ril_${i}_r1.fq.gz \
-2 /mnt/nfs/bioinfdata/home/NIOO/xinyap/unmapped_paired_RIL_reads/r2/ril_${i}_r2.fq.gz \
-f 4 -@ 20 -n \
-s /mnt/nfs/bioinfdata/home/NIOO/xinyap/unmapped_paired_RIL_reads/singleton/ril_${i}_singleton.fq.gz \
/mnt/nfs/bioinfdata/ngs2/ME2/raaijmakers_group/RIL2019_analysis/coassembly_mapping/allMP_${i}_*.bam ; done

