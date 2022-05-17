#A2
r1=$(ls /mnt/nfs/bioinfdata/home/NIOO/xinyap/analyses/unmapped_paired_RIL_reads/r1/* -m)
r2=$(ls -m /mnt/nfs/bioinfdata/home/NIOO/xinyap/analyses/unmapped_paired_RIL_reads/r2/*)
se=$(ls -m /mnt/nfs/bioinfdata/home/NIOO/xinyap/analyses/unmapped_paired_RIL_reads/singleton/*)

megahit -1 ${r1} -2 ${r2} -r ${se} \
-o /mnt/nfs/bioinfdata/home/NIOO/xinyap/data/unmapped_ril_assembly/megahit/ \
--k-list 27,33,55,77,99 -t 40

#A3
r1=$(ls /mnt/nfs/bioinfdata/home/NIOO/xinyap/data/unmapped_below20/r1/* -m)
r2=$(ls -m /mnt/nfs/bioinfdata/home/NIOO/xinyap/data/unmapped_below20/r2/*)
se=$(ls -m /mnt/nfs/bioinfdata/home/NIOO/xinyap/data/unmapped_below20/singleton/*)

megahit -1 ${r1} -2 ${r2} -r ${se} \
-o /mnt/nfs/bioinfdata/home/NIOO/xinyap/data/unmapped_below20/megahit/ \
--k-list 27,33,55,77,99 -t 40


