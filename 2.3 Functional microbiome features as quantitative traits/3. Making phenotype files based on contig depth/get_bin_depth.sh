for bin in /mnt/nfs/bioinfdata/home/NIOO/xinyap/data/bin_contig_name/*; do
	for f in /mnt/nfs/bioinfdata/ngs2/ME2/raaijmakers_group/RIL2019_analysis/coassembly_depths/*_contigs.csv; do
		n=${f%_contigs.csv}
		grep -wf $bin $f | awk -v var="${n##*/}" ' { total += $4 } END { print var, total/NR }' >> /mnt/nfs/bioinfdata/home/NIOO/xinyap/data/bin_depth/${bin##*/}
	done
done
