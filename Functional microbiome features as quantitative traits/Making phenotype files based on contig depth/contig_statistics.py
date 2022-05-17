#!/usr/bin/env python3
"""
Author: XinyaP
Script to calculate depth/coverage statistics from alignment bed file
"""

# import statements
from sys import argv
from pathlib import Path

# function
def get_cov(file_name):
    print("contig_id", "length", "total_depth", "avg_depth", "percentage_coverage", sep = "\t" )

    with file_name.open("r") as f:
        contigs = {}
        for line in f:
            col = line.split('\t')
            contig_id = col[0]
            #print(contig_id)
               #with open("/mnt/nfs/bioinfdata/ngs2/ME2/raaijmakers_group/RIL2019_analysis/assembly_evaluation/test/a3_rm_cov100.fa", "r") as f2:
                    #for line in f2:
                        #if contig_id in line:
                           # contig_length = line.split(' ')[3].strip().split('=')[1]
                            #contig_len = int(contig_length) 
            region = int(col[2]) - int(col[1])
            depth = int(col[3])
            if depth == 0:
                covered  = 0
            else:
                covered = region
            total_depth =covered * depth
            # print (covered)
            if contig_id not in contigs:
               contigs[contig_id] = {"length": region, "total_depth": total_depth,"avg_depth": 0, "percentage_coverage": covered}
            else:
               contigs[contig_id]["length"] += region
               contigs[contig_id]["total_depth"] += total_depth
               contigs[contig_id]["percentage_coverage"] += covered
              # contigs[contig_id]["avg_depth"] += total_depth/contig_len
              # contigs[contig_id]["percentage_coverage"] += covered/contig_len*100
       # print(contigs)
# write to file
        for contig in contigs:
            contig_len = contigs[contig]["length"]
           # print(contig_len)
            contigs[contig]["avg_depth"] = contigs[contig]["total_depth"]/contig_len
            contigs[contig]["percentage_coverage"] = contigs[contig]["percentage_coverage"]/contig_len*100

            contigs[contig]["avg_depth"] = round(contigs[contig]["avg_depth"], 2)
            contigs[contig]["percentage_coverage"] = round(contigs[contig]["percentage_coverage"], 2)
            stat = list(contigs[contig].values())
           # print(stat)
            output = '\t'.join(map(str, stat))
            print(contig+"\t"+output)

if __name__ == "__main__":
    get_cov(Path(argv[1]))
