#!/usr/bin/env python3
"""
Script by Alex Crits-Christoph, University of California, Berkeley
Adaptations by Elmar van der Wijk, Bioinformatics group, Wageningen University

changes:
    varBase --> var_base
    refBase --> ref_base
"""

import pysam
import sys
import pandas as pd
import argparse
from Bio import SeqIO
import numpy as np
from collections import defaultdict 

P2C = {'A':0, 'C':1, 'T':2, 'G':3}
C2P = {0:'A', 1:'C', 2:'T', 3:'G'}

def joint_snps(bams, snvs):
    snvs = pd.read_csv(snvs, sep="\t")
    snv_freqs = defaultdict(dict)
    for bam in args.bams.split(","):
        bam = bam.strip()

        samfile = pysam.AlignmentFile(bam)
        for index, row in snvs.iterrows():
            if row['allele_count'] == 2:
                var_count = 0
                con_count = 0

                for pileupcolumn in  samfile.pileup(row['scaffold'], start=row['position'], stop=row['position']+1, truncate = True, max_depth=100000,
                                                    stepper = 'nofilter', compute_baq= True,
                                                    ignore_orphans = True, ignore_overlaps = True,
                                                    min_base_quality = 30):
                    table = np.zeros(4, dtype=int)
                    for pileupread in pileupcolumn.pileups:
                        if not pileupread.is_del and not pileupread.is_refskip:
                            if pileupread.alignment.query_sequence[pileupread.query_position] == row['var_base']:
                                var_count += 1
                            elif pileupread.alignment.query_sequence[pileupread.query_position] == row['con_base']:
                                con_count += 1

                snv_freqs[(row['scaffold'], row['position'])][bam] = (var_count, con_count)

    print("scaffold\tposition\tvar_base\tcon_base", end="\t")
    for bam in args.bams.split(","):
        print(bam.strip(".bam") + "_var_base" + "\t" + bam.strip(".bam") + "_con_base", end="\t")

    print("")
    for index, row in snvs.iterrows():
        if row['allele_count'] == 2:
            print(row['scaffold'] + "\t" + str(row['position']) + "\t" + row['var_base'] + "\t" + row['con_base'], end="\t")
            for b in args.bams.split(","):
                print(str(snv_freqs[(row['scaffold'],row['position'])][b.strip()][0]) + "\t" + str(snv_freqs[(row['scaffold'],row['position'])][b.strip()][1]), end="\t")
            print()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Joint-snp calling",
             formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    # Required positional arguments
    parser.add_argument("bams", help="List of BAMs, comma separated")
    parser.add_argument("snvs", help="InStrain SNVs file.")
    args = parser.parse_args()
    run = joint_snps(args.bams, args.snvs)