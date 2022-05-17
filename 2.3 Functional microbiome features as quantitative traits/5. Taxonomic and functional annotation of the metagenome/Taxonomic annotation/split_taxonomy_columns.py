#!/usr/bin/env python3

import sys
from pathlib import Path
import pandas as pd
#import csv

def split_taxonomy(file):
    with file.open("r") as f:
        #next(f)
        tax_dict = {}
        for line in f:
            #print(line)
            line = line.strip()
            cols = line.split('\t')
            #print(cols)
            contig_id = cols[0]
            full_taxonomy = cols[1]
            #print(contig_id)
            #print(full_taxonomy)
            splited_taxonomy = full_taxonomy.split("|")
            k = splited_taxonomy[0][3:]
            p = splited_taxonomy[1][3:]
            c = splited_taxonomy[2][3:]
            o = splited_taxonomy[3][3:]
            f = splited_taxonomy[4][3:]
            g = splited_taxonomy[5][3:]
            s = splited_taxonomy[6][3:]
            tax_dict[contig_id] = {"kingdom": k, "phylum": p , "class": c, "order": o, "family": f, "genus": g, "species": s}
    return tax_dict


if __name__ == "__main__":
    output_taxonomy = split_taxonomy(Path(sys.argv[1]))
    #print(output_taxonomy)
    csv_filename = "10kb_taxonomy_splited.csv"
    df = pd.DataFrame(output_taxonomy)
    df.to_csv(csv_filename)
    #with open(csv_filename, 'w') as f:
     #   writer = csv.DictWriter(f, output_taxonomy.keys(), delimiter = "\t")
      #  writer.writeheader()
       # writer.writerow(output_taxonomy)
