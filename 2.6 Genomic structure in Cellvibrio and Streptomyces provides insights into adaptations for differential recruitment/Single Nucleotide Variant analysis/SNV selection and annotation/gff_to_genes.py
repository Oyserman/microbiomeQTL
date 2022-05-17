#!/usr/bin/env python3
"""
Script by Elmar van der Wijk, MSc thesis student at Bioinformatics Group, WUR

Parses out a GFF file (for example from Prokka output), to finish after the
##FASTA tag, and not include the scaffold information, so only the gene
annotations remains.

input:
    argv[1]: GFF file
    argv[2]: output filename

"""
from sys import argv


if __name__ == "__main__":
    with open(argv[2], "w") as outfile:
        with open(argv[1], "r") as infile:
            for line in infile.readlines():
                if not line.startswith(("##")):
                    outfile.write(line)
                else:
                    if line.startswith("##FASTA"):
                        break
                    else:
                        continue
