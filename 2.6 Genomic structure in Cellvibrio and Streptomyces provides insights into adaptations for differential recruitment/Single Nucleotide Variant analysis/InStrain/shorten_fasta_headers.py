#!/usr/bin/env python3
"""
Script by Elmar van der Wijk, MSc thesis student at Bioinformatics Group, WUR

Script to shorten the contig headers of a fasta file. The contigs that are
recognised should be made by either metaSPAdes or MEGAHIT and should have the
following format:
    
    >NODE_11842_length_12527_cov_6.351304           (metaSPAdes)
    >k99_19237814 flag=0 multi=5.7582 len=11259     (MEGAHIT)

input:
    argv[1]: name and path of fasta file
    argv[2]: name and path of output fasta file

output:
    The fasta headers will be pruned so they are shorter, while still
    maintaining their unique numbers.

    >NODE_11842     (metaSPAdes)
    >k99_19237814   (MEGAHIT)
"""
from sys import argv
import os.path


def yield_fasta(fasta_file):
    """Parse fasta file and yield each fasta header and sequence
    
    par fasta_header: str: path and name of fasta file.
    
    yield: str: fasta header without the ">", for each header found in the
        file
        list: sequence belonging to fasta header, with each line a list element
    """
    
    assert os.path.exists(fasta_file), \
        f"Fasta file: {fasta_file} does not exist, please check the path."
    
    with open(fasta_file, "r") as infile:
        sequence = []
        header = ""
        for line in infile.readlines():            
            line  = line.strip()
            if not line:
                continue
            elif line.startswith(">"):
                if line.split(">")[-1]:   # test if fasta has a name
                    if header:  # if header is already present, yield 
                        yield (header, sequence)
                         
                    header = line.split(">")[-1].strip()  # new header
                    sequence = []  # initialise new sequence
                else:
                    raise ValueError("Fasta header has no name")
            else:
                sequence += [line]
        yield (header, sequence)
        

def prune_header(header):
    """Prune fasta header to be shorter
    
    par header: str: the fasta headers that are recognised should be made by 
    either metaSPAdes or MEGAHIT and should have the following format:
        >NODE_11842_length_12527_cov_6.351304           (metaSPAdes)
        >k99_19237814 flag=0 multi=5.7582 len=11259     (MEGAHIT)
    
    return: str: shorter header:
        >NODE_11842     (metaSPAdes)
        >k99_19237814   (MEGAHIT)
    """
    header_new = ""
    if header.startswith("NODE"):
        header_new = header.split("_length_")[0]
    elif header.startswith("k99"):
        header_new = header.split(" ")[0]
    
    if header_new == "":
        raise ValueError(f"{header} is of not supported format.")
    else:
        return header_new


if __name__ == "__main__":
    
    print("-"*79)
    print("Creating new fasta file with pruned headers!")
    
    print(f"Writing to {argv[2]}")
    with open(argv[2], "w") as outf:
        for header, seq in yield_fasta(argv[1]):
            header_new = prune_header(header)
            outf.write(f">{header_new}\n")
            seq_str = '\n'.join(seq)
            outf.write(f"{seq_str}\n")
        
        
