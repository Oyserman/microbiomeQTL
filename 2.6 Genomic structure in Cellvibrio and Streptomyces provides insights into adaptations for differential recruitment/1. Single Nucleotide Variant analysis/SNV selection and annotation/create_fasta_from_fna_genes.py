#!/usr/bin/env python3
"""
Script by Elmar van der Wijk, MSc thesis student at Bioinformatics Group, WUR

Script to create a FASTA file from a genes.FNA file, created by Prodigal or similar
tools. The FASTA file can contain many genes, each having their own FASTA sequence.

input:
    argv[1]: text file with on the first line the name of the annotation, and
        on the second line a semi-colon separated list of Prokka gene IDs
    argv[2]: Prokka FNA file, as created by prokka_ffn_to_fna.py from a GFF an
        FFN file.
        format:
            >header # start # stop # strand (1 for forward, -1 for reverse) # attributes
    
            header: contig name with unique gene number X, e.g.
                NODE_456_cov_123_length_67890_X,
                unique Prokka gene identifier, e.g. BKPGNLEB_00001
                annotation, e.g. Aminomethyltransferase
    argv[3]: output filename
    
"""
from sys import argv
import os.path    

def parse_gene_file(gene_file):
    
    gene_list = []
    with open(gene_file, "r") as infile:
        i = 0
        for line in infile.readlines():
            line  = line.strip()
            if not line:
                continue
            else:
                if i == 0:
                    i += 1
                    annotation = line
                elif i == 1:
                    i += 1
                    gene_list = line.split(";")
                else:
                    print("WARNING: gene file contained more than 2 lines, "\
                          "please check")
    return [annotation, gene_list]


def parse_fna(fasta_file):
    """Parse fasta file and return each fasta header and sequence as dict
    
    par fasta_file: str: path and name of fasta file.
    
    return: dict: keys: str: fasta header without the ">", for each header found
                in the fasta file
            values: header: the full fasta header
                list: sequence belonging to fasta header, with each line a
                list element
    """
    
    assert os.path.exists(fasta_file), \
        f"Fasta file: {fasta_file} does not exist, please check the path."
    
    with open(fasta_file, "r") as infile:
        sequence = []
        header = ""
        prokka_gene_id = ""
        fasta_dict = {}
        for line in infile.readlines():            
            line  = line.strip()
            if not line:
                continue
            elif line.startswith(">"):
                if line.split(">")[-1]:   # test if fasta has a name
                    if prokka_gene_id:  # if header is already present, yield 
                        fasta_dict[prokka_gene_id] = [header, sequence]
                         
                    header = line.split(">")[-1].strip()  # new header
                    prokka_gene_id = line.split("#")[0].split(" ")[1]
                    sequence = []  # initialise new sequence
                else:
                    raise ValueError("Fasta header has no name")
            else:
                sequence += [line]
        fasta_dict[prokka_gene_id] = [header, sequence] # add final sequence
    return fasta_dict


def create_fasta_file(gene_list, fasta_fna_dict, output_name):
    
    new_fasta_dict = {}
    for gene_id in gene_list[1]:
        # header  = fasta_fna_dict[gene_id][0] # not used
        sequence = fasta_fna_dict[gene_id][1]
        new_fasta_dict[gene_id] = sequence
        
    with open(output_name, "w") as outfile:
        for gene_id, sequence in new_fasta_dict.items():
            outfile.write(f">{gene_id}\n")
            outfile.write("\n".join(sequence)+"\n")
                
    return output_name


if __name__ == "__main__":
    print("\n"+"-"*79)
    
    annotation_gene_list = parse_gene_file(argv[1])
    
    fna_dict = parse_fna(argv[2])
    
    print(f"\n Writing FASTA file to: \n{argv[3]}")
    
    create_fasta_file(annotation_gene_list, fna_dict, argv[3])
