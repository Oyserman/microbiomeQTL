#!/usr/bin/env python3
"""
Script by Elmar van der Wijk, MSc thesis student at Bioinformatics Group, WUR

Script to convert a Prokka FNA file, created using prokka_ffn_to_fna.py into a
tab-delimited list, containing information on the fasta headers:
    [contig and gene counter]    [prokka gene id]    [prokka annotations]

input:
    argv[1]: FNA file as created by prokka_ffn_to_fna.py from Prokka output
    argv[2]: output filename
"""


from sys import argv
import os.path


def parse_fna(fna_file):
    """Parse FNA file and return list of first three fasta header elements
    
    par fasta_header: str: path and name of fasta file.
    
    return: list: for each fasta header one element, consisting of 3 elements:
        [[contig and gene counter]    [prokka gene id]    [prokka annotations]]
    """
    
    assert os.path.exists(fna_file), \
        f"Fasta file: {fna_file} does not exist, please check the path."
    
    with open(fna_file, "r") as infile:
        header = ""
        prokka_gene_id = ""
        header_list = []
        for line in infile.readlines():            
            line  = line.strip()
            if not line:
                continue
            elif line.startswith(">"):
                if line.split(">")[-1]:   # test if fasta has a name
                
                    header = line[1:].strip()  # remove the > 
                    contig_id_counter = header.split("#")[0].split(" ")[0]
                    prokka_gene_id = header.split("#")[0].split(" ")[1]
                    annotation = " ".join(header.split("#")[0].split(" ")[2:-1])
                    
                    header_list += [[contig_id_counter, prokka_gene_id, annotation]]
                else:
                    raise ValueError("Fasta header has no name")

    return header_list


def write_table(header_list, output_filename):
    """Write tab-delimited file of FNA gene fasta header information
    
    par header_list: list:
        for each fasta header one element, consisting of 3 elements:
        [[contig and gene counter]    [prokka gene id]    [prokka annotations]]
    par output_filename: name and path of output file.    
    """
    
    with open(output_filename, "w") as outfile:
        for header in header_list:
            outfile.write(f"{header[0]}\t{header[1]}\t{header[2]}\n")
    return output_filename
    

if __name__ == "__main__": 
    print("\n"+"-"*79)
    print(f"\nReading FNA file:\n{argv[1]}")
    
    header_list = parse_fna(argv[1])
    
    print(f"\nWriting tab-delimited file:\n{argv[2]}")
    
    write_table(header_list, argv[2])
    
    
    