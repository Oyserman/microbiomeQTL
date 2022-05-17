#!/usr/bin/env python3
"""
Script by Elmar van der Wijk, MSc thesis student at Bioinformatics Group, WUR

Script to create from a Kraken2 output file a new scaffold-to-bin (stb) file.
This file will have for each fasta header in the reference genome a taxonomic
annotation. Kraken2 file is tab delimited:
    
    0               1       2       3       4       5       6       7
    contig_id       kingdom phylum  class   order   family  genus   species

input:
    argv[1]: Kraken2 output file containing all contig names belonging to a 
        feature of interest followed by taxonomic levels in tab-delimited
        manner.
    argv[2]: Output filename of the stb file. Will contain two tab-delimited
        rows: contig_id and taxonomy of interest.
    argv[3]: Taxonomic level of interest. See above. For example, if it is put
        to 7, the species level is kept.
    argv[4]: Bin file in fasta format containing all contigs of which the
        headers should be in the stb file.

output:
    File containing contig ID followed by the selected taxonomic level,
    tab-delimited.
"""

from sys import argv
import os.path


def parse_kraken_taxonomy(kraken_taxonomy_file):
    """Parse kraken2 output file and yield lists of contig names
    
    par kraken2_file: str: path and name of kraken2 output file containing
        contig names without ">" with newline-delimitation.
    
    return:
        taxonomy_dict: dict: key: contig_id, values: list of taxonomy
            annotations.
        taxonomy_levels: list: all taxonomy name annotations present in the
            first line of the file: kingdom, phylum, ... , species.
    """
    assert os.path.exists(kraken_taxonomy_file), \
        f"File: {kraken_taxonomy_file} does not exist, please check the path."
    
    taxonomy_dict = {}
    taxonomy_levels = []
    with open(kraken_taxonomy_file) as infile:
        for line in infile.readlines():
            line = line.strip()
            if line.startswith("contig_id"):
                taxonomy_levels = line.split("\t")
            elif not line:
                continue
            else:
                elements = line.split("\t")
                taxonomy_dict[elements[0]] = elements[1:]
    return taxonomy_dict, taxonomy_levels


def filter_taxonomy(taxonomy_dict, taxonomy_levels, filter_taxonomy):
    """Filter taxonomies to keep only taxonomy of interest
    
    Example: keep only genus taxonomy level.
    
    par taxonomy_dict: dict: key: contig_id, values: list of taxonomy
            annotations.
    par taxonomy_levels: list: all taxonomy name annotations present in the
            first line of the file: kingdom, phylum, ... , species.
    par filter_taxonomy: str: name of taxonomy to be kept. Must occur in
        taxonomy_levels.
    
    return: dict: key: contig_id, value: taxonomy level of interest value.
    """
    print(taxonomy_levels)
    taxonomy_indices = [i for i, elem in enumerate(taxonomy_levels)\
                        if elem == filter_taxonomy] 
    
    assert len(taxonomy_indices) == 1,\
        "Multiple taxonomic levels have the same name (kingdom, phylum, etc."\
            "\nPlease check taxonomy file."
    tax_index = taxonomy_indices[0] -1 # taxonomy_levels contains contig_id
    
    filtered_dict = {}
    
    for contig_id, values in taxonomy_dict.items():
        filtered_dict[contig_id] = values[tax_index]
    
    return filtered_dict


def yield_fasta_headers(fasta_file):
    """Parse fasta file and yield each fasta header without the ">"
    
    par fasta_header: str: path and name of fasta file.
    
    yield: dict:    index spades: Header if spades format, None otherwise.
                    index megahit: Header if megahit format, None otherwise.
    """
    
    assert os.path.exists(fasta_file), \
        f"File: {fasta_file} does not exist, please check the path."
    
    with open(fasta_file, "r") as infile:
        
        for line in infile.readlines():
            line  = line.strip()
            if not line:
                continue
            elif line.startswith(">"):
                header = line.split(">")[-1].strip()
                if header:   # test if fasta has a name
                    yield header       
            else:
                continue


def write_stb(filtered_tax_dict, outfilename): 
    """Write tab-delimited scaffold-to-bin file.
    
    par filtered_tax_dict: dict: key: contig_id, value: taxonomy level value.
    par outfilename: str: path and name of output file.
    
    return: outfilename
    """
    outfilename_core = ".".join(outfilename.split(".")[:-1])
    outfilename_new = outfilename
    ver = ""
    while os.path.exists(outfilename_new):
        if ver == "":
            print(f"\nOutput file {outfilename} existed, creating new file.")
            ver = 1
        else:
            ver += 1
        outfilename_new = f"{outfilename_core}_{ver}.stb"
    
    print(f"\nWriting scaffold-to-bin file to:\n{outfilename_new}\n")
      
    with open(outfilename_new, "w") as outf:
        for contig_id, tax in filtered_tax_dict.items():
            outf.write(contig_id+"\t"+tax+"\n")
        
    return outfilename
    
    

if __name__ == "__main__":
    print("-"*79)
    
    (taxonomy_dict, taxonomy_levels) = parse_kraken_taxonomy(argv[1])
    
    filtered_tax_dict = \
        filter_taxonomy(taxonomy_dict, taxonomy_levels, argv[3])
    
    print(f"\nNumber of contig IDs: {len(taxonomy_dict)}")
    print(f"Filtering on taxonomy level: {argv[3]}")
    
    bin_headers_list = yield_fasta_headers(argv[4])
    
    bin_filtered_tax_dict = {}
    for header in bin_headers_list:
        bin_filtered_tax_dict[header] = filtered_tax_dict[header]
    
    write_stb(bin_filtered_tax_dict, argv[2])
    
    
    