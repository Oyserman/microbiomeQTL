#!/usr/bin/env python3
"""
Script by Elmar van der Wijk, MSc thesis student at Bioinformatics Group, WUR

Script to write from a Kraken2 output file all contigs to a fasta file, based
on a reference fasta file. Kraken2 output file contains list of contig
names not in fasta format. Only the first element of a contig header of Kraken2
is taken.

Example:
    "example_fasta_header species name" becomes
    "> example_fasta_header"
    ATGGGACAGATACAAGATGACACAGATATGGCGACCGATG

input:
    argv[1]: Kraken2 output file containing all contig names belonging to a 
        feature of interest in newline-delimited manner.
    argv[2]: reference assembly genome, containing all contigs that comprise
        the Kraken2 output file.
    argv[3]: Output filename to contain fasta sequences of the contigs in both
        the bin file and kraken2 output. Ends with .fa or .fasta extension.
    argv[4]: minimum contig length, is used for filtering. Put to 0 for no
        filter.

output:
    File containing the unique contigs of the Kraken2 output file, in fasta 
    format, thus containing also the sequence and header.
"""
from sys import argv
import os.path
import copy


def parse_kraken2(kraken2_file):
    """Parse kraken2 output file and yield lists of contig names
    
    par kraken2_file: str: path and name of kraken2 output file containing
        contig names without ">" with newline-delimitation.
    
    yield: str: header, only first element if header contains spaces
    """
    assert os.path.exists(kraken2_file), \
        f"File: {kraken2_file} does not exist, please check the path."
    
    with open(kraken2_file) as infile:
        for line in infile.readlines():
            line = line.strip()
            if not line:
                continue
            else:
                header = line.split("\t")[0]
                header = header.split(" ")[0]
                yield header


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
        
        
def filter_contigs(contig_dict, minimum_length):
    """Filter contigs based on minimum length
    
    par contig_dict: dict: key: contig header, elements: sequence,
        sequence length
    par minimum_length: int: minimum length of the contigs after filtering.
    
    return: list of contigs that pass length requirement    
    """
    contig_thr_list = []
    contig_error_list = []
    contig_low_list = []
    
    for header, data in contig_dict.items():
        try:
            length = data[1]
        except:
            contig_error_list += [header]
        
        if length >= minimum_length:
            contig_thr_list += [header]
        else:
            contig_low_list += [header]
            
    print(f"\nNumber of contigs before filter: {len(contig_dict)}")
    print(f"Number of contigs containing errors: {len(contig_error_list)}")
    print(f"Number of contigs not passing threshold: {len(contig_low_list)}")
    print(f"Number of contigs after filter: {len(contig_thr_list)}")
    
    return contig_thr_list


def create_bin(contig_dict, out_filename):
    """Write a new bin.fasta file from headers and reference.fasta
    
    par contig_dict: dict: key: contig header, elements: sequence,
        sequence length
    par out_filename: str: name and path of output bin filename to contain all
        headers and corresponding sequences. Ends with .fa or .fasta extension.
        
    return: out_filename
    """
    ver = ""
    out_filename_core = "".join(out_filename.split(".fa")[:-1])
    if out_filename.split(".fa")[-1] == "sta":
        extension = "fasta"
    else:
        extension = "fa"
    out_filename_new = f"{out_filename_core}.{extension}"
    while os.path.exists(out_filename_new):
        if ver == "":
            print(f"\nOuput file already exists:"\
                  f"\n{out_filename}")
            ver = 0
        ver += 1
        out_filename_new = f"{out_filename_core}_{ver}.{extension}"
 
    with open(out_filename_new, "w") as outf:      
        for header, data in contig_dict.items():
            header = ">" + header
            seq = "\n".join(data[0])
            outf.write(f"{header}\n")
            outf.write(f"{seq}\n\n")
    return out_filename_new


if __name__ == "__main__":
    print(79*"-")
    print("[1/4] Parsing Kraken2 file")
    kraken_headers = list(parse_kraken2(argv[1]))
    print(kraken_headers[0:10])
    
    print("\n"+40*"-")
    print("[2/4] Retrieving fasta sequences from reference:")
    print(argv[2])
    
    kraken_header_dict = {}
    for ref_header, ref_seq in yield_fasta(argv[2]):
        ref_header_short = ref_header.split(" ")[0] #short is identical to kraken2
        if ref_header_short in kraken_headers:
            kraken_header_dict[ref_header] = [ref_seq, len("".join(ref_seq))]
    
    print("\n"+40*"-")
    print(f"[3/4] Filtering contigs on minimum size: {argv[4]}") 
    print(f"\n{argv[1]}")
    kraken_thr_list = filter_contigs(kraken_header_dict, int(argv[4]))
    
    kraken_thr_dict = {}
    for header in kraken_thr_list:
        kraken_thr_dict[header] = kraken_header_dict[header]
    
    print("\n"+40*"-")
    print("[4/4] Creating new bin file")     
    outfilename = create_bin(kraken_thr_dict, argv[3])
    print(f"\nNew bin file written to: {outfilename}")
    
    