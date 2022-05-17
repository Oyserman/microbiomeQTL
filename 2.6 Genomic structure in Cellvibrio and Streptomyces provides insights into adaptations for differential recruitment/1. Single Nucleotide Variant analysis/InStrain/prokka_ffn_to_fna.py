#!/usr/bin/env python3
"""
Script by Elmar van der Wijk, MSc thesis student at Bioinformatics Group, WUR

Script to convert a Prokka FFN file into a FNA file that has the same format as
a Prodigal FNA file. For each fasta sequence header information from the GFF
file will be parsed and added in format:
    
    >header # start # stop # strand (1 for forward, -1 for reverse) # attributes
    
    header: contig name with unique gene number X, e.g.
        NODE_456_cov_123_length_67890_X,
        unique Prokka gene identifier, e.g. BKPGNLEB_00001
        annotation, e.g. Aminomethyltransferase
    
input:
    argv[1]: Prokka FFN file
    argv[2]: Prokka GFF file, shortened by gff_to_genes.py
        seq_name
        source
        feature
        start
        end
        score
        strand
        frame
        attribute
    argv[3]: Contig list file: contains the long contig name of all contigs
        present in the FFN and GFF file.    
    argv[4]: output filename, including FNA naming

"""
from sys import argv
import os.path


def parse_ffn(ffn_file):
    """Parse fasta file and return each fasta header and sequence as dict
    
    par fasta_header: str: path and name of fasta file.
    
    return: dict: keys: str: sequence ID
            values:
                0: annotation name
                1: list: sequence belonging to fasta header, with each line a
                    list element
    """
    
    assert os.path.exists(ffn_file), \
        f"ffn file: {ffn_file} does not exist, please check the path."
    
    with open(ffn_file, "r") as infile:
        annotation = ""
        sequence = []
        header = ""
        fasta_dict = {}
        for line in infile.readlines():            
            line  = line.strip()
            if not line:
                continue
            elif line.startswith(">"):
                if line.split(">")[-1]:   # test if fasta has a name
                    if header:
                        fasta_dict[header] = [annotation, sequence]
                         
                    header_line = ">".join(line.split(">")[1:]).strip()  # new header
                    header = header_line.split(" ")[0]
                    annotation = " ".join(header_line.split(" ")[1:])
                    sequence = []  # initialise new sequence
                else:
                    raise ValueError("Fasta header has no name")
            else:
                sequence += [line]
        fasta_dict[header] = [annotation, sequence]    
    return fasta_dict


def parse_gff(gff_file):
    """Parse GFF file and create dict with sequence information
    
    par gff_file: GFF file containing no FASTA sequences or other formates,
        only the GFF part
    
    return gff_dict:
        keys: unique sequence IDs
        values: short contig ID, sequence start position, sequence end position
            strand (+/-/.), semicolon-separated list of attributes
    
    """
    gff_dict = {}
    with open(gff_file, "r") as infile:
        for line in infile.readlines():
            line = line.strip()
            if not line or line.startswith("#"): # capture lines to be ignored
                continue
            else:
                line_elements = line.split("\t")
                try: # TEST
                    contig_short = line_elements[0]
                    start = line_elements[3]
                    end = line_elements[4]
                    strand = line_elements[6]
                    attribute = line_elements[8]
                    seq_id = attribute.split(";")[0]
                    seq_id = seq_id.split("ID=")[-1]
                except IndexError: # TEST
                    print("ERROR")
                    print(line_elements)
                    break
                
                gff_dict[seq_id] = [contig_short, start, end, strand, attribute]
    return gff_dict


def parse_contig_list(contig_list_file):
    """Parse list of contigs and return dict of short and long contig IDs
    
    par contig_list_file: file containing a list of long contig IDs in a
        newline-delimited format. Contig names have to be in SPAdes format to
        be recognised.
        
    return: dict: short contig names as values and long contig names as keys.
        Short contig names are cut off before the _length_ part in the name.
    """
    
    contig_dict = {}
    with open(contig_list_file, "r") as infile:
        for line in infile.readlines():
            line = line.strip()
            if not line:
                continue
            else:
                contig_id_long = line
                contig_id_short = line.split("_length_")[0]
                contig_dict[contig_id_short] = contig_id_long
    return contig_dict
                
                
def create_fna(ffn_dict, gff_dict, long_contig_name_dict, fna_filename):
    """Create FNA file from ffn and GFF data
    
    FNA file will be created that has FASTA format with headers as followed:
        contig ID sequence ID annotation # start # end # strand # GFF attributes
    
    par ffn_dict:
        keys: unique sequence IDs
        values: annotation name and list of each line of the fasta sequence
    par gff_dict:
        keys: unique sequence IDs
        values: short contig ID, sequence start position, sequence end position
            strand (+/-/.), semicolon-separated list of attributes
    par long_contig_name_dict: short contig names as keys, and long contig
        names as values.
    par fna_filename: name of output FNA file:
    
    
    """
    
    contig_counter_dict = {}
    
    with open(fna_filename, "w") as outfile:
        for seq_id, data in ffn_dict.items():
            annotation = data[0]
            sequence = data[1]
            
            try: # test if key is valid
                gff_data = gff_dict[seq_id]
            except KeyError:
                print("ERROR")
                print(seq_id)
                print(annotation)
           
            contig_id_short = gff_data[0]
            start = gff_data[1]
            end = gff_data[2]
            strand = gff_data[3]
            attributes = gff_data[4]
            contig_id_long = long_contig_name_dict[contig_id_short]
            
            if strand == "+":
                strand = 1
            elif strand == "-":
                strand = -1
            else:
                strand = ""
            
            if not contig_id_long in contig_counter_dict.keys():
                contig_counter_dict[contig_id_long] = 0
            contig_counter_dict[contig_id_long] += 1
            counter = contig_counter_dict[contig_id_long]
            
            header = f">{contig_id_long}_{counter} {seq_id} {annotation} "\
                f"# {start} # {end} # {strand} # {attributes}\n"
            outfile.write(header)
            outfile.write("\n".join(sequence)+"\n")
    return fna_filename


if __name__ == "__main__":
    print("-"*79)
    print(f"\nExtracting fasta sequences from ffn file:\n{argv[1]}")
    ffn_dict = parse_ffn(argv[1])

    print(f"\nExtracting GFF file data:\n{argv[2]}\n")
    gff_dict = parse_gff(argv[2])
    
    print(f"\nExtracting long contig names from:\n{argv[3]}")
    contig_dict = parse_contig_list(argv[3])

    print(f"Writing new FNA file:\n{argv[4]}")
    
    create_fna(ffn_dict, gff_dict, contig_dict, argv[4])
     
    






    
    