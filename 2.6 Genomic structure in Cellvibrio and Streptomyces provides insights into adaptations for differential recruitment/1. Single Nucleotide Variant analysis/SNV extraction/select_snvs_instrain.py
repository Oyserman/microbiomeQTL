#!/usr/bin/env python3
"""
Script by Elmar van der Wijk, MSc thesis student at Bioinformatics Group, WUR

Script to parse out an inStrain SNV output file under
    Analysis.IS/output/Analysis.IS_SNVs.tsv
so it contains only SNVs that are from a SNV selection file. For that, see the
analsysis script for SNV selection in R:
    SNV_selection_MM_P_analysis.R (Elmar van der Wijk)
The list of SNVs that are of interest should be in the format:
    NODE_5985_length_18358_cov_180.368531_pos_97, which is:
    contig_id + _pos_ + SNP position on contig

input:
    argv[1]: selected SNV file
    argv[2]: inStrain output file to have SNV data extracted
    argv[3]: output file to contain extracted SNV data from argv[2], but only
        those that occur in argv[1]

output:
    output file to contain extracted SNV data from argv[2], but only
        those that occur in argv[1]
"""
from sys import argv


def parse_snv_pos_file(snv_pos_file):
    """Parse selected SNV file and return list of contig_id and position
    
    par snv_pos_file: path and name of selected SNV file, containing a list of
        all contig_ids with positions in the following format:
            [contig_id]_pos_[snv_position]
    
    return: list of [contig_id]_pos_[snv_position]
    """
    
    selected_snv_list = []
    with open(snv_pos_file, "r") as infile:
        for line in infile.readlines():
            line = line.strip()
            if not line:
                continue
            else:
                selected_snv_list = line
    return selected_snv_list


def parse_instrain_snv_file(instrain_snv_file):
    """Parse inStrain SNV output file and create dict per contig and SNV pos
    
    par instrain_snv_file: name and path of inStrain SNV output file. Should
        be in tab-delimited format with contig_id and SNV position as first and
        second columns, respectively.
    
    return: dict: key: [contig_id]_pos_[snv_position], value: full line of SNV
        file as produced by inStrain.
    """
    instrain_snv_dict = {}
    with open(instrain_snv_file, "r") as infile:
        for line in infile.readlines():
            line = line.strip()
            if not line:
                continue
            else:
                contig_id = line.split("\t")[0]
                snv_pos = line.split("\t")[1]
                contig_id_snv_pos = f"{contig_id}_pos_{snv_pos}"
                instrain_snv_dict[contig_id_snv_pos] = line
    return instrain_snv_dict


def create_selected_snv_file(selected_snv_pos_list, instrain_dict, outname):
    """Write inStrain SNV file containing only selected_SNVs
    
    par selected_snv_pos_list: list: contains [contig_id]_pos_[snv_position]
    par instrain_dict: dict: key: [contig_id]_pos_[snv_position], value: full
        line of SNV file as produced by inStrain.
    par outname: path and name of output file: a shorter version of the
        inStrain SNV file, containing only SNVs from the selected set.
    
    return: outname
    """
    with open(outname, "w") as outfile:
        for contig_id_pos in selected_snv_pos_list:
            outfile.write(instrain_dict[contig_id_pos]+"\n")
    return outname


if __name__ == "__main__":
    print("-"*79)
    print(f"\nReading SNV selection file:\n{argv[1]}")
    selected_snv_contig_pos_list = parse_snv_pos_file(argv[1])
    
    print(f"\nReading inStrain SNV output file:\n{argv[2]}")
    instrain_snv_dict = parse_instrain_snv_file(argv[2])
    
    print(f"Writing SNV file with selected SNVs to\n{argv[3]}")
    print(create_selected_snv_file(selected_snv_contig_pos_list,\
                                   instrain_snv_dict, argv[3])
    
    
    