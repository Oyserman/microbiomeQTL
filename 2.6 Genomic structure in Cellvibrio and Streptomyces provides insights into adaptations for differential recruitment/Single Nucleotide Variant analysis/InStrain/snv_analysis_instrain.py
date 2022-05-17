#!/usr/bin/env python
"""
Script by Elmar van der Wijk, MSc thesis student at Bioinformatics Group, WUR

Script to run InStrain on metagenome data for SNV analysis at different loci
of the assembled bacterial contis or genomes. Different samples are compared by
inStrain compare.
This pipeline only allows for comparison of BAM files, without extraction of
SNPs for the individual bam files making up the merged BAM file. For that, use
instrain_filter_reads_pipeline.py, which should be in the same folder as this
script, using the same input parameters.

2020/11/25: changed inStrain settings from default to MAPQ > 19

Input:
    argv[1]: path to the directory containing all the BAM files for analysis
    argv[2]: name and path of file to contain the names of the BAM files to
        be used by inStrain. The argv[1] folder will be read and the names of
        the BAM files are written to the argv[2] filename. If the file already
        exists, by for example manual creation, this file will be used instead.
    argv[3]: path and name of the reference .fasta file or annotated reference
        .fna file.
    argv[4]: output directory in which all the inStrain output will be put.
    argv[5]: number of parallel processes to run, more will run faster but will
        be more intensive in server usage.
    argv[6]: path and name of the gene annotation file. This file is based on
        the reference .fasta file and is created by for example Prodigal.
        If not present: "no_annotation".
    argv[7]: path and name of the scaffold-to-bin file (STB): Tab-delimited
        file indicating to which bin each fasta header in the reference genomee
        belongs.
        If not present: "no_stb".
     
    
    example usage:
        python3 snv_analysis_instrain.py [bam_folder/] [bam_file_selection.txt]
            [ref.fa / bin.x.fa] [inStrain_output/] 8 [ref.fa.genes.fna]
            [ref.fa.stb]
"""


from sys import argv
import subprocess
import os.path


def get_list_of_bam(data_dir, outfilename):
    """Write file of all BAM file paths in directory if not present
    
    par data_dir: str: path to folder containing all BAM files for 
        metaSNV analysis.
    par outfilename: str: path and name of file to write the list of BAM files
        to.
        
    return outfilename: str: path and name of file to which the list of BAM 
        files is written. If outfilename already exists the file is not written
        but only the name is returned.
        
     return cmd_out: str: list of newline-delimited BAM files and paths put in
         the file of outfilename.            
    """
    
    assert os.path.exists(data_dir), "BAM file folder does not exist,"\
        " please check path."
        
    if os.path.exists(outfilename): 
        with open(outfilename, "r") as df:
            bam_sample_list = []
            for line in df.readlines():
                line = line.strip()
                if not line:
                    continue
                elif line == "":
                    continue
                else:
                    bam_sample_list += [line]
        print("\nBAM list file already existed, was not changed:\n"\
              f"{outfilename}\n")      
        print("\nList of BAM files used for analysis:")
        for sample in bam_sample_list:
            print(sample)
        return outfilename, bam_sample_list
    
    
    # If file with all BAM files does not yet exist:  
    if data_dir.endswith("/"):
        data_dir = data_dir[:-1]
                   
    cmd = f"ls `pwd`/{data_dir}/*.bam"
    cmd_out = subprocess.check_output(cmd, shell=True)
    
    with open(outfilename, "w") as of:
        of.write(cmd_out.decode("utf-8"))
    
    print("\n> New BAM list file is written to:")
    print(f"{outfilename}")
    
    print("\nList of BAM files used for analysis:")
    print(f'\n{cmd_out.decode("utf-8")}')
    
    bam_sample_list = cmd_out.decode("utf-8").split("\n")[:-1]
    
    return outfilename, bam_sample_list
    

def call_instrain_profile(bam_file, fasta_ref, output_dir, processes, \
                          gene_annotation_file="no_annotation",\
                              stb_file="no_stb"):
    """Call inStrain profile to perform SNV calling
    
    par bam_file_name: str: path and name of BAM file to be analysed against
        reference genome. Name ends with ".bam".
    par fasta_ref: str: path and name of reference genome or bin in fasta
        format. The BAM file has to be mapped against this genome. Name ends
        with ".fa" or ".fasta".
    par output_dir: str: directory to contain the inStrain output directories
        of each individual analysis.
    par processes: int: number of parallel processes to run, more will run
        faster but will be more intensive in CPU usage.
    par gene_annotation_file: str: path and name of gene annotation file of the
        reference genome. Is for example created by Prodigal.
        If not present, put to "no_annotation".
    par stb_file: str: path and name of scaffold-to-bin file; a tab-delimited
        file with scaffold fasta header of the reference genome in the first
        column and the corresponding bin name in the second column.
        If not present, put to "no_stb".
    
    return: str: inStrain output directory path
    return: str: inStrain output
    """
    
    if not os.path.exists(output_dir):
        print("\ninStrain output directory does not exist, making new one:\n"\
              f"{output_dir}")
        cmd = f"mkdir {output_dir}"
        subprocess.check_call(cmd, shell=True)
    if not output_dir.endswith("/"):
        output_dir = output_dir + "/"

    bam_filename = bam_file.split(".bam")[0].split("/")[-1]
    ref_filename = fasta_ref.split(".fa")[0].split("/")[-1]
    
    sample_output_dir = f"{output_dir}{ref_filename}_{bam_filename}.IS"
    print(f"\nPutting inStrain output in:\n{sample_output_dir}")
    
    if gene_annotation_file == "no_annotation":
        gene_an_cmd = ""
    else:
        gene_an_cmd= f" -g {gene_annotation_file}"
    
    if stb_file == "no_stb":
        stb_cmd = ""
    else:
        stb_cmd = f" -s {stb_file}"
    
    cmd = f"inStrain profile {bam_file} {fasta_ref} -o {sample_output_dir}"\
        f" -p {processes}{gene_an_cmd}{stb_cmd}"
    cmd_out = subprocess.check_output(cmd, shell=True)
    
    return sample_output_dir, cmd_out


def call_instrain_compare(is_output_name_list, output_dir, processes):
    """Call inStrain compare to analyse/compare samples mapped to the same ref
    
    par is_output_list: list: list containing all paths to the inStrain output
        directories, usually having a dir.IS format.
    par output_dir: str: path and name of directory to contain the inStrain
        compare output. Directory will be names dir.IS.COMPARE    
    par processes: int: number of parallel processes to run, more will run
        faster but will be more intensive in CPU usage.
        
    return: inStrain compare output
    """
    
    for is_name in is_output_name_list:
        assert os.path.exists(is_name), f"inStrain output file {is_name} "\
            f"does not exist. Please check if inStrain profile ran correctly." 
    
    if not output_dir.endswith("/"):
        output_dir = output_dir + "/"
    output_dir_comp = output_dir+"IS.COMPARE"
    
    ver = ""
    while os.path.exists(output_dir_comp):
        print(f"\nOutput directory for inStrain compare already exists:\n"\
              f"{output_dir_comp}")
        if ver == "":
            ver = 1
        output_dir_comp = f"{output_dir}_{ver}/IS.COMPARE"
        ver += 1
    
    print(f"Writing inStrain compare output to:\n{output_dir_comp}")
    
    is_output_string = " ".join(is_output_name_list)
    cmd = f"inStrain compare -i {is_output_string} -o {output_dir_comp} "\
        f"-p {processes}"
    cmd_out = subprocess.check_output(cmd, shell=True)
    
    return cmd_out


if __name__ == "__main__":
    print("-"*79)
    
    try:
        gene_annotation_file = argv[6]
    except IndexError:
        gene_annotation_file = "no_annotation" 
        
        
    try:
        stb_file = argv[7]
    except IndexError:
        stb_file = "no_stb"
    
    (bam_list_file, bam_sample_list) = get_list_of_bam(argv[1], argv[2])
    
    is_output_name_list = []
    for b, bam_sample in enumerate(bam_sample_list):
        print("-"*79)
        print(f"\nAnalysing sample [{b+1}/{len(bam_sample_list)}] "\
              f"{bam_sample.split('/')[-1]}")
        (is_output_name, instrain_profile_output) = \
            call_instrain_profile(bam_sample, argv[3], argv[4], argv[5], \
                                  gene_annotation_file, stb_file)
        is_output_name_list += [is_output_name]
    
    assert len(is_output_name_list) > 0, \
        "No IS output was generated! Did everything run properly?"
    
    is_compare_out_dir = is_output_name_list[0].split("/")[:-1]
    is_compare_out_dir = "/".join(is_compare_out_dir)
    
    instrain_compare_output = call_instrain_compare(is_output_name_list, \
                                is_compare_out_dir, argv[5])
    
