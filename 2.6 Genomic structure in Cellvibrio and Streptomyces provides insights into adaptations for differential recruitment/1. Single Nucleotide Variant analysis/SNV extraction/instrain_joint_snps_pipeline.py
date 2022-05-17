#!/usr/bin/env python3
"""
Script by Elmar van der Wijk, MSc thesis student at Bioinformatics Group, WUR

Script to extract SNPs from a set of BAM files.
Python script filter_reads.py (Alex Crits-Christoph) filters reads based on
quality and ANI, using a reference fasta file. Reference fasta file will be 
created from initial reference fasta file containing all contigs and is 
filtered from the inStrain profile SNV output file, so only contigs that are
needed for the analysis are left. File naming occurs automatically.


Input:
    argv[1]: path to the directory containing all the BAM files for analysis
    argv[2]: name and path of file to contain the names of the BAM files to
        be used by inStrain. The argv[1] folder will be read and the names of
        the BAM files are written to the argv[2] filename. If the file already
        exists, by for example manual creation, this file will be used instead.
    argv[3]: path and name of the reference .fasta file or annotated reference
        .fna file. This will be used for filtering the reads.
    argv[4]: path and name of SNV file, as produced by inStrain profile or a 
        filtered version. Only the contigs present in this file will be kept in
        the reference file.
    argv[5]: path and name of the output directory. This will contain the tab-
        delimited file containing for each SNP for each metagenome the
        consensus SNP count and variation SNP count. This directory will also
        contain the newly created filtered fasta file, as well as the filtered
        BAM files used for the analysis.
    argv[6]: number of parallel processes to run, more will run faster but will
        be more intensive in server usage.
        
"""
from sys import argv
import subprocess
import os.path
import multiprocessing as mp

def snv_file_contigs(snv_file):
    """Parse contigs from an inStrain SNV output file
    
    par snv_file: path and name of inStrain SNV file. File is tab-delimited
        and contains contig and SNV position in first and second column.
        
    return: dict of unique contigs in SNV file with all SNV positions as values
    """
    contig_dict = {}
    with open(snv_file, "r") as infile:
        for line in infile.readlines():
            line = line.strip()
            if not line:
                continue
            elif line.startswith("scaffold"):
                continue
            else:
                contig_id = line.split("\t")[0]
                snv_pos = int(line.split("\t")[1])
                
                if not contig_id in contig_dict.keys(): 
                    contig_dict[contig_id] = []
                else:
                    contig_dict[line.split("\t")[0]] += [snv_pos]
    return contig_dict


def parse_fasta(fasta_file):
    """Parse fasta file and return each fasta header and sequence as dict
    
    par fasta_header: str: path and name of fasta file.
    
    return: dict: keys: str: fasta header without the ">", for each header found
                in the fasta file
            values: list: sequence belonging to fasta header, with each line a
                list element
    """
    
    assert os.path.exists(fasta_file), \
        f"Fasta file: {fasta_file} does not exist, please check the path."
    
    with open(fasta_file, "r") as infile:
        sequence = []
        header = ""
        fasta_dict = {}
        for line in infile.readlines():            
            line  = line.strip()
            if not line:
                continue
            elif line.startswith(">"):
                if line.split(">")[-1]:   # test if fasta has a name
                    if header:  # if header is already present, yield 
                        fasta_dict[header] = sequence
                         
                    header = line.split(">")[-1].strip()  # new header
                    sequence = []  # initialise new sequence
                else:
                    raise ValueError("Fasta header has no name")
            else:
                sequence += [line]
        return fasta_dict


def write_fasta_from_snv(fasta_contig_dict, snv_contig_dict, new_fasta):
    """Write new fasta file containing only contigs in SNV file
    
    par fasta_contig_dict: dict: keys: str: fasta header without the ">", for
                each header found in the fasta file
            values: list: sequence belonging to fasta header, with each line a
                list element
    par snv_contig_dict: dict of unique contigs in SNV file with all SNV
        positions as values
    par new_fasta: name and path of fasta file to be created. Will contain only
        contigs that are in the SNV file. Sequence is derived from fasta file.
        
    return: new_fasta
    
    """    
    with open(new_fasta, "w") as outfile:
        for contig_id, seq_list in fasta_contig_dict.items():
            if contig_id in snv_contig_dict.keys():
                outfile.write(">" + contig_id + "\n")
                outfile.write("\n".join(seq_list)+"\n")
            else:
                continue
    return new_fasta


def get_list_of_bam(data_dir, outfilename):
    """Write file of all BAM file paths in directory if not present
    
    par data_dir: str: path to folder containing all BAM files for 
        metaSNV analysis.
    par outfilename: str: path and name of file to write the list of BAM files
        to.
        
    return outfilename: str: path and name of file to which the list of BAM 
        files is written. If outfilename already exists the file is not written
        but only the name is returned.
        
     return cmd_out: str: list of BAM files and paths put in
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


def call_filter_reads(cmd):
    cmd_out = subprocess.check_output(cmd, shell=True)
    return cmd_out.decode("utf-8")


def pool_filter_reads(bam_sample_list, py_dir, fasta_ref, output_dir, \
                 processes, mismatch_threshold=0.95,\
                 min_mapq=-1, max_insert_length=3, min_insert_length=50):
    """Filter BAM files using filter_reads.py with specific cut-offs
    
    Call filter_reads.py (Alex Crits-Christoph) to filter BAM files using
    specific cut-offs. Log file is written to output directory.
    
    par bam_sample_list: list of paths and names of all BAM files to be
        filtered.
    par py_dir: str: path to directory containing the filter_reads.py script
    par fasta_ref: path and name of reference fasta file to which all BAM files
        are mapped.
    par output_dir: directory to which filtered BAM files will be written.
    par processes: int: number of parallel processes to run, more will run
        faster but will be more intensive in CPU usage.
    par mismatch_threshold: float: Minimum percentage identity of read pairs to
        consensus to use the reads, default = 0.95.
    par min_mapq: float: Minimum mapQ score of either read in a pair to use 
        that pair. Default = -1.
    par max_insert_length: float: Multiplier to determine maximum insert size
        between two reads, default is to use 3x median.
    par min_insert_length: Minimum insert size between two reads, default =50.

    return: filtered_bam_list: list of filtered BAM file paths and names
    return: cmd_out_list: list of outputs of command calling, for error testing
    """
    if not py_dir.endswith("/"):
        py_dir = py_dir + "/"
    if not output_dir.endswith("/"):
        output_dir= output_dir + "/"
    
    assert os.path.exists(f'{py_dir}filter_reads.py'),\
        "Python script 'filter_reads.py' is not present in current folder:\n"\
            f"{py_dir}"
    
    print(f"ANI threshold:\t\t{mismatch_threshold}\nMinimum mapQ:\t\t"\
          f"{min_mapq}\nMaximum insert length:\t\t{max_insert_length}x median"\
          f"\nMinimum insert length:\t\t{min_insert_length}")
    
    
    ref_core = fasta_ref.split("/")[-1]
    ref_core = ref_core.split(".fa")[0]
    
    filtered_bam_list = []
    cmd_list = []
    cmd_out_list = []
    for bam_file in bam_sample_list:
        bam_core = bam_file.split("/")[-1]
        bam_core = bam_core.split(".bam")[0]
        filtered_bam = f"{output_dir}{bam_core}_{ref_core}_filtered.bam"
        filtered_bam_list += [filtered_bam]
        
        if os.path.exists(filtered_bam):
            cmd_out_list += ["no run: file existed"]
            continue
        
        # original command
        cmd_list += [f"python3 {py_dir}filter_reads.py {bam_file} {fasta_ref} "\
            f"-m {mismatch_threshold} -q {min_mapq} -l {max_insert_length} "\
            f"-u {min_insert_length} -g {filtered_bam}"]
        
        # cmd_list += [f"python3 {py_dir}filter_reads.py {bam_file} {fasta_ref} "\
        # f"-g {filtered_bam}"] # no settings command
    
    pool = mp.Pool(processes)
    cmd_out_list = [pool.apply(call_filter_reads, args=([cmd])) for cmd in cmd_list]
    pool.close()
    
    with open(f"{output_dir}filter_reads_log.txt", "w") as outfile:
        for cmd in cmd_out_list:
            outfile.write(cmd+"\n")
    
    return filtered_bam_list, cmd_out_list
    

def samtools_index(filtered_bam_list, processes):
    """Index BAM files using samtools index
    
    par filtered_bam_list: paths and names of filtered BAM files.
    par processes: int: number of parallel processes to run, more will run
        faster but will be more intensive in CPU usage.
    
    return: cmd_out_list: list of outputs of command calling, for error testing
    return: indexed_bam_list: list of paths and names of indexed BAM files.
    """
    cmd_out_list = []
    indexed_bam_list = []
    for bam_file in filtered_bam_list:
        bam_core = bam_file.split(".bam")[0]
        bam_indexed_name = bam_core + ".bam.bai"
        indexed_bam_list += [bam_indexed_name]
        
        assert os.path.exists(bam_file),\
            f"BAM file: {bam_file} does not exists, did the filtering step "\
                f"succeed?"
        if os.path.exists(bam_indexed_name):
            print(f"Indexed BAM file already exists, "\
                  f"samtools index is not run:\n{bam_indexed_name}")
            cmd_out_list += ["not run: file existed"]
            continue
        
        cmd = f"samtools index {bam_file}"
        cmd_out = subprocess.check_output(cmd, shell=True)
        cmd_out_list += [cmd_out.decode("utf-8")]
        
    return indexed_bam_list, cmd_out_list


def run_joint_snps(joint_snps_path, bam_file, snv_file):
    """Extract SNPs from BAM file based on inStrain profile SNV file.
    
    Run joint_snps.py (by Alex Crits-Christoph) to extract SNP count for each
    SNP in the inStrain profile SNV output file. Consequently, parse output
    and return as table.
    
    par joint_snps_path: path of the joint_snps.py script
        (Alex Crits-Christoph)
    par bam_file: BAM file to have its SNPs extracted from SNV file.
    par snv_file: path and name of SNV file, as outputted by inStrain profile.
        All SNPs in this file will be extracted from the BAM files.
    
    return: output from joint_snps.py
    """
    cmd = f"python3 {joint_snps_path} {bam_file} {snv_file}"
    
    joint_snps_out_raw = subprocess.check_output(cmd, shell=True)\
        .decode("utf-8")
    
    js_table_out = []
    
    joint_snps_lines = joint_snps_out_raw.split("\n")
    for line in joint_snps_lines:
        line_elems = line.split("\t")
        js_table_out += [line_elems]
    
    return js_table_out


def pool_joint_snps(bam_list, snv_file, output_file, py_dir, processes):
    """Extract for BAM files SNP counts from inStrain SNV file.
    
    For each BAM file, joint_snps.py (by Alex Crits-Christoph) is run to
    extract SNP count for each SNP in the inStrain profile SNV output file.
    Output is concatenated in a single table. Uses parallel processing for each
    BAM file.
    
    par bam_list: list of all BAM files to have their SNPs extracted.
    par snv_file: path and name of SNV file, as outputted by inStrain profile.
        All SNPs in this file will be extracted from the BAM files.
    par output_file: path and name of output table.
    par py_dir: str: path to directory containing the joint_snps.py script
    par processes: int: number of parallel processes to run, more will run
        faster but will be more intensive in CPU usage.
        
    return: list containing row elements of final table as written to file.
    """
    if not py_dir.endswith("/"):
        py_dir = py_dir + "/"
    joint_snps_path = f'{py_dir}joint_snps.py'
    assert os.path.exists(joint_snps_path),\
        f"Python script 'joint_snps.py' is not present in current folder:\n"\
        f"{py_dir}"
    
    pool = mp.Pool(processes=processes)
    js_table_list = [pool.apply(run_joint_snps,\
                                args = (joint_snps_path, bam_file, snv_file))\
                     for bam_file in bam_list]
    pool.close()
    
    master_table = []
    for index, js_table in enumerate(js_table_list):
        for row, row_values in enumerate(js_table):
            if index == 0: # initialise table
                # print(row_values)
                master_table += [row_values[:4]]
                if row == 0:
                    master_table[row] += [row_values[4].split("/")[-1]]
                    master_table[row] += [row_values[5].split("/")[-1]]
                else:
                    master_table[row] += row_values[4:6]
            else:
                # add elem 4 and 5 of each row to each row of master
                if row == 0:
                    master_table[row] += [row_values[4].split("/")[-1]]
                    master_table[row] += [row_values[5].split("/")[-1]]
                else:
                    master_table[row] += row_values[4:6]
    
    with open(output_file, "w") as outfile:
        for row in master_table:
            line = "\t".join(row)
            outfile.write(f"{line}\n")
    return master_table
    


if __name__ == "__main__":
    print("-"*79)
    output_dir = argv[5]
    if not output_dir.endswith("/"):
        output_dir = output_dir + "/"
    if not os.path.exists(output_dir):
        subprocess.check_call(f"mkdir {output_dir}", shell=True)
    
    print(f"Output is put in:\n{output_dir}\n")
    
    print("-"*40)
    print(\
      "[1/5] Preparing and reference fasta based on SNP-containing contigs:\n")
    print(f"Fasta: {argv[3]}\nSNV file: {argv[4]}")
    
    fasta_contig_dict = parse_fasta(argv[3])
    snv_contig_dict = snv_file_contigs(argv[4])
    print(f"Number of contigs containing SNPs: {len(snv_contig_dict)}")
    
    snv_file_core = argv[3].split("/")[-1]
    snv_file_core = snv_file_core.split(".")[0]
    new_fasta = f"{output_dir}{snv_file_core}.fa"
    
    write_fasta_from_snv(fasta_contig_dict, snv_contig_dict, new_fasta)
    print(f"New fasta file is written to:\n{new_fasta}")
    
    
    print("-"*40)
    print("[2/5] Extracting BAM files from directory")
    
    bam_list_file, bam_sample_list = get_list_of_bam(argv[1], argv[2])
    
    
    print("-"*40)
    print("[3/5] Filtering BAM files\n")
    py_dir = "/".join(argv[0].split("/")[:-1])
    
    # original
    (filtered_bam_list, cmd_out_list) = \
        pool_filter_reads(bam_sample_list, py_dir, new_fasta,\
                          output_dir, int(argv[6]))
    # test
    # (filtered_bam_list, cmd_out_list) = \
    #     pool_filter_reads(bam_sample_list, py_dir, argv[3],\
    #                       output_dir, int(argv[6]), mismatch_threshold=0.94,\
    #                           min_insert_length=51)
    
    for i, cmd in enumerate(cmd_out_list):
        print(f"{filtered_bam_list[i]} {cmd}")
    
    
    print("-"*40)
    print("[4/5] Indexing BAM files: samtools index")
    
    (indexed_bam_list, samtool_index_cmd_out_list) = \
        samtools_index(filtered_bam_list, argv[5])
    # for i, cmd in enumerate(samtool_index_cmd_out_list):
    #     print(f"{indexed_bam_list[i]}: {cmd}")
    
    print("-"*40)
    print("[5/5] Extracting SNPs from individual BAM files using joint_snps.py")
    print(f"Number of parallel processes: {argv[6]}")
    joint_snps_output_file = f"{output_dir}joint_SNPs_output.tsv"
    
    master_table = pool_joint_snps(filtered_bam_list, argv[4],\
                                   joint_snps_output_file, py_dir,\
                                       int(argv[6]))
        
    print(f"Joint SNPs output file created: {joint_snps_output_file}")
    
    
    
    
    
    
    
    
    
