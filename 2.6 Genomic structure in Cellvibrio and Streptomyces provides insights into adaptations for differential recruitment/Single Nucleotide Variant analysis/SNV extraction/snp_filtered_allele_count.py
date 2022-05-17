#!/usr/bin/env python3
"""
Script by Elmar van der Wijk, MSc thesis student at Bioinformatics Group, WUR

Script to extract allelic SNP counts at SNP loci from BAM files. BAM files are
filtered according to inStrain standards and filtered BAM files are used for
allele count extraction.

Input:
    argv[1]: path to the directory containing all the BAM files for analysis
    argv[2]: name and path of file to contain the names of the BAM files to
        be used by inStrain. The argv[1] folder will be read and the names of
        the BAM files are written to the argv[2] filename. If the file already
        exists, by for example manual creation, this file will be used instead.
    argv[3]: path and name of the reference .fasta file or annotated reference
        .fna file. This will be used for filtering the reads.
    argv[4]: path and name of BED file: Contains a tab-separated list of
        contigs and SNP position. For each, depth all 4 nt alleles (A,C,G,T) is
        measured over the BAM files (argv[1&2])
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
import pysam


def expand_bed(bed_file_in, bed_file_out):
    """Expand each line of the bed file 4 times, once for each base
    
    par bed_file_in: bed file: tab-separated contig and SNP position.
    par bed_file_out: bed file, but for each line in the original file, 4 lines
        with a different base appended in tab-delimited manner.
    
    return: number of lines in bed file
    """
    line_counter = 0
    with open(bed_file_out, "w") as outfile:
        with open(bed_file_in, "r") as infile:
            for line in infile.readlines():
                line = line.strip()
                if not line:
                    continue
                else:
                    line_counter += 1
                    for base in ["A", "C", "G", "T"]:
                        outfile.write(f"{line}\t{base}\n")
    
    return line_counter                
    
    


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
            # cmd_out_list += ["no run: file existed"] # no false command should be added as this gives bugs later in the script
            continue
        
        cmd_list += [f"python3 {py_dir}filter_reads.py {bam_file} {fasta_ref} "\
            f"-m {mismatch_threshold} -q {min_mapq} -l {max_insert_length} "\
            f"-u {min_insert_length} -g {filtered_bam}"]
    
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
            # print(f"Indexed BAM file already exists, "\
            #       f"samtools index is not run:\n{bam_indexed_name}")
            cmd_out_list += [""]
            continue
        
        cmd = f"samtools index {bam_file}"
        cmd_out = subprocess.check_output(cmd, shell=True)
        cmd_out_list += [cmd_out.decode("utf-8")]
        
    return indexed_bam_list, cmd_out_list


def extract_snps(bam_bed_tuple):
    """Extract allele counts for each SNP
    
    par bam_bed_tuple: tuple:
        par bam_sample: filtered BAM file to have its SNP allele counts to be
            extracted.
        par bed_file_list: list of tuples: contig ID, SNP position, allele. There
            should be 4 alleles per position (ACGT).
    
    return: snp_count_list: list: allele count for each SNP allele of the bed
        file. First element is bam filename.
    """
    bam_sample = bam_bed_tuple[0]
    bed_file_list = bam_bed_tuple[1]
    
    print(bam_sample) # to see how it is running
    
    bam_core = bam_sample.split("/")[-1].split(".bam")[0]
    
    snp_count_list = [bam_core]
    sam_file = pysam.AlignmentFile(bam_sample)
    
    for snp_allele in bed_file_list:
        allele_count = 0
        for pileup_column in sam_file.pileup(snp_allele[0], start=snp_allele[1], stop=snp_allele[1]+1, truncate = True, max_depth=100000,
                                                        stepper = 'nofilter', compute_baq= True,
                                                        ignore_orphans = True, ignore_overlaps = True,
                                                        min_base_quality = 0):
                                
            for pileup_read in pileup_column.pileups:
                if not pileup_read.is_del and not pileup_read.is_refskip:
                    if pileup_read.alignment.query_sequence[pileup_read.query_position] == snp_allele[2]:
                        allele_count += 1
        snp_count_list += [allele_count]
    return snp_count_list


def extract_snps_parallel(bam_sample_list, bed_file_expanded, process):
    """Extract allele counts for each SNP
    
    par bam_sample_list: list of filtered BAM files to have their SNP allele
        counts to be extracted.
    par bed_file_expanded: bed file: tab-delimited, for each SNP, its position
        on contig and 4 alleles in order: ACGT, resulting in 4 lines per SNP.
    par process: int: number of parallel processes to run, more will run
        faster but will be more intensive in CPU usage.
    
    return: snp_total_count_list: list: for each row, the contig, SNP position
        and SNP allele, followed by the allele counts of each BAM file.
        First row is a header row, with contig, position, allele and BAM files
    """
    snp_total_count_list = [["contig", "position", "allele"]]
    
    bed_file_list = []
    
    with open(bed_file_expanded, "r") as bed_file:
        line_counter = 0
        for line in bed_file.readlines():
            line_counter += 1
            line = line.strip()
            if not line:
                continue
            else:
                snp_total_count_list += [[line]]
                contig = line.split("\t")[0]
                position = int(line.split("\t")[1])
                base = line.split("\t")[2]
                
                bed_file_list += [(contig, position, base)]    
    
    data_list = []
    for element in bam_sample_list:
        data_list += [(element, bed_file_list)]
    
    with mp.Pool(processes=process) as pool:
        snp_table = pool.map(extract_snps, iterable = data_list, chunksize = 1)
    
    for sample in snp_table:
        for i, element in enumerate(sample):
            snp_total_count_list[i] += [element]
    return snp_total_count_list


def write_count_table(output_file, bam_total_count_list, filter_bool):
    
    final_snp_counter = 0          
    with open(output_file, "w") as outfile:
        for i, row in enumerate(bam_total_count_list):
            if filter_bool and i != 0 and sum(row[1:]) == 0:
                continue
            else:
                line = "\t".join([str(n) for n in row])
                outfile.write(line + "\n")
                if i != 0:
                    final_snp_counter += 1
    return final_snp_counter
                    
            

if __name__ == "__main__":
    print("-"*79)
    output_dir = argv[5]
    if not output_dir.endswith("/"):
        output_dir = output_dir + "/"
    if not os.path.exists(output_dir):
        subprocess.check_call(f"mkdir {output_dir}", shell=True)
    
    print(f"Output is put in:\n{output_dir}\n")

    print("-"*40)
    print("[1/5] Extracting BAM files from directory")
    
    bam_list_file, bam_sample_list = get_list_of_bam(argv[1], argv[2])
    
    
    print("-"*40)
    print("[2/5] Filtering BAM files\n")
    py_dir = "/".join(argv[0].split("/")[:-1])
    
    (filtered_bam_list, cmd_out_list) = \
        pool_filter_reads(bam_sample_list, py_dir, argv[3],\
                          output_dir, int(argv[6]))

    
    # for i, cmd in enumerate(cmd_out_list):
    #     print(f"{filtered_bam_list[i]} {cmd}")
    
    
    print("-"*40)
    print("[3/5] Indexing BAM files: samtools index")
    
    (indexed_bam_list, samtool_index_cmd_out_list) = \
        samtools_index(filtered_bam_list, argv[5])

    print("\n" + "-"*40)
    print("[4/5] Preparing BED file\n")
    
    bed_file_in = argv[4]
    bed_in_core = bed_file_in.split(".")[0]
    bed_in_core = bed_in_core.split("/")[-1]
    bed_file_out = output_dir + f"{bed_in_core}_expanded.bed"
    
    bed_file_lines = expand_bed(bed_file_in, bed_file_out)
    print(f"Lines in bed file: {bed_file_lines}")
    print(f"Wrote expanded bed file to: {bed_file_out}")
    
    
    print("-"*40)
    print("[5/5] Extracting SNPs from individual BAM files - this might take a while")
    print(f"Number of parallel processes: {argv[6]}")
    
    snp_total_count_list = extract_snps_parallel(filtered_bam_list, bed_file_out, int(argv[6]))
    
    snp_table_file = f"{output_dir}{bed_in_core}_extracted_snp_allele_counts.tsv"
    
    snp_diff_count = write_count_table(snp_table_file, snp_total_count_list, True)
    
    
    print(f"Number of SNP alleles: {len(snp_total_count_list)-1}")
    print(f"Number of SNP alleles with any counts: {snp_diff_count}")
    print(f"Wrote output to {snp_table_file}")

    snp_table_file_full = f"{output_dir}{bed_in_core}_extracted_snp_allele_counts_full.tsv"
    
    snp_diff_count = write_count_table(snp_table_file_full, snp_total_count_list, False)
    
    print(f"Wrote full, unfiltered, output to {snp_table_file_full}")
    
    