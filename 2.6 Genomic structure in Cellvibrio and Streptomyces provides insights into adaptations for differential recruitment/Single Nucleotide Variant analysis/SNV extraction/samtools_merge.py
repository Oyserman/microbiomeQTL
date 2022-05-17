"""
Script by Elmar van der Wijk, MSc thesis student at Bioinformatics Group, WUR

Script to merge a list of .bam files into one .bam file

input:
    argv[1]: directory containing the .bam files
    argv[2]: file containing the paths and names of the .bam files to be
        merged. If the file does not exists, it is created and will contain all
        .bam files of the argv[1] directory.
    argv[3]: name and path of merged .bam file.
    argv[4]: number of parallel processes to run, more will run
        faster but will be more intensive in CPU usage.
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
              f"{outfilename}")      
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


def samtools_merge(bam_list_file, out_bam, processes):
    """Call samtools merge to merge BAM files in a new file
    
    par bam_list_file: str: path and name of file containing all the BAM files
        to be merged, in newline-delimited format.
    par out_bam: str: name and path of merged BAM file, consisting of all
        the BAM files in bam_list_file.
    par processes: int: number of parallel processes to run, more will run
        faster but will be more intensive in CPU usage.
    
    return out_bam: str: name and path of merged BAM file, consisting of all
        the BAM files in bam_list_file.
    return cmd_out: str: output of samtools merge command.    
    """
    if os.path.exists(out_bam):
        print("\nMerged BAM file already exists, no merge was performed:\n"\
              f"{out_bam}\nRemove .bam file if you want to merge.")
        return out_bam, "Samtools merge was not run."
    
    cmd = f"samtools merge {out_bam} -b {bam_list_file} --threads {processes}"
    cmd_out = subprocess.check_output(cmd, shell=True)
    cmd_out = cmd_out.decode("utf-8")
    print(f"\nSamtools output:\n{cmd_out}")
    print("\nSamtools merge was performed successfully, merged .bam file was "\
          f"created:\n{out_bam}")
    
    return out_bam, cmd_out


if __name__ == "__main__":
    print("-"*79)
    bam_list_file, bam_sample_list = get_list_of_bam(argv[1], argv[2])
    
    (out_bam, samtools_out) = samtools_merge(bam_list_file, argv[3], argv[4])
    
    
    