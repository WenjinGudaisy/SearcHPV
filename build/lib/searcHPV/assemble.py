from searcHPV.generate_assemble import *
import subprocess
import os

############
#Function: run pipeline for assemble
#fq1: raw sequencing data, fastq file
#fq2: raw sequencing data, fastq file
#out_dir: output directory for searcHPV
#virRef: virus reference genome
#window: the length of region searching for informative reads, default=300
def assemble(fq1, fq2, out_dir,virRef,gz,window):
    bam = f'{out_dir}/alignment/alignment.RG.indelre.mkdup.sort.bam'
    check_file(bam)
    assemble_out_dir = f'{out_dir}/assemble/'
    mkdir(assemble_out_dir)
    extract_read_name(bam,out_dir,virRef,window)
    script_read_seq = extract_read_seq(assemble_out_dir,fq1,fq2,gz)
    os.system(f'chmod +x {script_read_seq}')
    subprocess.call(script_read_seq)
    
    #preprocessForPear(assemble_out_dir)
    script_pear = pear(assemble_out_dir)
    os.system(f'chmod +x {script_pear}')
    subprocess.call(script_pear)

    preprocessForCap3(assemble_out_dir)
    script_cap3 = cap3(assemble_out_dir)
    os.system(f'chmod +x {script_cap3}')
    subprocess.call(script_cap3)

    return None
