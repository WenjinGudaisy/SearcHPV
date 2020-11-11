import sys
import os
import re
from general import *
if len(sys.argv) != 3:
    raise ValueError("Please input a configuration file and parameter for remove duplicate")


#config Dictionary formate:key1 = "Input File",key2 = "reference_file = ",value = full_path
        
#####################
##generate alignment bash file
#Function: generate alignment bash file for a sample
#Parameters: bash_file: output bash file
#ref: reference gene
#fq1: fastq1 file
#fq2: fastq2 file
#out_dir: outputPath
#sample: samle name
#samtools: full path of samtools
#bwa: full path of bwa
def generate_alignment_bash(bash_file,ref,fq1,fq2,out_dir,sample,samtools,bwa):
    with open(bash_file,'w') as output:
        output.write(f'''
{bwa} mem {ref} '<zcat {fq1}' '<zcat {fq2}' > {out_dir}/{sample}.sam
{samtools} view -bhS {out_dir}/{sample}.sam > {out_dir}/{sample}.bam
{samtools} sort -o {out_dir}/{sample}.sort {out_dir}/{sample}.bam
rm {out_dir}/{sample}.sam
echo \'alignment done\'
''')
    return None

 ################


 ################
 ##generate indel re alignment bash file
#Parameters: bash_file: output bash file
#ref: reference gene
#fq1: fastq1 file
#fq2: fastq2 file
#out_dir: outputPath
#sample: samle name
#samtools: full path of samtools
#bwa: full path of bwa
#picard: full path of picard
#gatk: full path of GATK
#java: full path of java
def generate_indel_alignment_bash(bash_File,ref,sample,java,picard,out_dir,samtools,gatk):
    with open(bash_File,'w') as output:
        output.write(f'''{java} -Xmx4g -jar {picard} \
AddOrReplaceReadGroups \
I={out_dir}/{sample}.bam \
O={out_dir}/{sample}.RG.bam \
RGID={sample} \
RGPL=illumina \
RGPU=NA \
RGSM={sample} \
RGLB={sample}
{samtools} sort -o {out_dir}/{sample}.RG.sort.bam {out_dir}/{sample}.RG.bam
{samtools} index {out_dir}/{sample}.RG.sort.bam
{java} -Xmx4g -jar {gatk} \
-T RealignerTargetCreator \
-R {ref} \
-I {out_dir}/{sample}.RG.sort.bam \
-o {out_dir}/{sample}.RG.intervals
{java} -Xmx4g -jar {gatk} \
-T IndelRealigner \
-R {ref} \
-I {out_dir}/{sample}.RG.sort.bam \
-targetIntervals {out_dir}/{sample}.RG.intervals \
-o {out_dir}/{sample}.RG.indelre.bam
{samtools} index {out_dir}/{sample}.RG.indelre.bam
echo \'indel alignment done\'''')



def generate_mkdup_bash(bash_File,sample,java,picard,out_dir):
    with open(bash_File,'a') as output:
        output.write(f'''
samtools sort -n -o {out_dir}/{sample}.RG.indelre.sortbyQ.bam {out_dir}/{sample}.RG.indelre.bam
{java} -Xmx4g -jar {picard} MarkDuplicates \
I={out_dir}/{sample}.RG.indelre.sortbyQ.bam \
O={out_dir}/{sample}.RG.indelre.mkdup.bam \
M={out_dir}/{sample}.RG.indelre.mkdup.txt \
TAGGING_POLICY=All ASSUME_SORT_ORDER=queryname
samtools sort -o {out_dir}/{sample}.RG.indelre.mkdup.sort.bam {out_dir}/{sample}.RG.indelre.mkdup.bam
{samtools} index {out_dir}/{sample}.RG.indelre.mkdup.sort.bam
echo \'indel alignment done\'''')
 
 

############################
##Process
##read parameters
config_file = sys.argv[1]
config_dic = read_config(config_file)
input_para = "Input File"
tools = "Full Path to third party tools"
output_para = "Output"
out_dir = config_dic[output_para]['output_directory']
script_dir = config_dic[output_para]['script_directory']
ref = config_dic[input_para]["hg_virus_reference"]
check_file(ref)
fq1 = config_dic[input_para]["fastq1"]
# check_file(fq1)
fq2 = config_dic[input_para]["fastq2"]
# check_file(fq2)
sample = config_dic[output_para]["output_sample_name"]
samtools = config_dic[tools]["samtools"]
bwa = config_dic[tools]["bwa"]
java = config_dic[tools]["java"]
picard = config_dic[tools]["picard"]
gatk = config_dic[tools]["gatk"]


##generate alignment bash
# print(out_dir)
if not os.path.exists(out_dir):
    os.system(f"mkdir {out_dir}")
out_dir = f'{out_dir}/alignment'
if not os.path.exists(out_dir):
    os.system(f"mkdir {out_dir}")

script_dir = f'{script_dir}/alignment'
if not os.path.exists(script_dir):
    os.system(f"mkdir {script_dir}")

alignment_File = script_dir + f"/{sample}.alignment.sh"
generate_alignment_bash(alignment_File,ref,fq1,fq2,out_dir,sample,samtools,bwa)
check_file(alignment_File)


##generate bash file for these two file
bash_File = script_dir + f"/{sample}.run.sh"

current_dir = os.getcwd()
if sys.argv[2] == "DontMk":
    ##generate indel alignment bash file
    indel_File = script_dir + f"/{sample}.indel.alignment.sh"
    # print(indel_File)
    generate_indel_alignment_bash(indel_File,ref,sample,java,picard,out_dir,samtools,gatk)
    with open(bash_File,'w') as output:
        output.write(f'''bash {alignment_File};
        bash {indel_File};''')

elif sys.argv[2] == "Mkdup":
    indel_File = script_dir + f"/{sample}.indel.alignment.sh"
    # print(indel_File)
    generate_indel_alignment_bash(indel_File,ref,sample,java,picard,out_dir,samtools,gatk)
    generate_mkdup_bash(indel_File,sample,java,picard,out_dir)

    with open(bash_File,'w') as output:
        output.write(f'''bash {alignment_File};
bash {indel_File};''')

elif sys.argv[2] == "OnlyMkdup":
    indel_File = script_dir + f"/{sample}.indel.alignment.sh"
    os.system(f'rm {indel_File}')
    generate_mkdup_bash(indel_File,sample,java,picard,out_dir)
    with open(bash_File,'w') as output:
        output.write(f'''bash {indel_File};''')






