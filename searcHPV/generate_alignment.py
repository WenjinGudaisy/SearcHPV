import sys
import os
import re
#config Dictionary formate:key1 = "Input File",key2 = "reference_file = ",value = full_path

#####################
##catenate reference genome
#Function: cat human genome + virus genome
#Parameters: humRef: human genome
#virRef: virus genome
#outputDir: output directory
def catRef(humRef, virRef, outputDir):
    os.system(f'cat {humRef} {virRef} > {outputDir}/hg_hpv.fa')
    newRef = f'{outputDir}/hg_hpv.fa'
    return(newRef)

#####################
##generate bash script for index reference genome
#Function: index customized reference
#bash_file: output bash script
#humRef: human genome
#virRef: virus genome
#newRef: customize genome
#outputDir: output directory
def indexRef(bash_file,ref):
    with open(bash_file,'w') as output:
        output.write(f'''bwa index {ref}
samtools faidx {ref}
picard CreateSequenceDictionary R={ref} O={ref.replace('.fa','.dict')}
''')
    return None




#####################
##generate alignment bash file
#Function: generate alignment bash file for a sample
#Parameters: bash_file: output bash file
#ref: reference gene
#fq1: fastq1 file
#fq2: fastq2 file
#out_dir: outputPath
#gz: if fastq file is in gz format: default = True
def generate_alignment_bash(bash_file,ref,fq1,fq2,out_dir,memory,thread,gz = True):
    with open(bash_file,'w') as output:
        if gz:
            output.write(f'''
    bwa mem -t {thread} {ref} '<zcat {fq1}' '<zcat {fq2}' > {out_dir}/alignment.sam
    samtools view -@ {thread} -bhS {out_dir}/alignment.sam > {out_dir}/alignment.bam
    samtools sort -@ {thread} {out_dir}/alignment.bam -o {out_dir}/alignment.sort.bam
    rm {out_dir}/alignment.sam
    echo \'alignment done\'
    ''')
        else:
            output.write(f'''
    bwa mem -t {thread} {ref} '<cat {fq1}' '<cat {fq2}' > {out_dir}/alignment.sam
    samtools view -@ {thread} -bhS {out_dir}/alignment.sam > {out_dir}/alignment.bam
    samtools sort -@ {thread} {out_dir}/alignment.bam -o {out_dir}/alignment.sort.bam
    rm {out_dir}/alignment.sam
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
def generate_indel_alignment_bash(bash_File,ref,out_dir,memory,thread):
    with open(bash_File,'w') as output:
        output.write(f'''picard \
AddOrReplaceReadGroups \
I={out_dir}/alignment.bam \
O={out_dir}/alignment.RG.bam \
RGID=sample \
RGPL=illumina \
RGPU=NA \
RGSM=sample \
RGLB=sample
samtools sort -@ {thread} {out_dir}/alignment.RG.bam -o {out_dir}/alignment.RG.sort.bam
samtools index -@ {thread} {out_dir}/alignment.RG.sort.bam
GenomeAnalysisTK \
-Xmx{memory} \
-T RealignerTargetCreator \
-R {ref} \
-I {out_dir}/alignment.RG.sort.bam \
-o {out_dir}/alignment.RG.intervals
GenomeAnalysisTK \
-Xmx{memory} \
-T IndelRealigner \
-R {ref} \
-I {out_dir}/alignment.RG.sort.bam \
-targetIntervals {out_dir}/alignment.RG.intervals \
-o {out_dir}/alignment.RG.indelre.bam
samtools index -@ {thread} {out_dir}/alignment.RG.indelre.bam
echo \'indel alignment done\'''')



def generate_mkdup_bash(bash_File,out_dir,thread):
    with open(bash_File,'w') as output:
        output.write(f'''
samtools sort -@ {thread} -n {out_dir}/alignment.RG.indelre.bam -o {out_dir}/alignment.RG.indelre.sortbyQ.bam
picard MarkDuplicates \
I={out_dir}/alignment.RG.indelre.sortbyQ.bam \
O={out_dir}/alignment.RG.indelre.mkdup.bam \
M={out_dir}/alignment.RG.indelre.mkdup.txt \
TAGGING_POLICY=All ASSUME_SORT_ORDER=queryname
samtools sort -@ {thread} {out_dir}/alignment.RG.indelre.mkdup.bam -o {out_dir}/alignment.RG.indelre.mkdup.sort.bam
samtools index -@ {thread} {out_dir}/alignment.RG.indelre.mkdup.sort.bam
echo \'indel alignment done\'''')


#delete intermediate bam file
def rm_inter_bam(bash_File,out_dir):
    
    with open(bash_File,'w') as output:
        output.write(f'''rm {out_dir}/alignment.RG.indelre.mkdup.bam
rm {out_dir}/alignment.RG.indelre.bam*
rm {out_dir}/alignment.RG.indelre.bai
rm {out_dir}/alignment.RG.sort.bam*
rm {out_dir}/alignment.RG.indelre.sortbyQ.bam 
rm {out_dir}/alignment.RG.bam
rm {out_dir}/alignment.bam
rm {out_dir}/alignment.sort.bam
''')



############################







