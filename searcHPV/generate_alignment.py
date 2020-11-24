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
def indexRef(bash_file,humRef,virRef,newRef,outputDir):
    with open(bash_file,'w') as output:
        output.write(f'''#bwa index {humRef}
#bwa index {virRef}
bwa index {newRef}
#samtools faidx {humRef}
#samtools faidx {virRef}
samtools faidx {newRef}
#java -Xmx4g -jar /home/wenjingu/tools/picard.jar \
CreateSequenceDictionary R={humRef} O={humRef.replace('.fa','.dict')}
#java -Xmx4g -jar /home/wenjingu/tools/picard.jar \
CreateSequenceDictionary R={virRef} O={virRef.replace('.fa','.dict')}
java -Xmx4g -jar /home/wenjingu/tools/picard.jar \
CreateSequenceDictionary R={newRef} O={newRef.replace('.fa','.dict')}
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
def generate_alignment_bash(bash_file,ref,fq1,fq2,out_dir):
    with open(bash_file,'w') as output:
        output.write(f'''
bwa mem {ref} '<zcat {fq1}' '<zcat {fq2}' > {out_dir}/alignment.sam
samtools view -bhS {out_dir}/alignment.sam > {out_dir}/alignment.bam
samtools sort -o {out_dir}/alignment.sort {out_dir}/alignment.bam
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
def generate_indel_alignment_bash(bash_File,ref,out_dir):
    with open(bash_File,'w') as output:
        output.write(f'''java -Xmx4g -jar /home/wenjingu/tools/picard.jar \
AddOrReplaceReadGroups \
I={out_dir}/alignment.bam \
O={out_dir}/alignment.RG.bam \
RGID=sample \
RGPL=illumina \
RGPU=NA \
RGSM=sample \
RGLB=sample
samtools sort -o {out_dir}/alignment.RG.sort.bam {out_dir}/alignment.RG.bam
samtools index {out_dir}/alignment.RG.sort.bam
java -Xmx4g -jar /sw/med/centos7/gatk/3.7/GenomeAnalysisTK.jar \
-T RealignerTargetCreator \
-R {ref} \
-I {out_dir}/alignment.RG.sort.bam \
-o {out_dir}/alignment.RG.intervals
java -Xmx4g -jar /sw/med/centos7/gatk/3.7/GenomeAnalysisTK.jar \
-T IndelRealigner \
-R {ref} \
-I {out_dir}/alignment.RG.sort.bam \
-targetIntervals {out_dir}/alignment.RG.intervals \
-o {out_dir}/alignment.RG.indelre.bam
samtools index {out_dir}/alignment.RG.indelre.bam
echo \'indel alignment done\'''')



def generate_mkdup_bash(bash_File,out_dir):
    with open(bash_File,'a') as output:
        output.write(f'''
samtools sort -n -o {out_dir}/alignment.RG.indelre.sortbyQ.bam {out_dir}/alignment.RG.indelre.bam
java -Xmx4g -jar /home/wenjingu/tools/picard.jar MarkDuplicates \
I={out_dir}/alignment.RG.indelre.sortbyQ.bam \
O={out_dir}/alignment.RG.indelre.mkdup.bam \
M={out_dir}/alignment.RG.indelre.mkdup.txt \
TAGGING_POLICY=All ASSUME_SORT_ORDER=queryname
samtools sort -o {out_dir}/alignment.RG.indelre.mkdup.sort.bam {out_dir}/alignment.RG.indelre.mkdup.bam
samtools index {out_dir}/alignment.RG.indelre.mkdup.sort.bam
echo \'indel alignment done\'''')
 
 

############################







