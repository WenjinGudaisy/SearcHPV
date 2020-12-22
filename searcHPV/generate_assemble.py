import sys
import os
import re
import pysam
from searcHPV.general import *

########################
#get informative reads (SP + PE) from bam file
#chrm:chromosome on human genome
#pos: start position
#bam: original bam file
#windowSize: default = 300, length of window for extracting reads
#ref_genome:name of virus chromosome
def get_reads(chrm,pos,bam,ref_genome,windowSize=300):
    samfile = pysam.AlignmentFile(bam, "rb")
    read_name_list = []
    #print (f'{chrm}:{pos}')
    for read in samfile.fetch(chrm,pos-int(windowSize),pos+int(windowSize)):
        read_name = read.query_name    
        read_ref = read.reference_name
        pair_ref = read.next_reference_name

        if read.is_duplicate is True:
            continue
        if read.is_qcfail is True:
            continue
        if read.is_unmapped is True:
            continue
        if read.is_secondary is True:
            continue
#pair end reads
        if pair_ref == ref_genome:
            read_name_list.append(read_name)
            
            
#split reads
        else:
            if read.has_tag('SA'):
                name = read.query_name
                SATag = read.get_tag('SA')
                nchrm = SATag.split(",")[0]
                if nchrm == ref_genome:
                    npos = SATag.split(",")[1]
                    ndir = SATag.split(",")[2]
                    ncigar = SATag.split(",")[3]
                    read_name_list.append(read_name)
                    
    read_name_set = set(read_name_list)

    samfile.close()
    return read_name_set

########################
#Function:
#extract informative reads name from original bam file
#bam: original bam file, alignment.RG.indelre.mkdup.sort.bam
#fusionRes: genome fusion result file from genome fusion step, all.filtered.clustered.result
#out_dir:out put directory for assemble
def extract_read_name(bam,fusionRes,out_dir,virRef):
    #####for this step, need to modify from Yifan's script
    #candidate_in = os.popen(f'grep \'{sample}\' {out_dir}/call_fusion/all.filtered.clustered.result')
    #############
    
    with open(fusionRes) as candidate_in:
        candidate_in = candidate_in.read().rstrip()

    #if there are fusion points
    if candidate_in != "":
        candidate_in = candidate_in.split(";")
#for each candidate site, collect reads with truncation and abberant pairs upstream and downstream windowbp and reads in HPV
    # sampleDir = f'/home/wenjingu/myScratch/ourResult/newVersion/targetedExomeResult/{sample}.{windowSize}'
    # try:
    #     subprocess.call(f'mkdir {sampleDir}')
    # except:
    #     print('dir exists')
    # try:
    #     subprocess.call(f'mkdir {sampleDir}/readName')
    # except:
    #     print('dir exists')
    #find the virus_chrm
    with open(virRef) as virRefFile:
        virus_chrm = virRefFile.readline().replace('>','').replace('\n','')
    for site in candidate_in:
        chrm = site.split(':')[0]
        pos = int(site.split(':')[1])
        #print(chrm,pos,bam,virus_chrm)
        read_list = get_reads(chrm=chrm,pos=pos,bam=bam,ref_genome =virus_chrm)
        outf_path = f'{out_dir}/{chrm}.{pos}/'
        mkdir(outf_path)

        outf = open(f'{outf_path}/readName.txt', 'w')
        
    #my dir
        for each in read_list:
            print(f'>{each}',file=outf)
        outf.close()
    return None

########################
#Function:generate bash script for extracting read sequence
#out_dir:output directory for assemble
#fq1:raw sequenction data, fastq file
#fq2:raw sequencing data, fastq file
#return: path of bash script
def extract_read_seq(out_dir,fq1,fq2):
    with open(f'{out_dir}/extractReadSequence.sh','w') as bash:
        bash.write('#!/bin/bash\n')
        for site in os.listdir(out_dir):
            if ".sh" not in site:
                readNamePath = f'{out_dir}/{site}/readName.txt'
                outputPath = f'{out_dir}/{site}/supportiveReads/'
                if not os.path.exists(outputPath):
                    os.system(f'mkdir -p {outputPath}')
                tab = '\"\\t\"'
                newLine = '\"\\n\"'
                bash.write(f'''zcat {fq1} | awk '{{if(NR%4!=0)ORS=" ";else ORS={newLine}}}1' | awk 'NR==FNR{{a[$1]=($1{tab}$2{newLine}$3{newLine}$4{newLine}$5)}}NR>FNR{{gsub(/>/,"@");if(a[$1]!={newLine})print a[$1]}}' - {readNamePath} | grep -v '^$' > {outputPath}/{site}.informativeReads.1.fq;
zcat {fq2} | awk '{{if(NR%4!=0)ORS=" ";else ORS={newLine}}}1' | awk 'NR==FNR{{a[$1]=($1{tab}$2{newLine}$3{newLine}$4{newLine}$5)}}NR>FNR{{gsub(/>/,"@");if(a[$1]!={newLine})print a[$1]}}' - {readNamePath} | grep -v '^$' > {outputPath}/{site}.informativeReads.2.fq;
''')
        bash.write('echo \'extract informative sequences done\'')
    return f'{out_dir}/extractReadSequence.sh'

########################
#Function:formalize fasta file for PEAR
#out_dir:output directory for assemble
def preprocessForPear(out_dir):
    listSites = os.listdir(out_dir)
    for site in listSites:
        if '.sh' not in site:
            outputPath = f'{out_dir}/{site}/supportiveReads/'
            with open(f'{outputPath}/{site}.informativeReads.1.fq','r') as inputFile:
                with open(f'{outputPath}/{site}.formalizedReads.1.fq','w') as outputFile:
                    i = 0
                    readList = inputFile.read().split('\n')
                    for each in readList:
                        if '@' in each and '+' in readList[i+4]:
                            rowSeq = readList[i+1] + readList[i+2]+readList[i+3]
                            rowInfo = readList[i+5] + readList[i+6]+readList[i+7]
                            outputFile.write(each+'\n'+rowSeq+'\n+\n'+rowInfo+'\n')
                        elif '@' in each:
                            rowSeq = readList[i+1] + readList[i+2]+readList[i+3]
                            rowInfo = readList[i+1] + readList[i+2]+readList[i+3]
                            outputFile.write(each+'\n'+rowSeq+'\n+\n'+rowInfo+'\n')
                        i += 1

            with open(f'{outputPath}/{site}.informativeReads.2.fq','r') as inputFile:
                with open(f'{outputPath}/{site}.formalizedReads.2.fq','w') as outputFile:
                    i = 0
                    readList = inputFile.read().split('\n')
                    for each in readList:
                        if '@' in each and '+' in readList[i+4]:
                            rowSeq = readList[i+1] + readList[i+2]+readList[i+3]
                            rowInfo = readList[i+5] + readList[i+6]+readList[i+7]
                            outputFile.write(each+'\n'+rowSeq+'\n+\n'+rowInfo+'\n')
                        elif '@' in each:
                            rowSeq = readList[i+1] + readList[i+2]+readList[i+3]
                            rowInfo = readList[i+1] + readList[i+2]+readList[i+3]
                            outputFile.write(each+'\n'+rowSeq+'\n+\n'+rowInfo+'\n')
                        i += 1  
    return None

########################
#Function:generate bash script for PEAR
#out_dir:output directory for assemble
#return:path of bash script
def pear(out_dir):
    listSites = os.listdir(out_dir)
    with open(f'{out_dir}/pear.sh','w') as output:
        output.write('#!/bin/bash\n')
        for site in listSites:
            if ".sh" not in site:
                outputPath = f'{out_dir}/{site}/pearOutput/'
                fqPath = f'{out_dir}/{site}/supportiveReads/'
                if not os.path.exists(outputPath):
                    os.mkdir(outputPath)
                output.write(f'''
    pear \
    -f {fqPath}/{site}.formalizedReads.1.fq \
    -r {fqPath}/{site}.formalizedReads.2.fq \
    -o {outputPath}/{site}''')
    return f'{out_dir}/pear.sh'

#############
#Function: formalize results from PEAR for CAP3
#out_dir:output directory for assemble
def preprocessForCap3(out_dir):
    listSites = os.listdir(out_dir)
    for site in listSites:
        if ".sh" not in site:
            pearOutputPath = f'{out_dir}/{site}/pearOutput/'

            #change unassemble reads to be fasta
            with open(f'{pearOutputPath}/{site}.unassembled.forward.fastq',encoding='utf8',errors='ignore') as fq:
                with open(f'{pearOutputPath}/{site}.1.fa','w') as fa:
                    #fq = fq.read().rstrip().decode('uft-16').split('\n')
                    fq = fq.read().rstrip().split('\n')
                    i = 0
                    for eachRow in fq:
                        if '@' in eachRow:
                            newName = eachRow.replace('@','>')
                            rowSeq = fq[i+1]
                            fa.write(newName + '.1' + '\n' + rowSeq + '\n')
                        i += 1
            with open(f'{pearOutputPath}/{site}.unassembled.reverse.fastq',encoding='utf8',errors='ignore') as fq:
                with open(f'{pearOutputPath}/{site}.2.fa','w') as fa:
                    #fq = fq.read().rstrip().decode('uft-16').split('\n')
                    fq = fq.read().rstrip().split('\n')
                    i = 0
                    for eachRow in fq:
                        if '@' in eachRow:
                            newName = eachRow.replace('@','>')
                            rowSeq = fq[i+1]
                            fa.write(newName + '.2' + '\n' + rowSeq + '\n')
                        i += 1
            #merge unassembled reads
            os.system(f'cat {pearOutputPath}/{site}.1.fa {pearOutputPath}/{site}.2.fa > {pearOutputPath}/{site}.unassembled.fa')

            #change assembled reads to be fasta
            with open(f'{pearOutputPath}/{site}.assembled.fastq',encoding='utf8',errors='ignore') as fq:
                with open(f'{pearOutputPath}/{site}.assembled.fa','w') as fa:
                    #fqList = fq.read().rstrip().decode('uft-16').split('\n')
                    fqList = fq.read().rstrip().split('\n')
                    i = 0
                    for eachRow in fqList:
                        if '@' in eachRow:
                            newName = eachRow.replace('@','>')
                            rowSeq = fqList[i+1]
                            fa.write(newName + '\n' + rowSeq + '\n')
                        i += 1
            #merge all reads
            os.system(f'cat {pearOutputPath}/{site}.assembled.fa {pearOutputPath}/{site}.unassembled.fa > {pearOutputPath}/{site}.all.fa')

    return None


##############
#Function: generate bash script for CAP3
#out_dir:output directory for assemble
#return: path of script for CAP3
def cap3(out_dir):
    listSites = os.listdir(out_dir)
    with open(f'{out_dir}/cap3.sh','w') as bashFile:
        bashFile.write('#!/bin/bash\n')
        for site in listSites:
            if ".sh" not in site:
                faPath = f'{out_dir}/{site}/pearOutput/'
                cap3OutputPath = f'{out_dir}/{site}/cap3Output/'
                mkdir(cap3OutputPath)
                bashFile.write(f'''
        cap3 {faPath}/{site}.all.fa > {cap3OutputPath}/{site}.CAP3
        ''')
        bashFile.write('echo \'CAP3 done\'')
    return f'{out_dir}/cap3.sh'