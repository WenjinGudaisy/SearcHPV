import os
import sys
import re
from searcHPV.general import *

#############
#Function: generate bash script for mapping contigs to human+hpv genome
#out_dir:out_dir for searcHPV
#humRef: human reference genome
#return: path of bash script 
def mapToRef(out_dir):
    newRef = f'{out_dir}/hg_hpv.fa'
    sitesPath = f'{out_dir}/assemble/'
    listSites = os.listdir(sitesPath)
    outputPath = f'{out_dir}/call_fusion_virus/'
    mkdir(outputPath)
    #if there are sites
    if listSites != ['']:
        with open( f'{outputPath}/alignContigsToRef.sh','w') as bashFile:
            bashFile.write('#!/bin/bash\n')
            for site in listSites:
                if ".sh" not in site:
                    contigPath = f'{sitesPath}/{site}/pearOutput/'            
                    sitePath = f'{out_dir}/call_fusion_virus/{site}/'
                    mkdir(sitePath)
                    bashFile.write( f'''bwa mem -R \'@RG\\tID:hpv\\tSM:hpv\\tLB:hpv\\tPL:ILLUMINA\' -M -t 8 \
{newRef} {contigPath}/{site}.all.fa.cap.contigs > {sitePath}/{site}.contig.sam;
samtools view -bhS {sitePath}/{site}.contig.sam > {sitePath}/{site}.contig.bam;
samtools sort {sitePath}/{site}.contig.bam {sitePath}/{site}.contig.sort;
samtools index {sitePath}/{site}.contig.sort.bam;
rm {sitePath}/{site}.contig.bam;\n''')
    return f'{outputPath}/alignContigsToRef.sh'




#############
#Function: generate bash script for mapping contigs to human genome
#out_dir:out_dir for searcHPV
#humRef: human reference genome
#return: path of bash script 
def mapToHgRef(out_dir,humRef):
    sitesPath = f'{out_dir}/assemble/'
    listSites = os.listdir(sitesPath)
    outputPath = f'{out_dir}/call_fusion_virus/'
    mkdir(outputPath)
    #if there are sites
    if listSites != ['']:
        with open( f'{outputPath}/alignContigsToGenome.sh','w') as bashFile:
            bashFile.write('#!/bin/bash\n')
            for site in listSites:
                if ".sh" not in site:
                    contigPath = f'{sitesPath}/{site}/pearOutput/'            
                    sitePath = f'{out_dir}/call_fusion_virus/{site}/'
                    mkdir(sitePath)
                    bashFile.write( f'''bwa mem -R \'@RG\\tID:hpv\\tSM:hpv\\tLB:hpv\\tPL:ILLUMINA\' -M -t 8 \
{humRef} {contigPath}/{site}.all.fa.cap.contigs > {sitePath}/{site}.contigToGenome.sam;
samtools view -bhS {sitePath}/{site}.contigToGenome.sam > {sitePath}/{site}.contigToGenome.bam;
samtools sort {sitePath}/{site}.contigToGenome.bam {sitePath}/{site}.contigToGenome.sort;
samtools index {sitePath}/{site}.contigToGenome.sort.bam;
rm {sitePath}/{site}.contigToGenome.bam;\n''')
    return f'{outputPath}/alignContigsToGenome.sh'

#############
#Function: generate bash script for mapping contigs to hpv genome
#out_dir:out_dir for searcHPV
#virRef: virus reference genome
#return: path of bash script 
def mapToVirRef(out_dir,virRef):
    sitesPath = f'{out_dir}/assemble/'
    listSites = os.listdir(sitesPath)
    outputPath = f'{out_dir}/call_fusion_virus/'
    mkdir(outputPath)
    #if there are sites
    if listSites != ['']:
        with open( f'{outputPath}/alignContigsToVirus.sh','w') as bashFile:
            bashFile.write('#!/bin/bash\n')
            for site in listSites:
                if ".sh" not in site:
                    contigPath = f'{sitesPath}/{site}/pearOutput/'
                    sitePath = f'{out_dir}/call_fusion_virus/{site}/'
                    bashFile.write( f'''bwa mem -R \'@RG\\tID:hpv\\tSM:hpv\\tLB:hpv\\tPL:ILLUMINA\' -M -t 8 \
{virRef} {contigPath}/{site}.all.fa.cap.contigs > {contigPath}/{site}.contigToVirus.sam;
samtools view -bhS {contigPath}/{site}.contigToVirus.sam >  {contigPath}/{site}.contigToVirus.bam;
samtools sort {contigPath}/{site}.contigToVirus.bam {sitePath}/{site}.contigToVirus.sort;
samtools index {sitePath}/{site}.contigToVirus.sort.bam;
samtools faidx {contigPath}/{site}.all.fa.cap.contigs;
rm {contigPath}/{site}.contigToVirus.bam;\n''')
    return f'{outputPath}/alignContigsToVirus.sh'

#############
#Function: count number of supprotive reads for each contig
#out_dir:out_dir for searcHPV
def count_sr(out_dir):
    dirName = f'{out_dir}/assemble/'
    sites = os.listdir(dirName)

    key = ""
    rowList = []
    for site in sites:
        if ".sh" not in site:
            sampleDic = {}
            cap3File = f'{dirName}/{site}/cap3Output/{site}.CAP3'
            with open(cap3File) as eachFile:
                for i in range(0, 5):
                    eachFile.readline()

                for row in eachFile:
                    if re.match('^\*',row):
                        
                        key = site+"."+\
                            row.rstrip().strip(' ').split(" ")[1]+row.rstrip().strip(' ').split(" ")[2]#contig number
                        rowList = []
                    else:
                        rowList.append(row)
                    sampleDic[key] = rowList
                    if re.match("^DETAILED",row):
                        break

        #         for items in sampleDic.items():
        #            print(items)

            with open(f'{out_dir}/call_fusion_virus/{site}/srNum.txt', "w") as sampleOutput:
                for key, value in sampleDic.items():
                    sampleOutput.write(key+"\t"+str(len(value)-1)+'\n')
    return None














