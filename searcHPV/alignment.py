from searcHPV.generate_alignment import *
import subprocess
from searcHPV.general import *
#####################
#pipeline of alignment
#Function:
#1. customize reference genome
#2. generate alignment bash script
#3. run alignment script
#fq1: paired-end raw sequencing data, fastq1
#fq2: paired-end raw sequencing data, fastq2
#humRef: human reference genome
#outputDir: output directory
#multi: if fastq file is in gz format: default = True
#index: if True, index the reference files; if False, not index the references files
#memory: memory size allocated
def alignment(fq1, fq2, humRef, virRef, outputDir, index, gz, memory,thread):
    #make output dir
    outputDir = os.path.abspath(outputDir)
    mkdir(outputDir)
    scriptDir = f'{outputDir}/alignment'
    mkdir(scriptDir)
    #catenate humRef and virRef
    ref = catRef(humRef,virRef, outputDir)

    #generate alignment bash

    alignmentFile = scriptDir + "/orignal.alignment.sh"
    indelFile = scriptDir + "/indel.alignment.sh"
    generate_alignment_bash(alignmentFile,ref,fq1,fq2,scriptDir,memory,thread,gz)
    generate_indel_alignment_bash(indelFile,ref,scriptDir, memory,thread)
    check_file(alignmentFile)
    check_file(indelFile)
    bashFile = scriptDir + f"/alignment.sh"

    #remove intermediate file
    rmInter = scriptDir + f"/rm_inter.sh"
    rm_inter_bam(rmInter,scriptDir)

    if index:
        #index humRef and virRef
        indexFileHum = scriptDir + "/index_hum.sh"
        indexRef(indexFileHum,humRef)
        indexFileVir = scriptDir + "/index_vir.sh"
        indexRef(indexFileVir,virRef)
        indexFile = scriptDir + "/index.sh"
        indexRef(indexFile,ref)
    else:
        indexFile = scriptDir + "/index.sh"
        indexRef(indexFile,ref)

    if index:
    ##generate mkdup alignment bash file
        mkdupFile = scriptDir + "/mkdup.alignment.sh"
        generate_mkdup_bash(mkdupFile,scriptDir,thread)

        with open(bashFile,'w') as output:
            output.write(f'''#!/bin/bash
bash {indexFileHum};
bash {indexFileVir};
bash {indexFile};
bash {alignmentFile};
bash {indelFile};
bash {mkdupFile};
bash {rmInter};''')
    else:
        ##generate mkdup alignment bash file
        mkdupFile = scriptDir + "/mkdup.alignment.sh"
        generate_mkdup_bash(mkdupFile,scriptDir,thread)

        with open(bashFile,'w') as output:
            output.write(f'''#!/bin/bash
bash {indexFile};
bash {alignmentFile};
bash {indelFile};
bash {mkdupFile}
bash {rmInter};''')


    #run scripts
    check_file(bashFile)
    os.system(f'chmod +x {bashFile}')
    subprocess.call(bashFile)
