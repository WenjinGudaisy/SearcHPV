from generate_alignment import *
import subprocess
from general import *


def alignment(fq1, fq2, humRef, virRef, outputDir, rmdup):
    #make output dir
    mkdir(outputDir)
    scriptDir = f'{outputDir}/alignment'
    mkdir(scriptDir)

    #catenate humRef and virRef
    ref = catRef(humRef,virRef, outputDir)

    #generate alignment bash

    alignmentFile = scriptDir + "/orignal.alignment.sh"
    indelFile = scriptDir + "/indel.alignment.sh"
    generate_alignment_bash(alignmentFile,ref,fq1,fq2,outputDir)
    generate_indel_alignment_bash(indelFile,ref,outputDir)
    check_file(alignmentFile)
    check_file(indelFile)

    #rmdup for alignment results
    
    bashFile = scriptDir + f"/alignment.sh"
    if not rmdup:
    ##generate indel alignment bash file
        with open(bashFile,'w') as output:
            output.write(f'''bash {alignmentFile};
bash {indelFile};''')
    else:
        generate_mkdup_bash(indelFile,outputDir)


        with open(bashFile,'w') as output:
            output.write(f'''bash {alignmentFile};
bash {indelFile};''')

    #run scripts
    check_file(bashFile)
    subprocess.call("alignment.sh")
