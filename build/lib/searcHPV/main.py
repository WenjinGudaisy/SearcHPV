import argparse 
from searcHPV.alignment import *
from searcHPV.genome_fusion import *
from searcHPV.assemble import *
from searcHPV.hpv_fusion import *

def main(): 

#pass the arguments from commline line
    parser = argparse.ArgumentParser(prog ='searcHPV', 
                                     description ='An HPV intgegration detection tool for targted capture sequencing') 
  
    parser.add_argument('-fastq1', type = str, dest= "fq1",
                        help ='sequencing data: fastq/fq.gz file', required=True,
                        )
    parser.add_argument('-fastq2', type = str, dest = 'fq2',
                        help ='sequencing data: fastq/fq.gz file', required=True,
                        )
    parser.add_argument('-humRef', type = str, dest = 'humRef',
                        help ='human reference genome: fasta file', required=True,
                        )
    parser.add_argument('-virRef', type = str, dest = 'virRef',
                        help ='HPV reference genome: fasta file', required=True,
                        )
    parser.add_argument('-window', type = int, default=300, dest = 'window',
                        help ='the length of region searching for informative reads, default=300', 
                        )
    parser.add_argument('-gz', action ='store_const', const = True,
                        default = False, dest ='gz', 
                        help ="if fastq files are in gz format")
    parser.add_argument('-output', type = str, default="./",dest = 'outputDir',
                        help ='output directory, default "./"', 
                        )
                          
    parser.add_argument('-alignment', action ='store_const', const = True,
                        default = False, dest ='alignment', 
                        help ="run the alignment step, step1")

    parser.add_argument('-genomeFusion', action ='store_const', const = True, 
                        default = False, dest ='genomeFusion', 
                        help ="call the genome fusion points, step2")

    parser.add_argument('-assemble', action ='store_const', const = True, 
                        default = False, dest ='assemble', 
                        help ="local assemble for each integration event, step3")

    parser.add_argument('-hpvFusion', action ='store_const', const = True, 
                        default = False, dest ='hpvFusion', 
                        help ="call the HPV fusion points, step4")
    



    args = parser.parse_args() 
    
  
    if args.alignment:
        alignment(fq1 = args.fq1, fq2 = args.fq2, humRef = args.humRef, virRef = args.virRef, outputDir = args.outputDir, gz = args.gz)
    elif args.genomeFusion: 
        genomeFusion(args.window,args.outputDir,args.virRef)
    elif args.assemble: 
        #check result from genomeFusion
        assemble(args.fq1, args.fq2, args.outputDir,args.virRef,args.gz)
    elif args.hpvFusion:
        #check result from assemble
        hpv_fusion(args.humRef,args.virRef,args.outputDir)
    else:
        alignment(fq1 = args.fq1, fq2 = args.fq2, humRef = args.humRef, virRef = args.virRef, outputDir = args.outputDir, gz = args.gz)
        genomeFusion(args.window,args.outputDir,args.virRef)
        assemble(args.fq1, args.fq2, args.outputDir,args.virRef,args.gz)
        hpv_fusion(args.humRef,args.virRef,args.outputDir)
  
    