import argparse 
from alignment import *
from genome_fusion import *
from assemble import *
from hpv_fusion import *

def main(): 
  
    parser = argparse.ArgumentParser(prog ='searcHPV', 
                                     description ='An HPV intgegration detection tool for targted capture sequencing') 
  
    parser.add_argument('-fastq1', type = str, nargs ='+', 
                        help ='sequencing data: fastq/fq.gz file', 
                        )
    parser.add_argument('-fastq2', type = str, nargs ='+', 
                        help ='sequencing data: fastq/fq.gz file', 
                        )
    parser.add_argument('-humRef', type = str, nargs ='+', 
                        help ='human reference genome: fasta file', 
                        )
    parser.add_argument('-virRef', type = str, nargs ='+', 
                        help ='HPV reference genome: fasta file', 
                        )
    parser.add_argument('-window', type = int, nargs ='+', default=300,
                        help ='the length of region searching for informative reads, default=300', 
                        )
    parser.add_argument('-output', type = str, nargs ='+', default="./",
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
        alinment()
    if args.genomeFusion: 
        #check result from alignment
        genomeFusion()
    if args.assemble: 
        #check result from genomeFusion
        assemble()
    if args.hpvFusion()
        #check result from assemble
        hpvFusion()
        
  
    