import os
import sys
from searcHPV.general import *
from searcHPV.generate_call_fusion import *

def genomeFusion(window,out_dir,virRef):

    bam = f'{out_dir}/alignment/alignment.RG.indelre.mkdup.sort.bam'
    check_file(bam)
    #make dir for output
    out_dir = f'{out_dir}/call_fusion'
    mkdir(out_dir)

    #find the virus_chrm
    with open(virRef) as virRefFile:
        virus_chrm = virRefFile.readline().replace('>','').replace('\n','')

    #identify fusion points on genome
    res = define_fusion(bam,virus_chrm,out_dir) #generate out_dir/genome_fusion.txt

    #filter and cluster fusion points
    ##sort result
    os.system(f'(head -n 1 {res} && tail -n +2 {res} | sort -k3,3rn) > {out_dir}/genome_fusion.sort.txt')


    #change format for cluster
    chrm_li = range(1,22)
    chrm_li = list(map(str,chrm_li))
    chrm_li+=['X','Y']

    with open(f'{out_dir}/all.result','w') as output:   
            fusion_li = []
            with open(f'{out_dir}/genome_fusion.sort.txt') as res:
                res.readline()
                for line in res.read().rstrip().split('\n'):
                    elements = line.rstrip().split('\t')
                    if elements != ['']:
                        #split read count
                        chrom = elements[0]
                        pos = elements[1]
                        count = elements[2]
                        #pair-end read count
                        pair_count = elements[3]
                        if chrom in chrm_li:
                            fusion_li.append(f'{elements[0]}:{elements[1]}:{count}:{pair_count}')
            to_print = ';'.join(fusion_li)
            output.write(f'{to_print}')

    ##cluster the events within 100bp from each other, maybe becasue of SVs or CNVs
    cluster_result(f'{out_dir}/all.result',f'{out_dir}/all.clustered.result')


    ##filter for sites with at least 2 split read  and 2 pairs of read support(high cutoff) and their summation greater than 5
    with open(f'{out_dir}/all.clustered.result') as inf:
        with open(f'{out_dir}/all.filtered.clustered.result','w') as outf:
            for line in inf.read().rstrip().split('\n'):
                elements = line.rstrip()
                if elements == "":
                    outf.write(elements)
                else:
                    pos_li = elements.split(';')
                    new_pos_li = []
                    for pos in pos_li:
                        single_count = int(pos.split(':')[2])
                        pair_count = int(pos.split(':')[3])
                        if single_count > 2 and pair_count >2:
                        #or (single_count > 5) or (pair_count > 5)
                            new_pos_li.append(pos + ':high')
                        elif single_count + pair_count >= 5:
                            new_pos_li.append(pos + ':low')
                            
                    new_pos_infor = ';'.join(new_pos_li)
                    outf.write(new_pos_infor)
    return None