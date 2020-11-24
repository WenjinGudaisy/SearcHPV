import pysam
import sys
from searcHPV.general import *
import re

##################
#find the end of each cigar string slot
#pos: start position of each cigar string slot
#cigar: cigar string returned by pysam

def find_end(pos,cigar):
    letters = re.split('\d',cigar)
    letters = list(filter(None, letters))
    numbers = re.split('\D',cigar)
    numbers = list(filter(None, numbers))
    numbers = list(map(int,numbers))
    length = 0
    if letters[0] == 'S' or letters[0] == 'H':
        return int(pos)
    elif letters[0] == 'M':
        for i,letter in enumerate(letters):
            if letter == 'M':
                length+=numbers[i]
            elif letter == 'D':
                length+=numbers[i]
            elif letter == 'I':
                length-=numbers[i]
            else:
                continue
        return int(pos)+length

#####################
#identify fusion points on human genome
#Function: identify fusion points on human genome by PE and SP reads
#bam: alignment.RG.indelre.mkdup.sort.bam
#virus_chrm: name of virus chromosome in virus reference genome
##out_dir:output path
#return: output file
def define_fusion(bam,virus_chrm,out_dir):
    out_dir = os.path.abspath(out_dir)
    paired_read_li = []
    clipped_read_li = []

    samfile = pysam.AlignmentFile(bam, 'rb')
    if samfile.count(virus_chrm) != 0:
        for read in samfile.fetch(virus_chrm):
            # some filter here
            if read.is_duplicate is True:
                continue
            if read.is_qcfail is True:
                continue
            if read.is_unmapped is True:
                continue
            if read.is_secondary is True:
                continue
            mapq = read.mapping_quality
            if mapq < 50:
                continue
            if read.next_reference_name != virus_chrm:
                paired_read_li.append(read)
            else:
                cigar = read.cigarstring
                if "H" in cigar:
                    clipped_read_li.append(read)
                    
                elif "S" in cigar:
                    clipped_read_li.append(read)
                    
                else:
                    continue



    next_pos_li = []
    for read in paired_read_li:
        next_pos_li.append(f"{read.next_reference_name}:{read.next_reference_start}")
    next_pos_li_sort = sorted(next_pos_li)

    from operator import itemgetter
    clipped_infor_li = []
    for read in clipped_read_li:
        if read.is_duplicate is True:
            continue
        if read.is_qcfail is True:
            continue
        if read.is_unmapped is True:
            continue
        if read.is_secondary is True:
            continue
        if read.has_tag('SA'):
            name = read.query_name
            SATag = read.get_tag('SA')
            nchrm = SATag.split(",")[0]
            npos = SATag.split(",")[1]
            ndir = SATag.split(",")[2]
            ncigar = SATag.split(",")[3]
            mapq = int(SATag.split(",")[4])
            if mapq >= 50:
                clipped_infor_li.append([nchrm,npos,ndir,ncigar,name])
                # if nchrm == '10':
                #     print(read,samfile.mate(read))

    #10:7041919:11:1;17:41153083:47:99;17:41153186:12:99;2:102871348:6:1;2:17687921:8:1;3:130887923:2:8;5:97050065:7:1;6:151124993:11:1;6:2255928:8:1;6:71772584:8:1;X:113760995:6:1
        else:
            continue
    sorted_clipped_infor_li = sorted(clipped_infor_li,key=itemgetter(0,1))
    #print(sorted_clipped_infor_li)

    pos_li = []
    for read in sorted_clipped_infor_li:
        chrm = read[0]
        if chrm == virus_chrm:
            continue
        else:
            end_pos = find_end(read[1],read[3])
            pos_li.append((chrm,end_pos))
            

    a = set(pos_li)
    #print(pos_li)
    b = []
    for i in a:
        count = pos_li.count(i)
        i += (str(count),)
        b.append(i)
    b = sorted(b)

    for i,pos_infor in enumerate(b):
        chrm = pos_infor[0]
        pos = int(pos_infor[1])
        samfile = pysam.AlignmentFile(bam, 'rb')
        paired_evidence_count = 0
        if pos < 150:
            start_pos = 1
        else:
            start_pos = pos -150
        if samfile.count(chrm,start_pos,pos+150) != 0:
            for read in samfile.fetch(chrm,start_pos,pos+150):
                # some filter here
                if read.is_duplicate is True:
                    continue
                if read.is_qcfail is True:
                    continue
                if read.is_unmapped is True:
                    continue
                if read.is_secondary is True:
                    continue
                mapq = read.mapping_quality
                if mapq < 50:
                    continue
                if read.next_reference_name == virus_chrm:
                    paired_evidence_count += 1
        b[i] += (str(paired_evidence_count),)
    
    with open(f'{out_dir}/genome_fusion.txt','w') as output:
        output.write('chrm\tpos\tsingle_evidence\tpaired_evidence\n')
        for i in b:
            output.write(f'{i[0]}\t{i[1]}\t{i[2]}\t{i[3]}\n')
    return f'{out_dir}/genome_fusion.txt'

##################
#cluster fusion points within certain base pair
#resultin: genome_fusion.txt
#outf: output file
#window: length of window for cluster, default = 100
def cluster_result(resultin,outf,window=100):
    with open(resultin) as resultin:
        with open(outf,'w') as output:
            out_li = []
            for line in resultin.read().rstrip().split('\n'):
                elements = line.rstrip()
                if elements == "":
                    out_li.append(elements)
                else:
                    pos_li = elements.split(';')
                    pos_li = sorted(pos_li)
                    new_pos_li = [pos_li[0]]
                    for i,pos in enumerate(pos_li[1:]):
                        old_infor = new_pos_li[-1].split(':')
                        old_chrm = old_infor[0]
                        old_p = int(old_infor[1])
                        old_count = int(old_infor[2])
                        infor = pos.split(':')
                        chrm = infor[0]
                        p = int(infor[1])
                        single_count = int(infor[2])
                        if chrm == old_chrm:
                            if abs(p-old_p) > window:
                                new_pos_li.append(':'.join(infor))
                            else:
                                if single_count > old_count:
                                    new_pos_li[-1] = ':'.join(infor)
                        else:
                            new_pos_li.append(':'.join(list(map(str,infor))))
                    
                    pos_infor = ';'.join(new_pos_li)   
                    output.write(pos_infor)
    return None

    



