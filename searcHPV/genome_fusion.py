import pysam
import sys
from general import *

config_file = sys.argv[1]
config_dic = read_config(config_file)
input_para = "Input File"
virus_chrm = config_dic[input_para]['virus_chromsome_name']


bam = sys.argv[2]

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
import re
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

print ('chrm\tpos\tsingle_evidence\tpaired_evidence')
for i in b:
    print (f'{i[0]}\t{i[1]}\t{i[2]}\t{i[3]}')

