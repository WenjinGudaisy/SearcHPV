import subprocess
import os
from searcHPV.generate_call_fusion_virus import *
from searcHPV.selection_contig_call_virus_insertion import *

##############
#Fuction: run pipeline for hpv fusion points calling
#humRef: human reference genome
#virRef: virus reference genome
#out_dir: output directory for seacHPV
#n: poly(n), n*d(A/T/C/G), will report low confidence if contig contains poly(n)
def virus_fusion(humRef,virRef,out_dir,n):
    script_map = mapToRef(out_dir)
    script_mapHg = mapToHgRef(out_dir,humRef)
    script_mapVir = mapToVirRef(out_dir,virRef)
    os.system(f'chmod +x {script_map}')
    subprocess.call(script_map)
    os.system(f'chmod +x {script_mapHg}')
    subprocess.call(script_mapHg)
    os.system(f'chmod +x {script_mapVir}')
    subprocess.call(script_mapVir)

    count_sr(out_dir)
    contigDic = read_mapping_info(out_dir)
    combinedContigDic = cal_virus_ins(contigDic,out_dir)
    selectedAllContig = select_contig(combinedContigDic)
    siteConfidence = siteConf(out_dir)
    # print(siteConfidence.keys())
    # print(selectedAllContig.keys())
    filteredSelectedContig = filter_res(selectedAllContig,siteConfidence)
    write_to_file(filteredSelectedContig,siteConfidence,out_dir,n)
