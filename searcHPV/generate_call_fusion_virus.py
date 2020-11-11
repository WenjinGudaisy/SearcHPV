import os
import sys
from general import *

if len(sys.argv) != 2:
    raise ValueError("Please input a configuration file")

config_file = sys.argv[1]
config_dic = read_config(config_file)
output_para = "Output"
input_para = "Input File"
out_dir = config_dic[output_para]['output_directory']
script_dir = config_dic[output_para]['script_directory']
sample = config_dic[output_para]["output_sample_name"]

bam = f'{out_dir}/alignment/{sample}.RG.indelre.mkdup.sort.bam'

out_dir = f'{out_dir}/call_fusion_virus'
if not os.path.exists(out_dir):
    os.system(f"mkdir {out_dir}")

script_dir = f'{script_dir}/call_fusion_virus'
if not os.path.exists(script_dir):
    os.system(f"mkdir {script_dir}")

current_dir = config_dic[input_para]['tool_path']
with open(f"{script_dir}/call_fusion_virus.sh",'w') as script:
    script.write(f'''#python {current_dir}/call_fusion_virus/mappingToGenome.py {config_file}
python {current_dir}/call_fusion_virus/mappingToHPV.py {config_file}
bash {script_dir}/alignContigsToGenome.sh
echo \'mapping to genome done\'
bash {script_dir}/alignContigsToHPV.sh
echo \'mapping to HPV done\'
python {current_dir}/call_fusion_virus/count_sr.py {config_file}
python {current_dir}/call_fusion_virus/selection_contig_call_virus_insertion.py {config_file}
echo \'call fusion for virus done\'''')
