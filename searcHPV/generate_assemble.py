import sys
import os
import re
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

script_dir = f'{script_dir}/assemble'
if not os.path.exists(script_dir):
    os.system(f"mkdir {script_dir}")
current_dir = config_dic[input_para]['tool_path'] + '/assemble'
with open(f"{script_dir}/assemble.sh",'w') as script:
    script.write(f'''python {current_dir}/extractReadName.py {config_file}
python {current_dir}/extractReadSequenceAWK.py {config_file}
bash  {script_dir}/extractReadSequence.sh
python {current_dir}/preprocessForPear.py {config_file}
python {current_dir}/pairEndMergingWithPear.py {config_file}
bash {script_dir}/pear.sh
python {current_dir}/preprocessForCap3.py {config_file}
python {current_dir}/cap3Assemble.py {config_file}
bash {script_dir}/cap3.sh
''')