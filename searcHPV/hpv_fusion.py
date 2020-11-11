from generate_call_fusion_virus import *
from selection_contig_call_virus_insertion import *

generate_scripts()

#run scripts
subprocess.call("hpv_fusion.sh")

select_contigs()

