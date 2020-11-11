from generate_alignment import generate_scripts
import subprocess

#generate_scripts
generate_scripts()


#run scripts
subprocess.call("alignment.sh")
