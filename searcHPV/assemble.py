from generate_assemble import generate_scripts
import subprocess

#generate_scripts
generate_scripts()


#run scripts
subprocess.call("assemble.sh")
