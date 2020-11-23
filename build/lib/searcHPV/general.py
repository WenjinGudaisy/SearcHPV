########
#Function: some general functions
import os
#check if the directory exist and if not creat the directory
def mkdir(dir):
    if not os.path.exists(dir):
        os.system(f"mkdir {dir}")
    return None
#check if the file exist
def check_file(file):
    if not os.path.isfile(file):
        raise ValueError(f"{file} does not exist")
