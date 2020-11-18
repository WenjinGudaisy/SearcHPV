########
#Function: some general functions
import os
#check if the directory exist and if not creat the directory
def mkdir(dir):
    if not os.path.exists(dir):
        os.system(f"mkdir {dir}")
