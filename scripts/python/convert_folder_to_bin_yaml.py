from os import listdir
from os.path import isfile, join
import sys
import yaml
path=snakemake.input[0]
#Get all files from a list
#We assume that the files have the following format: <name>.fa
#We will use <name> as an identifier in the created yaml

files = [f for f in listdir(path) if isfile(join(path, f))]
#print("Input files " + str(files))
bin_dir = {}
for f in files:
	bin_dir[f] =  f

bin_dir = {"samples": bin_dir}

with open("bin_list.yaml", "w") as f:
	yaml.dump(bin_dir, f)
   
