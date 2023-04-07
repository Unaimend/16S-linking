import os
print(snakemake.input.samples)
print(str(snakemake.input.read_path))

input = []
with open(str(snakemake.input.samples), "r") as f:
    input=f.read().splitlines()

print(input)

files = list(map(lambda x: str(snakemake.input.read_path) + "/"  + x, input))

print(files)

os.mkdir(files[0])



