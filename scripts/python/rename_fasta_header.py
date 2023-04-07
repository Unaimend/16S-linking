from os import listdir, mkdir
from os.path import isfile, join
import sys
import yaml
print(sys.version)
print(sys.executable)
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq


bins=snakemake.input

for bin in bins:
        print("Renaming " + bin)
        #read all sequences of a bin
        fasta_sequences = SeqIO.parse(open(bin),'fasta')
        # create a new empty sequence that has the current file name as id
        # this results in downstream program treating all contigs as a single one, i.e., 
        # a genome.

        # this might be BS a we should just prepend the file name
        newSeq = SeqIO.SeqRecord(Seq(""), description="")
        newSeq.id = Path(bin).name
        for fasta in fasta_sequences:
                s = str(newSeq.seq)
                s2 = str(fasta.seq)
                newSeq.seq = Seq(s + s2)
                count = SeqIO.write([newSeq], "renamed_bins/" + Path(bin).name, "fasta")
                


