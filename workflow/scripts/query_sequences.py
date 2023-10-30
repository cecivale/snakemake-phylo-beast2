# ------------------------------------------------------------------------------
#          ---        
#        / o o \    Snakemake workflow phylo-BEAST
#        V\ Y /V    Script to select sequences inside snakemake workflow
#    (\   / - \     
#     )) /    |     
#     ((/__) ||     Code by Ceci VA 
# ------------------------------------------------------------------------------

import os
import numpy as np
import pandas as pd

from Bio import SeqIO

def get_seq_format(seq_file):
    seq_name, seq_extension = os.path.splitext(seq_file)
    if seq_extension == ".fasta":
        return "fasta" 
    elif seq_extension in [".nex", ".nexus"]:
        return "nexus"
    else:
        print("Error: Only fasta and nexus files are recognized as valid sequence files.") 
        return

if __name__ == '__main__':
    ids = pd.read_csv(snakemake.input["ids"], sep = "\t", dtype = "str")

    with open(snakemake.output["sequences"], "w") as output:
        for seq_file in snakemake.input["sequences"]:
            seq_format = get_seq_format(seq_file)
            sequences = SeqIO.parse(seq_file, seq_format)
            for seq in sequences:
                id = seq.id
                if ids['seq_id'].str.contains(id).any():
                    seq.id = ids.query("seq_id=='" + id + "'")["sample_id"].item()
                    seq.description = ""
                    SeqIO.write(seq, output, 'fasta')



