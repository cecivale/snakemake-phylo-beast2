# ------------------------------------------------------------------------------
#          ---        
#        / o o \    Snakemake workflow phylo-BEAST
#        V\ Y /V    Script to output metadata file inside snakemake workflow
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

    if snakemake.params["metadata"] is not None:
        metadata = pd.read_csv(snakemake.params["metadata"], sep = "\t", dtype = "str")
    else:
        seq_format = get_seq_format(snakemake.params["sequences"])
        seq_ids = []
        for seq in SeqIO.parse(snakemake.params["sequences"], seq_format):
            seq_ids.append(seq.id)
        metadata = pd.DataFrame({"sample_id" : seq_ids, "seq_id" : seq_ids  })

        metadata.to_csv(snakemake.output["metadata"], index = False, sep="\t")




