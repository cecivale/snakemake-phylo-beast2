# ------------------------------------------------------------------------------
#          ---        
#        / o o \    Snakemake workflow phylo-BEAST
#        V\ Y /V    Script to query metadata inside snakemake workflow
#    (\   / - \     
#     )) /    |     
#     ((/__) ||     Code by Ceci VA 
# ------------------------------------------------------------------------------

import ast, os
import numpy as np
import pandas as pd

from Bio import SeqIO


def attr_dict_to_query(attr_dict):
    s = str(attr_dict)
    query = s.replace(
        " ", "").replace(
        "{", "").replace(
        "}", "").replace(
        "'", "").replace(
        ",", "' & ").replace(
        "From:", " >= '").replace(
        "To:", " <= '").replace(
        "Not:", " != '").replace(
        ":", " == '").replace(
        "\\n", ",") + "'"
    return query

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

    if snakemake.input["metadata"] == snakemake.input["sequences"]:
        seq_format = get_seq_format(snakemake.input["sequences"])
        seq_ids = []
        for seq in SeqIO.parse(snakemake.input["sequences"], seq_format):
            seq_ids.append(seq.id)
        metadata = pd.DataFrame({"sample_id" : seq_ids, "seq_id" : seq_ids  })
    else:
        metadata = pd.read_csv(snakemake.input["metadata"], sep = "\t", dtype = "str")
    
    if len(metadata) == 1: 
        ids = metadata # No metadata case
    else:
        if snakemake.params["select"] is not None:
            query = attr_dict_to_query(snakemake.params["select"])
            print(query)
            metadata = metadata.query(query)
        
        # assign sample id
        ids = pd.DataFrame()
        if "sample_id" in metadata:
            ids["sample_id"] = metadata["sample_id"]
        elif "seq_id" in metadata and "date" in metadata and snakemake.params["deme"] is not None:
            ids["sample_id"] = metadata["seq_id"] + "|" + snakemake.params["deme"] + "|" + metadata["date"]
        elif "seq_id" in metadata and "date" in metadata:
            ids["sample_id"] = metadata["seq_id"] + "|" + metadata["date"]
        elif "seq_id" in metadata:
            ids["sample_id"] = metadata["seq_id"]
        else:
            print("Error: No sequence or sample id in metadata.")
       
        if "date" in metadata:
            ids["date"] = metadata["date"]
        # ids[snakemake.params["seq_id"]] = metadata[[snakemake.params["seq_id"]]]
        if "seq_id" in metadata:
            ids["seq_id"] = metadata["seq_id"]
        if snakemake.params["deme"] is not None: 
            ids[["deme"]] = snakemake.params["deme"] 


    metadata.to_csv(snakemake.output["metadata"], index = False, sep="\t")
    ids.to_csv(snakemake.output["ids"], index = False, sep="\t")

