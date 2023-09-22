# ------------------------------------------------------------------------------
#          ---        
#        / o o \    Snakemake workflow phylo-BEAST
#        V\ Y /V    Script to query metadata inside snakemake workflow
#    (\   / - \     
#     )) /    |     
#     ((/__) ||     Code by Ceci VA 
# ------------------------------------------------------------------------------

import ast
import numpy as np
import pandas as pd


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


if __name__ == '__main__':

    metadata = pd.read_csv(snakemake.input["metadata"], sep = "\t", dtype = "str")

    if snakemake.params["select"] is not None:
        query = attr_dict_to_query(snakemake.params["select"])
        print(query)
        metadata = metadata.query(query)
    
    ids = pd.DataFrame()
    if "sample_id" in metadata:
        ids["sample_id"] = metadata["sample_id"]
    elif snakemake.params["deme"] is not None:
        ids["sample_id"] = metadata["seq_id"] + "|" + snakemake.params["deme"] + "|" + metadata["date"]
    else:
        ids["sample_id"] = metadata["seq_id"] + "|" + metadata["date"]

    ids["date"] = metadata["date"]
    ids[snakemake.params["seq_id"]] = metadata[[snakemake.params["seq_id"]]]

    if snakemake.params["deme"] is not None: 
        ids[["deme"]] = snakemake.params["deme"] 

    




    ids.to_csv(snakemake.output["ids"], index = False, sep="\t")

