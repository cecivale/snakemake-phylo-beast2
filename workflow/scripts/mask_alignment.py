# ------------------------------------------------------------------------------
#          ---        
#        / o o \    Project:  cov-euphydyn
#        V\ Y /V    Mask alignment. Script adapted from nextstrain ncov workflow.
#    (\   / - \     01 February 2023
#     )) /    |     
#     ((/__) ||     Code by Ceci VA 
# ------------------------------------------------------------------------------

import urllib, json, requests, os, sys, ast
import numpy as np
import pandas as pd

from Bio import SeqIO
from Bio.Seq import Seq
from io import StringIO

# Script and rule adapted from Nextrain ncov workflow

def mask_terminal_gaps(seq):
    L = len(seq)
    seq_trimmed = seq.lstrip('-')
    left_gaps = L - len(seq_trimmed)
    seq_trimmed = seq_trimmed.rstrip('-')
    right_gaps = L - len(seq_trimmed) - left_gaps
    return "N"*left_gaps + seq_trimmed + "N"*right_gaps


if __name__ == '__main__':

    begin_length = 0
    if snakemake.params["mask_from_beginning"]:
        begin_length = snakemake.params["mask_from_beginning"]
    end_length = 0
    if snakemake.params["mask_from_end"]:
        end_length = snakemake.params["mask_from_end"]

    with open(snakemake.input["alignment"], 'r') as sequences, open(snakemake.output["masked"], 'w') as outfile:
        records = SeqIO.parse(sequences, 'fasta')
        for record in records:
            seq = str(record.seq)
            if snakemake.params["mask_terminal_gaps"]:
                seq = mask_terminal_gaps(seq)

            start = "N" * begin_length
            middle = seq[begin_length:-end_length]
            end = "N" * end_length
            seq_list = list(start + middle + end)
            if snakemake.params["mask_sites"]:
                for site in snakemake.params["mask_sites"]:
                    if seq_list[site-1]!='-':
                        seq_list[site-1] = "N"
            record.seq = Seq("".join(seq_list))
            SeqIO.write(record, outfile, 'fasta')


