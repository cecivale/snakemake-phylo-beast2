
# ----------------------------------------------------------------------------------
#          ---        
#        / o o \    Snakemake workflow phylo-BEAST
#        V\ Y /V    Rules to load sequence data and other metadata for the analysis
#    (\   / - \     
#     )) /    |     
#     ((/__) ||     Code by Ceci VA 
# -----------------------------------------------------------------------------------


rule select_samples:
    message:
        """
        Load metadata file and select samples based on config.
        """
    input: 
        metadata = _get_metadata_file,
        sequences = _get_sequences_file
    output:
        metadata = "results/data/{dataset}/{prefix,.*}metadata.tsv",
        ids = "results/data/{dataset}/{prefix,.*}ids.tsv",
    params:
        select = _get_select_params,
        # seq_id = _get_seq_id,
        deme = _get_deme
    # log:
    #     "logs/select_samples_{dataset}" + "{prefix,.*}"[0:-1] +".txt"
    # conda:
    #     "envs/python-genetic-data.yaml"
    script: 
        "../scripts/query_metadata.py"
