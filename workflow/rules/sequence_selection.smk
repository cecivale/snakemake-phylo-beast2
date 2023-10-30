rule select_sequences:
    message:
        """
        Load the selected sequences for each dataset.
        """
    input:
        ids = _get_sequence_ids,
        sequences = _get_all_sequences_file
    output:
        sequences = "results/data/{dataset}/sequences{sufix,.*}.fasta"
    log:
        # "logs/load_seqs_{dataset}_{subsampling}.{dseed}.txt" # TODO add log, change run to shell
    # conda:
    #     "envs/python-genetic-data.yaml"
    script: 
        "../scripts/query_sequences.py"
