rule select_sequences:
    message:
        """
        Load the selected sequences for each dataset.
        """
    input:
        ids = _get_sequence_ids,
        sequences = _get_sequences_file
    output:
        sequences = "results/data/{dataset}/sequences{sufix,.*}.fasta"
    log:
        # "logs/load_seqs_{dataset}_{subsampling}.{dseed}.txt" # TODO add log, change run to shell
    # conda:
    #     "envs/python-genetic-data.yaml"
    shell:
        """
        python3 workflow/scripts/query_sequences.py \
            --ids_file {input.ids} \
            --output_file {output.sequences} \
            --drop_incomplete True 2>&1 | tee {log}
        """
