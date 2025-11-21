

# alignment
rule align:
    message:
        """
        Align with mafft
        """
    input:
        sequences = rules.select_sequences.output.sequences
    output:
        alignment = "results/data/{dataset}/aligned{sufix,.*}.fasta"
    params:
        to_align = lambda wildcards: _get_dataset_param("align", wildcards)
    log:
        "logs/align_{dataset}_{sufix,.*}.txt"
    # conda:
    #     "envs/mafft.yaml"
    shell:
        """
        if [ {params.to_align} == False ]; then
            scp {input.sequences} {output.alignment}
            echo 'Sequences already aligned, skipping alignment step.'
        fi
        #TODO align with mafft
        """
