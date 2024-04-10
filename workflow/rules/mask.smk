

rule mask:
    message:
        """
        Mask bases in alignment {input.alignment}
          - masking {params.mask_from_beginning} from beginning
          - masking {params.mask_from_end} from end
          - masking other sites: {params.mask_sites}
        Rule from nextstrain ncov workflow.
        """
    input:
        alignment = "results/data/{dataset}/aligned{sufix,.*}.fasta"
    output:
        masked =  "results/data/{dataset}/masked{sufix,.*}.fasta"
    log:
        "logs/mask_{dataset}{sufix,.*}.txt"
    benchmark:
        "benchmarks/mask_{dataset}{sufix,.*}}.txt"
    params:
        mask_from_beginning = config["mask"].get("mask_from_beginning"),
        mask_from_end = config["mask"].get("mask_from_end"),
        mask_sites = config["mask"].get("mask_sites"),
        mask_terminal_gaps = config["mask"].get("mask_terminal_gaps")
    # conda:
    #     "envs/python-genetic-data.yaml"
    script: 
        "../scripts/mask_alignment.py"
