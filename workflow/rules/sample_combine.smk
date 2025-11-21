rule combine_samples:
    message:
        """
        Combine sequences ids {input.ids} from demes.
        """
    input:
        ids = _get_ids_to_combine
    output:
        combined = "results/data/{dataset}/ids_combined{sufix,.*}.tsv" 
    log:
        # ("logs/combine_subsamples_{dataset}_{subsampling}.{dseed}.txt" if _is_structured
        # else "logs/combine_subsamples_{dataset}.txt")
    shell:
        """
        rm -f {output.combined}
        awk '(NR == 1) || (FNR > 1)' {input.ids} > {output.combined} 2>&1 | tee {log}
        """


rule combine_metadata:
    message:
        """
        Combine all metadata files
        """
    input:
        files = _get_all_metadata_files
    output:
        combined = "results/data/metadata_combined.tsv" 
    log:
        # ("logs/combine_subsamples_{dataset}_{subsampling}.{dseed}.txt" if _is_structured
        # else "logs/combine_subsamples_{dataset}.txt")
    shell:
        """
        rm -f {output.combined}
        awk '(NR == 1) || (FNR > 1)' {input.files} > {output.combined} 2>&1 | tee {log}
        """
