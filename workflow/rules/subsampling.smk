# ------------------------------------------------------------------------------
#          ---        
#        / o o \    Snakemake workflow phylo-BEAST
#        V\ Y /V    Rules for subsampling
#    (\   / - \     
#     )) /    |     
#     ((/__) ||     Code by Ceci VA 
# ------------------------------------------------------------------------------


rule subsample:
    message:
        """
        Subsample sequences 
        """
    input:
        ids = _get_ids_to_subsample,
        metadata = rules.load_metadata.output.metadata 
    output:
        ids = "results/data/{dataset}/{prefix,.*}ids_subsampled.{dseed}.tsv"
    params:
        include = lambda wildcards: _get_subsampling_param("include", wildcards), 
        exclude = lambda wildcards: _get_subsampling_param("exclude", wildcards),
        n = lambda wildcards: _get_subsampling_param("n", wildcards),
        method = lambda wildcards: _get_subsampling_param("method", wildcards),
        weights = lambda wildcards: _get_subsampling_param("weights", wildcards)
    log:
        "logs/subsample_{dataset}_{prefix,.*}.{dseed}.txt"
    # conda:
    #     "envs/r-genetic-data.yaml"
    script:
        "../scripts/subsample.R"

