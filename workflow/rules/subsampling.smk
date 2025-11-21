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
        metadata = "results/data/{dataset}/{prefix,.*}metadata.tsv",
        groups =  lambda wildcards: expand("results/data/{{dataset}}/{{prefix,.*}}groups.tsv", proxy=[]) if _get_subsampling_param("method", wildcards) == "groups" else [],
        weights =  lambda wildcards: expand("results/data/{{dataset}}/weights.tsv", proxy=[]) if _get_subsampling_param("method", wildcards) == "weights" else []
    output:
        ids = "results/data/{dataset}/{prefix,.*}ids_subsampled.{dseed}.tsv"
    params:
        include = lambda wildcards: _get_subsampling_param("include", wildcards), 
        exclude = lambda wildcards: _get_subsampling_param("exclude", wildcards),
        n = lambda wildcards: _get_subsampling_param("n", wildcards),
        method = lambda wildcards: _get_subsampling_param("method", wildcards),
        weights = lambda wildcards: _get_subsampling_param("weights", wildcards),
        grouping_var = lambda wildcards: _get_subsampling_param("grouping_var", wildcards),
    log:
        "logs/subsample_{dataset}_{prefix,.*}.{dseed}.txt"
    # conda:
    #     "envs/r-genetic-data.yaml"
    script:
        "../scripts/subsample.R"

