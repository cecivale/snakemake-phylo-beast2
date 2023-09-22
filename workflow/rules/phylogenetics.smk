# phylogenetics

rule iqtree:
    input:
        alignment = "results/data/{dataset}/aligned{sufix,.*}.fasta"
    output:
        tree =  "results/analysis/iqtree/{analysis}/{dataset}{sufix,.*}.treefile"
    params:
      outgroup = lambda wildcards: _get_analysis_param(wildcards, "iqtree", "outgroup"),
      outgroup_file = lambda wildcards: _get_analysis_param(wildcards, "iqtree", "outgroup_file"),
      alignment_outgroup =  "results/analysis/iqtree/{analysis}/temp/{dataset}{sufix,.*}.fasta",
      tree_args = lambda wildcards: _get_analysis_param(wildcards, "iqtree", "tree_args"),
      file_name = "results/analysis/iqtree/{analysis}/other/{dataset}{sufix,.*}",
      seed = _get_seed
    conda:
        "../envs/iqtree2.yaml"
    shell:
        """
        mkdir -p "results/analysis/iqtree/{wildcards.analysis}/temp/"
        mkdir -p "results/analysis/iqtree/{wildcards.analysis}/other/"
        cat {input.alignment} {params.outgroup_file} > {params.alignment_outgroup} 2>/dev/null  

        iqtree2 -s {params.alignment_outgroup}  \
        {params.tree_args} \
        -seed {params.seed} \
        -o '{params.outgroup}' \
        -pre {params.file_name} 
        mv {params.file_name}.treefile {output.tree}
        rm {params.alignment_outgroup}
        mv {params.file_name}.log logs/iqtree_{wildcards.dataset}{wildcards.sufix}.log
        """

# rule plot_tree_ann: #TODO add general metadata input
#     input:
#         tree = rules.iqtree.output.tree,
#         metadata =  rules.combine_metadata.output.combined
#     output:
#         fig = "results/analysis/iqtree/{analysis}/{dataset}{sufix,.*}.png",
#     params:
#         outgroup = lambda wildcards: _get_analysis_param(wildcards, "iqtree", "outgroup"),
#     script:
#         "../scripts/plot_mltree.R"

rule plot_iqtree: 
    input:
        tree = rules.iqtree.output.tree,
        ids = rules.combine_samples.output.combined
    output:
        fig = "results/analysis/iqtree/{analysis}/{dataset}{sufix,.*}.png",
    params:
        outgroup = lambda wildcards: _get_analysis_param(wildcards, "iqtree", "outgroup"),
    script:
        "../scripts/plot_tree.R"

# rule create_dates_file:
#     input:
#         metadata = rules.combine_metadata.output.combined
#     output:
#         dates = "results/{dataset}/data/dates_{subsampling}.{dseed}.csv",
#     run:
#         md = pd.read_csv(input.metadata, sep = "\t")
#         dates = pd.DataFrame() 
#         dates["strain"] = md["seq_name"]
#         dates["date"] = md["seq_name"].str.split("|", expand = True)[[2]]
#         dates.to_csv(output.dates, index = False, header = True, sep = ",")

rule treetime:
    input:
        metadata = rules.combine_samples.output.combined,
        alignment = "results/data/{dataset}/aligned{sufix,.*}.fasta",
        # tree = rules.iqtree.output.tree,
        tree = lambda wildcards: "results/analysis/iqtree/" + str(_get_analysis_param(wildcards, "treetime", "analysis_tree")) + "/{dataset}{sufix,.*}.treefile",
    output:
        tree = "results/analysis/treetime/{analysis}/{dataset}{sufix,.*}.nexus"
    params:
        clock_rate = lambda wildcards: _get_analysis_param(wildcards, "treetime", "clock_rate"),
        outgroup_file = lambda wildcards: _get_analysis_param(wildcards, "treetime", "outgroup_file"),
        alignment_outgroup =  "results/analysis/treetime/{analysis}/temp/{dataset}{sufix,.*}.fasta",
        outgroup = lambda wildcards: _get_analysis_param(wildcards, "treetime", "outgroup"),
    conda:
        "../envs/treetime.yaml"
    shell:
        """
        mkdir -p "results/analysis/treetime/{wildcards.analysis}/temp/"
        mkdir -p "results/analysis/treetime/{wildcards.analysis}/other/"
        cat {input.alignment} {params.outgroup_file} > {params.alignment_outgroup} 2>/dev/null  

        treetime --tree {input.tree} \
         --dates {input.metadata} \
         --aln {params.alignment_outgroup} \
         --name-column sample_id \
         --clock-rate {params.clock_rate} \
         --reroot '{params.outgroup}' \
         --max-iter 5 \
         --stochastic-resolve \
         --outdir results/analysis/treetime/{wildcards.analysis}/other/{wildcards.dataset}{wildcards.sufix}/
        
        mv results/analysis/treetime/{wildcards.analysis}/other/{wildcards.dataset}{wildcards.sufix}/timetree.nexus {output.tree}
        rm {params.alignment_outgroup}
        """

         
rule plot_tree: 
    input:
        tree = rules.treetime.output.tree,
        ids =  rules.combine_samples.output.combined
    output:
        fig = "results/analysis/treetime/{analysis}/{dataset}{sufix,.*}.png",
    params:
        outgroup = lambda wildcards: _get_analysis_param(wildcards, "treetime", "outgroup"),
    script:
        "../scripts/plot_tree.R"

