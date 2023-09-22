# phylodynamics

checkpoint beast:
#TODO change beast in _get_analysis_param to wildcard method
    input:
        alignment = "results/data/{dataset}/aligned{sufix,.*}.fasta",
        xml = lambda wildcards: _get_analysis_param(wildcards, "beast", "xml"),
        ids = lambda wildcards: _get_sequence_ids(wildcards)
    output:
        trace = "results/analysis/beast/{analysis}/chains/{dataset}{sufix,.*}.{chain}.r{i}.log",
        trees = "results/analysis/beast/{analysis}/chains/{dataset}{sufix,.*}.{chain}.r{i}.trees",
        is_converged = "results/analysis/beast/{analysis}/chains/is_converged_{dataset}{sufix,.*}.{chain}.r{i}.txt"
    params:
        beast_command = lambda wildcards: _get_analysis_param(wildcards, "beast", "command"), 
        action = lambda wildcards: _get_analysis_param(wildcards, "beast", "action"),
        user_check = lambda wildcards: _get_analysis_param(wildcards, "beast", "user_check"),,
        xml_params = lambda wildcards: str(_get_analysis_param(wildcards, "beast", "xml_params")).replace(":", "=").replace(
            "{", "\"").replace("}", "\"").replace(" ", "").replace("'", ""),
        # mrs = lambda wildcards: _get_mrs(wildcards),
        folder_name = "results/analysis/beast/{analysis}/chains",
        file_name = "{dataset}{sufix,.*}.{chain}.r{i}",
    log:
        "logs/beast_{analysis}_{dataset}_{sufix,.*}_{chain}.r{i}.txt"
    benchmark:
        "benchmarks/beast_{analysis}_{dataset}_{sufix,.*}_{chain}.r{i}.benchmark.txt"
    threads:
        lambda wildcards: _get_analysis_param(wildcards, "beast", "threads"),
    resources:
        runtime = lambda wildcards: _get_analysis_param(wildcards, "beast", "time"),
        mem_mb = lambda wildcards: _get_analysis_param(wildcards, "beast", "mem_mb")
    shell:
        """
        mkdir -p {params.folder_name}/running
        
        if [ {params.action} == "resume" ] || [ -f {output.is_converged} -a  $(< {output.is_converged})=="NO" ]; then
            ACTION="resume"
        else
            ACTION="overwrite"
        fi
        
        if [ {params.user_check} == True ]; then
            echo NO > {output.is_converged}
        else
            echo YES > {output.is_converged}
        fi
    
        {params.beast_command} \
            -D aligned_fasta={input.alignment} \
            -D {params.xml_params} \
            -D file_name={params.folder_name}/running/{params.file_name} \
            -seed {wildcards.chain} \
            -statefile "{params.folder_name}/{params.file_name}.state" \
            -$ACTION {input.xml} 2>&1 | tee -a {log}
        
        # We need the chains to be written in a different path than snakemake output so snakemake does not delete them if job fails
        # So we moved them once they are finished
        [ -f {params.folder_name}/running/{params.file_name}.log ] && mv {params.folder_name}/running/{params.file_name}.log {output.trace}
        [ -f {params.folder_name}/running/{params.file_name}.trees ] && mv {params.folder_name}/running/{params.file_name}.trees {output.trees} 

        """



runs = 0
def _is_converged(wildcards):
    # TODO implement automatic check of convergence based on ESS > 200
    global runs
    with checkpoints.beast.get(analysis = wildcards.analysis, dataset = wildcards.dataset, 
        sufix = wildcards.sufix, 
        chain = wildcards.chain, i = runs).output.is_converged.open() as f:
        s = f.read().strip()
        if s == "YES":
            return expand("results/analysis/beast/{{analysis}}/chains/{{dataset}}{{sufix,.*}}.{{chain}}.r{i}.{{output}}", i = range(0, runs+1))
        elif s == "NO":
            runs += 1
            checkpoints.beast.get(analysis = wildcards.analysis, dataset = wildcards.dataset, 
        sufix = wildcards.sufix, 
        chain = wildcards.chain, i = runs)


rule aggregate_runs:
    input:
        runs = _is_converged 
    output:
        combined_run = "results/analysis/beast/{analysis}/chains/{dataset}{sufix,.*}.{chain}.{output}",
    params:
        input_command = lambda wildcards, input: " -log ".join(input.runs)
    shell:
        """
        logcombiner -log {params.input_command} -o {output.combined_run} 2>&1 | tee -a {log}
        """
def _get_chains(wildcards):
    files = expand(
        "results/analysis/beast/{{analysis}}/chains/{{dataset}}{{sufix,.*}}.{chain}.{{output}}",
        chain = range(1, _get_analysis_param(wildcards, "beast", "chains") + 1))
    return files

rule combine_chains:
    message: 
        """
        Combine chain files: {input.chain_files} with LogCombiner.
        """
    input:
        chain_files = _get_chains  
    output:
        combined_chain = "results/analysis/beast/{analysis}/{dataset}{sufix,.*}.{output}",
    # log:
    #     "logs/combine_trace_{dataset}_{analysis}_{subsampling}.{dseed}.txt"
    # benchmark:
    #     "benchmarks/combine_trace_{dataset}_{analysis}_{subsampling}.{dseed}.benchmark.txt"
    params:
        burnin =  lambda wildcards: _get_analysis_param(wildcards, "beast", "burnin"),
        input_command = lambda wildcards, input: " -log ".join(input) 
    shell:
        """
        logcombiner -log {params.input_command} -o {output.combined_chain} -b {params.burnin}  2>&1 | tee -a {log}
        """

# rule mcc_tree:
#     message: 
#         """
#         Summarize trees to Maximum clade credibility tree with median node heights with TreeAnnotator.
#         """
#     input:
#         combined_trees = rules.combine_typedtrees.output.combined_typedtrees
#     output:
#         summary_tree = "results/{dataset}/beast/{analysis}/{subsampling}.{dseed}.comb.mcc.typed.node.tree"
#     log:
#         "logs/summarise_typedtrees_{dataset}_{analysis}_{subsampling}.{dseed}.txt"
#     benchmark:
#         "benchmarks/summarise_typedtree_{dataset}_{analysis}_{subsampling}.{dseed}.benchmark.txt"
#     params:
#         burnin = 0,
#         heights = "median"
#     shell:
#         """
#         treeannotator -heights {params.heights} -b {params.burnin} -lowMem {input.combined_trees} {output.summary_tree} 2>&1 | tee -a {log} 
#         """
