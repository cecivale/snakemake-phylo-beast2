# phylodynamics

checkpoint beast:
    input:
        # alignment = "results/data/{dataset}/masked{sufix,.*}.fasta",
        alignment = lambda wildcards: _get_alignment(wildcards), 
        xml = lambda wildcards: _get_analysis_param(wildcards, "xml"),
        ids = lambda wildcards: _get_sequence_ids(wildcards)
    output:
        trace = "results/analysis/{analysis}/chains/{dataset}{sufix,.*}.{chain}.r{i}.log",
        trees = "results/analysis/{analysis}/chains/{dataset}{sufix,.*}.{chain}.r{i}.trees",
        is_converged = "results/analysis/{analysis}/chains/is_converged_{dataset}{sufix,.*}.{chain}.r{i}.txt"
    params:
        beast_command = lambda wildcards: _get_analysis_param(wildcards, "command"), 
        action = lambda wildcards: _get_run_param(wildcards, "action"),
        user_check = lambda wildcards: _get_run_param(wildcards, "user_check"),
        xml_params = lambda wildcards: str(_get_analysis_param(wildcards, "xml_params")).replace(":", "=").replace(
            "{", "\"").replace("}", "\"").replace(" ", "").replace("'", ""),
        folder_name = "results/analysis/{analysis}/chains",
        file_name = "{dataset}{sufix,.*}.{chain}.r{i}",
        state_file = "{dataset}{sufix,.*}.{chain}.state",
        previous_file_name = lambda wildcards: wildcards.dataset + wildcards.sufix + "." + wildcards.chain + ".r" + str(int(wildcards.i) - 1) if wildcards.i != 0 else wildcards.dataset + wildcards.sufix + wildcards.chain + ".r{i}"
    log:
        "logs/beast_{analysis}_{dataset}_{sufix,.*}_{chain}.r{i}.txt"
    benchmark:
        "benchmarks/beast_{analysis}_{dataset}_{sufix,.*}_{chain}.r{i}.benchmark.txt"
    threads:
        lambda wildcards: _get_analysis_param(wildcards, "threads"),
    resources:
        runtime = lambda wildcards: _get_analysis_param(wildcards, "time"),
        mem_per_cpu = lambda wildcards: _get_analysis_param(wildcards,  "mem_mb")
    shell:
        """
        mkdir -p {params.folder_name}/running
        
        if [ {wildcards.i} == 0 ]; then

            rm -f results/analysis/{wildcards.analysis}/chains/is_converged_{wildcards.dataset}{wildcards.sufix}.{wildcards.chain}.r*
            
            if [ {params.action} == "resume" ]; then
                ACTION="resume"
            else
                ACTION="overwrite"
            fi
            
        else
            ACTION="resume"
            scp {params.folder_name}/{params.previous_file_name}.log {params.folder_name}/running/{params.file_name}.log 
            scp {params.folder_name}/{params.previous_file_name}.trees {params.folder_name}/running/{params.file_name}.trees

            rm {params.folder_name}/{params.previous_file_name}.log {params.folder_name}/{params.previous_file_name}.trees {params.folder_name}/is_converged_{params.previous_file_name}.txt
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
            -statefile "{params.folder_name}/{params.state_file}" \
            -$ACTION {input.xml} 2>&1 | tee -a {log}
        touch {params.folder_name}/running/{params.file_name}.log
        touch {params.folder_name}/running/{params.file_name}.trees

        # We need the chains to be written in a different path than snakemake output so snakemake does not delete them if job fails
        # So we moved them once they are finished
        [ -f {params.folder_name}/running/{params.file_name}.log ] && mv {params.folder_name}/running/{params.file_name}.log {output.trace}
        [ -f {params.folder_name}/running/{params.file_name}.trees ] && mv {params.folder_name}/running/{params.file_name}.trees {output.trees} 
        """
        
def _is_converged(wildcards):
    # ?TODO implement automatic check of convergence based on ESS > 200
    path = expand("results/analysis/{analysis}/chains/{dataset}{sufix}.{chain}.r", 
        analysis = wildcards.analysis,
        dataset = wildcards.dataset,
        sufix = wildcards.sufix,
        chain = wildcards.chain)
    trace = glob.glob(path[0] + "*.log")
    if len(trace) == 0:
        runs = 0
    else:
        start = ".r"
        end = ".log"
        runs = int(trace[0][trace[0].find(start)+len(start):trace[0].find(end)])

    with checkpoints.beast.get(analysis = wildcards.analysis, dataset = wildcards.dataset, 
        sufix = wildcards.sufix, 
        chain = wildcards.chain, i = runs).output.is_converged.open() as f:
        s = f.read().strip()
        print(s)
        if s == "YES":
            return "results/analysis/{analysis}/chains/{dataset}{sufix,.*}.{chain}.r" + str(runs) + ".{output}"
        elif s == "NO":
            runs += 1
            checkpoints.beast.get(analysis = wildcards.analysis, dataset = wildcards.dataset, 
        sufix = wildcards.sufix, 
        chain = wildcards.chain, i = runs)


rule aggregate_runs:
    input:
        run = _is_converged 
    output:
        chain = "results/analysis/{analysis}/chains/{dataset}{sufix,.*}.{chain}.{output}",
    shell:
        """
        mv {input.run} {output.chain} 2>&1 | tee -a {log}
        """

def _get_chains(wildcards):
    files = expand(
        "results/analysis/{{analysis}}/chains/{{dataset}}{{sufix,.*}}.{chain}.{{output}}",
        chain = range(1, _get_analysis_param(wildcards, "chains") + 1))
    return files

rule combine_chains:
    message: 
        """
        Combine chain files: {input.chain_files} with LogCombiner.
        """
    input:
        chain_files = _get_chains  
    output:
        combined_chain = "results/analysis/{analysis}/{dataset}{sufix}.{output}",
    # log:
    #     "logs/combine_trace_{dataset}_{analysis}_{subsampling}.{dseed}.txt"
    # benchmark:
    #     "benchmarks/combine_trace_{dataset}_{analysis}_{subsampling}.{dseed}.benchmark.txt"
    params:
        burnin =  lambda wildcards: _get_analysis_param(wildcards, "burnin"),
        input_command = lambda wildcards, input: " -log ".join(input) 
    shell:
        """
        logcombiner -log {params.input_command} -o {output.combined_chain} -b {params.burnin}  2>&1 | tee -a {log}
        """
rule downsample_trees:
    input:
        trees =  "results/analysis/{analysis}/{dataset}{sufix}.trees"
    output:
        downsampled_trees =  "results/analysis/{analysis}/{dataset}{sufix}.ds.trees"
    log:
        "logs/downsample_trees_{analysis}_{dataset}{sufix}.txt"
    params:
        command = config["logcombiner"].get("command"),
        resample = config["logcombiner"].get("resample"),
        burnin = lambda wildcards: _get_analysis_param(wildcards, "burnin"), 
    shell:
       """
        {params.command} -log {input.trees} -o {output.downsampled_trees} -b {params.burnin} -resample {params.resample}  2>&1 | tee {log}
        """

rule summarize_trees:
    message: 
        """
        Summarize trees to {wildcards.topo} tree with median node heights with TreeAnnotator.
        """
    input:
        trees = rules.downsample_trees.output.downsampled_trees
    output:
        summary_tree =  "results/analysis/{analysis}/{dataset}{sufix}.{topo}.tree"
    log:
        "logs/summarise_typedtrees_{dataset}_{analysis}_{sufix}_{topo}.txt"
    benchmark:
        "benchmarks/summarise_typedtree_{dataset}_{analysis}_{sufix}_{topo}.benchmark.txt"
    params:
        command = config["treeannotator"].get("command"),
        burnin = 0,
        heights = config["treeannotator"].get("heights"),
        topology = lambda wildcards: wildcards.topo
    shell:
        """
        {params.command} -topology {params.topology} -heights {params.heights} -b {params.burnin} -lowMem {input.trees} {output.summary_tree} 2>&1 | tee -a {log} 
        """

