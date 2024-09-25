# Common functions

def _is_structured(wildcards):
    return (config["datasets"][wildcards.dataset].get("structure") is not None)

def _is_subsampled(wildcards):
    return (config["datasets"][wildcards.dataset].get("replicates") is not None)

def _is_filtered(wildcards):
    return (config["datasets"][wildcards.dataset].get("filter") is not None)

def _is_filtered_deme(wildcards, deme):
    return (config["datasets"][wildcards.dataset][deme].get("filter") is not None)

def _get_ids_to_subsample(wildcards):
    if _is_filtered(wildcards):
        return "results/data/{dataset}/{prefix,.*}ids_filtered.tsv"
    else:
        return "results/data/{dataset}/{prefix,.*}ids.tsv"

def _get_ids_to_combine(wildcards):
    if _is_structured(wildcards):
        demes = config["datasets"][wildcards.dataset]["structure"].keys()
        files = []
        for deme in demes:
            if config["datasets"][wildcards.dataset]["structure"][deme].get("subsample") is not None:
                # files.append("results/data/{dataset}/" + deme + "/" + "ids_subsampled{sufix,*.}.tsv") 
                files.append("results/data/{dataset}/" + deme + "/" + "ids_subsampled{sufix}.tsv") 
            elif config["datasets"][wildcards.dataset]["structure"][deme].get("filter") is not None:
                files.append("results/data/{dataset}/" + deme + "/" + "ids_filtered.tsv") 
            else:
                files.append("results/data/{dataset}/" + deme + "/" + "ids.tsv") 
            # files.append("results/{{dataset}}/data/" + deme + "/" + "ids" + 
            #     ("_filtered" if (config["dataset"][wildcards.dataset]["structure"][deme].get("filter") is not None) else "") + ".tsv")
        return files
    return "results/data/{dataset}/ids.tsv"

def _get_alignment(wildcards):
    if _get_dataset_param("mask", wildcards):
        return "results/data/{dataset}/masked{sufix,.*}.fasta",
    else:
        return "results/data/{dataset}/aligned{sufix,.*}.fasta",

def  _get_sequence_ids(wildcards):
    data_dir = "results/data/" + wildcards.dataset + "/"
    if _is_subsampled(wildcards) and _is_structured(wildcards):
        return data_dir + "ids_combined" + wildcards.sufix + ".tsv"
    if _is_subsampled(wildcards):
        return data_dir + "ids_subsampled" + wildcards.sufix + ".tsv"
    if _is_structured(wildcards):
        return data_dir + "ids_combined.tsv"
    if _is_filtered(wildcards):
        return data_dir + "ids_filtered.tsv"
    else:
        return data_dir + "ids.tsv"


def  _get_metadata_file(wildcards):
    return _get_dataset_param("metadata", wildcards)

def  _get_all_metadata_files(wildcards):
    datasets = config["datasets"].keys()
    metadata_files = []
    for d in datasets:
        # TODO if structure
        demes = config["datasets"][d]["structure"].keys()
        for deme in demes:
            metadata = config["datasets"][d]["structure"][deme].get("metadata")
            if metadata not in metadata_files:
                metadata_files.append(metadata)
    return metadata_files

def  _get_sequences_file(wildcards):
    return _get_dataset_param("sequences", wildcards)

def  _get_all_sequences_file(wildcards):
    if _is_structured(wildcards):
        demes = config["datasets"][wildcards.dataset]["structure"].keys()
        sequences_files = []
        for deme in demes:
            sequences = config["datasets"][wildcards.dataset]["structure"][deme].get("sequences")
            if sequences not in sequences_files:
                sequences_files.append(sequences)
        return sequences_files
    else:
        return [_get_dataset_param("sequences", wildcards)]


def  _get_select_params(wildcards):
    return _get_dataset_param("select", wildcards)

def  _get_seq_id(wildcards):
    seq_id = _get_dataset_param("seq_id", wildcards) or config["lapis"]["seq_id"]
    return seq_id

def _get_deme(wildcards):
    if _is_structured(wildcards):
        return config["datasets"][wildcards.dataset]["structure"][wildcards.prefix[0:-1]].get("deme", 
            wildcards.prefix[0:-1]) 
    else:
        return None

def _get_seed(wildcards):
    # Returns data seed if several replicates or 1
    return wildcards.sufix[1:] or 1

def _get_analysis_param(wildcards, param):
    return config["analyses"][wildcards.analysis].get(param, config["beast"].get(param))

def _get_run_param(wildcards, param):
    return config["run"][wildcards.analysis].get(param, config["run"].get(param))

def  _get_dataset_param(param, wildcards):
    if _is_structured(wildcards) and wildcards.get("prefix") is not None:
        return config["datasets"][wildcards.dataset]["structure"][wildcards.prefix[0:-1]].get(param, config["datasets"][wildcards.dataset].get(param))
    else:
        return config["datasets"][wildcards.dataset].get(param)

def  _get_subsampling_param(param, wildcards):
    if _is_structured(wildcards):
        return config["datasets"][wildcards.dataset]["structure"][wildcards.prefix[0:-1]]["subsample"].get(param)
    else:
        return config["datasets"][wildcards.dataset]["subsample"].get(param)

def _get_mrs(wildcards):
    ids = pd.read_csv(_get_sequence_ids(wildcards), sep = '\t')
    dates = ids['sample_id'].str.split("|", expand=True)[2]
    return max(dates)

