from pathlib import Path

out_dir_path = Path(config["out_dir"])
draft_dir_path = Path(config["draft_dir"])
fastq_dir_path = Path(config["fastq_dir"])

log_dir_path = out_dir_path / config["log_dir"]
error_dir_path  = out_dir_path / config["error_dir"]
benchmark_dir_path =  out_dir_path / config["benchmark_dir"]
cluster_log_dir_path = Path(config["cluster_log_dir"])

if "species_list" not in config:
    config["species_list"] = [d.name for d in draft_dir_path.iterdir() if d.is_dir()]

def species_output_prefix_from_wildcards(wildcards):
    return out_dir_path / "{species}/{species}".format(species=wildcards.species)

def species_dir_from_wildcards(wildcards):
    return out_dir_path / "{species}".format(species=wildcards.species)

def species_fastq_dir_from_wildcards(wildcards):
    return out_dir_path / "{species}/fastq".format(species=wildcards.species)


localrules: all

rule all:
    input:
        expand(draft_dir_path / "{species}/{species}.draft.fasta.bwt", species=config["species_list"]),
        expand(draft_dir_path / "{species}/{species}.draft.fasta.fai", species=config["species_list"]),
        expand(out_dir_path / ("{species}/{species}_%s.txt" % config["restrictase"]), species=config["species_list"]),
        expand(out_dir_path / "{species}/merged_nodups.txt", species=config["species_list"])

include: "workflow/rules/JuicerPreprocessing.smk"

