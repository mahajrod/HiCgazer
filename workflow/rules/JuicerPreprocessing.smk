
rule restriction_sites:
    input:
        draft_dir_path / "{species}/{species}.draft.fasta"
    output:
        out_dir_path / ("{species}/{species}_%s.txt" % config["restrictase"])
    params:
        restrictase=config["restrictase"],
        output_prefix=species_output_prefix_from_wildcards # out_dir_path / "{species}/{species}"
    log:
        std=log_dir_path / "{species}/restriction_sites.log",
        cluster_log=cluster_log_dir_path / "{species}.restriction_sites.cluster.log",
        cluster_err=cluster_log_dir_path / "{species}.restriction_sites.cluster.err"
    benchmark:
        benchmark_dir_path / "{species}/restriction_sites.benchmark.txt"
    conda:
        "../../%s" % config["conda_config"]
    resources:
        cpus=config["restriction_sites_threads"],
        time=config["restriction_sites_time"],
        mem=config["restriction_sites_mem_mb"],
        partition=config["restriction_sites_partition"] if config["restriction_sites_partition"] else config["default_partition"]
    threads:
        config["restriction_sites_threads"]
    shell:
         "workflow/scripts/generate_site_positions.py {params.restrictase} {params.output_prefix} {input} > {log.std} 2>&1"

rule juicer:
    input:
        draft=draft_dir_path / "{species}/{species}.draft.fasta",
        draft_index=rules.bwa_index.output,
        restriction_sites=rules.restriction_sites.output,
        fastq_dir=(fastq_dir_path / "{species}").absolute()
    output:
        merged=out_dir_path / "{species}/merged_nodups.txt",
        chr_path=out_dir_path / "{species}/{species}.chr"
    params:
        restrictase=config["restrictase"],
        output_prefix=species_output_prefix_from_wildcards , #out_dir_path / "{species}/{species}",
        species_dir=species_dir_from_wildcards, #out_dir_path / "{species}",
        species_fastq_dir=species_fastq_dir_from_wildcards #out_dir_path / "{species}/fastq"
    log:
        std=log_dir_path / "{species}/juicer.log",
        cluster_log=cluster_log_dir_path / "{species}.juicer.cluster.log",
        cluster_err=cluster_log_dir_path / "{species}.juicer.cluster.err"
    benchmark:
        benchmark_dir_path / "{species}/juicer.benchmark.txt"
    conda:
        "../../%s" % config["conda_config"]
    resources:
        cpus=config["juicer_threads"],
        time=config["juicer_time"],
        mem=config["juicer_mem_mb"],
        partition=config["juicer_partition"] if config["juicer_partition"] else config["default_partition"]
    threads:
        config["juicer_threads"]
    shell:
         "ln -s {input.fastq_dir} {params.species_fastq_dir}; "
         " juicer.sh -q {resources.partition} -l {resources.partition}  -Q {resources.time} -L {resources.time} "
         " -t {threads} -g {wildcards.species}"
         " -d {params.species_dir} -g {params.output_prefix} "
         " -s {params.restrictase} -z {input.draft} -y {input.restriction_sites} "
         " -p {output.chr_path} > {log.std} 2>&1"




