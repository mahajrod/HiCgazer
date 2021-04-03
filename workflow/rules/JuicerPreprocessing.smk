rule bwa_index:
    input:
        draft_dir_path / "{species}/{species}.draft.fasta"
    output:
        draft_dir_path / "{species}/{species}.draft.fasta.bwt"
    log:
        std=log_dir_path / "{species}/bwa_index.log",
        cluster_log=cluster_log_dir_path / "{species}.bwa_index.cluster.log",
        cluster_err=cluster_log_dir_path / "{species}.bwa_index.cluster.err"
    benchmark:
        benchmark_dir_path / "{species}/bwa_index.benchmark.txt"
    conda:
        "../../%s" % config["conda_config"]
    resources:
        cpus=config["bwa_index_threads"],
        time=config["bwa_index_time"],
        mem=config["bwa_index_mem_mb"],
        partition=config["bwa_index_partition"] if config["bwa_index_partition"] else config["default_partition"]
    threads:
        config["bwa_index_threads"]
    shell:
        "bwa index {input} > {log.std} 2>&1"

rule ref_faidx:
    input:
        draft_dir_path / "{species}/{species}.draft.fasta"
    output:
        draft_dir_path / "{species}/{species}.draft.fasta.fai"
    log:
        std=log_dir_path / "{species}/faidx.log",
        cluster_log=cluster_log_dir_path / "{species}.faidx.cluster.log",
        cluster_err=cluster_log_dir_path / "{species}.faidx.cluster.err"
    benchmark:
        benchmark_dir_path / "{species}/faidx.benchmark.txt"
    conda:
        "../../%s" % config["conda_config"]
    resources:
        cpus=config["faidx_threads"],
        time=config["faidx_time"],
        mem=config["faidx_mem_mb"],
        partition=config["faidx_partition"] if config["faidx_partition"] else config["default_partition"]
    threads:
        config["faidx_threads"]
    shell:
         "samtools faidx {input} > {log.std} 2>&1"

rule restriction_sites:
    input:
        draft_dir_path / "{species}/{species}.draft.fasta"
    output:
        out_dir_path / ("{species}/{species}_%s.txt" % config["restrictase"])
    params:
        restrictase=config["restrictase"],
        output_prefix=out_dir_path / "{species}/{species}"
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
         "generate_site_positions.py {params.restrictase} {params.output_prefix} {input}> {log.std} 2>&1"

rule juicer:
    input:
        draft=draft_dir_path / "{species}/{species}.draft.fasta",
        draft_index=rules.bwa_index.output,
        restriction_sites=rules.restriction_sites.output,
        species_dir=directory(draft_dir_path / "{species}"),
        fastq_dir=(fastq_dir_path / "{species}").absolute()
    output:
        merged=out_dir_path / "{species}/merged_nodups.txt",
        chr_path=out_dir_path / "{species}/{species}.chr"
    params:
        restrictase=config["restrictase"],
        output_prefix=out_dir_path / "{species}/{species}",
        species_fastq_dir= out_dir_path / "{species}/fastq"
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
         "ln {input.fastq_dir} {params.species_fastq_dir}; juicer.sh -q {resources.partition}  -Q {resources.time} -L {resources.time} "
         " -t {threads} -g {wildcards.species}"
         " -d {input.species_dir} –g {params.output_prefix} "
         " –s {params.restrictase} –z {input.draft} –y {input.restriction_sites} "
         " –p {output.chr_path} > {log.std} 2>&1"
