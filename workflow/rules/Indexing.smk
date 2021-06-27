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
