rule bwa:
    input:
        index_bwt=draft_dir_path / "{species}/{species}.draft.fasta.bwt",
        forward_read=fastq_dir_path / "{species}/{chunk}_R1_{suffix}.fastq.gz",
        reverse_read=fastq_dir_path / "{species}/{chunk}_R2_{suffix}.fastq.gz" #TODO: add function to work with both compressed and not compressed files
    output:
        sam=out_dir_path / "{species}/{chunk}_{suffix}.sam",
    params:
        index_prefix=draft_dir_path / "{species}/{species}.draft.fasta"
    log:
        std=log_dir_path / "{species}/bwa.{chunk}_{suffix}.log",
        cluster_log=cluster_log_dir_path / "{species}.bwa.{chunk}_{suffix}.cluster.log",
        cluster_err=cluster_log_dir_path / "{species}.bwa.{chunk}_{suffix}.cluster.err"
    benchmark:
        benchmark_dir_path / "{species}/bwa.{chunk}_{suffix}.benchmark.txt"
    conda:
        "../../%s" % config["conda_config"]
    resources:
        cpus=config["bwa_threads"],
        time=config["bwa_time"],
        mem=config["bwa_mem_mb"],
        partition=config["bwa_partition"] if config["bwa_partition"] else config["default_partition"]
    threads:
        config["bwa_threads"]
    shell:
        "bwa mem -SP5M ${threads} ${params.index_prefix} <(zcat {input.forward_read}) <(zcat {input.reverse_read}) > {output} 2>{log.std}"

