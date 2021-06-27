#template
snakemake --profile profile/slurm/ --configfile config/default.yaml --config draft_dir= fastq_dir= --use-conda --printshellcmds --latency-wait 60

#example

snakemake --cores 32 --profile profile/slurm/ --configfile config/default.yaml --config default_partition=main --use-conda --printshellcmds --latency-wait 60

