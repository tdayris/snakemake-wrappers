rule link_fasta_for_vep:
    input:
        config["ref"]["fasta"],
    output:
        temp("resources/GRCh38.fasta"),
    threads: 1
    resources:
        mem_mb=get_1gb_per_attempt,
        time_min=get_15min_per_attempt,
        tmpdir=tmp,
    conda:
        str(workflow_source_dir / "envs" / "bash.yaml")
    log:
        "logs/vep/link_fasta.log",
    params:
        "--force --relative --symbolic --verbose",
    shell:
        "ln {params} {input} {output} > {log} 2>&1"