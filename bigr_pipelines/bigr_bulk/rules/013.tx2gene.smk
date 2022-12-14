rule tx2gene:
    input:
        "resources/{genome_build}.{genome_release}.gtf"
    output:
        tx2gene_small="resources/{genome_build}.{genome_release}.tx2gene_small.tsv",
    threads: 1
    resources:
        mem_mb=get_3gb_per_gb,
        time_min=get_5min_per_attempt,
        tmpdir="tmp",
    cache: True
    log:
        "logs/013.tx2gene/{genome_build}.{genome_release}.log"
    conda:
        str(workflow_source_dir / "envs" / "013.pandas.yaml")
    script:
        str(workflow_source_dir / "scripts" / "013.tx2gene.py")
    