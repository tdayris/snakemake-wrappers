rule star_fusion:
    input:
        junctions="star/{sample}/{sample}.{status}.Chimeric.out.junction",
        resource_lib=config["star_fusion"]["CTAT_resource_lib"],
    output:
        directory("star-fusions/{sample}/"),
    threads: config.get("max_threads", 20)
    resources:
        mem_mb=get_20gb_per_attempt,
        time_min=get_2h_per_attempt,
        tmpdir=tmp,
    retries: 1
    log:
        "logs/star-fusion/{sample}.log",
    params:
        extra=config["star_fusion"].get("extra", "--FusionInspector inspect"),
    conda:
        "../envs/fusions.yaml"
    shell:
        "STAR-Fusion "
        "--chimeric_junction {input.junctions} "
        "--genome_lib_dir {input.resource_lib} "
        "--CPU {threads} "
        "--output_dir {output} "
        "> {log} 2>&1 "
