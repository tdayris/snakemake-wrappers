rule star_fusion:
    input:
        junctions="star/{sample}/{sample}.Chimeric.out.junction",
        resource_lib=config["CTAT_resource_lib"]
    output:
        directory("star-fusions/{sample}/")
    threads: 20
    resources:
        mem_mb=get_20g_per_attempt,
        time_min=get_2h_per_attempt,
        tmpdir="tmp"
    log:
        "logs/star-fusion/{sample}.log"
    params:
        extra = "--FusionInspector inspect"
    env:
        "envs/fusions.yaml"
    shell:
        "STAR-Fusion "
        "--chimeric_junction {input.junctions} "
        "--genome_lib_dir {input.resource_lib} "
        "--CPU {threads} "
        "--output_dir {output} "
        "> {log} 2>&1 "
        