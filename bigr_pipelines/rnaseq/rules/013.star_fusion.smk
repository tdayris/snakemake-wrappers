"""
013.star_fusion
"""
rule 013_star_fusion:
    input:
        junctions="010.star/{sample}/chimera/{sample}.Chimeric.out.junction",
        resource_lib=config["star_fusion"]["CTAT_resource_lib"],
    output:
        directory("011.star-fusions/{sample}/"),
    threads: config.get("max_threads", 20)
    resources:
        mem_mb=get_20gb_per_attempt,
        time_min=get_2h_per_attempt,
        tmpdir="tmp",
    retries: 1
    log:
        "logs/011.star-fusion/{sample}.log",
    params:
        extra=config["star_fusion"].get("extra", "--FusionInspector inspect"),
    conda:
        str(workflow_source_dir / "envs" / "fusions.yaml")
    shell:
        "STAR-Fusion "
        "--chimeric_junction {input.junctions} "
        "--genome_lib_dir {input.resource_lib} "
        "--CPU {threads} "
        "--output_dir {output} "
        "> {log} 2>&1 "
