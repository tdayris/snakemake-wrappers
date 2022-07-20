rule seacr_callpeak:
    input:
        bg="bedtools/genomecov/{sample}.bedgraph",
    output:
        "seacr/{sample}.{mode}.bed",
    threads: 1
    resources:
        mem_mb=get_1gb_per_attempt,
        time_min=get_1h_per_attempt,
        tmpdir="tmp",
    log:
        "logs/seacr/{sample}.stringent.log",
    params:
        fdr=0.01,
        norm="norm",
        method=lambda wildcards: wildcards.mode,
    conda:
        "../envs/seacr.yaml"
    shell:
        "SEACR_1.3.sh {input.bg} {params.fdr} {params.norm} {wildcards.mode} seacr/{wildcards.sample} > {log} 2>&1"
