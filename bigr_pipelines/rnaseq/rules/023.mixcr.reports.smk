rule mixcr_export_clone:
    input:
        "mixcr/assembleContigs/clns/{sample}.clns",
    output:
        table=protected("data_output/Mixcr/{sample}.tsv"),
    threads: 1
    resources:
        mem_mb=get_2gb_per_attempt,
        time_min=get_1h_per_attempt,
        tmpdir="tmp",
    log:
        "logs/mixcr/export_clones/{sample}.log",
    params:
        extra=config["mixcr"].get("export_clones", "--preset full"),
    wrapper:
        "bio/mixcr/export"


rule mixcr_export_qc_align:
    input:
        "mixcr/extend/vdjca/{sample}.vdjca",
    output:
        protected("data_output/QC/{sample}.alignQc.pdf"),
    threads: 1
    resources:
        mem_mb=get_4gb_per_attempt,
        time_min=get_45min_per_attempt,
        tmpdir="tmp",
    log:
        "logs/mixcr/exportQc/{sample}.align.log",
    params:
        extra=config["mixcr"].get("exportqc_align", ""),
        mode="align",
    wrapper:
        "bio/mixcr/exportqc"


rule mixcr_export_qc_chain:
    input:
        "mixcr/assembleContigs/clns/{sample}.clns",
    output:
        protected("data_output/QC/{sample}.chainUsage.pdf"),
    threads: 1
    resources:
        mem_mb=get_4gb_per_attempt,
        time_min=get_45min_per_attempt,
        tmpdir="tmp",
    log:
        "logs/mixcr/exportQc/{sample}.align.log",
    params:
        extra=config["mixcr"].get("exportqc_chain", "--hide-non-functional"),
        mode="chainUsage",
    wrapper:
        "bio/mixcr/exportqc"
