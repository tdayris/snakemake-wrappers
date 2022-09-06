rule mixcr_extend:
    input:
        "mixcr/assemblePartial/clna2/{sample}.vdjca"
    output:
        report=temp("mixcr/extend/report/{sample}.txt"),
        json=temp("mixcr/extend/report/{sample}.log"),
        vdjca=temp("mixcr/extend/vdjca/{sample}.vdjca"),
    threads: config.get("max_threads", 20)
    resources:
        mem_mb=get_8gb_per_attempt,
        time_min=get_45min_per_attempt,
        tmpdir="tmp"
    params:
        extra=config["mixcr"].get("extend", "")
    log:
        "logs/mixcr/extend/{sample}.log"
    wrapper:
        "bio/mixcr/extend"


rule mixcr_assemble:
    input:
        "mixcr/extend/vdjca/{sample}.vdjca"
    output:
        report=temp("mixcr/assemble/report/{sample}.txt"),
        json=temp("mixcr/assemble/report/{sample}.json"),
        clna=temp("mixcr/assemble/clna/{sample}.clna"),
    threads: config.get("max_threads", 20)
    resources:
        mem_mb=get_10gb_per_attempt,
        time_min=get_45min_per_attempt,
        tmpdir="tmp"
    params:
        extra=config["mixcr"].get("assemble", "-OseparateByV=true -OseparateByJ=true")
    log:
        "logs/mixcr/extend/{sample}.log"
    wrapper:
        "bio/mixcr/extend"


rule mixcr_assemble_contigs:
    input:
        "mixcr/assemble/clna/{sample}.clna"
    output:
        report=temp("mixcr/assembleContigs/report/{sample}.txt"),
        json=temp("mixcr/assembleContigs/report/{sample}.json"),
        clns=temp("mixcr/assembleContigs/clns/{sample}.clns"),
    threads: config.get("max_threads", 20)
    resources:
        mem_mb=get_15gb_per_attempt,
        time_min=get_3h_per_attempt,
        tmpdir="tmp"
    log:
        "logs/mixcr/assembleContigs/{sample}.log"
    params:
        extra=config["mixcr"].get("assemble_contigs", "")
    wrapper:
        "bio/mixcr/assembleContigs"