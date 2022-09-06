rule mixcr_align:
    input:
        expand("fastp/trimmed/{sample}.{stream}.fastq", stream=streams, allow_missing=True)
    output:
        json=temp("mixcr/align/reports/{sample}.json"),
        report=temp("mixcr/align/reports/{sample}.txt"),
        vdjca=temp("mixcr/align/vdjca/{sample}.vdjca"),
    threads: config("max_threads", 20)
    resources:
        mem_mb=get_4gb_per_attempt,
        time_min=get_3h_per_attempt,
        tmpdir="tmp"
    params:
        extra=config["mixcr"].get("align", "-s hsa -p rna-seq -OvParameters.geneFeatureToAlign=VTranscriptWithout5UTRWithP -OallowPartialAlignments=true")
    log:
        "logs/mixcr/align/{sample}.log"
    wrapper:
        "bio/mixcr/align"


rule mixcr_assemble_partial_round_one:
    input:
        "mixcr/align/vdjca/{sample}.vdjca"
    output:
        report=temp("mixcr/assemblePartial/report1/{sample}.txt"),
        json=temp("mixcr/assemblePartial/report1/{sample}.json"),
        clna=temp("mixcr/assemblePartial/clna1/{sample}.clna"),
        vdjca=temp("mixcr/assemblePartial/clna1/{sample}.vdjca"),
    threads: 1
    resources:
        mem_mb=get_4gb_per_attempt,
        time_min=get_3h_per_attempt,
        tmpdir="tmp"
    params:
        extra=config["mixcr"].get("assemblePartial", "")
    log:
        "logs/mixcr/assemble/{sample}.round1.log"
    wrapper:
        "bio/mixcr/assemblepartial"


rule mixcr_assemble_partial_round_two:
    input:
        "mixcr/assemblePartial/clna1/{sample}.vdjca"
    output:
        report=temp("mixcr/assemblePartial/report2/{sample}.txt"),
        json=temp("mixcr/assemblePartial/report2/{sample}.json"),
        clna=temp("mixcr/assemblePartial/clna2/{sample}.clna"),
        vdjca=temp("mixcr/assemblePartial/clna2/{sample}.vdjca"),
    threads: 1
    resources:
        mem_mb=get_4gb_per_attempt,
        time_min=get_3h_per_attempt,
        tmpdir="tmp"
    params:
        extra=config["mixcr"].get("assemblePartial", "")
    log:
        "logs/mixcr/assemble/{sample}.round2.log"
    wrapper:
        "bio/mixcr/assemblepartial"
