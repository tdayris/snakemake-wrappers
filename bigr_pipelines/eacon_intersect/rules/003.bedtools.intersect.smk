rule bedtools_multi_intersect:
    input:
        bed=expand("bed/{sample}.bed", sample=sample_list),
        genome="",
    output: 
        "bedtools/multi_intersect.bed",
    threads: 1
    resources:
        mem_mb=lambda wildcards, input: input.size_mb + (1024 * 4),
        time_min=lambda wildcards, input: (input.size_mb / 1024) * 30,
        tmpdir="tmp"
    params:
        extra=lambda wildcards: f"-header -names {' '.join(sample_list)}"