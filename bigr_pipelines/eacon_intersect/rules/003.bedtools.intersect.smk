def acbs_file_list(wildcards):
    results = {"genome": config["ref"]["fasta_index"]}
    bedfiles = checkpoints.rsync_cbs.get(**wildcards).output[0]
    results["bed"] = expand(
        "bed/{sample}.bed",
        sample=glob_wildcards(os.path.join(bedfiles, "{sample}.Cut.cbs")).sample
    )
    return results



rule bedtools_multi_intersect:
    input:
        acbs_file_list
    output: 
        "bedtools/multi_intersect.bed",
    threads: 1
    resources:
        mem_mb=lambda wildcards, input: input.size_mb + (1024 * 4),
        time_min=lambda wildcards, input: (input.size_mb / 1024) * 30,
        tmpdir="tmp"
    params:
        extra="-header"
    conda:
        "../envs/bedtools.yaml"
    shell:
        "bedtools multiinter "
        "{params.extra} "
        "-g {input.genome} "
        "-i {input.bed} "
        "> {output} 2>&1"