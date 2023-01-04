def bed_file_list(wildcards):
    bedfiles = checkpoints.rsync_cbs.get(**wildcards).output[0]
    return expand(
        "bed/{sample}.bed",
        sample=glob_wildcards(os.path.join(bedfiles, "data_input/cbs_files/{sample}.Cut.cbs")).sample
    )



rule bedtools_multi_intersect:
    input:
        genome=config["ref"]["fasta_index"],
        bed="data_input/cbs_files"
    output: 
        "bedtools/multi_intersect.bed",
    threads: 1
    resources:
        mem_mb=lambda wildcards, input: input.size_mb + (1024 * 4),
        time_min=get_2h_per_attempt,
        tmpdir="tmp"
    params:
        extra="-header"
    conda:
        "../envs/bedtools.yaml"
    shell:
        "bedtools multiinter "
        "{params.extra} "
        "-g {input.genome} "
        "-i {input.bed}/*.bed "
        "> {output} 2>&1"
