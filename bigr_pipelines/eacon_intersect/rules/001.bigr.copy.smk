rule find_acbs_files:
    output:
        temp("acbs_list.txt")
    threads: 1
    resources:
        mem_mb=1024,
        time_min=lambda wildcards, attempt: attempt * 90,
        tmpdir="tmp"
    log:
        "logs/copy/find_acbs.log"
    params:
        extra='-type f -name "*Cut.csb"',
        exec_dir=config["OncoCytoDir"],
    shell:
        "find {params.exec_dir} {params.extra} > {output} 2> {log}"


checkpoint rsync_cbs:
    input:
        "acbs_list.txt"
    output:
        directory("data_input/cbs_files/")
    threads: 1
    resources:
        mem_mb=512,
        time_min=lambda wildcards, attempt: attempt * 90,
        tmpdir="tmp",
    conda:
        "../envs/rsync.yaml"
    log:
        "logs/copy/rsync_acbs.log"
    params:
        extra="--human-readable --partial --progress --verbose --checksum ",
    shell:
        "rsync {params.extra} --files-from {input} {output} > {log} 2>&1"