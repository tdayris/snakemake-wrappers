rule md5_hash:
    input:
        "{file}",
    output:
        protected("{file}.md5"),
    threads: 1
    resources:
        mem_mb=get_1gb_per_attempt,
        time_min=get_15min_per_attempt,
        tmpdir=tmp,
    retries: 1
    log:
        "logs/md5sum/{file}.log",
    params:
        extra="",
    conda:
        "../envs/bash.yaml"
    shell:
        "md5sum {params.extra} {input} > {output} 2> {log}"
