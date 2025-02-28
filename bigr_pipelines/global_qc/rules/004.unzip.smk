# Unzip Stats.json for multiqc inclusion
"""
004.unzip_stats
from
-> Entry job
by
-> 003.multiqc
"""


rule unzip_stats:
    input:
        stats
    output:
        temp("tmp/Stats.json"),
    threads: 1
    resources:
        mem_mb=get_768mb_per_attempt,
        time_min=get_15min_per_attempt,
        tmpdir="tmp",
    params:
        wd=f"-d {os.getcwd()}/tmp",
        extra="-o"
    log:
        "logs/003.unzip_stats.log",
    shell:
        "unzip {params.extra} {params.wd} {input} > {log} 2>&1"