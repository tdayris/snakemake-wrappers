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
        temp("Stats.json"),
    threads: 1
    resources:
        mem_mb=get_768mb_per_attempt,
        time_min=get_15min_per_attempt,
        tmpdir="tmp",
    log:
        "logs/003.unzip_stats.log",
    shell:
        'unzip -n -d "${{PWD}}" {input} > {log} 2>&1'