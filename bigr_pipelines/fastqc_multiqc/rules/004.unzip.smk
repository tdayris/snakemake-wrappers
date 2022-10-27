# Unzip Stats.json for multiqc inclusion
"""
004.unzip_stats
from
-> Entry job
by
-> 003.multiqc
"""


rule unzip_stats:
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
        'unzip -n -d "${{PWD}}" '
        "input/*/archive/*/unaligned/Stats/Stats.json.zip "
        "> {log} 2>&1"


"""
004.unzip_interop
from
-> Entry job
by
-> 003.multiqc
"""


rule uzip_interop:
    input:
        "InterOp.zip",
    output:
        directory("InterOp"),
    threads: 1
    resources:
        mem_mb=get_768mb_per_attempt,
        time_min=get_15min_per_attempt,
        tmpdir="tmp",
    log:
        "logs/003.unzip_interop.log",
    shell:
        'unzip -n -d "${{PWD}}" {input} > {log} 2>&1'


"""
004.unzip_runinfo
from
-> Entry job
by
-> 003.multiqc
"""


use rule uzip_interop as unzip_runinfo with:
    input:
        "RunInfo.xml.zip",
    output:
        "RunInfo.xml",
    log:
        "logs/003.unzip_runinfo.log",


"""
004.unzip_runparams
from
-> Entry job
by
-> 003.multiqc
"""


use rule uzip_interop as unzip_runparams with:
    input:
        "RunParameters.xml.zip",
    output:
        "RunParameters.xml",
    log:
        "logs/003.unzip_runparams.log",
