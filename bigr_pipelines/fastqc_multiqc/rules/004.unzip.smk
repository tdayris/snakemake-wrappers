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
    params:
        regex="input/*/archive/*/unaligned/Stats/Stats.json.zip"
    log:
        "logs/003.unzip_stats.log",
    shell:
        'unzip -n -d "${{PWD}}" '
        '{params.regex} '
        "> {log} 2>&1"


"""
004.unzip_interop
from
-> Entry job
by
-> 003.multiqc
"""


use rule unzip_stats as uzip_interop with:
    output:
        directory("InterOp"),
    params:
        regex="input/*/archive/uploadToKDIAnalysis/InterOp.zip"
    log:
        "logs/003.unzip_interop.log",


"""
004.unzip_runinfo
from
-> Entry job
by
-> 003.multiqc
"""


use rule unzip_stats as unzip_runinfo with:
    output:
        "RunInfo.xml",
    params:
        regex="input/*/archive/uploadToKDIAnalysis/RunInfo.xml.zip"
    log:
        "logs/003.unzip_runinfo.log",


"""
004.unzip_runparams
from
-> Entry job
by
-> 003.multiqc
"""


use rule unzip_stats as unzip_runparams with:
    output:
        "RunParameters.xml",
    params:
        regex="input/*/archive/uploadToKDIAnalysis/RunParameters.xml.zip"
    log:
        "logs/003.unzip_runparams.log",
