"""
004.unzip_interop
from
-> Entry job
by
-> 003.multiqc
"""


use rule unzip_stats as uzip_interop with:
    input:
        interop,
    output:
        temp(directory("tmp/InterOp")),
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
    input:
        runinfo,
    output:
        temp("tmp/RunInfo.xml"),
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
    input:
        runparams,
    output:
        temp("tmp/RunParameters.xml"),
    log:
        "logs/003.unzip_runparams.log",
