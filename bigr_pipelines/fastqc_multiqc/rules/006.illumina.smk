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
    input:
        runinfo,
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
    input:
        runparams,
    output:
        "RunParameters.xml",
    params:
        regex="input/*/archive/uploadToKDIAnalysis/RunParameters.xml.zip"
    log:
        "logs/003.unzip_runparams.log",
