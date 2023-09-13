"""
No command line provided
"""

rule awk_add_origin_brc_ctc:
    input:
        "vep/brc/{sample}_{version}_{manip}_{nb}.tsv"
    output:
        temp("vep/brc/{sample}_{version}_{manip}_{nb}.orig.tsv")
    threads: 1
    resources:
        mem_mb=512,
        time_min=10,
        tmpdir=tmp,
    log:
        "logs/vep/origin/{sample}_{version}_{manip}_{nb}.brc.log"
    params:
        begin='FS=OFS="\\t"',
        body=[
            ['NR == 1', 'print $0"\\tSample_Type\\tTool\\tCondition"', 'print $0"\\tBRC_CTC\\tBRC\\tCTC"']
        ]
    wrapper:
        "bio/awk"


rule awk_add_origin_brc_wbc:
    input:
        "vep/brc/{sample}_{version}_{manip}.tsv"
    output:
        temp("vep/brc/{sample}_{version}_{manip}.orig.tsv")
    threads: 1
    resources:
        mem_mb=512,
        time_min=10,
        tmpdir=tmp,
    log:
        "logs/vep/origin/{sample}_{version}_{manip}.brc.log"
    params:
        begin='FS=OFS="\\t"',
        body=[
            ['NR == 1', 'print $0"\\tSample_Type\\tTool\\tCondition"', 'print $0"\\tBRC_WBC\\tBRC\\tWBC"']
        ]
    wrapper:
        "bio/awk"


rule awk_add_origin_brc_baseline:
    input:
        "vep/brc/{sample}.baseline.tsv"
    output:
        temp("vep/brc/{sample}.baseline.orig.tsv")
    threads: 1
    resources:
        mem_mb=512,
        time_min=10,
        tmpdir=tmp,
    log:
        "logs/vep/origin/{sample}.brc.log"
    params:
        begin='FS=OFS="\\t"',
        body=[
            ['NR == 1', 'print $0"\\tSample_Type\\tTool\\tCondition"', 'print $0"\\tBRC_Baseline\\tBRC\\tBaseline"']
        ]
    wrapper:
        "bio/awk"