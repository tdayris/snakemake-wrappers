"""
No command line provided
"""

rule awk_add_origin_mutect_ctc:
    input:
        "vep/mutect/{sample}_{version}_{manip}_{nb}.tsv"
    output:
        temp("vep/mutect/{sample}_{version}_{manip}_{nb}.orig.tsv")
    threads: 1
    resources:
        mem_mb=512,
        time_min=10,
        tmpdir=tmp,
    log:
        "logs/vep/origin/{sample}_{version}_{manip}_{nb}.log"
    params:
        begin='FS=OFS="\\t"',
        body=[
            ['NR == 1', 'print $0"\\tSample_Type\\tTool\\tCondition"', 'print $0"\\tMutect_CTC\\tMutect\\tCTC"']
        ]
    wrapper:
        "bio/awk"


rule awk_add_origin_mutect_wbc:
    input:
        "vep/mutect/{sample}_{version}_{manip}.tsv"
    output:
        temp("vep/mutect/{sample}_{version}_{manip}.orig.tsv")
    threads: 1
    resources:
        mem_mb=512,
        time_min=10,
        tmpdir=tmp,
    log:
        "logs/vep/origin/{sample}.mutect2_wbc.log"
    params:
        begin='FS=OFS="\\t"',
        body=[
            ['NR == 1', 'print $0"\\tSample_Type\\tTool\\tCondition"', 'print $0"\\tMutect_WBC\\tMutect\\tWBC"']
        ]
    wrapper:
        "bio/awk"