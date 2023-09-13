rule awk_add_origin_hc_ctc:
    input:
        "vep/hc/{sample}_{version}_{manip}_{nb}.tsv"
    output:
        temp("vep/hc/{sample}_{version}_{manip}_{nb}.orig.tsv")
    threads: 1
    resources:
        mem_mb=512,
        time_min=10,
        tmpdir=tmp,
    log:
        "logs/vep/origin/{sample}_{version}_{manip}_{nb}.ctc.log"
    params:
        begin='FS=OFS="\\t"',
        body=[
            ['NR == 1', 'print $0"\\tSample_Type\\tTool\\tCondition"', 'print $0"\\tHC_CTC\\tHaplotypeCaller\\tCTC"']
        ]
    wrapper:
        "bio/awk"


rule awk_add_origin_hc_wbc:
    input:
        "vep/hc/{sample}_{version}_{manip}.tsv"
    output:
        temp("vep/hc/{sample}_{version}_{manip}.orig.tsv")
    threads: 1
    resources:
        mem_mb=512,
        time_min=10,
        tmpdir=tmp,
    log:
        "logs/vep/origin/{sample}_{version}_{manip}.wbc.log"
    params:
        begin='FS=OFS="\\t"',
        body=[
            ['NR == 1', 'print $0"\\tSample_Type\\tTool\\tCondition"', 'print $0"\\tHC_WBC\\tHaplotypeCaller\\tWBC"']
        ]
    wrapper:
        "bio/awk"


rule awk_add_origin_hc_ctc:
    input:
        "vep/hc/{sample}.baseline.tsv"
    output:
        temp("vep/hc/{sample}.baseline.orig.tsv")
    threads: 1
    resources:
        mem_mb=512,
        time_min=10,
        tmpdir=tmp,
    log:
        "logs/vep/origin/{sample}.baseline.log"
    params:
        begin='FS=OFS="\\t"',
        body=[
            ['NR == 1', 'print $0"\\tSample_Type\\tTool\\tCondition"', 'print $0"\\tHC_Baseline\\tHaplotypeCaller\\tBaseline"']
        ]
    wrapper:
        "bio/awk"