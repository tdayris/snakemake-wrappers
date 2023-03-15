rule add_origin_mutect_ctc:
    input:
        "vep/mutect/{sample}.tsv"
    output:
        temp("vep/mutect/{sample}.orig.tsv")
    threads: 1
    resources:
        mem_mb=512,
        time_min=10,
        tmpdir=tmp,
    log:
        "logs/vep/origin/{sample}.mutect.log"
    params:
        begin='FS=OFS="\\t"',
        body=[
            ['NR == 1', 'print $0"\\tSample_Type\\tTool\\tCondition"', 'print $0"\\tMutect_CTC\\tMutect\\tCTC"']
            # 'print $0"\\tMutect_CTC"'
        ]
    wrapper:
        "bio/awk"


rule add_origin_brc:
    input:
        "vep/bcr/{sample}.tsv"
    output:
        temp("vep/bcr/{sample}.orig.tsv")
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
            ['NR == 1', 'print $0"\\tSample_Type\\tTool\\tCondition"', 'print $0"\\tBRC_CTC\\tBRC\\tCTC"']
            # 'print $0"\\tBRC_CTC"'
        ]
    wrapper:
        "bio/awk"


rule add_origin_wbc_ctc:
    input:
        "vep/hc/{sample}.wbc.tsv"
    output:
        temp("vep/hc/{sample}.wbc.tmp.tsv")
    threads: 1
    resources:
        mem_mb=512,
        time_min=10,
        tmpdir=tmp,
    log:
        "logs/vep/origin/{sample}.wbc.log"
    group:
        "wbc_origin"
    params:
        begin='FS=OFS="\\t"',
        body=[
            ['NR == 1', 'print $0"\\tSample_Type\\tTool\\tCondition"', 'print $0"\\tHC_WBC\\tHaplotypeCaller\\tWBC"']
            # 'print $0"\\tHC_WBC"'
        ]
    wrapper:
        "bio/awk"


rule rename_sample_wbc_hc:
    input:
        "vep/hc/{sample}.wbc.tmp.tsv"
    output:
        temp("vep/hc/{sample}.wbc.orig.tsv")
    threads: 1
    resources:
        mem_mb=512,
        time_min=10,
        tmpdir=tmp,
    log:
        "logs/vep/rename/{sample}.wbc.log"
    group:
        "wbc_origin"
    params:
        sample=lambda wildcards: (
            "_".join(str(wildcards.sample).split("_")[:-1]) + "_WBC" 
            if re.search("_\d+$", str(wildcards.sample)) else str(wildcards.sample) + "_WBC"
        )
    conda:
        "../envs/bash.yaml"
    shell:
        "sed 's|{wildcards.sample}|{params.sample}|g' {input} > {output} 2> {log}"


rule add_origin_hc_ctc:
    input:
        "vep/hc/{sample}.ctc.tsv"
    output:
        temp("vep/hc/{sample}.ctc.orig.tsv")
    threads: 1
    resources:
        mem_mb=512,
        time_min=10,
        tmpdir=tmp,
    log:
        "logs/vep/origin/{sample}.ctc.log"
    params:
        begin='FS=OFS="\\t"',
        body=[
            ['NR == 1', 'print $0"\\tSample_Type\\tTool\\tCondition"', 'print $0"\\tHC_CTC\\tHaplotypeCaller\\tCTC"']
            # 'print $0"\\tHC_CTC"'
        ]
    wrapper:
        "bio/awk"


rule add_origin_hc_bseline:
    input:
        "vep/hc/{sample}.baseline.tsv"
    output:
        temp("vep/hc/{sample}.baseline.tmp.tsv")
    threads: 1
    resources:
        mem_mb=512,
        time_min=10,
        tmpdir=tmp,
    log:
        "logs/vep/origin/{sample}.baseline.log"
    group:
        "baseline_origin"
    params:
        begin='FS=OFS="\\t"',
        body=[
            ['NR == 1', 'print $0"\\tSample_Type\\tTool\\tCondition"', 'print $0"\\tHC_Germline\\tHaplotypeCaller\\tGermline"']
            # 'print $0"\\tHC_Germline"'
        ]
    wrapper:
        "bio/awk"


rule rename_sample_germline_hc:
    input:
        "vep/hc/{sample}.baseline.tmp.tsv"
    output:
        temp("vep/hc/{sample}.baseline.orig.tsv")
    threads: 1
    resources:
        mem_mb=512,
        time_min=10,
        tmpdir=tmp,
    log:
        "logs/vep/rename/{sample}.baseline.log"
    group:
        "baseline_origin"
    params:
        sample=lambda wildcards: (
            "_".join(str(wildcards.sample).split("_")[:-1]) + "_Germline_DNA" 
            if re.search("_\d+$", str(wildcards.sample)) else str(wildcards.sample) + "_Germline_DNA"
        )
    conda:
        "../envs/bash.yaml"
    shell:
        "sed 's|{wildcards.sample}|{params.sample}|g' {input} > {output} 2> {log}"


rule concat_to_bigtable:
    input:
        expand(
            "vep/{annot}/{sample}.orig.tsv",
            annot=["bcr", "mutect"],
            sample=samples_list,
        ),
        expand(
            "vep/hc/{sample}.wbc.orig.tsv",
            sample=samples_list,
        ),
        expand(
            "vep/hc/{sample}.ctc.orig.tsv",
            sample=samples_list,
        ),
        expand(
            "vep/hc/{sample}.baseline.orig.tsv",
            sample=samples_list,
        ),
    output:
        temp("bigtable/raw.tsv"),
    threads: 1
    resources:
        mem_mb=get_10gb_per_attempt,
        time_min=get_45min_per_attempt,
        tmpdir=tmp,
    log:
        "logs/bigtable/raw.log",
    conda:
        str(workflow_source_dir / "envs" / "bash.yaml")
    params:
        "'(NR == 1) || (FNR > 1)'",
    shell:
        "awk {params} {input} > {output} 2> {log}"


rule bigtable_header:
    input:
        "bigtable/raw.tsv",
    output:
        temp("bigtable/header.tsv"),
    threads: 1
    resources:
        mem_mb=get_1gb_per_attempt,
        time_min=get_15min_per_attempt,
        tmpdir=tmp,
    log:
        "logs/bigtable/header.log",
    params:
        "-n 1",
    conda:
        str(workflow_source_dir / "envs" / "bash.yaml")
    shell:
        "head {params} {input} > {output} 2> {log}"


rule bigtable_noheader:
    input:
        "bigtable/raw.tsv",
    output:
        temp("bigtable/noheader.tsv"),
    threads: 1
    resources:
        mem_mb=get_1gb_per_attempt,
        time_min=get_15min_per_attempt,
        tmpdir=tmp,
    log:
        "logs/bigtable/noheader.log",
    params:
        "'1d;s|.ctc.brc.vcf||g;s|.ctc.bcr.vcf||g;s|Baseline|Germline|g'",
    conda:
        str(workflow_source_dir / "envs" / "bash.yaml")
    shell:
        "sed {params} {input} > {output} 2> {log}"


rule bigtable_sort:
    input:
        "bigtable/noheader.tsv",
    output:
        temp("bigtable/sorted.tsv"),
    threads: 2
    resources:
        mem_mb=get_10gb_per_attempt,
        time_min=get_45min_per_attempt,
        tmpdir=tmp,
    log:
        "logs/bigtable/sort.log",
    params:
        "",
    conda:
        str(workflow_source_dir / "envs" / "bash.yaml")
    shell:
        "sort {params} {input} | uniq > {output}"


rule bigtable_output:
    input:
        header="bigtable/header.tsv",
        content="bigtable/sorted.tsv",
    output:
        protected("data_output/bigtable.tsv"),
    threads: 1
    resources:
        mem_mb=get_1gb_per_attempt,
        time_min=get_15min_per_attempt,
        tmpdir=tmp,
    log:
        "logs/bigtable/output.log",
    params:
        "",
    conda:
        str(workflow_source_dir / "envs" / "bash.yaml")
    shell:
        "cat {input.header} {input.content} > {output} 2> {log}"



rule remove_duplicates:
    input:
        "data_output/bigtable.tsv"
    output:
        temp("bigtable/bigtable.uniq.csv")
    threads: 3
    resources:
        mem_mb=get_1gb_per_attempt,
        time_min=get_15min_per_attempt,
        tmpdir=tmp,
    log:
        "logs/xsv/preset/bigtable.log"
    params:
        s="",
        u="",
        e="'s|,|;|g;s|\t|,|g'",
    conda:
        str(workflow_source_dir / "envs" / "bash.yaml")
    shell:
        "sed {params.e} {input} | xsv sort {params.s} | uniq {params.u} > {output} 2> {log}"


rule bigtable_annotated:
    input:
        bigtable="bigtable/bigtable.uniq.csv",
        label="Labels.csv"
    output:
        bittable="data_output/bigtable.annot.tsv"
    threads: 1
    resources:
        mem_mb=get_1gb_per_attempt,
        time_min=get_15min_per_attempt,
        tmpdir=tmp,
    log:
        "logs/bigtable/annot.log",
    params:
        "",
    conda:
        str(workflow_source_dir / "envs" / "bash.yaml")
    script:
        "../scripts/annotate_bigtable.py"