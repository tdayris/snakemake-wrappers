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
        body=['print $0"\\tMutect_CTC"']
    wrapper:
        "bio/awk"


rule add_origin_brc:
    input:
        "vep/brc/{sample}.tsv"
    output:
        temp("vep/brc/{sample}.orig.tsv")
    threads: 1
    resources:
        mem_mb=512,
        time_min=10,
        tmpdir=tmp,
    log:
        "logs/vep/origin/{sample}.brc.log"
    params:
        begin='FS=OFS="\\t"',
        body=['print $0"\\tBRC"']
    wrapper:
        "bio/awk"



rule add_origin_wbc_ctc:
    input:
        "vep/hc/{sample}.wbc.tsv"
    output:
        temp("vep/hc/{sample}.wbc.orig.tsv")
    threads: 1
    resources:
        mem_mb=512,
        time_min=10,
        tmpdir=tmp,
    log:
        "logs/vep/origin/{sample}.wbc.log"
    params:
        begin='FS=OFS="\\t"',
        body=['print $0"\\tHC_WBC"']
    wrapper:
        "bio/awk"


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
        body=['print $0"\\tHC_WBC"']
    wrapper:
        "bio/awk"


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
        pipe("bigtable/noheader.tsv"),
    threads: 1
    resources:
        mem_mb=get_1gb_per_attempt,
        time_min=get_15min_per_attempt,
        tmpdir=tmp,
    log:
        "logs/bigtable/noheader.log",
    params:
        "'1d;s|.ctc.brc.vcf||g;s|Baseline|Germline|g;'",
    conda:
        str(workflow_source_dir / "envs" / "bash.yaml")
    shell:
        "sed {params} {input} > {output} 2> {log}"


rule bigtable_sort:
    input:
        "bigtable/noheader.tsv",
    output:
        temp("bigtable/sorted.tsv"),
    threads: 1
    resources:
        mem_mb=get_10gb_per_attempt,
        time_min=get_45min_per_attempt,
        tmpdir=tmp,
    log:
        "logs/bigtable/sort.log",
    params:
        "-k1,1 -k2,2n",
    conda:
        str(workflow_source_dir / "envs" / "bash.yaml")
    shell:
        "sort {params} {input} > {output}"


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
