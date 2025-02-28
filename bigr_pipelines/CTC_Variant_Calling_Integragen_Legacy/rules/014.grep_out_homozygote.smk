"""
grep -v "./.:0,0:0:0,0,0" {patient}_CTC_all.cutadapt.sorted.rmmarkdup.hc.snps.filtered.PASS.{sample1_normal}.vcf > {patient}_CTC_all.cutadapt.sorted.rmmarkdup.hc.snps.filtered.PASS.{sample1_normal}.grepped.vcf
"""


rule grep_out_homozygote_bseline:
    input:
        "gatk/select_variants/baseline/{sample}.tmp.vcf",
    output:
        temp("gatk/select_variants/baseline/{sample}.vcf"),
    threads: 1
    resources:
        mem_mb=get_1gb_per_attempt,
        time_min=get_35min_per_attempt,
        tmpdir=tmp,
    log:
        "logs/grep/{sample}.baseline.log",
    conda:
        str(workflow_source_dir / "envs" / "bash.yaml")
    params:
        ' -v "./.:0,0:0:0,0,0"',
    shell:
        "grep {params} {input} > {output} 2> {log}"


rule grep_out_homozygote_wbc:
    input:
        "gatk/select_variants/wbc/{sample}_{version}_{manip}.tmp.vcf",
    output:
        temp("gatk/select_variants/wbc/{sample}_{version}_{manip}.vcf"),
    threads: 1
    resources:
        mem_mb=get_1gb_per_attempt,
        time_min=get_35min_per_attempt,
        tmpdir=tmp,
    log:
        "logs/grep/{sample}_{version}_{manip}.wbc.log",
    conda:
        str(workflow_source_dir / "envs" / "bash.yaml")
    params:
        ' -v "./.:0,0:0:0,0,0"',
    shell:
        "grep {params} {input} > {output} 2> {log}"


rule grep_out_homozygote_ctc:
    input:
        "gatk/select_variants/ctc/{sample}_{version}_{manip}_{nb}.tmp.vcf",
    output:
        temp("gatk/select_variants/ctc/{sample}_{version}_{manip}_{nb}.vcf"),
    threads: 1
    resources:
        mem_mb=get_1gb_per_attempt,
        time_min=get_35min_per_attempt,
        tmpdir=tmp,
    log:
        "logs/grep/{sample}_{version}_{manip}_{nb}.ctc.log",
    conda:
        str(workflow_source_dir / "envs" / "bash.yaml")
    params:
        ' -v "./.:0,0:0:0,0,0"',
    shell:
        "grep {params} {input} > {output} 2> {log}"