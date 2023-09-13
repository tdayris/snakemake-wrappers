"""
grep_pass (optional) : grep -P "^#" {patient}_CTC_all.cutadapt.sorted.rmmarkdup.hc.snps.filtered.vcf >{patient}_CTC_all.cutadapt.sorted.rmmarkdup.hc.snps.filtered.PASS.vcf ; grep -P "\tPASS\t" {patient}_CTC_all.cutadapt.sorted.rmmarkdup.hc.snps.filtered.vcf >>{patient}_CTC_all.cutadapt.sorted.rmmarkdup.hc.snps.filtered.PASS.vcf
"""


rule grep_pass_wbc_baseline:
    input:
        "gatk/variant_filtration/{sample}_{version}_{manip}.g.vcf",
    output:
        temp("grep/{sample}_{version}_{manip}.vcf"),
    threads: 1
    resources:
        mem_mb=get_1gb_per_attempt,
        time_min=get_35min_per_attempt,
        tmpdir=tmp,
    conda:
        str(workflow_source_dir / "envs" / "bash.yaml")
    log:
        "logs/grep/{sample}_{version}_{manip}.log",
    params:
        '-P "^#|\tPASS\t"',
    shell:
        "grep {params} {input} > {output} 2> {log}"


rule grep_pass_ctc_baseline:
    input:
        "gatk/variant_filtration/{sample}_{version}_{manip}_{nb}.g.vcf",
    output:
        temp("grep/{sample}_{version}_{manip}_{nb}.vcf"),
    threads: 1
    resources:
        mem_mb=get_1gb_per_attempt,
        time_min=get_35min_per_attempt,
        tmpdir=tmp,
    conda:
        str(workflow_source_dir / "envs" / "bash.yaml")
    log:
        "logs/grep/{sample}_{version}_{manip}_{nb}.log",
    params:
        '-P "^#|\tPASS\t"',
    shell:
        "grep {params} {input} > {output} 2> {log}"


rule grep_pass_baseline:
    input:
        "gatk/variant_filtration/{sample}.baseline.g.vcf",
    output:
        temp("grep/{sample}.baseline.vcf"),
    threads: 1
    resources:
        mem_mb=get_1gb_per_attempt,
        time_min=get_35min_per_attempt,
        tmpdir=tmp,
    conda:
        str(workflow_source_dir / "envs" / "bash.yaml")
    log:
        "logs/grep/{sample}.baseline.log",
    params:
        '-P "^#|\tPASS\t"',
    shell:
        "grep {params} {input} > {output} 2> {log}"