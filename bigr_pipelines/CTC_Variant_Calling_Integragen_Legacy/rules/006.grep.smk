"""
grep_pass (optional) : grep -P "^#" {patient}_CTC_all.cutadapt.sorted.rmmarkdup.hc.snps.filtered.vcf >{patient}_CTC_all.cutadapt.sorted.rmmarkdup.hc.snps.filtered.PASS.vcf ; grep -P "\tPASS\t" {patient}_CTC_all.cutadapt.sorted.rmmarkdup.hc.snps.filtered.vcf >>{patient}_CTC_all.cutadapt.sorted.rmmarkdup.hc.snps.filtered.PASS.vcf
"""


rule grep:
    input:
        "gatk/variant_filtration/baseline_wbc/{sample}.g.vcf",
    output:
        temp("grep/baseline_wbc/{sample}.vcf"),
    threads: 1
    resources:
        mem_mb=get_1gb_per_attempt,
        time_min=get_35min_per_attempt,
        tmpdir=tmp,
    conda:
        str(workflow_source_dir / "envs" / "bash.yaml")
    log:
        "logs/grep/{sample}.log",
    params:
        '-P "^#|\tPASS\t"',
    shell:
        "grep {params} {input} > {output} 2> {log}"
