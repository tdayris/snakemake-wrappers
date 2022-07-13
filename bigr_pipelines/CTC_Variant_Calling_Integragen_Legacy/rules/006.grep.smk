"""
grep_pass (optional) : grep -P "^#" {patient}_CTC_all.cutadapt.sorted.rmmarkdup.hc.snps.filtered.vcf >{patient}_CTC_all.cutadapt.sorted.rmmarkdup.hc.snps.filtered.PASS.vcf ; grep -P "\tPASS\t" {patient}_CTC_all.cutadapt.sorted.rmmarkdup.hc.snps.filtered.vcf >>{patient}_CTC_all.cutadapt.sorted.rmmarkdup.hc.snps.filtered.PASS.vcf
"""
rule grep:
    input:
        "gatk/variant_filtration/baseline_wbc/{sample}.g.vcf"
    output:
        temp("grep/baseline_wbc/{sample}.vcf")
    threads: 1
    resource:
        mem_mb = lambda wildcards, attempt: attempt * 512,
        time_min = lambda wildcards, attempt: attempt * 35,
        tmpdir = "tmp"
    log:
        "logs/grep/{sample}.log"
    params:
        '-P "^#|\tPASS\t"'
    shell:
        "grep {params} {input} > {output} 2> {log}"