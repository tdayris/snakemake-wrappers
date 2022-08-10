"""
grep -v "./.:0,0:0:0,0,0" {patient}_CTC_all.cutadapt.sorted.rmmarkdup.hc.snps.filtered.PASS.{sample1_normal}.vcf > {patient}_CTC_all.cutadapt.sorted.rmmarkdup.hc.snps.filtered.PASS.{sample1_normal}.grepped.vcf
"""


rule grep_out_homozygote:
    input:
        "gatk/select_variants/baseline/{sample}.tmp.vcf",
    output:
        temp("gatk/select_variants/baseline/{sample}.vcf"),
    threads: 1
    resources:
        mem_mb=get_1gb_per_attempt,
        time_min=get_35min_per_attempt,
        tmpdir="tmp",
    group:
        "retrieve_baseline"
    log:
        "logs/grep/{sample}.log",
    params:
        ' -v "./.:0,0:0:0,0,0"',
    shell:
        "grep {params} {input} > {output} 2> {log}"
