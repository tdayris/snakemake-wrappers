"""
grep -v "./.:0,0:0:0,0,0" {patient}_CTC_all.cutadapt.sorted.rmmarkdup.hc.snps.filtered.PASS.{sample1_normal}.vcf > {patient}_CTC_all.cutadapt.sorted.rmmarkdup.hc.snps.filtered.PASS.{sample1_normal}.grepped.vcf
"""
rule grep_out_homozygote:
    input:
        "gatk/select_variants/baseline/{sample}.tmp.vcf"
    output:
        temp("gatk/select_variants/baseline/{sample}.vcf")
    threads: 1
    resource:
        mem_mb = lambda wildcards, attempt: attempt * 512,
        time_min = lambda wildcards, attempt: attempt * 35,
        tmpdir = "tmp"
    group:
        "retrieve_baseline"
    log:
        "logs/grep/{sample}.log"
    params:
        ' -v "./.:0,0:0:0,0,0"'
    shell:
        "grep {params} {input} > {output} 2> {log}"