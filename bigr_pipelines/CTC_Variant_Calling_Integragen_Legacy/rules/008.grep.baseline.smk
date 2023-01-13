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
        tmpdir=tmp,
    group:
        "retrieve_baseline"
    log:
        "logs/grep/{sample}.log",
    conda:
        str(workflow_source_dir / "envs" / "bash.yaml")
    params:
        ' -v "./.:0,0:0:0,0,0"',
    shell:
        "grep {params} {input} > {output} 2> {log}"


rule zip_baseline_variants:
    input:
        "gatk/select_variants/baseline/{sample}.vcf"
    output:
        protected("data_output/Baseline/{sample}.vcf.gz")
    threads: 2
    resources:
        mem_mb=get_4gb_per_attempt,
        time_min=get_45min_per_attempt,
        tmpdir="tmp"
    log:
        "logs/bcftools/view/{sample}.baseline.log"
    params:
        extra=""
    wrapper:
        "bio/bcftools/view"


rule tabix_baseline_variants:
    input:
        "gatk/select_variants/Baseline/{sample}.vcf.gz"
    output:
        protected("data_output/Baseline/{sample}.vcf.gz.tbi")
    threads: 1
    resources:
        mem_mb=get_4gb_per_attempt,
        time_min=get_45min_per_attempt,
        tmpdir="tmp"
    log:
        "logs/tabix/{sample}.baseline.log"
    params:
        "-p vcf"
    wrapper:
        "bio/tabix/index"