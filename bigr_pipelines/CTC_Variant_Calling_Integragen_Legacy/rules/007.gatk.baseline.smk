"""
gatk 3.7 selectvariants (to separate each "normal" on independent VCFs) : java -Xmx8g -jar GenomeAnalysisTK.jar -T SelectVariants -R hg38_basechr.fa -V {patient}_CTC_all.cutadapt.sorted.rmmarkdup.hc.snps.filtered.PASS.vcf -o {patient}_CTC_all.cutadapt.sorted.rmmarkdup.hc.snps.filtered.PASS.{sample1_normal}.vcf -sn {sample1_normal}
"""
rule gatk_select_variants:
    input:
        vcf = "grep/baseline_wbc/{sample}.vcf",
        fasta = config["ref"]["fasta"]
    output:
        temp("gatk/select_variants/baseline/{sample}.tmp.vcf")
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 1024 * 10,
        time_min = lambda wildcards, attempt: attempt * 75,
        java_mem_gb = lambda wildcards, attempt: attempt * 9,
        tmpdir = "tmp"
    group:
        "retrieve_baseline"
    log:
        "logs/gatk/select_variants/baseline_wbc/{sample}.pass.log"
    params:
        "-sn {sample}.baseline"
    conda:
        "envs/conda/gatk3.yaml"
    shell:
        "java -Xmx{resources.java_mem_gb}GB "
        "-jar GenomeAnalysisTK.jar "
        "-T SelectVariants "
        "{params} "
        "-R {input.fasta} "
        "-V {input.vcf} "
        "-o {output} "
        "> {log} 2>&1"
    