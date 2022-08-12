"""
gatk 3.7 selectvariants (to separate each "normal" on independent VCFs) : java -Xmx8g -jar GenomeAnalysisTK.jar -T SelectVariants -R hg38_basechr.fa -V {patient}_CTC_all.cutadapt.sorted.rmmarkdup.hc.snps.filtered.PASS.vcf -o {patient}_CTC_all.cutadapt.sorted.rmmarkdup.hc.snps.filtered.PASS.{sample1_normal}.vcf -sn {sample1_normal}
"""


rule gatk_select_variants_wbc:
    input:
        vcf="grep/baseline_wbc/{sample}.vcf",
        fasta=config["ref"]["fasta"],
    output:
        temp("gatk/select_variants/baseline/{sample}.tmp.vcf"),
    threads: 1
    resources:
        mem_mb=get_10gb_per_attempt,
        time_min=get_45min_per_attempt,
        java_mem_gb=get_10gb_per_attempt,
        tmpdir="tmp",
    group:
        "retrieve_baseline"
    log:
        "logs/gatk/select_variants/baseline_wbc/{sample}.pass.log",
    params:
        extra="-sn {sample}.baseline",
        jar="/mnt/beegfs/userdata/t_dayris/GATK3.7/devs/GATK/GenomeAnalysisTK.jar"
    conda:
        str(workflow_source_dir / "envs" / "gatk.yaml")
    shell:
        "java -Xmx{resources.java_mem_gb}GB "
        "-jar {params.jar} "
        "-T SelectVariants "
        "{params.extra} "
        "-R {input.fasta} "
        "-V {input.vcf} "
        "-o {output} "
        "> {log} 2>&1"
