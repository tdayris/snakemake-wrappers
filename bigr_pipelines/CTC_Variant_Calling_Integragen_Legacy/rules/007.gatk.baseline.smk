"""
gatk 3.7 selectvariants (to separate each "normal" on independent VCFs) : java -Xmx8g -jar GenomeAnalysisTK.jar -T SelectVariants -R hg38_basechr.fa -V {patient}_CTC_all.cutadapt.sorted.rmmarkdup.hc.snps.filtered.PASS.vcf -o {patient}_CTC_all.cutadapt.sorted.rmmarkdup.hc.snps.filtered.PASS.{sample1_normal}.vcf -sn {sample1_normal}
"""


rule gatk_select_variants_baseline:
    input:
        vcf="grep/baseline_wbc/{sample}.vcf",
        fasta=config["ref"]["fasta"],
    output:
        temp("gatk/select_variants/baseline/{sample}.tmp.vcf"),
    threads: 1
    resources:
        mem_mb=4 * 1024,
        java_mem_gb=4 * 1024,
        time_min=get_6h_per_attempt,
        tmpdir=tmp,
    group:
        "retrieve_baseline"
    log:
        "logs/gatk/select_variants/baseline/{sample}.pass.log",
    params:
        extra="", #"-sn {sample}.baseline",
        jar="/mnt/beegfs/userdata/t_dayris/GATK3.7/devs/GATK/GenomeAnalysisTK.jar",
        tmp=tmp,
    conda:
        str(workflow_source_dir / "envs" / "gatk.yaml")
    shell:
        "SAMPLE_NAME=$(grep -P \"^#CHROM\" {input.vcf} | rev | cut -f 2 | rev); "
        "gatk "
        "-Xmx{resources.java_mem_gb}M "
        "-Djava.io.tmpdir=\"{params.tmp}\" "
        "-T SelectVariants "
        "{params.extra} "
        "-sn ${{SAMPLE_NAME}} "
        "-R {input.fasta} "
        "-V {input.vcf} "
        "-o {output} "
        "> {log} 2>&1"


rule gatk_select_variants_wbc:
    input:
        vcf="grep/baseline_wbc/{sample}.vcf",
        fasta=config["ref"]["fasta"],
    output:
        temp("gatk/select_variants/wbc/{sample}.tmp.vcf"),
    threads: 1
    resources:
        mem_mb=4 * 1024,
        java_mem_gb=4 * 1024,
        time_min=get_6h_per_attempt,
        tmpdir=tmp,
    group:
        "retrieve_wbc"
    log:
        "logs/gatk/select_variants/wbc/{sample}.pass.log",
    params:
        extra="",
        jar="/mnt/beegfs/userdata/t_dayris/GATK3.7/devs/GATK/GenomeAnalysisTK.jar",
        tmp=tmp,
    conda:
        str(workflow_source_dir / "envs" / "gatk.yaml")
    shell:
        "SAMPLE_NAME=$(grep -P \"^#CHROM\" {input.vcf} | rev | cut -f 1 | rev); "
        "gatk "
        "-Xmx{resources.java_mem_gb}M "
        "-Djava.io.tmpdir=\"{params.tmp}\" "
        "-T SelectVariants "
        "{params.extra} "
        "-sn ${{SAMPLE_NAME}} "
        "-R {input.fasta} "
        "-V {input.vcf} "
        "-o {output} "
        "> {log} 2>&1"