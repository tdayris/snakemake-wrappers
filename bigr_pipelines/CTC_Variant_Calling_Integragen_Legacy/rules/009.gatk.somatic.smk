"""
mutect 2.0 : java -Xmx8g -jar GenomeAnalysisTK.jar -T MuTect2 -I:tumor {sample_tumor}.cutadapt.sorted.bam -I:normal {sample_normal}.cutadapt.sorted.bam --output_mode EMIT_VARIANTS_ONLY -o {sample_tumor}vs{sample_normal}_mutect2.vcf.gz -R /mnt/beegfs/userdata/m_deloger/B20002_FRFA_01/DBs/bwa/hg38_basechr.fa --max_alt_alleles_in_normal_count 2 --max_alt_allele_in_normal_fraction 0.04 --maxReadsInRegionPerSample 100000
"""


rule mutect2:
    input:
        tumor="sambamba/markdup/{sample}.tumor.bam",
        normal="sambamba/markdup/{sample}.baseline.bam",
        wbc="sambamba/markdup/{sample}.wbc.bam",
        fasta=config["ref"]["fasta"],
    output:
        "gatk/mutect2/{sample}.vcf.gz",
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024 * 10,
        time_min=lambda wildcards, attempt: attempt * 75,
        java_mem_gb=lambda wildcards, attempt: attempt * 9,
        tmpdir="tmp",
    log:
        "logs/gatk/mutect2/{sample}.log",
    params:
        "--max_alt_alleles_in_normal_count 2 "
        "--max_alt_allele_in_normal_fraction 0.04 "
        "--maxReadsInRegionPerSample 100000 "
        "--output_mode EMIT_VARIANTS_ONLY ",
    conda:
        "envs/conda/gatk3.yaml"
    shell:
        "java -Xmx{resources.java_mem_gb}GB "
        "-jar GenomeAnalysisTK.jar "
        "{params} "
        "-T MuTect2 "
        "-I:tumor {input.tumor} "
        "-I:normal {input.normal} "
        "-o {output} "
        "-R {input.fasta} "
        "> {log} 2>&1 "
