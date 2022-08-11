"""
mutect 2.0 : java -Xmx8g -jar GenomeAnalysisTK.jar -T MuTect2 -I:tumor {sample_tumor}.cutadapt.sorted.bam -I:normal {sample_normal}.cutadapt.sorted.bam --output_mode EMIT_VARIANTS_ONLY -o {sample_tumor}vs{sample_normal}_mutect2.vcf.gz -R /mnt/beegfs/userdata/m_deloger/B20002_FRFA_01/DBs/bwa/hg38_basechr.fa --max_alt_alleles_in_normal_count 2 --max_alt_allele_in_normal_fraction 0.04 --maxReadsInRegionPerSample 100000
"""


rule mutect2:
    input:
        unpack(get_trio),
    output:
        "gatk/mutect2/{sample}.vcf.gz",
    threads: 1
    resources:
        mem_mb=get_10gb_per_attempt,
        java_mem_gb=get_10gb_per_attempt,
        time_min=get_2h_per_attempt,
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
