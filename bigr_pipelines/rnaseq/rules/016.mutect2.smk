"""
This rule calls germline variants with GATK Mutect2
"""


rule mutect2_germline:
    input:
        fasta=config["reference"]["genome"],
        fasta_index=config["reference"]["genome_index"],
        fasta_dict=config["reference"]["genome_dict"],
        map="gatk/splitncigarreads/{sample}.bam",
        map_index="gatk/splitncigarreads/{sample}.bam.bai",
        germline=config["reference"]["af_only"],
        germline_tbi=config["reference"]["af_only"],
        intervals=config["reference"]["capturekit_bed"],
    output:
        vcf=temp("mutect2/call/{sample}.vcf.gz"),
        f1r2=temp("mutect2/f1r2/{sample}.tar.gz"),
    threads: config.get("max_threads", 20)
    resources:
        time_min=get_5h_per_attempt,
        mem_mb=get_10gb_per_attempt,
        tmpdir="tmp",
    retries: 1
    params:
        extra=config["gatk"].get(
            "mutect2",
        (
        "--max-reads-per-alignment-start 0 "
                "--disable-read-filter MateOnSameContigOrNoMappedMateReadFilter "
            ),
        ),
    log:
        "logs/gatk/mutect2/call/{sample}.log",
    wrapper:
        "bio/gatk/mutect"
