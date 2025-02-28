# This rule calls germline variants with GATK Mutect2
"""
016.mutect2_germline
from
-> 010.gatk_split_n_cigar_reads
by
-> 017.learn_read_orientation_model
-> 018.mutect2_filter
"""


rule mutect2_germline:
    input:
        fasta=config["reference"]["genome"],
        fasta_index=config["reference"]["genome_index"],
        fasta_dict=config["reference"]["genome_dict"],
        map="010.gatk/splitncigarreads/{sample}.bam",
        map_index="010.gatk/splitncigarreads/{sample}.bam.bai",
        germline=config["reference"]["af_only"],
        germline_tbi=config["reference"]["af_only"],
        intervals=config["reference"]["capturekit_bed"],
    output:
        vcf=temp("016.mutect2/call/{sample}.vcf.gz"),
        f1r2=temp("016.mutect2/f1r2/{sample}.tar.gz"),
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
        "logs/016.gatk/mutect2/call/{sample}.log",
    wrapper:
        "bio/gatk/mutect"
