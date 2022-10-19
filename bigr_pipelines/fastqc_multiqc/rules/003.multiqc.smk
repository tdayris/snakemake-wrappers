# QC report aggregation
"""
003.multiqc
from
-> 001.fastqc
-> 001.fastq_screen
by
-> End job
"""
rule multiqc:
    input:
        fqc_zip=expand(
            "fastqc/{sample}_fastqc.zip",
            sample=design["Sample_id"],
        ),
        fqc_html=expand(
            "fastqc/{sample}.html",
            sample=design["Sample_id"],
        ),
        txt=expand(
            "fastq_screen/{sample}.fastq_screen.txt",
            sample=design["Sample_id"],
        ),
        png=expand(
            "fastq_screen/{sample}.fastq_screen.png",
            sample=design["Sample_id"],
        )
    output:
        "multiqc/multiqc.html",
        directory("multiqc/multiqc_data")
    threads: 1
    resources:
        mem_mb=get_2gb_and_6gb_per_attempt,
        time_min=get_1h_per_attempt,
        tmpdir="tmp"
    params:
        "--flat"
    log:
        "logs/003.multiqc.log"
    wrapper:
        "bio/multiqc"



# Additional behaviour for demultiplexing automaton
# QC report aggregation
"""
003.multiqc
from
-> 001.fastqc
-> 001.fastq_screen
-> 003.unzip_stats
by
-> End job
"""
use rule multiqc as irods_complient with:
    input:
        fqc_zip=expand(
            "fastqc/{sample}_fastqc.zip",
            sample=design["Sample_id"],
        ),
        fqc_html=expand(
            "fastqc/{sample}.html",
            sample=design["Sample_id"],
        ),
        txt=expand(
            "fastq_screen/{sample}.fastq_screen.txt",
            sample=design["Sample_id"],
        ),
        png=expand(
            "fastq_screen/{sample}.fastq_screen.png",
            sample=design["Sample_id"],
        ),
        bcl_json="Stats.json"
    output:
        "output/multiqc.html",
        directory("output/multiqc_data")
    group:
        "stats_inclusion"


# Unzip Stats.json for multiqc inclusion
"""
003.unzip_stats
from
-> Entry job
by
-> 003.multiqc
"""
rule unzip_stats:
    output:
        temp("Stats.json")
    threads: 1
    resources:
        mem_mb=get_768mb_per_attempt,
        time_min=get_15min_per_attempt,
        tmpdir="tmp"
    log:
        "logs/003.unzipping.log"
    group:
        "stats_inclusion"
    shell:
        'unzip -n -d "${{PWD}}" '
        'input/*/archive/*/unaligned/Stats/Stats.json.zip '
        '> {log} 2>&1'