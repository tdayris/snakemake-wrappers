# QC report aggregation
"""
005.multiqc
from
-> 002.fastqc
-> 002.fastq_screen
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
        # txt=expand(
        #     "fastq_screen/{sample}.fastq_screen.txt",
        #     sample=design["Sample_id"],
        # ),
        # png=expand(
        #     "fastq_screen/{sample}.fastq_screen.png",
        #     sample=design["Sample_id"],
        # ),
    output:
        "data_output/multiqc.html",
        directory("data_output/multiqc_data"),
    threads: 1
    resources:
        mem_mb=get_2gb_and_6gb_per_attempt,
        time_min=get_1h_per_attempt,
        tmpdir="tmp",
    params:
        "--flat",
    log:
        "logs/003.multiqc.log",
    wrapper:
        "bio/multiqc"


# Additional behaviour for demultiplexing automaton
# QC report aggregation
"""
005.multiqc
from
-> 002.fastqc
-> 002.fastq_screen
-> 004.unzip_stats
-> 004.unzip_runparams
-> 004.unzip_runinfo
-> 004.unzip_interop
by
-> End job
"""


rule multiqc_stats:
    input:
        **demux_input(stats, interop, runinfo, runparams)
    output:
        "output/multiqc.html",
        directory("output/multiqc_data"),
    threads: 1
    resources:
        mem_mb=get_2gb_and_6gb_per_attempt,
        time_min=get_1h_per_attempt,
        tmpdir="tmp",
    params:
        "--flat",
    log:
        "logs/003.multiqc.log",
    wrapper:
        "bio/BiGR/multiqc_illumina"
