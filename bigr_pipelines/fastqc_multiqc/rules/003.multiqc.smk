
rule multiqc:
    input:
        fqc_zip=expand(
            "fastqc/{sample}_{stream}_fastqc.zip",
            sample=design["Sample_id"],
            stream=["1", "2"]
        ),
        fqc_html=expand(
            "fastqc/{sample}.{stream}.html",
            sample=design["Sample_id"],
            stream=["1", "2"]
        ),
        txt=expand(
            "fastq_screen/{sample}.{stream}.fastq_screen.txt",
            sample=design["Sample_id"],
            stream=["1", "2"]
        ),
        png=expand(
            "fastq_screen/{sample}.{stream}.fastq_screen.png",
            sample=design["Sample_id"],
            stream=["1", "2"]
        )
    output:
        "multiqc/multiqc.html",
        directory("multiqc/multiqc_data")
    message:
        "Gathering all quality reports in {output}"
    threads: 1
    resources:
        mem_mb=lambda wildcard, attempt: attempt * 2048,
        time_min=lambda wildcard, attempt: attempt * 50,
        tmpdir="tmp"
    params:
        "--flat"
    log:
        "logs/multiqc.log"
    wrapper:
        str(worflow_source_dir /  "bio"/ "multiqc")



# Additional behaviour for demultiplexing automaton
use rule multiqc as irods_complient with:
    input:
        fqc_zip=expand(
            "fastqc/{sample}_{stream}_fastqc.zip",
            sample=design["Sample_id"],
            stream=["1", "2"]
        ),
        fqc_html=expand(
            "fastqc/{sample}.{stream}.html",
            sample=design["Sample_id"],
            stream=["1", "2"]
        ),
        txt=expand(
            "fastq_screen/{sample}.{stream}.fastq_screen.txt",
            sample=design["Sample_id"],
            stream=["1", "2"]
        ),
        png=expand(
            "fastq_screen/{sample}.{stream}.fastq_screen.png",
            sample=design["Sample_id"],
            stream=["1", "2"]
        ),
        bcl_json="Stats.json"
    output:
        "output/multiqc.html",
        directory("output/multiqc_data")
    group:
        "stats_inclusion"


rule unzip_stats:
    output:
        temp("Stats.json")
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024,
        time_min=lambda wildcards, attempt: attempt * 15,
        tmpdir="tmp"
    log:
        "logs/unzipping.log"
    group:
        "stats_inclusion"
    shell:
        'unzip -n -d "${{PWD}}" '
        'input/*/archive/*/unaligned/Stats/Stats.json.zip '
        '> {log} 2>&1'