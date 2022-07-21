rule multiqc:
    input:
        html=expand(
            "fastp/html/pe/{sample}_{status}.fastp.html",
            sample=design["Sample_id"],
            status=["normal", "tumor"],
        ),
        json=expand(
            "fastp/json/pe/{sample}_{status}.fastp.json",
            sample=design["Sample_id"],
            status=["normal", "tumor"],
        ),
        sambamba_metrics=expand(
            "sambamba/markdup/{sample}_{status}.bam",
            sample=design["Sample_id"],
            status=["normal", "tumor"],
        ),
        fastq_screen=expand(
            "fastq_screen/{sample}.{stream}.{status}.fastq_screen.{ext}",
            sample=design["Sample_id"],
            stream=["1", "2"],
            ext=["txt", "png"],
            status=["normal", "tumor"],
        ),
        picard_summary=expand(
            "picard/alignment_summary/{sample}_{status}.summary.txt",
            sample=design["Sample_id"],
            status=["normal", "tumor"],
        ),
        snpeff=expand(
            "snpeff_snpsift/snpeff/csvstats/{sample}.genes.txt",
            sample=design["Sample_id"],
        ),
    output:
        report(
            "multiqc/variant_calling_somatic.html",
            caption="../common/reports/multiqc.rst",
            category="Quality Controls",
        ),
    message:
        "Aggregating quality reports from SnpEff"
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: min(attempt * 1536, 10240),
        time_min=lambda wildcards, attempt: attempt * 35,
        tmpdir="tmp",
    log:
        "logs/multiqc.log",
    wrapper:
        "bio/multiqc"
