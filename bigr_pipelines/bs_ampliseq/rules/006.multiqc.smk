rule multiqc:
    input:
        mapping=expand("bismark/align/{sample}_PE_report.txt", sample=sample_list),
        extract=expand("bismark/meth/{sample}_splitting_report.txt", sample=sample_list),
        nuc=expand("bismark/align/{sample}.nucleotide_stats.txt", sample=sample_list),
        m_bias=expand("bismark/meth/{sample}.M-bias.txt", sample=sample_list),
        fastp=expand(
            "fastp/{ext}/{sample}.fastp.{ext}", sample=sample_list, ext=["html", "json"]
        ),
    output:
        "data_output/MultiQC/Bismark.html",
    params:
        config.get("multiqc", "--force"),
    log:
        "logs/multiqc.log",
    wrapper:
        "bio/multiqc"
