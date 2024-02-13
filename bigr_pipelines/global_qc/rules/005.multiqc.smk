# QC report aggregation
"""
005.multiqc
from
-> 002.fastqc
-> 002.fastq_screen
-> 009.concat_cellranger_multi_RNA.smk
-> 009.concat_cellranger_atac_ATAC.smk
-> 009.concat_cellranger_arc_RNA_ATAC.smk
by
-> End job
"""
    
def CR_html_input(wildcards):
    if isScRNAData or isScATACData or isScRNAATACData:
        all_samples = sum([list(dic_DATA[i].keys()) for i in list(dic_DATA.keys())], [])
        return ["cellranger/web_summary/" + s + "_web_summary_mqc.html" for s in all_samples ]
    else:
        return []
    
def CR_csv_input(wildcards):
    res=list()
    if isScRNAData :
        res.append("cellranger/csv_concat/CellRanger_RNA_summary_mqc.csv")
    if isScATACData :
        res.append("cellranger/csv_concat/CellRanger_ATAC_summary_mqc.csv")
    if isScRNAATACData :
        res.append("cellranger/csv_concat/CellRanger_RNA_ATAC_summary_mqc.csv")
    return res

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
        CR_csv_res=CR_csv_input,
        #CR_html_res=CR_html_input,
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
        "--flat --force --module fastp --module fastq_screen ",
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
        "--flat --force --module fastp --module fastq_screen",
    log:
        "logs/003.multiqc.log",
    wrapper:
        # "bio/BiGR/multiqc_illumina"
        "bio/multiqc"
