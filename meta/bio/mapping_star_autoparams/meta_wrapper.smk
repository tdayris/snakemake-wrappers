import yaml
import typing

def get_xsv_columns(stranding: str) -> str:
    """provide xsv select arguments to extract read counts from star output"""
    extra: str = "--delimiter $'\t' --no-header"
    if stranding == "unstranded":
        return f"{extra} 1,2"
    if stranding == "foreward":
        return f"{extra} 1,3"
    if stranding == "reverse"
        return f"{extra} 1,4"
    raise ValueError(
        f"Unexepcted {stranding=}, "
        "expected one of ['unstranded', 'foreward', 'reverse']"
    )


def get_star_params(genome_statistics: str) -> str:
    """provide star align arguments to align reads with given genome"""
    extra: str = str(
        "--outFilterType 'BySJout' "
        "--outFilterMultimapNmax 20 "
        "--alignSJoverhangMin 8 "
        "--alignSJDBoverhangMin 1 "
        "--outFilterMismatchNmax 999 "
        "--outFilterMismatchNoverReadLmax 0.04 "
        "--alignMatesGapMax 1000000 "
        "--outSAMattributes 'All' "
        "--twopassMode 'Basic' "
        "--outSAMtype 'BAM' 'Unsorted' "
    )
    with open(genome_statistics, "r") as yaml_genome_stream:
        genome_data: list[typing.Any] = yaml.safe_load(yaml_genome_stream)
        intron_min: int = genome_data[0]
    #    --alignIntronMin 20 minimum intron length
    # --alignIntronMax 1000000 maximum intron length

    return extra


rule star_index:
    """build genome index for star mapping"""
    input:
        fasta="<genome_sequence>",
        fasta_index="<genome_index>",
    output:
        directory("<reference>/<species>.<build>.<release>/star_index"),
    log:
        "<log>/star_index/<species>.<build>.<release>.log",
    benchmark:
        "<benchmark>/star_index/<species>.<build>.<release>.tsv"
    threads: 20
    params:
        extra="",
        sjdbOverhang=100,
    wrapper:
        "v3.3.7/bio/star/index"
    

rule star_align:
    """align reads over indexed genome with star"""
    input:
        fq1=["<upstream_reads>"],
        fq2=["<downstream_reads>"],
        idx="<reference>/<species>.<build>.<release>/star_index",
        statistics="<genome_agat_statistics>"
    output:
        aln=temp("<tmp>/star_align/<species>.<build>.<release>.{sample}/{sample}.bam"),
        log=temp("<tmp>/star_align/<species>.<build>.<release>.{sample}/Log.out"),
        log_progress=temp("<tmp>/star_align/<species>.<build>.<release>.{sample}/Log.progress.out"),
        log_final=temp("<tmp>/star_align/<species>.<build>.<release>.{sample}/Log.final.out"),
        reads_per_gene=temp("<tmp>/star_align/<species>.<build>.<release>.{sample}/ReadsPerGene.out.tab"),
        sj=temp("<tmp>/star_align/<species>.<build>.<release>.{sample}/SJ.out.tab"),
    log:
        "<log>/star_align/{sample}.<species>.<build>.<release>.log",
    benchmark:
        "<benchmark>/star_align/{sample}.<species>.<build>.<release>.tsv"
    threads: 20
    params:
        extra="",
    wrapper:
        "v9.4.2/bio/star/align"


rule sambamba_sort_reads:
    """sort star mapped reads to gain disk space and work faster"""
    input:
        "<tmp>/star_align/<species>.<build>.<release>.{sample}/{sample}.bam",
    output:
        temp("<tmp>/sambamba_sort_reads/<species>.<build>.<release>.{sample}.bam"),
    log:
        "<log>/sambamba_sort_reads/{sample}.<species>.<build>.<release>.log",
    benchmark:
        "<benchmark>/sambamba_sort_reads/{sample}.<species>.<build>.<release>.tsv",
    threads: 20
    params:
        extra="",
    wrapper:
        "v3.11.0/bio/sambamba/sort"


rule sambamba_index_sorted_reads:
    """access faster to sorted reads chunks"""
    input:
        "<tmp>/sambamba_sort_reads/<species>.<build>.<release>.{sample}.bam",
    output:
        temp("<tmp>/sambamba_sort_reads/<species>.<build>.<release>.{sample}.bam.bai"),
    log:
        "<log>/sambamba_index_sorted_reads/{sample}.<species>.<build>.<release>.log",
    benchmark:
        "<benchmark>/sambamba_index_sorted_reads/{sample}.<species>.<build>.<release>.tsv",
    threads: 3
    params:
        extra="",
    wrapper:
        "v6.1.0/bio/sambamba/index"


rule sambamba_mark_duplicates:
    """mark or remove artefactual duplicates in mapping"""
    input:
        "<tmp>/sambamba_sort_reads/<species>.<build>.<release>.{sample}.bam",
    output:
        temp("<tmp>/sambamba_mark_duplicates/<species>.<build>.<release>.{sample}.bam"),
    log:
        "<log>/sambamba_mark_duplicates/<species>.<build>.<release>.{sample}.log",
    benchmark:
        "<benchmark>/sambamba_mark_duplicates/<species>.<build>.<release>.{sample}.tsv",
    threads: 20
    params:
        extra="--remove-duplicates",
    wrapper:
        "v6.1.0/bio/sambamba/markdup"


rule samtools_cram_aligned_reads:
     """compress mapped reads to save disk space"""
    input:
        "<tmp>/sambamba_mark_duplicates/<species>.<build>.<release>.{sample}.bam",
    output:
        "<results>/<species>.<build>.<release>/mapping_star/{sample}.cram",
    log:
        "<log>/samtools_cram_aligned_reads/{sample}.star.<species>.<build>.<release>.log",
    benchmark:
        "<benchmark>/samtools_cram_aligned_reads/{sample}.star.<species>.<build>.<release>.tsv",
    threads: 8
    params:
        extra="",
    wrapper:
        "v9.14.0/bio/samtools/view"


rule samtools_index_cram_reads:
    """access faster to cramed chunks"""
    input:
        "<results>/<species>.<build>.<release>/mapping_star/{sample}.cram",
    output:
        "<results>/<species>.<build>.<release>/mapping_star/{sample}.cram.crai",
    log:
        "<logs>/samtools_index_cram_reads/{sample}.<species>.<build>.<release>.log",
    benchmark:
        "<benchmark>/samtools_index_cram_reads/{sample}.<species>.<build>.<release>.tsv",
    threads: 3
    params:
        extra="",
    wrapper:
        "v9.14.0/bio/samtools/index"


rule xsv_extract_read_counts:
    """extracts unstranded read counts over genes"""
    input:
        table="<tmp>/star_align/<species>.<build>.<release>/{sample}/ReadsPerGene.out.tab",
    output:
        temp("<tmp>/xsv_extract_read_counts/<species>.<build>.<release>.{sample}.csv"),
    log:
        "<log>/xsv_extract_read_counts/{sample}.{stranding}.<species>.<build>.<release>.log",
    benchmark:
        "<benchmark>/xsv_extract_read_counts/{sample}.{stranding}.<species>.<build>.<release>.tsv",
    threads: 1
    params:
        subcommand="select",
        extra=lambda wildcards: get_xsv_columns(str(wildcards.stranding)),
    wrapper:
        "v3.4.0/utils/xsv"


use rule xsv_extract_read_counts as xsv_add_header with:
    """provide sample name to count table"""
    input:
        "<tmp>/xsv_extract_read_counts/<species>.<build>.<release>.{sample}.csv",
    output:
        temp("<results>/<species>.<build>.<release>/mapping_star/{sample}.{stranding}.read_counts.csv"),
    log:
        "<log>/xsv_add_header/<species>.<build>.<release>.{sample}.log",
    benchmark:
        "<benchmark>/xsv_add_header/<species>.<build>.<release>.{sample}.tsv"
    threads: 1
    params:
        subcommand="cat row",
        extra=lambda wildcards: f'<( echo "gene_name,{wildcards.sample}" )',


use rule xsv_extract_read_counts as xsv_aggregate_counts with:
    """provide a single table for all read counts"""
    input:
        expand(
            "<results>/<species>.<build>.<release>/mapping_star/{sample}.{stranding}.read_counts.csv",
            stranding={"foreward", "reverse", "unstranded"},
            sample=config["samples"],
        ),
    output:
        temp("<results>/<species>.<build>.<release>/mapping_star/aggregated_{stranding}_read_counts.csv"),
    log:
        "<log>/xsv_aggregate_counts/<species>.<build>.<release>.log",
    benchmark:
        "<benchmark>/xsv_aggregate_counts/<species>.<build>.<release>.tsv"
    threads: 1
    params:
        subcommand="cat columns",
        extra="",
    

rule bioconvert_star_csv_counts_to_xls:
    """let coworkers access raw counts with ease"""
    input:
        "<results>/<species>.<build>.<release>/mapping_star/aggregated_{stranding}_read_counts.csv",
    output:
        "<results>/<species>.<build>.<release>/mapping_star/{sample}.{stranding}.read_counts.xls",
    log:
        "<log>/bioconvert_star_csv_counts_to_xls/{stranding}.<species>.<build>.<release>.log",
    benchmark:
        "<benchmark>/bioconvert_star_csv_counts_to_xls/stranding}.<species>.<build>.<release>.tsv",
    threads: 1
    params:
        converter="csv2xls",
        extra="",
    wrapper:
        "v9.15.0/bio/bioconvert"


rule compress_read_counts:
    """save disk space"""
    input:
        "<results>/<species>.<build>.<release>/mapping_star/aggregated_{stranding}_read_counts.csv",
    output:
        "<results>/<species>.<build>.<release>/mapping_star/aggregated_{stranding}_read_counts.csv.gz",
    log:
        "<log>/compress_read_counts/<species>.<build>.<release>.{stranding}.log",
    benchmark:
        "<benchmark>/compress_read_counts/<species>.<build>.<release>.{stranding}.tsv"
    threads: 7
    params:
        extra="--compression-level 9",
    wrapper:
        "v9.15.0/utils/crabz"
