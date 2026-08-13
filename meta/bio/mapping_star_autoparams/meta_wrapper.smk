import io
import typing
import yaml


def get_xsv_columns(stranding: str) -> str:
    """provide xsv select arguments to extract read counts from star output"""
    extra: str = "--delimiter $'\t' --no-header"
    if stranding == "unstranded":
        return f"{extra} 1,2"
    if stranding == "foreward":
        return f"{extra} 1,3"
    if stranding == "reverse":
        return f"{extra} 1,4"
    raise ValueError(
        f"Unexepcted {stranding=}, "
        "expected one of {'unstranded', 'foreward', 'reverse'}"
    )


def get_star_params(
    genome_statistics: io.TextIOWrapper,
    chimeric: bool = False,
) -> str:
    """provide star align arguments to align reads with given genome"""
    extra: str = str(
        " --outFilterType 'BySJout'"
        " --outFilterMultimapNmax 20"
        " --alignSJoverhangMin 8"
        " --alignSJDBoverhangMin 1"
        " --outFilterMismatchNmax 999"
        " --outFilterMismatchNoverReadLmax 0.04"
        " --twopassMode 'Basic'"
        " --outSAMtype BAM Unsorted"
    )
    genome_data: list[typing.Any] = yaml.safe_load(genome_statistics)
    try:
        intron_min: int = genome_data[0]["transcript"]["without_isoforms"][
            "Shortest intron into exon part (bp)"
        ]
    except KeyError:
        intron_min = 21
    try:
        intron_max: int = genome_data[0]["transcript"]["without_isoforms"][
            "Longest intron into exon part (bp)"
        ]
        matex_gap_max: int = intron_max + 300
    except KeyError:
        intron_max = 1000000
        matex_gap_max = 1000000

    extra += str(
        f" --alignIntronMin {intron_min}"
        f" --alignIntronMax {intron_max}"
        f" --alignMatesGapMax {matex_gap_max}"
    )

    if chimeric is True:
        raise NotImplementedError("Not chimeric alignment implemented. Sorry.")
    else:
        extra += " --outSAMattributes 'All'"

    print(extra)
    return extra


rule star_index:
    """build genome index for star mapping"""
    input:
        fasta="<genome_sequence>",
        fasta_index="<genome_index>",
        gtf="<genome_annotation>",
    output:
        directory("<resources>/<species>.<build>.<release>/star_index"),
    log:
        "<logs>/star_index/<species>.<build>.<release>.log",
    benchmark:
        "<benchmarks>/star_index/<species>.<build>.<release>.tsv"
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
        idx="<resources>/<species>.<build>.<release>/star_index",
        statistics="<genome_agat_statistics>",
    output:
        aln=temp("<temp>/star_align/<species>.<build>.<release>.{sample}/{sample}.bam"),
        log=temp("<temp>/star_align/<species>.<build>.<release>.{sample}/Log.out"),
        log_progress=temp(
            "<temp>/star_align/<species>.<build>.<release>.{sample}/Log.progress.out"
        ),
        log_final=temp(
            "<temp>/star_align/<species>.<build>.<release>.{sample}/Log.final.out"
        ),
        sj=temp("<temp>/star_align/<species>.<build>.<release>.{sample}/SJ.out.tab"),
    log:
        "<logs>/star_align/{sample}.<species>.<build>.<release>.log",
    benchmark:
        "<benchmarks>/star_align/{sample}.<species>.<build>.<release>.tsv"
    threads: 20
    params:
        extra=parse_input(input.statistics, parser=get_star_params),
    wrapper:
        "v9.4.2/bio/star/align"


rule sambamba_sort_reads:
    """sort star mapped reads to gain disk space and work faster"""
    input:
        "<temp>/star_align/<species>.<build>.<release>.{sample}/{sample}.bam",
    output:
        temp("<temp>/sambamba_sort_reads/<species>.<build>.<release>.{sample}.bam"),
    log:
        "<logs>/sambamba_sort_reads/{sample}.<species>.<build>.<release>.log",
    benchmark:
        "<benchmarks>/sambamba_sort/{sample}.<species>.<build>.<release>.star_reads.tsv"
    threads: 20
    params:
        extra="",
    wrapper:
        "v3.11.0/bio/sambamba/sort"


rule sambamba_index_sorted_reads:
    """access faster to sorted reads chunks"""
    input:
        "<temp>/sambamba_sort_reads/<species>.<build>.<release>.{sample}.bam",
    output:
        temp("<temp>/sambamba_sort_reads/<species>.<build>.<release>.{sample}.bam.bai"),
    log:
        "<logs>/sambamba_index_sorted_reads/{sample}.<species>.<build>.<release>.log",
    benchmark:
        "<benchmarks>/sambamba_index/{sample}.<species>.<build>.<release>.sorted_star_reads.tsv"
    threads: 3
    params:
        extra="",
    wrapper:
        "v6.1.0/bio/sambamba/index"


rule sambamba_mark_duplicates:
    """mark or remove artefactual duplicates in mapping"""
    input:
        "<temp>/sambamba_sort_reads/<species>.<build>.<release>.{sample}.bam",
    output:
        temp("<temp>/sambamba_mark_duplicates/<species>.<build>.<release>.{sample}.bam"),
    log:
        "<logs>/sambamba_mark_duplicates/<species>.<build>.<release>.{sample}.log",
    benchmark:
        "<benchmarks>/sambamba_mark_duplicates/<species>.<build>.<release>.{sample}.tsv"
    threads: 20
    params:
        extra="--remove-duplicates",
    wrapper:
        "v6.1.0/bio/sambamba/markdup"


use rule sambamba_index_sorted_reads as sambamba_index_deduplicated_reads with:
    input:
        "<temp>/sambamba_mark_duplicates/<species>.<build>.<release>.{sample}.bam",
    output:
        temp(
            "<temp>/sambamba_mark_duplicates/<species>.<build>.<release>.{sample}.bam.bai"
        ),
    log:
        "<logs>/sambamba_index_deduplicated_reads/<species>.<build>.<release>.{sample}.log",
    benchmark:
        "<benchmarks>/sambamba_index/<species>.<build>.<release>.{sample}.deduplicated_reads.tsv"


rule samtools_cram_aligned_reads:
    """compress mapped reads to save disk space"""
    input:
        "<temp>/sambamba_mark_duplicates/<species>.<build>.<release>.{sample}.bam",
        "<temp>/sambamba_mark_duplicates/<species>.<build>.<release>.{sample}.bam.bai",
        ref="<genome_sequence>",
        ref_idx="<genome_index>",
    output:
        "<results>/<species>.<build>.<release>/mapping_star/{sample}.cram",
    log:
        "<logs>/samtools_cram_aligned_reads/{sample}.star.<species>.<build>.<release>.log",
    benchmark:
        "<benchmarks>/samtools_view/{sample}.star.<species>.<build>.<release>.cram_deduplicated_reads.tsv"
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
        "<benchmarks>/samtools_index_cram_reads/{sample}.<species>.<build>.<release>.tsv"
    threads: 3
    params:
        extra="",
    wrapper:
        "v9.14.0/bio/samtools/index"
