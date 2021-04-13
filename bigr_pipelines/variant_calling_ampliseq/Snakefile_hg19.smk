import sys

sys.path.append("../common/python")

from file_manager import read_design
from files_linker import link_fq

from snakemake.utils import min_version
min_version("6.0")

container: "docker://continuumio/miniconda3:4.4.10"

default_config_variant_calling_ampliseq = {
    "design" = read_design("design.tsv"),
    "threads": 20,
    "ref": {
        "fasta": "/mnt/beegfs/database/bioinfo/Index_DB/Fasta/Gencode/GRCH37/release_19/DNA/GRCh37.p13.genome.fa",
        "dbsnp": "/mnt/beegfs/database/bioinfo/Index_DB/dbSNP/hg19/144_20150605/All_20150605.vcf.gz",
        "gwascat": "/mnt/beegfs/database/bioinfo/Index_DB/GWASCatalog/gwas_catalog_v1.0.2-studies_r2020-05-03.tsv",
        "snpeff": "/mnt/beegfs/database/bioinfo/Index_DB/SnpEff/GRCh37.75/",
        "gmt": "/mnt/beegfs/database/bioinfo/Index_DB/MSigDB/GeneSets.gmt",
        "kaviar": "/mnt/beegfs/database/bioinfo/Index_DB/Kaviar/HG19/Kaviar-160204-Public/vcfs/Kaviar-160204-Public-hg19-trim.vcf.gz",
        "dbsnp": "/mnt/beegfs/database/bioinfo/Index_DB/dbSNP/hg19/144_20150605/All_20150605.vcf.gz",
        "cosmic": "/mnt/beegfs/database/bioinfo/COSMIC/73_20150629/CosmicCodingMuts.vcf"
    },
    "fastp": {
        "fastp_extra": ""
    }
}

fastq_links = link_fq(
    condig["design"].Sample_id,
    condig["design"].Upstream_fastq,
    condig["design"].Downstream_fastq
)

try:
    if config == dict():
        config = default_config_variant_calling_ampliseq
except NameError:
    config = default_config_variant_calling_ampliseq


module bwa_fixmate:
    snakefile: "../../meta/bio/bwa_fixmate"
    configfile: {
        "threads": config["threads"],
        "genome": config["ref"]["fasta"]
    }


module gatk_bqsr:
    snakefile: "../../meta/bio/gatk_bqsr"
    configfile: {
        "threads": config["threads"],
        "genome": config["ref"]["fasta"],
        "dbsnp": config["ref"]["dbsnp"]
    }


module varscan2_calling:
    snakefile: "../../meta/bio/varscan2_calling"
    configfile: {
        "genome": config["ref"]["fasta"]
    }


rule all:
    input:
        expand(
            "snpeff/gwascat/{sample}.vcf.gz",
            sample=config["design"]["Sample_id"]
        ),
        expand(
            "snpeff/gwascat/{sample}.vcf.gz.tbi",
            sample=config["design"]["Sample_id"]
        ),

#################################
### FINAL VCF FILE INDEXATION ###
#################################

rule tabix_index:
    input:
        "snpeff/gwascat/{sample}.vcf.gz"
    output:
        "snpeff/gwascat/{sample}.vcf.gz.tbi"
    message:
        "Indexing {wildcards.sample} final annotated VCF with tabix"
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1020,
        time_min=lambda wildcards, attempt: attempt * 45
    log:
        "logs/pbgzip/post_gwascat/{sample}.log"
    wrapper:
        "/bio/tabix"


rule compress_pbgzip:
    input:
        "snpeff/gwascat/{sample}.vcf"
    output:
        "snpeff/gwascat/{sample}.vcf.gz"
    message:
        "Compressing {wildcards.sample} final annotated VCF with pbgzip"
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1020,
        time_min=lambda wildcards, attempt: attempt * 45
    log:
        "logs/pbgzip/post_gwascat/{sample}.log"
    wrapper:
        "/bio/compress/pbgzip"

###########################
### VCF FILE ANNOTATION ###
###########################

rule snpsift_gwascat:
    input:
        call = "snpeff/cosmic/{sample}.vcf",
        gwascat = config["ref"]["gwascat"]
    output:
        call = temp("snpeff/gwascat/{sample}.vcf")
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1020 + 4096,
        time_min=lambda wildcards, attempt: attempt * 45
    wrapper:
        "/bio/snpsift/gwascat"


rule snpsift_cosmic:
    input:
        call="snpeff/dbsnp/{sample}.vcf",
        database=config["ref"]["cosmic"]
    output:
        call=temp("snpeff/cosmic/{sample}.vcf")
    log:
        "logs/snpsift/cosmic/{sample}.log"
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1020 + 4096,
        time_min=lambda wildcards, attempt: attempt * 45
    wrapper:
        "/bio/snpsift/annotate"


rule snpsift_dbsnp:
    input:
        call="snpeff/kaviar/{sample}.vcf",
        database=config["ref"]["dbsnp"]
    output:
        call=temp("snpeff/dbsnp/{sample}.vcf")
    log:
        "logs/snpsift/dbsnp/{sample}.log"
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1020 + 4096,
        time_min=lambda wildcards, attempt: attempt * 45
    wrapper:
        "/bio/snpsift/annotate"


rule snpsift_kaviar:
    input:
        call="snpsift/gmt/{sample}.vcf",
        database=config["ref"]["kaviar"]
    output:
        call=temp("snpeff/kaviar/{sample}.vcf")
    log:
        "logs/snpsift/kaviar/{sample}.log"
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1020 + 4096,
        time_min=lambda wildcards, attempt: attempt * 45
    wrapper:
        "/bio/snpsift/annotate"


rule snpsift_gmt:
    input:
        call = "snpeff/{sample}.vcf",
        gmt = config["ref"]["gmt"]
    output:
        call = temp("snpsift/gmt/{sample}.vcf")
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1020 + 4096,
        time_min=lambda wildcards, attempt: attempt * 45
    wrapper:
        "/bio/snpsift/genesets"


rule snpeff_annotate:
    input:
        calls="snpsift/vartype/{sample}.vcf",
        db=config["ref"]["snpeff"]
    output:
        calls=temp("snpeff/{sample}.vcf"),
        stats="snpeff/{sample}.html",
        csvstats="snpeff/{sample}.csv"
    log:
        "logs/snpeff/{sample}.log"
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1020 + 4096,
        time_min=lambda wildcards, attempt: attempt * 45
    wrapper:
        "/bio/snpeff/annotate"


rule snpsift_vartype:
    input:
        vcf="bcftools/{sample}.vcf.gz",
        vcf_tbi="bcftools/{sample}.vcf.gz.tbi"
    output:
        vcf=temp("snpsift/vartype/{sample}.vcf")
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1020 + 4096,
        time_min=lambda wildcards, attempt: attempt * 45
    log:
        "logs/snpsift/varType/{sample}.log"
    wrapper:
        "/bio/snpsift/varType"


###############################
### VARIANT CALLING VARSCAN ###
###############################

use rule * from varscan2_calling as varscan2_calling_*

use rule varscan2_calling_samtools_mpilup with:
    input:
        bam="gatk/recal_bam/{sample}.bam",
        reference_genome=config['ref']['fasta'],
        reference_genome_idx=get_fasta_index_from_genome_path(config['ref']['fasta']),


##############################
### GATK BAM RECALIBRATION ###
##############################

use rule gatk_apply_baserecalibrator from gatk_bqsr as gatk_bqsr_apply_baserecalibrator with:
    input:
        bam="samtools/sort/{sample}.bam",
        bam_index="samtools/sort/{sample}.bam.bai",
        ref=config['ref']['fasta'],
        ref_idx=get_fasta_index_from_genome_path(config['ref']['fasta']),
        ref_dict=get_fasta_dict_from_genome_path(config['ref']['fasta']),
        recal_table="gatk/recal_data_table/{sample}.grp"


use rule gatk_compute_baserecalibration_table from gatk_bqsr as gatk_bqsr_compute_baserecalibration_table_bigr_hg19 with:
    input:
        bam="samtools/sort/{sample}.bam",
        bam_index="samtools/sort/{sample}.bam.bai",
        ref=config['ref']['fasta'],
        ref_idx=get_fasta_index_from_genome_path(config['ref']['fasta']),
        ref_dict=get_fasta_dict_from_genome_path(config['ref']['fasta']),
        known=config['ref']['dbsnp'],
        known_idx=get_vcf_tbi_from_db_path(config['ref']['dbsnp'])


###################
### BWA MAPPING ###
###################

use rule * from bwa_fixmate as bwa_fixmate_*

use rule bwa_fixmate_bwa_mem as bwa_fixmate_bwa_mem_bigr with:
    input:
        reads = expand(
            "fastp/trimmed/pe/{sample}.{stream}.fastq",
            "fastp/trimmed/pe/{sample}.{stream}.fastq",
            stream=["1", "2"],
            allow_missing=True
        )


############################
### FASTP FASTQ CLEANING ###
############################

rule fastp_clean:
    input:
        sample=expand(
            "reads/{sample}.{stream}.fq.gz",
            stream=["1", "2"],
            allow_missing=True
        ),
    output:
        trimmed=expand(
            "fastp/trimmed/pe/{sample}.{stream}.fastq",
            "fastp/trimmed/pe/{sample}.{stream}.fastq",
            stream=["1", "2"],
            allow_missing=True
        ),
        html="fastp/html/pe/{sample}.fastp.html",
        json="fastp/json/pe/{sample}.fastp.json"
    message: "Cleaning {wildcards.sample} with Fastp"
    threads: 1
    resources:
        mem_mb=lambda wildcard, attempt: min(attempt * 4096, 15360),
        time_min=lambda wildcard, attempt: attempt * 45
    params:
        adapters=config.get("fastp_adapters", None),
        extra=config.get("fastp_extra", "")
    log:
        "logs/fastp/{sample}.log"
    wrapper:
        "/bio/fastp"


#################################################
### Gather files from iRODS or mounting point ###
#################################################

rule bigr_copy:
    output:
        "reads/{sample}.{stream}.fq.gz"
    message:
        "Gathering {wildcards.sample} fastq file ({wildcards.stream})"
    threads: 1
    resources:
        mem_mb=lambda wildcard, attempt: min(attempt * 1024, 2048),
        time_min=lambda wildcard, attempt: attempt * 45
    params:
        input=lambda wildcards, output: fastq_links[output]
    log:
        "logs/bigr_copy/{sample}.{stream}.log"
    wrapper:
        "bio/bigr/copy"
