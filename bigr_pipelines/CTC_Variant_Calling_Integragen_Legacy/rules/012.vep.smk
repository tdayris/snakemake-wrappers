"""
ensembl VEP 87.0 refseq (on each normal VCF/TSV) : singularity run -B /mnt/beegfs:/mnt/beegfs /mnt/beegfs/software/vep/87/ensembl-vep_release_87.0.sif vep -i {patient}_CTC_all.cutadapt.sorted.rmmarkdup.hc.snps.filtered.PASS.{sample1_normal}.vcf -o {patient}_CTC_all.cutadapt.sorted.rmmarkdup.hc.snps.filtered.PASS.{sample1_normal}.annotated.vcf --species homo_sapiens --assembly GRCh38 --cache --dir_cache /mnt/beegfs/database/bioinfo/vep/87/ --everything --offline --fasta /mnt/beegfs/database/bioinfo/vep/87/homo_sapiens/87_GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa --vcf --refseq
"""


rule ensembl_vep:
    input:
        vcf="gatk/mutect2/{sample}.vcf",
        cache=config["ref"]["vep"],
        fasta=config["ref"]["fasta"],
    output:
        vcf=temp("vep/annotate/{sample}.vcf"),
    threads: 1
    resources:
        mem_mb=get_10gb_per_attempt,
        time_min=get_45min_per_attempt,
        tmpdir="tmp",
    log:
        "logs/vep/{sample}.log",
    params:
        "--species homo_sapiens "
        "--assembly GRCh38 "
        "--cache "
        "--everything "
        "--offline "
        "--vcf "
        "--refseq ",
    container:
        "/mnt/beegfs/software/vep/87/ensembl-vep_release_87.0.sif"
    shell:
        "vep {params} "
        "-i {input.vcf} "
        "-o {output.vcf} "
        "--dir_cache {input.cache} "
        "--fasta {input.fasta} "
        "> {log} 2>&1 "


rule ensemblvep_hc:
    input:
        cancer_genes=config.get("cancer_genes", "Cancer.genes.cleaned.txt"),
        vcfs=["vep/annotate/{sample}.vcf"],
    output:
        vcf=temp("vep/hc/{sample}.vcf"),
    threads: 1
    resources:
        mem_mb=get_10gb_per_attempt,
        time_min=get_45min_per_attempt,
        tmpdir="tmp",
    log:
        "logs/vep_hc.log",
    params:
        organism=config.get("vep_db", "hg38"),
    container:
        "/mnt/beegfs/software/vep/87/ensembl-vep_release_87.0.sif"
    script:
        "scripts/ensmblVEP_hc.R"


rule ensemblvep_bcr:
    input:
        cancer_genes=config.get("cancer_genes", "Cancer.genes.cleaned.txt"),
        vcfs=["vep/hc/{sample}.vcf"],
    output:
        vcf=temp("vep/mutect/{sample}.vcf"),
    threads: 1
    resources:
        mem_mb=get_10gb_per_attempt,
        time_min=get_45min_per_attempt,
        tmpdir="tmp",
    log:
        "logs/vep_hc.log",
    params:
        organism=config.get("vep_db", "hg38"),
    container:
        "/mnt/beegfs/software/vep/87/ensembl-vep_release_87.0.sif"
    script:
        "scripts/ensmblVEP_mutect.R"


rule ensemblvep_bcr:
    input:
        cancer_genes=config.get("cancer_genes", "Cancer.genes.cleaned.txt"),
        vcfs=["vep/mutect/{sample}.vcf"],
    output:
        vcf=temp("vep/bcr/{sample}.vcf"),
    threads: 1
    resources:
        mem_mb=get_10gb_per_attempt,
        time_min=get_45min_per_attempt,
        tmpdir="tmp",
    log:
        "logs/vep_hc.log",
    params:
        organism=config.get("vep_db", "hg38"),
    container:
        "/mnt/beegfs/software/vep/87/ensembl-vep_release_87.0.sif"
    script:
        "scripts/ensmblVEP_bcr.R"
