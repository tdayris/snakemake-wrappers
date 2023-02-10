### Some tools requires a VEP annotation instead of
### SnpEff. However, SnpEff provides annotation that
### are missing in VEP and vice versa.

rule get_vep_cache_human_99:
    output:
        directory("resources/vep/homo_sapiens.GRCh38.99")
    threads: 1
    resources:
        mem_mb=get_2go_per_attempt,
        time_min=get_5h_per_attempt,
        tmpdir=tmp
    params:
        species="homo_sapiens",
        build="GRCh38",
        release="99",
    log:
        "logs/vep/cache/homo_sapiens.GRCh38.99.log"
    cache: "omit-software"
    wrapper:
        "bio/vep/cache"


rule get_vep_cache_mouse_99:
    output:
        directory("resources/vep/mus_musculus.GRCm38.99")
    threads: 1
    resources:
        mem_mb=get_2go_per_attempt,
        time_min=get_5h_per_attempt,
        tmpdir=tmp
    params:
        species="mus_musculus",
        build="GRCm38",
        release="99",
    log:
        "logs/vep/cache/mus_musculus.GRCm38.99.log"
    cache: "omit-software"
    wrapper:
        "bio/vep/cache"


rule vep_annotate:
    input:
        calls="snpsift/fixed/{sample}.vcf",
        bam="sambamba/markdup/{sample}_tumor.bam",
        bam_index="sambamba/markdup/{sample}_tumor.bam.bai",
        cache="resources/vep/{species}.{build}.{release}",
        fasta=config["reference"]["fasta"],
        fai=config["reference"]["fasta_index"],
        gtf=config["reference"]["gtf"],
        gtf_csi=config["reference"]["gtf_csi"],
    ouptut:
        calls=temp("vep/annotate/{sample}.vcf"),
        stats=temp("vep/stats/{sample}_vep.txt"),
    threads:
        config.get("max_threads", 1)
    resources:
        mem_mb=get_8go_per_attempt,
        time_min=get_3h_per_attempt,
        tmpdir=tmp,
    log:
        "logs/vep/annotate/{sample}.log"
    params:
        extra="--vcf_info_field CSQ --hgvs --hgvsg --force_overwrite --offline --cache --format vcf --vcf --stats_text"
        cache_version="99"
        species="GRCh38" if config["reference"].get("ncbi_build", "GRCh38").lower() == "grch38" else "GRCm38"
        build=config["reference"].get("ncbi_build", "GRCh38")
    conda:
        "../envs/vep.yaml"
    shell:
        "vep {params.extra} "
        "--fork {threads} "
        "--fasta {input.fasta} "
        "--gtf {input.gtf} "
        "--dir_cache {input.cache} "
        "--cache_version {params.cache_version} "
        "--species {params.species} "
        "--assembly {params.build}"
        "--output_file {output.calls} "
        "--stats_file {output.stats} "
        "--input_file {input.calls} "
        "--bam {input.bam} "
        "> {log} 2>&1 "


rule vcf_report:
    input:
        calls="vep/annotate/{sample}.vcf.gz",
        bam="sambamba/markdup/{sample}_tumor.bam",
        fasta=config["reference"]["fasta"],
        fai=config["reference"]["fasta_index"],
    output:
        directory("data_output/HTML/{sample}_report")
    threads: config.get("max_threads", 1)
    resources:
        mem_mb=get_2go_per_attempt,
        time_min=get_1h_per_attempt,
        tmpdir=tmp,
    log:
        "logs/rbt/vcf_report/{sample}.log"
    params:
        extra="--annotation-field CSQ "
    conda:
        "../envs/rbt.yaml"
    shell:
        "rbt vcf-report "
        "{params.extra} "
        "--threads {threads} "
        "{input.fasta} "
        "{output}"