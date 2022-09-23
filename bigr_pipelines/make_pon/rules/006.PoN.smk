rule create_genomics_db:
    input:
        aln=expand(
            "samtools/filter/{sample}.bam",
            sample=design.index
        ),
        intervals=config[genome_id]["bed"],
        genomicsdb=expand(
            "mutect2/filter/{sample}.vcf.gz",
            sample=design.index
        ),
        ref=config[genome_id]["fasta"],
        ref_idx=fai_file,
        ref_dict=dict_file
    output:
        genomicsdb=directory("data_output/PoN.gdb")
    resources:
        mem_mb=get_4gb_per_attempt,
        time_min=get_4h_per_attempt,
        tmpdir="tmp"
    log:
        "logs/create_genomics_db.log"
    params:
        extra=""
    env:
        "envs/gatk.yaml"
    script:
        str(workflow_source_dir / "006.PoN.GenomicsDBImport.py")
    

rule creat_somatic_pon:
    input:
        ref=config[genome_id]["fasta"],
        ref_idx=fai_file,
        ref_dict=dict_file
        aln=expand(
            "samtools/filter/{sample}.bam",
            sample=design.index
        ),
        intervals=config[genome_id]["bed"],
        gdb=directory("PoN.gdb")
    output:
        pon=protected("data_output/PoN.vcf.gz")
    resources:
        mem_mb=get_20gb_per_attempt,
        time_min=get_8h_per_attempt,
        tmpdir="tmp"
    params:
        extra="--create-output-variant-index true --create-output-variant-md5 true"
    log:
        "logs/gatk/create_pon.log"
    env:
        "envs/gatk.yaml"
    script:
        str(workflow_source_dir / "006.PoN.CreateSomaticPanelOfNormals.py")
