"""
Compute Itegrated Genome Size from the capture kit bed
"""
rule estimate_igs:
    input:
        bed = config["reference"]["capture_kit_bed"]
    output:
        yaml = temp("igs.yaml")
    threads: 1
    resources:
        mem_mb=get_5gb_per_attempt,
        time_min=get_15min_per_attempt,
        tmpdir="tmp",
    retries: 1
    log:
        "logs/igs_estimation.log"
    wrapper:
        "bio/tmb/igs_estimation"



"""
Compute statistics over somatic mutations' coverage
"""
rule extract_somatic_mutations:
    input:
        vcf="data_output/VCF/{sample}.vcf.gz",
        vcf_tbi="data_output/VCF/{sample}.vcf.gz.tbi",
    output:
        yaml = temp("tmb/{sample}.yaml")
    threads: 1
    resources:
        mem_mb=get_1gb_per_attempt,
        time_min=get_15min_per_attempt,
        tmpdir="tmp",
    retries: 1
    log:
        "logs/extract_somatic_mutations/{sample}.log"
    params:
        filter_in = config.get("filter_in", []),
        filter_out = config.get("filter_out", []),
        min_coverage = config["tmb"].get("min_coverage", 10),
        allele_depth = config["tmb"].get("allele_depth_keyname", "AD")
    wrapper:
        "bio/tmb/extract_somatic"


"""
Compute Tumor Molecular Burden from previous indexes and stats
"""
rule compute_tmb:
    input:
        igs = temp("igs.yaml"),
        samples = expand("tmb/{sample}.yaml", sample=sample_list)
    output:
        tsv = protected("data_output/TMB.tsv")
    threads: 1
    resources:
        mem_mb=get_2gb_per_attempt,
        time_min=get_15min_per_attempt,
        tmpdir="tmp",
    retries: 2
    log:
        "logs/tmb.log"
    params:
        high_threshold = config["tmb"].get("tmb_highness_threshold", 20)
    wrapper:
        "bio/tmb/compute_tmb"