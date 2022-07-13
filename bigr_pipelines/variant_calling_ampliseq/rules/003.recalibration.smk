gatk_bqsr_meta_config = {
    "threads": config["threads"], 
    "genome": config["ref"]["fasta"], 
    "dbsnp": config["ref"]["dbsnp"]
}


module gatk_bqsr_meta:
    snakefile: "../../../meta/bio/gatk_bqsr/test/Snakefile"
    config: gatk_bqsr_meta_config


use rule gatk_apply_baserecalibrator from gatk_bqsr_meta with:
    input:
        bam="samtools/filter/{sample}.bam",
        bam_index=get_bai("samtools/filter/{sample}.bam"),
        ref=config['ref']['fasta'],
        ref_idx=get_fai(config['ref']['fasta']),
        ref_dict=get_dict(config['ref']['fasta']),
        recal_table="gatk/recal_data_table/{sample}.grp"


use rule gatk_compute_baserecalibration_table from gatk_bqsr_meta with:
    input:
        bam="samtools/filter/{sample}.bam",
        bam_index=get_bai("samtools/filter/{sample}.bam"),
        ref=config['ref']['fasta'],
        ref_idx=get_fai(config['ref']['fasta']),
        ref_dict=get_dict(config['ref']['fasta']),
        known=config['ref']['dbsnp'],
        known_idx=get_tbi(config['ref']['dbsnp'])
