varscan2_germline_meta_config = {
    "genome": config["ref"]["fasta"], 
    "bed": config["ref"]["capture_kit_bed"]
}


module varscan2_germline_meta:
    snakefile: "../../../meta/bio/varscan2_germline/test/Snakefile"
    config: varscan2_germline_meta_config


use rule * from varscan2_germline_meta as *


use rule samtools_mpilup from varscan2_germline_meta with:
    input:
        bam="gatk/recal_bam/{sample}.bam",
        reference_genome=config['ref']['fasta'],
        reference_genome_idx=get_fai(config['ref']['fasta']),
