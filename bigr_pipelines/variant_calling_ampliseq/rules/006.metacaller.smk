# Merge calls from Varscan2 and Mutect2.

metacaller_germline_config = {
    "genome": config["ref"]["fasta"], 
    "bed": config["ref"]["capture_kit_bed"]
}

module metacaller_germline_meta:
    snakefile: "../../../meta/bio/meta_caller_germline/test/Snakefile"
    config: metacaller_germline_config


use rule * from metacaller_germline_meta