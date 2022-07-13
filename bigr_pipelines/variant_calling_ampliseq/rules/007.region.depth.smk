# Compute mean region size according to a given bin

coverage_variant_region_config = {
    "bed": config["ref"]["capture_kit_bed"],
    "threads": config.get("threads", 20),
    "bin_size": config.get(
        "deeptools", {"bin_size": "10"}
    ).get("bin_size", "10")
}


module coverage_variant_region_meta:
    snakefile: "../../../meta/bio/coverage_variant_region/test/Snakefile"
    config: coverage_variant_region_config


use rule * from coverage_variant_region_meta