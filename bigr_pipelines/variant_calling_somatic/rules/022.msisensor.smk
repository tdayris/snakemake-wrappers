msi_sensor_config = {
    "fasta": config["ref"]["fasta"],
    "bed": config["ref"]["capture_kit_bed"],
    "msi_scan_extra": config["msisensor_pro"].get("scan", ""),
    "msi_pro_extra": config["msisensor_pro"].get("msi", ""),
}


module missensor_pro_meta:
    snakefile:
        str(
            workflow_source_dir
            / ".."
            / ".."
            / "meta"
            / "bio"
            / "msi_sensor_pro"
            / "test"
            / "Snakefile"
        )
    config:
        msi_sensor_config


use rule * from missensor_pro_meta
