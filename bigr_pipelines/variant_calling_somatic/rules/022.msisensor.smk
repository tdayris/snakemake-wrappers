rule msisensor_pro_msi:
    input:
        microsat=config["reference"]["msi_scan"],
        normal="sambamba/markdup/{sample}_normal.bam",
        normal_idx="sambamba/markdup/{sample}_normal.bam.bai",
        tumor="sambamba/markdup/{sample}_tumor.bam",
        tumor_idx="sambamba/markdup/{sample}_tumor.bam.bai",
        bed=config["reference"]["bed"]
    output:
        temp(
            multiext(
                "msisensor/{sample}/{sample}",
                ".msi",
                ".msi_dis",
                ".msi_germline",
                ".msi_somatic"
            )
        )
    threads: config.get("max_threads", 20)
    resources:
        mem_mb=get_6gb_per_attempt,
        time_min=get_45min_per_attempt,
        tmpdir="tmp"
    log:
        "logs/msisensor_pro/msi/{sample}.log"
    params:
        extra = config["msisensor"].get("extra", "")
    wrapper:
        "bio/msisensor_pro/msi"
