rule download_ensembl_variations:
    output:
        "<resources>/{species}.{build}.{release}/{build}.{type}.vcf.gz",
    log:
        "<logs>/download_ensembl_variations/{species}.{build}.{release}.{type}.log",
    cache: "omit-software"
    params:
        url=config.get("url", "ftp://ftp.ensembl.org/pub"),
    wrapper:
        "v9.4.1/bio/reference/ensembl-variation"
