rule mixcr_individual_postanalysis:
    input:
        samples=lambda wildcards: expand(
            "mixcr/assembleContigs/clns/{sample}.clns",
            sample=samples_per_prefixes[wildcards.comparison]
        ),
        metadata="mixcr/post_analysis/{comparison}/metadata.tsv"
    output:
        json=temp("mixcr/post_analysis/{comparison}/individual_post_analysis.json.gz"),
        tsv=temp("mixcr/post_analysis/{comparison}/individual_post_analysis.tsv"),
    threads: 1
    resources:
        mem_mb=get_4gb_per_attempt,
        time_min=get_3h_per_attempt,
        tmpdir="tmp",
    log:
        "logs/mixcr/postanalysis/{sample}.individual.log"
    params:
        extra=config["mixcr"].get("individual_postanalysis", "--default-downsampling read --default-weight-function read"),
        subcommand="individual"
    wrapper:
        "bio/mixcr/postanalysis"


rule mixcr_overlap_postanalysis:
    input:
        samples=lambda wildcards: expand(
            "mixcr/assembleContigs/clns/{sample}.clns",
            sample=samples_per_prefixes[wildcards.comparison]
        ),
        metadata="mixcr/post_analysis/{comparison}/metadata.tsv"
    output:
        json=temp("mixcr/post_analysis/{comparison}/overlap_post_analysis.json.gz"),
        tsv=temp("mixcr/post_analysis/{comparison}/overlap_post_analysis.tsv"),
    threads: 1
    resources:
        mem_mb=get_4gb_per_attempt,
        time_min=get_3h_per_attempt,
        tmpdir="tmp",
    log:
        "logs/mixcr/postanalysis/{sample}.overlap.log"
    params:
        extra=config["mixcr"].get("overlap_postanalysis", "--default-downsampling read --default-weight-function read"),
        subcommand="overlap"
    wrapper:
        "bio/mixcr/postanalysis"