rule bedsimilarity:
    input:
        left="bed/{sample1}.bed",
        right="bed/{sample2}.bed",
    output:
        pipe("similarity/{sample1}_{sample2}.txt")
    threads: 1
    resources:
        mem_mb=lambda wildcards, input: input.size_mb + (1024 * 4),
        time_min=lambda wildcards, input: int(input.size_mb / 1024) * 30,
        tmpdir="tmp"
    params:
        left_overlap=config["jaccard"].get("a_overlaps_b", "1E-9"),
        right_overlap=config["jaccard"].get("b_overlaps_a", "1E-9"),
        minimum_overlap_percent=config["jaccard"].get("minimum_overlap_percent", 0.7),
        extra="-e"
    log:
        "logs/bedtools_jaccard/{sample1}_{sample2}.log"
    conda:
        "../envs/bedtools.yaml"
    shell:
        "bedtools jaccard "
        "{params.extra} "
        "-f {params.left_overlap} "
        "-F {params.right_overlap} "
        "-r {minimum_overlap_percent} "
        "-a {left} -b {right} "
        "> {output} 2> {log}"


rule annotate_jaccard:
    input:
        "similarity/{sample1}_{sample2}.txt",
    output:
        temp("similarity/{sample1}_{sample2}.annotated.txt"),
    threads: 1
    resources:
        mem_mb=128,
        time_min=2,
        tmpdir="tmp",
    log:
        "logs/annotate_jaccard/{sample1}_{sample2}.log"
    script:
        results = [wildcards.sample1, wildcards.sample2]
        with open(snakemake.input) as similarity:
            _ = next(similarity)
            results += next(similarity).split("\t")
        
        with open(snakemake.output) as similarity:
            similarity.write()