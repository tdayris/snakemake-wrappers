# Format term / gene relation from TSV
"""
029.term2gene_TERMS
from
-> Entry Job
by
-> 029.enricher_TERMS
"""


rule term2gene_TERMS_PPI:
    input:
        lambda wildcards: config["clusterprofiler"]["ppi"][wildcards.database],
    output:
        temp("026.clusterprofiler/databases/{database}.ENSEMBLPROT.term2gene.tsv"),
    threads: 1
    resources:
        mem_mb=get_1gb_per_attempt,
        time_min=get_10min_per_attempt,
        tmpdir="tmp",
    params:
        begin='FS=OFS="\t"',
        body=["print $3 FS $1"],
    group:
        "prepare_terms_ppi"
    log:
        "logs/029.awk/prepare_terms2gene/{database}.log",
    wrapper:
        "bio/awk"


# Expand terms in human readable format for TSV files
"""
029.terms2name_TERMS
from
-> Entry Job
by
-> 029.enricher_TERMS
"""


rule terms2name_TERMS_PPI:
    input:
        lambda wildcards: config["clusterprofiler"]["ppi"][wildcards.database],
    output:
        temp("026.clusterprofiler/databases/{database}.ENSEMBLPROT.term2name.tsv"),
    threads: 1
    resources:
        mem_mb=get_1gb_per_attempt,
        time_min=get_10min_per_attempt,
        tmpdir="tmp",
    params:
        begin='FS=OFS="\t"',
        body=["print $3 FS $4"],
    group:
        "prepare_terms_ppi"
    log:
        "logs/029.awk/prepare_terms2name/{database}.log",
    wrapper:
        "bio/awk"


"""
027.term2name_GMT
from
-> 027.term2name_GMT
-> 026.expand_rank_list
by
-> End job
"""


use rule enricher_GMT as enricher_TSV with:
    input:
        gene="026.clusterprofiler/gene_lists/{keytype}/{comparison}.tsv",
        term2gene="026.clusterprofiler/databases/{database}.{keytype}.term2gene.tsv",
        term2name="026.clusterprofiler/databases/{database}.{keytype}.term2name.tsv",
