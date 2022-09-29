# Perform dot plots on DGE results

"""
030.dotplot
from
-> 029.enricher_TERMS
-> 028.enrichDO
-> 028.enrichDGN
-> 028.enrichDO
-> 027.enrich_GMT
by
-> End job
"""


rule dotplot:
    input:
        rds="027.enrich/{database}/{comparison}/enrich.{database}.{keytype}.RDS",
    output:
        png=protected(
            "data_output/GSEA/{comparison}/{database}.{keytype}/dotplot.enrich.png"
        ),
    threads: 1
    resources:
        time_min=get_15min_per_attempt,
        mem_mb=get_4gb_per_attempt,
        tmpdir="tmp",
    params:
        dotplot_extra=config["clusterprofiler"].get("dotplot_extra", ""),
    log:
        "logs/030.clusterprofiler/dotplot/enrich.{db}.{comparison}.{keytype}.log",
    wrapper:
        "bio/clusterProfiler/dotplot"


# Perform bar plots on DGE results

"""
030.barplot
from
-> 029.enricher_TERMS
-> 028.enrichDO
-> 028.enrichDGN
-> 028.enrichDO
-> 027.enrich_GMT
by
-> End job
"""


rule barplot:
    input:
        rds="027.enrich/{database}/{comparison}/enrich.{database}.{keytype}.RDS",
    output:
        png=protected(
            "data_output/GSEA/{comparison}/{database}.{keytype}/dotplot.enrich.png"
        ),
    threads: 1
    resources:
        time_min=get_15min_per_attempt,
        mem_mb=get_4gb_per_attempt,
        tmpdir="tmp",
    params:
        barplot_extra=config["clusterprofiler"].get("barplot_extra", ""),
    log:
        "logs/030.clusterprofiler/barplot/{method}.{db}.{comparison}.{keytype}.log",
    wrapper:
        "bio/clusterProfiler/barplot"


# Perform dot plots on DGE results

"""
030.upsetplot
from
-> 029.enricher_TERMS
-> 028.enrichDO
-> 028.enrichDGN
-> 028.enrichDO
-> 027.enrich_GMT
by
-> End job
"""


rule upsetplot:
    input:
        rds="027.enrich/{database}/{comparison}/enrich.{database}.{keytype}.RDS",
    output:
        png=protected(
            "data_output/GSEA/{comparison}/{database}.{keytype}/upsetplot.enrich.png"
        ),
    threads: 1
    resources:
        time_min=get_15min_per_attempt,
        mem_mb=get_4gb_per_attempt,
        tmpdir="tmp",
    params:
        upsetplot_extra=config.get("upsetplot_extra", "n = 5"),
    log:
        "logs/upsetplot/{method}.{db}.{comparison}.{keytype}.log",
    wrapper:
        "bio/clusterProfiler/upsetplot"
