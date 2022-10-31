# Perform dot plots on DGE results

"""
030.enrichplot_enricher
from
-> 029.enricher_TERMS
-> 028.enrichDO
-> 028.enrichDGN
-> 028.enrichDO
-> 027.enrich_GMT
by
-> End job
"""


rule enrichplot_enricher:
    input:
        rds="027.enrich/{database}.{keytype}/{comparison}/enrich.{database}.{keytype}.RDS",
    output:
        barplot=protected(
            "data_output/GSEA/{comparison}/{database}.{keytype}/barplot.enrich.png"
        ),
        dotplot=protected(
            "data_output/GSEA/{comparison}/{database}.{keytype}/dotplot.enrich.png"
        ),
        cnetplot=protected(
            "data_output/GSEA/{comparison}/{database}.{keytype}/cnetplot.enrich.png"
        ),
        heatplot=protected(
            "data_output/GSEA/{comparison}/{database}.{keytype}/heatplot.enrich.png"
        ),
        upsetplot=protected(
            "data_output/GSEA/{comparison}/{database}.{keytype}/upsetplot.enrich.png"
        ),
        pmcplot=protected(
            "data_output/GSEA/{comparison}/{database}.{keytype}/pmcplot.enrich.png"
        ),
    threads: 1
    resources:
        time_min=get_35min_per_attempt,
        mem_mb=get_4gb_per_attempt,
        tmpdir="tmp",
    params:
        enricher_extra="pvalueCutoff = 0.1, qvalueCutoff = 0.1",
        barplot_extra="showCategory = 5",
        dotplot_extra="showCategory = 5",
        cnetplot_extra="showCategory = 5",
        heatplot_extra="showCategory = 5",
        upsetplot_extra="n = 5",
        pmcplot_extra="period=2012:2022",
        png_extra="height = 768, width = 1024, units = 'px', type = 'cairo'",
    log:
        "logs/030.clusterprofiler/dotplot/enrich.{database}.{comparison}.{keytype}.log",
    wrapper:
        "bio/clusterprofiler/enrichplot"