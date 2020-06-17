#!/usr/bin/python3.8
# conding: utf-8

"""
Set up a pcaExplorer launcher script
"""

script = f"""#!/usr/bin/R

# This script is produced automatically, please check-it before any execution

base::library(
    package = "DESeq2",
    quietly = TRUE
);

base::library(
    package = "pcaExplorer",
    quietly = TRUE
);
base::message(
    'Libraries loaded. Now loading input datasets.'
);

dds <- base::readRDS(
    file = '{snakemake.input["dds"]}'
);
dst <- base::readRDS(
    file = '{snakemake.input["dst"]}'
);
coldata <- utils::read.table(
    file = '{snakemake.input["coldata"]}'
);
annotation <- base::readRDS(
    file = '{snakemake.input["annotation"]}'
);
pca2go <- base::readRDS(
    file = '{snakemake.input["pca2go"]}'
);
base::message(
    'Input files loaded. Now launching pcaExplorer.'
);

pcaExplorer::pcaExplorer(
    dds = dds,
    dst = dst,
    coldata = coldata,
    pca2go = pca2go,
    annotation = annotation
);
"""


with open(snakemake.output["script"], "w") as out_script:
    out_script.write(script)
