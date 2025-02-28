#Sources:

Find these pipelines on Flamingo at: /mnt/beegfs/pipelines/unofficial-snakemake-wrappers/1.1.0/ 

In doubt, use the latest available version.

# Pipelines available in Gustave Roussy only:

1. [RNA-seq](https://github.com/tdayris/snakemake-wrappers/tree/Unofficial/bigr_pipelines#rna-seq-bulk) bulk
1. [WES](https://github.com/tdayris/snakemake-wrappers/tree/Unofficial/bigr_pipelines#somatic-wes), for normal/tumor pairs

See below for more information

## RNA-Seq bulk

[Documentation summary](https://github.com/tdayris/snakemake-wrappers/blob/Unofficial/bigr_pipelines/rnaseq/readme.md):

1. [TLDR](https://github.com/tdayris/snakemake-wrappers/tree/Unofficial/bigr_pipelines/rnaseq#tldr-run-this-pipeline)
1. [design.tsv file](https://github.com/tdayris/snakemake-wrappers/tree/Unofficial/bigr_pipelines/rnaseq#design-file)
1. [config.yaml file](https://github.com/tdayris/snakemake-wrappers/tree/Unofficial/bigr_pipelines/rnaseq#configyaml)
1. [Classical use](https://github.com/tdayris/snakemake-wrappers/tree/Unofficial/bigr_pipelines/rnaseq#classical-use)
1. [Quality controls](https://github.com/tdayris/snakemake-wrappers/tree/Unofficial/bigr_pipelines/rnaseq#quality-controls)
1. [Quantification](https://github.com/tdayris/snakemake-wrappers/tree/Unofficial/bigr_pipelines/rnaseq#quantification)
1. [Differential Gene Expression](https://github.com/tdayris/snakemake-wrappers/tree/Unofficial/bigr_pipelines/rnaseq#differential-gene-expression)
1. [Aggregate factors](https://github.com/tdayris/snakemake-wrappers/tree/Unofficial/bigr_pipelines/rnaseq#aggregate-factors)
1. [Ignore factors](https://github.com/tdayris/snakemake-wrappers/tree/Unofficial/bigr_pipelines/rnaseq#ignore-factors)
1. [Perform only a subset of DGE](https://github.com/tdayris/snakemake-wrappers/tree/Unofficial/bigr_pipelines/rnaseq#perform-only-a-subset-of-dge)
1. [Immune cell population deconvolution]() Under construction
1. [Fusions](https://github.com/tdayris/snakemake-wrappers/tree/Unofficial/bigr_pipelines/rnaseq#fusions) Under construction
1. [Variant Calling](https://github.com/tdayris/snakemake-wrappers/tree/Unofficial/bigr_pipelines/rnaseq#variant-calling) Under construction


## Somatic WES

[Documentation summary](https://github.com/tdayris/snakemake-wrappers/tree/Unofficial/bigr_pipelines/wes_somatic):

1. [TLDR](https://github.com/tdayris/snakemake-wrappers/tree/Unofficial/bigr_pipelines/wes_somatic#tldr-run-this-pipeline)
1. [design.tsv file](https://github.com/tdayris/snakemake-wrappers/tree/Unofficial/bigr_pipelines/wes_somatic#design-file)
1. [config.yaml file](https://github.com/tdayris/snakemake-wrappers/tree/Unofficial/bigr_pipelines/wes_somatic#config-file)
1. [Classical use](https://github.com/tdayris/snakemake-wrappers/tree/Unofficial/bigr_pipelines/wes_somatic#classical-use)
1. [Quality controls](https://github.com/tdayris/snakemake-wrappers/tree/Unofficial/bigr_pipelines/wes_somatic#quality-controls)
1. [Mapping only](https://github.com/tdayris/snakemake-wrappers/tree/Unofficial/bigr_pipelines/wes_somatic#mapping-only)
1. [Variant calling](https://github.com/tdayris/snakemake-wrappers/tree/Unofficial/bigr_pipelines/wes_somatic#variant-calling)
1. [CNV calling](https://github.com/tdayris/snakemake-wrappers/tree/Unofficial/bigr_pipelines/wes_somatic#cnv-calling)
1. [Tumor Molecular Burden](https://github.com/tdayris/snakemake-wrappers/tree/Unofficial/bigr_pipelines/wes_somatic#tumor-molecular-burden)
1. [Micro satellites instability](https://github.com/tdayris/snakemake-wrappers/tree/Unofficial/bigr_pipelines/wes_somatic#miscosatellites-instability)
1. [Fusions](https://github.com/tdayris/snakemake-wrappers/tree/Unofficial/bigr_pipelines/wes_somatic#fusions) Under construction

# The Snakemake Wrapper Repository

The Snakemake Wrapper Repository is a collection of reusable wrappers that allow to quickly use popular command line tools
from [Snakemake](https://snakemake.readthedocs.io) rules and workflows.

Visit https://snakemake-wrappers.readthedocs.io for more information.
