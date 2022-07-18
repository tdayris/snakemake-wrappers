from snakemake.utils import min_version
from pathlib import Path
from yaml import dump
min_version("7.0")

import sys

worflow_source_dir = Path(snakemake.workflow.srcdir("."))
common = str(worflow_source_dir / ".." / "common" / "python")
sys.path.append(common)

from dataframes import *
from file_manager import *
from files_linker import *
from graphics import *
from write_yaml import *
from messages import message


#################
### Preambule ###
#################

logging.basicConfig(
    filename="snakemake.dge_deseq2.log",
    filemode="w",
    level=logging.DEBUG
)

default_config = read_yaml(worflow_source_dir / "config.hg38.yaml")
configfile: get_config(default_config)

try:
    design = pandas.read_csv("design.tsv", sep="\t", header=0, index_col=0)
    design["Sample_id"] = design.index.tolist()
    #design.set_index("Sample_id", inplace=True)
except FileNotFoundError:
    logging.error(
        """A design file is required for this pipeline. It is a TSV with
        the following columns:

        1. Sample_id (case matters): Name of your sample, unique and composed
           with a least 1 letter (no sample should have numerical names only,
           it would make R fail while parsing sample names with DESeq2)
        2. Upstream_file (case matters): Path to the file, it can be
           an absolute path, a relative path, or a iRODS url.
        3. Downstream_file (case matters): Path to the file, it can be
           an absolute path, a relative path, or a iRODS url.
        4. XXXX: A name of your choice, unique and understandable. It will be
           used as comparison name within DESeq2 and graphs. It contains levels
           for each single sample. Do not use only integers or floats for your
           level name: R and DESeq2 behaves stangely with them.
        5. YYYY: A name of your choice, unique and understandable. It will be
           used as comparison name within DESeq2 and graphs. It contains levels
           for each single sample. Do not use only integers or floats for your
           level name: R and DESeq2 behaves stangely with them.
        Etc, etc. You can have any other condition name. Name them as you want,
        these names must be unique and understandable. It will be used as
        comparison name within DESeq2 and graphs. It contains levels for each
        single sample. Do not use only integers or floats for your level name:
        R and DESeq2 behaves stangely with them.
        """
    )

fastq_links = link_fq(
    design.Sample_id,
    design.Upstream_file,
    design.Downstream_file
)

# A list that holds all comparisons expected for this snakemake pipeline
comparison_levels = list(yield_comps(
    complete_design=design,
    aggregate=config["design"].get("aggregate_col"),
    remove=config["design"].get("remove_col"),
    contains=config["design"].get("include_only")
))

# Stored as a list for futrther re-use
output_prefixes = [
    f"DGE_considering_factor_{factor}_comparing_{test}_vs_{ref}"
    for factor, test, ref in comparison_levels
]
# print(output_prefixes)

# An iterator that holds all samples involved in the comparisons
# listed above
samples_iterator = yield_samples(
    complete_design=design.copy(),
    aggregate=config["design"].get("aggregate_col"),
    remove=config["design"].get("remove_col")
)

expected_pcas = [
    f"figures/DGE_considering_factor_{factor}_comparing_{test}_vs_{ref}/pca/pca_{factor}_{axes}_{elipse}.png"
    for (factor, test, ref) in comparison_levels
    for axes in ["ax_1_ax_2", "ax_2_ax_3"] # , "ax_3_ax_4"]
    for elipse in ["with_elipse", "without_elipse"]
]

condition_dict = {
    f"DGE_considering_factor_{factor}_comparing_{test}_vs_{ref}": relation_condition_sample(design.copy(), factor, test, ref)
    for factor, test, ref in comparison_levels
}

samples_per_prefixes = dict(zip(output_prefixes, condition_dict))
samples_per_prefixes = {
    prefix: list(condition_dict[prefix].keys())
    for prefix in output_prefixes
}
logging.debug(samples_per_prefixes)


############################
### Wilcards constraints ###
############################

wildcard_constraints:
    comparison=r"|".join(output_prefixes),
    factor=r"|".join(map(str, [i[0] for i in comparison_levels])),
    test=r"|".join(map(str, [i[1] for i in comparison_levels])),
    ref=r"|".join(map(str, [i[2] for i in comparison_levels])),
    axes=r"|".join(["ax_1_ax_2", "ax_2_ax_3", "ax_3_ax_4"]),
    elipse=r"|".join(["with_elipse", "without_elipse"])


###################
### Target rule ###
###################

rule target_dge_deseq:
    input:
        multiqc=expand(
            "results/{comparison}/MultiQC.{comparison}.html",
            comparison=output_prefixes
        ),
        gseaapp=expand(
            "results/{comparison}/deseq2_{subset}_{comparison}.tsv",
            comparison=output_prefixes,
            subset=["complete_results", "sorted_on_fold_change", "sorted_on_pval"]
        ),
        csv_report=expand(
            "results/{comparison}/html_table_deseq2_{subset}.tar.bz2",
            comparison=output_prefixes,
            subset=["complete_results", "sorted_on_fold_change", "sorted_on_pval"]
        ),
        deseq2_wald=expand(
            "deseq2/{comparison}/wald.{comparison}.RDS",
            comparison=output_prefixes
        ),
        pcas=expected_pcas,
        general_pcas=expand(
            "figures/pca/general.pca_{axes}.png",
            axes=["PC1_PC2", "PC2_PC1"]
        ),
        counts_with_deseq2=expand(
            "results/{comparison}/deseq2_{content}_with_counts_{comparison}.tsv",
            comparison=output_prefixes,
            content=["complete_results", "sorted_on_pval", "sorted_on_fold_change"]
        ),
        gene_plots=expand(
            "results/{comparison}/gene_plots/{gene}.png",
            comparison=output_prefixes,
            gene=config.get("genes_of_interest", ["ENSG00000141510"])
        )
        #consensus=expand(
        #    "consensusclusterplus/{comparison}",
        #    comparison=output_prefixes
        #)


##############################
### DESeq2 post processing ###
##############################


include: "rules/003.tximport_deseq2.smk"


###########################
### tximprot and DESeq2 ###
###########################

include: "rules/002.tximport_deseq2.smk"


# #############################
# ### Salmon quantification ###
# #############################


include: "rules/001.salmon_quant_pipeline.smk"
