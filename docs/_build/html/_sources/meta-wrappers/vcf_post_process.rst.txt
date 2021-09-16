.. _`vcf_post_process`:

VCF_POST_PROCESS
================

Format gzipped VCF to TSV with renamed columns and human-friendly content when possible.


Example
-------

This meta-wrapper can be used by integrating the following into your workflow:

.. code-block:: python

    import sys
    from pathlib import Path

    worflow_source_dir = Path(next(iter(workflow.get_sources()))).absolute().parent
    common = str(worflow_source_dir / "../../../../bigr_pipelines/common/python")
    sys.path.append(common)

    from file_manager import *

    default_config_vcf_post_process = {
        "ncbi_build": "GRCh38",
        "center": "GustaveRoussy",
        "annotation_tag": "ANN=",
        "sample_list": list(),
        "genome": "/path/to/ref.fasta"
    }

    try:
        if config == dict():
            config = default_config_vcf_post_process
    except NameError:
        config = default_config_vcf_post_process



    """
    Big concatenation for maftools
    """
    rule concat_all_mafs:
        input:
            expand(
                "maf/maftools/{sample}.maf",
                sample=config["sample_list"]
            )
        output:
            "maf/maftools/complete.maf"
        message:
            "Merging separates maf files in a single cohort maf"
        threads: 1
        resources:
            mem_mb=lambda wildcards, attempt: attempt * 1024,
            time_min=lambda wildcards, attempt: attempt * 15,
            tmpdir="tmp"
        log:
            "logs/maftools/concat.log"
        shell:
            "cat {input} > {output} 2> {log}"


    """
    Assure MAFtools compatibility and human readability
    """
    rule rename_snpsift_maf_cols:
        input:
            tsv="maf/extracted/{sample}.tsv"
        output:
            tsv="maf/maftools/{sample}.maf"
        message:
            "Renaming columns to fit MAFtools requirement in {wildcards.sample}"
        threads: 1
        resources:
            mem_mb=lambda wildcards, attempt: attempt * 1024,
            time_min=lambda wildcards, attempt: attempt * 15,
            tmpdir="tmp"
        log:
            "logs/maftools/rename/{sample}.log"
        params:
            add_cols=True,
            ncbi_build=config.get("NBCI_build", "GRCh38"),
            center=config.get("center", "GustaveRoussy"),
            Tumor_Sample_Barcode=lambda wildcards: f"{wildcards.sample}_tumor",
            Matched_Norm_Sample_Barcode=lambda wildcards: f"{wildcards.sample}_normal"
        wrapper:
            "bio/BiGR/rename_snpsift_maf_cols"


    """
    Extracting all INFO/FORMAT data, the list is built from vcf header
    """
    rule extract_all_fields:
        input:
            call="maf/splitted/{sample}.vcf"
        output:
            tsv=temp("maf/extracted/{sample}.tsv")
        message:
            "Extracting variant annotations for {wildcards.sample}"
        threads: 2
        resources:
            mem_mb=lambda wildcards, attempt: attempt * 10240,
            time_min=lambda wildcards, attempt: attempt * 25,
            tmpdir="tmp"
        log:
            "logs/snpsift/extract_all_fields/{sample}.log"
        params:
            annotation_tag=config.get("annotation_tag", "ANN=")
        wrapper:
            "bio/snpsift/extractAllFields"


    """
    Split annotation since it may lead to errors in MAFtools and/or in a result
    file more difficult to read by a human.
    """
    rule split_vcf_features:
        input:
            call="gatk/variant_filtration/{sample}.vcf"
        output:
            call=temp("maf/splitted/{sample}.vcf")
        message:
            "Splitting variant annotations for {wildcards.sample}"
        threads: 1
        resources:
            mem_mb=lambda wildcards, attempt: attempt * 1024,
            time_min=lambda wildcards, attempt: attempt * 15,
            tmpdir="tmp"
        log:
            "logs/split_vcf_features/{sample}.log"
        params:
            annotation_tag="ANN="
        wrapper:
            "bio/BiGR/split_vcf_features"


    """
    Add filter tags
    """
    rule gatk_variant_filtration:
        input:
            vcf="snpsift/annotate_corrected/{sample}.vcf.gz",
            vcf_tbi=get_tbi("snpsift/annotate_corrected/{sample}.vcf.gz"),
            ref=config["genome"]
        output:
            vcf="gatk/variant_filtration/{sample}.vcf"
        message:
            "Filtering VCF for {wildcards.sample}"
        threads: 1
        resources:
            mem_mb=lambda wildcards, attempt: attempt * 10240,
            time_min=lambda wildcards, attempt: attempt * 25,
            tmpdir="tmp"
        log:
            "logs/gatk/variant_filtration/{sample}.log"
        params:
            filters={
                "Depth60X": "DP > 59",
                "VAF10pct": "AF >= 0.1",
                "GeneralPopulation": "dbNSFP_ExAC_Adj_AF <= 0.001",
                "VAF5pct": "AF >= 0.05",
            }
        wrapper:
            "bio/gatk/variantfiltration"


    rule fix_SnpSift_Annotate:
        input:
            call="snpsift/dbnsfp/{sample}.vcf.gz",
            tbi=get_tbi("snpsift/dbnsfp/{sample}.vcf.gz")
        output:
            call=temp("snpsift/annotate_corrected/{sample}.vcf")
        message:
            "Correcting annotation type error in {wildcards.sample}"
        threads: 4
        resources:
            mem_mb=lambda wildcards, attempt: attempt * 1024,
            time_min=lambda wildcards, attempt: attempt * 15,
            tmpdir="tmp"
        params:
            remove_ends="'s/END=\([0-9]\+,\?\)\+//g'",
            remove_empty="'s/;;/;/g'",
            remove_eofield="'s/;\t/\t/g'"
        log:
            "logs/snpsift/annotate_corrected/{sample}.log"
        shell:
            "gunzip -c {input.call} | "
            "sed {params.remove_ends} | "
            "sed {params.remove_empty} | "
            "sed {params.remove_eofield} "
            "> {output.call} 2> {log}"

Note that input, output and log file paths can be chosen freely, as long as the dependencies between the rules remain as listed here.
For additional parameters in each individual wrapper, please refer to their corresponding documentation (see links below).

When running with

.. code-block:: bash

    snakemake --use-conda

the software dependencies will be automatically deployed into an isolated environment before execution.



Used wrappers
---------------------

The following individual wrappers are used in this meta-wrapper:


* :ref:`bio/BiGR/split_vcf_features`

* :ref:`bio/snpsift/extractAllFields`

* :ref:`bio/BiGR/rename_snpsift_maf_cols`


Please refer to each wrapper in above list for additional configuration parameters and information about the executed code.






Notes
-----

VCF header will be used to extract information, no need to list them all.




Authors
-------


* Thibault Dayris

