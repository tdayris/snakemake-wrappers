import logging

from typing import Dict, List, Union

def get_gentrome(wildcards) -> Dict[str, str]:
    """
    Unpack function for resources.salmon.smk
    """
    genome_build = wildcards.genome_build
    genome_release = wildcards.genome_release
    config_genome = confi.get("resources", {}).get("dna_sequence_fasta")
    config_trnome = config.get("resources", {}).get("cdna_sequence_fasta")
    if config_genome and config_trnome:
        logging.debug("Using user-defined cDNA and DNA sequences")
        return {"transcriptome": config_trnome, "genome": config_genome}

    logging.debug("Using downloaded cDNA and DNA sequences")
    return {
        "transcriptome": f"resources/{genome_build}.{genome_release}.cdna.fasta",
        "genome": f"resources/{genome_build}.{genome_release}.dna.fasta",
    }


def get_bigr(wildcards, output) -> Union[List[str], str]:
    """
    Unpack function for bigr copy. The input depends on the
    output name.
    """
    sample = str(wildcards.sample_stream)
    stream = None
    if sample.endswith((".1", ".2")):
        if sample[-1] else "Downstream_file":
            stream = "Downstream_file"
        else:
            stream = "Upstream_file"
        sample = sample[:-2]

    return design.loc[sample, stream]


def get_salmon_index(wildcards) -> List[str]:
    """
    Unpack function for quantification.salmon.smk
    """
    genome_build = wildcards.genome_build
    genome_release = wildcards.genome_release
    salmon_index = config.get("index", {}).get("salmon")
    genome_gtf = config.get("reference", {}).get("gtf")
    if salmon_index and gtf:
        logging.debug("Using user-defined salmon index and GTF")
        return {"index": salmon_index, "gtf": genome_gtf}

    logging.debug("Re-building salmon index")
    return {
        "index": multiext(
            "salmon/{genome_build}.{genome_release}_index/",
            "complete_ref_lens.bin",
            "ctable.bin",
            "ctg_offsets.bin",
            "duplicate_clusters.tsv",
            "info.json",
            "mphf.bin",
            "pos.bin",
            "pre_indexing.log",
            "rank.bin",
            "refAccumLengths.bin",
            "ref_indexing.log",
            "reflengths.bin",
            "refseq.bin",
            "seq.bin",
            "versionInfo.json",
        ),
        "gtf": "resources/GRCh38.108.gtf",
    }


def get_salmon_fastq(wildcards) -> Dict[str, str]:
    """
    Unpack function for quantification.salmon.smk
    """
    if "Downstream_file" in design.columns.tolist():
        logging.debug(
            f"Salmon quant uses pair-ended mode for {wildcards.sample}"
        )
        return {
            "r1": f"020.trimming/fastp/pe/{wildcards.sample}.1.fastq",
            "r2": f"020.trimming/fastp/pe/{wildcards.sample}.2.fastq",
        }

    logging.debug(
        f"Salmon quant uses single-ended mode for {wildcards.sample}"
    )
    return {"r": f"020.trimming/fastp/se/{wildcards.sample}.fastq"}