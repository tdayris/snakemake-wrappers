#!/usr/bin/python3.8
# conding: utf-8

"""
Rename variant effects from SnpEff to default MafTools
"""

__author__ = "Thibault Dayris"
__copyright__ = "Copyright 2021, Thibault Dayris"
__email__ = "thibault.dayris@gustaveroussy.fr"
__license__ = "MIT"

import logging
import pandas
import numpy
# import snakemake

from snakemake.shell import shell

logging.basicConfig(
    filename=snakemake.log[0],
    filemode="w",
    level=logging.DEBUG
)

# log = snakemake.log_fmt_shell(stdout=False, stderr=True)

accepted_terms = [
    "Frame_Shift_Del",
    "Frame_Shift_Ins",
    "In_Frame_Del",
    "In_Frame_Ins",
    "Missense_Mutation",
    "Nonsense_Mutation",
    "Splice_Site",
    "Synonymous_Variant",
    "5_prime_UTR_variant",
    "3_prime_UTR_variant"
]

def translate(snpeff_term: str) -> str:
    """
    Guess the best MAFTools term from SnpEff/SnpSift very
    precise annotation
    """
    snpeff_term = snpeff_term.split("&")

    nonsense = [
        "stop_retained_variant",
        "stop_lost",
        "start_lost",
        "stop_gained",
        "stop_retained_variant",
    ]
    if any(i in nonsense for i in snpeff_term):
        return "Nonsense_Mutation"

    frame_ins = [
        "inframe_insertion",
    ]
    if any(i in frame_ins for i in snpeff_term):
        return "In_Frame_Ins"

    frame_del = [
        "inframe_deletion"
    ]
    if any(i in frame_del for i in snpeff_term):
        return "In_Frame_Del"

    frame_shift_ins = [
        "duplication", # Duplication of a large chromosome segment (over 1% or 1,000,000 bases)
        "disruptive_inframe_insertion",
        "frameshift_variant",
    ]
    if any(i in frame_shift_ins for i in snpeff_term):
        return "Frame_Shift_Ins"

    frame_shift_del = [
        "chromosome", # A large parte (over 1%) of the chromosome was deleted.
        "disruptive_inframe_deletion",
        "exon_loss_variant",
        "feature_ablation",

    ]
    if any(i in frame_shift_del for i in snpeff_term):
        return "Frame_Shift_Del"

    splice = [
        "splice_acceptor_variant",
        "splice_donor_variant",
        "splice_region_variant",
        "exon_loss",
        "exon_loss_variant",
    ]
    if any(i in splice for i in snpeff_term):
        return "Splice_Site"

    p5 = [
        "5_prime_UTR_premature_start_codon_gain_variant",
        "start_retained",
        "5_prime_UTR_variant",
        "5_prime_UTR_truncation"
    ]
    if any(i in p5 for i in snpeff_term):
        return "5_prime_UTR_variant"

    p3 = [
        "3_prime_UTR_variant",
        "3_prime_UTR_truncation"
    ]
    if any(i in p3 for i in snpeff_term):
        return "3_prime_UTR_variant"

    missense = [
        "coding_sequence_variant",
        "inversion", # Inversion of a large chromosome segment (over 1% or 1,000,000 bases).
        "exon_variant",
        "gene_variant",
        "gene_fusion",
        "bidirectional_gene_fusion",
        "rearranged_at_DNA_level",
        'missense_variant',
        'initiator_codon_variant',
        "structural_interaction_variant",
        "rare_amino_acid_variant",
        "transcript_variant",
    ]
    if any(i in missense for i in snpeff_term):
        return "Missense_Mutation"
    return "Synonymous_Variant"


logging.debug("Reading MAF file")
df = pandas.read_csv(snakemake.input["maf"], sep="\t", header=0, index_col=None)
logging.debug("Maf file loaded")
df["Variant_Classification"] = [
    translate(v) for v in df["Variant_Classification"]
]
logging.debug("Variant_Classification translated")
df.to_csv(snakemake.output["maf"], sep="\t", index=False)
logging.debug("Variant_Classification translated")

# for idx, (key, value) in enumerate(translation_dict.items()):
#     replacement = "'s/{value}/{key}/g'"
#     if idx == 0:
#         shell("sed {replacement} {snakemake.input.maf} > {snakemake.input.maf}.tmp {log}")
#     else:
#         shell("sed -i {replacement} {snakemake.input.maf}.tmp {log}")
#
# log = snakemake.log_fmt_shell(stdout=True, stderr=True)
# shell("mv -v {snakemake.input.maf}.tmp {snakemake.output.maf} {log}")