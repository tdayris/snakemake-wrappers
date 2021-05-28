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

log = snakemake.log_fmt_shell(stdout=False, stderr=True)

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


translation_dict = {
    "Synonymous_Variant": "chromosome_number_variation",
    "Splice_Site": "exon_loss_variant",
    "Missense_Mutation": "frameshift_variant",
    "Nonsense_Mutation": "stop_gained",
    "Missense_Mutation": "stop_lost",
    "Nonsense_Mutation": "start_lost",
    "Splice_Site": "splice_acceptor_variant",
    "Splice_Site": "splice_donor_variant",
    "Missense_Mutation": "rare_amino_acid_variant",
    "Missense_Mutation": "missense_variant",
    "Frame_Shift_Ins": "disruptive_inframe_insertion",
    "In_Frame_Ins": "conservative_inframe_insertion",
    "Frame_Shift_Del": "disruptive_inframe_deletion",
    "In_Frame_Del": "conservative_inframe_deletion",
    "5_prime_UTR_variant": "5_prime_UTR_truncation+exon_loss_variant",
    "3_prime_UTR_variant": "3_prime_UTR_truncation+exon_loss",
    "Splice_Site": "splice_branch_variant",
    "Splice_Site": "splice_region_variant",
    "Nonsense_Mutation": "stop_retained_variant",
    "Missense_Mutation": "initiator_codon_variant",
    "Synonymous_Variant": "synonymous_variant",
    "Missense_Mutation": "initiator_codon_variant+non_canonical_start_codon",
    "Nonsense_Mutation": 'stop_retained_variant',
    "Missense_Mutation": "coding_sequence_variant",
    "5_prime_UTR_variant": "5_prime_UTR_variant",
    "3_prime_UTR_variant": "3_prime_UTR_variant",
    "5_prime_UTR_variant": "5_prime_UTR_premature_start_codon_gain_variant",
    "Synonymous_Variant": "upstream_gene_variant",
    "Synonymous_Variant": "downstream_gene_variant",
    "Synonymous_Variant": "TF_binding_site_variant",
    "Synonymous_Variant": "regulatory_region_variant",
    "Synonymous_Variant": "miRNA",
    "Synonymous_Variant": "custom",
    "Missense_Mutation": "sequence_feature",
    "Splice_Site": "conserved_intron_variant",
    "Synonymous_Variant": "intron_variant",
    "Synonymous_Variant": "intragenic_variant",
    "Splice_Site": "conserved_intergenic_variant",
    "Synonymous_Variant": "intergenic_region",
    "5_prime_UTR_variant": "coding_sequence_variant",
    "Synonymous_Variant": "non_coding_exon_variant",
    "Synonymous_Variant": "nc_transcript_variant",
    "5_prime_UTR_variant": "gene_variant",
    "Synonymous_Variant": "chromosome"
}

for idx, (key, value) in enumerate(translation_dict.items()):
    replacement = "'s/\t{value}\t/\t{key}\t/g'"
    if idx == 0:
        shell("sed {replacement} {snakemake.input.maf} > {snakemake.input.maf}.tmp {log}")
    else:
        shell("sed -i {replacement} {snakemake.input.maf}.tmp {log}")

log = snakemake.log_fmt_shell(stdout=True, stderr=True)
shell("mv -v {snakemake.input.maf}.tmp {snakemake.output.maf} {log}")
