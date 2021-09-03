"""Snakemake wrapper for SnpSift extractFields"""

__author__ = "Thibault Dayris"
__copyright__ = "Copyright 2020, Dayris Thibault"
__email__ = "thibault.dayris@gustaveroussy.fr"
__license__ = "MIT"

import re

from snakemake.shell import shell
from snakemake_wrapper_utils.java import get_java_opts

java_opts = get_java_opts(snakemake)

log = snakemake.log_fmt_shell(stdout=False, stderr=True)
extra = snakemake.params.get("extra", "")
min_threads = 1

# Uncompression shall be done according to user-defined input
incall = snakemake.input["call"]
if snakemake.input["call"].endswith("bcf"):
    min_threads += 1
    incall = "bcftools view {}".format(incall)
elif snakemake.input["call"].endswith("gz"):
    min_threads += 1
    incall = "gunzip -c {}".format(incall)
else:
    incall = "cat {}".format(incall)
    min_threads += 1


# Each (un)compression step raises the threads requirements
# if snakemake.threads < min_threads:
#     raise ValueError(
#         "At least {} threads required, {} provided".format(
#             min_threads, snakemake.threads
#         )
#     )

fields = ["'CHROM'", "'POS'", "'ID'", "'REF'", "'ALT'", "'FILTER'"]
help = [
    "CHROM: Chromosome name",
    "POS: Position over the chromosome",
    "ID: Variant identifier"
]
with open(snakemake.input.call, "r") as vcf_stream:
    ignore_format = snakemake.params.get("ignore_format", False)
    regex_general = r'##(\w+)=(.+$)'
    regex_id = r"<ID=([^,]+),(.*)>"
    for line in vcf_stream:
        if not line.startswith("##"):
            break

        htype, hcontent = re.findall(regex_general, line)[0]
        if htype.lower() ==  "format" and ignore_format:
            hid, hinfo = re.findall(regex_id, line)[0]
            fields.append("'{}[{}]'".format(htype, hid))
        elif htype.lower() == "info":
            hid, hinfo = re.findall(regex_id, line)[0]
            if hid.lower() == "ann":
                fields += [
                    "'ANN[*].ALLELE'",
                    "'ANN[*].EFFECT'",
                    "'ANN[*].IMPACT'",
                    "'ANN[*].GENE'",
                    "'ANN[*].GENEID'",
                    "'ANN[*].FEATURE'",
                    "'ANN[*].FEATUREID'",
                    "'ANN[*].BIOTYPE'",
                    "'ANN[*].RANK'",
                    "'ANN[*].HGVS_C'",
                    "'ANN[*].HGVS_P'",
                    "'ANN[*].CDNA_POS'",
                    "'ANN[*].CDNA_LEN'",
                    "'ANN[*].CDS_POS'",
                    "'ANN[*].CDS_LEN'",
                    "'ANN[*].AA_POS'",
                    "'ANN[*].AA_LEN'",
                    "'ANN[*].DISTANCE'",
                    "'ANN[*].ERRORS'"
                ]
            elif hid.lower() == "ref":
                fields += [
                    "'EFF[*].EFFECT'",
                    "'EFF[*].IMPACT'",
                    "'EFF[*].FUNCLASS'",
                    "'EFF[*].CODON'",
                    "'EFF[*].AA'",
                    "'EFF[*].AA_LEN'",
                    "'EFF[*].GENE'",
                    "'EFF[*].BIOTYPE'",
                    "'EFF[*].CODING'",
                    "'EFF[*].TRID'",
                    "'EFF[*].RANK'"
                ]
            elif hid.lower() == "lof":
                fields += [
                    "'LOF[*].GENE'",
                    "'LOF[*].GENEID'",
                    "'LOF[*].NUMTR'",
                    "'LOF[*].PERC'"
                ]
            elif hid.lower() == "nmd":
                fields += [
                    "'NMD[*].GENE'",
                    "'NMD[*].GENEID'",
                    "'NMD[*].NUMTR'",
                    "'NMD[*].PERC'"
                ]
            else:
                fields.append(hid)


if isinstance(fields, list):
    fields = " ".join(fields)

shell(
    "{incall} | "
    "SnpSift extractFields"  # Tool and its subcommand
    " {java_opts} {extra}"  # Extra parameters
    " - "  # Path to input vcf file
    " {fields}"  # The fields to extract
    " > {snakemake.output.tsv}"  # Path to output file
    " {log}"  # Logging behaviour
)
