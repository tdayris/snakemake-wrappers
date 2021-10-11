"""Snakemake wrapper which copies FORMAT to INFO"""

__author__ = "Thibault Dayris"
__copyright__ = "Copyright 2021, Dayris Thibault"
__email__ = "thibault.dayris@gustaveroussy.fr"
__license__ = "MIT"

import datetime
import logging

def create_header(sample_name, field_name):
    return f"""##INFO=<ID=FORMAT_{sample_name}_{field_name},Number=.,Type=String,Description="Copy of the value of the field {field_name} from FORMAT concerning sample/tool {sample_name}">"""

colnames = None
version = 1.0
name = "vcf_format_to_info"
url = f"github.com/tdayris/snakemake-wrappers/tree/Unofficial/bio/BiGR/{name}/wrapper.py"
header = f'##BiGRCommandLine=<ID={name},CommandLine="{url}",Version={version},Date={datetime.date.today()}>\n'
format_headers = []

logging.info("Looking for all possible format fields among all samples")
with open(snakemake.input.call, "r") as instream:
    for line in instream:
        if line.startswith("##"):
            continue

        if line.startswith("#"):
            colnames = line[1:-1].split("\t")
            continue

        formats = line.split("\t")[8]
        samples = colnames[9:]
        for sample in samples:
            for format in formats.split(":"):
                format_headers.append(create_header(sample, format))

header += "\n".join(set(format_headers)) + "\n"
logging.debug(header)
logging.info("Formats were built, now annotating VCF")

with (open(snakemake.input.call, "r") as instream,
      open(snakemake.output.call, "w") as outstream):

    for line in instream:
        if line.startswith("##"):
            outstream.write(line)
            continue
        if line.startswith("#"):
            colnames = line[1:-1].split("\t")
            outstream.write(header)
            outstream.write(line)
        else:
            chomp = line[:-1].split("\t")

            for sample in colnames[9:]:
                if chomp[7] == ".":
                    chomp[7] = ""
                else:
                    chomp[7] += ";"

                chomp[7] += ";".join(
                    "FORMAT_{}_{}={}".format(sample, k, v) for k, v in zip(
                        chomp[8].split(":"), chomp[9].split(":")
                    )
                )
            outstream.write("\t".join(chomp) + "\n")
