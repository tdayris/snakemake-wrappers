"""Snakemake wrapper which copies FORMAT to INFO"""

__author__ = "Thibault Dayris"
__copyright__ = "Copyright 2021, Dayris Thibault"
__email__ = "thibault.dayris@gustaveroussy.fr"
__license__ = "MIT"

import re

with (open(snakemake.input.call, "r") as instream,
      open(snakemake.output.call, "w") as outstream):
    regex_general = r'##(\w+)=(.+$)'
    regex_id = r"<ID=([^,]+),(.*)>"
    for line in instream:
        if line.startswith("##"):
            htype, hcontent = re.findall(regex_general, line)[0]
            if htype.lower() ==  "format":
                hid, hinfo = re.findall(regex_id, line)[0]
                outstream.write("##INFO=<ID=FORMAT_{},{}>\n".format(
                    hid, hinfo
                ))
            outstream.write(line)
            continue
        if line.startswith("#"):
            outstream.write(line)
        else:
            chomp = line[:-1].split("\t")

            if chomp[7] == ".":
                chomp[7] = ""
            else:
                chomp[7] += ";"

            chomp[7] += ";".join(
                "FORMAT_{}={}".format(k, v) for k, v in zip(
                    chomp[8].split(":"), chomp[9].split(":")
                )
            )
            outstream.write("\t".join(chomp) + "\n")
