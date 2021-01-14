#!/usr/bin/python3.8
# -*- coding: utf-8 -*-

"""Snakemake wrapper for bash sed"""

__author__ = "Thibault Dayris"
__copyright__ = "Copyright 2019, Thibault Dayris"
__email__ = "thibault.dayris@gustaveroussy.fr"
__license__ = "MIT"


from snakemake.shell import shell

extra = snakemake.params.get("extra", "")
if "-i" in extra:
    raise ValueError("Do not use -i with this wrapper.")

temp_outfile = snakemake.params.get("tmp", "sed.temp")
# fs = field_separator
fs = snakemake.params.get("field_separator", "\t")
intermediar = []

# Case of a simple csv to tsv conversion
if snakemake.params.get("csv_to_tsv", False) is True:
    intermediar.append("'s/,/\t/g'")
    field_separator = "\t"

# Case of a path removal
if snakemake.params.get("basename", False) is True:
    intermediar.append("'s/{fs}[^{fs}]*\//{fs}/g'")


if snakemake.params.get("wrap_lines") is not None:
    intermediar.append(
        "-E ':x {N ; s/\n/ /g ; s/(.{80,80})/\1\n/ ; /\n/!bx ; P ; D}'"
    )


if snakemake.params.get("split_on_commas") is not None:
    intermediar.append(
        "-E ':x {N ; s/\n/ /g ; s/,/,\n/ ; /\n/!bx ; s/^ *// ; P ; D}'"
    )


if snakemake.params.get("join_on_backslashes") is not None:
    intermediar.append(
        "-e ':x /\\$/ { N; s/\\\n//g ; bx }'"
    )


if snakemake.params.get("join_on_spaces") is not None:
    intermediar.append(
        "-E ':a ; $!N ; s/\n\s+/ / ; ta ; P ; D'"
    )


# Case of line removal on several criteria
if (remove_list := snakemake.params.get("remove_list")) is not None:
    for element_to_remove in remove_list:
        intermediar.append(f"'/{element_to_remove}/d'")
if (remove_lines := snakemake.params.get("remove_lines")) is not None:
    for line_to_remove in remove_lines:
        intermediar.append(f"'{line_to_remove}d'")

# Case of series of replacements
if (replace_dict := snakemake.params.get("replace_dict")) is not None:
    for replace, replacement in replace_dict.items():
        intermediar.append(f"'s/{replace}/{replacement}/g'")

# Other cases
if (regex := snakemake.params.get("regex")) is not None:
    if isinstance(regex, str):
        intermediar.append(regex)
    elif isinstance(regex, list):
        intermediar += regex

# Case user did not provide anything
if len(intermediar) == 0:
    raise ValueError("Could not find any sed regex to apply.")


# Apply sed one after each other.
sed_iterator = iter(intermediar)
sed_value = next(sed_iterator, None)
log = snakemake.log_fmt_shell(stdout=False, stderr=True)
shell("sed {extra} {sed_value} {snakemake.input} > {temp_outfile} {log}")

sed_value = next(sed_iterator, None)
log = snakemake.log_fmt_shell(stdout=True, stderr=True, append=True)
while sed_value is not None:
    shell("sed -i {extra} {sed_value} {temp_outfile} {log}")
    sed_value = next(sed_iterator, None)

shell("mv -v {temp_outfile} {snakemake.output[0]} {log}")
