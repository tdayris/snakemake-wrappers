"""Snakemake wrapper for agat conversion scripts"""

__author__ = "Thibault Dayris"
__copyright__ = "Copyright 2023, Thibault Dayris"
__email__ = "thibault.dayris@gustaveroussy.fr"
__license__ = "MIT"

import os.path

from snakemake.shell import shell
from tempfile import TemporaryDirectory

log = snakemake.log_fmt_shell(stdout=True, stderr=True, append=True)
extra = snakemake.params.get("extra", "")


# Acquire agat script name
script = str(snakemake.params["script"])
if not script.lower().strip().endswith(".pl"):
    script += ".pl"


# Acquire input files
input_cmd = ""
for argname, file_path in snakemake.input.items():
    # Deal with long/short options, since some agat scripts
    # only accepts 'input' as long option and input.input is
    # a protected name in Snakemake.
    dash = "-" if len(str(argname)) == 1 else "--"
    if isinstance(file_path, list):
        input_cmd += " ".join([
            f" {dash}{argname} {fp} "
            for fp in file_path
        ])
    else:
        input_cmd += f" {dash}{argname} {file_path} "

# Acquire output file(s)
# output_cmd = ""
# if script == "agat_convert_sp_gff2zff.pl":
#     output_cmd += f""
#     prefix = os.path.commonprefix(list(map(str, snakemake.output)))[:-1]
#     output_cmd += f" --output {prefix} "
# elif script == "agat_sp_extract_attributes.pl":
#     ext = os.path.splitext(snakemake.output[0])[-1]
#     prefix = os.path.commonprefix(list(map(str , snakemake.output)))[:-1]
#     output_cmd += f" --output {prefix}{ext} "

# else:



with TemporaryDirectory() as tempdir:
    output_cmd = ""
    # The following script do not let user choose output file names
    if script in ["agat_convert_sp_gff2zff.pl", "agat_sp_extract_attributes.pl"]:
        output_cmd += f" --output {tempdir}/outfile "
    else:
        for argname, file_path in snakemake.output.items():
            # Deal with long/short options, since some agat scripts
            # only accepts 'output' as long option and output.output is
            # a protected name in Snakemake.
            dash = "-" if len(str(argname)) == 1 else "--"
            output_cmd += f" {dash}{argname} {file_path} "

    shell("{script} {extra} {input_cmd} {output_cmd} {log}")
    shell("ls -lrth {tempdir} {log}")
    shell("ls -lrth {log}")
    shell("pwd {log}")

    # Forwarding output files for script that do not
    # let user choose output file name(s)
    if script == "agat_convert_sp_gff2zff.pl":
        annot = snakemake.output.get("annotation")
        if annot:
            shell(f"mv --verbose {tempdir}/outfile.ann {annot} {log}")
        
        fasta = snakemake.output.get("fasta")
        if fasta:
            shell(f"mv --verbose {tempdir}/outfile.dna {fasta} {log}")

    elif script == "agat_sp_extract_attributes.pl":
        extra_args = iter(extra.split(" "))
        extra = next(extra_args, None)   
        
        # Acquire all expected output files
        fields = []
        while extra is not None:
            if extra in ["--attribute", "--att", "-a"]:
                fields = next(extra_args, None).split(",")
                break

            next(extra_args, None)
        
        # This expects user listed output file(s) in the same
        # order as the attributes listed in `extra`
        for field, outfile in zip(fields, snakemake.output):
            shell(f"mv --verbose {tempdir}/outfile_{field} {outfile} {log}")
            
        
