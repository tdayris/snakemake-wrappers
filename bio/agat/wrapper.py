"""Snakemake wrapper for agat conversion scripts"""

__author__ = "Thibault Dayris"
__copyright__ = "Copyright 2023, Thibault Dayris"
__email__ = "thibault.dayris@gustaveroussy.fr"
__license__ = "MIT"

import shlex

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


with TemporaryDirectory() as tempdir:
    output_cmd = ""
    # The following scripts do not let user choose output file names
    if script in ["agat_sp_extract_attributes.pl"]:
        output_cmd += f" --output {tempdir}/outfile "
    elif script in ["agat_sp_filter_by_ORF_size.pl"]:
        output_cmd += f" --output {tempdir}/outfile.gff "
    else:
        for argname, file_path in snakemake.output.items():
            # Deal with long/short options, since some agat scripts
            # only accepts 'output' as long option and output.output is
            # a protected name in Snakemake.
            argname = "output" if script in ["agat_convert_sp_gff2zff.pl"] else argname
            dash = "-" if len(str(argname)) == 1 else "--"
            output_cmd += f" {dash}{argname} {file_path} "

    shell("{script} {extra} {input_cmd} {output_cmd} {log}")
    shell("echo 'ls {tempdir}...'")
    shell("ls -lrth {tempdir} {log}")
    shell("echo 'ls local...'")
    shell("ls -lrth {log}")
    shell("echo 'tree {tempdir}...'")
    shell("tree {tempdir} {log}")
    shell("echo 'we are at...'")
    shell("pwd {log}")

    # Forwarding output files for script that do not
    # let user choose output file name(s)
    if script == "agat_sp_extract_attributes.pl":
        extra_args = iter(shlex.split(extra))
        extra = next(extra_args, None)   
        
        # Acquire all expected output files
        fields = []
        while extra is not None:
            if extra in ["--attribute", "--att", "-a"]:
                if "=" in extra:
                    fields = extra.split("=")[-1]
                else:
                    fields = next(extra_args, None).split(",")
                break

            next(extra_args, None)
        
        # This expects user listed output file(s) in the same
        # order as the attributes listed in `extra`
        for field, outfile in zip(fields, snakemake.output):
            shell(f"mv --verbose {tempdir}/outfile_{field} {outfile} {log}")

    elif script == "agat_sp_filter_by_ORF_size.pl":
        filters = ["sup", "NOT_sup"]
        prot_length = 100

        # Acquiering protein size chosen by user (if any)
        extra_args = iter(shlex.split(extra))
        extra = next(extra_args, None)

        while extra is not None:
            if extra in ["--size", "-s"]:
                if "=" in extra:
                    prot_length = extra.split("=")[-1]
                else:
                    prot_length = next(extra_args, None).split(",")
                break

            next(extra_args, None)
        
        # Move expected output files as uses intends
        passing = snakemake.output.get("passing")
        if passing:
            shell("mv --verbose {tempdir}/outfile_{filters[0]}{prot_length}.gff {passing} {log}")
        
        not_passing = snakemake.output.get("not_passing")
        if not_passing:
            shell("mv --verbose {tempdir}/outfile_{filters[1]}{prot_length}.gff {not_passing} {log}")
            
        
