"""Snakemake wrapper for agat conversion scripts"""

__author__ = "Thibault Dayris"
__copyright__ = "Copyright 2023, Thibault Dayris"
__email__ = "thibault.dayris@gustaveroussy.fr"
__license__ = "MIT"

import shlex

from snakemake.shell import shell
from tempfile import TemporaryDirectory

log = snakemake.log_fmt_shell(stdout=True, stderr=True, append=True)


def move_additional_file(snakemake_key, source, log=log) -> None:
    """
    Move the given target to the given location if target exists.
    Raises an error destination is not reachable.
    """
    destination = snakemake.output.get(snakemake_key)
    if destination:
        shell(f"mv --verbose {source} {destination} {log}")


extra = snakemake.params.get("extra", "")
if any(str(outfile).endswith(".pdf") for outfile in snakemake.output):
    extra += "--plot"


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

    # Deal with some agat function allowing multiple input files
    # with argument repetition.
    if isinstance(file_path, list):
        input_cmd += " ".join([f" {dash}{argname} {fp} " for fp in file_path])
    else:
        input_cmd += f" {dash}{argname} {file_path} "


with TemporaryDirectory() as tempdir:
    output_cmd = ""
    outfile_exceptions = [
        "agat_sp_extract_attributes.pl",
        "agat_sp_fix_longest_ORF.pl",
        "agat_sp_manage_UTRs.pl",
        "agat_sp_manage_functional_annotation.pl",
        "agat_sp_manage_introns.pl",
        "agat_sp_prokka_fix_fragmented_gene_annotations.pl",
    ]
    # The following scripts do not let user choose output file names
    if script in outfile_exceptions:
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

    if script == "agat_sp_prokka_fix_fragmented_gene_annotations.pl":
        if "fasta" in snakemake.output.keys():
            extra += " --frags "
        if "gff" in snakemake.output.keys():
            extra += " --pseudo "

    shell("{script} {extra} {input_cmd} {output_cmd} {log}")
    shell("echo 'ls {tempdir}...'")
    shell("ls -lrth --color  {tempdir} {log}")
    shell("echo 'ls local...'")
    shell("ls -lrth --color ${{PWD}} {log}")
    shell("echo 'tree {tempdir}...'")
    shell("tree -Rf {tempdir} {log}")
    shell("tree -Rf {log}")
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
            move_additional_file(outfile, f"{tempdir}/outfile_{field}")

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
        move_additional_file(
            "passing", f"{tempdir}/outfile_{filters[0]}{prot_length}.gf"
        )
        move_additional_file(
            "not_passing", f"{tempdir}/outfile_{filters[1]}{prot_length}.gff"
        )

    elif script == "agat_sp_fix_longest_ORF.pl":
        move_additional_file("modified", f"{tempdir}/outfile-only_modified.gff")
        move_additional_file("intact", f"{tempdir}/outfile-intact.gff")
        move_additional_file("both", f"{tempdir}/outfile-all.gff")
        move_additional_file("report", f"{tempdir}/outfile-report.txt")
    elif script == "agat_sp_manage_UTRs.pl":
        move_additional_file("over_or_equals", "outfile/*_overORequal*.gff")
        move_additional_file("under", "outfile/*_under*.gff")
        move_additional_file("report", "outfile/report.txt")
    elif script == "agat_sp_manage_functional_annotation.pl":
        move_additional_file("gff", f"{tempdir}/outfile/*.gff")
        move_additional_file("cdd", f"{tempdir}/outfile/CDD.txt")
        move_additional_file(
            "duplicates", f"{tempdir}/outfile/duplicatedNameFromBlast.txt"
        )
        move_additional_file("error", f"{tempdir}/outfile/error.txt")
        move_additional_file("go", f"{tempdir}/outfile/GO.txt")
        move_additional_file("interpro", f"{tempdir}/outfile/InterPro.txt")
        move_additional_file("mobi", f"{tempdir}/outfile/MobiDBLite.txt")
        move_additional_file("panther", f"{tempdir}/outfile/PANTHER.txt")
        move_additional_file("report", f"{tempdir}/outfile/report.txt")
        move_additional_file("family", f"{tempdir}/outfile/SUPERFAMILY.txt")
    elif script == "agat_sp_manage_introns.pl":
        move_additional_file("report", f"{tempdir}/outfile/report.txt")
        move_additional_file("cds", f"{tempdir}/outfile/intronPlot_cds.pdf")
        move_additional_file("exon", f"{tempdir}/outfile/intronPlot_exon.pdf")
        move_additional_file(
            "five_prime_utr", f"{tempdir}/outfile/intronPlot_five_prime_utr.pdf"
        )
        move_additional_file(
            "three_prime_utr", f"{tempdir}/outfile/intronPlot_three_prime_utr.pdf"
        )
    elif script == "agat_sp_prokka_fix_fragmented_gene_annotations.pl":
        move_additional_file("report", f"{tempdir}/outfile/report.txt")
        move_additional_file("fasta", f"{tempdir}/outfile/*.fa")
        move_additional_file("gff", f"{tempdir}/outfile/*.gff")
    
