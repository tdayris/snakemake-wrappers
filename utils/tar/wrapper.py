# coding: utf-8

"""Snakemake wrapper for tar"""

__author__ = "Thibault Dayris"
__copyright__ = "Copyright 2025, Thibault Dayris"
__license__ = "MIT"

from getpass import getuser
from grp import getgrgid
from os import getuid, getgid
from os.path import splitext
from snakemake.shell import shell
from snakemake.io import AnnotatedString
from snakemake.ioutils import subpath
from snakemake_wrapper_utils.snakemake import get_format
from shlex import quote

extra = snakemake.params.get("extra", "")
log = snakemake.log_fmt_shell(stdout=True, stderr=True, append=True)

if snakemake.params.get("explicit_owner"):
    user_i = getuid()
    user_n = getuser()
    extra += f" --owner='{user_n}:{user_i}'"

if snakemake.params.get("explicit_group"):
    grp_i = getgid()
    grp_n = getgrgid(grp_i).gr_name
    extra += f" --group='{grp_n}:{grp_i}'"

archive = snakemake.output.get("archive") or snakemake.input.get("archive")
if archive:
    extra += f" --file={quote(str(archive))}"

    try:
        out_format = get_format(
            path=archive,
            ignore_compression=False,
        )
        compress_fmt = ("bzip2", "gzip", "lzma", "zstd", "lzop")
        if out_format in compress_fmt:
            extra += f" --{out_format} --no-auto-compress"
    except ValueError:
        pass

incremental_backup = snakemake.input.get("incremental_backup")
if incremental_backup:
    extra += f" --listed-incremental={quote(str(incremental_backup))}"

files_from = snakemake.input.get("files_list")
if files_from:
    extra += f" --files-from={quote(str(files_from))}"

exclude_from = snakemake.input.get("exclusion_list")
if exclude_from:
    extra += f" --exclude-from={quote(str(exclude_from))}"

infiles = snakemake.input.get("data", "")
outfiles = snakemake.output.get("data", "")
outdir = snakemake.output.get("outdir", "")
if isinstance(outfiles, list):
    for o in outfiles:
        parent = subpath(str(o), parent=True)
        shell("mkdir --parent --verbose {parent:q} {log}")
elif outdir:
    shell("mkdir --parent --verbose {outdir:q} {log}")
    outfiles = f" --directory={quote(outdir)}"
else:
    parent = subpath(str(outfiles), parent=True)
    shell("mkdir --parent --verbose {parent:q} {log}")




shell("tar {extra} {infiles} {outfiles} {log}")
