# coding: utf-8

"""Snakemake wrappers used to download and deploy librarian sources"""

__author__ = "Thibault Dayris"
__copyright__ = "Copyright 2025, Thibault Dayris"
__license__ = "MIT"

from snakemake.shell import shell
from shlex import quote
from tempfile import TemporaryDirectory

extra = snakemake.params.get("extra", "")
log = snakemake.log_fmt_shell(stdout=True, stderr=True, append=True)


url = snakemake.params.get(
    "url",
    "https://github.com/kartva/Librarian/releases/download/v1.1.1/librarian.tar.gz",
)

shell(
    "( mkdir --parent --verbose {snakemake.output:q} && "
    "  wget --quiet {url:q} && "
    "  tar -xvzf 'librarian.tar.gz' -C {snakemake.output:q} ) {log} "
)
