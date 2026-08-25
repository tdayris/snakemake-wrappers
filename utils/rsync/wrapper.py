# coding: utf-8

"""Snakemake wrappers for rsync"""

__author__ = "Thibault Dayris"
__copyright__ = "Copyright 2026, Thibault Dayris"
__license__ = "MIT"

from os.path import isdir, islink
from snakemake.shell import shell
from snakemake_wrapper_utils.snakemake import is_arg
from shlex import quote
from tempfile import TemporaryDirectory

extra = snakemake.params.get("extra", "")
log = snakemake.log_fmt_shell(stdout=False, stderr=True)

# Optional IO files
exclude_from = snakemake.input.get("exclude")
if exclude_from:
    extra += f" --exclude-from={quote(exclude_from)}"

include_from = snakemake.input.get("include")
if include_from:
    extra += f" --include-from={quote(include_from)}"

log_file = snakemake.output.get("log")
if log_file:
    extra += f" --log-file={quote(log_file)}"

password = snakemake.input.get("password")
if password:
    extra += f" --password-file={quote(password)}"

early_input = snakemake.input.get("early_input")
if early_input:
    extra += f" --early-input={quote(early_input)}"

batch_update = snakemake.output.get("batch")
if batch_update:
    extra += f" --write-batch={quote(batch_update)}"

wo_batch_update = snakemake.output.get("wo_batch")
if wo_batch_update:
    extra += f" --only-write-batch={quote(wo_batch_update)}"

r_batch_update = snakemake.output.get("read_batch")
if r_batch_update:
    extra += f" --read-batch={quote(r_batch_update)}"

config = snakemake.input.get("config")
if config:
    extra += f" --config={quote(config)}"

compression = is_arg("-z", extra) or is_arg("--compress-threads", extra)
if compression:
    extra += f" --compress-threads={snakemake.threads - 1}"

files_from = snakemake.input.get("from")
src = snakemake.input.get("src", "")
if files_from and src:
    raise ValueError("Please provide either input.src OR input.from, not both")
elif not (files_from or src):
    raise ValueError("Please provide either input.src OR input.from")
elif files_from:
    extra += f" --files-from={quote(files_from)}"

with TemporaryDirectory() as tempdir:
    shell("rsync {extra} --temp-dir={tempdir:q} {src} {snakemake.output.dst} {log}")
