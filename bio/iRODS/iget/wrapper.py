#!/usr/bin/python3.8
# conding: utf-8

"""
A script that handles iRODS copy.

WARNING: Input file has to be defined in params, as it does not exists outside
iRODS environment. Therefore, snakemake would raise a FileNotFoundError
(false positive in best cases).

snakemake.RemoteProvider.iRODS is not used since authentication cannot be
securely transmitted within the wrapper.
"""

__author__ = "Thibault Dayris"
__copyright__ = "Copyright 2020, Thibault Dayris"
__email__ = "thibault.dayris@gustaveroussy.fr"
__license__ = "MIT"


from snakemake.shell import shell

# Prepare logging
log = snakemake.log_fmt_shell(stdout=True, stderr=True)
extra = snakemake.params.get("extra", "-VK")

if not isinstance(snakemake.params['input'], str):
    raise TypeError(
        "A single file is expected in this wrapper. Multiple copies should be "
        "performed in multiple call of this pipeline."
    )

shell(
    " iget "  # Copy command from iRODS
    " {extra} "  # Extra parameters
    " -N {snakemake.threads} "  # Maximum number of threads
    " {snakemake.params['input']} "  # Input collection
    " {snakemake.output[0]} "  # Output path
    " {log} "  # Logging behavior
)
