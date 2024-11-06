# coding: utf-8

"""Snakemake wrapper designed to make syntethic tests datasets for Sopa"""

__author__ = "Dayris Thibault"
__copyright__ = "Copyright 2024, Dayris Thibault"
__email__ = "thibault.dayris@gustaveroussy.fr"
__license__ = "MIT"

from sopa.utils.data import blobs, uniform

# Generate test dataset
if snakemake.params.use_blob:
    dataset = blobs(**snakemake.params.get("extra", {}))
else:
    dataset = uniform(**snakemake.params.get("extra", {}))


# Save results on disk
dataset.write(
    file_path=snakemake.output[0],
    **snakemake.params.get("write", {}),
)

