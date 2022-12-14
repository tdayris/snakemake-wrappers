import datetime
import os
import logging

from snakemake.utils import min_version

##################################
### Initiate logging behaviour ###
##################################
logging.basicConfig(
    filename=f"snakemake.bigr_bulk.{datetime.timestamp(datetime.now())}.log",
    filemode="w", 
    level=logging.DEBUG
)


################################################
### Check versions and environment variables ###
################################################
min_version("7.18")
bigr_bulk_version = "0.0.1"

envvars:
    "SNAKEMAKE_OUTPUT_CACHE"