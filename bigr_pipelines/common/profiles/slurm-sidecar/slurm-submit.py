#!/usr/bin/env python3
"""
Snakemake SLURM submit script.
"""
import json
import logging
import os

import requests
from snakemake.utils import read_job_properties

import slurm_utils

logger = logging.getLogger(__name__)

SIDECAR_VARS = os.environ.get("SNAKEMAKE_CLUSTER_SIDECAR_VARS", None)
DEBUG = bool(int(os.environ.get("SNAKEMAKE_SLURM_DEBUG", "0")))

if DEBUG:
    logging.basicConfig(level=logging.DEBUG)
    logger.setLevel(logging.DEBUG)


def register_with_sidecar(jobid):
    if SIDECAR_VARS is None:
        return
    sidecar_vars = json.loads(SIDECAR_VARS)
    url = "http://localhost:%d/job/register/%s" % (sidecar_vars["server_port"], jobid)
    logger.debug("POST to %s", url)
    headers = {"Authorization": "Bearer %s" % sidecar_vars["server_secret"]}
    requests.post(url, headers=headers)


# cookiecutter arguments
# SBATCH_DEFAULTS = "--account ${USER}"
CLUSTER = f"--user={os.environ.get('USER')}"
# CLUSTER_CONFIG = "flamingo.yaml"

RESOURCE_MAPPING = {
    "time": ("time", "runtime", "walltime", "time_min"),
    "mem": ("mem", "mem_mb", "ram", "memory"),
    "mem-per-cpu": ("mem-per-cpu", "mem_per_cpu", "mem_per_thread"),
    "nodes": ("nodes", "nnodes"),
    "partition": ("partition", "queue"),
}

# parse job
jobscript = slurm_utils.parse_jobscript()
job_properties = read_job_properties(jobscript)

sbatch_options = {}
# cluster_config = slurm_utils.load_cluster_config(CLUSTER_CONFIG)

# 1) sbatch default arguments and cluster
# Ignore sbatch defaults
# sbatch_options.update(slurm_utils.parse_sbatch_defaults(SBATCH_DEFAULTS))
# Ignore cluster config
# sbatch_options.update(slurm_utils.parse_sbatch_defaults(CLUSTER))

# 2) cluster_config defaults
# Ignore cluster config
# sbatch_options.update(cluster_config["__default__"])

# 3) Convert resources (no unit conversion!) and threads
# Assume all submissions are already in proper unit
# sbatch_options.update(slurm_utils.convert_job_properties(job_properties, RESOURCE_MAPPING))

# 4) cluster_config for particular rule
# Ignore cluster config
# sbatch_options.update(cluster_config.get(job_properties.get("rule"), {}))

# 5) cluster_config options
# Ignore cluster config
# sbatch_options.update(job_properties.get("cluster", {}))

# 6) Format pattern in snakemake style
sbatch_options = slurm_utils.format_values(sbatch_options, job_properties)
if "partition" not in sbatch_options.keys():
    sbatch_options["partition"] = slurm_utils.set_partition(sbatch_options.get("time_min", 0))

if not "output" in sbatch_options.keys():
    sbatch_options["output"] = "logs/slurm/slurm-%x-%j-%N.out"

if not "error" in sbatch_options.keys():
    sbatch_options["error"] = "logs/slurm/slurm-%x-%j-%N.err"

if not "mail_type" in sbatch_options.keys():
    sbatch_options["mail_type"] = "END,FAIL"
    sbatch_options["mail_user"] = "thibault.dayris@gustaveroussy.fr"

# ensure sbatch output dirs exist
for o in ("output", "error"):
    slurm_utils.ensure_dirs_exist(sbatch_options[o]) if o in sbatch_options else None

# submit job and echo id back to Snakemake (must be the only stdout)
jobid = slurm_utils.submit_job(jobscript, **sbatch_options)
logger.debug("Registering %s with sidecar...", jobid)
register_with_sidecar(jobid)
logger.debug("... done registering with sidecar")
print(jobid)
