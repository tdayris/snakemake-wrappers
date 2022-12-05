"""
This Snakefile deals with IO from iRODS
"""

# Copy files on BiGR Flamingo
# iRODS paths are accepted
"""
001.bigr_copy
from
-> Entry job
by
-> End job
"""


rule bigr_copy:
    output:
        "{output_file}",
    threads: 1
    resources:
        mem_mb=get_1gb_per_attempt,
        time_min=get_45min_per_attempt,
    retries: 2
    params:
        input=lambda wildcards: io_dict[wildcards.output_file],
        irods_extra=irods_extra
    log:
        "logs/001.bigr_copy/{output_file}.log",
    wrapper:
        "bio/BiGR/copy"