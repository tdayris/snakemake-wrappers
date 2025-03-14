__author__ = "Johannes Köster"
__copyright__ = "Copyright 2019, Johannes Köster"
__email__ = "johannes.koester@uni-due.de"
__license__ = "MIT"

import subprocess as sp
import sys
from itertools import product
from snakemake.shell import shell
from snakemake.logging import logger

species = snakemake.params.species.lower()
release = int(snakemake.params.release)
build = snakemake.params.build

branch = ""
if release >= 81 and build == "GRCh37":
    # use the special grch37 branch for new releases
    branch = "grch37/"
elif snakemake.params.get("branch"):
    branch = snakemake.params.branch + "/"

log = snakemake.log_fmt_shell(stdout=False, stderr=True)

spec = ("{build}" if int(release) > 75 else "{build}.{release}").format(
    build=build, release=release
)

suffixes = ""
datatype = snakemake.params.get("datatype", "")
chromosome = snakemake.params.get("chromosome", "")
if datatype == "dna":
    if chromosome:
        suffixes = [f"dna.chromosome.{chrom}.fa.gz" for chrom in chromosome]
    else:
        suffixes = ["dna.primary_assembly.fa.gz", "dna.toplevel.fa.gz"]
elif datatype == "cdna":
    suffixes = ["cdna.all.fa.gz"]
elif datatype == "cds":
    suffixes = ["cds.all.fa.gz"]
elif datatype == "ncrna":
    suffixes = ["ncrna.fa.gz"]
elif datatype == "pep":
    suffixes = ["pep.all.fa.gz"]
else:
    raise ValueError("invalid datatype, must be one of dna, cdna, cds, ncrna, pep")

if chromosome:
    if not datatype == "dna":
        raise ValueError(
            "Invalid datatype. To select individual chromosomes, the datatype must be dna."
        )

url = snakemake.params.get("url", "ftp://ftp.ensembl.org/pub")
spec = spec.format(build=build, release=release)
url_prefix = f"{url}/{branch}release-{release}/fasta/{species}/{datatype}/{species.capitalize()}.{spec}"

for suffix in suffixes:
    success = False
    url = f"{url_prefix}.{suffix}"

    try:
        shell("curl -sSf {url} > /dev/null 2> /dev/null")
    except sp.CalledProcessError:
        if chromosome:
            logger.error(
                f"Unable to download the requested chromosome sequence from Ensembl at: {url_prefix}.{suffix}.",
            )
            break
        else:
            continue

    shell("(curl -L {url} | gzip -d >> {snakemake.output[0]}) {log}")
    success = True
    if not chromosome:
        break

if not success:
    if not chromosome:
        if len(suffixes) > 1:
            url = f"{url_prefix}.[{'|'.join(suffixes)}]"
        else:
            url = f"{url_prefix}.{suffixes[0]}"
        logger.error(
            f"Unable to download the requested reference sequence data from Ensembl at: {url}.",
        )
    logger.error(
        "Please check whether above URL is currently available (might be a temporal server issue). \n"
        "Apart from that, did you check that this combination of species, build, and release is actually provided?",
    )
    exit(1)
