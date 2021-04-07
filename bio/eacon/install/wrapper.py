#!/bin/python3.8

# This wrapper downloads and install optional EaCoN packages

__author__ = "Thibault Dayris"
__copyright__ = "Copyright 2020, Thibault Dayris"
__email__ = "thibault.dayris@gustaveroussy.fr"
__license__ = "MIT"

from snakemake.utils import makedirs
from snakemake.shell import shell
log = snakemake.log_fmt_shell(stdout=True, stderr=True, append=True)

out_dir = snakemake.output[0]
makedirs(out_dir)

shell('R --vanilla -e "options(unzip = \'internal\'); Sys.setenv(TAR = \'$(which tar)\'); devtools::install_github(\'gustaveroussy/apt.oncoscan.2.4.0\'); " {log}')
shell('R --vanilla -e \'options(unzip = "internal"); Sys.setenv(TAR = "$(which tar)"); devtools::install_github("gustaveroussy/apt.cytoscan.2.4.0"); \' {log}')

if snakemake.params.get("OncoScan", True) is True:
    shell('wget "https://nextcloud.gustaveroussy.fr/s/jrWWJDdpJbaR82o/download" -O {out_dir}/OncoScan.na33.r4_0.1.0.tar.gz {log}')
    shell('R  --vanilla CMD INSTALL --clean --data-compress=gzip {out_dir}/OncoScan.na33.r4_0.1.0.tar.gz {log}')

    shell('wget "https://nextcloud.gustaveroussy.fr/s/ZpnYYwmKPzaeHWj/download" -O {out_dir}/OncoScan.na36.r1_0.1.0.tar.gz {log}')
    shell('R  --vanilla CMD INSTALL --clean --data-compress=gzip {out_dir}/OncoScan.na36.r1_0.1.0.tar.gz {log}')

if snakemake.params.get("OncoScanCNV", True) is True:
    shell('wget "https://nextcloud.gustaveroussy.fr/s/BanRj6fAn4HFAA5/download" -O {out_dir}/OncoScanCNV.na33.r2_0.1.0.tar.gz {log}')
    shell('R  --vanilla CMD INSTALL --clean --data-compress=gzip {out_dir}/OncoScanCNV.na33.r2_0.1.0.tar.gz {log}')

    shell('wget "https://nextcloud.gustaveroussy.fr/s/MQ9LwiZAHxnzJ2D/download" -O {out_dir}/OncoScanCNV.na36.r1_0.1.0.tar.gz {log}')
    shell('R  --vanilla CMD INSTALL --clean --data-compress=gzip {out_dir}/OncoScanCNV.na36.r1_0.1.0.tar.gz {log}')


if snakemake.params.get("CytoScan750K", True) is True:
    shell('wget "https://nextcloud.gustaveroussy.fr/s/zx6iwPNKxX798Zg/download" -O {out_dir}/CytoScan750K.Array.na33.r4_0.1.0.tar.gz {log}')
    shell('R  --vanilla CMD INSTALL --clean --data-compress=gzip {out_dir}/CytoScan750K.Array.na33.r4_0.1.0.tar.gz {log}')

    shell('wget "https://nextcloud.gustaveroussy.fr/s/riXGCQNBENdkQSM/download" -O {out_dir}/CytoScan750K.Array.na36.r1_0.1.0.tar.gz {log}')
    shell('R  --vanilla CMD INSTALL --clean --data-compress=gzip {out_dir}/CytoScan750K.Array.na36.r1_0.1.0.tar.gz {log}')

if snakemake.params.get("CytoScanHD", True) is True:
    shell('wget "https://nextcloud.gustaveroussy.fr/s/ip294gysJdcccYm/download" -O {out_dir}/CytoScanHD.Array.na33.r4_0.1.0.tar.gz {log}')
    shell('R  --vanilla CMD INSTALL --clean --data-compress=gzip {out_dir}/CytoScanHD.Array.na33.r4_0.1.0.tar.gz {log}')

    shell('wget "https://nextcloud.gustaveroussy.fr/s/SjRmBFreAee9mqD/download" -O {out_dir}/CytoScanHD.Array.na36.r1_0.1.0.tar.gz {log}')
    shell('R  --vanilla CMD INSTALL --clean --data-compress=gzip {out_dir}/CytoScanHD.Array.na36.r1_0.1.0.tar.gz {log}')

if snakemake.params.get("genomewide", False) is True:
    shell('R --vanilla CMD \'devtools::install_github("gustaveroussy/apt.snp6.1.20.0");\' {log}')
    shell('wget "https://nextcloud.gustaveroussy.fr/s/46iyPjPPjFsni5S/download" -O {out_dir}/GenomeWideSNP.6.na35.r1_0.1.0.tar.gz {log}')
    shell('R  --vanilla CMD INSTALL --clean --data-compress=gzip {out_dir}/GenomeWideSNP.6.na35.r1_0.1.0.tar.gz {log}')

if snakemake.params.get("norm", True) is True:
    shell('wget "https://nextcloud.gustaveroussy.fr/s/Zc7JR3QaAk6rFBi/download" -O {out_dir}/rcnorm_0.1.5.tar.gz {log}')
    shell('R  --vanilla CMD INSTALL --clean --data-compress=gzip {out_dir}/rcnorm_0.1.5.tar.gz {log}')

    shell('wget "https://nextcloud.gustaveroussy.fr/s/gfLnN8xrndFHMdM/download" -O {out_dir}/affy.CN.norm.data_0.1.2.tar.gz {log}')
    shell('R --vanilla CMD INSTALL --clean --data-compress=gzip {out_dir}/affy.CN.norm.data_0.1.2.tar.gz {log}')

if snakemake.params.get("EaCoN_dev", False) is True:
    shell('R --vanilla -e \'options(unzip = "internal"); Sys.setenv(TAR = \"$(which tar)\"); devtools::install_github("gustaveroussy/EaCoN", force=TRUE); \' {log}')

if snakemake.params.get("EaCoN_chromosomes", False) is True:
    shell('R --vanilla -e \'options(unzip = "internal"); Sys.setenv(TAR = \"$(which tar)\"); devtools::install_github("aoumess/chromosomes", force=TRUE);\' {log}')
