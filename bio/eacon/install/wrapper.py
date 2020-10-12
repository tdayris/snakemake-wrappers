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
    shell(f'wget "https://partage.gustaveroussy.fr/pydio_public/b88fb8?dl=true&file=/OncoScan.na33.r4_0.1.0.tar.gz" -O {out_dir}/OncoScan.na33.r4_0.1.0.tar.gz {{log}}')
    shell(f'R  --vanilla CMD INSTALL --clean --data-compress=gzip {out_dir}/OncoScan.na33.r4_0.1.0.tar.gz {log}')

    shell(f'wget "https://partage.gustaveroussy.fr/pydio_public/582a03?dl=true&file=/OncoScan.na36.r1_0.1.0.tar.gz" -O {out_dir}/OncoScan.na36.r1_0.1.0.tar.gz {{log}}')
    shell(f'R  --vanilla CMD INSTALL --clean --data-compress=gzip {out_dir}/OncoScan.na36.r1_0.1.0.tar.gz {log}')

if snakemake.params.get("OncoScanCNV", True) is True:
    shell(f'wget "https://partage.gustaveroussy.fr/pydio_public/cd59c8?dl=true&file=/OncoScanCNV.na33.r2_0.1.0.tar.gz" -O {out_dir}/OncoScanCNV.na33.r2_0.1.0.tar.gz {{log}}')
    shell(f'R  --vanilla CMD INSTALL --clean --data-compress=gzip {out_dir}/OncoScanCNV.na33.r2_0.1.0.tar.gz {log}')

    shell(f'wget "https://partage.gustaveroussy.fr/pydio_public/41d8af?dl=true&file=/OncoScanCNV.na36.r1_0.1.0.tar.gz" -O {out_dir}/OncoScanCNV.na36.r1_0.1.0.tar.gz {{log}}')
    shell(f'R  --vanilla CMD INSTALL --clean --data-compress=gzip {out_dir}/OncoScanCNV.na36.r1_0.1.0.tar.gz {log}')


if snakemake.params.get("CytoScan750K", True) is True:
    shell(f'wget "https://partage.gustaveroussy.fr/pydio_public/74d4cf?dl=true&file=/CytoScan750K.Array.na33.r4_0.1.0.tar.gz" -O {out_dir}/CytoScan750K.Array.na33.r4_0.1.0.tar.gz {{log}}')
    shell(f'R  --vanilla CMD INSTALL --clean --data-compress=gzip {out_dir}/CytoScan750K.Array.na33.r4_0.1.0.tar.gz {log}')

    shell(f'wget "https://partage.gustaveroussy.fr/pydio_public/656d13?dl=true&file=/CytoScan750K.Array.na36.r1_0.1.0.tar.gz" -O {out_dir}/CytoScan750K.Array.na36.r1_0.1.0.tar.gz {{log}}')
    shell(f'R  --vanilla CMD INSTALL --clean --data-compress=gzip {out_dir}/CytoScan750K.Array.na36.r1_0.1.0.tar.gz {log}')

if snakemake.params.get("CytoScanHD", True) is True:
    shell(f'wget "https://partage.gustaveroussy.fr/pydio_public/bc4b54?dl=true&file=/CytoScanHD.Array.na33.r4_0.1.0.tar.gz" -O {out_dir}/CytoScanHD.Array.na33.r4_0.1.0.tar.gz {{log}}')
    shell(f'R  --vanilla CMD INSTALL --clean --data-compress=gzip {out_dir}/CytoScanHD.Array.na33.r4_0.1.0.tar.gz {log}')

    shell(f'wget "https://partage.gustaveroussy.fr/pydio_public/24b026?dl=true&file=/CytoScanHD.Array.na36.r1_0.1.0.tar.gz" -O {out_dir}/CytoScanHD.Array.na36.r1_0.1.0.tar.gz {{log}}')
    shell(f'R  --vanilla CMD INSTALL --clean --data-compress=gzip {out_dir}/CytoScanHD.Array.na36.r1_0.1.0.tar.gz {log}')

if snakemake.params.get("norm", True) is True:
    shell(f'wget "https://partage.gustaveroussy.fr/pydio_public/e6fe22?dl=true&file=/rcnorm_0.1.5.tar.gz" -O {out_dir}/rcnorm_0.1.5.tar.gz {{log}}')
    shell(f'R  --vanilla CMD INSTALL --clean --data-compress=gzip {out_dir}/rcnorm_0.1.5.tar.gz {log}')

    shell(f'wget "https://partage.gustaveroussy.fr/pydio_public/083305?dl=true&file=/affy.CN.norm.data_0.1.2.tar.gz" -O {out_dir}/affy.CN.norm.data_0.1.2.tar.gz {{log}}')
    shell(f'R --vanilla CMD INSTALL --clean --data-compress=gzip {out_dir}/affy.CN.norm.data_0.1.2.tar.gz {log}')

if snakemake.params.get("EaCoN_dev", False) is True:
    shell('R --vanilla -e \'options(unzip = "internal"); Sys.setenv(TAR = \"$(which tar)\"); devtools::install_github("gustaveroussy/EaCoN", force=TRUE); \' {log}')
