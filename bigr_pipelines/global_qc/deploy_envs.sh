# this script allow to download singularity images containing CellRanger for the global_qc pipeline
# using: bash deploy_envs.sh

cd $(dirname $0)
wget https://zenodo.org/records/10652954/files/cellranger_arc_v2.0.2.simg?download=1 -O ./envs/cellranger_arc_v2.0.2.simg
wget https://zenodo.org/records/10652954/files/cellranger_atac_v2.1.0.simg?download=1 -O ./envs/cellranger_atac_v2.1.0.simg
wget https://zenodo.org/records/10652954/files/cellranger_v7.2.0.simg?download=1 -O ./envs/cellranger_v7.2.0.simg