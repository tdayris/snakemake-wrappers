#!/usr/bin/bash
set -euo 'pipefail'
shopt -s 'nullglob'

echo "INFO: sourcing functions"

in_container() {
    echo "INFO: Checking wether we're in a container"
    [ -f /.dockerenv ] && return 0
    grep -qsE '(docker|containerd|kubepods|podman|lxc|bwrap)' /proc/1/cgroup 2>/dev/null
    grep -qs 'container' /proc/self/mountinfo 2>/dev/null
}


apt_update_image() {
    echo "INFO: Updating image with apt"
    apt update 
    apt upgrade -y
    apt clean 
    apt autoremove -y 
}

apt_install_apptainer() {
    echo "INFO: Installing apptainer from .deb files"
    apt install -y build-essential wget bash bzip2 ca-certificates gnupg2 squashfs-tools git 
    apt install -y "/opt/apptainer_1.5.4-trixie+_amd64.deb" "/opt/apptainer-suid_1.5.4-trixie+_amd64.deb" 
    apt install --fix-broken 
    apt clean 
    apt autoremove -y 
}

install_conda() {
    echo "INFO: Installing conda from Miniforge shell script"
    bash "/opt/Miniforge3-Linux-x86_64.sh" -bfc 
    eval "$(/root/miniforge3/bin/conda 'shell.bash' 'hook')" 
    echo 'eval "$(/root/miniforge3/bin/conda shell.bash hook)" ; ' >> "/etc/bash.bashrc" 
    conda clean --all --yes 
    . "/etc/bash.bashrc"
}

conda_install_snakemake() {
    local TARGET_ENV_PIN="${1}"
    echo "INFO: Installing ${TARGET_ENV_PIN} environment"
    mkdir --parent --verbose "/root/miniforge3/envs/snakemake"
    conda create -p "/root/miniforge3/envs/snakemake" --file "${TARGET_ENV_PIN}"
    echo "export CONDARC='/etc/conda/condarc'; " >> "/etc/bash.bashrc"
    . "/etc/bash.bashrc"
    conda clean --all --yes 
    # Ensure conda is available in bashrc
    # and run conda init + conda activate on startup
    echo "conda activate /root/miniforge3/envs/snakemake; " >> "/etc/bash.bashrc" 
    conda activate "/root/miniforge3/envs/snakemake"
}

pip_install_snakemake_for_google() {
    echo "INFO: Installing snakemake with pip for google-batch and cloud storage"
    # Install snakemake in '.snakemake/pip-deployments' for persistance 
    python3 -m pip install --root-user-action 'ignore' \
        --target "${PWD}/.snakemake/pip-deployments" \
        'snakemake' \
        'pandas==2.3.3' \
        'snakemake-storage-plugin-gcs==1.1.4' \
        'snakemake-storage-plugin-http==0.3.1' \
        'snakemake-storage-plugin-fs==1.1.3' \
        'snakemake-software-deployment-plugin-conda==0.8.4' \
        'snakemake-software-deployment-plugin-container==0.7.5' \
        'snakemake-executor-plugin-googlebatch==0.5.1' 
    python3 -m pip cache purge 
}

pip_install_snakemake_for_slurm() {
    echo "INFO: Installing snakemake for slurm submissions"
    # Install snakemake in '.snakemake/pip-deployments' for persistance 
    python3 -m pip install --root-user-action 'ignore' \
        --target "${PWD}/.snakemake/pip-deployments" \
        'snakemake' \
        'pandas==2.3.3' \
        'snakemake-storage-plugin-http==0.3.1' \
        'snakemake-storage-plugin-fs==1.1.3' \
        'snakemake-software-deployment-plugin-conda==0.8.4' \
        'snakemake-software-deployment-plugin-container==0.7.5' \
        'snakemake-executor-plugin-slurm==2.8.0'
    python3 -m pip cache purge 
}

pip_install_snakemake_local() {
    echo "INFO: Installing snakemake for lical use"
    # Install snakemake in '.snakemake/pip-deployments' for persistance 
    python3 -m pip install --root-user-action 'ignore' \
        --target "${PWD}/.snakemake/pip-deployments" \
        'snakemake' \
        'pandas==2.3.3' \
        'snakemake-software-deployment-plugin-conda==0.8.4' \
        'snakemake-software-deployment-plugin-container==0.7.5'
    python3 -m pip cache purge 
}


fix_bashrc() {
    echo "INFO: Fixing and sourcing bashrc"
    # Fix PS1 syntax error on Dockers
    sed -i 's|^\[|#\[|' "/etc/bash.bashrc" 

    # Run bashrc to ensure it has no error
    cat "/etc/bash.bashrc" 
    chmod u+x "/etc/bash.bashrc"
    chmod u+x "/opt/init_bashrc.sh"
    . "/etc/bash.bashrc" 
}

add_links() {
    echo "INFO: Adding conda links"
    # Build essential links
    ln -sfv "/root/miniforge3/bin/conda" "/usr/local/bin/conda" 
    ln -sfv "/root/miniforge3/bin/mamba" "/usr/local/bin/mamba" 
}

clean_image() {
    echo "INFO: Cleaning deb and shell"
    # Remove artifacts
    rm -fv "/opt/apptainer_1.5.4-trixie+_amd64.deb" \
        "/opt/apptainer-suid_1.5.4-trixie+_amd64.deb" \
        "/opt/Miniforge3-Linux-x86_64.sh"
}

usage() {
    cat <<EOF
Usage: $(basename "$0") [OPTIONS]

Options:
  -f, --force              Skip pre-run checks
  -m, --mode <MODE>        Execution mode (default: local)
                           MODE is one of: local, google, slurm
  -h, --help               Show this help message and exit
EOF
}

die() {
    echo "ERROR: $*" >&2
    exit 1
}


run_checks() {
    if in_container; then
        echo "INFO: I'm in container"
    else
        die "I might not be in container. Exiting."
    fi
}

FORCE=false
MODE="local"

while [[ $# -gt 0 ]]; do
    case "$1" in
        -f|--force)
            FORCE=true
            shift
            ;;
        -m|--mode)
            [[ $# -lt 2 ]] && die "--mode requires a value (local, google, slurm)"
            MODE="$2"
            shift 2
            ;;
        --mode=*)
            MODE="${1#*=}"
            shift
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        --)
            shift
            ARGS+=("$@")
            break
            ;;
        -*)
            die "Unknown option: $1 (see --help)"
            ;;
        *)
            ARGS+=("$1")
            shift
            ;;
    esac
done

if [[ "$FORCE" == true ]]; then
    echo "INFO: --force given, skipping checks."
else
    run_checks
fi

echo "INFO: Starting (mode: $MODE, force: $FORCE)"

apt_update_image
apt_install_apptainer
install_conda

case "$MODE" in
    local)
        echo "INFO: Running locally..."
        conda_install_snakemake "/opt/environment_local.linux-64.pin.txt"
        # pip_install_snakemake_local
        ;;
    google)
        echo "INFO: Running in google mode..."
        conda_install_snakemake "/opt/environment_google.linux-64.pin.txt"
        # pip_install_snakemake_for_google
        ;;
    slurm)
        echo "INFO: Running in Slurm mode..."
        conda_install_snakemake "/opt/environment_slurm.linux-64.pin.txt"
        # pip_install_snakemake_for_slurm
        ;;
esac

fix_bashrc
add_links
clean_image

echo "INFO: Process over."
