#!/usr/bin/bash
set -euo 'pipefail'
shopt -s 'nullglob'

build_docker_images() {
    echo -e "\nDOCKER\n"
    echo "Login if needed:"
    echo "docker login ghcr.io -u tdayris"

    echo "Build docker image from dockerfile ${VERSION}-${MODE}"
    docker build --build-arg MODE="${MODE}" --debug --progress 'plain' \
        -t snakemake-conda-apptainer:${VERSION}-${MODE} .
    echo "Define release tag"
    docker tag snakemake-conda-apptainer:${VERSION}-${MODE} \
        ghcr.io/tdayris/snakemake-conda-apptainer:${VERSION}-${MODE}
}

usage() {
    cat <<EOF
Usage $(basename "${0}") [OPTIONS]

Options:
  -u, --update-conda-envs   Use snakedeploy to update conda envs
  -p, --pin-conda-envs      Use snakedeploy to pin conda envs 
  -c, --clear-docker-cache  Clean docker cache to make space
  -t, --test-docker-image   Run the docker image with 'snakemake --help'
  --push-docker-image       Push built image to github
  -v, --version             Image version (default: 0.15)
  -m, --mode                Mode is one of: local, google, slurm (default: local)
  -h, --help                Show this message and exit
EOF
}

die() {
    echo "ERROR: $*" >&2
    exit 1
}

VERSION="0.15"
MODE="local"
UPDATE_CONDA_ENVS=false
PIN_CONDA_ENVS=false
CLEAR_DOCKER_CACHE=false
PUSH_DOCKER_IMAGE=false
TEST_DOCKER_IMAGE=false

while [[ $# -gt 0 ]]; do
    case "$1" in
        -u|--update-conda-envs)
            UPDATE_CONDA_ENVS=true
            shift
            ;;
        -p|--pin-conda-envs)
            PIN_CONDA_ENVS=true
            shift
            ;;
        -c|--clear-docker-cache)
            CLEAR_DOCKER_CACHE=true
            shift
            ;;
        -t|--test-docker-image)
            TEST_DOCKER_IMAGE=true
            shift
            ;;
        --push-docker-image)
            PUSH_DOCKER_IMAGE=true
            shift
            ;;
        -m|--mode)
            [[ $# -lt 2 ]] && die "--mode requires a value (local, batch, slurm)"
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
            die "Unknown option: $1 (see --help)"
            ;;
    esac
done

case "$MODE" in
    local|google|slurm)
        ;;
    *)
        die "Invalid mode '$MODE'. Must be one of: local, google, slurm"
        ;;
esac

if [[ "${CLEAR_DOCKER_CACHE}" == true ]]; then
    echo "INFO cleaning docker cache, images and all."
    docker builder prune --all #--force
fi

if [[ "${UPDATE_CONDA_ENVS}" == true ]]; then
    echo "INFO: Updating conda env for ${MODE}"
    snakedeploy update-conda-envs "environment_${MODE}.yaml"
fi

if [[ "${PIN_CONDA_ENVS}" == true ]]; then
    echo "INFO: Pinning conda env for ${MODE}"
    rm -v "environment_${MODE}.linux-64.pin.txt"
    snakedeploy pin-conda-envs "environment_${MODE}.yaml"
fi

build_docker_images "${MODE}"

if [[ "${PUSH_DOCKER_IMAGE}" == true ]]; then
    echo "INFO: Push ${MODE} image to github"
    docker push "ghcr.io/tdayris/snakemake-conda-apptainer:${VERSION}-${MODE}"
fi

if [[ "${TEST_DOCKER_IMAGE}" == true ]]; then
    echo "INFO: Testing docker image"
    docker run -it "snakemake-conda-apptainer:${VERSION}-${MODE}" snakemake --help
fi

echo "INFO: Process over."
