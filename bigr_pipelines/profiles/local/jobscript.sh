#!/bin/sh
# properties = {properties}

set -euo 'pipefail'
shopt -s 'nullglob'

date
hostname
whoami
lshw -json

{exec_job}
