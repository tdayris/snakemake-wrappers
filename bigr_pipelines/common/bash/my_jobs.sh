#!/usr/bin/env bash
set -e

squeue_format='%.10i %.9P %.35j %.8u %.8T %.10M %.10l %.6D %.3C %.10m %R'
sed_colors='s/RUNNING/\033[1;36mRUNNING\033[0m/g;s/PENDING/\033[1;33mPENDING\033[0m/g;s/COMPLETED/\033[1;32mCOMPLETED\033[0m/g;'

squeue -o "${squeue_format}" -S j -u "${USER}" | sed "${sed_colors}"