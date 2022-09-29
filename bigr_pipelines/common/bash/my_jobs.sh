#!/usr/bin/env bash
set -e

squeue_format='%.10i %.9P %.35j %.8u %.8T %.10M %.10l %.6D %.3C %.10m %R'
sed_colors='s/RUNNING/\x1b[1;36mRUNNING\x1b[0m/g;s/PENDING/\x1b[1;33mPENDING\x1b[0m/g;s/COMPLETED/\x1b[1;32mCOMPLETED\x1b[0m/g;s/shortq/\x1b[1;36mshortq\x1b[0m/g;s/mediumq/\x1b[1;33mmediumq\x1b[0m/g;s/longq/\x1b[1;33mlongq\x1b[0m/g;s/verylongq/\x1b[41mverylongq\x1b[0m/g'

squeue -o "${squeue_format}" -S j -u "${USER}" | sed "${sed_colors}"