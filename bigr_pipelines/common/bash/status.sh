#!/usr/bin/env bash

# This script looks for a given snakemake repository and searches for the
# current pipeline status (returns a string and an error code)
# Returns:
# ERROR with exit code 1
# ON_GOING with exit code 0
# DONE with exit code 0
# UNKNOWN with exit code 2
# DIR_NOT_FOUND with exit code 3

# Change

if [ ! -d "${1}" ] ; then
    echo "DIR_NOT_FOUND"
    # exit 3
fi

cd "${1}" || exit 3

if [ -f "DONE" ] ; then
    echo "DONE"
    # exit 0
fi

if [ -f "ERROR" ] ; then
    echo "ERROR"
    # exit 1
fi

if [ -f "ON_GOING" ] ; then
    echo "ON_GOING"
    # exit 0
fi

if [[ $(grep -iP "(Wildcard|Workflow|Syntax)Error" *.log --quiet) -eq 0 ]]; then
    echo "ERROR"
    # exit 1
fi

echo "UNKNOWN"
# exit 2

exit 0