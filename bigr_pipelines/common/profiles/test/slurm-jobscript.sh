#!/bin/bash
# properties = {properties}

hostname
date

conda info --envs

/usr/bin/time -v -p bash -c "{exec_job}"
