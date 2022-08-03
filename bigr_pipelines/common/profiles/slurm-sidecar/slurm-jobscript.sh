#!/bin/bash
# properties = {properties}

hostname
date

/usr/bin/time -v -p bash -c "{exec_job}"
