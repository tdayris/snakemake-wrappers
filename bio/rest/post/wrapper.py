#!/usr/bin/env python3
# coding: utf-8

"""Snakemake wrapper to post rest requests"""

import json
import requests


files = {str(k), str(v) for k, v in snakemake.input.items()}
extra = snakemake.params.get("extra", {})
address = snakemake.params["address"]

request = requests.post(address, files=files, **extra)
request_result = request.json()

with open(snakemake.output["json"], "w") as request_json:
    json.dump(request_result, request_json)