#!/usr/bin/python

__author__ = "Thibault Dayris"
__copyright__ = "Copyright 2021, Thibault Dayris"
__email__ = "thibault.dayris@gustaveroussy.fr"
__license__ = "MIT"

"""
This script contains functions to read and write yaml files from various object 
types in python. Depending on the needs in the upcomming pipelines.
"""

import yaml

from pathlib import Path
from typing import Any, Union


def write_yaml(output_yaml: Union[str, Path], data: dict[str, Any]) -> None:
    """
    Save given dictionnary as Yaml-formatted text file
    """
    if isinstance(str, output_yaml):
        write_yaml_from_str(output_yaml, data)
    else:
        write_yaml_from_path(output_yaml, data)


def write_yaml_from_str(output_yaml: str, data: dict[str, Any]) -> None:
    """
    Save given dictionnary as Yaml-formatted text file, 
    the output_yaml is a string.
    """
    with open(output_yaml, "w") as outyaml:
        yaml.dump(data, outyaml, default_flow_style=False)


def write_yaml_from_path(output_yaml: Path, data: dict[str, Any]) -> None:
    """
    Save given dictionnary as Yaml-formatted text file,
    the output_yaml is a Path
    """
    with output_yaml.open("w") as outyaml:
        yaml.dump(data, outyaml, default_flow_style=False)