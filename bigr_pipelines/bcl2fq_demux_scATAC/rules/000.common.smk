
#####################################
# Check Snakemake version           #
# Import search-and-buil functions  #
# for config and design             #
#####################################

# Official libraries
import os
import functools

from snakemake.utils import min_version
from pathlib import Path
from yaml import dump
from typing import Any, Dict, List

min_version("7.5")

import sys

workflow_source_dir = Path(snakemake.workflow.srcdir(".."))
common = str(workflow_source_dir / ".." / "common" / "python")
sys.path.append(common)h-and-buil functions  #


from snakemake.utils import min_version
from pathlib import Path
from yaml import dump
from typing import Any, Dict, List

