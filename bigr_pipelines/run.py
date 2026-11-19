# coding: utf-8

"""This script helps user with their pipeline launch"""

import argparse
import os
import os.path
import pathlib
import rich.console
import rich.syntax
import rich_argparse
import shlex
import shutil
import sys
import typing

from loguru import logger
from snakemake.shell import shell

def get_pipelines(directory: pathlib.Path) -> typing.Generator[str, None, None]:
    """Return a set of available pipelines"""
    for child in directory.iterdir():
        snakefile: pathlib.Path = child / "Snakefile"
        if child.is_dir() and str(child.name) != "profiles" and snakefile.exists():
            yield str(child.name)

def get_profiles(directory: pathlib.Path) -> typing.Generator[str, None, None]:
    """Return a set of available profiles"""
    for child in directory.iterdir():
        expected_config = child / "config.yaml"
        if child.is_dir() and expected_config.exists():
            yield str(child.name)


def main() -> None:
    launcher_path: pathlib.PosixPath = pathlib.Path(os.path.realpath(__file__))
    pipelines: tuple[str] = tuple(get_pipelines(launcher_path.parent))
    profiles: tuple[str] = tuple(get_profiles(launcher_path.parent / "profiles"))
    working_dir: pathlib.PosixPath = pathlib.Path(".")

    parser = argparse.ArgumentParser(
        formatter_class=rich_argparse.ArgumentDefaultsRichHelpFormatter,
        description=__doc__,
    )
    parser.add_argument(
        "-n", "--pipeline-name",
        type=str,
        help="Name of a pipeline to run",
        choices=pipelines,
        default=pipelines[0],
    )
    parser.add_argument(
        "-p", "--profile",
        type=str,
        help="Name of a profile to use",
        choices=profiles,
        default="local",
    )
    parser.add_argument(
        "-c", "--config-file",
        type=str,
        default="config/config.yaml",
        help="Path to configuration file",
    )
    parser.add_argument(
        "-d", "--dry-run",
        action="store_true",
        default=False,
        help="Don't run anything, just show final command line and check for input",
    )
    parser.add_argument(
        "-r", "--repo",
        action="store_true",
        default=False,
        help="Add path to this directory as a key to the configuration",
    )
    parser.add_argument(
        "-v", "--verbosity",
        type=str,
        choices={"debug", "info"},
        default="info",
        help="Logging level",
    )
    args = parser.parse_args()
    logger.remove()
    logger.add(sys.stderr, level=args.verbosity.upper())
    logger.debug("Command line arguments parsed")
    extra: str = ""

    snakefile: pathlib.Path = launcher_path.parent / args.pipeline_name / "Snakefile"
    if snakefile.exists():
        extra += f" --snakefile '{snakefile.resolve()}'"
    else:
        raise FileNotFoundError(f"Could not find {snakefile=}")

    profile: pathlib.Path = launcher_path.parent / "profiles" / args.profile
    workflow_profile: pathlib.Path = snakefile.parent / "workflow-profiles" / args.profile
    if profile.exists():
        extra += f" --profile '{profile.resolve()}'"
    else:
        raise FileNotFoundError(f"Could not find {profile=}")

    if workflow_profile.exists():
        extra += f" --workflow-profile '{workflow_profile.resolve()}'"
    else:
        logger.warning(
            f"Could not find workflow profile {workflow_profile=}. This could "
            "be normal and expected."
        )

    if str(args.profile) in ("flamingo", "slurm"):
        partitions: pathlib.Path = profile / "partitions.yaml"
        if partitions.exists():
            extra += f" --slurm-partition-config '{partitions.resolve()}'"


    config: pathlib.PosixPath = pathlib.Path(args.config_file)
    if config.exists():
        extra += f" --configfile '{config.resolve()}'"
    else:
        raise FileNotFoundError(f"Could not find {config=}")
    logger.info("All required files found")

    irods_home: str = pathlib.Path(os.path.expanduser("~/.irods/irods_environment.json"))
    if irods_home.exists():
        extra += f" --storage-irods-home '{irods_home.resolve()}'"
    else:
        logger.warning("No irods default environment found.")

    if args.repo:
        repo_dir = launcher_path.parent.parent
        extra += f" --config repo='{repo_dir.resolve()}'"


    cmd: str = str(f"snakemake {extra} ")
    logger.debug("Command line built, running Snakemake:")
    rich.console.Console().print(rich.syntax.Syntax(cmd.replace("--", "\\\n  --"), "bash", line_numbers=False, word_wrap=True,))
    if not args.dry_run:
        shell(cmd)
    else:
        print(cmd)
        logger.debug("Dry-run mode. Stopping here.")


if __name__ == "__main__":
    main()
