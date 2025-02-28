#!/usr/bin/python3.11
# -*- coding: utf-8 -*-

"""
Snakemake wrapper for bash copy within the IGR's Flamingo cluster
written in actor oriented programming scheme with pykka
"""

__author__ = "Thibault Dayris"
__copyright__ = "Copyright 2022, Thibault Dayris"
__email__ = "thibault.dayris@gustaveroussy.fr"
__license__ = "MIT"

import os
import os.path as op
import pykka
import time

from enum import StrEnum
from tempfile import TemporaryDirectory
from typing import List, Optional, Tuple, Union

# from snakemake.shell import shell

Files = Union[str, List[str]]

NextFilePairs = object()
ManagerType = pykka.ThreadingActor

class CopyStatus(StrEnum):
    TO_COPY = "Not copied yet"
    TO_CONCAT = "Not concatenated yet"
    ONGOING_COPY = "On going copy"
    ONGOING_CONCAT = "On going concatenation"
    DONE = "Done"


class FilePairs(object):
    def __init__(self, source, destination, status) -> None:
        self.source = source
        self.destination = destination
        self.status = status

    def __str__(self) -> str:
        return f"FilePairs[{self.status.value}]: {self.source} -> {self.destination}"


class FileConcat(object):
    def __init__(self, destination, status, *sources) -> None:
        self.sources = sources
        self.destination = destination
        self.status = status

    def __str__(self) -> str:
        return f"FileConcat[{self.status.value} / {'Ready to go' if self.ready() else 'Waiting for data'}] -> {self.destination}"

    def ready(self) -> bool:
        return all(source.status == CopyStatus.DONE for source in self.sources)


class WorkToDo(object):
    def __init__(self) -> None:
        self.copylist = []

    def add(self, files: Union[FilePairs, FileConcat]) -> None:
        self.copylist.append(files)

    def __iter__(self):
        return iter(self.copylist)



class Manager(pykka.ThreadingActor):
    def __init__(self, work) -> None:
        super().__init__()
        self.work_to_do = work
        self.index = -1


    def add(self, data) -> None:
        print(str(data))
        self.work_to_do.add(data)


    def __next__(self):
        self.index += 1
        return self.work_to_do[self.index]


    def __str__(self):
        ratio = [0, 0]
        for data in self.work_to_do:
            ratio[0] += 1 if data.status == CopyStatus.DONE else 0
            ratio[1] += 1
        return f"{ratio[0]} / {ratio[1]} ({(ratio[0]/ratio[1])*100} %)"


    def on_receive(self, message) -> Optional[List[Union[int, FilePairs]]]:
        if isinstance(message, dict):
            if "," in message["source"]:
                # Then it's as list of files to concat
                intermediary_copies = []
                for source in message["source"].split(","):
                    with TemporaryDirectory() as tmpdir:
                        intermediary_copy = FilePairs(
                            source=source, 
                            destination=tmpdir, 
                            status=CopyStatus.TO_COPY
                        )
                        intermediary_copies.append(intermediary_copy)
                        self.add(intermediary_copy)

                self.add(
                    FileConcat(
                        destination=message["destination"], 
                        status=CopyStatus.TO_CONCAT, 
                        *intermediary_copies
                    )
                )

            elif isinstance(message["source"], list):
                # Then there are multiple files that goes to a single directory
                if not op.exists(message["destination"]):
                    os.makedirs(message["destination"], exist_ok=True)
                for source in message["source"]:
                    self.add(
                        FilePairs(
                            source=source,
                            destination=message["destination"],
                            status=CopyStatus.TO_COPY
                        )
                    )
                    
            else:
                # Then it's a simple copy.
                self.add(
                    FilePairs(
                        source=message["source"], 
                        destination=message["destination"], 
                        status=CopyStatus.TO_COPY
                    )
                )

        elif message is NextFilePairs:
            if self.ready:
                return next(self.work_to_do, None)
            else:
                raise ValueError("I'm not ready to work.")



class Copy(pykka.ThreadingActor):
    def __init__(
        self, 
        manager: ManagerType,
        log: str,
        cold_storage: Tuple[str] = ("/mnt/isilon", "/mnt/archivage", "/mnt/install"), 
        irods_prefix: Tuple[str] = ("/odin/kdi/dataset/"),
        iget_template: str="iget -vKr {} {} >> {} 2>&1", 
        rsync_template: str="rsync -cvrhP {} {} >> {} 2>&1", 
        symlink_template: str="ln -sfr {} {} >> {} 2>&1",
        cat_template: str="cat {} > {} 2>> {}"
    ) -> None:
        super().__init__()
        self.manager = manager
        self.iget_template = iget_template
        self.rsync_template = rsync_template
        self.symlink_template = symlink_template
        self.cold_storage = cold_storage
        self.irods_prefix = irods_prefix
        self.cat_template = cat_template
        self.log = log


    def copy(self) -> None:
        file_to_copy = self.manager.ask(NextFilePairs)
        if file_to_copy is None:
            time.sleep(60 * 5)
        elif file_to_copy is FilePairs:
            file_to_copy.status = CopyStatus.ONGOING_COPY
            if file_to_copy.source.startswith(self.irods_prefix):
                shell(self.iget_template.format(
                    file_to_copy.source, file_to_copy.destination, self.log
                ))
            elif file_to_copy.source.startswith(self.cold_storage):
                shell(self.rsync_template.format(
                    file_to_copy.source, file_to_copy.destination, self.log
                ))
            else:
                shell(self.symlink_template.format(
                    file_to_copy.source, file_to_copy.destination, self.log
                ))
        elif file_to_copy is FileConcat:
            file_to_copy.status = CopyStatus.ONGOING_CONCAT
            if file_to_copy.ready():
                shell(self.cat_template.format(
                    file_to_copy.source, file_to_copy.destination, self.log
                ))
        else:
            raise ValueError(f"I don't understand: {file_to_copy}")

        file_to_copy.status = CopyStatus.DONE


def shell(string: str) -> None:
    print(f"shell: {string}")