import os
import sys

sys.path.insert(0, os.path.dirname(__file__))

import _git as _git
import _version
from hatchling.metadata.plugin.interface import MetadataHookInterface


def _get_subpackages(directory):
    subpackages = []
    for entry in os.scandir(directory):
        if entry.is_dir() and entry.name.startswith("wool"):
            subpackages.append(entry.name)
    return subpackages


class MetadataHook(MetadataHookInterface):
    PLUGIN_NAME = "wool-metadata"

    def update(self, metadata):
        version = _version.PythonicVersion.parse.git()
        metadata["version"] = str(version)
