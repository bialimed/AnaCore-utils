#!/usr/bin/env python
import os
import re
from setuptools import setup


def get_version():
    version = None
    notes_filepath = "RELEASES_NOTES.md"
    if os.path.exists(notes_filepath):
        with open(notes_filepath) as FH:
            first_line = FH.readline()
            version = re.search(r"^\#\s+.+\s+(.+)\s+\[", first_line).groups()[0]  # Example: "# v2.5.0 [DEV]"
    return version


def load_scripts(path):
    scripts = []
    for filename in os.listdir(path):
        filepath = os.path.join(path, filename)
        if os.path.isdir(filepath):
            if filename != "test":
                load_scripts(filepath)
        else:
            if filename.endswith(".py") and not filename.startswith("__"):
                scripts.append(filepath)
    print(scripts)
    return scripts


setup(
    version=get_version(),
    scripts=load_scripts("bin")
)
