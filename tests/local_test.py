"""A collection of heavier tests to run locally."""

import shutil
import subprocess
from pathlib import Path

import pytest


@pytest.fixture(scope="module")
def module_path():
    """Parent directory of the project."""
    return Path(__file__).parent.parent



def test_example_config(module_path, tmp_path):
    """The example file should run fine."""

    shutil.copytree(module_path / "workflow", tmp_path / "workflow")
    assert subprocess.run("snakemake")


def test_example_europe(module_path):
    """The European example should run file."""
