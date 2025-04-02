"""A collection of heavier tests to run locally.

IMPORTANT: it's recommended to avoid running all tests at the same time!
You'll likely run out of memory.
"""

import shutil
import subprocess
from pathlib import Path

import pytest


def copy_workflow(workflow_path, dest_path):
    """Copy workflow-relevant folders to a temporary location."""
    for folder in ["workflow", "config"]:
        shutil.copytree(workflow_path / folder, dest_path / folder)


@pytest.fixture(scope="module")
def module_path():
    """Parent directory of the project."""
    return Path(__file__).parent.parent


def test_config_example(module_path, tmp_path):
    """The example file should result in a successful run."""
    copy_workflow(module_path, tmp_path)
    result_file = "results/shapes.parquet"
    config_file = Path(tmp_path / "config/config.yaml")
    subprocess.run(
        f"snakemake --cores 3 --configfile={config_file} {result_file}",
        shell=True,
        check=True,
        cwd=tmp_path,
    )
    assert Path(tmp_path / result_file).exists()

def test_europe_example(module_path, tmp_path):
    """Run a heavy workflow building shapes for European models."""
    copy_workflow(module_path, tmp_path)
    result_file = "results/shapes.parquet"
    config_file = Path(tmp_path / "config/europe_example.yaml")
    subprocess.run(
        f"snakemake --cores 3 --configfile={config_file} {result_file}",
        shell=True,
        check=True,
        cwd=tmp_path,
    )
    assert Path(tmp_path / result_file).exists()
