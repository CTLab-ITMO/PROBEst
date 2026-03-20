# MIT License
#
# Copyright (c) 2025 CTLab-ITMO
#
# Authors: Daniil Smutin, Aleksandr Serdiukov, Vitalii Dravgelis, Artem Ivanov,
# Aleksei Zabashta, Sergey Muravyov, and the CTLab-ITMO university team.
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

"""Integration test: general FISH dataset with JIT FASTA bases, AI + visualization."""

import shutil
import subprocess
import sys
from pathlib import Path

import pytest

PROJECT_ROOT = Path(__file__).resolve().parent.parent


def _require_blast_tools():
    for cmd in ("makeblastdb", "blastn"):
        if shutil.which(cmd) is None:
            pytest.skip(f"{cmd} not on PATH")


@pytest.mark.integration
def test_integration_fish(tmp_path):
    """Same CLI as: …/test.fna, fasta_base true/false dirs, -N 3, --visualize True --AI True.

    True-base genomes must be present in git as ``*.fna.gz`` under ``fasta_base/true_base``
    (plain ``*.fna`` is gitignored); false bases already use gz in the repo.
    """
    _require_blast_tools()

    output = tmp_path / "output"
    cmd = [
        sys.executable,
        str(PROJECT_ROOT / "pipeline.py"),
        "-i",
        str(PROJECT_ROOT / "data/test/general/test.fna"),
        "-o",
        str(output),
        "-tb",
        str(PROJECT_ROOT / "data/test/general/fasta_base/true_base"),
        "-fb",
        str(PROJECT_ROOT / "data/test/general/fasta_base/false_base_1"),
        str(PROJECT_ROOT / "data/test/general/fasta_base/false_base_2"),
        "-a",
        "FISH",
        "--PRIMER_PICK_PRIMER",
        "5",
        "--PRIMER_NUM_RETURN",
        "5",
        "-N",
        "3",
        "--visualize",
        "True",
        "--AI",
        "True",
    ]

    result = subprocess.run(
        cmd,
        cwd=PROJECT_ROOT,
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        pytest.fail(
            f"pipeline.py exited {result.returncode}\n"
            f"--- stdout ---\n{result.stdout}\n--- stderr ---\n{result.stderr}"
        )

    assert (output / "modeling_results.tsv").is_file()
    assert (output / "output_dedegenerated.fa").is_file()
    viz_dir = output / "visualizations"
    assert viz_dir.is_dir()
    assert any(viz_dir.glob("*_visualization.png"))
