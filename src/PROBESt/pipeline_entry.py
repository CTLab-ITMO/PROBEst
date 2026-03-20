"""Console entry for the main PROBESt pipeline (``pipeline.py`` at repository root)."""

import os
import runpy
import sys


def main() -> None:
    try:
        import pipeline as _pipeline
    except ModuleNotFoundError:
        _pipeline = None

    if _pipeline is not None:
        _pipeline.main()
        return

    here = os.path.dirname(os.path.abspath(__file__))
    root = os.path.normpath(os.path.join(here, os.pardir, os.pardir))
    pipeline_py = os.path.join(root, "pipeline.py")
    if os.path.isfile(pipeline_py):
        runpy.run_path(pipeline_py, run_name="__main__")
        return

    print(
        "PROBESt: could not load the pipeline (install the package with pip or "
        "run from a source tree containing pipeline.py).",
        file=sys.stderr,
    )
    raise SystemExit(1)
