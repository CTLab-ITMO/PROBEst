# Extraction submodule

This folder contains a **file-driven, multi-pass extraction pipeline** built on **Outlines** (JSON-guided generation) and **Ollama** (local model serving). The pipeline reads configuration, prompts, and JSON Schemas from disk, runs a configured sequence of passes, writes versioned artifacts to an output directory (never overwriting), and can optionally persist results + timing metadata into an SQLite database.

## What it does

- **Pass-based extraction**: runs a configurable sequence of passes (e.g., `A_core`, `B_index`, `C_sequences`, ...) using Outlines JSON schema guidance.
- **Artifacts on disk**: writes raw text, pretty JSON, and error logs for each pass without overwriting prior runs.
- **Final stitching + validation**: stitches per-pass outputs into a final “FULL” object and can validate it against a “full schema” if configured.
- **SQLite optional**: can insert stitched results into SQLite via `hyb_db.insert_article_object(...)`.
- **Perf sidecars + continuation**: each JSON artifact can have a `*.perf.json` sidecar; the same metrics can also be mirrored into SQLite (`pipeline_artifacts`) and used for “resume” mode.

## Repository layout

The pipeline expects a “project directory” that contains:

- `config/pipeline.json` (main configuration)
- `passes/<pass_name>/schema.json` and `passes/<pass_name>/prompt.txt` (per-pass assets)
- `passes/common.txt` (shared prompt prefix, optional)
- `schema/json/article.json` (full schema for final validation, optional)
- input directory with source files (configured in `pipeline.json`)

The config shown in `config/pipeline.json` includes keys such as:
- `model_names`, `ollama_base_url`, `ollama_parameters`, `timeout_s`
- `input_dir`, `out_dir`, `article_glob`
- `pre_passes`, `construct_single_experiment_passes`, `passes`

## Installation

Python dependencies (minimum set used by the pipeline):

```bash
pip install -r requirements.txt
```

Or use the conda/mamba to initialize environment from `environment.yml`.

You also need:
- **Ollama** running locally (or reachable over HTTP), matching `ollama_base_url` in `config/pipeline.json`.

Optional:
- If `db_path` is set in config, SQLite will be used and schema will be auto-created.

### Environment variables

- `OPEN_BUTTON_TOKEN` (optional): if set, it is passed as a Bearer token in Ollama client headers.

## How to run

### 1) Configure `config/pipeline.json`

Edit paths to match your machine. In the attached example, `input_dir` is set to an absolute path and `article_glob` uses a recursive pattern.

Key fields you typically tune:
- `model_names`: list of Ollama model identifiers to run.
- `ollama_parameters`: e.g. `num_ctx`, `num_predict`, `temperature`, `seed`.
- `timeout_s`, `ollama_base_url`
- `out_dir`, `db_path`

### 2) Run the pipeline

#### CLI

From the repository root (or anywhere, as long as you pass the correct project directory):

```bash
python extraction/pipeline_filedriven.py extraction --fresh
```

- `project_dir` is the folder containing `config/`, `passes/`, etc.
- omit `--fresh` to enable continuation/resume behavior.

#### Python

```python
from extraction.pipeline_filedriven import run_project
run_project("extraction", fresh=False)
```

## Outputs

Artifacts are written under `out_dir` (from `pipeline.json`).

The pipeline writes, per pass and per model/article:
- raw text: `*.txt`
- JSON outputs: `*.json`
- log JSON: `*.log.json`
- error logs: `logs/*.log`
- perf sidecars: `*.perf.json` (one per emitted JSON artifact)

Perf sidecars include timestamps, wallclock duration, and (when Ollama reports it) token counts.

## Continuation / resume mode

When `db_path` is configured, the pipeline can skip already completed work:

- default `fresh=False`: for each `(model_name, article_name)`, if a successful `pass_name="FULL"` is recorded in the DB, the article can be skipped.
- `--fresh`: disables skipping and forces re-processing.

Implementation note:
- completion is tracked in `pipeline_artifacts` and queried via `hyb_db.get_completed_passes(...)`.

## Database schema (optional)

If `db_path` is set, `hyb_db` auto-creates tables and views and inserts:
- stitched article objects (`insert_article_object`)
- artifact-level perf bookkeeping (`pipeline_artifacts`)

## Overall design (short)

- **Config-first**: a project is a directory of config + prompts + schemas, making experiments easy to reproduce and version-control.
- **Multi-pass extraction**: each pass targets a specific sub-problem and produces a structured JSON artifact.
- **Immutable artifacts**: outputs are timestamped and never overwritten, enabling auditing and comparisons across runs.
- **Optional persistence**: results and metrics can be stored in SQLite for analysis and “resume” behavior.
