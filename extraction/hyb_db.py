# -*- coding: utf-8 -*-
# MIT License
#
# Copyright (c) 2025 CTLab-ITMO
#
# Authors: Aleksandr Serdiukov, Vitalii Dravgelis, Daniil Smutin, Artem Ivanov,
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


import outlines
from outlines.types import JsonSchema
import ollama
import re
from typing import Optional, Tuple, Dict, Any, List
import json
from pathlib import Path
from tqdm import tqdm
#from __future__ import annotations
import sqlite3
from contextlib import contextmanager, closing
from datetime import datetime, timezone
from loguru import logger
from ollama import chat, ChatResponse
from json_repair import repair_json
import os, sys
from jsonschema import Draft202012Validator

# -*- coding: utf-8 -*-
"""
SQLite dataset builder for hybridization-article extractions.

Public API:
    init_db(db_path)

    # HYBRIDIZATION ARTICLE / EXPERIMENT STRUCTURE
    insert_article_object(db_path, article_obj, model_name, article_name)
    insert_seqdesc_object(...)

    # PERFORMANCE / PIPELINE METRICS
    insert_perf_event(db_path, event_dict)
    insert_perf_events(db_path, [event_dict, ...])

    # NEW (ADDED FOR PIPELINE CONTINUATION + SIDE-CAR PERF TRACKING)
    insert_pipeline_artifact(db_path, artifact_dict)
    get_pipeline_artifacts_for_article(db_path, model_name, article_name)
    get_completed_passes(db_path, model_name, article_name)

Features:
- Auto-initializes schema (tables, indexes, views).
- Preserves every run (no overwrites).
- Normalizes sense/antisense & prime markers.
- Guards against non-oligo "probes" (skips probe insertion but keeps experiment).
- Includes Ollama-style helper tools with Google docstrings.
- perf_events table for timings/tokens of every step/question.
- NEW: pipeline_artifacts table for per-pass / per-file sidecar metrics and
       continuation bookkeeping. Each JSON artifact the pipeline writes on disk
       (per pass, per article, per model) can have a "sidecar" JSON with timing
       and token usage. We mirror that data into pipeline_artifacts so that:
         * downstream QC / benchmarking code can query timings and tokens
         * the pipeline can resume/continue work by checking which passes for a
           given (model_name, article_name) have already succeeded.
  The pipeline will:
    - emit a sidecar JSON next to every produced .json/.log.json/etc. file
      containing timing, token counts, and file paths
    - call insert_pipeline_artifact(...) with the same metadata
  Because this module always calls _ensure_schema() on connect, the new table
  will be created automatically in older existing DBs without migration steps.

Tables overview
---------------
articles / runs / raw_payloads / experiments / ... :
    Structured hybridization experiment data (final stitched objects).

seqdesc_* :
    Per-sequence descriptors from sequence descriptor passes.

perf_events :
    Fine-grained timing and token usage for any granular step/question.

pipeline_artifacts :
    Coarse-grained artifact-level bookkeeping used for:
        - sidecar perf metrics per produced JSON artifact
        - continuation / resume logic (which pipeline passes already finished)
"""

import json
import re
import sqlite3
from contextlib import contextmanager
from datetime import datetime, timezone
from typing import Any, Dict, Optional, Tuple, List


@contextmanager
def _db(db_path: str):
    """Context manager for SQLite connection with FK + WAL enabled."""
    conn = sqlite3.connect(db_path)
    try:
        conn.execute("PRAGMA foreign_keys = ON;")
        conn.execute("PRAGMA journal_mode = WAL;")
        yield conn
        conn.commit()
    finally:
        conn.close()


# ----------------------------- Schema DDL ----------------------------- #

_TABLES_AND_INDEXES_SQL = """
CREATE TABLE IF NOT EXISTS articles (
    id                  INTEGER PRIMARY KEY,
    doi                 TEXT NOT NULL UNIQUE,
    latest_article_name TEXT,
    latest_abstract     TEXT,
    latest_topic        TEXT,
    created_at          TEXT NOT NULL DEFAULT (datetime('now'))
);

CREATE TABLE IF NOT EXISTS runs (
    id           INTEGER PRIMARY KEY,
    article_id   INTEGER NOT NULL,
    model_name   TEXT    NOT NULL,
    article_name TEXT,
    branch       TEXT    NOT NULL CHECK (branch IN ('experiments','no_sequences')),
    created_at   TEXT    NOT NULL,
    FOREIGN KEY (article_id) REFERENCES articles(id) ON DELETE CASCADE
);
CREATE INDEX IF NOT EXISTS idx_runs_article ON runs(article_id);
CREATE INDEX IF NOT EXISTS idx_runs_created ON runs(created_at);
CREATE INDEX IF NOT EXISTS idx_runs_model   ON runs(model_name);

CREATE TABLE IF NOT EXISTS raw_payloads (
    run_id  INTEGER PRIMARY KEY,
    json    TEXT NOT NULL,
    FOREIGN KEY (run_id) REFERENCES runs(id) ON DELETE CASCADE
);

CREATE TABLE IF NOT EXISTS experiments (
    id                          INTEGER PRIMARY KEY,
    run_id                      INTEGER NOT NULL,
    id_exp                      TEXT    NOT NULL,
    type                        TEXT,
    description                 TEXT    NOT NULL,
    raw_description             TEXT,
    organism                    TEXT,
    technology                  TEXT,
    annealing_qualitative       INTEGER,  -- NULL/0/1
    rna_impurities_qualitative  INTEGER,  -- NULL/0/1
    UNIQUE (run_id, id_exp),
    FOREIGN KEY (run_id) REFERENCES runs(id) ON DELETE CASCADE
);
CREATE INDEX IF NOT EXISTS idx_experiments_run   ON experiments(run_id);
CREATE INDEX IF NOT EXISTS idx_experiments_idexp ON experiments(id_exp);

CREATE TABLE IF NOT EXISTS oligos (
    id                      INTEGER PRIMARY KEY,
    raw                     TEXT NOT NULL,
    sequence                TEXT,
    length_bases            INTEGER,
    prime_prefix            INTEGER CHECK (prime_prefix IN (3,5)),
    five_prime_label        TEXT,
    three_prime_label       TEXT,
    sense_antisense         TEXT CHECK (sense_antisense IN ('sense','antisense')),
    provenance_source_type  TEXT,
    provenance_page         INTEGER,
    provenance_section      TEXT,
    provenance_quote        TEXT,
    provenance_notes        TEXT
);
CREATE INDEX IF NOT EXISTS idx_oligos_seq ON oligos(sequence);

CREATE TABLE IF NOT EXISTS probes (
    id               INTEGER PRIMARY KEY,
    experiment_id    INTEGER NOT NULL,
    name             TEXT    NOT NULL,
    amplicon_id      TEXT,
    oligo_id         INTEGER NOT NULL,
    fluorophore      TEXT,
    quencher         TEXT,
    sense_antisense  TEXT CHECK (sense_antisense IN ('sense','antisense')),
    notes            TEXT,
    FOREIGN KEY (experiment_id) REFERENCES experiments(id) ON DELETE CASCADE,
    FOREIGN KEY (oligo_id)     REFERENCES oligos(id)
);
CREATE INDEX IF NOT EXISTS idx_probes_name ON probes(name);
CREATE INDEX IF NOT EXISTS idx_probes_exp  ON probes(experiment_id);

CREATE TABLE IF NOT EXISTS target_sequences (
    id             INTEGER PRIMARY KEY,
    experiment_id  INTEGER NOT NULL,
    oligo_id       INTEGER NOT NULL,
    FOREIGN KEY (experiment_id) REFERENCES experiments(id) ON DELETE CASCADE,
    FOREIGN KEY (oligo_id)     REFERENCES oligos(id)
);

CREATE TABLE IF NOT EXISTS primer_pairs (
    id               INTEGER PRIMARY KEY,
    experiment_id    INTEGER NOT NULL,
    forward_oligo_id INTEGER NOT NULL,
    reverse_oligo_id INTEGER NOT NULL,
    FOREIGN KEY (experiment_id)    REFERENCES experiments(id) ON DELETE CASCADE,
    FOREIGN KEY (forward_oligo_id) REFERENCES oligos(id),
    FOREIGN KEY (reverse_oligo_id) REFERENCES oligos(id)
);
CREATE INDEX IF NOT EXISTS idx_primers_exp ON primer_pairs(experiment_id);

CREATE TABLE IF NOT EXISTS related_sequences (
    id             INTEGER PRIMARY KEY,
    experiment_id  INTEGER NOT NULL,
    oligo_id       INTEGER NOT NULL,
    description    TEXT,
    FOREIGN KEY (experiment_id) REFERENCES experiments(id) ON DELETE CASCADE,
    FOREIGN KEY (oligo_id)     REFERENCES oligos(id)
);
CREATE INDEX IF NOT EXISTS idx_relseqs_exp ON related_sequences(experiment_id);

CREATE TABLE IF NOT EXISTS outcomes (
    id                 INTEGER PRIMARY KEY,
    experiment_id      INTEGER NOT NULL,
    outcome            INTEGER,   -- NULL/0/1
    comparative_notes  TEXT,
    FOREIGN KEY (experiment_id) REFERENCES experiments(id) ON DELETE CASCADE
);
CREATE INDEX IF NOT EXISTS idx_outcomes_exp ON outcomes(experiment_id);

CREATE TABLE IF NOT EXISTS measurements (
    id                      INTEGER PRIMARY KEY,
    experiment_id           INTEGER NOT NULL,
    key                     TEXT    NOT NULL,
    raw                     TEXT    NOT NULL,
    value                   REAL,
    unit                    TEXT,
    si_value                REAL,
    si_unit                 TEXT,
    assumptions             TEXT,
    provenance_source_type  TEXT,
    provenance_page         INTEGER,
    provenance_section      TEXT,
    provenance_quote        TEXT,
    provenance_notes        TEXT,
    FOREIGN KEY (experiment_id) REFERENCES experiments(id) ON DELETE CASCADE
);
CREATE INDEX IF NOT EXISTS idx_measurements_exp_key ON measurements(experiment_id, key);

CREATE TABLE IF NOT EXISTS pairings (
    id                        INTEGER PRIMARY KEY,
    experiment_id             INTEGER NOT NULL,
    paired_with_probe_name    TEXT,
    relationship              TEXT,
    FOREIGN KEY (experiment_id) REFERENCES experiments(id) ON DELETE CASCADE
);
CREATE INDEX IF NOT EXISTS idx_pairings_exp ON pairings(experiment_id);

CREATE TABLE IF NOT EXISTS extraction_report_entries (
    id            INTEGER PRIMARY KEY,
    run_id        INTEGER NOT NULL,
    experiment_id INTEGER,
    kind          TEXT NOT NULL CHECK (kind IN ('missing','uncertain')),
    json_pointer  TEXT NOT NULL,
    notes         TEXT,
    FOREIGN KEY (run_id)        REFERENCES runs(id) ON DELETE CASCADE,
    FOREIGN KEY (experiment_id) REFERENCES experiments(id) ON DELETE CASCADE
);
CREATE INDEX IF NOT EXISTS idx_report_run_kind ON extraction_report_entries(run_id, kind);

CREATE TABLE IF NOT EXISTS no_sequences_explanations (
    id          INTEGER PRIMARY KEY,
    run_id      INTEGER NOT NULL,
    explanation TEXT    NOT NULL,
    FOREIGN KEY (run_id) REFERENCES runs(id) ON DELETE CASCADE
);
CREATE INDEX IF NOT EXISTS idx_no_seq_run ON no_sequences_explanations(run_id);

/* generic performance/timing/token metrics for all steps and questions */
CREATE TABLE IF NOT EXISTS perf_events (
    id                   INTEGER PRIMARY KEY,
    namespace            TEXT NOT NULL CHECK (namespace IN ('pre_pass','pass','query','construct','stitch','db_insert','other')),
    article_name         TEXT,
    model_name           TEXT,
    article_doi          TEXT,
    pass_name            TEXT,
    sequence_key         TEXT,
    question_param       TEXT,
    started_at           TEXT,
    finished_at          TEXT,
    duration_ms          REAL,
    prompt_tokens        INTEGER,
    completion_tokens    INTEGER,
    total_tokens         INTEGER,
    tokens_per_sec       REAL,
    sidecar_path         TEXT,
    notes                TEXT
);
CREATE INDEX IF NOT EXISTS idx_perf_ns_article_model ON perf_events(namespace, article_name, model_name);

/* NEW TABLE:
   pipeline_artifacts captures artifact-level metadata (per produced JSON file).
   It is designed for:
     - performance sidecar ingestion (duration, token counts, sidecar paths)
     - resume/continuation bookkeeping (which passes are already done)
   The UNIQUE constraint prevents exact duplicate rows for the same (model,article,pass,file),
   but allows multiple historical attempts over time if artifact_path differs
   (timestamps are embedded in filenames in the pipeline).
*/
CREATE TABLE IF NOT EXISTS pipeline_artifacts (
    id                   INTEGER PRIMARY KEY,
    model_name           TEXT NOT NULL,
    article_name         TEXT NOT NULL,
    pass_name            TEXT NOT NULL,
    artifact_path        TEXT NOT NULL,
    sidecar_path         TEXT,
    started_at           TEXT,
    finished_at          TEXT,
    duration_ms          REAL,
    prompt_tokens        INTEGER,
    completion_tokens    INTEGER,
    total_tokens         INTEGER,
    tokens_per_sec       REAL,
    success              INTEGER,  -- NULL/0/1
    notes                TEXT,
    UNIQUE (model_name, article_name, pass_name, artifact_path)
);
CREATE INDEX IF NOT EXISTS idx_pa_lookup   ON pipeline_artifacts(model_name, article_name, pass_name);
CREATE INDEX IF NOT EXISTS idx_pa_finished ON pipeline_artifacts(finished_at);
"""

_VIEWS_SQL = """
CREATE VIEW IF NOT EXISTS view_experiments_flat AS
SELECT
    a.doi                               AS doi,
    r.model_name                        AS model_name,
    r.article_name                      AS article_name,
    r.created_at                        AS run_created_at,
    e.id                                AS experiment_id,
    e.id_exp                            AS id_exp,
    e.type                              AS exp_type,
    e.description                       AS exp_description,
    e.organism                          AS organism,
    e.technology                        AS technology,
    p.name                              AS probe_name,
    p.amplicon_id                       AS amplicon_id,
    p.fluorophore                       AS probe_fluorophore,
    p.quencher                          AS probe_quencher,
    po.sequence                         AS probe_sequence,
    po.five_prime_label                 AS probe_5p_label,
    po.three_prime_label                AS probe_3p_label,
    tgo.sequence                        AS target_sequence,
    o.outcome                           AS outcome_bool,
    o.comparative_notes                 AS outcome_notes
FROM experiments e
JOIN runs r           ON r.id = e.run_id
JOIN articles a       ON a.id = r.article_id
LEFT JOIN probes p    ON p.experiment_id = e.id
LEFT JOIN oligos po   ON po.id = p.oligo_id
LEFT JOIN target_sequences ts ON ts.experiment_id = e.id
LEFT JOIN oligos tgo   ON tgo.id = ts.oligo_id
LEFT JOIN outcomes o   ON o.experiment_id = e.id;

CREATE VIEW IF NOT EXISTS view_measurements_flat AS
SELECT
    a.doi,
    r.model_name,
    r.article_name,
    r.created_at        AS run_created_at,
    e.id                AS experiment_id,
    e.id_exp,
    m.key,
    m.raw,
    m.value,
    m.unit,
    m.si_value,
    m.si_unit,
    m.assumptions
FROM measurements m
JOIN experiments e ON e.id = m.experiment_id
JOIN runs r       ON r.id = e.run_id
JOIN articles a   ON a.id = r.article_id;
"""


def _ensure_schema(conn: sqlite3.Connection) -> None:
    """Create all tables, indexes, and views if they don't exist."""
    cur = conn.cursor()
    cur.executescript(_TABLES_AND_INDEXES_SQL)
    cur.executescript(_VIEWS_SQL)
    conn.commit()


# ----------------------------- Public API ----------------------------- #

def init_db(db_path: str) -> None:
    """Create (if not exists) the SQLite database schema, indices, and views.

    Args:
      db_path: Path to the SQLite file. Created if it doesn't exist.
    """
    with _db(db_path) as conn:
        _ensure_schema(conn)


def insert_article_object(db_path: str, article_obj: Dict[str, Any],
                          model_name: str, article_name: Optional[str]) -> int:
    """Insert a schema-conformant JSON object into the SQLite DB.

    Auto-creates the DB schema if missing. Preserves every run.

    Args:
      db_path: SQLite file path.
      article_obj: Dict that conforms to the Hybridization Article schema.
      model_name: Model identifier (e.g., 'Qwen2.5-Instruct-1M:14b').
      article_name: Name/key for the source file processed.

    Returns:
      run_id (int) for this insertion.
    """
    with _db(db_path) as conn:
        _ensure_schema(conn)
        cur = conn.cursor()

        doi = article_obj.get("doi", "unknown")
        if not doi:
            raise ValueError("Input must contain a top-level 'doi' string.")

        has_experiments = isinstance(article_obj.get("experiments"), list)
        branch = "experiments" if has_experiments else "no_sequences"

        article_id = _get_or_create_article(
            cur,
            doi=doi,
            article_name=article_name or article_obj.get("article_name"),
            abstract=article_obj.get("abstract"),
            topic=article_obj.get("topic"),
        )

        run_id = _create_run(cur, article_id, model_name, article_name, branch, raw_json=article_obj)

        # Top-level extraction report (if any)
        _insert_extraction_report(cur, run_id, article_obj.get("extraction_report"), experiment_id=None)

        if branch == "no_sequences":
            explanation = article_obj.get("explanation_why_does_not_this_article_have_any_hybridization_probes_sequences") or ""
            cur.execute(
                "INSERT INTO no_sequences_explanations (run_id, explanation) VALUES (?, ?)",
                (run_id, explanation),
            )
            return run_id

        # ---- experiments branch ----
        for exp in (article_obj.get("experiments") or []):
            id_exp = exp.get("id_exp")
            desc = exp.get("description") or ""
            cur.execute(
                """
                INSERT INTO experiments
                    (run_id, id_exp, type, description, raw_description,
                     organism, technology, annealing_qualitative, rna_impurities_qualitative)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
                """,
                (
                    run_id,
                    id_exp,
                    exp.get("type"),
                    desc,
                    exp.get("raw_description"),
                    (exp.get("metadata") or {}).get("organism"),
                    (exp.get("metadata") or {}).get("technology"),
                    _to_int_bool(((exp.get("metadata") or {}).get("annealing") or {}).get("qualitative")),
                    _to_int_bool(((exp.get("metadata") or {}).get("rna_impurities") or {}).get("qualitative")),
                ),
            )
            experiment_id = cur.lastrowid

            # Per-experiment extraction report
            _insert_extraction_report(cur, run_id, exp.get("extraction_report"), experiment_id=experiment_id)

            # Sequences
            seqs = exp.get("sequences") or {}

            # Validate the probe looks like a real oligo
            probe = seqs.get("probe") or {}
            if not _has_real_probe(probe):
                # Record and skip probe insertion, but keep experiment row and any metadata/measurements/outcomes
                _insert_extraction_report(
                    cur, run_id,
                    {"missing": ["/experiments/*/sequences/probe/oligo/sequence"],
                     "notes": "Rejected probable non-oligo probe (no bases/labels/length)."},
                    experiment_id=experiment_id,
                )
            else:
                # Target (optional)
                tgt = seqs.get("target_sequence")
                if isinstance(tgt, dict) and (tgt.get("raw") is not None):
                    tgt_oligo_id = _insert_oligo(cur, tgt)
                    cur.execute(
                        "INSERT INTO target_sequences (experiment_id, oligo_id) VALUES (?, ?)",
                        (experiment_id, tgt_oligo_id),
                    )

                # Probe (required by schema; normalized before insert)
                probe_oligo = probe.get("oligo") or {}
                probe_oligo_id = _insert_oligo(cur, probe_oligo)
                sa = _coerce_sa(probe.get("sense_antisense"), probe.get("name"))
                cur.execute(
                    """
                    INSERT INTO probes
                        (experiment_id, name, amplicon_id, oligo_id, fluorophore, quencher, sense_antisense, notes)
                    VALUES (?, ?, ?, ?, ?, ?, ?, ?)
                    """,
                    (
                        experiment_id,
                        probe.get("name"),
                        probe.get("amplicon_id"),
                        probe_oligo_id,
                        probe.get("fluorophore"),
                        probe.get("quencher"),
                        sa,
                        probe.get("notes"),
                    ),
                )

                # Primers (optional)
                primers = seqs.get("primer_sequences")
                if isinstance(primers, dict):
                    fwd = primers.get("forward") or {}
                    rev = primers.get("reverse") or {}
                    fwd_id = _insert_oligo(cur, fwd)
                    rev_id = _insert_oligo(cur, rev)
                    cur.execute(
                        "INSERT INTO primer_pairs (experiment_id, forward_oligo_id, reverse_oligo_id) VALUES (?, ?, ?)",
                        (experiment_id, fwd_id, rev_id),
                    )

                # Related sequences (0..N)
                for rs in (seqs.get("related_sequences") or []):
                    r_oligo = rs.get("related_sequence")
                    if isinstance(r_oligo, dict) and (r_oligo.get("raw") is not None):
                        r_oligo_id = _insert_oligo(cur, r_oligo)
                        cur.execute(
                            "INSERT INTO related_sequences (experiment_id, oligo_id, description) VALUES (?, ?, ?)",
                            (experiment_id, r_oligo_id, rs.get("description")),
                        )

            # Measurements (experiment_properties + metadata)
            exprops = exp.get("experiment_properties") or {}
            concs = (exprops.get("concentrations") or {})
            _insert_measurement(cur, experiment_id, "experiment_properties.concentrations.dna_rna_concentration",
                                concs.get("dna_rna_concentration"))
            _insert_measurement(cur, experiment_id, "experiment_properties.concentrations.concentration_SI",
                                concs.get("concentration_SI"))

            params = (exprops.get("parameters_SI") or {})
            for key in ("temperature", "Tris", "Na", "K", "Mg", "DMSO"):
                _insert_measurement(cur, experiment_id, f"experiment_properties.parameters_SI.{key}", params.get(key))

            meta = exp.get("metadata") or {}
            _insert_measurement(cur, experiment_id, "metadata.pH", meta.get("pH"))
            ann = meta.get("annealing") or {}
            _insert_measurement(cur, experiment_id, "metadata.annealing.quantitative", ann.get("quantitative"))
            rimp = meta.get("rna_impurities") or {}
            _insert_measurement(cur, experiment_id, "metadata.rna_impurities.quantitative", rimp.get("quantitative"))

            # Outcome
            out = exp.get("outcome") or {}
            cur.execute(
                "INSERT INTO outcomes (experiment_id, outcome, comparative_notes) VALUES (?, ?, ?)",
                (experiment_id, _to_int_bool(out.get("outcome")), out.get("comparative_notes")),
            )
            _insert_measurement(cur, experiment_id, "outcome.fluorescence", out.get("fluorescence"))

            # Pairing (optional)
            pair = exp.get("pairing") or {}
            if pair.get("paired_with_probe_name") or pair.get("relationship"):
                cur.execute(
                    "INSERT INTO pairings (experiment_id, paired_with_probe_name, relationship) VALUES (?, ?, ?)",
                    (experiment_id, pair.get("paired_with_probe_name"), pair.get("relationship")),
                )

        return run_id


# NEW: perf events API -------------------------------------------------- #

def _event_defaults(ev: Dict[str, Any]) -> Dict[str, Any]:
    d = dict(ev or {})
    for k in ("namespace","article_name","model_name","article_doi","pass_name",
              "sequence_key","question_param","started_at","finished_at",
              "duration_ms","prompt_tokens","completion_tokens","total_tokens",
              "tokens_per_sec","sidecar_path","notes"):
        d.setdefault(k, None)
    return d

def insert_perf_event(db_path: str, event: Dict[str, Any]) -> int:
    """Insert a single performance/timing event row.

    Typical usage:
        insert_perf_event(db_path, {
            "namespace": "pass",
            "article_name": "paper123",
            "model_name": "my/model:7b",
            "article_doi": "10.xxxx/yyy",
            "pass_name": "A_core",
            "sequence_key": None,
            "question_param": None,
            "started_at": "...",
            "finished_at": "...",
            "duration_ms": 1234.5,
            "prompt_tokens": 4567,
            "completion_tokens": 890,
            "total_tokens": 5457,
            "tokens_per_sec": 12.34,
            "sidecar_path": "/path/to/file.sidecar.json",
            "notes": "ok"
        })
    """
    with _db(db_path) as conn:
        _ensure_schema(conn)
        cur = conn.cursor()
        e = _event_defaults(event)
        cur.execute(
            """
            INSERT INTO perf_events (
                namespace, article_name, model_name, article_doi, pass_name,
                sequence_key, question_param, started_at, finished_at, duration_ms,
                prompt_tokens, completion_tokens, total_tokens, tokens_per_sec,
                sidecar_path, notes
            )
            VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)
            """,
            (
                e["namespace"], e["article_name"], e["model_name"], e["article_doi"], e["pass_name"],
                e["sequence_key"], e["question_param"], e["started_at"], e["finished_at"], e["duration_ms"],
                e["prompt_tokens"], e["completion_tokens"], e["total_tokens"], e["tokens_per_sec"],
                e["sidecar_path"], e["notes"]
            ),
        )
        return cur.lastrowid

def insert_perf_events(db_path: str, events: List[Dict[str, Any]]) -> List[int]:
    """Bulk insert multiple performance events."""
    ids: List[int] = []
    if not events:
        return ids
    with _db(db_path) as conn:
        _ensure_schema(conn)
        cur = conn.cursor()
        for ev in events:
            e = _event_defaults(ev)
            cur.execute(
                """
                INSERT INTO perf_events (
                    namespace, article_name, model_name, article_doi, pass_name,
                    sequence_key, question_param, started_at, finished_at, duration_ms,
                    prompt_tokens, completion_tokens, total_tokens, tokens_per_sec,
                    sidecar_path, notes
                )
                VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)
                """,
                (
                    e["namespace"], e["article_name"], e["model_name"], e["article_doi"], e["pass_name"],
                    e["sequence_key"], e["question_param"], e["started_at"], e["finished_at"], e["duration_ms"],
                    e["prompt_tokens"], e["completion_tokens"], e["total_tokens"], e["tokens_per_sec"],
                    e["sidecar_path"], e["notes"]
                ),
            )
            ids.append(cur.lastrowid)
    return ids


# NEW: pipeline_artifacts API ------------------------------------------ #

def _artifact_defaults(rec: Dict[str, Any]) -> Dict[str, Any]:
    """Normalize/complete an artifact record dict before DB insert.

    Expected keys in `rec`:
        model_name         : str
        article_name       : str
        pass_name          : str  (e.g. 'A_core', 'SeqDesc-OPTIM', 'FULL')
        artifact_path      : str  (absolute or project-rel path to main JSON)
        sidecar_path       : str|None (path to sidecar .perf.json)
        started_at         : str|None (ISO8601 UTC)
        finished_at        : str|None (ISO8601 UTC)
        duration_ms        : float|None
        prompt_tokens      : int|None
        completion_tokens  : int|None
        total_tokens       : int|None
        tokens_per_sec     : float|None
        success            : bool|int|None
        notes              : str|None
    """
    d = dict(rec or {})
    for k in (
        "model_name", "article_name", "pass_name", "artifact_path",
        "sidecar_path", "started_at", "finished_at", "duration_ms",
        "prompt_tokens", "completion_tokens", "total_tokens", "tokens_per_sec",
        "success", "notes"
    ):
        d.setdefault(k, None)

    # Convert success -> int bool (1/0/NULL)
    d["success"] = _to_int_bool(d.get("success"))
    return d


def insert_pipeline_artifact(db_path: str, artifact: Dict[str, Any]) -> int:
    """Insert a single pipeline_artifacts row.

    This captures per-pass / per-article / per-model artifact bookkeeping
    (timings + token usage for the JSON file) as well as marking a pass
    as 'successfully finished' for continuation.

    NOTE:
    We do a plain INSERT. The table has a UNIQUE constraint on
    (model_name, article_name, pass_name, artifact_path), so the pipeline
    should generate unique artifact_path names (it already includes a
    timestamped suffix in filenames). We intentionally do NOT overwrite.

    Returns:
      row_id (int).
    """
    with _db(db_path) as conn:
        _ensure_schema(conn)
        cur = conn.cursor()
        r = _artifact_defaults(artifact)
        cur.execute(
            """
            INSERT INTO pipeline_artifacts (
                model_name, article_name, pass_name, artifact_path, sidecar_path,
                started_at, finished_at, duration_ms,
                prompt_tokens, completion_tokens, total_tokens, tokens_per_sec,
                success, notes
            )
            VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?)
            """,
            (
                r["model_name"],
                r["article_name"],
                r["pass_name"],
                r["artifact_path"],
                r["sidecar_path"],
                r["started_at"],
                r["finished_at"],
                r["duration_ms"],
                r["prompt_tokens"],
                r["completion_tokens"],
                r["total_tokens"],
                r["tokens_per_sec"],
                r["success"],
                r["notes"],
            ),
        )
        return cur.lastrowid


def get_pipeline_artifacts_for_article(
    db_path: str,
    model_name: str,
    article_name: str,
) -> List[Dict[str, Any]]:
    """Fetch all recorded artifacts for (model_name, article_name).

    This is useful for:
      - debugging
      - external QC scripts
      - continuation logic in the pipeline (to see what's already done)

    Returns:
      A list (possibly empty) of dicts, newest-first by finished_at (NULL last).
    """
    with _db(db_path) as conn:
        _ensure_schema(conn)
        cur = conn.cursor()
        cur.execute(
            """
            SELECT
                id,
                model_name,
                article_name,
                pass_name,
                artifact_path,
                sidecar_path,
                started_at,
                finished_at,
                duration_ms,
                prompt_tokens,
                completion_tokens,
                total_tokens,
                tokens_per_sec,
                success,
                notes
            FROM pipeline_artifacts
            WHERE model_name = ?
              AND article_name = ?
            ORDER BY
                CASE WHEN finished_at IS NULL THEN 1 ELSE 0 END,
                finished_at DESC
            """,
            (model_name, article_name),
        )
        rows = cur.fetchall()

    out_rows: List[Dict[str, Any]] = []
    for row in rows:
        (
            rid, mname, aname, pass_name, artifact_path, sidecar_path,
            started_at, finished_at, duration_ms,
            prompt_tokens, completion_tokens, total_tokens, tokens_per_sec,
            success, notes
        ) = row
        out_rows.append(
            {
                "id": rid,
                "model_name": mname,
                "article_name": aname,
                "pass_name": pass_name,
                "artifact_path": artifact_path,
                "sidecar_path": sidecar_path,
                "started_at": started_at,
                "finished_at": finished_at,
                "duration_ms": duration_ms,
                "prompt_tokens": prompt_tokens,
                "completion_tokens": completion_tokens,
                "total_tokens": total_tokens,
                "tokens_per_sec": tokens_per_sec,
                "success": bool(success) if success is not None else None,
                "notes": notes,
            }
        )
    return out_rows


def get_completed_passes(
    db_path: str,
    model_name: str,
    article_name: str,
) -> Dict[str, Dict[str, Any]]:
    """Return a summary of which logical passes have already completed
    successfully for (model_name, article_name).

    The pipeline can use this to implement continuation / resume:
      - If a pass_name appears here with success True, the pipeline MAY skip
        regenerating that pass for that (model, article), unless --fresh was
        requested or the pipeline.json layout changed and the user wants that
        pass rerun anyway.

    Returns:
      dict:
        {
           "A_core": {
               "finished_at": "...",
               "artifact_path": ".../paper__A_core__model__timestamp.json",
               "sidecar_path": ".../paper__A_core__model__timestamp.perf.json",
               "duration_ms": 1234.5,
               "total_tokens": 9876,
               ...
           },
           ...
        }

      If multiple rows exist for the same pass_name we pick the most recent row
      with success == 1 (by finished_at DESC). If none are successful for a
      pass, that pass won't appear in the dict.
    """
    artifacts = get_pipeline_artifacts_for_article(
        db_path=db_path, model_name=model_name, article_name=article_name
    )

    best: Dict[str, Dict[str, Any]] = {}
    for row in artifacts:
        pname = row["pass_name"]
        if not row.get("success"):
            continue
        prev = best.get(pname)
        if prev is None:
            best[pname] = row
            continue
        # pick newer finished_at
        prev_finished = prev.get("finished_at")
        cur_finished = row.get("finished_at")
        # if prev_finished is None but cur_finished not None -> prefer current
        # if both not None -> compare lexicographically (ISO8601 so OK)
        take_current = False
        if prev_finished is None and cur_finished is not None:
            take_current = True
        elif prev_finished is not None and cur_finished is not None:
            if str(cur_finished) > str(prev_finished):
                take_current = True
        elif prev_finished is None and cur_finished is None:
            # keep first, arbitrary
            take_current = False
        # else prev has finished_at but current doesn't — keep prev.
        if take_current:
            best[pname] = row

    # Only expose a shallow summary for convenience
    summarized: Dict[str, Dict[str, Any]] = {}
    for pname, row in best.items():
        summarized[pname] = {
            "finished_at": row.get("finished_at"),
            "artifact_path": row.get("artifact_path"),
            "sidecar_path": row.get("sidecar_path"),
            "duration_ms": row.get("duration_ms"),
            "prompt_tokens": row.get("prompt_tokens"),
            "completion_tokens": row.get("completion_tokens"),
            "total_tokens": row.get("total_tokens"),
            "tokens_per_sec": row.get("tokens_per_sec"),
            "notes": row.get("notes"),
        }
    return summarized


# ----------------------------- Ollama-style helper tools ----------------------------- #

def to_si(value: Optional[float], unit: Optional[str]) -> Tuple[Optional[float], Optional[str]]:
    """Convert a numeric value and unit to SI.

    Supports common units from hybridization papers:
    - Temperature: °C -> K (K = °C + 273.15), K stays K.
    - Concentration: M, mM, µM/um, nM -> mol/m^3 (1 mM = 1 mol/m^3).
    - Percent: % -> dimensionless fraction (value/100).

    Args:
      value: The numeric value parsed from the article, or None if unknown.
      unit: The unit string as written in the article (e.g., '°C', 'C', 'mM', '%'), or None.

    Returns:
      A pair (si_value, si_unit):
        - si_value: The value converted to SI, or None if not convertible.
        - si_unit: The SI unit string ('K', 'mol/m^3', 'dimensionless'), or None if not convertible.

    Examples:
      >>> to_si(25, '°C')
      (298.15, 'K')
      >>> to_si(2, 'mM')
      (2.0, 'mol/m^3')
      >>> to_si(10, '%')
      (0.1, 'dimensionless')
    """
    if value is None or unit is None:
        return None, None

    u = unit.strip().lower().replace('µ', 'u')
    # Temperature
    if u in {'c', '°c', 'deg c', 'celsius'}:
        return value + 273.15, 'K'
    if u in {'k', 'kelvin'}:
        return value, 'K'

    # Concentration (to mol/m^3)
    if u in {'m', 'mol/l'}:
        return value * 1000.0, 'mol/m^3'
    if u in {'mm', 'mmol/l', 'mmol', 'mm'}:  # 'mm' for mM as often OCR'd
        return value * 1.0, 'mol/m^3'
    if u in {'um', 'umol/l', 'µm', 'µmol/l', 'micromolar'}:
        return value * 1e-3, 'mol/m^3'
    if u in {'nm', 'nmol/l', 'nanomolar'}:
        return value * 1e-6, 'mol/m^3'

    # Percent
    if u in {'%', 'percent', 'perc'}:
        return value / 100.0, 'dimensionless'

    return None, None


OLIGO_RE = re.compile(r"""
^\s*
(?:(?P<prime>(?:5|3)(?:['′’]|0|O)?)\s*-\s*)?
(?:(?P<prefix>(?:[A-Za-z0-9+]+-)+))?
(?P<seq>[ACGUTRYSWKMBDHVN]+)
(?:(?P<suffix>(?:-[A-Za-z0-9+]+)+))?
(?:\s*\(\s*(?P<len>\d+)\s*(?:b|bp)\s*\)\s*)?
\s*$
""", re.X)

def parse_oligo(raw: Optional[str]) -> Dict[str, Any]:
    """Parse a decorated oligo string into schema-ready parts.

    Accepts OCR-prone patterns like "50-FAM-...-BHQ2 (27 b)" and normalizes:
    - prime_prefix: 5 or 3 when 5′/3′ (includes 50/5O variants)
    - sequence: IUPAC bases (uppercase)
    - length_bases: integer if present
    - labels: all labels in order; five_prime_label and three_prime_label are the first/last, respectively

    Args:
      raw: The exact oligo string from the article (may include labels and length), or None.

    Returns:
      A dict matching the 'decoratedOligo' shape (minus provenance):
        {
          "raw": str or None,
          "sequence": str or None,
          "length_bases": int or None,
          "prime_prefix": 5|3|None,
          "five_prime_label": str or None,
          "three_prime_label": str or None,
          "labels": List[str],
          "sense_antisense": None
        }
    """
    result: Dict[str, Any] = {
        "raw": raw,
        "sequence": None,
        "length_bases": None,
        "prime_prefix": None,
        "five_prime_label": None,
        "three_prime_label": None,
        "labels": [],
        "sense_antisense": None
    }
    if not raw:
        return result

    m = OLIGO_RE.match(raw)
    if not m:
        return result

    prime = m.group('prime')
    if prime:
        result["prime_prefix"] = 5 if prime.startswith('5') else 3

    seq = m.group('seq')
    if seq:
        result["sequence"] = seq.upper()

    if m.group('len'):
        result["length_bases"] = int(m.group('len'))

    labels: List[str] = []
    if m.group('prefix'):
        labels += [x for x in m.group('prefix').split('-') if x]
    if m.group('suffix'):
        labels += [x for x in m.group('suffix').split('-') if x]
    result["labels"] = labels
    if labels:
        result["five_prime_label"] = labels[0]
        result["three_prime_label"] = labels[-1]

    return result


def make_measurement(raw: Optional[str],
                     value: Optional[float] = None,
                     unit: Optional[str] = None) -> Dict[str, Any]:
    """Build a 'measurement' object with SI conversion.

    Convenience helper to populate the schema's measurement type while keeping the raw text.

    Args:
      raw: The raw textual measurement from the article (e.g., '58 °C', '2 mM', '10%').
      value: Parsed numeric value, if available.
      unit: Parsed unit string as written in the article (e.g., '°C', 'mM', '%').

    Returns:
      A dict with keys: raw, value, unit, si_value, si_unit, assumptions.
      Unknown or unsupported units yield si_value/si_unit = None.
    """
    si_value, si_unit = to_si(value, unit) if (value is not None and unit is not None) else (None, None)
    return {
        "raw": raw or "",
        "value": value,
        "unit": unit,
        "si_value": si_value,
        "si_unit": si_unit,
        "assumptions": None
    }


# ----------------------------- Normalization / validation helpers ----------------------------- #

_SA_MAP = {
    's': 'sense',
    'sense': 'sense',
    'as': 'antisense',
    'antisense': 'antisense',
    '+': 'sense',
    '-': 'antisense',
    'forward': 'sense',
    'reverse': 'antisense',
}
_SA_NAME_RE = re.compile(r"\)\s*(as|s)\s*$", re.IGNORECASE)

def _detect_sa_from_name(probe_name: Optional[str]) -> Optional[str]:
    """Infer sense/antisense from a trailing '(...)s' or '(...)as' in the probe name.

    Args:
      probe_name: Probe name (e.g., 'N3-FAM(27)s').

    Returns:
      'sense', 'antisense', or None if not inferable.
    """
    if not probe_name:
        return None
    m = _SA_NAME_RE.search(probe_name.strip())
    if not m:
        return None
    g = m.group(1).lower()
    return 'antisense' if g == 'as' else 'sense'


def _coerce_sa(value: Optional[str], probe_name: Optional[str] = None) -> Optional[str]:
    """Coerce various encodings to 'sense'/'antisense'/None.

    Args:
      value: A string like 's', 'as', 'sense', 'antisense', '+', '-', or None.
      probe_name: Fallback context to infer sense/antisense from the name suffix.

    Returns:
      'sense', 'antisense', or None.
    """
    if value is None or (isinstance(value, str) and not value.strip()):
        return _detect_sa_from_name(probe_name)
    v = str(value).strip().lower()
    if v in _SA_MAP:
        return _SA_MAP[v]
    return _detect_sa_from_name(probe_name)


def _coerce_prime_prefix(value: Any) -> Optional[int]:
    """Clamp prime prefix to {3, 5} or None.

    Handles OCR-like strings such as '5', "5'", '50', '5O', '5′'.

    Args:
      value: Raw input for prime prefix.

    Returns:
      3, 5, or None.
    """
    if value is None:
        return None
    s = str(value).strip()
    if s.startswith('5'):
        return 5
    if s.startswith('3'):
        return 3
    return None


def _has_real_probe(probe: Dict[str, Any]) -> bool:
    """Heuristic gate: reject obviously non-oligo 'probes'.

    Accepts a probe only if at least one of these holds:
      - >= 6 IUPAC bases appear in oligo.sequence or oligo.raw
      - a known label is present (FAM/ROX/Cy5/BHQ1/BHQ2/RTQ1) in labels/five/three
      - length_bases is present

    Args:
      probe: The 'probe' dict from the schema.

    Returns:
      True if looks like a real oligo; False otherwise.
    """
    if not isinstance(probe, dict):
        return False
    oligo = probe.get("oligo") or {}
    raw = (oligo.get("raw") or "")
    seq = (oligo.get("sequence") or "")
    has_bases = bool(re.search(r"[ACGUTRYSWKMBDHVN]{6,}", (seq or raw).upper()))
    has_label = any(bool(oligo.get(k)) for k in ("five_prime_label", "three_prime_label")) \
                or bool(oligo.get("labels"))
    has_length = bool(oligo.get("length_bases"))
    return has_bases or has_label or has_length


# ----------------------------- DB helpers ----------------------------- #

def _utcnow_iso() -> str:
    """UTC timestamp in ISO8601 format."""
    return datetime.now(timezone.utc).isoformat()


def _get_or_create_article(cur: sqlite3.Cursor, doi: str,
                           article_name: Optional[str],
                           abstract: Optional[str],
                           topic: Optional[str]) -> int:
    """Fetch article.id by DOI, creating the row if needed (and refreshing metadata)."""
    cur.execute("SELECT id FROM articles WHERE doi = ?", (doi,))
    row = cur.fetchone()
    if row:
        article_id = row[0]
        cur.execute(
            """
            UPDATE articles
               SET latest_article_name = COALESCE(?, latest_article_name),
                   latest_abstract     = COALESCE(?, latest_abstract),
                   latest_topic        = COALESCE(?, latest_topic)
             WHERE id = ?
            """,
            (article_name, abstract, topic, article_id),
        )
        return article_id
    cur.execute(
        """
        INSERT INTO articles (doi, latest_article_name, latest_abstract, latest_topic, created_at)
        VALUES (?, ?, ?, ?, ?)
        """,
        (doi, article_name, abstract, topic, _utcnow_iso()),
    )
    return cur.lastrowid


def _create_run(cur: sqlite3.Cursor, article_id: int, model_name: str,
                article_name: Optional[str], branch: str,
                raw_json: Dict[str, Any]) -> int:
    """Create a run row and persist the raw JSON payload."""
    cur.execute(
        """
        INSERT INTO runs (article_id, model_name, article_name, branch, created_at)
        VALUES (?, ?, ?, ?, ?)
        """,
        (article_id, model_name, article_name, branch, _utcnow_iso()),
    )
    run_id = cur.lastrowid
    cur.execute("INSERT INTO raw_payloads (run_id, json) VALUES (?, ?)",
                (run_id, json.dumps(raw_json, ensure_ascii=False)))
    return run_id


def _insert_provenance_cols(entity: Dict[str, Any]) -> Tuple[Optional[str], Optional[int], Optional[str], Optional[str], Optional[str]]:
    """Extract provenance fields with safe defaults."""
    prov = entity.get("provenance") or {}
    return (
        prov.get("source_type"),
        prov.get("page"),
        prov.get("section"),
        prov.get("quote"),
        prov.get("notes"),
    )


def _insert_oligo(cur: sqlite3.Cursor, oligo: Dict[str, Any]) -> int:
    """Insert an oligo row after normalizing prime_prefix and sense/antisense."""
    cleaned = dict(oligo or {})
    cleaned["prime_prefix"] = _coerce_prime_prefix(cleaned.get("prime_prefix"))
    cleaned["sense_antisense"] = _coerce_sa(cleaned.get("sense_antisense"))

    ps, pg, sc, qu, no = _insert_provenance_cols(cleaned)
    cur.execute(
        """
        INSERT INTO oligos
            (raw, sequence, length_bases, prime_prefix,
             five_prime_label, three_prime_label, sense_antisense,
             provenance_source_type, provenance_page, provenance_section,
             provenance_quote, provenance_notes)
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """,
        (
            cleaned.get("raw", ""),
            cleaned.get("sequence"),
            cleaned.get("length_bases"),
            cleaned.get("prime_prefix"),
            cleaned.get("five_prime_label"),
            cleaned.get("three_prime_label"),
            cleaned.get("sense_antisense"),
            ps, pg, sc, qu, no,
        ),
    )
    return cur.lastrowid


def _insert_measurement(cur: sqlite3.Cursor, experiment_id: int, key: str, m: Optional[Dict[str, Any]]) -> None:
    """Insert a measurement if present."""
    if not m:
        return
    ps, pg, sc, qu, no = _insert_provenance_cols(m)
    cur.execute(
        """
        INSERT INTO measurements
            (experiment_id, key, raw, value, unit, si_value, si_unit, assumptions,
             provenance_source_type, provenance_page, provenance_section, provenance_quote, provenance_notes)
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """,
        (
            experiment_id,
            key,
            (m.get("raw") or ""),
            m.get("value"),
            m.get("unit"),
            m.get("si_value"),
            m.get("si_unit"),
            m.get("assumptions"),
            ps, pg, sc, qu, no,
        ),
    )


def _insert_extraction_report(cur: sqlite3.Cursor, run_id: int,
                              report: Optional[Dict[str, Any]],
                              experiment_id: Optional[int] = None) -> None:
    """Insert extraction report entries (missing/uncertain pointers)."""
    if not report:
        return
    for kind in ("missing", "uncertain"):
        for ptr in report.get(kind, []) or []:
            cur.execute(
                """
                INSERT INTO extraction_report_entries (run_id, experiment_id, kind, json_pointer, notes)
                VALUES (?, ?, ?, ?, ?)
                """,
                (run_id, experiment_id, kind, str(ptr), report.get("notes")),
            )


def _to_int_bool(val: Optional[bool]) -> Optional[int]:
    """Convert Python bool/None -> 1/0/NULL for SQLite."""
    if val is None:
        return None
    return 1 if bool(val) else 0

# ──────────────────────────────────────────────────────────────────────
# Sequence-descriptors DB (no collision; separate "seqdesc_*" namespace)
# ──────────────────────────────────────────────────────────────────────

def _extract_doi_from_text(text: str) -> Optional[str]:
    """Heuristic DOI extractor from article text (fallback)."""
    if not text:
        return None
    m = re.search(r"\b10\.\d{4,9}/[^\s\"'<>]+", text, flags=re.I)
    return m.group(0).rstrip(".,);]") if m else None

def _ensure_seqdesc_schema(conn: sqlite3.Connection) -> None:
    """Create the seqdesc_* schema if it does not exist."""
    conn.execute("PRAGMA foreign_keys = ON;")
    # Runs table
    conn.execute("""
        CREATE TABLE IF NOT EXISTS seqdesc_runs (
            id            INTEGER PRIMARY KEY AUTOINCREMENT,
            created_at    TEXT NOT NULL,
            model_name    TEXT NOT NULL,
            article_name  TEXT NOT NULL,
            doi           TEXT,
            source_path   TEXT,
            raw_json      TEXT NOT NULL
        );
    """)
    conn.execute("CREATE INDEX IF NOT EXISTS idx_seqdesc_runs_article ON seqdesc_runs(article_name);")
    conn.execute("CREATE INDEX IF NOT EXISTS idx_seqdesc_runs_doi     ON seqdesc_runs(doi);")
    conn.execute("CREATE INDEX IF NOT EXISTS idx_seqdesc_runs_model   ON seqdesc_runs(model_name);")

    # Sequences table (one row per sequence key in the run)
    conn.execute("""
        CREATE TABLE IF NOT EXISTS seqdesc_sequences (
            id                              INTEGER PRIMARY KEY AUTOINCREMENT,
            run_id                          INTEGER NOT NULL REFERENCES seqdesc_runs(id) ON DELETE CASCADE,
            sequence_key                    TEXT NOT NULL,   -- the dict key (probe string as found)
            is_seq                          INTEGER,         -- NULL/0/1
            sequence_full                   TEXT,
            sequence_normalized             TEXT,
            sequence_expanded               TEXT,
            sequence_backbone               TEXT,
            sequence_backbone_expanded      TEXT,
            fluorophore                     TEXT,
            quencher                        TEXT,
            target_raw                      TEXT,
            target_normalized               TEXT,
            primers_forward                 TEXT,
            primers_reverse                 TEXT,
            pH                              REAL,
            annealing_raw                   TEXT,
            T_value                         REAL,
            T_unit                          TEXT,
            Tris_value                      REAL,
            Tris_unit                       TEXT,
            Na_value                        REAL,
            Na_unit                         TEXT,
            K_value                         REAL,
            K_unit                          TEXT,
            Mg_value                        REAL,
            Mg_unit                         TEXT,
            DMSO_value                      REAL,
            DMSO_unit                       TEXT,
            outcome                         INTEGER,         -- NULL/0/1
            raw_json                        TEXT NOT NULL
        );
    """)
    conn.execute("CREATE INDEX IF NOT EXISTS idx_seqdesc_sequences_run ON seqdesc_sequences(run_id);")

    # Modifications table (0..N per sequence)
    conn.execute("""
        CREATE TABLE IF NOT EXISTS seqdesc_modifications (
            id                  INTEGER PRIMARY KEY AUTOINCREMENT,
            sequence_id         INTEGER NOT NULL REFERENCES seqdesc_sequences(id) ON DELETE CASCADE,
            modification_position INTEGER,
            modification_type   TEXT,
            modification_description TEXT
        );
    """)
    conn.execute("CREATE INDEX IF NOT EXISTS idx_seqdesc_mods_seq ON seqdesc_modifications(sequence_id);")

    # Helpful views (namespaced)
    conn.execute("""
        CREATE VIEW IF NOT EXISTS seqdesc_v_sequences AS
        SELECT
            r.id           AS run_id,
            r.created_at   AS run_created_at,
            r.model_name   AS model_name,
            r.article_name AS article_name,
            r.doi          AS doi,
            s.*
        FROM seqdesc_sequences s
        JOIN seqdesc_runs r ON r.id = s.run_id;
    """)
    conn.execute("""
        CREATE VIEW IF NOT EXISTS seqdesc_v_modifications AS
        SELECT
            s.run_id,
            s.id          AS sequence_id,
            s.sequence_key,
            m.modification_position,
            m.modification_type,
            m.modification_description
        FROM seqdesc_modifications m
        JOIN seqdesc_sequences s ON s.id = m.sequence_id;
    """)

def _coerce_bool_to_int(x: Any) -> Optional[int]:
    if x is None:
        return None
    if isinstance(x, bool):
        return 1 if x else 0
    # sometimes LLMs send "true"/"false"
    xs = str(x).strip().lower()
    if xs in {"true", "1", "yes"}:
        return 1
    if xs in {"false", "0", "no"}:
        return 0
    return None

def _coerce_float(x: Any) -> Optional[float]:
    try:
        return float(x) if x is not None else None
    except Exception:
        return None

def _extract_measure(obj: Any) -> Tuple[Optional[float], Optional[str]]:
    """obj like {"value": 50, "unit": "mM"} or None -> (50.0, 'mM')"""
    if isinstance(obj, dict):
        return _coerce_float(obj.get("value")), (obj.get("unit") if obj.get("unit") is not None else None)
    return None, None

def insert_seqdesc_object(
    db_path: Path | str,
    *,
    article_name: str,
    doi: Optional[str],
    model_name: str,
    sequence_descriptors: List[Tuple[str, Dict[str, Any]]],
    source_path: Optional[Path] = None,
) -> int:
    """Insert one 'run' of sequence descriptors and return run_id.

    The payload shape:
      {
        "<sequence_key>": {
          "is_seq": bool|None,
          "sequence_full": str|None,
          ...
          "modifications": [{"modification_position": int, "modification_type": str, "modification_description": str}, ...],
          "primers": {"forward": str|None, "reverse": str|None},
          "T": {"value": float, "unit": str}|None,
          ...
        },
        ...
      }
    """
    created_at = datetime.utcnow().strftime("%Y-%m-%dT%H:%M:%SZ")
    raw_json = json.dumps(sequence_descriptors, ensure_ascii=False)

    db_path = Path(db_path)
    db_path.parent.mkdir(parents=True, exist_ok=True)

    with closing(sqlite3.connect(str(db_path))) as conn:
        conn.execute("PRAGMA journal_mode = WAL;")
        _ensure_seqdesc_schema(conn)

        with conn:  # transaction
            cur = conn.cursor()
            cur.execute(
                """
                INSERT INTO seqdesc_runs(created_at, model_name, article_name, doi, source_path, raw_json)
                VALUES (?, ?, ?, ?, ?, ?)
                """,
                (
                    created_at,
                    model_name,
                    article_name,
                    doi,
                    str(source_path) if source_path else None,
                    raw_json,
                ),
            )
            run_id = cur.lastrowid

            for seq_key, payload in (sequence_descriptors or []):
                # For very sparse entries, payload can be {} — guard everything
                payload = payload or {}

                is_seq = _coerce_bool_to_int(payload.get("is_seq"))
                seq_full = payload.get("sequence_full")
                seq_norm = payload.get("sequence_normalized")
                seq_exp = payload.get("sequence_expanded")
                seq_bb = payload.get("sequence_backbone")
                seq_bb_exp = payload.get("sequence_backbone_expanded")
                fluor = payload.get("fluorophore")
                quen = payload.get("quencher")
                target_raw = payload.get("target_raw")
                target_norm = payload.get("target_normalized")

                primers = payload.get("primers") or {}
                primers_forward = primers.get("forward")
                primers_reverse = primers.get("reverse")

                pH_val = _coerce_float(payload.get("pH"))
                anneal_raw = payload.get("annealing_raw")

                T_val, T_unit = _extract_measure(payload.get("T"))
                Tris_val, Tris_unit = _extract_measure(payload.get("Tris"))
                Na_val, Na_unit = _extract_measure(payload.get("Na"))
                K_val, K_unit = _extract_measure(payload.get("K"))
                Mg_val, Mg_unit = _extract_measure(payload.get("Mg"))
                DMSO_val, DMSO_unit = _extract_measure(payload.get("DMSO"))

                outcome = _coerce_bool_to_int(payload.get("outcome"))

                cur.execute(
                    """
                    INSERT INTO seqdesc_sequences(
                        run_id, sequence_key, is_seq,
                        sequence_full, sequence_normalized, sequence_expanded,
                        sequence_backbone, sequence_backbone_expanded,
                        fluorophore, quencher,
                        target_raw, target_normalized,
                        primers_forward, primers_reverse,
                        pH, annealing_raw,
                        T_value, T_unit,
                        Tris_value, Tris_unit,
                        Na_value, Na_unit,
                        K_value, K_unit,
                        Mg_value, Mg_unit,
                        DMSO_value, DMSO_unit,
                        outcome, raw_json
                    )
                    VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                    """,
                    (
                        run_id, seq_key, is_seq,
                        seq_full, seq_norm, seq_exp,
                        seq_bb, seq_bb_exp,
                        fluor, quen,
                        target_raw, target_norm,
                        primers_forward, primers_reverse,
                        pH_val, anneal_raw,
                        T_val, T_unit,
                        Tris_val, Tris_unit,
                        Na_val, Na_unit,
                        K_val, K_unit,
                        Mg_val, Mg_unit,
                        DMSO_val, DMSO_unit,
                        outcome, json.dumps(payload, ensure_ascii=False),
                    ),
                )
                sequence_id = cur.lastrowid

                # Modifications (array of objects)
                for m in payload.get("modifications") or []:
                    if not isinstance(m, dict):
                        continue
                    cur.execute(
                        """
                        INSERT INTO seqdesc_modifications(
                            sequence_id, modification_position, modification_type, modification_description
                        ) VALUES (?,?,?,?)
                        """,
                        (
                            sequence_id,
                            m.get("modification_position"),
                            m.get("modification_type"),
                            m.get("modification_description"),
                        ),
                    )

        return run_id
