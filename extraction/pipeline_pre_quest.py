# pipeline_filedriven.py
# -*- coding: utf-8 -*-
"""
File-driven multi-pass extractor with Outlines + Ollama.

- Reads config, prompts, and schemas from disk (Git-friendly).
- Runs A..F passes (configurable) with Outlines JSON-guided generation.
- Saves raw text (*.txt), pretty JSON (*.json), and errors (*.log), never overwriting.
- Stitches pass outputs into a full object, validates against full schema (if provided),
  and optionally inserts into SQLite via hyb_db.insert_article_object.

Requirements:
  pip install outlines ollama jsonschema tqdm

Usage (script):
  from pipeline_filedriven import run_project
  run_project("your_project_dir")

The project_dir must contain (by default):
  config/pipeline.json
  passes/<pass_name>/{schema.json,prompt.txt}
  schemas/full.json
  inputs/*.txt
"""

import json
import logging
import re
import os, sys
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import ollama
import outlines
from jsonschema import Draft202012Validator
from outlines.types import JsonSchema
from tqdm import tqdm

API_TOKEN = os.getenv("OPEN_BUTTON_TOKEN", None)


# ──────────────────────────────────────────────────────────────────────
# Config models
# ──────────────────────────────────────────────────────────────────────


@dataclass
class PassConfig:
    """Single extraction pass config loaded from pipeline.json."""

    name: str  # e.g., "A_core"
    schema_path: Path  # path to JSON Schema file
    prompt_path: Path  # path to the prompt .txt file
    timeout: Optional[int]


@dataclass
class PipelineConfig:
    """Pipeline config loaded from config/pipeline.json."""

    model_names: List[str]
    ollama_parameters: Dict[str, Any]
    ollama_base_url: str
    timeout_s: Optional[int]
    input_dir: Path
    out_dir: Path
    full_schema_path: Optional[Path]
    common_prompt_path: Optional[Path]
    construct_single_experiment_passes: List[PassConfig]
    db_path: Optional[Path]
    article_glob: str
    pre_passes: List[PassConfig]
    passes: List[PassConfig]


def model_name_encode(model_name: str) -> str:
    return model_name.replace("/", "_").replace("\\", "_").replace(":", "_")


def load_pipeline_config(project_dir: Path) -> PipelineConfig:
    """Load pipeline.json and construct a PipelineConfig.

    Expected JSON structure in config/pipeline.json:
    {
      "model_name": "myaniu/qwen2.5-1m:7b",
      "num_ctx": 131072,
      "num_predict": 65536,
      "timeout_s": 1800,
      "input_dir": "inputs",
      "out_dir": "out",
      "full_schema_path": "schemas/full.json",
      "db_path": "out/massive.sqlite",   // or null to skip DB
      "article_glob": "*.txt",
      "passes": [
        {"name": "A_core",      "schema": "passes/A_core/schema.json",      "prompt": "passes/A_core/prompt.txt"},
        {"name": "B_index",     "schema": "passes/B_index/schema.json",     "prompt": "passes/B_index/prompt.txt"},
        {"name": "C_sequences", "schema": "passes/C_sequences/schema.json", "prompt": "passes/C_sequences/prompt.txt"},
        {"name": "D_parameters","schema": "passes/D_parameters/schema.json","prompt": "passes/D_parameters/prompt.txt"},
        {"name": "E_outcomes",  "schema": "passes/E_outcomes/schema.json",  "prompt": "passes/E_outcomes/prompt.txt"},
        {"name": "F_pairings",  "schema": "passes/F_pairings/schema.json",  "prompt": "passes/F_pairings/prompt.txt"}
      ]
    }
    """
    cfg_path = project_dir / "config" / "pipeline.json"
    data = json.loads(cfg_path.read_text(encoding="utf-8"))

    def _opt_path(p) -> Optional[Path]:
        return (project_dir / p) if p else None

    pre_passes: List[PassConfig] = []
    for p in data["pre_passes"]:
        pre_passes.append(
            PassConfig(
                name=p["name"],
                schema_path=project_dir / p["schema"],
                prompt_path=project_dir / p["prompt"],
                timeout=p.get("timeout", None),
            )
        )

    construct_single_experiment_passes = []
    for p in data["construct_single_experiment_passes"]:
        construct_single_experiment_passes.append(
            PassConfig(
                name=p["name"],
                schema_path=project_dir / p["schema"],
                prompt_path=project_dir / p["prompt"],
                timeout=p.get("timeout", None),
            )
        )

    passes: List[PassConfig] = []
    for p in data["passes"]:
        passes.append(
            PassConfig(
                name=p["name"],
                schema_path=project_dir / p["schema"],
                prompt_path=project_dir / p["prompt"],
                timeout=p.get("timeout", None),
            )
        )

    return PipelineConfig(
        model_names=list(data.get("model_names", [])),
        ollama_parameters=dict(data.get("ollama_parameters", {})),
        ollama_base_url=str(data.get("ollama_base_url", None)),
        timeout_s=int(data.get("timeout_s", None)),
        input_dir=project_dir / data.get("input_dir", "inputs"),
        out_dir=project_dir / data.get("out_dir", "out"),
        full_schema_path=_opt_path(data.get("full_schema_path")),
        common_prompt_path=_opt_path(data.get("common_prompt_path")),
        construct_single_experiment_passes=construct_single_experiment_passes,
        db_path=_opt_path(data.get("db_path")),
        article_glob=data.get("article_glob", "*.txt"),
        pre_passes=pre_passes,
        passes=passes,
    )


# ──────────────────────────────────────────────────────────────────────
# Logging
# ──────────────────────────────────────────────────────────────────────


class TqdmLoggingHandler(logging.Handler):
    def emit(self, record):
        try:
            msg = self.format(record)
            tqdm.write(msg)
            self.flush()
        except (KeyboardInterrupt, SystemExit):
            raise
        except:
            self.handleError(record)


def _make_logger(log_dir: Path) -> logging.Logger:
    log_dir.mkdir(parents=True, exist_ok=True)
    logger = logging.getLogger("pipeline_filedriven")
    logger.setLevel(logging.INFO)
    logger.handlers.clear()
     

    #ch = logging.StreamHandler(sys.stdout)
    ch = TqdmLoggingHandler()
    ch.setLevel(logging.INFO)
    ch.setFormatter(logging.Formatter("%(asctime)s | %(levelname)s | %(message)s"))
    logger.addHandler(ch)

    fh = logging.FileHandler(log_dir / "pipeline.log", encoding="utf-8")
    fh.setLevel(logging.INFO)
    fh.setFormatter(logging.Formatter("%(asctime)s | %(levelname)s | %(message)s"))
    logger.addHandler(fh)
    return logger


# ──────────────────────────────────────────────────────────────────────
# Tools (Ollama helpers) — Google-style docstrings
# ──────────────────────────────────────────────────────────────────────


def to_si(
    value: Optional[float], unit: Optional[str]
) -> Tuple[Optional[float], Optional[str]]:
    """Convert a numeric value and unit to SI.

    Supports temperature and common concentrations.

    Args:
      value: Parsed numeric value or None.
      unit: Unit string as written (e.g., '°C', 'mM', 'µM', 'nM', '%', 'K').

    Returns:
      A pair (si_value, si_unit), or (None, None) if unknown.
    """
    if value is None or unit is None:
        return None, None
    u = unit.strip().lower().replace("µ", "u")
    if u in {"c", "°c", "deg c", "celsius"}:
        return value + 273.15, "K"
    if u in {"k", "kelvin"}:
        return value, "K"
    if u in {"m", "mol/l"}:
        return value * 1000.0, "mol/m^3"
    if u in {"mm", "mmol/l", "mmol", "mm"}:
        return value * 1.0, "mol/m^3"
    if u in {"um", "umol/l", "µm", "µmol/l", "micromolar"}:
        return value * 1e-3, "mol/m^3"
    if u in {"nm", "nmol/l", "nanomolar"}:
        return value * 1e-6, "mol/m^3"
    if u in {"%", "percent", "perc"}:
        return value / 100.0, "dimensionless"
    return None, None


OLIGO_RE = re.compile(
    r"^\s*(?:(?P<prime>(?:5|3)(?:['′’]|0|O)?)\s*-\s*)?(?:(?P<prefix>(?:[A-Za-z0-9+]+-)+))?"
    r"(?P<seq>[ACGUTRYSWKMBDHVN]+)(?:(?P<suffix>(?:-[A-Za-z0-9+]+)+))?"
    r"(?:\s*\(\s*(?P<len>\d+)\s*(?:b|bp)\s*\)\s*)?\s*$",
    re.X,
)


def parse_oligo(raw: Optional[str]) -> Dict[str, Any]:
    """Parse a decorated oligo string into structured parts.

    Args:
      raw: The exact oligo string from the article (may include labels and length).

    Returns:
      A dict with keys: raw, sequence, length_bases, prime_prefix,
      five_prime_label, three_prime_label, labels, sense_antisense (None).
    """
    result = {
        "raw": raw,
        "sequence": None,
        "length_bases": None,
        "prime_prefix": None,
        "five_prime_label": None,
        "three_prime_label": None,
        "labels": [],
        "sense_antisense": None,
    }
    if not raw:
        return result
    m = OLIGO_RE.match(raw)
    if not m:
        return result
    prime = m.group("prime")
    if prime:
        result["prime_prefix"] = 5 if prime.startswith("5") else 3
    seq = m.group("seq")
    if seq:
        result["sequence"] = seq.upper()
    if m.group("len"):
        result["length_bases"] = int(m.group("len"))
    labels: List[str] = []
    if m.group("prefix"):
        labels += [x for x in m.group("prefix").split("-") if x]
    if m.group("suffix"):
        labels += [x for x in m.group("suffix").split("-") if x]
    result["labels"] = labels
    if labels:
        result["five_prime_label"] = labels[0]
        result["three_prime_label"] = labels[-1]
    return result


def make_measurement(
    raw: Optional[str], value: Optional[float] = None, unit: Optional[str] = None
) -> Dict[str, Any]:
    """Build a 'measurement' object with SI conversion.

    Args:
      raw: Raw textual measurement (e.g., '58 °C', '2 mM', '10%').
      value: Parsed numeric value, if available.
      unit: Unit string as written.

    Returns:
      A dict with keys: raw, value, unit, si_value, si_unit, assumptions (None).
    """
    si_value, si_unit = (
        to_si(value, unit) if (value is not None and unit is not None) else (None, None)
    )
    return {
        "raw": raw or "",
        "value": value,
        "unit": unit,
        "si_value": si_value,
        "si_unit": si_unit,
        "assumptions": None,
    }


# ──────────────────────────────────────────────────────────────────────
# JSON helpers
# ──────────────────────────────────────────────────────────────────────


def repair_json(text: str) -> str:
    """Best-effort JSON repair for streamed outputs."""
    start = text.find("{")
    end = text.rfind("}")
    if start == -1 or end == -1 or end <= start:
        return text
    candidate = text[start : end + 1]
    try:
        json.loads(candidate)
        return candidate
    except Exception:
        candidate = re.sub(r",\s*([}\]])", r"\1", candidate)
        json.loads(candidate)
        return candidate


# ──────────────────────────────────────────────────────────────────────
# Outlines runner
# ──────────────────────────────────────────────────────────────────────


def _now_stamp() -> str:
    return datetime.utcnow().strftime("%Y%m%d_%H%M%S_%f")


def think_generate(
    model: outlines.models.ollama.Ollama,
    model_input: outlines.inputs.Chat | str | list,
    logger: logging.Logger,
    output_type: Optional[Any] = None,
    think: bool = True,
    **kwargs: Any,
) -> str:
    if think:
        try:
            logger.debug("Trying thinking mode")
            response = model.generate(
                model_input=model_input, output_type=output_type, think=True, **kwargs
            )
            return response
        except ollama.ResponseError:
            logger.warning(
                f"Seems that model {model.model_name} does not support thinking."
            )

    logger.debug("Trying non-thinking mode")
    response = model.generate(
        model_input=model_input, output_type=output_type, think=False, **kwargs
    )

    return response


def run_single_pass(
    model: Any,
    article_text: str,
    pass_cfg: PassConfig,
    out_base: Path,
    article_stem: str,
    tools: List[Any],
    logger: logging.Logger,
    ollama_parameters: Dict[str, Any],
    model_name: str,
) -> Dict[str, Any]:
    """Run one pass (schema+prompt from files), save raw+json+log, return object."""
    txt_dir = out_base / "txt"
    json_dir = out_base / "json"
    log_dir = out_base / "logs"
    for d in (txt_dir, json_dir, log_dir):
        d.mkdir(parents=True, exist_ok=True)

    js = JsonSchema(pass_cfg.schema_path.read_text(encoding="utf-8"))
    validator = Draft202012Validator(json.loads(js.schema))
    prompt = pass_cfg.prompt_path.read_text(encoding="utf-8")

    stamp = _now_stamp()
    raw_txt_path = (
        txt_dir
        / f"{article_stem}__{pass_cfg.name}__{model_name_encode(model_name)}__{stamp}.txt"
    )
    json_out_path = (
        json_dir
        / f"{article_stem}__{pass_cfg.name}__{model_name_encode(model_name)}__{stamp}.json"
    )
    err_log_path = (
        log_dir
        / f"{article_stem}__{pass_cfg.name}__{model_name_encode(model_name)}__{stamp}.log"
    )

    logger.info(f"[{pass_cfg.name}:{model_name}] generating …")
    response = ""
    try:
        response = think_generate(
            model=model,
            model_input=prompt
            + "\n"
            + "And here is the article text you must base your answer on:\n\n<article>\n"
            + article_text
            + "\n<\\article>\n",
            output_type=js,
            options=ollama_parameters,
            logger=logger,
            keep_alive="30s",
        )
    except Exception as e:
        logger.exception(f"[{pass_cfg.name}:{model_name}] stream error")
        err_log_path.write_text(f"STREAM ERROR:\n{e}\n", encoding="utf-8")
        raise

    raw_txt_path.write_text(response, encoding="utf-8")

    try:
        fixed = repair_json(response)
        obj = json.loads(fixed)
    except Exception as e:
        logger.exception(f"[{pass_cfg.name}:{model_name}] JSON parse error")
        err_log_path.write_text(
            f"JSON ERROR:\n{e}\nRAW:\n{response}\n", encoding="utf-8"
        )
        raise

    errors = sorted(validator.iter_errors(obj), key=lambda er: er.path)
    if errors:
        msg = "\n".join(str(e) for e in errors)
        logger.error(f"[{pass_cfg.name}:{model_name}] validation errors:\n{msg}")
        err_log_path.write_text(
            f"VALIDATION ERRORS:\n{msg}\nJSON:\n{json.dumps(obj, indent=2)}",
            encoding="utf-8",
        )
    else:
        logger.info(f"[{pass_cfg.name}:{model_name}] validation OK")
        logger.info(f"[{pass_cfg.name}] validation OK [{model_name}]")

    json_out_path.write_text(
        json.dumps(obj, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    return obj


def run_construct_single_experiment_pass(
    model: Any,
    article_text: str,
    sequence: str,
    sequence_id: int,
    pass_cfg: PassConfig,
    out_base: Path,
    article_stem: str,
    tools: List[Any],
    logger: logging.Logger,
    ollama_parameters: Dict[str, Any],
    model_name: str,
) -> Dict[str, Any]:
    """Run one pass (schema+prompt from files), save raw+json+log, return object."""
    txt_dir = out_base / "txt"
    json_dir = out_base / "json"
    log_dir = out_base / "logs"
    for d in (txt_dir, json_dir, log_dir):
        d.mkdir(parents=True, exist_ok=True)

    js = JsonSchema(pass_cfg.schema_path.read_text(encoding="utf-8"))
    validator = Draft202012Validator(json.loads(js.schema))
    prompt = pass_cfg.prompt_path.read_text(encoding="utf-8")

    stamp = _now_stamp()
    raw_txt_path = (
        txt_dir
        / f"{article_stem}__{pass_cfg.name}__{sequence_id}__{model_name_encode(model_name)}__{stamp}.txt"
    )
    json_out_path = (
        json_dir
        / f"{article_stem}__{pass_cfg.name}__{sequence_id}__{model_name_encode(model_name)}__{stamp}.json"
    )
    err_log_path = (
        log_dir
        / f"{article_stem}__{pass_cfg.name}__{sequence_id}__{model_name_encode(model_name)}__{stamp}.log"
    )

    logger.info(f"[{pass_cfg.name}:{model_name}] generating …")
    response = ""
    try:
        response = think_generate(
            model=model,
            model_input=outlines.inputs.Chat(
                [
                    {
                        "role": "system",
                        "content": prompt
                        + "\n"
                        + "And here is the article text you must base your answer on:\n\n<article>\n"
                        + article_text
                        + "\n<\\article>\n",
                    },
                    {
                        "role": "user",
                        "content": "Let's describe a single nucleotide sequence!",
                    },
                    {
                        "role": "assistant",
                        "content": "Sure! Let's describe one! But before we start, could you please tell me in which format you would like me to provide you an answer?",
                    },
                    {
                        "role": "user",
                        "content": "Great question! I would like your answer to satisfy the following JSON schema:\n```json"
                        + js.schema
                        + "\n```\n\nIs it OK?",
                    },
                    {
                        "role": "assistant",
                        "content": "Absolutely! Now please provide the nucleotide sequence you want me to describe in terms of tthe hybridization experiment design and I will provide you its description strictly following your provided JSON schema!",
                    },
                    {
                        "role": "user",
                        "content": sequence,
                    },
                ]
            ),
            logger=logger,
            output_type=js,
            options=ollama_parameters,
            keep_alive="30s",
        )

    except Exception as e:
        logger.exception(f"[{pass_cfg.name}:{model_name}] stream error")
        err_log_path.write_text(f"STREAM ERROR:\n{e}\n", encoding="utf-8")
        raise

    raw_txt_path.write_text(response, encoding="utf-8")

    try:
        fixed = repair_json(response)
        obj = json.loads(fixed)
    except Exception as e:
        logger.exception(f"[{pass_cfg.name}:{model_name}] JSON parse error")
        err_log_path.write_text(
            f"JSON ERROR:\n{e}\nRAW:\n{response}\n", encoding="utf-8"
        )
        raise

    errors = sorted(validator.iter_errors(obj), key=lambda er: er.path)
    if errors:
        msg = "\n".join(str(e) for e in errors)
        logger.error(f"[{pass_cfg.name}:{model_name}] validation errors:\n{msg}")
        err_log_path.write_text(
            f"VALIDATION ERRORS:\n{msg}\nJSON:\n{json.dumps(obj, indent=2)}",
            encoding="utf-8",
        )
    else:
        logger.info(f"[{pass_cfg.name}:{model_name}] validation OK")
        logger.info(f"[{pass_cfg.name}] validation OK [{model_name}]")

    json_out_path.write_text(
        json.dumps(obj, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    return obj


def run_query_model(
    model: Any,
    article_text: str,
    sequences: List[str],
    out_base: Path,
    article_stem: str,
    common_prompt_path: Path,
    logger: logging.Logger,
    ollama_parameters: Dict[str, Any],
    model_name: str,
    tqdm_position: int = 0,
) -> Dict[str, Any]:
    """Run one pass (schema+prompt from files), save raw+json+log, return object."""
    pass_name = "query_chat"
    txt_dir = out_base / "txt"
    json_dir = out_base / "json"
    log_dir = out_base / "logs"
    for d in (txt_dir, json_dir, log_dir):
        d.mkdir(parents=True, exist_ok=True)

    prompt = common_prompt_path.read_text(encoding="utf-8")

    stamp = _now_stamp()
    raw_txt_path = (
        txt_dir
        / f"{article_stem}__{pass_name}__{model_name_encode(model_name)}__{stamp}.txt"
    )
    json_log_path = (
        json_dir
        / f"{article_stem}__{pass_name}__{model_name_encode(model_name)}__{stamp}.log.json"
    )
    json_out_path = (
        json_dir
        / f"{article_stem}__{pass_name}__{model_name_encode(model_name)}__{stamp}.json"
    )
    err_log_path = (
        log_dir
        / f"{article_stem}__{pass_name}__{model_name_encode(model_name)}__{stamp}.log"
    )

    logger.info(f"[{pass_name}:{model_name}] generating …")

    def ask_with_schema(chat_messages: outlines.inputs.Chat, schema: JsonSchema):
        response = ""
        try:
            response = think_generate(
                model=model,
                model_input=chat_messages,
                output_type=schema,
                options=ollama_parameters,
                logger=logger,
                keep_alive="30s",
            )
        except Exception as e:
            logger.exception(f"[{pass_name}:{model_name}] stream error")
            err_log_path.write_text(f"STREAM ERROR:\n{e}\n", encoding="utf-8")
            raise

        with open(raw_txt_path, mode="at", encoding="utf-8") as f:
            f.write(f"> {chat_messages.messages[-1]}\n< ")
            f.write(response)
            f.write("\n\n")

        try:
            fixed = repair_json(response)
            obj = json.loads(fixed)
        except Exception as e:
            logger.exception(f"[{pass_name}:{model_name}] JSON parse error")
            err_log_path.write_text(
                f"JSON ERROR:\n{e}\nRAW:\n{response}\n", encoding="utf-8"
            )
            raise

        return obj, response

    base_chat = outlines.inputs.Chat(
        [
            {
                "role": "system",
                "content": prompt
                + "\n"
                + "And here is the article text you must base your answer on:\n\n<article>\n"
                + article_text
                + "\n<\\article>\n",
            }
        ]
    )
    answers = []

    try:

        def parse_sequence(seq: str, base_chat: outlines.inputs.Chat):
            chat = outlines.inputs.Chat(base_chat.messages)
            questions_to_schema: List[Tuple[str, str, Dict[str, Any]]] = [
                (
                    "is_seq",
                    "Check the whole article text. Is your picked sequence really a probe sequence or a part of probe sequence in this article text? Put true here if and only if this sequence is being described and presented as a hybridization probe. If that's a random abbreviation or nucleotide-looking string which is not a hybridization probe or otherwise not a hybridization probe, put false here.",
                    {"type": "boolean"},
                ),
                (
                    "sequence_full",
                    "Provide this sequence fully as a probe sequence in IUPAC-normalized format: from 5' to 3' end, with fluorophore and quencher. Use capital Latin letters, digits and dashes, you may also use parentheses and apostrophy. Put null here if not applicable.",
                    {
                        "type": ["string", "null"],
                        "minLength": 5,
                        "maxLength": 150,
                    },
                ),
                (
                    "sequence_normalized",
                    "Provide this probe sequence in IUPAC-normalized format: from 5' to 3' end, with fluorophore and quencher. Use capital Latin letters, digits and dashes, you may also use parentheses and apostrophy. Put null here if not applicable.",
                    {
                        "type": ["string", "null"],
                        "minLength": 5,
                        "maxLength": 150,
                        "pattern": r"^5'-([a-zA-Z0-9(_)'-]*-)?([a-zA-Z0-9()']*?[ACGUTRYSWKMBDHVN]{5,}[a-zA-Z0-9()']*?)(-[a-zA-Z0-9(_)'-]*)?-3'$",
                    },
                ),
                (
                    "sequence_expanded",
                    "Provide this probe sequence in expanded IUPAC format (with all repeats expanded and no parentheses in the probe sequence backbone body): from 5' to 3' end, with fluorophore and quencher. Use capital Latin letters, digits and dashes, you may also use parentheses and apostrophy. Put null here if not applicable.",
                    {
                        "type": ["string", "null"],
                        "minLength": 5,
                        "maxLength": 150,
                        "pattern": r"^5'-([a-zA-Z0-9_'-]*-)?([a-zA-Z0-9']*?[ACGUTRYSWKMBDHVN]{5,}[a-zA-Z0-9']*?)(-[a-zA-Z0-9_'-]*)?-3'$",
                    },
                ),
                (
                    "sequence_backbone",
                    "Now provide only the probe sequence body from 5' to 3', without any fluorophores, modifications and quenchers. Use capital Latin letters, digits and dashes, you may also use parentheses and apostrophy. Put null here if not applicable.",
                    {
                        "type": ["string", "null"],
                        "minLength": 5,
                        "maxLength": 150,
                        "pattern": r"^5'-([ACGUTRYSWKMBDHVN0-9()]{5,})-3'$",
                    },
                ),
                (
                    "sequence_backbone_expanded",
                    "Now provide only the expanded probe sequence body from 5' to 3' with all repeats expanded, without any fluorophores, modifications and quenchers. Use capital Latin letters, digits, dashes and apostrophy. Only the expanded backbone of probe sequence body. Put null here if not applicable.",
                    {
                        "type": ["string", "null"],
                        "minLength": 5,
                        "maxLength": 150,
                        "pattern": r"^5'-([ACGUTRYSWKMBDHVN0-9]{5,})-3'$",
                    },
                ),
                (
                    "fluorophore",
                    "Provide the fluorophore of this probe. Use capital Latin letters, digits and dashes, you may also use an apostrophy. Put null here if not applicable or not present in the text of the article.",
                    {
                        "type": ["string", "null"],
                        "minLength": 3,
                        "maxLength": 150,
                        "pattern": r"^[A-Z0-9']{3,}$",
                    },
                ),
                (
                    "quencher",
                    "Provide the quencher of this probe. Use capital Latin letters, digits and dashes, you may also use an apostrophy. Put null here if not applicable or not present in the text of the article.",
                    {
                        "type": ["string", "null"],
                        "minLength": 3,
                        "maxLength": 150,
                        "pattern": r"^[A-Z0-9']{3,}$",
                    },
                ),
                (
                    "modifications",
                    "Now provide the modifications of the probe sequence as an array, where each element is a modification and its position in 5'-3' direction. Use Latin letters, digits and dashes, you may also use parentheses and apostrophy. Provide an empty array if not present in the article text.",
                    {
                        "type": "array",
                        "minItems": 0,
                        "maxItems": 150,
                        "items": {
                            "type": "object",
                            "additionalProperties": False,
                            "required": [
                                "modification_position",
                                "modification_type",
                                "modification_description",
                            ],
                            "properties": {
                                "modification_position": {
                                    "type": "integer",
                                    "minimum": 1,
                                },
                                "modification_type": {
                                    "type": "string",
                                    "maxLength": 100,
                                    "minLength": 1,
                                },
                                "modification_description": {
                                    "type": "string",
                                    "minLength": 1,
                                    "maxLength": 150,
                                },
                            },
                        },
                    },
                ),
                (
                    "target_raw",
                    "Describe the target to which this probe was designed to hybridize.",
                    {"type": "string", "minLength": 5, "maxLength": 250},
                ),
                (
                    "target_normalized",
                    "Now provide the target sequence to which this probe should hybridize, from 5' to 3'. Use capital Latin letters, digits and dashes, you may also use parentheses and apostrophy. Put null here if not applicable or if the exact sequence is not present in the article text.",
                    {
                        "type": ["string", "null"],
                        "minLength": 5,
                        "maxLength": 150,
                        "pattern": r"^(5')?([A-Z0-9_()'-]*)[-]?([ACGUTRYSWKMBDHVN0-9()]{5,})[-]?([A-Z0-9_()-]*)(3')?$",
                    },
                ),
                (
                    "primers",
                    "Describe the primer sequences in IUPAC-normalized format, each from 5' to 3' end. Use capital Latin letters, digits and dashes, parentheses and apostrophy. Put null to the primer if it is not present in the article text.",
                    {
                        "type": "object",
                        "additionalProperties": False,
                        "required": ["forward", "reverse"],
                        "properties": {
                            "forward": {
                                "type": ["string", "null"],
                                "minLength": 5,
                                "maxLength": 150,
                                "pattern": r"^(5')?([A-Z0-9_()'-]*)[-]?([ACGUTRYSWKMBDHVN0-9()]{5,})[-]?([A-Z0-9_()-]*)(3')?$",
                            },
                            "reverse": {
                                "type": ["string", "null"],
                                "minLength": 5,
                                "maxLength": 150,
                                "pattern": r"^(5')?([A-Z0-9_()'-]*)[-]?([ACGUTRYSWKMBDHVN0-9()]{5,})[-]?([A-Z0-9_()-]*)(3')?$",
                            },
                        },
                    },
                ),
                (
                    "pH",
                    "Describe the pH in this experiment. Only put null here if this information is not present in the article text and can't be inferred from the whole article text.",
                    {"type": ["number", "null"]},
                ),
                (
                    "annealing_raw",
                    "Describe the annealing in this experiment. Provide the raw description string. If that's can't be inferred from the whole article text, explain why.",
                    {"type": ["string"], "minLength": 10, "maxLength": 250},
                ),
                (
                    "T",
                    "Describe the melting temperature in this experiment and provide the measurement unit. Only put null here if this information is not present in the article text and can't be inferred from the whole article text.",
                    {
                        "type": ["object", "null"],
                        "additionalProperties": False,
                        "required": ["value", "unit"],
                        "properties": {
                            "value": {"type": "number"},
                            "unit": {"type": "string", "minLength": 1, "maxLength": 10},
                        },
                    },
                ),
                (
                    "Tris",
                    "Describe the amount of Tris in this experiment and provide the measurement unit. Only put null here if this information is not present in the article text and can't be inferred from the whole article text.",
                    {
                        "type": ["object", "null"],
                        "additionalProperties": False,
                        "required": ["value", "unit"],
                        "properties": {
                            "value": {"type": "number"},
                            "unit": {"type": "string", "minLength": 1, "maxLength": 10},
                        },
                    },
                ),
                (
                    "Na",
                    "Describe the amount of Na (Sodium) in this experiment and provide the measurement unit. Only put null here if this information is not present in the article text and can't be inferred from the whole article text.",
                    {
                        "type": ["object", "null"],
                        "additionalProperties": False,
                        "required": ["value", "unit"],
                        "properties": {
                            "value": {"type": "number"},
                            "unit": {"type": "string", "minLength": 1, "maxLength": 10},
                        },
                    },
                ),
                (
                    "K",
                    "Describe the amount of K (Potassium) in this experiment and provide the measurement unit. Only put null here if this information is not present in the article text and can't be inferred from the whole article text.",
                    {
                        "type": ["object", "null"],
                        "additionalProperties": False,
                        "required": ["value", "unit"],
                        "properties": {
                            "value": {"type": "number"},
                            "unit": {"type": "string", "minLength": 1, "maxLength": 10},
                        },
                    },
                ),
                (
                    "Mg",
                    "Describe the amount of Mg (Magnesium) in this experiment and provide the measurement unit. Only put null here if this information is not present in the article text and can't be inferred from the whole article text.",
                    {
                        "type": ["object", "null"],
                        "additionalProperties": False,
                        "required": ["value", "unit"],
                        "properties": {
                            "value": {"type": "number"},
                            "unit": {"type": "string", "minLength": 1, "maxLength": 10},
                        },
                    },
                ),
                (
                    "DMSO",
                    "Describe the amount of DMSO in this experiment and provide the measurement unit. Only put null here if this information is not present in the article text and can't be inferred from the whole article text.",
                    {
                        "type": ["object", "null"],
                        "additionalProperties": False,
                        "required": ["value", "unit"],
                        "properties": {
                            "value": {"type": "number"},
                            "unit": {"type": "string", "minLength": 1, "maxLength": 10},
                        },
                    },
                ),
                (
                    "outcome",
                    "Describe the outcome of this hybridization experiment based on the article text. Put true in case of successful hybridization of this probe to target, put false in case of unsuccessful and put null if this information is not present in the article.",
                    {"type": ["boolean", "null"]},
                ),
            ]

            seq_desc: Dict[str, Any] = dict()

            for param, query, schema in tqdm(
                questions_to_schema,
                desc="Questions to the sequence",
                position=tqdm_position + 1,
                leave=False,
            ):
                try:
                    chat.add_user_message(
                        query
                        + "\nAnd here is the schema yout answer has to follow:\n```json\n"
                        + json.dumps(schema)
                        + "```\n"
                    )
                    response, raw = ask_with_schema(
                        chat_messages=chat, schema=JsonSchema(schema)
                    )
                    answers.append({"seq": seq, "param": param, "response": response})
                    seq_desc[param] = response
                    chat.add_assistant_message(raw)
                except Exception as e:
                    logger.exception(
                        f"Exception on sequence {seq} during query: {query}", e
                    )
            return seq_desc

        described_sequences: Dict[str, Dict[str, Any]] = dict()
        for seq in tqdm(
            sequences, desc="Found sequences", position=tqdm_position, leave=False
        ):
            base_chat_with_sequence = outlines.inputs.Chat(base_chat.messages)
            base_chat_with_sequence.add_user_message(
                "Let's pick and analyze a single probe sequence from the article text. Provide the probe sequence which we will describe in all the following messages."
            )
            base_chat_with_sequence.add_assistant_message(seq)
            base_chat_with_sequence.add_user_message(
                f"Great choice! Let's analyze nucleotidic sequence {seq} for the rest of this chat!"
            )
            try:
                sequence_descriptor = parse_sequence(
                    seq, base_chat=base_chat_with_sequence
                )
                described_sequences[seq] = sequence_descriptor
                answers.append(
                    {"sequence": seq, "sequence_descriptor": sequence_descriptor}
                )
            except Exception as e:
                logger.exception(
                    f'[{pass_name}:{model_name}] Sequence "{seq}" error computing descriptor'
                )
                err_log_path.write_text(
                    f'[{pass_name}:{model_name}] Sequence "{seq}" error computing descriptor',
                    encoding="utf-8",
                )
    finally:
        json_log_path.write_text(
            json.dumps(answers, indent=2, ensure_ascii=False), encoding="utf-8"
        )
        json_out_path.write_text(
            json.dumps(described_sequences, indent=2, ensure_ascii=False),
            encoding="utf-8",
        )
    return described_sequences


# ──────────────────────────────────────────────────────────────────────
# Stitcher (to your full object)
# ──────────────────────────────────────────────────────────────────────


def _merge_reports(*reports: Optional[Dict[str, Any]]) -> Dict[str, Any]:
    out = {"missing": [], "uncertain": [], "notes": None}
    notes = []
    for r in reports:
        if not r:
            continue
        out["missing"].extend(r.get("missing") or [])
        out["uncertain"].extend(r.get("uncertain") or [])
        if r.get("notes"):
            notes.append(str(r["notes"]))
    out["missing"] = list(dict.fromkeys(out["missing"]))
    out["uncertain"] = list(dict.fromkeys(out["uncertain"]))
    out["notes"] = " | ".join(notes) if notes else None
    return out


def _to_si(
    value: Optional[float], unit: Optional[str]
) -> Tuple[Optional[float], Optional[str]]:
    return to_si(value, unit)


def _to_measurement_full(m_lite: Optional[Dict[str, Any]]) -> Optional[Dict[str, Any]]:
    if not m_lite:
        return None
    raw = m_lite.get("raw") or ""
    value = m_lite.get("value")
    unit = m_lite.get("unit")
    si_value, si_unit = (
        _to_si(value, unit)
        if (value is not None and unit is not None)
        else (None, None)
    )
    return {
        "raw": raw,
        "value": value,
        "unit": unit,
        "si_value": si_value,
        "si_unit": si_unit,
        "assumptions": None,
        "provenance": {
            "source_type": "unknown",
            "page": None,
            "section": None,
            "quote": None,
            "notes": None,
        },
    }


def _detect_sa_from_name(name: Optional[str]) -> Optional[str]:
    if not name:
        return None
    n = name.strip().lower()
    if n.endswith(")as"):
        return "antisense"
    if n.endswith(")s"):
        return "sense"
    return None


def _coerce_sa(value: Optional[str], name: Optional[str]) -> Optional[str]:
    m = {
        "s": "sense",
        "as": "antisense",
        "sense": "sense",
        "antisense": "antisense",
        "+": "sense",
        "-": "antisense",
    }
    if value is None or (isinstance(value, str) and not value.strip()):
        return _detect_sa_from_name(name)
    return m.get(str(value).strip().lower(), _detect_sa_from_name(name))


def _to_oligo_full(ol: Optional[Dict[str, Any]]) -> Optional[Dict[str, Any]]:
    if not ol:
        return None
    return {
        "raw": ol.get("raw") or "",
        "sequence": ol.get("sequence"),
        "length_bases": ol.get("length_bases"),
        "prime_prefix": None,
        "five_prime_label": ol.get("five_prime_label"),
        "three_prime_label": ol.get("three_prime_label"),
        "labels": [],
        "sense_antisense": ol.get("sense_antisense"),
        "provenance": {
            "source_type": "unknown",
            "page": None,
            "section": None,
            "quote": None,
            "notes": None,
        },
    }


def stitch_full(
    A_core: Dict[str, Any],
    B_index: Dict[str, Any],
    C_sequences: Dict[str, Any],
    D_parameters: Dict[str, Any],
    E_outcomes: Dict[str, Any],
    F_pairings: Dict[str, Any],
) -> Dict[str, Any]:
    core = {
        "doi": A_core.get("doi"),
        "abstract": A_core.get("abstract"),
        "topic": A_core.get("topic"),
    }
    E: Dict[str, Dict[str, Any]] = {}
    for e in B_index.get("experiments") or []:
        E[e["id_exp"]] = {
            "id_exp": e["id_exp"],
            "raw_description": e.get("raw_description"),
            "type": e.get("type"),
            "description": e.get("description"),
            "metadata": {},
            "sequences": {},
            "experiment_properties": {},
            "outcome": {},
            "pairing": {},
            "extraction_report": {"missing": [], "uncertain": [], "notes": None},
        }

    for item in C_sequences.get("items") or []:
        ie = item["id_exp"]
        if ie not in E:
            continue
        prb = item.get("probe") or {}
        seqs = {}
        seqs["probe"] = {
            "name": prb.get("name"),
            "amplicon_id": prb.get("amplicon_id"),
            "oligo": _to_oligo_full(prb.get("oligo")),
            "fluorophore": prb.get("fluorophore"),
            "quencher": prb.get("quencher"),
            "sense_antisense": _coerce_sa(prb.get("sense_antisense"), prb.get("name")),
            "notes": prb.get("notes"),
        }
        tgt = item.get("target_sequence")
        seqs["target_sequence"] = _to_oligo_full(tgt) if tgt is not None else None
        pr = item.get("primer_sequences")
        if isinstance(pr, dict):
            seqs["primer_sequences"] = {
                "forward": _to_oligo_full(pr.get("forward")),
                "reverse": _to_oligo_full(pr.get("reverse")),
            }
        else:
            seqs["primer_sequences"] = None
        rels = []
        for rs in item.get("related_sequences") or []:
            rels.append(
                {
                    "related_sequence": _to_oligo_full(rs.get("related_sequence")),
                    "description": rs.get("description"),
                }
            )
        seqs["related_sequences"] = rels
        E[ie]["sequences"] = seqs

    for item in D_parameters.get("items") or []:
        ie = item["id_exp"]
        if ie not in E:
            continue
        MD: Dict[str, Any] = {}
        _md = item.get("metadata") or {}
        MD["organism"] = _md.get("organism")
        MD["technology"] = _md.get("technology")
        ann = _md.get("annealing")
        if ann is None:
            MD["annealing"] = None
        elif isinstance(ann, dict):
            MD["annealing"] = {
                "quantitative": _to_measurement_full(ann.get("quantitative")),
                "qualitative": ann.get("qualitative"),
            }
        else:
            MD["annealing"] = None
        MD["pH"] = _to_measurement_full(_md.get("pH"))
        ri = _md.get("rna_impurities")
        if ri is None:
            MD["rna_impurities"] = None
        elif isinstance(ri, dict):
            MD["rna_impurities"] = {
                "quantitative": _to_measurement_full(ri.get("quantitative")),
                "qualitative": ri.get("qualitative"),
            }
        else:
            MD["rna_impurities"] = None

        EP: Dict[str, Any] = {}
        concs = (item.get("experiment_properties") or {}).get("concentrations") or {}
        EP["concentrations"] = {
            "dna_rna_concentration": _to_measurement_full(
                concs.get("dna_rna_concentration")
            ),
            "concentration_SI": _to_measurement_full(concs.get("concentration_SI")),
        }
        pars = (item.get("experiment_properties") or {}).get("parameters_SI") or {}
        EP["parameters_SI"] = {
            "temperature": _to_measurement_full(pars.get("temperature")),
            "Tris": _to_measurement_full(pars.get("Tris")),
            "Na": _to_measurement_full(pars.get("Na")),
            "K": _to_measurement_full(pars.get("K")),
            "Mg": _to_measurement_full(pars.get("Mg")),
            "DMSO": _to_measurement_full(pars.get("DMSO")),
        }
        E[ie]["metadata"] = MD
        E[ie]["experiment_properties"] = EP

    for item in E_outcomes.get("items") or []:
        ie = item["id_exp"]
        if ie not in E:
            continue
        E[ie]["outcome"] = {
            "outcome": item.get("outcome"),
            "fluorescence": _to_measurement_full(item.get("fluorescence")),
            "comparative_notes": item.get("comparative_notes"),
        }

    for item in F_pairings.get("items") or []:
        ie = item["id_exp"]
        if ie not in E:
            continue
        E[ie]["pairing"] = {
            "paired_with_probe_name": item.get("paired_with_probe_name"),
            "relationship": item.get("relationship"),
        }

    full_report = _merge_reports(
        A_core.get("extraction_report"),
        B_index.get("extraction_report"),
        C_sequences.get("extraction_report"),
        D_parameters.get("extraction_report"),
        E_outcomes.get("extraction_report"),
        F_pairings.get("extraction_report"),
    )

    return {
        "doi": core["doi"],
        "abstract": core["abstract"],
        "topic": core["topic"],
        "experiments": list(E.values()),
        "extraction_report": full_report,
    }


def _deep_merge_keep_left(a, b):
    """Shallow-friendly deep merge: keep a's non-null scalars; use b if a is None.
    - Dicts: recurse.
    - Lists: concatenate (no dedup).
    - Scalars: prefer a unless a is None/empty, then b.
    """
    if a is None:
        return b
    if b is None:
        return a
    if isinstance(a, dict) and isinstance(b, dict):
        out = dict(a)
        for k, bv in b.items():
            av = out.get(k)
            out[k] = _deep_merge_keep_left(av, bv) if k in out else bv
        return out
    if isinstance(a, list) and isinstance(b, list):
        return a + b
    # prefer a unless it's falsy and b is truthy
    return a if a not in (None, "", []) else b


def aggregate_c_outputs(outputs: Dict[str, Dict[str, Any]]) -> Dict[str, Any]:
    """Build a consolidated C_sequences object from any of: C_sequences, C1_probe_core, C2_target_primers, C3_related."""
    # Start with single-pass C if present
    base = outputs.get("C_sequences") or {
        "items": [],
        "extraction_report": {"missing": [], "uncertain": [], "notes": None},
    }

    # Build item index by id_exp from base
    items_map: Dict[str, Dict[str, Any]] = {}
    for it in base.get("items", []):
        if not isinstance(it, dict) or "id_exp" not in it:
            continue
        items_map[it["id_exp"]] = dict(it)

    def _merge_from(pass_name: str, fields: List[str]):
        obj = outputs.get(pass_name)
        if not isinstance(obj, dict):
            return
        for it in obj.get("items", []):
            if not isinstance(it, dict) or "id_exp" not in it:
                continue
            ie = it["id_exp"]
            tgt = items_map.setdefault(ie, {"id_exp": ie})
            for f in fields:
                if f in it:
                    tgt[f] = _deep_merge_keep_left(tgt.get(f), it[f])

        # merge extraction report
        br = base.get("extraction_report") or {
            "missing": [],
            "uncertain": [],
            "notes": None,
        }
        er = obj.get("extraction_report") or {
            "missing": [],
            "uncertain": [],
            "notes": None,
        }
        br["missing"] = list((br.get("missing") or []) + (er.get("missing") or []))
        br["uncertain"] = list(
            (br.get("uncertain") or []) + (er.get("uncertain") or [])
        )
        br_notes = [n for n in [br.get("notes"), er.get("notes")] if n]
        br["notes"] = " | ".join(br_notes) if br_notes else None
        base["extraction_report"] = br

    # Merge micro-passes over base (C1, C2, C3). The field names match your schemas.
    _merge_from("C1_probe_core", ["probe"])
    _merge_from("C2_target_primers", ["target_sequence", "primer_sequences"])
    _merge_from("C3_related", ["related_sequences"])

    # Produce consolidated list
    merged_items = list(items_map.values())
    # Normalize: ensure all top-level keys exist for stitch_full
    for it in merged_items:
        it.setdefault("probe", None)
        it.setdefault("target_sequence", None)
        it.setdefault("primer_sequences", None)
        it.setdefault("related_sequences", [])

    return {
        "items": merged_items,
        "extraction_report": base.get("extraction_report")
        or {"missing": [], "uncertain": [], "notes": None},
    }


# ──────────────────────────────────────────────────────────────────────
# Project runner
# ──────────────────────────────────────────────────────────────────────


def run_project(project_dir: str | Path) -> None:
    """Run the pipeline as configured by files under project_dir."""
    project_dir = Path(project_dir)
    cfg = load_pipeline_config(project_dir)

    out_base = cfg.out_dir
    out_base.mkdir(parents=True, exist_ok=True)
    logger = _make_logger(out_base / "logs")

    headers = dict()
    if API_TOKEN is not None:
        headers["Authorization"] = f"Bearer {API_TOKEN}"

    # Ollama client + Outlines model
    client = ollama.Client(
        host=cfg.ollama_base_url, timeout=cfg.timeout_s, headers=headers
    )

    ollama_models = client.list()
    for model_name in tqdm(cfg.model_names, desc="LLM Models", position=0, leave=False):
        model = outlines.from_ollama(client, model_name)
        tools = [to_si, parse_oligo, make_measurement]

        # Optional full-schema validator
        full_validator = None
        if cfg.full_schema_path and cfg.full_schema_path.exists():
            try:
                full_schema_text = cfg.full_schema_path.read_text(encoding="utf-8")
                full_validator = Draft202012Validator(json.loads(full_schema_text))
                logger.info("Loaded full schema for final validation.")
            except Exception:
                logger.exception(
                    "Failed to load/parse full schema; proceeding without final validation."
                )

        logger.info(f"Article glob: {cfg.article_glob}")

        # Iterate input articles
        files = sorted(
            cfg.input_dir.glob(cfg.article_glob), key=lambda s: str(s).upper()
        )
        logger.info(f"Files: {files}")

        for art_path in tqdm(files, desc="Articles", position=1, leave=False):
            article_name = art_path.stem
            logger.info(f"=== {article_name} : {model_name} ===")
            article_text = art_path.read_text(encoding="utf-8")

            # Run configured pre-passes
            outputs: Dict[str, Dict[str, Any]] = {}
            for p in tqdm(
                cfg.pre_passes,
                desc=f"{article_name} pre-passes",
                position=2,
                leave=False,
            ):
                try:
                    outputs[p.name] = run_single_pass(
                        model=model,
                        article_text=article_text,
                        pass_cfg=p,
                        out_base=out_base,
                        article_stem=article_name,
                        tools=tools,
                        logger=logger,
                        ollama_parameters=cfg.ollama_parameters,
                        model_name=model_name,
                    )
                except Exception:
                    logger.exception(
                        f"Pass failed: {p.name} : {article_name} : {model_name}"
                    )

            all_found_sequences = list(
                sorted(
                    set(
                        set(outputs.get("SeqPrompt_strict", [])).union(
                            outputs.get("SeqPrompt", [])
                        )
                    )
                )
            )
            all_found_sequences_str = ", ".join(all_found_sequences)
            logger.info("Pre-passes done, found sequences: " + all_found_sequences_str)

            sequence_descriptors = run_query_model(
                model=model,
                article_text=article_text,
                sequences=all_found_sequences,
                out_base=out_base,
                article_stem=article_name,
                common_prompt_path=cfg.common_prompt_path,
                ollama_parameters=cfg.ollama_parameters,
                logger=logger,
                model_name=model_name,
                tqdm_position=3,
            )

            stamp = _now_stamp()
            full_dir = out_base / "json_full"
            full_dir.mkdir(parents=True, exist_ok=True)
            full_seq_desc_path = (
                full_dir
                / f"{article_name}_{model_name_encode(model_name)}__SeqDesc-FULL__{stamp}.json"
            )
            full_seq_desc_path.write_text(
                json.dumps(sequence_descriptors, indent=2, ensure_ascii=False),
                encoding="utf-8",
            )

            for i, seq in enumerate(
                tqdm(
                    all_found_sequences,
                    desc=f"{article_name}: sequences construction",
                    leave=False,
                    position=3,
                )
            ):
                for construct_pass in tqdm(
                    cfg.construct_single_experiment_passes,
                    desc="Construction schemas",
                    leave=False,
                ):
                    try:
                        run_construct_single_experiment_pass(
                            model=model,
                            article_text=article_text,
                            sequence=seq,
                            sequence_id=i,
                            pass_cfg=construct_pass,
                            out_base=out_base,
                            article_stem=article_name,
                            tools=tools,
                            logger=logger,
                            ollama_parameters=cfg.ollama_parameters,
                            model_name=model_name,
                        )
                    except Exception:
                        logger.exception(
                            f"Pass failed: {p.name} : {article_name} : {model_name}"
                        )

            # for p in tqdm(cfg.passes, desc=f"{article_name} passes", leave=False):
            #     try:
            #         outputs[p.name] = run_single_pass(
            #             model=model,
            #             article_text=article_text,
            #             pass_cfg=p,
            #             out_base=out_base,
            #             article_stem=article_name,
            #             tools=tools,
            #             logger=logger,
            #             ollama_parameters=cfg.ollama_parameters,
            #             model_name=model_name,
            #         )
            #     except Exception:
            #         logger.exception(f"Pass failed: {p.name} : {article_name} : {model_name}")

            # # Stitch only if the expected pass names are present
            # try:
            #     A = outputs.get("A_core", {})
            #     B = outputs.get("B_index", {})
            #     # C = outputs.get("C_sequences", {})
            #     C = aggregate_c_outputs(outputs)
            #     D = outputs.get("D_parameters", {})
            #     E = outputs.get("E_outcomes", {})
            #     F = outputs.get("F_pairings", {})
            #     full_obj = stitch_full(A, B, C, D, E, F)

            #     # Final validation
            #     if full_validator:
            #         errs = sorted(full_validator.iter_errors(full_obj), key=lambda e: e.path)
            #         if errs:
            #             logger.error(f"[FULL] validation errors for {article_name} : {model_name}:\n" + "\n".join(str(e) for e in errs))
            #         else:
            #             logger.info(f"[FULL] validation OK for {article_name} : {model_name}")

            #     # Save full object (timestamped)
            #     stamp = _now_stamp()
            #     full_dir = out_base / "json_full"
            #     full_dir.mkdir(parents=True, exist_ok=True)
            #     full_path = full_dir / f"{article_name}_{model_name_encode(model_name)}__FULL__{stamp}.json"
            #     full_path.write_text(json.dumps(full_obj, indent=2, ensure_ascii=False), encoding="utf-8")
            #     logger.info(f"[FULL] wrote {full_path.name} {article_name} : {model_name}")

            #     # Optional DB insert
            #     if cfg.db_path:
            #         try:
            #             from hyb_db import insert_article_object  # your earlier module
            #             run_id = insert_article_object(
            #                 db_path=str(cfg.db_path),
            #                 article_obj=full_obj,
            #                 model_name=model_name,
            #                 article_name=article_name,
            #             )
            #             logger.info(f"[DB] inserted run_id={run_id} for {article_name} : {model_name}")
            #         except Exception:
            #             logger.exception("[DB] insertion failed")
            # except Exception:
            #     logger.exception(f"[FULL] stitching failed for {article_name} : {model_name}")


# Optional CLI hook (project_dir arg)
if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python pipeline_filedriven.py <project_dir>")
        sys.exit(1)
    run_project(sys.argv[1])
