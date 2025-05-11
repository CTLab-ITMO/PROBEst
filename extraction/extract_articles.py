#!/usr/bin/env python3
import os
import re
import json
import tempfile
import threading
import argparse
import logging
import queue
import warnings
import pdfplumber
import pytesseract
import requests
import torch
import sqlite3
from datetime import datetime
from PIL import Image
from tqdm import tqdm
from dicttoxml import dicttoxml
from lxml import etree
from transformers import (
    AutoTokenizer,
    AutoModelForCausalLM,
    TextIteratorStreamer,
    GenerationConfig,
)
from rich.logging import RichHandler

# ── ENVIRONMENT & WARNINGS ─────────────────────────────────────────────────────
os.environ["XFORMERS_ENABLE_SDW"] = "0"  # disable sliding-window attention
warnings.filterwarnings("ignore", "Sliding Window Attention")

# ── Configuration ──────────────────────────────────────────────────────────────
MODEL_NAME    = "Qwen/Qwen2.5-7B-Instruct-1M"
PROMPT_FILE   = "extraction_prompt.txt"
INPUT_DIR     = "articles"
RAW_OUT_DIR   = os.path.join("outputs", "raw")
JSON_OUT_DIR  = os.path.join("outputs", "json")
DTD_DIR       = "schema"
DTD_FILES     = {k: os.path.join(DTD_DIR, f) for k, f in {"hybridization_experiments": "experiment.dtd", "probe_samples": "probe.dtd"}.items()}
DB_PATH       = "outputs/extraction.db"

OLLAMA_URL    = "http://localhost:11434/api/generate"
OLLAMA_MODEL  = "myaniu/qwen2.5-1m:14b"

# ── Logging setup ───────────────────────────────────────────────────────────────
logging.basicConfig(
    level="INFO",
    format="%(asctime)s • %(levelname)s • %(message)s",
    handlers=[RichHandler(rich_tracebacks=True)]
)
log = logging.getLogger("extractor")

# ── Helpers ───────────────────────────────────────────────────────────────────
def ensure_dirs():
    for d in (RAW_OUT_DIR, JSON_OUT_DIR, os.path.dirname(DB_PATH)):
        os.makedirs(d, exist_ok=True)

# PDF text extraction with OCR fallback
def extract_text_from_pdf(path: str) -> str:
    text = ""
    with pdfplumber.open(path) as pdf:
        for page in pdf.pages:
            txt = page.extract_text() or ""
            if txt.strip(): text += txt + "\n"
            else:
                img = page.to_image(resolution=300).original
                text += pytesseract.image_to_string(img) + "\n"
    return text

def json_to_etree(js, parent=None):
    """
    Recursively convert JSON to an etree so that
    lists under key 'experiment' become repeated <experiment> tags, etc.
    """
    if parent is None:
        parent = etree.Element("root")
    for k, v in js.items():
        if isinstance(v, dict):
            node = etree.SubElement(parent, k)
            json_to_etree(v, node)
        elif isinstance(v, list):
            # produce one child per list-item, using the singular form of k
            singular = k.rstrip('s')  # crude singularization
            for item in v:
                node = etree.SubElement(parent, singular)
                if isinstance(item, (dict, list)):
                    json_to_etree(item, node)
                else:
                    node.text = str(item)
        else:
            node = etree.SubElement(parent, k)
            node.text = str(v)
    return parent

def item_func(parent_name):
    # when parent is 'experiment', name each list item 'experiment'
    if parent_name == 'experiment':
        return 'experiment'
    if parent_name == 'concentration':
        return 'concentration'
    if parent_name == 'parameter':
        return 'parameter'
    if parent_name == 'numeric_parameter':
        return 'numeric_parameter'
    # fallback for all others
    return parent_name.rstrip('s')


def validate_dtd(root: etree.Element) -> bool:
    dtd_path = DTD_FILES[root.tag]
    dtd = etree.DTD(open(dtd_path, 'rb'))
    if not dtd.validate(root):
        errs = "\n".join(str(e) for e in dtd.error_log.filter_from_errors())
        log.error(f"DTD validation failed ({dtd_path}):\n{errs}")
        return False
    return True


# JSON -> XML -> DTD validation
def validate_json_via_dtd(json_obj: dict) -> bool:
    if 'article_data' not in json_obj:
        log.error("Missing article_data in JSON object")
        return False

    if 'hybridization_experiments' not in json_obj['article_data']:
        log.error("Missing hybridization_experiments in the JSON object")
        return False
    
    root_he = json_to_etree(json_obj['article_data']['hybridization_experiments'])
    root_he.tag = 'hybridization_experiments'

    he = validate_dtd(root_he)

    if not he:
        return False
    
    if 'probe_samples' not in json_obj['article_data']:
        log.warning("No probe samples were found in the article?")
        return True

    root_ps = json_to_etree(json_obj['article_data']['probe_samples'])
    root_ps.tag = 'probe_samples'

    return validate_dtd(root_ps)
    


# Ollama API streaming
def call_ollama(prompt: str):
    headers = {"Content-Type": "application/json"}
    payload = {"model": OLLAMA_MODEL, "prompt": prompt, "max_tokens":2048, "stream":True}
    with requests.post(OLLAMA_URL, json=payload, headers=headers, stream=True) as resp:
        resp.raise_for_status()
        buffer = b""
        for chunk in resp.iter_content(chunk_size=None):
            if chunk is None: continue
            buffer += chunk
            # split on newline
            while b'\n' in buffer:
                line, buffer = buffer.split(b'\n', 1)
                if not line: continue
                data = json.loads(line.decode('utf-8'))
                yield data.get('response', '')
                if data.get('done', False): return

# ── SQLite DB setup ────────────────────────────────────────────────────────────
def init_db(conn):
    c = conn.cursor()
    # Articles
    c.execute('''
    CREATE TABLE IF NOT EXISTS articles (
      id INTEGER PRIMARY KEY,
      name TEXT UNIQUE,
      doi TEXT,
      processed_at TEXT
    )''')
    # Experiments
    c.execute('''
    CREATE TABLE IF NOT EXISTS experiments (
      id INTEGER PRIMARY KEY AUTOINCREMENT,
      article_id INTEGER,
      id_exp INTEGER,
      raw_description TEXT,
      type TEXT,
      organism TEXT,
      rna_impurities TEXT,
      annealing TEXT,
      ph REAL,
      FOREIGN KEY(article_id) REFERENCES articles(id)
    )''')
    # Concentrations
    c.execute('''
    CREATE TABLE IF NOT EXISTS concentrations (
      exp_id INTEGER,
      dna_rna_concentration TEXT,
      id_probe INTEGER,
      concentration_article TEXT,
      concentration_si REAL,
      FOREIGN KEY(exp_id) REFERENCES experiments(id)
    )''')
    # String params
    c.execute('''
    CREATE TABLE IF NOT EXISTS string_parameters (
      exp_id INTEGER,
      temperature TEXT,
      tris TEXT,
      na TEXT,
      k TEXT,
      mg TEXT,
      dmso TEXT,
      FOREIGN KEY(exp_id) REFERENCES experiments(id)
    )''')
    # Numeric params
    c.execute('''
    CREATE TABLE IF NOT EXISTS numeric_parameters (
      exp_id INTEGER,
      temperature REAL,
      tris REAL,
      na REAL,
      k REAL,
      mg REAL,
      dmso REAL,
      FOREIGN KEY(exp_id) REFERENCES experiments(id)
    )''')
    # Probe samples
    c.execute('''
    CREATE TABLE IF NOT EXISTS probe_groups (
      id INTEGER PRIMARY KEY AUTOINCREMENT,
      article_id INTEGER,
      id_group INTEGER,
      FOREIGN KEY(article_id) REFERENCES articles(id)
    )''')
    c.execute('''
    CREATE TABLE IF NOT EXISTS results (
      group_id INTEGER,
      id_exp INTEGER,
      outcome TEXT,
      fluorescence REAL,
      FOREIGN KEY(group_id) REFERENCES probe_groups(id)
    )''')
    c.execute('''
    CREATE TABLE IF NOT EXISTS related_sequences (
      group_id INTEGER,
      sequence TEXT,
      FOREIGN KEY(group_id) REFERENCES probe_groups(id)
    )''')
    c.execute('''
    CREATE TABLE IF NOT EXISTS probes (
      group_id INTEGER,
      id_probe INTEGER,
      probe_sequence TEXT,
      FOREIGN KEY(group_id) REFERENCES probe_groups(id)
    )''')
    c.execute('''
    CREATE TABLE IF NOT EXISTS modifications (
      group_id INTEGER,
      id_probe INTEGER,
      mod_pos INTEGER,
      mod_type TEXT,
      FOREIGN KEY(group_id) REFERENCES probe_groups(id)
    )''')
    conn.commit()

# Insert data into DB
def insert_into_db(conn, article_name, data, force=False):
    data = data['article_data']
    c = conn.cursor()

    # 1) ARTICLE: insert or skip
    c.execute("SELECT id FROM articles WHERE name = ?", (article_name,))
    existing = c.fetchone()
    if existing:
        if not force:
            log.info("Skipping existing article %s", article_name)
            return
        # delete cascade (assumes FK with ON DELETE CASCADE or manual cleanup)
        log.info("Reprocessing article %s (force)", article_name)
        c.execute("DELETE FROM articles WHERE id = ?", (existing[0],))
        conn.commit()

    doi = data['hybridization_experiments']['doi']
    c.execute(
        "INSERT INTO articles (name, doi, processed_at) VALUES (?, ?, ?)",
        (article_name, doi, datetime.utcnow().isoformat())
    )
    article_id = c.lastrowid

    # 2) HYBRIDIZATION EXPERIMENTS
    for exp in data['hybridization_experiments']['experiment']:
        c.execute(
            """INSERT INTO experiments
               (article_id, id_exp, raw_description, type, organism,
                rna_impurities, annealing, ph)
               VALUES (?, ?, ?, ?, ?, ?, ?, ?)""",
            (
                article_id,
                exp['id_exp'],
                exp['raw_description'],
                exp['type'],
                exp['organism'],
                str(exp['rna_impurities']),
                str(exp['annealing']),
                exp['ph']
            )
        )
        exp_row = c.lastrowid

        # concentrations (one row per concentration)
        for conc in exp['concentrations']['concentration']:
            c.execute(
                "INSERT INTO concentrations VALUES (?, ?, ?, ?, ?)",
                (
                    exp_row,
                    conc['dna_rna_concentration'],
                    conc['id_probe'],
                    conc['concentration_article'],
                    conc['concentration_si']
                )
            )

        # string_parameters
        for sp in exp['string_parameters']['parameter']:
            c.execute(
                "INSERT INTO string_parameters VALUES (?, ?, ?, ?, ?, ?, ?)",
                (
                    exp_row,
                    sp['temperature'],
                    sp['tris'],
                    sp['na'],
                    sp['k'],
                    sp['mg'],
                    sp['dmso']
                )
            )

        # numeric_parameters
        for np_ in exp['numeric_parameters']['numeric_parameter']:
            c.execute(
                "INSERT INTO numeric_parameters VALUES (?, ?, ?, ?, ?, ?, ?)",
                (
                    exp_row,
                    np_['temperature'],
                    np_['tris'],
                    np_['na'],
                    np_['k'],
                    np_['mg'],
                    np_['dmso']
                )
            )

    # 3) PROBE SAMPLES
    # There may be multiple sample_group entries
    for grp in data['probe_samples']['sample_group']:
        c.execute(
            "INSERT INTO probe_groups (article_id, id_group) VALUES (?, ?)",
            (article_id, grp['id_group'])
        )
        group_row = c.lastrowid

        # result rows
        for res in grp['results']['result']:
            c.execute(
                "INSERT INTO results VALUES (?, ?, ?, ?)",
                (
                    group_row,
                    res['id_exp'],
                    str(res['outcome']),
                    res['fluorescence']
                )
            )

        # related sequences at group level
        for seq in grp['sequences']['related_sequences']['related_sequence']:
            c.execute(
                "INSERT INTO related_sequences VALUES (?, ?)",
                (group_row, seq)
            )

        # probes (and their modifications + per-probe related_sequences)
        for pr in grp['sequences']['probe_sequences']['probes']['probe']:
            c.execute(
                "INSERT INTO probes VALUES (?, ?, ?)",
                (group_row, pr['id_probe'], pr['probe_sequence'])
            )
            # per-probe related sequences (if you also need to record these separately)
            for r in pr.get('related_sequences', {}):
                c.execute(
                    "INSERT INTO related_sequences VALUES (?, ?)",
                    (group_row, r)
                )
            # modifications
            for mod in pr.get('modification', []):
                c.execute(
                    "INSERT INTO modifications VALUES (?, ?, ?, ?)",
                    (
                        group_row,
                        pr['id_probe'],
                        mod['modification_pos'],
                        mod['modification_type']
                    )
                )

    conn.commit()

# Main workflow

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", choices=["local","ollama"], default="local")
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()

    log.info("Preparing directories and DB…")
    ensure_dirs()
    conn = sqlite3.connect(DB_PATH)
    init_db(conn)

    with open(PROMPT_FILE, 'r') as f:
        instruction = f.read().strip()

    if args.mode == 'local':
        log.info("Loading local model %s…", MODEL_NAME)
        device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
        tokenizer = AutoTokenizer.from_pretrained(MODEL_NAME, trust_remote_code=True)
        model = AutoModelForCausalLM.from_pretrained(MODEL_NAME, trust_remote_code=True).to(device)
        model.eval()
        tokenizer.pad_token = tokenizer.eos_token
        model.config.attn_config = None
        def local_generate(prompt):
            tokens = tokenizer(prompt, return_tensors="pt", padding=True, truncation=True).to(device)
            streamer = TextIteratorStreamer(tokenizer, skip_prompt=True, timeout=10.0)
            thread = threading.Thread(
                target=model.generate,
                kwargs={
                    'input_ids': tokens.input_ids,
                    'attention_mask': tokens.attention_mask,
                    'generation_config': GenerationConfig(max_new_tokens=2048),
                    'streamer': streamer,
                    'do_sample': False
                }, daemon=True)
            thread.start()
            q = streamer.text_queue
            while thread.is_alive() or not q.empty():
                try: yield q.get(timeout=1.0)
                except queue.Empty: continue
        infer_fn = local_generate
    else:
        log.info("Using Ollama at %s…", OLLAMA_URL)
        infer_fn = call_ollama

    all_articles = queue.Queue() 
    for f_article in sorted(set(os.listdir(INPUT_DIR))):
        all_articles.put(f_article)

    ok_articles: int = 0
    cached_articles: int = 0
    failed_articles: int = 0
    
    while not all_articles.empty():
        fname = all_articles.get()
        if not fname.lower().endswith('.pdf'): continue
        base = os.path.splitext(fname)[0]
        raw_path = os.path.join(RAW_OUT_DIR, f"{base}.txt")
        json_path = os.path.join(JSON_OUT_DIR, f"{base}.json")

        log.info("Parsing article %s", base)

        # If JSON exists: validate and DB insert (or skip)
        if os.path.isfile(json_path):
            with open(json_path,'r') as f: parsed = json.load(f)
            # validate
            if validate_json_via_dtd(parsed):
                insert_into_db(conn, base, parsed, force=args.force)
                cached_articles += 1
                log.info("Valid JSON cache was found for article %s", base)
                continue
            else:
                log.info("JSON cache was found for article %s, but failed validation", base)
        # If raw exists but no JSON: parse & validate
        if os.path.isfile(raw_path):
            raw_str = open(raw_path).read().strip()
            parsed = json.loads(raw_str)
            # save JSON
            if validate_json_via_dtd(parsed):
                with open(json_path,'w', encoding='utf-8') as f: json.dump(parsed,f,indent=2)
                insert_into_db(conn, base, parsed, force=args.force)
                cached_articles += 1
                log.info("Valid RAW cache was found for article %s", base)
                continue
            else:
                log.info("RAW cache was found for article %s, but failed validation", base)

        # ELSE: run full extraction
        log.info("▶ Processing %s…", fname)
        full_text = extract_text_from_pdf(os.path.join(INPUT_DIR,fname))
        prompt = f"{instruction}\n\n{full_text}"
        # stream LLM
        out_accum = ''
        for chunk in infer_fn(prompt):
            print(chunk,end='',flush=True); out_accum+=chunk
            if out_accum.count('```')>=2: break
        print()
        m = re.search(r"```json\n(.*?)\n```", out_accum, re.S)
        if not m: log.error("No JSON delimiters for %s", fname); continue
        raw_json = m.group(1).strip()
        with open(raw_path,'w') as f: f.write(raw_json)
        parsed = json.loads(raw_json)
        # validate + insert
        if validate_json_via_dtd(parsed):
            with open(json_path,'w',encoding='utf-8') as f: json.dump(parsed,f,indent=2)
            insert_into_db(conn, base, parsed, force=args.force)
            ok_articles += 1
            log.info("Article %s was parsed successfully", base)
        else:
            failed_articles += 1
            all_articles.put(fname)
            log.info("Article %s was not successfully parsed and re-enqueued", base)

    conn.close()
    log.info("%d new articles were successfully parsed, %d articles had valid cache, %d article extraction attempts failed", ok_articles, cached_articles, failed_articles)
    log.info("✅ All done.")

if __name__=='__main__': 
    main()
