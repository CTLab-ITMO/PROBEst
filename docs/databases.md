## Databases Overview

PROBESt ships two main data sources for probe curation: (1) ProbeBase-derived tables parsed and optionally re-aligned to genomes, and (2) an article-extracted database built from PDF/PMC papers via the extraction pipeline.

### ProbeBase-derived datasets
- Location: `data/databases/open/`
  - `probeBase.csv`: raw scraped ProbeBase table.
  - `probeBase_false.csv`: negative/decoy entries for specificity checks.
  - `test_ML_database.csv`: numeric feature table for AI filtration training/validation.
- Parsing scripts (under `scripts/databases/`):
  - `probeBase.py`: fetch + parse a single ProbeBase page to a uniform row.
  - `probeBase_parse.py`: end-to-end processor that pivots the raw CSV, extracts species/genus from taxonomy, optionally downloads genomes (`data/genomes`), runs BLAST against each genus, and writes enriched rows with flanking sequences.
    - Key steps:
      1) Read raw CSV → pivot to wide format (`id`, `Sequence`, `Taxonomy`).
      2) Extract species and genera (`extract_species_from_taxonomy` / `clean_taxon`).
      3) Optionally download reference genomes via `ncbi-genome-download`; reuse cached genomes in `data/genomes`.
      4) For each probe/genera: create BLAST DB, search probe, capture identity/length/mismatches/evalue/bitscore and ±15bp flanks.
      5) Save enriched CSV (default: `data/databases/open/probeBase_genome_results.csv` when run as `__main__`).
- Supporting scripts:
  - `probeBase_wide.py`, `probeBase_analyse.R`, `probeBase_wide.R`, `generate_noisy_probes.py`, `LLM_benchmarkig.R` (analysis/augmentation helpers).

### Article-extracted database
- Pipeline location: `extraction/`
  - Driver: `extract_articles.py`
  - Prompt: `extraction/extraction_prompt.txt`
  - Schemas/DTDs: `extraction/schema/experiment.dtd`, `extraction/schema/probe.dtd`
  - Outputs: `extraction/outputs/`
    - `extraction.db` / `_extraction.db`: SQLite with structured tables (articles, experiments, concentrations, probes, primers, targets).
    - `raw/*.txt`: raw text extracted from PDFs (OCR fallback via `pytesseract`).
    - `json/*.json`: structured LLM-extracted records per article.
- High-level flow (`extract_articles.py`):
  1) Read PDF/PMC content (default `INPUT_DIR=/mnt/Models/articles2/`), extract text with `pdfplumber` + OCR fallback.
  2) Prompt LLM (Qwen/Ollama defaults) using `extraction_prompt.txt` to produce structured JSON.
  3) Validate JSON against DTDs by converting to XML (`experiment.dtd`, `probe.dtd`).
  4) Persist results into SQLite (`outputs/extraction.db`) and save raw/text/JSON artifacts.
- Related data hints:
  - `data/articles/pmc_result.txt`, `pmid-nucleotide-set.txt`: search result seed lists.
  - `data/articles/pdf_checked/`: manually verified PDFs (if present).

### Using these datasets in the pipeline
- ProbeBase-derived tables can seed true/false BLAST databases or provide labeled examples for AI filtration (e.g., `test_ML_database.csv`).
- The article-extracted SQLite/JSON exports can be transformed into FASTA/TSV inputs for downstream probe design; schemas ensure each probe is tied to experimental context.

### Regeneration quick start
- Parse ProbeBase to genomes:
  ```bash
  python scripts/databases/probeBase_parse.py \
    data/databases/open/probeBase.csv \
    data/databases/open/probeBase_genome_results.csv
  ```
- Run article extraction (GPU recommended):
  ```bash
  cd extraction
  python extract_articles.py --help  # see options
  ```


