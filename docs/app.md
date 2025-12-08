## Web App (`app/app.py`)

Flask UI for wrapping `pipeline.py`. Users upload input FASTA plus BLAST databases (true/off-target), the app builds databases, runs the pipeline, and returns top probes and a downloadable results zip.

### Endpoints
- `GET /`: Renders `templates/index.html`.
- `POST /process`: Main entry. Validates uploads, builds BLAST DBs, runs pipeline, returns JSON.
- `GET /download`: Zips the session output directory and streams it back.

### Upload expectations (via `/process`)
- Input FASTA: single file or archive (`.zip`, `.tar.gz`, `.tgz`, `.tar`) containing FASTA. First found FASTA is used.
- True database: required archive of FASTA(s); converted to BLAST DB via `scripts/generator/prep_db.sh`.
- False databases: ≥1 archives of FASTA(s); each becomes its own BLAST DB.
- All uploads land in `app/uploads/` under a UUID session prefix; results in `app/results/<session_id>/`.

### Pipeline invocation
- Runs `python ../pipeline.py` with:
  - Required: `-i <input.fa> -tb <true_db> -fb <false_db...> -c <contig_table> -o <session_result_dir>`
  - Optional form fields mapped to args: `threads (-t)`, `algorithm (-a)`, `iterations (-N)`, `top (-T)`, `mutation_rate (-M)`, `indel_rate (-I)`, `set_size (-S)`, `append (-A)`, Primer3 knobs (`PRIMER_*`), BLAST knobs (`word_size`, `reward`, `penalty`, `gapopen`, `gapextend`, `evalue`), probe_check knobs (`max_mismatch`, `multimap_max`, `negative_max`, `min_ident`).
- `prep_db.sh` is made executable and run per archive to generate BLAST DBs plus contig tables under the session.

### Result parsing & response
- Reads `<session>/output.fa` and extracts probe headers (`>H{hits}_{name}`) to list probes and top 5 by hit count.
- `/process` JSON response: `{success, top_probes, total_probes}` or `{error}`.
- `/download` returns a zip of the entire session output directory (FASTA, TSVs, logs).

### Runtime limits and configs
- Max upload size: 500 MB (`MAX_CONTENT_LENGTH`).
- Upload/result roots configurable in `app.config['UPLOAD_FOLDER']` and `RESULTS_FOLDER`; created on startup.
- Secret key generated per run (session cookies only used to track session/output paths).

### Running the app
```bash
cd app
export FLASK_APP=app.py
python app.py  # dev server on 0.0.0.0:5000, debug=True
```

### Notes
- Archives must contain FASTA; otherwise `/process` returns 400.
- If `output.fa` is missing after pipeline run, `/process` returns 500 with an error message.
- Currently picks the first FASTA in an input archive; extend `find_fasta_files`/merge logic if you need multi-FASTA support.


