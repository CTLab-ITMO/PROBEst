## AI Filtration Module

The AI module scores BLASTn positive hits to drop likely off-target probes before multimapping filtration. It runs automatically inside `pipeline.py` when `--AI` is enabled (default).

### Runtime usage in `pipeline.py`
- Source: `PROBESt/AI.py`
- Entry: `apply_ai_filtration_to_blast_file(...)`
- Inputs: `positive_hits.tsv`, trained model (`data/AIL.pt` by default), the iteration `merged.fa`
- Flow:
  1) Read and extend BLAST hits with numeric features
  2) Load TorchClassifier (architecture auto-detected if not passed)
  3) Predict per-hit `ai_score` and per-probe `sum_score`
  4) Keep probes with `sum_score <= threshold` (median by default)
  5) Rewrite `positive_hits.tsv` for downstream probe_check
- CLI flags:
  - `--AI / --no-AI` toggle
  - `--ai_threshold_method {median,mean}` (if exposed via args; default median)

### Feature set
Numeric columns expected (missing ones are zero-filled): `Formamide`, `GCcontent`, `Lengthnt`, `evalue`, `mismatches`, `length`, `bitscore`, `identity`, `score`, `hairpin_prob`, `dimer_DNA`, `dimer_DNA_flank`, `dimer_probe`, `dimer_probe_DNA`.

### Model loading
- Wrapper: `load_torch_classifier(model_path, input_size, model_architecture=None)`
- Auto-detects architecture from checkpoint if absent.
- Supported architectures (from `models_registry`): `GAILDiscriminator`, `GAILWide`, `GAILDeep`, `GAILNarrow`, `GAILWithDropout`, `GAILWideDeep`, `GAILWideDropout`, `GAILWideBatchNorm`, `GAILWideExtra`, `GAILWideBalanced`.

### Training script (`scripts/generator/ML_filtration.py`)
Purpose: train and compare classifiers for BLAST-hit filtration and export the best Torch model.

Key steps:
- Load training CSV (`data/databases/open/test_ML_database.csv`), option to tokenize sequences with `tokenize_table`.
- Split into train/val/test; drop non-numeric columns except `type`.
- Train candidates:
  - Torch models via `TorchClassifier`: `ShallowNet`, `WideNet`, `ResidualNet`, `GAILDiscriminator`, `TabTransformer`.
  - Baselines: `DeepNeuralNetworkModel`, `LogisticRegressionModel`.
- Validate on held-out data; plot ROC/metrics and learning curves.
- Architecture search focused on GAIL variants (`GAILDeep`, `GAILWide`, `GAILNarrow`, `GAILWideExtra`, etc.).
- Select best model by F1; save checkpoint (weights, scaler, learning rate, pos_weight, architecture name) to `data/GAIL.pt` (can be renamed to `AIL.pt` for pipeline use).
- Exports predictions on test split and comparison plots (`tests_outs/`).

### Using a trained model in the pipeline
1) Place the chosen checkpoint at `data/AIL.pt` (or set `model_path` when calling `apply_ai_filtration_to_blast_file`).
2) Ensure `input_size` matches training (default 14 numeric features).
3) Optionally set `model_architecture`; otherwise it is inferred.

### Quick start: retrain AI
```bash
python scripts/generator/ML_filtration.py
# Outputs: data/GAIL.pt, tests_outs/* plots, test_predictions.csv
```

### Notes and limits
- Skips AI stage if `positive_hits.tsv` or model is missing/empty.
- Thresholding is median-based by default; mean is supported.
- Uses class weighting (`weight_pos`) to handle imbalance; baked into saved checkpoint.
