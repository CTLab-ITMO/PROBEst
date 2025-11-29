# PROBESt Web Application

A web-based interface for the PROBESt probe generation and optimization tool.

## Features

- **User-friendly web interface** for all PROBESt parameters
- **File upload** for input FASTA, BLAST databases, and contig tables
- **Real-time processing** with progress indicators
- **Results visualization** showing top 5 best probes
- **Download results** as a ZIP file containing all output files

## Installation

1. Install Flask and required dependencies:
```bash
pip install flask werkzeug
```

2. Ensure all PROBESt dependencies are installed (see main README.md)

3. Make sure `primer3_core` and `blastn` are available in your PATH

## Running the Application

From the project root directory:

```bash
cd app
python app.py
```

The application will be available at `http://localhost:5000`

## Usage

### Required Inputs

1. **Input FASTA File**: Input FASTA file for probe generation (.fna recommended)
   - Can be a single FASTA file (.fa, .fasta, .fna) or an archive (.zip, .tar.gz) containing FASTA files
   
2. **True Base Archive**: Archive containing FASTA files for primer adjustment
   - Must be a .zip or .tar.gz archive
   - Archive should contain one or more FASTA files (.fa, .fasta, .fna, optionally gzipped)
   - The archive will be automatically preprocessed using `prep_db.sh` to create a BLAST database
   
3. **False Base Archive(s)**: Archive(s) containing FASTA files for non-specific testing
   - Must be .zip or .tar.gz archives
   - Can upload multiple archives
   - Each archive will be automatically preprocessed using `prep_db.sh` to create BLAST databases
   
**Note**: The contig table is automatically generated from the true base archive during preprocessing.

### Optional Parameters

The web interface provides access to all optional parameters:

- **Basic Parameters**: Threads, Algorithm type
- **Evolutionary Algorithm**: Iterations, mutation rates, set sizes
- **Primer3 Parameters**: Primer size ranges, product size
- **BLAST Parameters**: Word size, scoring parameters
- **Probe Check Parameters**: Mismatch thresholds, identity requirements

### Results

After processing:

1. **Top 5 Probes** are displayed on screen with:
   - Probe name
   - Hit count (number of matches)
   - Sequence length
   - Full sequence

2. **Download Results** button provides a ZIP file containing:
   - `output.fa`: Final probe sequences
   - `stats.csv`: Iteration statistics
   - All intermediate files from the pipeline

## File Structure

```
app/
├── app.py              # Flask application
├── templates/
│   └── index.html     # Main web interface
├── static/
│   ├── style.css      # Styling
│   └── script.js      # Client-side JavaScript
├── uploads/           # Temporary file uploads (created automatically)
└── results/            # Processing results (created automatically)
```

## Notes

- Uploaded files are stored temporarily and cleaned up after processing
- Processing may take several minutes depending on input size and parameters
- The application uses session-based file management for security


