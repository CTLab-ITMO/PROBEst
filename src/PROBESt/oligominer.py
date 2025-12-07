from Bio import SeqIO
import os
import subprocess
import re
import glob


def fastq_to_fasta(fastq_file: str, fasta_file: str) -> None:
    """
    Converts a FASTQ file to FASTA format.
    
    Args:
        fastq_file (str): Path to input FASTQ file.
        fasta_file (str): Path to output FASTA file.
    """
    with open(fasta_file, 'w') as out_fasta:
        for record in SeqIO.parse(fastq_file, "fastq"):
            # Write as FASTA: header and sequence
            out_fasta.write(f">{record.id}\n{record.seq}\n")


def initial_set_generation(args, out_dir: str) -> None:
    """
    Generates an initial set of probes using OligoMiner and writes them to a FASTA file.

    This function performs the following steps:
    1. Processes each sequence in the input FASTA file separately.
    2. Runs blockParse.py to generate FASTQ files with candidate probes.
    3. Converts FASTQ files to FASTA format and merges them.
    4. Optionally appends additional sequences to the FASTA file.

    Args:
        args: Parsed command-line arguments.
        out_dir (str): Path to the output directory.
    """
    input_fasta = os.path.join(out_dir, "input.fa")
    
    # Create temporary directory for OligoMiner processing
    oligominer_tmp = os.path.join(out_dir, "oligominer_tmp")
    os.makedirs(oligominer_tmp, exist_ok=True)
    
    # Get OligoMiner path
    oligominer_path = args.oligominer_path
    if not oligominer_path:
        raise ValueError("OligoMiner path not specified. Please set --oligominer_path.")
    
    # Determine Python interpreter for OligoMiner (Python 2 required)
    python_cmd = args.oligominer_python
    
    # Check if OLIGOMINER_PYTHON environment variable is set (from conda activation script)
    if not python_cmd:
        python_cmd = os.environ.get('OLIGOMINER_PYTHON', None)
    
    if not python_cmd:
        # Try to auto-detect Python 2
        python2_found = False
        for py_cmd in ['python2.7', 'python2', 'python']:
            result = subprocess.run(f"which {py_cmd}", shell=True, capture_output=True)
            if result.returncode == 0:
                # Check if it's actually Python 2
                version_check = subprocess.run(
                    f"{py_cmd} -c 'import sys; exit(0 if sys.version_info[0] == 2 else 1)'",
                    shell=True,
                    capture_output=True
                )
                if version_check.returncode == 0:
                    python_cmd = py_cmd
                    python2_found = True
                    break
        
        if not python2_found:
            raise RuntimeError(
                "Python 2 not found for OligoMiner. OligoMiner requires Python 2. "
                "Please install Python 2 (e.g., 'python2' or 'python2.7') or specify it using "
                "--oligominer_python argument. The installation script creates a separate "
                "conda environment 'probest_oligominer' with Python 2.7 and Biopython."
            )
    
    # Verify the specified Python is actually Python 2
    version_check = subprocess.run(
        f"{python_cmd} -c 'import sys; exit(0 if sys.version_info[0] == 2 else 1)'",
        shell=True,
        capture_output=True
    )
    if version_check.returncode != 0:
        raise RuntimeError(
            f"Specified Python interpreter '{python_cmd}' is not Python 2. "
            "OligoMiner requires Python 2. Please specify a Python 2 interpreter using "
            "--oligominer_python (e.g., 'python2' or 'python2.7')."
        )
    
    # Check if Biopython is installed in Python 2 environment
    biopython_check = subprocess.run(
        f"{python_cmd} -c 'from Bio.SeqUtils import MeltingTemp'",
        shell=True,
        capture_output=True
    )
    
    if biopython_check.returncode != 0:
        error_msg = biopython_check.stderr.decode('utf-8', errors='ignore') if biopython_check.stderr else "Unknown error"
        install_cmd = f"{python_cmd} -m pip install biopython==1.76"
        
        raise RuntimeError(
            f"Biopython is not installed in the Python 2 environment ({python_cmd}). "
            "OligoMiner requires Biopython 1.76 for Python 2. Please install it using:\n"
            f"  {install_cmd}\n"
            f"Error: {error_msg[:200]}"
        )
    
    print(f"Using Python interpreter for OligoMiner: {python_cmd}")
    
    # Process each sequence in the input FASTA file
    sequences = list(SeqIO.parse(input_fasta, "fasta"))
    all_fastq_files = []
    
    # Process each sequence separately
    for record in sequences:
        # Create a temporary FASTA file for this sequence
        seq_fasta = os.path.join(oligominer_tmp, f"{record.id}.fasta")
        with open(seq_fasta, "w") as f:
            SeqIO.write(record, f, "fasta")
        
        # Extract prefix for naming
        prefix = re.sub(r'[^a-zA-Z0-9_]', '_', record.id)
        prefix = prefix[:50]  # Limit length
        
        # Step 1: Run blockParse.py
        # blockParse.py generates a FASTQ file with candidate probes
        # Note: blockParse.py writes output to current working directory using filename stem
        # It splits on '.' and takes first part, so we need to handle this carefully
        blockparse_log = os.path.join(oligominer_tmp, f'{prefix}_blockparse.log')
        
        # Get the base name of the input file (without path) for output filename
        seq_fasta_basename = os.path.basename(seq_fasta)
        # blockParse uses split('.')[0] to get the stem, so we need to match that
        # For "AE006914.1_cds_AAL02539.1_1.fasta", split('.')[0] = "AE006914"
        seq_fasta_stem = seq_fasta_basename.split('.')[0]
        # blockParse will create {stem}.fastq in the current working directory
        # Since we run with cwd=oligominer_tmp, the output will be there
        expected_fastq = os.path.join(oligominer_tmp, f'{seq_fasta_stem}.fastq')
        
        # Run blockParse.py in the tmp directory so output goes there
        # Redirect both stdout and stderr to log file to suppress verbose output
        # Use relative path for input file since we're changing to that directory
        blockparse_cmd = [
            python_cmd,
            os.path.join(oligominer_path, 'blockParse.py'),
            '-l', str(args.oligominer_probe_length),
            '-L', str(args.oligominer_probe_length),  # Set max length same as min for fixed length
            '-T', str(args.oligominer_temperature),
            '-f', seq_fasta_basename  # Use basename since we're running in oligominer_tmp
        ]
        # Redirect output to log file to suppress verbose messages
        with open(blockparse_log, 'w') as log_file:
            result = subprocess.run(
                blockparse_cmd,
                cwd=oligominer_tmp,  # Change to tmp directory so output goes there
                stdout=log_file,
                stderr=subprocess.STDOUT,  # Redirect stderr to stdout
                text=True
            )
        
        # Check for errors in blockParse and whether probes were found
        probe_found = False
        probe_count = 0
        
        if os.path.exists(blockparse_log):
            with open(blockparse_log, 'r') as f:
                log_content = f.read()
                
                # Check for Python version errors
                if 'SyntaxError' in log_content or 'print' in log_content.lower():
                    raise RuntimeError(
                        f"OligoMiner blockParse.py failed with Python version error for {record.id}. "
                        f"OligoMiner requires Python 2. Please set --oligominer_python to a Python 2 interpreter "
                        f"(e.g., 'python2' or 'python2.7'). Error: {log_content[:500]}"
                    )
                
                # Check if probes were identified
                if 'candidate probes identified' in log_content:
                    # Extract number of probes
                    match = re.search(r'(\d+) candidate probes identified', log_content)
                    if match:
                        probe_count = int(match.group(1))
                        if probe_count > 0:
                            probe_found = True
                elif 'No candidate probes discovered' in log_content:
                    probe_found = False
                    probe_count = 0
        
        # If blockParse failed or no probes found, skip this sequence
        if result.returncode != 0 or not probe_found or probe_count == 0:
            continue  # Skip silently - no probes found for this sequence
        
        # Check if FASTQ was generated
        fastq_file = expected_fastq
        if not os.path.exists(fastq_file):
            # Probes were found but FASTQ not created - try to find it with different patterns
            # Search for any .fastq files in the tmp directory that might match
            possible_fastq = glob.glob(os.path.join(oligominer_tmp, '*.fastq'))
            if possible_fastq:
                # Use the most recently created one (likely from this run)
                fastq_file = max(possible_fastq, key=os.path.getmtime)
            else:
                # Still not found - skip this sequence
                continue
        
        # Add to list of FASTQ files to convert
        if os.path.exists(fastq_file) and os.path.getsize(fastq_file) > 0:
            all_fastq_files.append(fastq_file)
    
    print("OligoMiner done")
    
    # Step 2: Convert all FASTQ files to FASTA and merge them
    output_fasta_path = os.path.join(out_dir, "output.fa")
    with open(output_fasta_path, 'w') as out_fasta:
        for fastq_file in all_fastq_files:
            try:
                for record in SeqIO.parse(fastq_file, "fastq"):
                    # Write as FASTA: header and sequence
                    out_fasta.write(f">{record.id}\n{record.seq}\n")
            except Exception as e:
                # Skip files that can't be parsed
                continue
    
    # Append additional sequences if provided
    if hasattr(args, 'add_set') and args.add_set:
        add_fasta_cmd = f"cat {' '.join(args.add_set)} >> {output_fasta_path}"
        subprocess.run(add_fasta_cmd, shell=True, executable="/bin/bash")
