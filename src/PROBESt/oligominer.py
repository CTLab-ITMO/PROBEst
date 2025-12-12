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
    1. Runs blockParse.py on the entire input FASTA file to generate a FASTQ file with candidate probes.
    2. Converts the FASTQ file to FASTA format.
    3. Optionally appends additional sequences to the FASTA file.

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
    
    # Step 1: Run blockParse.py on the entire input FASTA file
    # blockParse.py generates a FASTQ file with candidate probes
    # It writes output to current working directory using filename stem (split('.')[0])
    # blockParse uses split('.')[0] on the input file path to get the stem
    # For "/path/to/input.fa", split('.')[0] gives "/path/to/input"
    # Then it writes "{outName}.fastq" where outName is that path
    # Since it writes to current working directory, it creates ".fastq" in the main output directory
    blockparse_log = os.path.join(oligominer_tmp, 'blockparse.log')
    
    # The FASTQ file is created in the main output directory (args.output), not in out_dir
    # Check for .fastq in the main output directory
    main_output_dir = args.output
    expected_fastq = os.path.join(main_output_dir, '.fastq')
    
    # Build blockParse command - only add -l/-L and -T if provided
    # Use absolute path for input file to avoid path issues
    input_fasta_abs = os.path.abspath(input_fasta)
    blockparse_cmd = [
        python_cmd,
        os.path.join(oligominer_path, 'blockParse.py'),
        '-f', input_fasta_abs  # Use absolute path to input file
    ]
    
    # Add probe length arguments only if provided
    if args.oligominer_probe_length is not None:
        blockparse_cmd.extend(['-l', str(args.oligominer_probe_length)])
        blockparse_cmd.extend(['-L', str(args.oligominer_probe_length)])  # Set max length same as min for fixed length
    
    # Add temperature argument only if provided
    if args.oligominer_temperature is not None:
        blockparse_cmd.extend(['-T', str(args.oligominer_temperature)])
    
    # Redirect output to log file to suppress verbose messages
    # Run from main output directory so FASTQ is created there
    with open(blockparse_log, 'w') as log_file:
        result = subprocess.run(
            blockparse_cmd,
            cwd=main_output_dir,  # Run from main output directory so FASTQ is created there
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
                    f"OligoMiner blockParse.py failed with Python version error. "
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
    
    # Check if FASTQ was generated
    # blockParse writes to current working directory (which is now main_output_dir)
    # It creates .fastq in the main output directory
    fastq_file = expected_fastq
    
    # Final check
    if not os.path.exists(fastq_file) or os.path.getsize(fastq_file) == 0:
        if not probe_found or probe_count == 0:
            raise RuntimeError(
                "OligoMiner blockParse.py found no candidate probes. "
                "Try adjusting probe generation parameters or using a different initial generator."
            )
        else:
            raise RuntimeError(
                f"OligoMiner blockParse.py found {probe_count} probes but FASTQ file was not found. "
                f"Expected at: {expected_fastq}. "
                "Check the log file for details."
            )
    
    print("OligoMiner done")
    
    # Step 2: Convert FASTQ file to FASTA
    output_fasta_path = os.path.join(out_dir, "output.fa")
    with open(output_fasta_path, 'w') as out_fasta:
        try:
            for record in SeqIO.parse(fastq_file, "fastq"):
                # Write as FASTA: header and sequence
                out_fasta.write(f">{record.id}\n{record.seq}\n")
        except Exception as e:
            raise RuntimeError(
                f"Failed to convert FASTQ to FASTA: {e}. "
                f"FASTQ file: {fastq_file}"
            )
    
    # Append additional sequences if provided
    if hasattr(args, 'add_set') and args.add_set:
        add_fasta_cmd = f"cat {' '.join(args.add_set)} >> {output_fasta_path}"
        subprocess.run(add_fasta_cmd, shell=True, executable="/bin/bash")
