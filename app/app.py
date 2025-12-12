# MIT License
#
# Copyright (c) 2025 CTLab-ITMO
#
# Authors: Daniil Smutin, Aleksandr Serdiukov, Vitalii Dravgelis, Artem Ivanov,
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


import os
import sys
import re
import zipfile
import tarfile
import subprocess
import glob
from flask import Flask, render_template, request, jsonify, send_file, session
from werkzeug.utils import secure_filename
import uuid

# Add parent directory to path to import PROBESt modules
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

app = Flask(__name__)
app.secret_key = os.urandom(24)
app.config['MAX_CONTENT_LENGTH'] = 500 * 1024 * 1024  # 500MB max file size
app.config['UPLOAD_FOLDER'] = os.path.join(os.path.dirname(__file__), 'uploads')
app.config['RESULTS_FOLDER'] = os.path.join(os.path.dirname(__file__), 'results')

# Create necessary directories
os.makedirs(app.config['UPLOAD_FOLDER'], exist_ok=True)
os.makedirs(app.config['RESULTS_FOLDER'], exist_ok=True)

def parse_fasta_output(fasta_path):
    """Parse the output.fa file and extract probe information."""
    probes = []
    try:
        with open(fasta_path, 'r') as f:
            current_name = None
            current_seq = None
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    if current_name and current_seq:
                        # Extract hit count from header (format: >H{hit_count}_{name})
                        match = re.match(r'>H(\d+)_(.+)', current_name)
                        if match:
                            hit_count = int(match.group(1))
                            probe_name = match.group(2)
                            probes.append({
                                'name': probe_name,
                                'sequence': current_seq,
                                'hits': hit_count
                            })
                    current_name = line[1:]  # Remove '>'
                    current_seq = ''
                else:
                    if current_seq is not None:
                        current_seq += line
            # Handle last probe
            if current_name and current_seq:
                match = re.match(r'>H(\d+)_(.+)', current_name)
                if match:
                    hit_count = int(match.group(1))
                    probe_name = match.group(2)
                    probes.append({
                        'name': probe_name,
                        'sequence': current_seq,
                        'hits': hit_count
                    })
    except Exception as e:
        print(f"Error parsing FASTA: {e}")
    return probes

def get_top_probes(probes, n=5):
    """Get top N probes sorted by hit count."""
    sorted_probes = sorted(probes, key=lambda x: x['hits'], reverse=True)
    return sorted_probes[:n]

def extract_archive(archive_path, extract_dir):
    """Extract .zip or .tar.gz archive to directory."""
    os.makedirs(extract_dir, exist_ok=True)
    
    if archive_path.endswith('.zip'):
        with zipfile.ZipFile(archive_path, 'r') as zip_ref:
            zip_ref.extractall(extract_dir)
    elif archive_path.endswith('.tar.gz') or archive_path.endswith('.tgz'):
        with tarfile.open(archive_path, 'r:gz') as tar_ref:
            tar_ref.extractall(extract_dir)
    elif archive_path.endswith('.tar'):
        with tarfile.open(archive_path, 'r') as tar_ref:
            tar_ref.extractall(extract_dir)
    else:
        raise ValueError(f"Unsupported archive format: {archive_path}")

def find_fasta_files(directory):
    """Find all FASTA files in a directory recursively."""
    fasta_files = []
    patterns = ['*.fa', '*.fasta', '*.fna', '*.fa.gz', '*.fasta.gz', '*.fna.gz']
    
    for pattern in patterns:
        # Search recursively
        fasta_files.extend(glob.glob(os.path.join(directory, '**', pattern), recursive=True))
        # Search in current directory
        fasta_files.extend(glob.glob(os.path.join(directory, pattern), recursive=False))
    
    # Remove duplicates and sort
    fasta_files = sorted(list(set(fasta_files)))
    return fasta_files

def run_prep_db(fasta_files, database_name, contig_table_path, tmp_dir, script_path):
    """Run prep_db.sh to create BLAST database from FASTA files."""
    if not fasta_files:
        raise ValueError("No FASTA files found in archive")
    
    # Get absolute paths
    fasta_files = [os.path.abspath(f) for f in fasta_files]
    database_name = os.path.abspath(database_name)
    contig_table_path = os.path.abspath(contig_table_path)
    tmp_dir = os.path.abspath(tmp_dir)
    script_path = os.path.abspath(script_path)
    
    # Build command
    cmd = ['bash', script_path, 
           '-n', database_name,
           '-c', contig_table_path,
           '-t', tmp_dir] + fasta_files
    
    result = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
        cwd=os.path.dirname(script_path)
    )
    
    if result.returncode != 0:
        raise RuntimeError(f"prep_db.sh failed: {result.stderr}")
    
    return database_name

def run_pipeline(args_dict, session_id):
    """Run the pipeline with given arguments."""
    try:
        # Create temporary directory for this session
        temp_dir = os.path.join(app.config['RESULTS_FOLDER'], session_id)
        os.makedirs(temp_dir, exist_ok=True)
        
        # Build command
        cmd = ['python', os.path.join(os.path.dirname(__file__), '..', 'pipeline.py')]
        
        # Add required arguments
        cmd.extend(['-i', args_dict['input']])
        cmd.extend(['-tb', args_dict['true_base']])
        cmd.extend(['-fb'] + args_dict['false_base'])
        cmd.extend(['-c', args_dict['contig_table']])
        cmd.extend(['-o', temp_dir])
        
        # Add optional arguments
        if args_dict.get('threads'):
            cmd.extend(['-t', str(args_dict['threads'])])
        if args_dict.get('algorithm'):
            cmd.extend(['-a', args_dict['algorithm']])
        if args_dict.get('iterations'):
            cmd.extend(['-N', str(args_dict['iterations'])])
        if args_dict.get('top'):
            cmd.extend(['-T', str(args_dict['top'])])
        if args_dict.get('mutation_rate'):
            cmd.extend(['-M', str(args_dict['mutation_rate'])])
        if args_dict.get('indel_rate'):
            cmd.extend(['-I', str(args_dict['indel_rate'])])
        if args_dict.get('set_size'):
            cmd.extend(['-S', str(args_dict['set_size'])])
        if args_dict.get('append') is not None:
            cmd.extend(['-A', str(args_dict['append'])])
        
        # Primer3 arguments
        if args_dict.get('PRIMER_PICK_PRIMER'):
            cmd.extend(['--PRIMER_PICK_PRIMER', str(args_dict['PRIMER_PICK_PRIMER'])])
        if args_dict.get('PRIMER_NUM_RETURN'):
            cmd.extend(['--PRIMER_NUM_RETURN', str(args_dict['PRIMER_NUM_RETURN'])])
        if args_dict.get('PRIMER_OPT_SIZE'):
            cmd.extend(['--PRIMER_OPT_SIZE', str(args_dict['PRIMER_OPT_SIZE'])])
        if args_dict.get('PRIMER_MIN_SIZE'):
            cmd.extend(['--PRIMER_MIN_SIZE', str(args_dict['PRIMER_MIN_SIZE'])])
        if args_dict.get('PRIMER_MAX_SIZE'):
            cmd.extend(['--PRIMER_MAX_SIZE', str(args_dict['PRIMER_MAX_SIZE'])])
        if args_dict.get('PRIMER_PRODUCT_SIZE_RANGE'):
            cmd.extend(['--PRIMER_PRODUCT_SIZE_RANGE', args_dict['PRIMER_PRODUCT_SIZE_RANGE']])
        
        # BLAST arguments
        if args_dict.get('word_size'):
            cmd.extend(['--word_size', str(args_dict['word_size'])])
        if args_dict.get('reward'):
            cmd.extend(['--reward', str(args_dict['reward'])])
        if args_dict.get('penalty'):
            cmd.extend(['--penalty', str(args_dict['penalty'])])
        if args_dict.get('gapopen'):
            cmd.extend(['--gapopen', str(args_dict['gapopen'])])
        if args_dict.get('gapextend'):
            cmd.extend(['--gapextend', str(args_dict['gapextend'])])
        if args_dict.get('evalue'):
            cmd.extend(['--evalue', str(args_dict['evalue'])])
        
        # Probe check arguments
        if args_dict.get('max_mismatch'):
            cmd.extend(['--max_mismatch', str(args_dict['max_mismatch'])])
        if args_dict.get('multimap_max'):
            cmd.extend(['--multimap_max', str(args_dict['multimap_max'])])
        if args_dict.get('negative_max'):
            cmd.extend(['--negative_max', str(args_dict['negative_max'])])
        if args_dict.get('min_ident'):
            cmd.extend(['--min_ident', str(args_dict['min_ident'])])
        
        # Execute pipeline
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            cwd=os.path.join(os.path.dirname(__file__), '..')
        )
        
        if result.returncode != 0:
            return {
                'success': False,
                'error': result.stderr or result.stdout
            }
        
        # Parse results
        output_fasta = os.path.join(temp_dir, 'output.fa')
        if not os.path.exists(output_fasta):
            return {
                'success': False,
                'error': 'Output file not generated'
            }
        
        probes = parse_fasta_output(output_fasta)
        top_probes = get_top_probes(probes, 5)
        
        return {
            'success': True,
            'probes': probes,
            'top_probes': top_probes,
            'output_dir': temp_dir
        }
        
    except Exception as e:
        return {
            'success': False,
            'error': str(e)
        }

@app.route('/')
def index():
    return render_template('index.html')

@app.route('/process', methods=['POST'])
def process():
    try:
        session_id = str(uuid.uuid4())
        session['session_id'] = session_id
        
        # Handle file uploads
        args_dict = {}
        
        # Get prep_db.sh script path
        prep_db_script = os.path.join(os.path.dirname(__file__), '..', 'scripts', 'generator', 'prep_db.sh')
        if not os.path.exists(prep_db_script):
            return jsonify({'error': 'prep_db.sh script not found'}), 500
        
        # Make script executable
        os.chmod(prep_db_script, 0o755)
        
        # Required file uploads - Input FASTA (can be in archive)
        if 'input_file' not in request.files or not request.files['input_file'].filename:
            return jsonify({'error': 'Input FASTA file is required'}), 400
        
        input_file = request.files['input_file']
        input_filename = secure_filename(input_file.filename)
        input_archive_path = os.path.join(app.config['UPLOAD_FOLDER'], session_id + '_input_archive')
        input_file.save(input_archive_path)
        
        # Extract input if it's an archive
        input_extract_dir = os.path.join(app.config['UPLOAD_FOLDER'], session_id + '_input_extracted')
        if input_filename.endswith(('.zip', '.tar.gz', '.tgz', '.tar')):
            extract_archive(input_archive_path, input_extract_dir)
            # Find FASTA files in extracted directory
            fasta_files = find_fasta_files(input_extract_dir)
            if not fasta_files:
                return jsonify({'error': 'No FASTA files found in input archive'}), 400
            # Use first FASTA file or merge if multiple (for now, use first)
            args_dict['input'] = fasta_files[0]
        else:
            # Assume it's a FASTA file
            args_dict['input'] = input_archive_path
        
        # Target database database - extract FASTA and create BLAST database
        if 'true_base_file' not in request.files or not request.files['true_base_file'].filename:
            return jsonify({'error': 'Target database database is required'}), 400
        
        true_base_file = request.files['true_base_file']
        true_base_filename = secure_filename(true_base_file.filename)
        
        if not true_base_filename.endswith(('.zip', '.tar.gz', '.tgz', '.tar')):
            return jsonify({'error': 'Target database must be a .zip or .tar.gz archive containing FASTA files'}), 400
        
        true_base_archive_path = os.path.join(app.config['UPLOAD_FOLDER'], session_id + '_true_base_archive')
        true_base_file.save(true_base_archive_path)
        
        # Extract archive
        true_base_extract_dir = os.path.join(app.config['UPLOAD_FOLDER'], session_id + '_true_base_extracted')
        extract_archive(true_base_archive_path, true_base_extract_dir)
        
        # Find FASTA files
        true_base_fasta_files = find_fasta_files(true_base_extract_dir)
        if not true_base_fasta_files:
            return jsonify({'error': 'No FASTA files found in Target database archive'}), 400
        
        # Create BLAST database using prep_db.sh
        true_base_db_path = os.path.join(app.config['UPLOAD_FOLDER'], session_id + '_true_base_db')
        contig_table_path = os.path.join(app.config['UPLOAD_FOLDER'], session_id + '_contig_table.tsv')
        true_base_tmp_dir = os.path.join(app.config['UPLOAD_FOLDER'], session_id + '_true_base_tmp')
        
        run_prep_db(true_base_fasta_files, true_base_db_path, contig_table_path, 
                   true_base_tmp_dir, prep_db_script)
        
        args_dict['true_base'] = true_base_db_path
        args_dict['contig_table'] = contig_table_path
        
        # Offtarget database databases - extract FASTA and create BLAST databases
        false_base_files = request.files.getlist('false_base_files')
        if not false_base_files or not any(f.filename for f in false_base_files):
            return jsonify({'error': 'At least one offtarget database database is required'}), 400
        
        args_dict['false_base'] = []
        for i, false_base_file in enumerate(false_base_files):
            if false_base_file.filename:
                false_base_filename = secure_filename(false_base_file.filename)
                
                if not false_base_filename.endswith(('.zip', '.tar.gz', '.tgz', '.tar')):
                    return jsonify({'error': f'Offtarget database {i+1} must be a .zip or .tar.gz archive containing FASTA files'}), 400
                
                false_base_archive_path = os.path.join(app.config['UPLOAD_FOLDER'], 
                                                      session_id + f'_false_base_{i}_archive')
                false_base_file.save(false_base_archive_path)
                
                # Extract archive
                false_base_extract_dir = os.path.join(app.config['UPLOAD_FOLDER'], 
                                                     session_id + f'_false_base_{i}_extracted')
                extract_archive(false_base_archive_path, false_base_extract_dir)
                
                # Find FASTA files
                false_base_fasta_files = find_fasta_files(false_base_extract_dir)
                if not false_base_fasta_files:
                    return jsonify({'error': f'No FASTA files found in offtarget database {i+1} archive'}), 400
                
                # Create BLAST database using prep_db.sh
                false_base_db_path = os.path.join(app.config['UPLOAD_FOLDER'], 
                                                  session_id + f'_false_base_{i}_db')
                false_base_tmp_dir = os.path.join(app.config['UPLOAD_FOLDER'], 
                                                  session_id + f'_false_base_{i}_tmp')
                
                run_prep_db(false_base_fasta_files, false_base_db_path, 
                           os.path.join(app.config['UPLOAD_FOLDER'], session_id + f'_false_base_{i}_contig.tsv'),
                           false_base_tmp_dir, prep_db_script)
                
                args_dict['false_base'].append(false_base_db_path)
        
        # Optional parameters
        args_dict['threads'] = request.form.get('threads', '1')
        args_dict['algorithm'] = request.form.get('algorithm', 'FISH')
        args_dict['iterations'] = request.form.get('iterations', '10')
        args_dict['top'] = request.form.get('top', '10')
        args_dict['mutation_rate'] = request.form.get('mutation_rate', '0.05')
        args_dict['indel_rate'] = request.form.get('indel_rate', '0.05')
        args_dict['set_size'] = request.form.get('set_size', '10')
        args_dict['append'] = request.form.get('append', 'True')
        
        # Primer3 parameters
        args_dict['PRIMER_PICK_PRIMER'] = request.form.get('PRIMER_PICK_PRIMER', '10')
        args_dict['PRIMER_NUM_RETURN'] = request.form.get('PRIMER_NUM_RETURN', '10')
        args_dict['PRIMER_OPT_SIZE'] = request.form.get('PRIMER_OPT_SIZE', '25')
        args_dict['PRIMER_MIN_SIZE'] = request.form.get('PRIMER_MIN_SIZE', '15')
        args_dict['PRIMER_MAX_SIZE'] = request.form.get('PRIMER_MAX_SIZE', '30')
        args_dict['PRIMER_PRODUCT_SIZE_RANGE'] = request.form.get('PRIMER_PRODUCT_SIZE_RANGE', '100-1000')
        
        # BLAST parameters
        args_dict['word_size'] = request.form.get('word_size', '7')
        args_dict['reward'] = request.form.get('reward', '3')
        args_dict['penalty'] = request.form.get('penalty', '-3')
        args_dict['gapopen'] = request.form.get('gapopen', '6')
        args_dict['gapextend'] = request.form.get('gapextend', '3')
        args_dict['evalue'] = request.form.get('evalue', '1')
        
        # Probe check parameters
        args_dict['max_mismatch'] = request.form.get('max_mismatch', '5')
        args_dict['multimap_max'] = request.form.get('multimap_max', '1')
        args_dict['negative_max'] = request.form.get('negative_max', '0')
        args_dict['min_ident'] = request.form.get('min_ident', '70')
        
        # Run pipeline
        result = run_pipeline(args_dict, session_id)
        
        if not result['success']:
            return jsonify({'error': result['error']}), 500
        
        # Store results in session
        session['top_probes'] = result['top_probes']
        session['output_dir'] = result['output_dir']
        
        return jsonify({
            'success': True,
            'top_probes': result['top_probes'],
            'total_probes': len(result['probes'])
        })
        
    except Exception as e:
        return jsonify({'error': str(e)}), 500

@app.route('/download')
def download():
    try:
        session_id = session.get('session_id')
        if not session_id:
            return jsonify({'error': 'No session found'}), 400
        
        output_dir = session.get('output_dir')
        if not output_dir or not os.path.exists(output_dir):
            return jsonify({'error': 'Results not found'}), 404
        
        # Create zip file
        zip_path = os.path.join(app.config['RESULTS_FOLDER'], f'{session_id}_results.zip')
        with zipfile.ZipFile(zip_path, 'w', zipfile.ZIP_DEFLATED) as zipf:
            for root, dirs, files in os.walk(output_dir):
                for file in files:
                    file_path = os.path.join(root, file)
                    arcname = os.path.relpath(file_path, output_dir)
                    zipf.write(file_path, arcname)
        
        return send_file(zip_path, as_attachment=True, 
                        download_name=f'PROBEst_results_{session_id}.zip')
        
    except Exception as e:
        return jsonify({'error': str(e)}), 500

if __name__ == '__main__':
    app.run(debug=True, host='0.0.0.0', port=5000)

