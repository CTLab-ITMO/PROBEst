#!/usr/bin/env python3
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


"""
PROBESt Database Search Application

This Flask app provides a web interface for searching nucleotide probes
in the combined database. Users can search by target or sequence,
and view probe properties along with 2.5D and 3D structure visualizations.
"""

import os
import sys
import csv
import re
from pathlib import Path
from flask import Flask, render_template, request, jsonify
import json

# Add parent directory to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

app = Flask(__name__)
app.secret_key = os.urandom(24)

# Database path
BASE_DIR = Path(__file__).parent.parent
DATABASE_PATH = BASE_DIR / "data" / "databases" / "full_database.csv"


def load_database():
    """Load the probe database from CSV file."""
    probes = []
    if not DATABASE_PATH.exists():
        return probes
    
    try:
        with open(DATABASE_PATH, 'r', encoding='utf-8') as f:
            reader = csv.DictReader(f)
            for row in reader:
                # Only include rows with valid sequences
                seq = row.get('Sequence', '').strip()
                if seq and seq != 'NA' and len(seq) >= 3:
                    probes.append(row)
    except Exception as e:
        print(f"Error loading database: {e}")
    
    return probes


def clean_sequence(seq: str) -> str:
    """Clean and normalize a sequence string."""
    if not seq:
        return ""
    seq = str(seq).strip().upper()
    # Remove whitespace
    seq = re.sub(r'\s+', '', seq)
    return seq


def search_by_sequence(probes, query_seq):
    """Search probes by sequence (exact or partial match)."""
    query_seq = clean_sequence(query_seq)
    if not query_seq:
        return []
    
    results = []
    for probe in probes:
        probe_seq = clean_sequence(probe.get('Sequence', ''))
        if query_seq in probe_seq or probe_seq in query_seq:
            results.append(probe)
    
    return results


def search_by_target(probes, query_target):
    """Search probes by target (name, specificity, taxonomy, etc.)."""
    query_target = query_target.lower().strip()
    if not query_target:
        return []
    
    results = []
    search_fields = ['Name', 'Specificity', 'Target rRNA', 'Taxonomy', 
                     'Category', 'Accession no.']
    
    for probe in probes:
        for field in search_fields:
            value = str(probe.get(field, '')).lower()
            if query_target in value:
                results.append(probe)
                break
    
    return results


def format_probe_for_display(probe):
    """Format probe data for display, removing NA values and organizing fields."""
    # Define important fields to display
    important_fields = [
        'id', 'Name', 'Sequence', 'Length [nt]', 'Accession no.',
        'Category', 'Specificity', 'Target rRNA', 'Taxonomy',
        'Formamide [%]', 'G+C content [%]', 'Position', 'References',
        'Remarks', 'Test', 'Used in Microarray', 'Source',
        'fluorophore', 'quencher', 'sense_antisense'
    ]
    
    formatted = {}
    for field in important_fields:
        value = probe.get(field, '')
        if value and value != 'NA' and str(value).strip():
            formatted[field] = str(value).strip()
    
    # Add any other non-NA fields
    for key, value in probe.items():
        if key not in formatted and value and value != 'NA' and str(value).strip():
            formatted[key] = str(value).strip()
    
    return formatted


def check_rnafold_available():
    """Check if RNAfold command is available."""
    try:
        import subprocess
        result = subprocess.run(
            ['which', 'RNAfold'],
            capture_output=True,
            text=True,
            timeout=2
        )
        if result.returncode == 0:
            return True
        # Try alternative: check if RNAfold is in PATH by trying to run it
        result = subprocess.run(
            ['RNAfold', '--version'],
            capture_output=True,
            text=True,
            timeout=2
        )
        return result.returncode == 0
    except Exception:
        return False


def predict_rna_structure(sequence):
    """Predict RNA secondary structure using RNAfold."""
    try:
        import subprocess
        
        # Check if RNAfold is available
        if not check_rnafold_available():
            return None, "RNAfold command not found in PATH"
        
        # Run RNAfold to get secondary structure
        input_data = f">probe\n{sequence}\n"
        result = subprocess.run(
            ['RNAfold', '--noPS', '--noLP'],
            input=input_data,
            capture_output=True,
            text=True,
            timeout=30
        )
        
        if result.returncode != 0:
            return None, f"RNAfold failed: {result.stderr}"
        
        # Parse RNAfold output
        # Format: >probe\nsequence\nstructure (energy)
        lines = [line.strip() for line in result.stdout.strip().split('\n') if line.strip()]
        
        if len(lines) >= 2:
            # Skip the header line (">probe"), get sequence and structure
            # Lines format: sequence\nstructure (energy)
            sequence_line = lines[0] if not lines[0].startswith('>') else lines[1]
            structure_line = lines[1] if not lines[0].startswith('>') else lines[2]
            
            # Extract structure (dot-bracket notation) and energy
            # Structure line format: "structure ( -12.34)"
            parts = structure_line.split()
            dot_bracket = parts[0] if parts else structure_line
            
            energy = None
            # Look for energy in parentheses
            import re
            energy_match = re.search(r'\(([-+]?\d+\.?\d*)\)', structure_line)
            if energy_match:
                energy = float(energy_match.group(1))
            
            return {
                'sequence': sequence_line,
                'structure': dot_bracket,
                'energy': energy,
                'raw_output': result.stdout
            }, None
        else:
            return None, "Unexpected RNAfold output format"
            
    except subprocess.TimeoutExpired:
        return None, "RNAfold timed out"
    except FileNotFoundError:
        return None, "RNAfold command not found. Make sure ViennaRNA is installed and in PATH"
    except Exception as e:
        return None, f"Error running RNAfold: {str(e)}"


def generate_rnaglib_visualization(sequence):
    """
    Generate RNA structure visualization using rnaglib.
    Returns data for 2.5D and 3D visualization.
    
    Note: rnaglib typically works with PDB structures. For sequences,
    we predict structure first using RNAfold, then can use rnaglib for analysis.
    """
    try:
        # Clean sequence
        seq = clean_sequence(sequence)
        if not seq or len(seq) < 3:
            return None
        
        # Convert DNA to RNA if needed (replace T with U)
        rna_seq = seq.replace('T', 'U')
        
        # Check for rnaglib
        rnaglib_available = False
        try:
            import rnaglib
            rnaglib_available = True
        except ImportError:
            pass
        
        # Try to predict structure using RNAfold
        structure_data, error_msg = predict_rna_structure(rna_seq)
        
        if structure_data:
            # Successfully predicted structure
            result = {
                'sequence': structure_data.get('sequence', rna_seq),
                'length': len(rna_seq),
                'has_structure': True,
                'structure': structure_data['structure'],
                'energy': structure_data.get('energy'),
                'method': 'RNAfold'
            }
            
            # Try to generate rnaglib graph if available
            if rnaglib_available:
                try:
                    import rnaglib
                    import networkx as nx
                    
                    # Try different methods to convert dot-bracket to graph
                    graph = None
                    try:
                        # Method 1: Try graph_from_dotbracket if available
                        from rnaglib.utils import graph_from_dotbracket
                        graph = graph_from_dotbracket(
                            structure_data['structure'],
                            structure_data.get('sequence', rna_seq)
                        )
                    except (ImportError, AttributeError):
                        try:
                            # Method 2: Try using rnaglib's graph construction
                            from rnaglib.representations import GraphConstructor
                            # This is a placeholder - actual implementation depends on rnaglib API
                            pass
                        except (ImportError, AttributeError):
                            pass
                    
                    # If we have a graph, export it
                    if graph:
                        result['rnaglib_available'] = True
                        result['method'] = 'RNAfold + rnaglib'
                        result['graph_nodes'] = graph.number_of_nodes()
                        result['graph_edges'] = graph.number_of_edges()
                        
                        # Export graph data for visualization
                        # Convert to JSON-serializable format
                        nodes = []
                        edges = []
                        for node in graph.nodes(data=True):
                            nodes.append({
                                'id': node[0],
                                'data': node[1] if len(node) > 1 else {}
                            })
                        for edge in graph.edges(data=True):
                            edges.append({
                                'source': edge[0],
                                'target': edge[1],
                                'data': edge[2] if len(edge) > 2 else {}
                            })
                        
                        result['graph_data'] = {
                            'nodes': nodes,
                            'edges': edges
                        }
                except Exception as e:
                    print(f"Error generating rnaglib graph: {e}")
                    result['rnaglib_available'] = True
                    result['rnaglib_error'] = str(e)
                    result['method'] = 'RNAfold + rnaglib (graph generation failed)'
            else:
                result['rnaglib_available'] = False
                result['note'] = 'Structure predicted. Install rnaglib for advanced visualization.'
            
            return result
        else:
            # Structure prediction failed
            result = {
                'sequence': rna_seq,
                'length': len(rna_seq),
                'has_structure': False,
                'error': error_msg or 'Structure prediction failed'
            }
            
            if rnaglib_available:
                result['rnaglib_available'] = True
                result['note'] = f'rnaglib available but structure prediction failed: {error_msg}'
            else:
                result['rnaglib_available'] = False
                result['note'] = f'Structure prediction failed: {error_msg}. Install rnaglib and ViennaRNA for full functionality.'
            
            return result
            
    except Exception as e:
        import traceback
        print(f"Error generating visualization: {e}")
        print(traceback.format_exc())
        return {
            'sequence': clean_sequence(sequence).replace('T', 'U'),
            'length': len(clean_sequence(sequence)),
            'has_structure': False,
            'error': str(e)
        }


# Load database on startup
PROBE_DATABASE = load_database()
print(f"Loaded {len(PROBE_DATABASE)} probes from database")


@app.route('/')
def index():
    """Main page with search interface."""
    return render_template('database.html')


@app.route('/search', methods=['POST'])
def search():
    """Search for probes by sequence or target."""
    try:
        data = request.get_json()
        query = data.get('query', '').strip()
        search_type = data.get('type', 'sequence')  # 'sequence' or 'target'
        
        if not query:
            return jsonify({'error': 'Query is required'}), 400
        
        # Perform search
        if search_type == 'sequence':
            results = search_by_sequence(PROBE_DATABASE, query)
        else:
            results = search_by_target(PROBE_DATABASE, query)
        
        # Limit results to top 50
        results = results[:50]
        
        # Format results
        formatted_results = []
        for probe in results:
            formatted = format_probe_for_display(probe)
            formatted_results.append(formatted)
        
        return jsonify({
            'success': True,
            'count': len(formatted_results),
            'results': formatted_results
        })
        
    except Exception as e:
        return jsonify({'error': str(e)}), 500


@app.route('/probe/<probe_id>', methods=['GET'])
def get_probe(probe_id):
    """Get detailed information about a specific probe."""
    try:
        # Search for probe by ID or sequence
        probe = None
        for p in PROBE_DATABASE:
            if p.get('id') == probe_id or p.get('Sequence', '').strip() == probe_id:
                probe = p
                break
        
        if not probe:
            return jsonify({'error': 'Probe not found'}), 404
        
        formatted = format_probe_for_display(probe)
        sequence = formatted.get('Sequence', '')
        
        # Generate visualization
        visualization = generate_rnaglib_visualization(sequence)
        
        return jsonify({
            'success': True,
            'probe': formatted,
            'visualization': visualization
        })
        
    except Exception as e:
        return jsonify({'error': str(e)}), 500


@app.route('/visualize', methods=['POST'])
def visualize():
    """Generate visualization for a sequence."""
    try:
        data = request.get_json()
        sequence = data.get('sequence', '').strip()
        
        if not sequence:
            return jsonify({'error': 'Sequence is required'}), 400
        
        visualization = generate_rnaglib_visualization(sequence)
        
        if not visualization:
            return jsonify({'error': 'Failed to generate visualization'}), 500
        
        return jsonify({
            'success': True,
            'visualization': visualization
        })
        
    except Exception as e:
        return jsonify({'error': str(e)}), 500


if __name__ == '__main__':
    app.run(debug=True, host='0.0.0.0', port=5001)
