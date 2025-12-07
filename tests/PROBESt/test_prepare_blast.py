"""Tests for prepare_blast module."""

import pytest
import os
import sys
import tempfile
import shutil
from pathlib import Path

# Add the project root directory to the Python path
project_root = Path(__file__).parent.parent.parent
sys.path.insert(0, str(project_root))
sys.path.insert(0, str(project_root / 'src'))

from src.PROBESt.prepare_blast import (
    is_fasta_directory,
    get_fasta_files,
    prepare_blast_database,
    prepare_bases_if_needed,
    FASTA_EXTENSIONS
)


@pytest.fixture
def temp_fasta_dir(tmp_path):
    """Create a temporary directory with FASTA files."""
    fasta_dir = tmp_path / "fasta_files"
    fasta_dir.mkdir()
    
    # Create sample FASTA files
    fasta1 = fasta_dir / "seq1.fna"
    fasta1.write_text(">seq1_header\nATGCATGCATGCATGCATGC\n>seq1_header2\nCGTACGTACGTACGTACGTA\n")
    
    fasta2 = fasta_dir / "seq2.fa"
    fasta2.write_text(">seq2_header\nGGGGGAAAAATTTTTCCCCC\n")
    
    return fasta_dir


@pytest.fixture
def temp_empty_dir(tmp_path):
    """Create a temporary empty directory."""
    empty_dir = tmp_path / "empty_dir"
    empty_dir.mkdir()
    return empty_dir


@pytest.fixture
def temp_non_fasta_dir(tmp_path):
    """Create a temporary directory with non-FASTA files."""
    non_fasta_dir = tmp_path / "non_fasta"
    non_fasta_dir.mkdir()
    
    txt_file = non_fasta_dir / "file.txt"
    txt_file.write_text("Not a FASTA file")
    
    return non_fasta_dir


class TestIsFastaDirectory:
    """Tests for is_fasta_directory function."""
    
    def test_returns_true_for_directory_with_fasta_files(self, temp_fasta_dir):
        """Should return True for directory containing FASTA files."""
        assert is_fasta_directory(str(temp_fasta_dir)) is True
    
    def test_returns_false_for_empty_directory(self, temp_empty_dir):
        """Should return False for empty directory."""
        assert is_fasta_directory(str(temp_empty_dir)) is False
    
    def test_returns_false_for_directory_without_fasta_files(self, temp_non_fasta_dir):
        """Should return False for directory without FASTA files."""
        assert is_fasta_directory(str(temp_non_fasta_dir)) is False
    
    def test_returns_false_for_file_path(self, temp_fasta_dir):
        """Should return False for a file path."""
        fasta_file = temp_fasta_dir / "seq1.fna"
        assert is_fasta_directory(str(fasta_file)) is False
    
    def test_returns_false_for_nonexistent_path(self):
        """Should return False for nonexistent path."""
        assert is_fasta_directory("/nonexistent/path") is False
    
    def test_detects_gzipped_fasta_files(self, tmp_path):
        """Should detect directories with gzipped FASTA files."""
        gz_dir = tmp_path / "gzipped"
        gz_dir.mkdir()
        
        # Create a .fna.gz file (just touch it - we only check extension)
        gz_file = gz_dir / "seq.fna.gz"
        gz_file.write_bytes(b"fake gzip content")
        
        assert is_fasta_directory(str(gz_dir)) is True


class TestGetFastaFiles:
    """Tests for get_fasta_files function."""
    
    def test_returns_all_fasta_files(self, temp_fasta_dir):
        """Should return all FASTA files in directory."""
        fasta_files = get_fasta_files(str(temp_fasta_dir))
        assert len(fasta_files) == 2
        
        basenames = {os.path.basename(f) for f in fasta_files}
        assert "seq1.fna" in basenames
        assert "seq2.fa" in basenames
    
    def test_returns_empty_list_for_empty_directory(self, temp_empty_dir):
        """Should return empty list for empty directory."""
        fasta_files = get_fasta_files(str(temp_empty_dir))
        assert fasta_files == []
    
    def test_returns_sorted_list(self, temp_fasta_dir):
        """Should return sorted list of files."""
        fasta_files = get_fasta_files(str(temp_fasta_dir))
        assert fasta_files == sorted(fasta_files)
    
    def test_includes_all_supported_extensions(self, tmp_path):
        """Should find files with all supported extensions."""
        ext_dir = tmp_path / "extensions"
        ext_dir.mkdir()
        
        # Create files with different extensions
        for ext in ['.fa', '.fasta', '.fna']:
            f = ext_dir / f"file{ext}"
            f.write_text(f">header\nATGC\n")
        
        fasta_files = get_fasta_files(str(ext_dir))
        assert len(fasta_files) == 3


class TestPrepareBlastDatabase:
    """Tests for prepare_blast_database function."""
    
    def test_raises_error_for_empty_directory(self, temp_empty_dir, tmp_path):
        """Should raise ValueError for directory with no FASTA files."""
        with pytest.raises(ValueError, match="No FASTA files found"):
            prepare_blast_database(
                fasta_dir=str(temp_empty_dir),
                output_db_path=str(tmp_path / "db"),
                contig_table_path=str(tmp_path / "contigs.tsv")
            )
    
    def test_raises_error_for_missing_script(self, temp_fasta_dir, tmp_path):
        """Should raise FileNotFoundError if prep_db.sh is not found."""
        with pytest.raises(FileNotFoundError, match="prep_db.sh script not found"):
            prepare_blast_database(
                fasta_dir=str(temp_fasta_dir),
                output_db_path=str(tmp_path / "db"),
                contig_table_path=str(tmp_path / "contigs.tsv"),
                script_path="/nonexistent/path"
            )


class TestPrepareBases:
    """Tests for prepare_bases_if_needed function."""
    
    def test_returns_unchanged_paths_for_non_directories(self, tmp_path):
        """Should return unchanged paths when inputs are not FASTA directories."""
        true_base = "/path/to/blast/db"
        false_bases = ["/path/to/false1", "/path/to/false2"]
        
        result_true, result_false = prepare_bases_if_needed(
            true_base=true_base,
            false_bases=false_bases,
            output_dir=str(tmp_path),
            contig_table_path=str(tmp_path / "contigs.tsv")
        )
        
        assert result_true == true_base
        assert result_false == false_bases
    
    def test_creates_blast_db_directory(self, tmp_path):
        """Should create .blast_db directory in output."""
        true_base = "/nonexistent/path"
        false_bases = []
        
        prepare_bases_if_needed(
            true_base=true_base,
            false_bases=false_bases,
            output_dir=str(tmp_path),
            contig_table_path=str(tmp_path / "contigs.tsv")
        )
        
        blast_db_dir = tmp_path / ".blast_db"
        assert blast_db_dir.exists()


class TestFastaExtensions:
    """Tests for FASTA_EXTENSIONS constant."""
    
    def test_contains_common_extensions(self):
        """Should contain common FASTA file extensions."""
        assert '.fa' in FASTA_EXTENSIONS
        assert '.fasta' in FASTA_EXTENSIONS
        assert '.fna' in FASTA_EXTENSIONS
    
    def test_contains_gzipped_extensions(self):
        """Should contain gzipped FASTA extensions."""
        assert '.fa.gz' in FASTA_EXTENSIONS
        assert '.fasta.gz' in FASTA_EXTENSIONS
        assert '.fna.gz' in FASTA_EXTENSIONS


# Integration tests (require makeblastdb and prep_db.sh)
@pytest.mark.skipif(
    not shutil.which("makeblastdb"),
    reason="makeblastdb not available"
)
class TestIntegration:
    """Integration tests that require external tools."""
    
    def test_full_database_creation(self, temp_fasta_dir, tmp_path):
        """Should create a complete BLAST database from FASTA directory."""
        output_db = tmp_path / "test_db"
        contig_table = tmp_path / "contigs.tsv"
        
        # Get the actual script path
        script_path = project_root / "scripts" / "generator"
        
        if not (script_path / "prep_db.sh").exists():
            pytest.skip("prep_db.sh not found")
        
        result = prepare_blast_database(
            fasta_dir=str(temp_fasta_dir),
            output_db_path=str(output_db),
            contig_table_path=str(contig_table),
            script_path=str(script_path)
        )
        
        assert result == str(output_db)
        # Check that BLAST database files were created
        assert any(f.name.startswith("test_db") for f in tmp_path.iterdir())
