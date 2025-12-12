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


import unittest
import os
import subprocess
import shutil

prep_db="./prep_db.sh"

class TestMergeFastaAndCreateBlastDB(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        """Set up test environment."""
        cls.test_dir = "test_data"
        os.makedirs(cls.test_dir, exist_ok=True)

        # Create sample FASTA files
        cls.fasta_files = [
            os.path.join(cls.test_dir, "file1.fa"),
            os.path.join(cls.test_dir, "file2.fa.gz"),
        ]
        with open(cls.fasta_files[0], "w") as f:
            f.write(">seq1\nATGC\n>seq2\nCGTA\n")
        with open(cls.fasta_files[1], "w") as f:
            f.write(">seq3\nTTAA\n>seq4\nGGCC\n")
        # Compress the second file to simulate a gzipped FASTA
        subprocess.run(["gzip", cls.fasta_files[1]])

    @classmethod
    def tearDownClass(cls):
        """Clean up test environment."""
        shutil.rmtree(cls.test_dir)

    def test_merge_fasta_and_create_blastdb(self):
        """Test the prep_db.sh script."""
        # Define paths for output files
        database_name = os.path.join(self.test_dir, "test_db")
        contig_name = os.path.join(self.test_dir, "contig_names.txt")
        tmp_dir = os.path.join(self.test_dir, "tmp")

        # Run the script
        cmd = [
            prep_db,
            "-n", database_name,
            "-c", contig_name,
            "-t", tmp_dir,
            self.fasta_files[0],
            self.fasta_files[1] + ".gz",  # Use the gzipped file
        ]
        result = subprocess.run(cmd, capture_output=True, text=True)

        # Check if the script executed successfully
        self.assertEqual(result.returncode, 0, f"Script failed with error: {result.stderr}")

        # Verify the BLAST database was created
        self.assertTrue(os.path.exists(database_name + ".nhr"), "BLAST database not created")

        # Verify the contig names file was created
        self.assertTrue(os.path.exists(contig_name), "Contig names file not created")

        # Verify the content of the contig names file
        with open(contig_name, "r") as f:
            contig_content = f.read()
        self.assertIn("file1\tseq1", contig_content)
        self.assertIn("file2\tseq3", contig_content)

        # Verify the merged FASTA file was created
        merged_fasta = os.path.join(tmp_dir, "merged.fa")
        self.assertTrue(os.path.exists(merged_fasta), "Merged FASTA file not created")

        # Verify the content of the merged FASTA file
        with open(merged_fasta, "r") as f:
            merged_content = f.read()
        self.assertIn(">seq1", merged_content)
        self.assertIn(">seq3", merged_content)

    def test_script_with_missing_arguments(self):
        """Test the script with missing arguments."""
        # Run the script without required arguments
        cmd = ["./merge_fasta_and_create_blastdb.sh"]
        result = subprocess.run(cmd, capture_output=True, text=True)

        # Check if the script failed as expected
        self.assertNotEqual(result.returncode, 0, "Script should fail with missing arguments")
        self.assertIn("error", result.stderr.lower(), "Expected error message not found")

if __name__ == "__main__":
    unittest.main()