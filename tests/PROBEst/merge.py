import unittest
import os
import subprocess
from src.PROBESt.merge import merge

class TestMergeFunctions(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        """Set up test environment."""
        cls.test_dir = "test_data"
        os.makedirs(cls.test_dir, exist_ok=True)

        # Create a sample input FASTA file for primer merging
        cls.input_fasta_primer = os.path.join(cls.test_dir, "input_primer.fasta")
        with open(cls.input_fasta_primer, "w") as f:
            f.write(">seq1_LEFT\nATGC\n>seq1_RIGHT\nCGTA\n>seq2_LEFT\nTTAA\n>seq2_RIGHT\nGGCC\n")

        # Create a sample input FASTA file for FISH merging
        cls.input_fasta_fish = os.path.join(cls.test_dir, "input_fish.fasta")
        with open(cls.input_fasta_fish, "w") as f:
            f.write(">seq1\nATGC\n>seq2\nCGTA\n")

        # Define output and temporary files
        cls.output_primer = os.path.join(cls.test_dir, "output_primer.fasta")
        cls.output_fish = os.path.join(cls.test_dir, "output_fish.fasta")
        cls.tmp_file = os.path.join(cls.test_dir, "tmp_table.txt")

        # Define the path to the bash scripts directory
        cls.script_path = "scripts/"

    @classmethod
    def tearDownClass(cls):
        """Clean up test environment."""
        os.remove(cls.input_fasta_primer)
        os.remove(cls.input_fasta_fish)
        if os.path.exists(cls.output_primer):
            os.remove(cls.output_primer)
        if os.path.exists(cls.output_fish):
            os.remove(cls.output_fish)
        if os.path.exists(cls.tmp_file):
            os.remove(cls.tmp_file)
        os.rmdir(cls.test_dir)

    def test_merge_primer_algorithm(self):
        """Test the merge function with the 'primer' algorithm."""
        # Define input parameters
        algo = "primer"
        input_file = self.input_fasta_primer
        output_file = self.output_primer
        tmp_file = self.tmp_file
        N = 5  # Number of 'N' characters to use as a separator
        script_path = self.script_path

        # Run the merge function
        merge(algo, input_file, output_file, tmp_file, N, script_path)

        # Verify the output file was created
        self.assertTrue(os.path.exists(output_file), "Output file not created")

        # Verify the content of the output file
        with open(output_file, "r") as f:
            output_content = f.read()
        expected_output = ">seq1\nATGCNNNNNCGTA\n>seq2\nTTAANNNNNGGC\n"
        self.assertEqual(output_content, expected_output, "Output content does not match expected result")

    def test_merge_fish_algorithm(self):
        """Test the merge function with the 'FISH' algorithm."""
        # Define input parameters
        algo = "FISH"
        input_file = self.input_fasta_fish
        output_file = self.output_fish
        tmp_file = self.tmp_file
        N = 5  # Number of 'N' characters (not used in FISH algorithm)
        script_path = self.script_path

        # Run the merge function
        merge(algo, input_file, output_file, tmp_file, N, script_path)

        # Verify the output file was created
        self.assertTrue(os.path.exists(output_file), "Output file not created")

        # Verify the content of the output file
        with open(output_file, "r") as f:
            output_content = f.read()
        expected_output = ">seq1\nATGC\n>seq2\nCGTA\n"
        self.assertEqual(output_content, expected_output, "Output content does not match expected result")

    def test_merge_invalid_algorithm(self):
        """Test the merge function with an invalid algorithm."""
        # Define input parameters
        algo = "invalid_algorithm"
        input_file = self.input_fasta_primer
        output_file = self.output_primer
        tmp_file = self.tmp_file
        N = 5
        script_path = self.script_path

        # Verify that the function raises a KeyError for an invalid algorithm
        with self.assertRaises(KeyError, msg="Expected KeyError for invalid algorithm"):
            merge(algo, input_file, output_file, tmp_file, N, script_path)

if __name__ == "__main__":
    unittest.main()