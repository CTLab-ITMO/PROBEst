import unittest
import os
import subprocess

class TestFasta2TableScript(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        """Set up test environment."""
        cls.test_dir = "test_data"
        os.makedirs(cls.test_dir, exist_ok=True)

        # Create sample input FASTA file
        cls.input_fasta = os.path.join(cls.test_dir, "input.fasta")
        cls.output_table = os.path.join(cls.test_dir, "output_table.txt")

        with open(cls.input_fasta, "w") as f:
            f.write(">header1\nATGC\n>header2\nCGTA\n>header3\nTTAA\n")

    @classmethod
    def tearDownClass(cls):
        """Clean up test environment."""
        os.remove(cls.input_fasta)
        if os.path.exists(cls.output_table):
            os.remove(cls.output_table)
        os.rmdir(cls.test_dir)

    def test_fasta2table_script(self):
        """Test the fasta2table.sh script with sample data."""
        # Run the script
        cmd = [
            "./fasta2table.sh",
            "-i", self.input_fasta,
            "-o", self.output_table,
        ]
        result = subprocess.run(cmd, capture_output=True, text=True)

        # Check if the script executed successfully
        self.assertEqual(result.returncode, 0, f"Script failed with error: {result.stderr}")

        # Verify the output file was created
        self.assertTrue(os.path.exists(self.output_table), "Output file not created")

        # Verify the content of the output file
        with open(self.output_table, "r") as f:
            output_content = f.read()
        expected_output = ">header1\tATGC\n>header2\tCGTA\n>header3\tTTAA\n"
        self.assertEqual(output_content, expected_output, "Output content does not match expected result")

    def test_script_with_missing_arguments(self):
        """Test the script with missing arguments."""
        # Run the script without required arguments
        cmd = ["./fasta2table.sh"]
        result = subprocess.run(cmd, capture_output=True, text=True)

        # Check if the script failed as expected
        self.assertNotEqual(result.returncode, 0, "Script should fail with missing arguments")
        self.assertIn("error", result.stderr.lower(), "Expected error message not found")

if __name__ == "__main__":
    unittest.main()