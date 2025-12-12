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