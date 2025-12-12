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

probe_check = "./probe_check.sh"

class TestProbeCheckScript(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        """Set up test environment."""
        cls.test_dir = "test_data"
        os.makedirs(cls.test_dir, exist_ok=True)

        # Create sample input files
        cls.true_file = os.path.join(cls.test_dir, "true_mappings.txt")
        cls.contig_table = os.path.join(cls.test_dir, "contig_table.txt")
        cls.filter_file1 = os.path.join(cls.test_dir, "false_mappings_1.txt")
        cls.filter_file2 = os.path.join(cls.test_dir, "false_mappings_2.txt")
        cls.output_file = os.path.join(cls.test_dir, "filtered_probes.txt")
        cls.python_filter_script = os.path.join(cls.test_dir, "filter.py")

        # Create sample true mappings file
        with open(cls.true_file, "w") as f:
            f.write("probe1\tgenome1\t90\t100\t0\t95\t1\n")
            f.write("probe2\tgenome2\t80\t90\t2\t92\t2\n")
            f.write("probe3\tgenome3\t70\t80\t5\t85\t3\n")

        # Create sample contig table
        with open(cls.contig_table, "w") as f:
            f.write("contig1\tgenome1\n")
            f.write("contig2\tgenome2\n")
            f.write("contig3\tgenome3\n")

        # Create sample false mappings files
        with open(cls.filter_file1, "w") as f:
            f.write("probe1\tgenome1\t90\t100\t0\t95\t1\n")
            f.write("probe4\tgenome4\t80\t90\t2\t92\t2\n")

        with open(cls.filter_file2, "w") as f:
            f.write("probe2\tgenome2\t80\t90\t2\t92\t2\n")
            f.write("probe5\tgenome5\t70\t80\t5\t85\t3\n")

        # Create a dummy Python filter script
        with open(cls.python_filter_script, "w") as f:
            f.write("import sys\n")
            f.write("input_file = sys.argv[1]\n")
            f.write("remove_file = sys.argv[2]\n")
            f.write("output_file = sys.argv[3]\n")
            f.write("max_out = int(sys.argv[4])\n")
            f.write("with open(input_file, 'r') as f_in, open(remove_file, 'r') as f_rm, open(output_file, 'w') as f_out:\n")
            f.write("    remove_probes = set(line.strip() for line in f_rm)\n")
            f.write("    for line in f_in:\n")
            f.write("        probe = line.split()[0]\n")
            f.write("        if probe not in remove_probes:\n")
            f.write("            f_out.write(line)\n")

    @classmethod
    def tearDownClass(cls):
        """Clean up test environment."""
        shutil.rmtree(cls.test_dir)

    def test_probe_check_script(self):
        """Test the probe_check.sh script with sample data."""
        # Run the script
        cmd = [
            "./probe_check.sh",
            "-t", self.true_file,
            "-o", self.output_file,
            "-p", self.python_filter_script,
            "-d", self.contig_table,
            "-e", "5",
            "-i", "90",
            "-a", "1",
            "-b", "1",
            self.filter_file1,
            self.filter_file2,
        ]
        result = subprocess.run(cmd, capture_output=True, text=True)

        # Check if the script executed successfully
        self.assertEqual(result.returncode, 0, f"Script failed with error: {result.stderr}")

        # Verify the output file was created
        self.assertTrue(os.path.exists(self.output_file), "Output file not created")

        # Verify the content of the output file
        with open(self.output_file, "r") as f:
            output_content = f.read()
        self.assertIn("probe1", output_content)
        self.assertNotIn("probe2", output_content)
        self.assertNotIn("probe3", output_content)

    def test_script_with_missing_arguments(self):
        """Test the script with missing arguments."""
        # Run the script without required arguments
        cmd = [probe_check]
        result = subprocess.run(cmd, capture_output=True, text=True)

        # Check if the script failed as expected
        self.assertNotEqual(result.returncode, 0, "Script should fail with missing arguments")
        self.assertIn("error", result.stderr.lower(), "Expected error message not found")

if __name__ == "__main__":
    unittest.main()