import unittest
import os
from Bio import SeqIO
from src.PROBESt.primer3 import primer_template, parse_primer3_output


class TestPrimer3Functions(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        """Set up test environment."""
        cls.test_dir = "test_data"
        os.makedirs(cls.test_dir, exist_ok=True)

        # Create a sample FASTA file
        cls.fasta_file = os.path.join(
            cls.test_dir, "test_sequence.exon.mod.fna")
        with open(cls.fasta_file, "w") as f:
            f.write(">seq1\nATGCATGCATGC\n>seq2\nCGTAACGTAACG\n")

        # Create a sample Primer3 output file
        cls.primer3_output_file = os.path.join(
            cls.test_dir, "primer3_output.txt")
        with open(cls.primer3_output_file, "w") as f:
            f.write("""SEQUENCE_ID=test_sequence_seq1
PRIMER_LEFT_0_SEQUENCE=ATGC
PRIMER_RIGHT_0_SEQUENCE=CGTA
SEQUENCE_ID=test_sequence_seq2
PRIMER_LEFT_0_SEQUENCE=TTAA
PRIMER_RIGHT_0_SEQUENCE=GGCC
""")

    @classmethod
    def tearDownClass(cls):
        """Clean up test environment."""
        os.remove(cls.fasta_file)
        os.remove(cls.primer3_output_file)
        os.rmdir(cls.test_dir)

    def test_primer_template(self):
        """Test the primer_template function."""
        # Define input parameters
        PRIMER_PICK_PRIMER = 1
        PRIMER_OPT_SIZE = 20
        PRIMER_MIN_SIZE = 18
        PRIMER_MAX_SIZE = 22
        PRIMER_PRODUCT_SIZE_RANGE = "100-300"
        PRIMER_NUM_RETURN = 1

        # Generate the primer template
        result = primer_template(
            self.fasta_file,
            PRIMER_PICK_PRIMER,
            PRIMER_OPT_SIZE,
            PRIMER_MIN_SIZE,
            PRIMER_MAX_SIZE,
            PRIMER_PRODUCT_SIZE_RANGE,
            PRIMER_NUM_RETURN
        )

        # Expected output
        expected_output = """SEQUENCE_ID=test_sequence_seq1
SEQUENCE_TEMPLATE=ATGCATGCATGC
PRIMER_TASK=generic
PRIMER_PICK_LEFT_PRIMER=1
PRIMER_PICK_RIGHT_PRIMER=1
PRIMER_PICK_INTERNAL_OLIGO=0
PRIMER_OPT_SIZE=20
PRIMER_MIN_SIZE=18
PRIMER_MAX_SIZE=22
PRIMER_PRODUCT_SIZE_RANGE=100-300
PRIMER_NUM_RETURN=1
PRIMER_EXPLAIN_FLAG=1
=
SEQUENCE_ID=test_sequence_seq2
SEQUENCE_TEMPLATE=CGTAACGTAACG
PRIMER_TASK=generic
PRIMER_PICK_LEFT_PRIMER=1
PRIMER_PICK_RIGHT_PRIMER=1
PRIMER_PICK_INTERNAL_OLIGO=0
PRIMER_OPT_SIZE=20
PRIMER_MIN_SIZE=18
PRIMER_MAX_SIZE=22
PRIMER_PRODUCT_SIZE_RANGE=100-300
PRIMER_NUM_RETURN=1
PRIMER_EXPLAIN_FLAG=1
=
"""

        # Verify the output matches the expected result
        self.assertEqual(result, expected_output,
                         "Primer template output does not match expected result")

    def test_parse_primer3_output(self):
        """Test the parse_primer3_output function."""
        # Parse the Primer3 output file
        result = parse_primer3_output(self.primer3_output_file)

        # Expected output
        expected_output = [
            ("test_sequence_seq1", "0", "LEFT", "ATGC"),
            ("test_sequence_seq1", "0", "RIGHT", "CGTA"),
            ("test_sequence_seq2", "0", "LEFT", "TTAA"),
            ("test_sequence_seq2", "0", "RIGHT", "GGCC"),
        ]

        # Verify the parsed output matches the expected result
        self.assertEqual(result, expected_output,
                         "Parsed Primer3 output does not match expected result")


if __name__ == "__main__":
    unittest.main()
