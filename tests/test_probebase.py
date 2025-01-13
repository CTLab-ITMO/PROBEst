import unittest
import pandas as pd
from scripts import parse_probebase_page

class TestProbebase(unittest.TestCase):
    def test_responce(self):
        """
        Test page with responce problem
        """
        data = "https://probebase.csb.univie.ac.at/pb_report/probe"
        resp = parse_probebase_page(data)
        self.assertIsInstance(resp, Warning)

    def test_responce_empty(self):
        """
        Test page with empty table
        """
        data = "https://probebase.csb.univie.ac.at/pb_report/probe/1"
        resp = parse_probebase_page(data)
        self.assertIsInstance(resp, Warning)

    def test_responce_empty(self):
        """
        Test page with content
        """
        data = "https://probebase.csb.univie.ac.at/pb_report/probe/2"
        resp = parse_probebase_page(data)
        self.assertIsInstance(resp, pd.Series)

if __name__ == '__main__':
    unittest.main()