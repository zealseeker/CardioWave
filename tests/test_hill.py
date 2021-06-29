import unittest
import warnings
import pandas as pd
from cdwave import hillcurve


class TestHillCurve(unittest.TestCase):
    cs = [-2, -1, 0, 1, 2]
    vs = [100, 90, 50, 10, 0]
    data = pd.DataFrame({'concentration': cs, 'response': vs})
    def test_tcpl_hill(self):
        hill = hillcurve.TCPLHill(self.data['concentration'], self.data['response'], logit=False)
        self.assertEqual(len(hill.popt), 4)

    def test_tcpl_gainloss(self):
        hill = hillcurve.TCPLGainLoss(self.data['concentration'], self.data['response'], logit=False)
        self.assertEqual(len(hill.popt), 7)

    def test_tcpl_plain(self):
        hill = hillcurve.TCPLPlain(self.data['concentration'], self.data['response'], logit=False)
        self.assertEqual(len(hill.popt), 1)
