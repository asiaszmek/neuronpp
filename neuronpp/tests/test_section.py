import re
import unittest
from neuron import h

from neuronpp.cells.cell import Cell


class TestSection(unittest.TestCase):
    def setUp(self):
        self.cell = Cell(name="cell")

    def tearDown(self):
        self.cell.remove_immediate_from_neuron()

        l = len(list(h.allsec()))
        if len(list(h.allsec())) != 0:
            raise RuntimeError("Not all section have been removed after teardown. "
                               "Sections left: %s" % l)

    def test_default_name_convetion(self):
        """
        If the name exists it will add a number starting from 2 to the proposed name.
        """
        pat = re.compile("[1-9]")

        for i in range(10):
            self.cell.add_sec("dend")

        for i, dend in enumerate(self.cell.filter_secs(name="dend")):
            res = pat.search(dend.name)
            if i == 0:
                self.assertIsNone(res)
            else:
                self.assertEqual(i+1, int(dend.name.replace("dend", "")))
            dend.remove_immediate_from_neuron()

    def test_regex_search(self):
        for i in range(10):
            self.cell.add_sec("dend")

        for i, dend in enumerate(self.cell.filter_secs("regex:dend[1-9]")):
            self.assertEqual(i+2, int(dend.name.replace("dend", "")))

        for d in self.cell.secs:
            d.remove_immediate_from_neuron()


if __name__ == '__main__':
    unittest.main()

