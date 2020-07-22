import unittest
import waLBerla as wlb
from waLBerla import field
import numpy as np


class BasicDatatypesTest(unittest.TestCase):

    def test_numpyConversion(self):
        f = wlb.field.create([3, 3, 3], ghostLayers=2)
        a = np.asarray(f.buffer())
        self.assertEqual(a.shape, (3, 3, 3, 1))

        f = field.create([3, 3, 3], ghostLayers=2)
        f.bufferWithGhostLayers()
        b = np.asarray(f.buffer(True))
        self.assertEqual(b.shape, (7, 7, 7, 1))

    def test_swapDataPointers(self):
        f1 = wlb.field.createField([30] * 3, float)
        f2 = wlb.field.createField([30] * 3, float)

        nf1 = np.asarray(f1.buffer())
        nf2 = np.asarray(f2.buffer())

        nf1[:] = 1
        nf2[:] = 2

        del nf2

        f1.swapDataPointers(f2)

        del f2

        # Should free part1
        del nf1

        # should free part2
        del f1
