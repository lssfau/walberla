import unittest
import numpy as np
import waLBerla as wlb
from waLBerla import field, createUniformBlockGrid


class FieldModuleTest(unittest.TestCase):

    def testFieldAsBlockData(self):
        blocks = createUniformBlockGrid(blocks=(1, 1, 1), cellsPerBlock=(3, 2, 2), periodic=(True, False, False))
        field.addToStorage(blocks, 'myField', np.float64, fSize=3, ghostLayers=0, initValue=0.0)
        my_field = wlb.field.toArray(blocks[0]['myField'])
        self.assertEqual(my_field[0, 0, 0, 0], 0)
        my_field[0, 0, 0, 0] = 42.0
        self.assertEqual(my_field[0, 0, 0, 0], 42.0)

        self.assertRaises(IndexError, my_field.__getitem__, (3, 0, 0))

    def testNumpyConversionWithoutGhostLayers(self):
        blocks = createUniformBlockGrid(blocks=(1, 1, 1), cellsPerBlock=(1, 2, 3), periodic=(True, False, False))
        field.addToStorage(blocks, 'f1', np.float64, fSize=4, ghostLayers=0, initValue=2.0)
        field.addToStorage(blocks, 'f2', np.float64, fSize=5, ghostLayers=0, initValue=2.0)

        f1np = field.toArray(blocks[0]['f1'])
        f2np = field.toArray(blocks[0]['f2'])
        self.assertEqual(f1np[0, 0, 0, 0], 2.0)
        self.assertEqual(f1np.shape, (1, 2, 3, 4))
        self.assertEqual(f2np.shape, (1, 2, 3, 5))

    def testGhostLayerExtraction(self):
        size = (10, 5, 4)
        gl = 3
        blocks = createUniformBlockGrid(blocks=(1, 1, 1), cellsPerBlock=size, periodic=(True, False, False))
        field.addToStorage(blocks, 'f', np.float64, fSize=3, ghostLayers=gl, initValue=0.0)

        f = blocks[0]['f']

        view1 = field.toArray(f, with_ghost_layers=True)
        self.assertEqual(view1[:, :, :, 0].shape, tuple([s + 2 * gl for s in size]))

        view2 = field.toArray(f, with_ghost_layers=False)
        self.assertEqual(view2[:, :, :, 0].shape, tuple(size))

        view3 = field.toArray(f, with_ghost_layers=2)
        self.assertEqual(view3[:, :, :, 0].shape, tuple([s + 2 * 2 for s in size]))

        view4 = field.toArray(f, with_ghost_layers=[2, False, True])
        self.assertEqual(view4[:, :, :, 0].shape, tuple([size[0] + 2 * 2, size[1] + 2 * 0, size[2] + 2 * gl]))


if __name__ == '__main__':
    unittest.main()
