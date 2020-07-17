import unittest
from waLBerla import field, createUniformBlockGrid


class FieldModuleTest(unittest.TestCase):

    def testFieldAsBlockData(self):
        blocks = createUniformBlockGrid(cells=(3, 2, 2), periodic=(1, 0, 0))
        field.addToStorage(blocks, 'myField', float, fSize=3, ghostLayers=0, initValue=0.0)
        myField = blocks[0]['myField']
        self.assertEqual(myField[0, 0, 0, 0], 0)
        myField[0, 0, 0, 0] = 42.0
        self.assertEqual(myField[0, 0, 0, 0], 42.0)

        self.assertRaises(IndexError, myField.__getitem__, (3, 0, 0))

    def testNumpyConversionWithoutGhostLayers(self):
        f1 = field.createField([1, 2, 3, 4], float, 2, field.zyxf)
        f2 = field.createField([1, 2, 3, 5], float, 4, field.zyxf)
        f1np = field.toArray(f1)
        f2np = field.toArray(f2)
        self.assertEqual(f1np[0, 0, 0, 0], 0)
        self.assertEqual(f1np.shape, (1, 2, 3, 4))
        self.assertEqual(f2np.shape, (1, 2, 3, 5))

        f1np[0, 0, 0, 0] = 1
        f2np[0, 0, 0, 0] = 2
        self.assertEqual(f1[0, 0, 0, 0], 1)
        self.assertEqual(f2[0, 0, 0, 0], 2)

    def testNumpyConversionWithGhostLayers(self):
        f = field.createField([1, 2, 3, 1], float, 2, field.zyxf)
        fnp = field.toArray(f, withGhostLayers=True)

        self.assertEqual(fnp[0, 0, 0, 0], 0)
        self.assertEqual(fnp.shape, (1 + 4, 2 + 4, 3 + 4, 1))
        fnp[0, 0, 0, 0] = 42
        self.assertEqual(f[-2, -2, -2, 0], 42)

    def testGhostLayerExtraction(self):
        size = [10, 5, 4]
        gl = 3
        f = field.createField(size, float, ghostLayers=gl)

        view1 = field.toArray(f, withGhostLayers=True)
        self.assertEqual(view1[:, :, :, 0].shape, tuple([s + 2 * gl for s in size]))

        view2 = field.toArray(f, withGhostLayers=False)
        self.assertEqual(view2[:, :, :, 0].shape, tuple(size))

        view3 = field.toArray(f, withGhostLayers=2)
        self.assertEqual(view3[:, :, :, 0].shape, tuple([s + 2 * 2 for s in size]))

        view4 = field.toArray(f, withGhostLayers=[2, False, True])
        self.assertEqual(view4[:, :, :, 0].shape, tuple([size[0] + 2 * 2, size[1] + 2 * 0, size[2] + 2 * gl]))


if __name__ == '__main__':
    unittest.main()
