import unittest
from waLBerla import field, createUniformBlockGrid


class BlockforestModuleTest(unittest.TestCase):

    def testMemoryManagement1(self):
        """Testing correct reference counting of block data"""
        blocks = createUniformBlockGrid(cells=(2, 2, 2))
        field.addToStorage(blocks, "TestField", float)
        f = blocks[0]['TestField']
        stridesBefore = f.strides
        del blocks
        # create another block structure - this has triggered segfault
        # when previous blockstructure was already freed
        blocks = createUniformBlockGrid(cells=(2, 2, 2))  # noqa: F841

        # The first block structure must exist here, since we hold a reference to block data
        # if it would have been deleted already f.strides should lead to segfault or invalid values
        self.assertEqual(stridesBefore, f.strides)

    def testMemoryManagement2(self):
        """Testing correct reference counting of block data
           Holding only a numpy array pointing to a waLBerla field should still hold the blockstructure alive"""
        blocks = createUniformBlockGrid(cells=(2, 2, 2))
        field.addToStorage(blocks, "TestField", float)
        npf = field.toArray(blocks[0]['TestField'])
        npf[:, :, :, :] = 42.0
        del blocks
        # create another block structure - this has triggered segfault
        # when previous blockstructure was already freed
        blocks = createUniformBlockGrid(cells=(2, 2, 2))  # noqa: F841
        self.assertEqual(npf[0, 0, 0, 0], 42.0)


if __name__ == '__main__':
    unittest.main()
