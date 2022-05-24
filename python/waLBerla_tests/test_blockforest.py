import unittest
import numpy as np
import waLBerla as wlb
from waLBerla import field, createUniformBlockGrid, AABB

field_data_type = np.float64 if wlb.build_info.build_with_double_accuracy else np.float32


class BlockforestModuleTest(unittest.TestCase):

    def testMemoryManagement1(self):
        """Testing correct reference counting of block data"""
        blocks = createUniformBlockGrid(blocks=(1, 1, 1), cellsPerBlock=(2, 2, 2))
        field.addToStorage(blocks, "TestField", field_data_type)
        f = blocks[0]["TestField"]
        strides_before = f.strides
        del blocks
        # create another block structure - this has triggered segfault
        # when previous blockstructure was already freed
        blocks = createUniformBlockGrid(blocks=(1, 1, 1), cellsPerBlock=(2, 2, 2))  # noqa: F841

        # The first block structure must exist here, since we hold a reference to block data
        # if it would have been deleted already f.strides should lead to segfault or invalid values
        self.assertEqual(strides_before, f.strides)

    def testMemoryManagement2(self):
        """Testing correct reference counting of block data
           Holding only a numpy array pointing to a waLBerla field should still hold the blockstructure alive"""
        blocks = createUniformBlockGrid(blocks=(1, 1, 1), cellsPerBlock=(2, 2, 2))
        field.addToStorage(blocks, "TestField", field_data_type)
        npf = field.toArray(blocks[0]["TestField"])
        npf[:, :, :] = 42.0
        del blocks
        # create another block structure - this has triggered segfault
        # when previous blockstructure was already freed
        blocks = createUniformBlockGrid(blocks=(1, 1, 1), cellsPerBlock=(2, 2, 2))  # noqa: F841
        self.assertEqual(npf[0, 0, 0], 42.0)

    def testMemoryManagement3(self):
        """Same as testMemoryManagement2, but with iterators"""
        blocks = createUniformBlockGrid(blocks=(1, 1, 1), cellsPerBlock=(2, 2, 2))
        field.addToStorage(blocks, "TestField", field_data_type)
        for block in blocks:
            for name in block.fieldNames:
                if name == "TestField":
                    f = block[name]
                    npf = field.toArray(f)
        npf[:, :, :] = 42.0
        del blocks, block, name, f
        blocks = createUniformBlockGrid(blocks=(1, 1, 1), cellsPerBlock=(2, 2, 2))  # noqa: F841
        self.assertEqual(npf[0, 0, 0], 42.0)

    def testExceptions(self):
        """Check that the right exceptions are thrown when nonexistent or non-convertible fields are accessed"""
        blocks = createUniformBlockGrid(blocks=(1, 1, 1), cellsPerBlock=(2, 2, 2))
        with self.assertRaises(ValueError) as cm:
            blocks[0]["cell bounding box"]
        self.assertEqual(str(cm.exception), "This blockdata is not accessible from Python")
        with self.assertRaises(IndexError) as cm:
            blocks[0]["nonexistent"]
        self.assertEqual(str(cm.exception), "No blockdata with the given name found")

    def testGeneralFunctionality(self):
        blocks = createUniformBlockGrid(blocks=(1, 1, 1), cellsPerBlock=(2, 2, 2))
        self.assertEqual(blocks.getNumberOfLevels(), 1)

        aabb = blocks.getDomain
        aabb2 = AABB(1.0, 1.0, 1.0, 1.2, 1.2, 1.2)

        self.assertEqual(aabb.min, (0.0, 0.0, 0.0))
        self.assertEqual(aabb.max, (2.0, 2.0, 2.0))
        self.assertEqual(aabb.size, (2.0, 2.0, 2.0))
        self.assertEqual(aabb.empty, False)
        self.assertEqual(aabb.volume, 8.0)
        self.assertEqual(aabb.center, (1.0, 1.0, 1.0))
        self.assertEqual(aabb.contains(aabb2), True)
        self.assertEqual(aabb2.contains(aabb), False)
        self.assertEqual(aabb2.contains((1.2, 1.2, 1.2)), False)
        self.assertEqual(aabb2.contains((1.1, 1.1, 1.1)), True)


if __name__ == '__main__':
    unittest.main()
