import unittest
import waLBerla as wlb


class CoreTest(unittest.TestCase):

    def test_CellInterval(self):
        ci1 = wlb.CellInterval(0, 0, 0, 5, 5, 5)
        ci2 = wlb.CellInterval([0] * 3, [5] * 3)
        self.assertEqual(ci1, ci2, "Equality comparison of CellIntervals failed.")
        self.assertFalse(ci1 != ci2, "Inequality check for CellIntervals wrong ")

        self.assertEqual(ci1.min, (0, 0, 0), "CellInterval min wrong")
        self.assertEqual(ci1.max, (5, 5, 5), "CellInterval max wrong")

        self.assertFalse(ci1.empty())

        ci1.intersect(ci2)
        self.assertTrue(ci1.contains(ci2))

        ci2.expand(1)
        self.assertFalse(ci1.contains(ci2))

    def test_AABB(self):
        aabb1 = wlb.AABB(0, 0, 0, 5, 5, 5)
        aabb2 = wlb.AABB([0] * 3, [5] * 3)
        self.assertEqual(aabb1, aabb2)


if __name__ == '__main__':
    unittest.main()
