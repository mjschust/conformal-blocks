import unittest
import conformal_blocks.cbbundle as cbd


class MyTestCase(unittest.TestCase):

    def test_A2N6(self):
        liealg = cbd.TypeALieAlgebra(2)
        cbb = cbd.ConformalBlocksBundle(liealg, [(2, 0), (0, 1), (2, 0), (1, 0), (0, 0), (0, 2)], 2)
        f_curve = [(1,), (5,), (3, 6), (2, 4)]
        self.assertEqual(0, cbb.intersect_F_curve(f_curve))

        cbb = cbd.ConformalBlocksBundle(liealg, [(0, 2), (0, 2), (0, 0), (2, 0), (0, 2), (0, 2)], 2)
        f_curve = [(4,), (2,), (3, 6), (1, 5)]
        self.assertEqual(2, cbb.intersect_F_curve(f_curve))

        cbb = cbd.ConformalBlocksBundle(liealg, [(2, 0), (1, 1), (2, 0), (0, 1), (0, 1), (2, 0)], 2)
        f_curve = [(2,), (3, 6), (5,), (1, 4)]
        self.assertEqual(0, cbb.intersect_F_curve(f_curve))

        cbb = cbd.ConformalBlocksBundle(liealg, [(0, 2), (0, 0), (0, 0), (0, 0), (0, 2), (0, 2)], 2)
        f_curve = [(2,), (3, 4, 5), (1,), (6,)]
        self.assertEqual(0, cbb.intersect_F_curve(f_curve))

        cbb = cbd.ConformalBlocksBundle(liealg, [(1, 0), (2, 0), (1, 1), (0, 2), (0, 0), (1, 0)], 2)
        f_curve = [(1, 3), (5,), (2, 4), (6,)]
        self.assertEqual(0, cbb.intersect_F_curve(f_curve))

        cbb = cbd.ConformalBlocksBundle(liealg, [(0, 1), (0, 1), (2, 0), (0, 2), (1, 0), (0, 2)], 2)
        f_curve = [(3, 6), (1,), (4,), (2, 5)]
        self.assertEqual(0, cbb.intersect_F_curve(f_curve))

        cbb = cbd.ConformalBlocksBundle(liealg, [(0, 0), (2, 0), (2, 0), (2, 0), (0, 1), (1, 1)], 2)
        f_curve = [(4,), (1, 5, 6), (2,), (3,)]
        self.assertEqual(0, cbb.intersect_F_curve(f_curve))

        cbb = cbd.ConformalBlocksBundle(liealg, [(2, 0), (0, 1), (0, 0), (2, 0), (1, 1), (0, 1)], 2)
        f_curve = [(4,), (2,), (1, 5), (3, 6)]
        self.assertEqual(0, cbb.intersect_F_curve(f_curve))

        cbb = cbd.ConformalBlocksBundle(liealg, [(0, 1), (0, 1), (1, 1), (0, 2), (1, 0), (0, 0)], 2)
        f_curve = [(2,), (4, 5), (3, 6), (1,)]
        self.assertEqual(0, cbb.intersect_F_curve(f_curve))

        cbb = cbd.ConformalBlocksBundle(liealg, [(0, 0), (1, 1), (1, 1), (1, 1), (0, 2), (0, 0)], 2)
        f_curve = [(1, 4), (2, 5), (6,), (3,)]
        self.assertEqual(0, cbb.intersect_F_curve(f_curve))

        cbb = cbd.ConformalBlocksBundle(liealg, [(1, 0), (1, 1), (1, 0), (1, 0), (1, 0), (0, 2)], 2)
        f_curve = [(1,), (5,), (4, 6), (2, 3)]
        self.assertEqual(0, cbb.intersect_F_curve(f_curve))

        cbb = cbd.ConformalBlocksBundle(liealg, [(0, 0), (1, 0), (0, 1), (0, 0), (1, 0), (1, 1)], 2)
        f_curve = [(1, 4, 6), (5,), (2,), (3,)]
        self.assertEqual(0, cbb.intersect_F_curve(f_curve))

        cbb = cbd.ConformalBlocksBundle(liealg, [(0, 0), (2, 0), (1, 0), (1, 0), (1, 1), (0, 1)], 2)
        f_curve = [(4,), (2, 3), (1, 6), (5,)]
        self.assertEqual(1, cbb.intersect_F_curve(f_curve))

        cbb = cbd.ConformalBlocksBundle(liealg, [(2, 0), (0, 2), (2, 0), (1, 1), (0, 2), (0, 0)], 2)
        f_curve = [(3, 6), (1,), (2, 5), (4,)]
        self.assertEqual(0, cbb.intersect_F_curve(f_curve))

        cbb = cbd.ConformalBlocksBundle(liealg, [(1, 1), (0, 2), (2, 0), (1, 0), (0, 1), (1, 1)], 2)
        f_curve = [(5,), (1,), (4, 6), (2, 3)]
        self.assertEqual(0, cbb.intersect_F_curve(f_curve))

        cbb = cbd.ConformalBlocksBundle(liealg, [(2, 0), (0, 0), (0, 2), (0, 0), (0, 0), (0, 2)], 2)
        f_curve = [(3,), (2, 4), (5, 6), (1,)]
        self.assertEqual(0, cbb.intersect_F_curve(f_curve))

    def test_A3N6(self):
        liealg = cbd.TypeALieAlgebra(3)
        cbb = cbd.ConformalBlocksBundle(liealg, [(0, 2, 0), (0, 1, 1), (0, 0, 1), (1, 1, 0), (0, 0, 1), (0, 1, 1)], 2)
        f_curve = [(6,), (1,), (2, 3, 4), (5,)]
        self.assertEqual(0, cbb.intersect_F_curve(f_curve))

        cbb = cbd.ConformalBlocksBundle(liealg, [(2, 0, 0), (0, 0, 2), (0, 0, 1), (2, 0, 0), (1, 0, 0), (0, 0, 1)], 2)
        f_curve = [(5, 6), (2, 3), (4,), (1,)]
        self.assertEqual(0, cbb.intersect_F_curve(f_curve))

        cbb = cbd.ConformalBlocksBundle(liealg, [(0, 1, 1), (1, 0, 1), (0, 1, 0), (1, 1, 0), (0, 0, 2), (0, 0, 1)], 2)
        f_curve = [(2, 3), (4,), (1, 6), (5,)]
        self.assertEqual(0, cbb.intersect_F_curve(f_curve))

        cbb = cbd.ConformalBlocksBundle(liealg, [(1, 0, 1), (0, 1, 1), (1, 0, 0), (1, 1, 0), (0, 0, 2), (0, 2, 0)], 2)
        f_curve = [(4,), (5,), (1, 3), (2, 6)]
        self.assertEqual(0, cbb.intersect_F_curve(f_curve))

        cbb = cbd.ConformalBlocksBundle(liealg, [(0, 0, 0), (0, 1, 1), (0, 0, 2), (0, 1, 1), (1, 0, 1), (0, 1, 1)], 2)
        f_curve = [(6,), (1,), (3, 4, 5), (2,)]
        self.assertEqual(0, cbb.intersect_F_curve(f_curve))

        cbb = cbd.ConformalBlocksBundle(liealg, [(0, 1, 0), (0, 2, 0), (0, 0, 2), (1, 0, 1), (0, 1, 1), (0, 0, 1)], 2)
        f_curve = [(1, 3), (2,), (5,), (4, 6)]
        self.assertEqual(3, cbb.intersect_F_curve(f_curve))

        cbb = cbd.ConformalBlocksBundle(liealg, [(2, 0, 0), (0, 1, 0), (1, 0, 1), (1, 0, 1), (1, 0, 0), (0, 0, 1)], 2)
        f_curve = [(3,), (4,), (1, 2), (5, 6)]
        self.assertEqual(3, cbb.intersect_F_curve(f_curve))

        cbb = cbd.ConformalBlocksBundle(liealg, [(0, 0, 1), (0, 1, 1), (1, 0, 0), (0, 0, 0), (0, 0, 0), (0, 2, 0)], 2)
        f_curve = [(1, 4), (3,), (2, 5), (6,)]
        self.assertEqual(0, cbb.intersect_F_curve(f_curve))

        cbb = cbd.ConformalBlocksBundle(liealg, [(1, 1, 0), (0, 1, 1), (2, 0, 0), (0, 0, 0), (1, 1, 0), (0, 2, 0)], 2)
        f_curve = [(2, 4), (1,), (3, 5), (6,)]
        self.assertEqual(0, cbb.intersect_F_curve(f_curve))

        cbb = cbd.ConformalBlocksBundle(liealg, [(0, 0, 1), (0, 0, 0), (0, 0, 0), (0, 0, 2), (0, 1, 1), (0, 2, 0)], 2)
        f_curve = [(3, 5), (1, 2), (6,), (4,)]
        self.assertEqual(0, cbb.intersect_F_curve(f_curve))

    def test_B3N5(self):
        # Level 1
        liealg = cbd.TypeBLieAlgebra(3)
        cbb = cbd.ConformalBlocksBundle(liealg, [(1, 0, 0), (0, 0, 1), (1, 0, 0), (0, 0, 1), (1, 0, 0)], 1)
        f_curve = [(3,), (1,), (4,), (2, 5)]
        self.assertEqual(1, cbb.intersect_F_curve(f_curve))

        cbb = cbd.ConformalBlocksBundle(liealg, [(0, 0, 1), (0, 0, 0), (0, 0, 0), (0, 0, 1), (0, 0, 1)], 1)
        f_curve = [(5,), (2,), (1, 3), (4,)]
        self.assertEqual(0, cbb.intersect_F_curve(f_curve))

        cbb = cbd.ConformalBlocksBundle(liealg, [(0, 0, 1), (0, 0, 1), (0, 0, 1), (1, 0, 0), (0, 0, 1)], 1)
        f_curve = [(3,), (1,), (4, 5), (2,)]
        self.assertEqual(2, cbb.intersect_F_curve(f_curve))

        cbb = cbd.ConformalBlocksBundle(liealg, [(0, 0, 0), (0, 0, 0), (1, 0, 0), (0, 0, 1), (0, 0, 1)], 1)
        f_curve = [(2,), (1, 5), (3,), (4,)]
        self.assertEqual(0, cbb.intersect_F_curve(f_curve))

        cbb = cbd.ConformalBlocksBundle(liealg, [(0, 0, 0), (0, 0, 1), (0, 0, 1), (0, 0, 0), (0, 0, 1)], 1)
        f_curve = [(5,), (1, 2), (3,), (4,)]
        self.assertEqual(0, cbb.intersect_F_curve(f_curve))

        cbb = cbd.ConformalBlocksBundle(liealg, [(1, 0, 0), (0, 0, 0), (0, 0, 0), (0, 0, 1), (1, 0, 0)], 1)
        f_curve = [(1,), (2, 5), (4,), (3,)]
        self.assertEqual(0, cbb.intersect_F_curve(f_curve))

        cbb = cbd.ConformalBlocksBundle(liealg, [(0, 0, 0), (1, 0, 0), (0, 0, 1), (0, 0, 0), (1, 0, 0)], 1)
        f_curve = [(4,), (5,), (2,), (1, 3)]
        self.assertEqual(0, cbb.intersect_F_curve(f_curve))

        cbb = cbd.ConformalBlocksBundle(liealg, [(1, 0, 0), (0, 0, 0), (1, 0, 0), (1, 0, 0), (0, 0, 1)], 1)
        f_curve = [(4,), (2,), (3, 5), (1,)]
        self.assertEqual(0, cbb.intersect_F_curve(f_curve))

        cbb = cbd.ConformalBlocksBundle(liealg, [(0, 0, 0), (0, 0, 0), (1, 0, 0), (0, 0, 0), (0, 0, 0)], 1)
        f_curve = [(5,), (1,), (3,), (2, 4)]
        self.assertEqual(0, cbb.intersect_F_curve(f_curve))

        cbb = cbd.ConformalBlocksBundle(liealg, [(0, 0, 0), (0, 0, 0), (1, 0, 0), (0, 0, 0), (0, 0, 1)], 1)
        f_curve = [(3,), (2,), (4, 5), (1,)]
        self.assertEqual(0, cbb.intersect_F_curve(f_curve))

        # Level 2
        cbb = cbd.ConformalBlocksBundle(liealg, [(0, 0, 2), (0, 0, 2), (0, 0, 1), (0, 0, 2), (1, 0, 1)], 2)
        f_curve = [(2,), (4,), (1, 5), (3,)]
        self.assertEqual(4, cbb.intersect_F_curve(f_curve))

        cbb = cbd.ConformalBlocksBundle(liealg, [(1, 0, 0), (1, 0, 0), (0, 0, 0), (0, 0, 1), (1, 0, 0)], 2)
        f_curve = [(5,), (2,), (1, 4), (3,)]
        self.assertEqual(0, cbb.intersect_F_curve(f_curve))

        cbb = cbd.ConformalBlocksBundle(liealg, [(0, 1, 0), (0, 1, 0), (1, 0, 0), (0, 0, 1), (2, 0, 0)], 2)
        f_curve = [(5,), (4,), (3,), (1, 2)]
        self.assertEqual(0, cbb.intersect_F_curve(f_curve))

        cbb = cbd.ConformalBlocksBundle(liealg, [(2, 0, 0), (0, 1, 0), (1, 0, 0), (2, 0, 0), (0, 0, 1)], 2)
        f_curve = [(3, 5), (4,), (1,), (2,)]
        self.assertEqual(0, cbb.intersect_F_curve(f_curve))

        cbb = cbd.ConformalBlocksBundle(liealg, [(1, 0, 1), (2, 0, 0), (0, 0, 0), (2, 0, 0), (0, 0, 0)], 2)
        f_curve = [(3,), (4,), (1, 2), (5,)]
        self.assertEqual(0, cbb.intersect_F_curve(f_curve))

        cbb = cbd.ConformalBlocksBundle(liealg, [(0, 0, 2), (0, 0, 0), (1, 0, 0), (0, 0, 1), (2, 0, 0)], 2)
        f_curve = [(3, 5), (1,), (2,), (4,)]
        self.assertEqual(0, cbb.intersect_F_curve(f_curve))

        cbb = cbd.ConformalBlocksBundle(liealg, [(0, 1, 0), (0, 0, 0), (0, 1, 0), (0, 0, 1), (1, 0, 0)], 2)
        f_curve = [(2,), (1,), (3, 5), (4,)]
        self.assertEqual(0, cbb.intersect_F_curve(f_curve))

        cbb = cbd.ConformalBlocksBundle(liealg, [(0, 0, 1), (0, 0, 1), (0, 0, 2), (2, 0, 0), (0, 0, 2)], 2)
        f_curve = [(4, 5), (1,), (3,), (2,)]
        self.assertEqual(2, cbb.intersect_F_curve(f_curve))

        cbb = cbd.ConformalBlocksBundle(liealg, [(0, 1, 0), (0, 1, 0), (0, 1, 0), (0, 0, 1), (2, 0, 0)], 2)
        f_curve = [(2, 3), (5,), (4,), (1,)]
        self.assertEqual(0, cbb.intersect_F_curve(f_curve))

        cbb = cbd.ConformalBlocksBundle(liealg, [(1, 0, 1), (1, 0, 0), (0, 0, 1), (1, 0, 1), (0, 0, 1)], 2)
        f_curve = [(5,), (1, 3), (4,), (2,)]
        self.assertEqual(3, cbb.intersect_F_curve(f_curve))

    def test_C2N6(self):
        liealg = cbd.TypeCLieAlgebra(2)
        cbb = cbd.ConformalBlocksBundle(liealg, [(1, 1), (2, 0), (0, 0), (2, 0), (1, 0), (0, 1)], 2)
        f_curve = [(6,), (2, 4), (3,), (1, 5)]
        self.assertEqual(0, cbb.intersect_F_curve(f_curve))

        cbb = cbd.ConformalBlocksBundle(liealg, [(0, 1), (2, 0), (0, 2), (0, 1), (1, 0), (0, 2)], 2)
        f_curve = [(3,), (4,), (1, 2, 6), (5,)]
        self.assertEqual(0, cbb.intersect_F_curve(f_curve))

        cbb = cbd.ConformalBlocksBundle(liealg, [(1, 1), (0, 0), (0, 2), (0, 2), (0, 2), (0, 1)], 2)
        f_curve = [(1, 4), (6,), (3, 5), (2,)]
        self.assertEqual(0, cbb.intersect_F_curve(f_curve))

        cbb = cbd.ConformalBlocksBundle(liealg, [(1, 1), (2, 0), (0, 1), (0, 0), (1, 0), (0, 0)], 2)
        f_curve = [(2, 3), (6,), (1, 5), (4,)]
        self.assertEqual(0, cbb.intersect_F_curve(f_curve))

        cbb = cbd.ConformalBlocksBundle(liealg, [(1, 0), (2, 0), (1, 0), (0, 1), (1, 0), (1, 0)], 2)
        f_curve = [(2,), (4, 6), (1,), (3, 5)]
        self.assertEqual(3, cbb.intersect_F_curve(f_curve))

        cbb = cbd.ConformalBlocksBundle(liealg, [(1, 0), (0, 0), (0, 1), (0, 1), (1, 1), (0, 0)], 2)
        f_curve = [(1,), (4, 5, 6), (3,), (2,)]
        self.assertEqual(0, cbb.intersect_F_curve(f_curve))

        cbb = cbd.ConformalBlocksBundle(liealg, [(1, 0), (1, 1), (0, 1), (1, 0), (0, 2), (0, 1)], 2)
        f_curve = [(3, 5, 6), (2,), (1,), (4,)]
        self.assertEqual(0, cbb.intersect_F_curve(f_curve))

        cbb = cbd.ConformalBlocksBundle(liealg, [(1, 0), (0, 2), (1, 0), (0, 0), (1, 1), (0, 2)], 2)
        f_curve = [(1,), (5,), (3, 4), (2, 6)]
        self.assertEqual(0, cbb.intersect_F_curve(f_curve))

        cbb = cbd.ConformalBlocksBundle(liealg, [(0, 2), (1, 0), (0, 0), (1, 1), (0, 2), (2, 0)], 2)
        f_curve = [(5,), (2,), (1,), (3, 4, 6)]
        self.assertEqual(1, cbb.intersect_F_curve(f_curve))

        cbb = cbd.ConformalBlocksBundle(liealg, [(0, 0), (1, 1), (0, 1), (1, 1), (1, 1), (0, 0)], 2)
        f_curve = [(1,), (2,), (3,), (4, 5, 6)]
        self.assertEqual(0, cbb.intersect_F_curve(f_curve))

    def test_D4N5(self):
        # Level 1
        liealg = cbd.TypeDLieAlgebra(4)
        cbb = cbd.ConformalBlocksBundle(liealg, [(0, 0, 0, 1), (1, 0, 0, 0), (0, 0, 0, 1), (0, 0, 0, 1), (0, 0, 0, 1)],
                                        1)
        f_curve = [(4,), (3,), (1, 2), (5,)]
        self.assertEqual(0, cbb.intersect_F_curve(f_curve))

        cbb = cbd.ConformalBlocksBundle(liealg, [(0, 0, 0, 0), (1, 0, 0, 0), (0, 0, 1, 0), (0, 0, 1, 0), (0, 0, 0, 1)],
                                        1)
        f_curve = [(4,), (3, 5), (1,), (2,)]
        self.assertEqual(0, cbb.intersect_F_curve(f_curve))

        cbb = cbd.ConformalBlocksBundle(liealg, [(0, 0, 0, 1), (0, 0, 1, 0), (0, 0, 0, 0), (1, 0, 0, 0), (0, 0, 0, 1)],
                                        1)
        f_curve = [(3, 4), (1,), (2,), (5,)]
        self.assertEqual(0, cbb.intersect_F_curve(f_curve))

        cbb = cbd.ConformalBlocksBundle(liealg, [(0, 0, 0, 0), (1, 0, 0, 0), (0, 0, 1, 0), (0, 0, 0, 0), (1, 0, 0, 0)],
                                        1)
        f_curve = [(2, 5), (1,), (4,), (3,)]
        self.assertEqual(0, cbb.intersect_F_curve(f_curve))

        cbb = cbd.ConformalBlocksBundle(liealg, [(0, 0, 1, 0), (0, 0, 1, 0), (0, 0, 1, 0), (1, 0, 0, 0), (0, 0, 0, 0)],
                                        1)
        f_curve = [(3, 5), (4,), (1,), (2,)]
        self.assertEqual(0, cbb.intersect_F_curve(f_curve))

        cbb = cbd.ConformalBlocksBundle(liealg, [(1, 0, 0, 0), (0, 0, 0, 1), (0, 0, 0, 1), (1, 0, 0, 0), (0, 0, 1, 0)],
                                        1)
        f_curve = [(3,), (1, 5), (4,), (2,)]
        self.assertEqual(0, cbb.intersect_F_curve(f_curve))

        cbb = cbd.ConformalBlocksBundle(liealg, [(0, 0, 0, 0), (0, 0, 0, 1), (1, 0, 0, 0), (0, 0, 0, 0), (0, 0, 1, 0)],
                                        1)
        f_curve = [(2,), (3,), (5,), (1, 4)]
        self.assertEqual(0, cbb.intersect_F_curve(f_curve))

        cbb = cbd.ConformalBlocksBundle(liealg, [(0, 0, 0, 0), (0, 0, 0, 1), (0, 0, 0, 1), (0, 0, 0, 1), (0, 0, 0, 1)],
                                        1)
        f_curve = [(5,), (1, 2), (4,), (3,)]
        self.assertEqual(2, cbb.intersect_F_curve(f_curve))

        cbb = cbd.ConformalBlocksBundle(liealg, [(0, 0, 1, 0), (0, 0, 1, 0), (0, 0, 0, 0), (0, 0, 0, 0), (0, 0, 1, 0)],
                                        1)
        f_curve = [(3,), (5,), (4,), (1, 2)]
        self.assertEqual(0, cbb.intersect_F_curve(f_curve))

        cbb = cbd.ConformalBlocksBundle(liealg, [(0, 0, 0, 1), (0, 0, 0, 1), (0, 0, 0, 1), (0, 0, 0, 1), (1, 0, 0, 0)],
                                        1)
        f_curve = [(2,), (3,), (1,), (4, 5)]
        self.assertEqual(0, cbb.intersect_F_curve(f_curve))

        # Level 2
        cbb = cbd.ConformalBlocksBundle(liealg, [(1, 0, 1, 0), (1, 0, 0, 0), (0, 0, 0, 1), (1, 0, 1, 0), (1, 0, 1, 0)],
                                        2)
        f_curve = [(3,), (1,), (4,), (2, 5)]
        self.assertEqual(0, cbb.intersect_F_curve(f_curve))

        cbb = cbd.ConformalBlocksBundle(liealg, [(1, 0, 0, 0), (0, 0, 0, 2), (0, 0, 0, 0), (0, 1, 0, 0), (0, 0, 1, 1)],
                                        2)
        f_curve = [(1, 2), (5,), (3,), (4,)]
        self.assertEqual(0, cbb.intersect_F_curve(f_curve))

        cbb = cbd.ConformalBlocksBundle(liealg, [(0, 1, 0, 0), (1, 0, 0, 1), (1, 0, 1, 0), (1, 0, 0, 1), (1, 0, 0, 1)],
                                        2)
        f_curve = [(5,), (1, 2), (3,), (4,)]
        self.assertEqual(0, cbb.intersect_F_curve(f_curve))

        cbb = cbd.ConformalBlocksBundle(liealg, [(1, 0, 0, 0), (0, 0, 0, 2), (0, 0, 1, 1), (0, 0, 1, 0), (1, 0, 0, 1)],
                                        2)
        f_curve = [(2,), (5,), (1,), (3, 4)]
        self.assertEqual(2, cbb.intersect_F_curve(f_curve))

        cbb = cbd.ConformalBlocksBundle(liealg, [(0, 1, 0, 0), (2, 0, 0, 0), (1, 0, 0, 0), (1, 0, 1, 0), (0, 0, 1, 0)],
                                        2)
        f_curve = [(2, 5), (3,), (4,), (1,)]
        self.assertEqual(2, cbb.intersect_F_curve(f_curve))

        cbb = cbd.ConformalBlocksBundle(liealg, [(0, 0, 0, 0), (2, 0, 0, 0), (1, 0, 1, 0), (0, 0, 1, 1), (0, 0, 2, 0)],
                                        2)
        f_curve = [(3,), (1,), (5,), (2, 4)]
        self.assertEqual(0, cbb.intersect_F_curve(f_curve))

        cbb = cbd.ConformalBlocksBundle(liealg, [(0, 0, 1, 1), (0, 0, 0, 1), (1, 0, 0, 0), (1, 0, 1, 0), (0, 1, 0, 0)],
                                        2)
        f_curve = [(5,), (1,), (3,), (2, 4)]
        self.assertEqual(3, cbb.intersect_F_curve(f_curve))

        cbb = cbd.ConformalBlocksBundle(liealg, [(1, 0, 0, 1), (0, 0, 0, 1), (1, 0, 0, 1), (0, 0, 1, 1), (0, 0, 2, 0)],
                                        2)
        f_curve = [(2,), (1,), (3, 4), (5,)]
        self.assertEqual(0, cbb.intersect_F_curve(f_curve))

        cbb = cbd.ConformalBlocksBundle(liealg, [(1, 0, 0, 0), (2, 0, 0, 0), (0, 0, 0, 1), (2, 0, 0, 0), (0, 0, 1, 0)],
                                        2)
        f_curve = [(2, 3), (4,), (5,), (1,)]
        self.assertEqual(1, cbb.intersect_F_curve(f_curve))

        cbb = cbd.ConformalBlocksBundle(liealg, [(0, 0, 0, 2), (1, 0, 1, 0), (0, 0, 1, 1), (0, 0, 1, 0), (0, 0, 1, 0)],
                                        2)
        f_curve = [(4,), (5,), (2, 3), (1,)]
        self.assertEqual(0, cbb.intersect_F_curve(f_curve))

if __name__ == '__main__':
    unittest.main()
