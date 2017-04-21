import unittest
import conformal_blocks.cbbundle as cbd


class MyTestCase(unittest.TestCase):

    def test_C2N6(self):
        #Level 1
        liealg = cbd.TypeCLieAlgebra(2)
        cbb = cbd.ConformalBlocksBundle(liealg, [(1, 0), (0, 0), (0, 0), (1, 0), (1, 0), (0, 1)], 1)
        self.assertEqual(0, cbb.get_rank(), "Rank incorrect: ([(1, 0), (0, 0), (0, 0), (1, 0), (1, 0), (0, 1)], 1)")

        cbb = cbd.ConformalBlocksBundle(liealg, [(0, 1), (0, 1), (0, 1), (0, 1), (0, 1), (1, 0)], 1)
        self.assertEqual(0, cbb.get_rank(), "Rank incorrect: ([(0, 1), (0, 1), (0, 1), (0, 1), (0, 1), (1, 0)], 1)")

        cbb = cbd.ConformalBlocksBundle(liealg, [(0, 0), (1, 0), (0, 0), (0, 1), (0, 0), (1, 0)], 1)
        self.assertEqual(1, cbb.get_rank(), "Rank incorrect: ([(0, 0), (1, 0), (0, 0), (0, 1), (0, 0), (1, 0)], 1)")

        cbb = cbd.ConformalBlocksBundle(liealg, [(0, 1), (0, 0), (0, 0), (0, 1), (1, 0), (0, 1)], 1)
        self.assertEqual(0, cbb.get_rank(), "Rank incorrect: ([(0, 1), (0, 0), (0, 0), (0, 1), (1, 0), (0, 1)], 1)")

        cbb = cbd.ConformalBlocksBundle(liealg, [(0, 0), (0, 0), (1, 0), (0, 1), (0, 0), (0, 1)], 1)
        self.assertEqual(0, cbb.get_rank(), "Rank incorrect: ([(0, 0), (0, 0), (1, 0), (0, 1), (0, 0), (0, 1)], 1)")

        cbb = cbd.ConformalBlocksBundle(liealg, [(1, 0), (0, 0), (0, 1), (0, 1), (0, 1), (0, 0)], 1)
        self.assertEqual(0, cbb.get_rank(), "Rank incorrect: ([(1, 0), (0, 0), (0, 1), (0, 1), (0, 1), (0, 0)], 1)")

        cbb = cbd.ConformalBlocksBundle(liealg, [(0, 0), (0, 1), (1, 0), (0, 0), (0, 1), (0, 1)], 1)
        self.assertEqual(0, cbb.get_rank(), "Rank incorrect: ([(0, 0), (0, 1), (1, 0), (0, 0), (0, 1), (0, 1)], 1)")

        cbb = cbd.ConformalBlocksBundle(liealg, [(1, 0), (1, 0), (0, 1), (0, 1), (0, 0), (0, 0)], 1)
        self.assertEqual(1, cbb.get_rank(), "Rank incorrect: ([(1, 0), (1, 0), (0, 1), (0, 1), (0, 0), (0, 0)], 1)")

        cbb = cbd.ConformalBlocksBundle(liealg, [(0, 1), (1, 0), (1, 0), (0, 0), (0, 1), (1, 0)], 1)
        self.assertEqual(0, cbb.get_rank(), "Rank incorrect: ([(0, 1), (1, 0), (1, 0), (0, 0), (0, 1), (1, 0)], 1)")

        cbb = cbd.ConformalBlocksBundle(liealg, [(1, 0), (0, 0), (1, 0), (1, 0), (1, 0), (0, 0)], 1)
        self.assertEqual(2, cbb.get_rank(), "Rank incorrect: ([(1, 0), (0, 0), (1, 0), (1, 0), (1, 0), (0, 0)], 1)")

    def test_A2N6(self):
        #Level 2
        liealg = cbd.TypeALieAlgebra(2)
        cbb = cbd.ConformalBlocksBundle(liealg, [(1, 1), (2, 0), (0, 1), (2, 0), (0, 0), (0, 2)], 2)
        self.assertEqual(0, cbb.get_rank(), "Rank incorrect: ([(1, 1), (2, 0), (0, 1), (2, 0), (0, 0), (0, 2)], 2)")

        cbb = cbd.ConformalBlocksBundle(liealg, [(0, 1), (0, 2), (0, 1), (1, 0), (0, 1), (0, 2)], 2)
        self.assertEqual(2, cbb.get_rank(), "Rank incorrect: ([(0, 1), (0, 2), (0, 1), (1, 0), (0, 1), (0, 2)], 2)")

        cbb = cbd.ConformalBlocksBundle(liealg, [(2, 0), (0, 2), (0, 1), (1, 1), (1, 1), (1, 1)], 2)
        self.assertEqual(0, cbb.get_rank(), "Rank incorrect: ([(2, 0), (0, 2), (0, 1), (1, 1), (1, 1), (1, 1)], 2)")

        cbb = cbd.ConformalBlocksBundle(liealg, [(0, 2), (2, 0), (0, 2), (2, 0), (2, 0), (1, 1)], 2)
        self.assertEqual(0, cbb.get_rank(), "Rank incorrect: ([(0, 2), (2, 0), (0, 2), (2, 0), (2, 0), (1, 1)], 2)")

        cbb = cbd.ConformalBlocksBundle(liealg, [(1, 1), (2, 0), (2, 0), (0, 2), (1, 0), (1, 0)], 2)
        self.assertEqual(0, cbb.get_rank(), "Rank incorrect: ([(1, 1), (2, 0), (2, 0), (0, 2), (1, 0), (1, 0)], 2)")

        cbb = cbd.ConformalBlocksBundle(liealg, [(0, 0), (1, 1), (1, 0), (0, 1), (0, 2), (0, 0)], 2)
        self.assertEqual(0, cbb.get_rank(), "Rank incorrect: ([(0, 0), (1, 1), (1, 0), (0, 1), (0, 2), (0, 0)], 2)")

        cbb = cbd.ConformalBlocksBundle(liealg, [(1, 0), (2, 0), (0, 0), (1, 0), (1, 1), (2, 0)], 2)
        self.assertEqual(1, cbb.get_rank(), "Rank incorrect: ([(1, 0), (2, 0), (0, 0), (1, 0), (1, 1), (2, 0)], 2)")

        cbb = cbd.ConformalBlocksBundle(liealg, [(0, 2), (0, 2), (0, 2), (0, 1), (0, 1), (0, 0)], 2)
        self.assertEqual(0, cbb.get_rank(), "Rank incorrect: ([(0, 2), (0, 2), (0, 2), (0, 1), (0, 1), (0, 0)], 2)")

        cbb = cbd.ConformalBlocksBundle(liealg, [(1, 1), (0, 1), (0, 0), (0, 1), (1, 0), (1, 0)], 2)
        self.assertEqual(3, cbb.get_rank(), "Rank incorrect: ([(1, 1), (0, 1), (0, 0), (0, 1), (1, 0), (1, 0)], 2)")

        cbb = cbd.ConformalBlocksBundle(liealg, [(1, 0), (1, 0), (2, 0), (0, 1), (1, 1), (1, 0)], 2)
        self.assertEqual(0, cbb.get_rank(), "Rank incorrect: ([(1, 0), (1, 0), (2, 0), (0, 1), (1, 1), (1, 0)], 2)")

        #Level 3
        liealg = cbd.TypeALieAlgebra(2)
        cbb = cbd.ConformalBlocksBundle(liealg, [(0, 0), (2, 0), (1, 0), (1, 1), (1, 0), (0, 1)], 3)
        self.assertEqual(4, cbb.get_rank(), "Rank incorrect: ([(0, 0), (2, 0), (1, 0), (1, 1), (1, 0), (0, 1)], 3)")

        cbb = cbd.ConformalBlocksBundle(liealg, [(1, 0), (1, 0), (2, 1), (0, 2), (1, 1), (3, 0)], 3)
        self.assertEqual(0, cbb.get_rank(), "Rank incorrect: ([(1, 0), (1, 0), (2, 1), (0, 2), (1, 1), (3, 0)], 3)")

        cbb = cbd.ConformalBlocksBundle(liealg, [(1, 2), (1, 2), (0, 1), (0, 3), (0, 1), (1, 1)], 3)
        self.assertEqual(0, cbb.get_rank(), "Rank incorrect: ([(1, 2), (1, 2), (0, 1), (0, 3), (0, 1), (1, 1)], 3)")

        cbb = cbd.ConformalBlocksBundle(liealg, [(1, 1), (0, 2), (2, 0), (2, 1), (1, 2), (0, 2)], 3)
        self.assertEqual(0, cbb.get_rank(), "Rank incorrect: ([(1, 1), (0, 2), (2, 0), (2, 1), (1, 2), (0, 2)], 3)")

        cbb = cbd.ConformalBlocksBundle(liealg, [(0, 3), (1, 0), (1, 2), (3, 0), (0, 3), (0, 2)], 3)
        self.assertEqual(0, cbb.get_rank(), "Rank incorrect: ([(0, 3), (1, 0), (1, 2), (3, 0), (0, 3), (0, 2)], 3)")

        cbb = cbd.ConformalBlocksBundle(liealg, [(1, 2), (0, 1), (1, 0), (1, 2), (0, 3), (0, 2)], 3)
        self.assertEqual(0, cbb.get_rank(), "Rank incorrect: ([(1, 2), (0, 1), (1, 0), (1, 2), (0, 3), (0, 2)], 3)")

        cbb = cbd.ConformalBlocksBundle(liealg, [(0, 1), (1, 2), (2, 0), (2, 1), (0, 3), (0, 1)], 3)
        self.assertEqual(3, cbb.get_rank(), "Rank incorrect: ([(0, 1), (1, 2), (2, 0), (2, 1), (0, 3), (0, 1)], 3)")

        cbb = cbd.ConformalBlocksBundle(liealg, [(1, 1), (1, 1), (0, 2), (0, 2), (0, 0), (1, 0)], 3)
        self.assertEqual(6, cbb.get_rank(), "Rank incorrect: ([(1, 1), (1, 1), (0, 2), (0, 2), (0, 0), (1, 0)], 3)")

        cbb = cbd.ConformalBlocksBundle(liealg, [(0, 1), (3, 0), (3, 0), (0, 0), (0, 2), (1, 0)], 3)
        self.assertEqual(0, cbb.get_rank(), "Rank incorrect: ([(0, 1), (3, 0), (3, 0), (0, 0), (0, 2), (1, 0)], 3)")

        cbb = cbd.ConformalBlocksBundle(liealg, [(0, 2), (2, 0), (3, 0), (1, 2), (0, 1), (2, 1)], 3)
        self.assertEqual(0, cbb.get_rank(), "Rank incorrect: ([(0, 2), (2, 0), (3, 0), (1, 2), (0, 1), (2, 1)], 3)")

    def test_A1N5(self):
        #Level 1
        liealg = cbd.TypeALieAlgebra(1)
        cbb = cbd.ConformalBlocksBundle(liealg, [(1,), (0,), (0,), (1,), (0,)], 1)
        self.assertEqual(1, cbb.get_rank(), "Rank incorrect: ([(1,), (0,), (0,), (1,), (0,)], 1)")

        cbb = cbd.ConformalBlocksBundle(liealg, [(1,), (0,), (1,), (0,), (0,)], 1)
        self.assertEqual(1, cbb.get_rank(), "Rank incorrect: ([(1,), (0,), (1,), (0,), (0,)], 1)")

        cbb = cbd.ConformalBlocksBundle(liealg, [(1,), (1,), (0,), (0,), (1,)], 1)
        self.assertEqual(0, cbb.get_rank(), "Rank incorrect: ([(1,), (1,), (0,), (0,), (1,)], 1)")

        cbb = cbd.ConformalBlocksBundle(liealg, [(0,), (0,), (1,), (0,), (0,)], 1)
        self.assertEqual(0, cbb.get_rank(), "Rank incorrect: ([(0,), (0,), (1,), (0,), (0,)], 1)")

        cbb = cbd.ConformalBlocksBundle(liealg, [(0,), (0,), (0,), (1,), (1,)], 1)
        self.assertEqual(1, cbb.get_rank(), "Rank incorrect: ([(0,), (0,), (0,), (1,), (1,)], 1)")

        cbb = cbd.ConformalBlocksBundle(liealg, [(1,), (0,), (0,), (1,), (0,)], 1)
        self.assertEqual(1, cbb.get_rank(), "Rank incorrect: ([(1,), (0,), (0,), (1,), (0,)], 1)")

        cbb = cbd.ConformalBlocksBundle(liealg, [(0,), (0,), (0,), (0,), (1,)], 1)
        self.assertEqual(0, cbb.get_rank(), "Rank incorrect: ([(0,), (0,), (0,), (0,), (1,)], 1)")

        cbb = cbd.ConformalBlocksBundle(liealg, [(1,), (1,), (0,), (0,), (1,)], 1)
        self.assertEqual(0, cbb.get_rank(), "Rank incorrect: ([(1,), (1,), (0,), (0,), (1,)], 1)")

        cbb = cbd.ConformalBlocksBundle(liealg, [(0,), (1,), (0,), (0,), (0,)], 1)
        self.assertEqual(0, cbb.get_rank(), "Rank incorrect: ([(0,), (1,), (0,), (0,), (0,)], 1)")

        cbb = cbd.ConformalBlocksBundle(liealg, [(0,), (0,), (1,), (1,), (0,)], 1)
        self.assertEqual(1, cbb.get_rank(), "Rank incorrect: ([(0,), (0,), (1,), (1,), (0,)], 1)")

        #Level 2
        liealg = cbd.TypeALieAlgebra(1)
        cbb = cbd.ConformalBlocksBundle(liealg, [(2,), (0,), (0,), (2,), (0,)], 2)
        self.assertEqual(1, cbb.get_rank(), "Rank incorrect: ([(2,), (0,), (0,), (2,), (0,)], 2)")

        cbb = cbd.ConformalBlocksBundle(liealg, [(1,), (1,), (1,), (1,), (2,)], 2)
        self.assertEqual(2, cbb.get_rank(), "Rank incorrect: ([(1,), (1,), (1,), (1,), (2,)], 2)")

        cbb = cbd.ConformalBlocksBundle(liealg, [(1,), (0,), (2,), (0,), (2,)], 2)
        self.assertEqual(0, cbb.get_rank(), "Rank incorrect: ([(1,), (0,), (2,), (0,), (2,)], 2)")

        cbb = cbd.ConformalBlocksBundle(liealg, [(0,), (2,), (2,), (1,), (2,)], 2)
        self.assertEqual(0, cbb.get_rank(), "Rank incorrect: ([(0,), (2,), (2,), (1,), (2,)], 2)")

        cbb = cbd.ConformalBlocksBundle(liealg, [(0,), (2,), (1,), (1,), (1,)], 2)
        self.assertEqual(0, cbb.get_rank(), "Rank incorrect: ([(0,), (2,), (1,), (1,), (1,)], 2)")

        cbb = cbd.ConformalBlocksBundle(liealg, [(2,), (2,), (1,), (2,), (1,)], 2)
        self.assertEqual(1, cbb.get_rank(), "Rank incorrect: ([(2,), (2,), (1,), (2,), (1,)], 2)")

        cbb = cbd.ConformalBlocksBundle(liealg, [(2,), (2,), (2,), (1,), (0,)], 2)
        self.assertEqual(0, cbb.get_rank(), "Rank incorrect: ([(2,), (2,), (2,), (1,), (0,)], 2)")

        cbb = cbd.ConformalBlocksBundle(liealg, [(2,), (1,), (1,), (0,), (2,)], 2)
        self.assertEqual(1, cbb.get_rank(), "Rank incorrect: ([(2,), (1,), (1,), (0,), (2,)], 2)")

        cbb = cbd.ConformalBlocksBundle(liealg, [(2,), (1,), (1,), (1,), (0,)], 2)
        self.assertEqual(0, cbb.get_rank(), "Rank incorrect: ([(2,), (1,), (1,), (1,), (0,)], 2)")

        cbb = cbd.ConformalBlocksBundle(liealg, [(0,), (0,), (2,), (2,), (2,)], 2)
        self.assertEqual(0, cbb.get_rank(), "Rank incorrect: ([(0,), (0,), (2,), (2,), (2,)], 2)")

        #Level 3
        liealg = cbd.TypeALieAlgebra(1)
        cbb = cbd.ConformalBlocksBundle(liealg, [(3,), (1,), (1,), (2,), (1,)], 3)
        self.assertEqual(2, cbb.get_rank(), "Rank incorrect: ([(3,), (1,), (1,), (2,), (1,)], 3)")

        cbb = cbd.ConformalBlocksBundle(liealg, [(0,), (0,), (0,), (3,), (2,)], 3)
        self.assertEqual(0, cbb.get_rank(), "Rank incorrect: ([(0,), (0,), (0,), (3,), (2,)], 3)")

        cbb = cbd.ConformalBlocksBundle(liealg, [(1,), (1,), (2,), (3,), (1,)], 3)
        self.assertEqual(2, cbb.get_rank(), "Rank incorrect: ([(1,), (1,), (2,), (3,), (1,)], 3)")

        cbb = cbd.ConformalBlocksBundle(liealg, [(0,), (2,), (1,), (2,), (1,)], 3)
        self.assertEqual(2, cbb.get_rank(), "Rank incorrect: ([(0,), (2,), (1,), (2,), (1,)], 3)")

        cbb = cbd.ConformalBlocksBundle(liealg, [(2,), (0,), (3,), (2,), (2,)], 3)
        self.assertEqual(0, cbb.get_rank(), "Rank incorrect: ([(2,), (0,), (3,), (2,), (2,)], 3)")

        cbb = cbd.ConformalBlocksBundle(liealg, [(0,), (1,), (0,), (1,), (2,)], 3)
        self.assertEqual(1, cbb.get_rank(), "Rank incorrect: ([(0,), (1,), (0,), (1,), (2,)], 3)")

        cbb = cbd.ConformalBlocksBundle(liealg, [(0,), (1,), (2,), (0,), (1,)], 3)
        self.assertEqual(1, cbb.get_rank(), "Rank incorrect: ([(0,), (1,), (2,), (0,), (1,)], 3)")

        cbb = cbd.ConformalBlocksBundle(liealg, [(0,), (3,), (3,), (1,), (3,)], 3)
        self.assertEqual(0, cbb.get_rank(), "Rank incorrect: ([(0,), (3,), (3,), (1,), (3,)], 3)")

        cbb = cbd.ConformalBlocksBundle(liealg, [(1,), (2,), (3,), (2,), (0,)], 3)
        self.assertEqual(1, cbb.get_rank(), "Rank incorrect: ([(1,), (2,), (3,), (2,), (0,)], 3)")

        cbb = cbd.ConformalBlocksBundle(liealg, [(0,), (3,), (3,), (1,), (2,)], 3)
        self.assertEqual(0, cbb.get_rank(), "Rank incorrect: ([(0,), (3,), (3,), (1,), (2,)], 3)")



if __name__ == '__main__':
    unittest.main()
