from __future__ import division
import unittest
import conformal_blocks.cbbundle as cbd
from fractions import Fraction


class MyTestCase(unittest.TestCase):
    def test_basic(self):
        liealg = cbd.TypeBLieAlgebra(3)
        self.assertEqual([1, 0, 0], liealg._convert_funds_to_epsilons((1, 0, 0)))
        self.assertEqual([1, 1, 0], liealg._convert_funds_to_epsilons((0, 1, 0)))
        self.assertEqual([Fraction(1,2), Fraction(1,2), Fraction(1,2)], liealg._convert_funds_to_epsilons((0, 0, 1)))
        self.assertEqual([2, 1, 0], liealg._convert_funds_to_epsilons((1, 1, 0)))

        self.assertEqual((1, 0, 0), liealg._convert_epsilons_to_funds((1, 0, 0)))
        self.assertEqual((0, 1, 0), liealg._convert_epsilons_to_funds((1, 1, 0)))
        self.assertEqual((1, -1, 2), liealg._convert_epsilons_to_funds((1, 0, 1)))
        self.assertEqual((2, 2, 0), liealg._convert_epsilons_to_funds((4, 2, 0)))

        self.assertEqual(2, liealg.killing_form((1,0,1), (0,1,0)))

        self.assertEqual(5, liealg.dual_coxeter())

        self.assertEqual(5, liealg.get_level((2,1,1)))

        roots = [(1,0,0), (0,1,0), (0,0,1), (1,1,0), (0,1,1), (1,1,1), (0,1,2), (1,1,2), (1,2,2)]
        for rt in liealg.get_positive_roots():
            self.assertTrue(tuple(rt.root_coords) in roots)
            roots.remove(tuple(rt.root_coords))
        self.assertTrue(len(roots) == 0)

        self.assertEqual((2, -1, 0), liealg._convert_roots_to_funds((1, 0, 0)))
        self.assertEqual((-1, 2, -2), liealg._convert_roots_to_funds((0, 1, 0)))
        self.assertEqual((0, -1, 2), liealg._convert_roots_to_funds((0, 0, 1)))
        self.assertEqual((1, -1, 2), liealg._convert_roots_to_funds((1, 1, 2)))

        wts = [(0,0,0),(0,0,1),(0,0,2),(0,0,3),(0,0,4),(0,1,0),(0,1,1),(0,1,2),(0,2,0),(1,0,0),(1,0,1),(1,0,2),
            (1,0,3),(1,1,0),(1,1,1),(2,0,0),(2,0,1),(2,0,2),(2,1,0),(3,0,0),(3,0,1),(4,0,0)]
        for wt in liealg.get_weights(4):
            self.assertTrue(wt in wts)
            wts.remove(wt)
        self.assertTrue(len(wts) == 0)

    def test_reflection(self):
        liealg = cbd.TypeBLieAlgebra(2)

        self.assertEqual(((3,2), -1), liealg.reflect_to_alcove_with_parity((3,4), 6))
        self.assertEqual(((1, 4), -1), liealg.reflect_to_alcove_with_parity((1, 6), 6))

    def test_fusion(self):
        liealg = cbd.TypeBLieAlgebra(2, exact=False)
        V = liealg.fusion((1,1), (1,2), 3)
        self.assertTrue(V[(0, 1)] == 1)
        self.assertTrue(V[(0, 3)] == 1)
        self.assertTrue(V[(2, 1)] == 1)
        self.assertTrue(V[(1, 1)] == 2)






if __name__ == '__main__':
    unittest.main()
