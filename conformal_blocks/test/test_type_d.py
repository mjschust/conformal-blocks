from __future__ import division
import unittest
import conformal_blocks.cbbundle as cbd

class MyTestCase(unittest.TestCase):
    def test_basic(self):
        liealg = cbd.TypeDLieAlgebra(3, exact=False)
        self.assertEqual([1, 0, 0], liealg._convert_funds_to_epsilons((1, 0, 0)))
        self.assertEqual([1/2, 1/2, -1/2], liealg._convert_funds_to_epsilons((0, 1, 0)))
        self.assertEqual([1/2, 1/2, 1/2], liealg._convert_funds_to_epsilons((0, 0, 1)))
        self.assertEqual([5/2, 3/2, -1/2], liealg._convert_funds_to_epsilons((1, 2, 1)))

        self.assertEqual((1, 0, 0), liealg._convert_epsilons_to_funds([1, 0, 0]))
        self.assertEqual((0, 1, 0), liealg._convert_epsilons_to_funds([1/2, 1/2, -1/2]))
        self.assertEqual((0, 0, 1), liealg._convert_epsilons_to_funds([1/2, 1/2, 1/2]))
        self.assertEqual((1, 2, 1), liealg._convert_epsilons_to_funds([5/2, 3/2, -1/2]))

        self.assertEqual((2, -1, -1), liealg._convert_roots_to_funds((1, 0, 0)))
        self.assertEqual((-1, 2, 0), liealg._convert_roots_to_funds((0, 1, 0)))
        self.assertEqual((-1, 0, 2), liealg._convert_roots_to_funds((0, 0, 1)))
        self.assertEqual((0, 1, 1), liealg._convert_roots_to_funds((1, 1, 1)))

        self.assertEqual(5 / 2, liealg.killing_form((1, 2, 1), (1, 0, 0)))
        self.assertEqual(1 / 4, liealg.killing_form((0, 1, 0), (0, 0, 1)))

        self.assertEqual(4, liealg.dual_coxeter())

        self.assertEqual(1, liealg.get_level((1, 0, 0)))
        self.assertEqual(1, liealg.get_level((0, 1, 0)))
        self.assertEqual(1, liealg.get_level((0, 0, 1)))
        self.assertEqual(3, liealg.get_level((1, 1, 1)))

        self.assertEqual((1, 0, 2), liealg.get_dual_weight((1, 2, 0)))

        roots = [(1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 1, 0), (1, 0, 1), (1, 1, 1)]
        for rt in liealg.get_positive_roots():
            self.assertTrue(tuple(rt.root_coords) in roots)
            roots.remove(tuple(rt.root_coords))
        self.assertTrue(len(roots) == 0)

        wts = [(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1)]
        for wt in liealg.get_weights(1):
            self.assertTrue(wt in wts)
            wts.remove(wt)
        self.assertTrue(len(roots) == 0)

        wts = [(0,0,0), (1,0,0), (2,0,0), (0,1,0), (1,1,0), (0,0,1), (1,0,1), (0,1,1), (0,2,0), (0,0,2)]
        for wt in liealg.get_weights(2):
            self.assertTrue(wt in wts)
            wts.remove(wt)
        self.assertTrue(len(roots) == 0)

        self.assertEqual((0, 0, 0), liealg.reflect_to_chamber((0, 0, 0)))
        self.assertEqual((1, 0, 1), liealg.reflect_to_chamber((1, 0, 1)))
        self.assertEqual((0, 1, 0), liealg.reflect_to_chamber((-1, 0, 1)))
        self.assertEqual((0, 1, 1), liealg.reflect_to_chamber((1, 0, -2)))
        self.assertEqual((1, 1, 0), liealg.reflect_to_chamber((1, -1, -2)))

        self.assertEqual(((0, 0, 0), 1), liealg.reflect_to_chamber_with_parity((0, 0, 0)))
        self.assertEqual(((1, 0, 1), 1), liealg.reflect_to_chamber_with_parity((1, 0, 1)))
        self.assertEqual((0, 1, 0), liealg.reflect_to_chamber((-1, 0, 1)))
        print(liealg._convert_funds_to_epsilons((-1,0,1)))
        self.assertEqual((0, 1, 1), liealg.reflect_to_chamber((1, 0, -2)))
        self.assertEqual((1, 1, 0), liealg.reflect_to_chamber((1, -1, -2)))


        liealg = cbd.TypeDLieAlgebra(4, exact=False)
        self.assertEqual([1, 0, 0, 0], liealg._convert_funds_to_epsilons((1, 0, 0, 0)))
        self.assertEqual([1, 1, 0, 0], liealg._convert_funds_to_epsilons((0, 1, 0, 0)))
        self.assertEqual([1 / 2, 1 / 2, 1 / 2, -1 / 2], liealg._convert_funds_to_epsilons((0, 0, 1, 0)))
        self.assertEqual([1 / 2, 1 / 2, 1 / 2, 1 / 2], liealg._convert_funds_to_epsilons((0, 0, 0, 1)))
        self.assertEqual([4, 3, 1, 1], liealg._convert_funds_to_epsilons((1, 2, 0, 2)))

        self.assertEqual((1, 0, 0, 0), liealg._convert_epsilons_to_funds([1, 0, 0, 0]))
        self.assertEqual((0, 1, 0, 0), liealg._convert_epsilons_to_funds([1, 1, 0, 0]))
        self.assertEqual((0, 0, 1, 0), liealg._convert_epsilons_to_funds([1 / 2, 1 / 2, 1 / 2, -1 / 2]))
        self.assertEqual((0, 0, 0, 1), liealg._convert_epsilons_to_funds([1 / 2, 1 / 2, 1 / 2, 1 / 2]))
        self.assertEqual((1, 2, 0, 2), liealg._convert_epsilons_to_funds([4, 3, 1, 1]))

        self.assertEqual((2, -1, 0, 0), liealg._convert_roots_to_funds((1, 0, 0, 0)))
        self.assertEqual((-1, 2, -1, -1), liealg._convert_roots_to_funds((0, 1, 0, 0)))
        self.assertEqual((0, -1, 2, 0), liealg._convert_roots_to_funds((0, 0, 1, 0)))
        self.assertEqual((0, -1, 0, 2), liealg._convert_roots_to_funds((0, 0, 0, 1)))
        self.assertEqual((0, 1, 0, 0), liealg._convert_roots_to_funds((1, 2, 1, 1)))

        self.assertEqual(1, liealg.killing_form((1, 0, 0, 0), (0, 1, 0, 0)))
        self.assertEqual(1 / 2, liealg.killing_form((0, 0, 1, 0), (0, 0, 0, 1)))
        self.assertEqual(4.5, liealg.killing_form((0, 0, 0, 1), (1, 2, 0, 2)))

        self.assertEqual(6, liealg.dual_coxeter())

        self.assertEqual(1, liealg.get_level((1, 0, 0, 0)))
        self.assertEqual(2, liealg.get_level((0, 1, 0, 0)))
        self.assertEqual(1, liealg.get_level((0, 0, 1, 0)))
        self.assertEqual(1, liealg.get_level((0, 0, 0, 1)))
        self.assertEqual(5, liealg.get_level((1, 1, 1, 1)))

        self.assertEqual((1, 0, 2, 0), liealg.get_dual_weight((1, 0, 2, 0)))

        roots = [(1, 0, 0, 0), (0, 1, 0, 0), (0, 0, 1, 0), (0, 0, 0, 1), (1, 1, 0, 0), (0, 1, 1, 0),
                 (0, 1, 0, 1), (1, 1, 1, 0), (1, 1, 0, 1), (0, 1, 1, 1), (1, 1, 1, 1), (1, 2, 1, 1)]
        for rt in liealg.get_positive_roots():
            self.assertTrue(tuple(rt.root_coords) in roots)
            roots.remove(tuple(rt.root_coords))
        self.assertTrue(len(roots) == 0)

        wts = [(0, 0, 0, 0), (1, 0, 0, 0), (0, 0, 1, 0), (0, 0, 0, 1)]
        for wt in liealg.get_weights(1):
            self.assertTrue(wt in wts)
            wts.remove(wt)
        self.assertTrue(len(roots) == 0)

        wts = [(0, 0, 0, 0), (0, 1, 0, 0), (1, 0, 0, 0), (2, 0, 0, 0), (0, 0, 1, 0), (1, 0, 1, 0),
               (0, 0, 0, 1), (1, 0, 0, 1), (0, 0, 1, 1), (0, 0, 2, 0), (0, 0, 0, 2)]
        for wt in liealg.get_weights(2):
            self.assertTrue(wt in wts)
            wts.remove(wt)
        self.assertTrue(len(roots) == 0)

        self.assertEqual((0, 0, 0, 0), liealg.reflect_to_chamber((0, 0, 0, 0)))
        self.assertEqual((1, 0, 1, 0), liealg.reflect_to_chamber((1, 0, 1, 0)))
        self.assertEqual((0, 0, 0, 1), liealg.reflect_to_chamber((-1, 0, 1, 0)))
        self.assertEqual((0, 0, 1, 1), liealg.reflect_to_chamber((1, 0, -2, 0)))
        self.assertEqual((0, 1, 1, 0), liealg.reflect_to_chamber((1, -1, -2, 1)))



if __name__ == '__main__':
    unittest.main()
