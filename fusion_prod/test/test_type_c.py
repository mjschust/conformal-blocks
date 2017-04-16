from __future__ import division
import unittest
import fusion_prod.cbd as cbd


class MyTestCase(unittest.TestCase):
    def test_basic(self):
        liealg = cbd.TypeCLieAlgebra(3)
        self.assertEqual([1, 0, 0], liealg._convert_funds_to_epsilons((1, 0, 0)))
        self.assertEqual([1, 1, 0], liealg._convert_funds_to_epsilons((0, 1, 0)))
        self.assertEqual([1, 1, 1], liealg._convert_funds_to_epsilons((0, 0, 1)))
        self.assertEqual([2, 1, 0], liealg._convert_funds_to_epsilons((1, 1, 0)))

        self.assertEqual((1, 0, 0), liealg._convert_epsilons_to_funds((1, 0, 0)))
        self.assertEqual((0, 1, 0), liealg._convert_epsilons_to_funds((1, 1, 0)))
        self.assertEqual((1, -1, 1), liealg._convert_epsilons_to_funds((1, 0, 1)))
        self.assertEqual((2, 2, 0), liealg._convert_epsilons_to_funds((4, 2, 0)))

        self.assertEqual(1.5, liealg.killing_form((1,0,1), (0,1,0)))

        self.assertEqual(4, liealg.dual_coxeter())

        self.assertEqual(4, liealg.get_level((2,1,1)))

        roots = [(1,0,0), (0,1,0), (0,0,1), (1,1,0), (0,1,1), (1,1,1), (0,2,1), (1,2,1), (2,2,1)]
        for rt in liealg.get_positive_roots():
            self.assertTrue(tuple(rt.root_coords) in roots)

        self.assertEqual([2, -1, 0], liealg._convert_roots_to_funds((1, 0, 0)))
        self.assertEqual([-1, 2, -1], liealg._convert_roots_to_funds((0, 1, 0)))
        self.assertEqual([0, -2, 2], liealg._convert_roots_to_funds((0, 0, 1)))
        self.assertEqual([0, 1, 0], liealg._convert_roots_to_funds((1, 2, 1)))




    def test_reflection(self):
        liealg = cbd.TypeCLieAlgebra(3)

        print(liealg.reflect_to_chamber_with_parity(liealg._convert_epsilons_to_funds([-1,1,4])))

        print(liealg.reflect_to_alcove_with_parity(liealg._convert_epsilons_to_funds([-1, 1, 4]), 2))

    def test_orbit(self):
        liealg = cbd.TypeCLieAlgebra(3)

        for wt in liealg.get_orbit_iter((0,0,0)):
            print(liealg._convert_funds_to_epsilons(wt))

    def test_tensor(self):
        liealg = cbd.TypeCLieAlgebra(3)

        print(liealg.fusion((1,1,0), (0,1,0), 2))



if __name__ == '__main__':
    unittest.main()
