from __future__ import division
import unittest
import fusion_prod.cbbundle as cbd


class MyTestCase(unittest.TestCase):
    def test_basic(self):
        liealg = cbd.TypeCLieAlgebra(1)
        wt = tuple([1])
        self.assertItemsEqual([1], liealg._convert_funds_to_epsilons(wt))
        rt = cbd._Root(liealg, [1])
        self.assertItemsEqual([2], rt)
        self.assertItemsEqual([2], liealg._convert_funds_to_epsilons(rt))
        self.assertEqual(1, liealg.killing_form(wt, rt))
        rt_list = liealg.get_positive_roots()


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

        self.assertEqual((2, -1, 0), liealg._convert_roots_to_funds((1, 0, 0)))
        self.assertEqual((-1, 2, -1), liealg._convert_roots_to_funds((0, 1, 0)))
        self.assertEqual((0, -2, 2), liealg._convert_roots_to_funds((0, 0, 1)))
        self.assertEqual((0, 1, 0), liealg._convert_roots_to_funds((1, 2, 1)))




    def test_reflection(self):
        liealg = cbd.TypeCLieAlgebra(1)
        wt = tuple([-1])
        self.assertItemsEqual([1], liealg.reflect_to_chamber(wt))

        liealg = cbd.TypeCLieAlgebra(2)



    def test_orbit(self):
        liealg = cbd.TypeCLieAlgebra(1)
        orbit = {(0,)}
        for wt in liealg.get_orbit_iter((0,)):
            self.assertTrue(wt in orbit)
            orbit.remove(wt)
        self.assertTrue(len(orbit) == 0)

        orbit_dict = {(1,), (-1,)}
        for wt in liealg.get_orbit_iter((-1,)):
            self.assertTrue(wt in orbit_dict)
        self.assertTrue(len(orbit) == 0)

        liealg = cbd.TypeCLieAlgebra(2)
        orbit = {(0, 0)}
        for wt in liealg.get_orbit_iter((0,0)):
            self.assertTrue(wt in orbit)
            orbit.remove(wt)
        self.assertTrue(len(orbit) == 0)

        orbit = {(1, 0), (-1, 0), (-1, 1), (1, -1)}
        for wt in liealg.get_orbit_iter((1,0)):
            self.assertTrue(wt in orbit)
            orbit.remove(wt)
        self.assertTrue(len(orbit) == 0)

        orbit = {(0, 1), (2, -1), (-2, 1), (0, -1)}
        for wt in liealg.get_orbit_iter((0, 1)):
            self.assertTrue(wt in orbit)
            orbit.remove(wt)
        self.assertTrue(len(orbit) == 0)

        orbit = {(1, 2), (-1, 3), (5, -2), (5, -3), (-5, 3), (-5, 2), (1, -3), (-1, -2)}
        for wt in liealg.get_orbit_iter((1, 2)):
            self.assertTrue(wt in orbit)
            orbit.remove(wt)
        self.assertTrue(len(orbit) == 0)

        liealg = cbd.TypeCLieAlgebra(3)
        orbit = {(0, 0, 0)}
        for wt in liealg.get_orbit_iter((0, 0, 0)):
            self.assertTrue(wt in orbit)
            orbit.remove(wt)
        self.assertTrue(len(orbit) == 0)

        orbit = {(1, 0, 0), (-1, 1, 0), (0, -1, 1), (0, 1, -1), (1, -1, 0), (-1, 0, 0)}
        for wt in liealg.get_orbit_iter((1, 0, 0)):
            self.assertTrue(wt in orbit)
            orbit.remove(wt)
        self.assertTrue(len(orbit) == 0)

        orbit = {(0, 1, 0), (1, -1, 1), (-1, 0, 1), (1, 1, -1), (-1, 2, -1), (2, -1, 0), (1,-2, 1), (-2, 1, 0),
                 (-1,-1, 1), (1, 0,-1), (-1, 1,-1), (0,-1, 0)}
        for wt in liealg.get_orbit_iter((0, 1, 0)):
            self.assertTrue(wt in orbit)
            orbit.remove(wt)
        self.assertTrue(len(orbit) == 0)

        orbit = {(0, 0, 1), (0, 2,-1), (2,-2, 1), (-2, 0, 1), (2, 0,-1), (-2, 2,-1), (0,-2, 1), (0, 0,-1)}
        for wt in liealg.get_orbit_iter((0, 0, 1)):
            self.assertTrue(wt in orbit)
            orbit.remove(wt)
        self.assertTrue(len(orbit) == 0)


    def test_tensor(self):
        liealg = cbd.TypeCLieAlgebra(1)
        decomp = liealg.tensor((0,), (0,))
        dec_wt = tuple([0])
        self.assertTrue(dec_wt in decomp, "Tensor decomp incorrect")
        self.assertEqual(1, decomp[dec_wt], "Tensor decomp incorrect")

        decomp = liealg.tensor((0,), (1,))
        dec_wt = tuple([1])
        self.assertTrue(dec_wt in decomp, "Tensor decomp incorrect")
        self.assertEqual(1, decomp[dec_wt], "Tensor decomp incorrect")

        decomp = liealg.tensor((1,), (1,))
        dec_wt = tuple([2])
        self.assertTrue(dec_wt in decomp, "Tensor decomp incorrect")
        self.assertEqual(1, decomp[dec_wt], "Tensor decomp incorrect")
        dec_wt = tuple([0])
        self.assertTrue(dec_wt in decomp, "Tensor decomp incorrect")
        self.assertEqual(1, decomp[dec_wt], "Tensor decomp incorrect")

        decomp = liealg.tensor((2,), (1,))
        dec_wt = tuple([3])
        self.assertTrue(dec_wt in decomp, "Tensor decomp incorrect")
        self.assertEqual(1, decomp[dec_wt], "Tensor decomp incorrect")
        dec_wt = tuple([1])
        self.assertTrue(dec_wt in decomp, "Tensor decomp incorrect")
        self.assertEqual(1, decomp[dec_wt], "Tensor decomp incorrect")

        decomp = liealg.tensor((5,), (2,))
        dec_wt = tuple([7])
        self.assertTrue(dec_wt in decomp, "Tensor decomp incorrect")
        self.assertEqual(1, decomp[dec_wt], "Tensor decomp incorrect")
        dec_wt = tuple([5])
        self.assertTrue(dec_wt in decomp, "Tensor decomp incorrect")
        self.assertEqual(1, decomp[dec_wt], "Tensor decomp incorrect")
        dec_wt = tuple([3])
        self.assertTrue(dec_wt in decomp, "Tensor decomp incorrect")
        self.assertEqual(1, decomp[dec_wt], "Tensor decomp incorrect")

        liealg = cbd.TypeCLieAlgebra(2)

        print(liealg.tensor((2, 0), (2, 0)))
        print(liealg.tensor((2, 0), (2, 0)))




if __name__ == '__main__':
    unittest.main()
