'''
Created on Nov 2, 2016

@author: mjschust
'''
import unittest
import conformal_blocks.cbbundle as cbd

class Test(unittest.TestCase):


    def test_sl4_basic(self):
        
        liealg = cbd.TypeALieAlgebra(3)
        self.assertEqual([1, 0, 0, 0], liealg._convert_funds_to_epsilons((1, 0, 0)))
        self.assertEqual([1, 1, 0, 0], liealg._convert_funds_to_epsilons((0, 1, 0)))
        self.assertEqual([1, 1, 1, 0], liealg._convert_funds_to_epsilons((0, 0, 1)))
        self.assertEqual([3, 3, 1, 0], liealg._convert_funds_to_epsilons((0, 2, 1)))

        self.assertEqual((1, 0, 0), liealg._convert_epsilons_to_funds([1, 0, 0, 0]))
        self.assertEqual((0, 1, 0), liealg._convert_epsilons_to_funds([1, 1, 0, 0]))
        self.assertEqual((0, 0, 1), liealg._convert_epsilons_to_funds([1, 1, 1, 0]))
        self.assertEqual((1, 0, 0), liealg._convert_epsilons_to_funds([1, 0, 0, 0]))

        rt = cbd._Root(liealg, [1, 0, 0])
        self.assertEqual((2, -1, 0), rt, "Fundamental coords incorrect")
        rt = cbd._Root(liealg, [0, 1, 0])
        self.assertEqual((-1, 2, -1), rt, "Fundamental coords incorrect")
        rt = cbd._Root(liealg, [0, 0, 1])
        self.assertEqual((0, -1, 2), rt, "Fundamental coords incorrect")
        rt = cbd._Root(liealg, [1, 1, 0])
        self.assertEqual((1,1,-1), rt, "Fundamental coords incorrect")
        self.assertEqual([1,0,-1,0], liealg._convert_funds_to_epsilons(rt), "Epsilon coords incorrect")

        wt = tuple([0, 2, 1])
        self.assertEqual(2, liealg.killing_form(wt, rt), "Killing product incorrect")
        self.assertItemsEqual([1,0,1], liealg.reflect_to_chamber(rt), "Reflection incorrect")

        self.assertEqual(4, liealg.dual_coxeter())

        roots = [(1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 1, 0), (0, 1, 1), (1, 1, 1)]
        for rt in liealg.get_positive_roots():
            self.assertTrue(tuple(rt.root_coords) in roots)
            self.assertEqual(2, liealg.length_squared(rt))
            roots.remove(tuple(rt.root_coords))
        self.assertTrue(len(roots) == 0)

        wt = liealg._convert_epsilons_to_funds([0, 0, 0, 1])
        self.assertEqual(((1, 0, 0), -1), liealg.reflect_to_chamber_with_parity(wt))

        wt = liealg._convert_epsilons_to_funds([2, -1, 1, 0])
        self.assertEqual(((1, 1, 1), 1), liealg.reflect_to_chamber_with_parity(wt))

    def test_sl4_orbit(self):

        liealg = cbd.TypeALieAlgebra(3)
        orbit = {(0, 0, 0)}
        for wt in liealg.get_orbit_iter((0, 0, 0)):
            self.assertTrue(wt in orbit)
            orbit.remove(wt)
        self.assertTrue(len(orbit) == 0)

        orbit = {(1, 0, 0), (-1, 1, 0), (0, -1, 1), (0, 0, -1)}
        for wt in liealg.get_orbit_iter((1, 0, 0)):
            self.assertTrue(wt in orbit)
            orbit.remove(wt)
        self.assertTrue(len(orbit) == 0)

        orbit = {(0, 1, 0), (1,-1, 1), (-1, 0, 1), (1, 0,-1), (-1, 1,-1), (0,-1, 0)}
        for wt in liealg.get_orbit_iter((0, 1, 0)):
            self.assertTrue(wt in orbit)
            orbit.remove(wt)
        self.assertTrue(len(orbit) == 0)

        orbit = {(0, 0, 1), (1, -1, 0), (0, 1, -1), (-1, 0, 0)}
        for wt in liealg.get_orbit_iter((0, 0, 1)):
            self.assertTrue(wt in orbit)
            orbit.remove(wt)
        self.assertTrue(len(orbit) == 0)

        orbit = {(1, 0, 1), (-1, 1, 1), (1, 1,-1), (0,-1, 2), (-1, 2,-1), (2,-1, 0), (0, 1,-2),
                 (1,-2, 1), (-2, 1, 0), (1,-1,-1), (-1,-1, 1), (-1, 0,-1)}
        for wt in liealg.get_orbit_iter((1, 0, 1)):
            self.assertTrue(wt in orbit)
            orbit.remove(wt)
        self.assertTrue(len(orbit) == 0)



        
    def test_sl2_weights(self):
        
        liealg = cbd.TypeALieAlgebra(1)
        wt = tuple([1])
        self.assertItemsEqual([1,0], liealg._convert_funds_to_epsilons(wt))
        rt = cbd._Root(liealg, [1])
        self.assertItemsEqual([2], rt)
        self.assertItemsEqual([2,0], liealg._convert_funds_to_epsilons(rt))
        self.assertEqual(1, liealg.killing_form(wt, rt),)
        wt2 = tuple([-1])
        self.assertItemsEqual([1], liealg.reflect_to_chamber(wt))
        rt_list = liealg.get_positive_roots()

        wt = tuple([0])
        orbit_wt = wt
        orbit_dict = {orbit_wt}
        for wt2 in liealg.get_orbit_iter(wt):
            self.assertTrue(wt2 in orbit_dict)

        wt = tuple([-1])
        orbit_wt1 = tuple([-1])
        orbit_wt2 = tuple([1])
        orbit_dict = {orbit_wt1, orbit_wt2}
        for wt2 in liealg.get_orbit_iter(wt):
            self.assertTrue(wt2 in orbit_dict)


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.test_SL4Weights']
    unittest.main()