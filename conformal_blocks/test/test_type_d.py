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
        self.assertEqual((0, 0, 1, 1), liealg.reflect_to_chamber((-1, 2, -2, 0)))
        self.assertEqual((0, 1, 1, 0), liealg.reflect_to_chamber((1, -1, -2, 1)))

    def test_orbit(self):
        liealg = cbd.TypeDLieAlgebra(3, exact=False)

        orbit = {(0, 0, 0)}
        for wt in liealg.get_orbit_iter((0, 0, 0)):
            self.assertTrue(wt in orbit)
            orbit.remove(wt)
        self.assertTrue(len(orbit) == 0)

        orbit = {(1, 0, 0), (-1, 1, 1), (0, -1, 1), (0, 1, -1), (1, -1, -1), (-1, 0, 0)}
        for wt in liealg.get_orbit_iter((1, 0, 0)):
            self.assertTrue(wt in orbit)
            orbit.remove(wt)
        self.assertTrue(len(orbit) == 0)

        orbit = {(0, 1, 0), (1, -1, 0), (-1, 0, 1), (0, 0, -1)}
        for wt in liealg.get_orbit_iter((0, 1, 0)):
            self.assertTrue(wt in orbit)
            orbit.remove(wt)
        self.assertTrue(len(orbit) == 0)

        orbit = {(0, 0, 1), (1, 0, -1), (-1, 1, 0), (0, -1, 0)}
        for wt in liealg.get_orbit_iter((0, 0, 1)):
            self.assertTrue(wt in orbit)
            orbit.remove(wt)
        self.assertTrue(len(orbit) == 0)

        orbit = {(1,1,0),(-1,2,1),(2,-1,0),(1,-2,1),(0,2,-1),(-2,1,2),(-1,-1,2),(2,-2,-1),(0,1,-2),(1,-1,-2),(-2,0,1),(-1,0,-1)}
        for wt in liealg.get_orbit_iter((1, 1, 0)):
            self.assertTrue(wt in orbit, msg=str(wt))
            orbit.remove(wt)
        self.assertTrue(len(orbit) == 0)

        orbit = {(1,0,1),(-1,1,2),(2,0,-1),(0,-1,2),(1,1,-2),(-2,2,1),(2,-1,-2),(-1,2,-1),(0,-2,1),(-2,1,0),(1,-2,-1),(-1,-1,0)}
        for wt in liealg.get_orbit_iter((1, 0, 1)):
            self.assertTrue(wt in orbit, msg=str(wt))
            orbit.remove(wt)
        self.assertTrue(len(orbit) == 0)

        orbit = {(0,1,1),(1,-1,1),(1,1,-1),(-1,0,2),(2,-1,-1),(-1,2,0),(1,0,-2),(-2,1,1),(1,-2,0),(-1,1,-1),(-1,-1,1),(0,-1,-1)}
        for wt in liealg.get_orbit_iter((0, 1, 1)):
            self.assertTrue(wt in orbit, msg=str(wt))
            orbit.remove(wt)
        self.assertTrue(len(orbit) == 0)

        orbit = {(1,1,1),(-1,2,2),(2,-1,1),(2,1,-1),(1,-2,2),(1,2,-2),(-2,1,3),(3,-1,-1),(-2,3,1),(-1,-1,3),(3,-2,-2),(-1,3,-1),(1,1,-3),(-3,2,2),(1,-3,1),(2,-1,-3),(-3,1,1),(2,-3,-1),(-1,2,-2),(-1,-2,2),(-2,1,-1),(-2,-1,1),(1,-2,-2),(-1,-1,-1)}
        for wt in liealg.get_orbit_iter((1, 1, 1)):
            self.assertTrue(wt in orbit, msg=str(wt))
            orbit.remove(wt)
        self.assertTrue(len(orbit) == 0)

        liealg = cbd.TypeDLieAlgebra(4, exact=False)

        orbit = {(0, 0, 0, 0)}
        for wt in liealg.get_orbit_iter((0, 0, 0, 0)):
            self.assertTrue(wt in orbit, msg=str(wt))
            orbit.remove(wt)
        self.assertTrue(len(orbit) == 0)

        orbit = {(1,0,0,0),(-1,1,0,0),(0,-1,1,1),(0,0,-1,1),(0,0,1,-1),(0,1,-1,-1),(1,-1,0,0),(-1,0,0,0)}
        for wt in liealg.get_orbit_iter((1, 0, 0, 0)):
            self.assertTrue(wt in orbit, msg=str(wt))
            orbit.remove(wt)
        self.assertTrue(len(orbit) == 0)

        orbit = {(0,1,0,0),(1,-1,1,1),(-1,0,1,1),(1,0,-1,1),(1,0,1,-1),(-1,1,-1,1),(-1,1,1,-1),(1,1,-1,-1),(0,-1,0,2),(-1,2,-1,-1),(0,-1,2,0),(2,-1,0,0),(0,1,0,-2),(1,-2,1,1),(0,1,-2,0),(-2,1,0,0),(1,-1,1,-1),(-1,-1,1,1),(1,-1,-1,1),(-1,0,1,-1),(1,0,-1,-1),(-1,0,-1,1),(-1,1,-1,-1),(0,-1,0,0)}
        for wt in liealg.get_orbit_iter((0, 1, 0, 0)):
            self.assertTrue(wt in orbit, msg=str(wt))
            orbit.remove(wt)
        self.assertTrue(len(orbit) == 0)

        orbit = {(0,0,1,0),(0,1,-1,0),(1,-1,0,1),(-1,0,0,1),(1,0,0,-1),(-1,1,0,-1),(0,-1,1,0),(0,0,-1,0)}
        for wt in liealg.get_orbit_iter((0, 0, 1, 0)):
            self.assertTrue(wt in orbit, msg=str(wt))
            orbit.remove(wt)
        self.assertTrue(len(orbit) == 0)

        orbit = {(0,0,0,1),(0,1,0,-1),(1,-1,1,0),(-1,0,1,0),(1,0,-1,0),(-1,1,-1,0),(0,-1,0,1),(0,0,0,-1)}
        for wt in liealg.get_orbit_iter((0, 0, 0, 1)):
            self.assertTrue(wt in orbit, msg=str(wt))
            orbit.remove(wt)
        self.assertTrue(len(orbit) == 0)

        orbit = {(1,1,0,0),(-1,2,0,0),(2,-1,1,1),(1,-2,2,2),(-2,1,1,1),(2,0,-1,1),(2,0,1,-1),(-1,-1,2,2),(1,0,-2,2),(1,0,2,-2),(-2,2,-1,1),(-2,2,1,-1),(2,1,-1,-1),(-1,1,-2,2),(-1,1,2,-2),(1,2,-2,-2),(0,-2,1,3),(-2,3,-1,-1),(0,-2,3,1),(3,-1,0,0),(0,-1,-1,3),(-1,3,-2,-2),(0,-1,3,-1),(3,-2,0,0),(0,1,1,-3),(1,-3,2,2),(0,1,-3,1),(-3,2,0,0),(0,2,-1,-3),(2,-3,1,1),(0,2,-3,-1),(-3,1,0,0),(1,-1,2,-2),(-1,-2,2,2),(1,-1,-2,2),(2,-2,1,-1),(-2,-1,1,1),(2,-2,-1,1),(-1,0,2,-2),(1,1,-2,-2),(-1,0,-2,2),(-2,0,1,-1),(2,-1,-1,-1),(-2,0,-1,1),(-1,2,-2,-2),(-2,1,-1,-1),(1,-2,0,0),(-1,-1,0,0)}
        for wt in liealg.get_orbit_iter((1, 1, 0, 0)):
            self.assertTrue(wt in orbit, msg=str(wt))
            orbit.remove(wt)
        self.assertTrue(len(orbit) == 0)

        orbit = {(1,0,1,0),(-1,1,1,0),(1,1,-1,0),(0,-1,2,1),(-1,2,-1,0),(2,-1,0,1),(0,1,-2,1),(0,0,2,-1),(1,-2,1,2),(-2,1,0,1),(2,0,0,-1),(1,-1,-1,2),(0,2,-2,-1),(-1,-1,1,2),(1,0,1,-2),(-2,2,0,-1),(-1,0,-1,2),(1,1,-1,-2),(2,-2,0,1),(-1,1,1,-2),(0,-2,2,1),(-1,2,-1,-2),(2,-1,0,-1),(-2,0,0,1),(0,-1,2,-1),(0,0,-2,1),(1,-2,1,0),(-2,1,0,-1),(0,1,-2,-1),(-1,-1,1,0),(1,-1,-1,0),(-1,0,-1,0)}
        for wt in liealg.get_orbit_iter((1, 0, 1, 0)):
            self.assertTrue(wt in orbit, msg=str(wt))
            orbit.remove(wt)
        self.assertTrue(len(orbit) == 0)

        orbit = {(1,0,0,1),(-1,1,0,1),(1,1,0,-1),(0,-1,1,2),(-1,2,0,-1),(2,-1,1,0),(0,0,-1,2),(0,1,1,-2),(1,-2,2,1),(-2,1,1,0),(2,0,-1,0),(0,2,-1,-2),(1,-1,2,-1),(-1,-1,2,1),(1,0,-2,1),(-2,2,-1,0),(2,-2,1,0),(-1,0,2,-1),(1,1,-2,-1),(-1,1,-2,1),(0,-2,1,2),(-2,0,1,0),(2,-1,-1,0),(-1,2,-2,-1),(0,-1,-1,2),(0,0,1,-2),(-2,1,-1,0),(1,-2,0,1),(0,1,-1,-2),(-1,-1,0,1),(1,-1,0,-1),(-1,0,0,-1)}
        for wt in liealg.get_orbit_iter((1, 0, 0, 1)):
            self.assertTrue(wt in orbit, msg=str(wt))
            orbit.remove(wt)
        self.assertTrue(len(orbit) == 0)

        orbit = {(0,1,1,0),(1,-1,2,1),(0,2,-1,0),(-1,0,2,1),(1,1,-2,1),(1,0,2,-1),(2,-2,1,2),(-1,2,-2,1),(-1,1,2,-1),(2,-1,-1,2),(1,2,-2,-1),(-2,0,1,2),(2,0,1,-2),(1,-2,0,3),(-1,3,-2,-1),(0,-1,3,0),(-2,1,-1,2),(2,1,-1,-2),(3,-2,0,1),(-2,2,1,-2),(-1,-1,0,3),(1,1,0,-3),(2,-3,1,2),(0,2,-3,0),(-2,3,-1,-2),(3,-1,0,-1),(-3,1,0,1),(0,-2,3,0),(-1,2,0,-3),(2,-1,1,-2),(-2,-1,1,2),(2,-2,-1,2),(1,-3,2,1),(-3,2,0,-1),(0,1,-3,0),(1,-2,2,-1),(-2,1,1,-2),(2,0,-1,-2),(-2,0,-1,2),(-1,-2,2,1),(1,-1,-2,1),(-1,-1,2,-1),(1,0,-2,-1),(-2,2,-1,-2),(-1,0,-2,1),(-1,1,-2,-1),(0,-2,1,0),(0,-1,-1,0)}
        for wt in liealg.get_orbit_iter((0, 1, 1, 0)):
            self.assertTrue(wt in orbit, msg=str(wt))
            orbit.remove(wt)
        self.assertTrue(len(orbit) == 0)

        orbit = {(0,1,0,1),(1,-1,1,2),(0,2,0,-1),(-1,0,1,2),(1,0,-1,2),(1,1,1,-2),(2,-2,2,1),(-1,1,-1,2),(-1,2,1,-2),(1,2,-1,-2),(2,-1,2,-1),(-2,0,2,1),(2,0,-2,1),(0,-1,0,3),(-1,3,-1,-2),(1,-2,3,0),(3,-2,1,0),(-2,1,2,-1),(2,1,-2,-1),(-2,2,-2,1),(0,2,0,-3),(2,-3,2,1),(-1,-1,3,0),(1,1,-3,0),(-3,1,1,0),(3,-1,-1,0),(-2,3,-2,-1),(0,-2,0,3),(2,-2,2,-1),(-2,-1,2,1),(2,-1,-2,1),(-1,2,-3,0),(-3,2,-1,0),(1,-3,1,2),(0,1,0,-3),(-2,0,2,-1),(2,0,-2,-1),(-2,1,-2,1),(1,-2,-1,2),(-1,-2,1,2),(1,-1,1,-2),(-2,2,-2,-1),(-1,-1,-1,2),(1,0,-1,-2),(-1,0,1,-2),(0,-2,0,1),(-1,1,-1,-2),(0,-1,0,-1)}
        for wt in liealg.get_orbit_iter((0, 1, 0, 1)):
            self.assertTrue(wt in orbit, msg=str(wt))
            orbit.remove(wt)
        self.assertTrue(len(orbit) == 0)

        orbit = {(0,0,1,1),(0,1,-1,1),(0,1,1,-1),(1,-1,0,2),(0,2,-1,-1),(1,-1,2,0),(-1,0,0,2),(1,1,0,-2),(2,-2,1,1),(-1,0,2,0),(1,1,-2,0),(-1,2,0,-2),(2,-1,1,-1),(-2,0,1,1),(2,-1,-1,1),(-1,2,-2,0),(1,-2,2,0),(-2,1,1,-1),(2,0,-1,-1),(-2,1,-1,1),(1,-2,0,2),(-1,-1,2,0),(1,0,-2,0),(-2,2,-1,-1),(-1,-1,0,2),(1,0,0,-2),(-1,1,-2,0),(0,-2,1,1),(-1,1,0,-2),(0,-1,-1,1),(0,-1,1,-1),(0,0,-1,-1)}
        for wt in liealg.get_orbit_iter((0, 0, 1, 1)):
            self.assertTrue(wt in orbit, msg=str(wt))
            orbit.remove(wt)
        self.assertTrue(len(orbit) == 0)

        orbit = {(1,1,1,0),(-1,2,1,0),(2,-1,2,1),(1,2,-1,0),(1,-2,3,2),(-1,3,-1,0),(-2,1,2,1),(2,1,-2,1),(2,0,2,-1),(3,-2,1,2),(-1,-1,3,2),(1,1,-3,2),(1,0,3,-2),(2,-3,2,3),(-2,3,-2,1),(-2,2,2,-1),(3,-1,-1,2),(2,2,-2,-1),(-3,1,1,2),(3,0,1,-2),(-1,2,-3,2),(-1,1,3,-2),(2,-1,-2,3),(1,3,-3,-2),(-2,-1,2,3),(2,0,2,-3),(1,-3,1,4),(-2,4,-2,-1),(0,-2,4,1),(-3,2,-1,2),(3,1,-1,-2),(4,-2,0,1),(-3,3,1,-2),(1,-2,-1,4),(-1,4,-3,-2),(0,-1,4,-1),(-2,1,-2,3),(2,2,-2,-3),(4,-3,0,1),(-2,2,2,-3),(-1,-2,1,4),(1,1,1,-4),(2,-4,2,3),(0,2,-4,1),(-3,4,-1,-2),(4,-1,0,-1),(-4,2,0,1),(0,-3,4,1),(-1,-1,-1,4),(1,2,-1,-4),(3,-4,1,2),(0,3,-4,-1),(-2,4,-2,-3),(4,-2,0,-1),(-4,1,0,1),(0,-2,4,-1),(-1,2,1,-4),(2,-1,2,-3),(-2,-2,2,3),(2,-2,-2,3),(1,-4,3,2),(-4,3,0,-1),(0,1,-4,1),(-1,3,-1,-4),(3,-2,1,-2),(-3,-1,1,2),(3,-3,-1,2),(2,-4,2,1),(-4,2,0,-1),(0,2,-4,-1),(1,-2,3,-2),(-2,1,2,-3),(2,1,-2,-3),(-2,0,-2,3),(-1,-3,3,2),(1,-1,-3,2),(2,-3,2,-1),(-3,1,1,-2),(3,-1,-1,-2),(-3,0,-1,2),(-2,-2,2,1),(2,-2,-2,1),(-1,-1,3,-2),(1,1,-3,-2),(-2,3,-2,-3),(-1,0,-3,2),(-2,-1,2,-1),(2,-1,-2,-1),(-3,2,-1,-2),(-2,0,-2,1),(-1,2,-3,-2),(1,-3,1,0),(-2,1,-2,-1),(-1,-2,1,0),(1,-2,-1,0),(-1,-1,-1,0)}
        for wt in liealg.get_orbit_iter((1, 1, 1, 0)):
            self.assertTrue(wt in orbit, msg=str(wt))
            orbit.remove(wt)
        self.assertTrue(len(orbit) == 0)

        orbit = {(1,1,0,1),(-1,2,0,1),(2,-1,1,2),(1,2,0,-1),(1,-2,2,3),(-1,3,0,-1),(-2,1,1,2),(2,0,-1,2),(2,1,1,-2),(3,-2,2,1),(-1,-1,2,3),(1,0,-2,3),(1,1,2,-3),(2,-3,3,2),(-2,2,-1,2),(-2,3,1,-2),(2,2,-1,-2),(3,-1,2,-1),(-3,1,2,1),(3,0,-2,1),(-1,1,-2,3),(-1,2,2,-3),(1,3,-2,-3),(2,-1,3,-2),(-2,-1,3,2),(2,0,-3,2),(0,-2,1,4),(-2,4,-1,-2),(1,-3,4,1),(4,-2,1,0),(-3,2,2,-1),(3,1,-2,-1),(-3,3,-2,1),(0,-1,-1,4),(-1,4,-2,-3),(1,-2,4,-1),(4,-3,1,0),(-2,1,3,-2),(2,2,-3,-2),(-2,2,-3,2),(0,2,1,-4),(2,-4,3,2),(-1,-2,4,1),(1,1,-4,1),(-4,2,1,0),(4,-1,-1,0),(-3,4,-2,-1),(0,-3,1,4),(0,3,-1,-4),(3,-4,2,1),(-1,-1,4,-1),(1,2,-4,-1),(-4,1,1,0),(4,-2,-1,0),(-2,4,-3,-2),(0,-2,-1,4),(2,-2,3,-2),(-2,-2,3,2),(2,-1,-3,2),(-1,2,-4,1),(-4,3,-1,0),(1,-4,2,3),(0,1,1,-4),(3,-3,2,-1),(-3,-1,2,1),(3,-2,-2,1),(-1,3,-4,-1),(-4,2,-1,0),(2,-4,1,2),(0,2,-1,-4),(-2,0,3,-2),(2,1,-3,-2),(-2,1,-3,2),(1,-2,-2,3),(-1,-3,2,3),(1,-1,2,-3),(-3,0,2,-1),(3,-1,-2,-1),(-3,1,-2,1),(2,-3,-1,2),(-2,-2,1,2),(2,-2,1,-2),(-2,3,-3,-2),(-1,-1,-2,3),(1,1,-2,-3),(-1,0,2,-3),(-3,2,-2,-1),(-2,-1,-1,2),(2,-1,-1,-2),(-2,0,1,-2),(1,-3,0,1),(-1,2,-2,-3),(-1,-2,0,1),(-2,1,-1,-2),(1,-2,0,-1),(-1,-1,0,-1)}
        for wt in liealg.get_orbit_iter((1, 1, 0, 1)):
            self.assertTrue(wt in orbit, msg=str(wt))
            orbit.remove(wt)
        self.assertTrue(len(orbit) == 0)

        orbit = {(1,0,1,1),(-1,1,1,1),(1,1,-1,1),(1,1,1,-1),(0,-1,2,2),(-1,2,-1,1),(-1,2,1,-1),(2,-1,0,2),(1,2,-1,-1),(2,-1,2,0),(0,1,-2,2),(0,1,2,-2),(1,-2,1,3),(-1,3,-1,-1),(1,-2,3,1),(-2,1,0,2),(2,1,0,-2),(3,-2,1,1),(-2,1,2,0),(2,1,-2,0),(1,-1,-1,3),(0,3,-2,-2),(1,-1,3,-1),(-1,-1,1,3),(1,1,1,-3),(2,-3,2,2),(-1,-1,3,1),(1,1,-3,1),(-2,3,0,-2),(3,-1,1,-1),(-3,1,1,1),(3,-1,-1,1),(-2,3,-2,0),(-1,0,-1,3),(1,2,-1,-3),(3,-3,1,1),(-1,0,3,-1),(1,2,-3,-1),(-1,2,1,-3),(2,-1,2,-2),(-2,-1,2,2),(2,-1,-2,2),(-1,2,-3,1),(1,-3,3,1),(-3,2,1,-1),(3,0,-1,-1),(-3,2,-1,1),(1,-3,1,3),(-1,3,-1,-3),(3,-2,1,-1),(-3,0,1,1),(3,-2,-1,1),(-1,3,-3,-1),(1,-2,3,-1),(-2,1,2,-2),(2,1,-2,-2),(-2,1,-2,2),(1,-2,-1,3),(-1,-2,3,1),(1,0,-3,1),(-3,3,-1,-1),(-1,-2,1,3),(1,0,1,-3),(2,-3,2,0),(-3,1,1,-1),(3,-1,-1,-1),(-3,1,-1,1),(2,-3,0,2),(-1,-1,3,-1),(1,1,-3,-1),(-2,3,-2,-2),(-1,-1,-1,3),(1,1,-1,-3),(-1,1,-3,1),(0,-3,2,2),(-1,1,1,-3),(-2,-1,2,0),(2,-1,-2,0),(-3,2,-1,-1),(-2,-1,0,2),(2,-1,0,-2),(-1,2,-3,-1),(1,-3,1,1),(-1,2,-1,-3),(0,-1,-2,2),(0,-1,2,-2),(-2,1,-2,0),(-1,-2,1,1),(-2,1,0,-2),(1,-2,-1,1),(1,-2,1,-1),(0,1,-2,-2),(-1,-1,-1,1),(-1,-1,1,-1),(1,-1,-1,-1),(-1,0,-1,-1)}
        for wt in liealg.get_orbit_iter((1, 0, 1, 1)):
            self.assertTrue(wt in orbit, msg=str(wt))
            orbit.remove(wt)
        self.assertTrue(len(orbit) == 0)

        orbit = {(1,1,1,1),(-1,2,1,1),(2,-1,2,2),(1,2,-1,1),(1,2,1,-1),(1,-2,3,3),(-1,3,-1,1),(-1,3,1,-1),(-2,1,2,2),(2,1,-2,2),(2,1,2,-2),(3,-2,1,3),(1,3,-1,-1),(3,-2,3,1),(-1,-1,3,3),(1,1,-3,3),(1,1,3,-3),(2,-3,2,4),(-1,4,-1,-1),(2,-3,4,2),(-2,3,-2,2),(-2,3,2,-2),(3,-1,-1,3),(2,3,-2,-2),(3,-1,3,-1),(-3,1,1,3),(3,1,1,-3),(4,-3,2,2),(-3,1,3,1),(3,1,-3,1),(-1,2,-3,3),(-1,2,3,-3),(2,-1,-2,4),(1,4,-3,-3),(2,-1,4,-2),(-2,-1,2,4),(2,1,2,-4),(3,-4,3,3),(-2,-1,4,2),(2,1,-4,2),(1,-3,1,5),(-2,5,-2,-2),(1,-3,5,1),(-3,2,-1,3),(3,2,-1,-3),(5,-3,1,1),(-3,2,3,-1),(3,2,-3,-1),(-3,4,1,-3),(4,-1,2,-2),(-4,1,2,2),(4,-1,-2,2),(-3,4,-3,1),(1,-2,-1,5),(-1,5,-3,-3),(1,-2,5,-1),(-2,1,-2,4),(2,3,-2,-4),(5,-4,1,1),(-2,1,4,-2),(2,3,-4,-2),(-2,3,2,-4),(3,-1,3,-3),(-3,-1,3,3),(3,-1,-3,3),(-2,3,-4,2),(-1,-2,1,5),(1,2,1,-5),(3,-5,3,3),(-1,-2,5,1),(1,2,-5,1),(-3,5,-1,-3),(5,-2,1,-1),(-5,2,1,1),(5,-2,-1,1),(-3,5,-3,-1),(1,-4,5,1),(-4,3,2,-2),(4,1,-2,-2),(-4,3,-2,2),(1,-4,1,5),(-1,-1,-1,5),(1,3,-1,-5),(4,-5,2,2),(-1,-1,5,-1),(1,3,-5,-1),(-2,5,-2,-4),(5,-3,1,-1),(-5,1,1,1),(5,-3,-1,1),(-2,5,-4,-2),(1,-3,5,-1),(-3,2,3,-3),(3,2,-3,-3),(-3,2,-3,3),(1,-3,-1,5),(-1,3,1,-5),(3,-2,3,-3),(-3,-2,3,3),(3,-2,-3,3),(-1,3,-5,1),(2,-5,4,2),(-5,3,1,-1),(5,-1,-1,-1),(-5,3,-1,1),(2,-5,2,4),(-1,-3,5,1),(1,1,-5,1),(-4,5,-2,-2),(-1,-3,1,5),(1,1,1,-5),(-1,4,-1,-5),(4,-3,2,-2),(-4,-1,2,2),(4,-3,-2,2),(-1,4,-5,-1),(3,-5,3,1),(-5,2,1,-1),(5,-2,-1,-1),(-5,2,-1,1),(3,-5,1,3),(-1,-2,5,-1),(1,2,-5,-1),(-3,5,-3,-3),(-1,-2,-1,5),(1,2,-1,-5),(2,-3,4,-2),(-3,1,3,-3),(3,1,-3,-3),(-3,1,-3,3),(2,-3,-2,4),(-2,-3,4,2),(2,-1,-4,2),(-5,4,-1,-1),(-2,-3,2,4),(2,-1,2,-4),(-1,2,-5,1),(1,-5,3,3),(-1,2,1,-5),(3,-4,3,-1),(-4,1,2,-2),(4,-1,-2,-2),(-4,1,-2,2),(3,-4,-1,3),(-3,-2,3,1),(3,-2,-3,1),(-5,3,-1,-1),(-3,-2,1,3),(3,-2,1,-3),(-1,3,-5,-1),(2,-5,2,2),(-1,3,-1,-5),(-2,-1,4,-2),(2,1,-4,-2),(-3,4,-3,-3),(-2,-1,-2,4),(2,1,-2,-4),(-2,1,-4,2),(-1,-4,3,3),(-2,1,2,-4),(1,-2,-3,3),(1,-2,3,-3),(-3,-1,3,-1),(3,-1,-3,-1),(-4,3,-2,-2),(-3,-1,-1,3),(3,-1,-1,-3),(-3,1,-3,1),(-2,-3,2,2),(-3,1,1,-3),(2,-3,-2,2),(2,-3,2,-2),(-2,3,-4,-2),(1,-4,1,1),(-2,3,-2,-4),(-1,-1,-3,3),(-1,-1,3,-3),(1,1,-3,-3),(-3,2,-3,-1),(-1,-3,1,1),(-3,2,-1,-3),(-2,-1,-2,2),(-2,-1,2,-2),(2,-1,-2,-2),(1,-3,-1,1),(1,-3,1,-1),(-1,2,-3,-3),(-1,-2,-1,1),(-1,-2,1,-1),(-2,1,-2,-2),(1,-2,-1,-1),(-1,-1,-1,-1)}
        for wt in liealg.get_orbit_iter((1, 1, 1, 1)):
            self.assertTrue(wt in orbit, msg=str(wt))
            orbit.remove(wt)
        self.assertTrue(len(orbit) == 0)


    def test_tensor(self):
        liealg = cbd.TypeDLieAlgebra(3, exact=False)

        liealg = cbd.TypeDLieAlgebra(4, exact=False)

        product = {(0,0,0,0): 1 ,(0,1,0,0): 1 ,(2,0,0,0): 1}
        for wt, mult in liealg.tensor((1,0,0,0), (1,0,0,0)).iteritems():
            self.assertTrue(wt in product, msg="Weight not in product: " + str(wt))
            self.assertTrue(product[wt] == mult, msg="Weight has wrong multiplicity: " + str(wt))
            product.pop(wt)
        self.assertTrue(len(product) == 0)

        product = {(0,0,1,2): 2,(0,0,3,0): 1,(0,0,3,2): 1,(0,1,1,0): 2,(0,1,1,2): 2,(0,1,3,0): 1,(0,2,1,0): 2,(1,0,0,1): 1,(1,0,0,3): 1,(1,0,2,1): 3,(1,1,0,1): 3,(1,1,2,1): 1,(1,2,0,1): 1,(2,0,1,0): 1,(2,0,1,2): 1,(2,1,1,0): 1}
        for wt, mult in liealg.tensor((1,0,1,0), (0,1,1,1)).iteritems():
            self.assertTrue(wt in product, msg="Weight not in product: " + str(wt))
            self.assertTrue(product[wt] == mult, msg="Weight has wrong multiplicity: " + str(wt))
            product.pop(wt)
        self.assertTrue(len(product) == 0)


if __name__ == '__main__':
    unittest.main()
