'''
Created on Nov 19, 2016

@author: mjschust
'''
from __future__ import division
import unittest
import fusion_prod.cbd as cbd

class Test(unittest.TestCase):
    
    def testSL2OrbitIter(self):
        liealg = cbd.TypeALieAlgebra(1)
        wt = cbd.Weight(liealg, [0])
        orbit_wt = wt
        orbit_dict = {orbit_wt}
        for wt2 in liealg.get_orbit_iter(wt):
            self.assertTrue(wt2 in orbit_dict, "Orbit incorrect")
            
        wt = cbd.Weight(liealg, [-1])
        orbit_wt1 = cbd.Weight(liealg, [-1])
        orbit_wt2 = cbd.Weight(liealg, [1])
        orbit_dict = {orbit_wt1, orbit_wt2}
        for wt2 in liealg.get_orbit_iter(wt):
            self.assertTrue(wt2 in orbit_dict, "Orbit incorrect")
    
    def testSL2Tensor(self):
        liealg = cbd.TypeALieAlgebra(1)
        decomp = liealg.tensor([0], [0])
        dec_wt = cbd.Weight(liealg, [0])
        self.assertTrue(dec_wt in decomp, "Tensor decomp incorrect")
        self.assertEqual(1, decomp[dec_wt], "Tensor decomp incorrect")

        decomp = liealg.tensor([0], [1])
        dec_wt = cbd.Weight(liealg, [1])
        self.assertTrue(dec_wt in decomp, "Tensor decomp incorrect")
        self.assertEqual(1, decomp[dec_wt], "Tensor decomp incorrect")

        decomp = liealg.tensor([1], [1])
        dec_wt = cbd.Weight(liealg, [2])
        self.assertTrue(dec_wt in decomp, "Tensor decomp incorrect")
        self.assertEqual(1, decomp[dec_wt], "Tensor decomp incorrect")
        dec_wt = cbd.Weight(liealg, [0])
        self.assertTrue(dec_wt in decomp, "Tensor decomp incorrect")
        self.assertEqual(1, decomp[dec_wt], "Tensor decomp incorrect")

        decomp = liealg.tensor([2], [1])
        dec_wt = cbd.Weight(liealg, [3])
        self.assertTrue(dec_wt in decomp, "Tensor decomp incorrect")
        self.assertEqual(1, decomp[dec_wt], "Tensor decomp incorrect")
        dec_wt = cbd.Weight(liealg, [1])
        self.assertTrue(dec_wt in decomp, "Tensor decomp incorrect")
        self.assertEqual(1, decomp[dec_wt], "Tensor decomp incorrect")

        decomp = liealg.tensor([5], [2])
        dec_wt = cbd.Weight(liealg, [7])
        self.assertTrue(dec_wt in decomp, "Tensor decomp incorrect")
        self.assertEqual(1, decomp[dec_wt], "Tensor decomp incorrect")
        dec_wt = cbd.Weight(liealg, [5])
        self.assertTrue(dec_wt in decomp, "Tensor decomp incorrect")
        self.assertEqual(1, decomp[dec_wt], "Tensor decomp incorrect")
        dec_wt = cbd.Weight(liealg, [3])
        self.assertTrue(dec_wt in decomp, "Tensor decomp incorrect")
        self.assertEqual(1, decomp[dec_wt], "Tensor decomp incorrect")
    
    def testSL2Fusion(self):
        liealg = cbd.TypeALieAlgebra(1)
        decomp = liealg.fusion([0], [0],1)
        dec_wt = cbd.Weight(liealg, [0])
        self.assertTrue(dec_wt in decomp, "Tensor decomp incorrect")
        self.assertEqual(1, decomp[dec_wt], "Tensor decomp incorrect")

        decomp = liealg.fusion([0], [1], 1)
        dec_wt = cbd.Weight(liealg, [1])
        self.assertTrue(dec_wt in decomp, "Tensor decomp incorrect")
        self.assertEqual(1, decomp[dec_wt], "Tensor decomp incorrect")

        decomp = liealg.fusion([1], [1],1)
        dec_wt = cbd.Weight(liealg, [2])
        self.assertFalse(dec_wt in decomp, "Tensor decomp incorrect")
        dec_wt = cbd.Weight(liealg, [0])
        self.assertTrue(dec_wt in decomp, "Tensor decomp incorrect")
        self.assertEqual(1, decomp[dec_wt], "Tensor decomp incorrect")
        
        decomp = liealg.fusion([1], [1],2)
        dec_wt = cbd.Weight(liealg, [2])
        self.assertTrue(dec_wt in decomp, "Tensor decomp incorrect")
        self.assertEqual(1, decomp[dec_wt], "Tensor decomp incorrect")
        dec_wt = cbd.Weight(liealg, [0])
        self.assertTrue(dec_wt in decomp, "Tensor decomp incorrect")
        self.assertEqual(1, decomp[dec_wt], "Tensor decomp incorrect")

        decomp = liealg.fusion([5], [2],7)
        dec_wt = cbd.Weight(liealg, [7])
        self.assertTrue(dec_wt in decomp, "Tensor decomp incorrect")
        self.assertEqual(1, decomp[dec_wt], "Tensor decomp incorrect")
        dec_wt = cbd.Weight(liealg, [5])
        self.assertTrue(dec_wt in decomp, "Tensor decomp incorrect")
        self.assertEqual(1, decomp[dec_wt], "Tensor decomp incorrect")
        dec_wt = cbd.Weight(liealg, [3])
        self.assertTrue(dec_wt in decomp, "Tensor decomp incorrect")
        self.assertEqual(1, decomp[dec_wt], "Tensor decomp incorrect")
        
        decomp = liealg.fusion([5], [2],6)
        dec_wt = cbd.Weight(liealg, [7])
        self.assertFalse(dec_wt in decomp, "Tensor decomp incorrect")
        dec_wt = cbd.Weight(liealg, [5])
        self.assertTrue(dec_wt in decomp, "Tensor decomp incorrect")
        self.assertEqual(1, decomp[dec_wt], "Tensor decomp incorrect")
        dec_wt = cbd.Weight(liealg, [3])
        self.assertTrue(dec_wt in decomp, "Tensor decomp incorrect")
        self.assertEqual(1, decomp[dec_wt], "Tensor decomp incorrect")
        
        decomp = liealg.fusion([5], [2],5)
        dec_wt = cbd.Weight(liealg, [7])
        self.assertFalse(dec_wt in decomp, "Tensor decomp incorrect")
        dec_wt = cbd.Weight(liealg, [5])
        self.assertFalse(dec_wt in decomp and decomp[dec_wt] > 0, "Tensor decomp incorrect")
        dec_wt = cbd.Weight(liealg, [3])
        self.assertTrue(dec_wt in decomp, "Tensor decomp incorrect")
        self.assertEqual(1, decomp[dec_wt], "Tensor decomp incorrect")

    def testSL2MultiFusion(self):
        liealg = cbd.TypeALieAlgebra(1)
        wt1 = cbd.Weight(liealg, [0])
        wt2 = cbd.Weight(liealg, [1])
        decomp = liealg.multi_fusion([wt1, wt2, wt2], 1)
        dec_wt = cbd.Weight(liealg, [2])
        self.assertFalse(dec_wt in decomp, "Tensor decomp incorrect")
        dec_wt = cbd.Weight(liealg, [0])
        self.assertTrue(dec_wt in decomp, "Tensor decomp incorrect")
        self.assertEqual(1, decomp[dec_wt], "Tensor decomp incorrect")

        decomp = liealg.multi_fusion([wt2, wt2, wt2], 1)
        dec_wt = cbd.Weight(liealg, [1])
        self.assertTrue(dec_wt in decomp, "Tensor decomp incorrect")
        self.assertEqual(1, decomp[dec_wt], "Tensor decomp incorrect")

        wt1 = cbd.Weight(liealg, [4])
        wt2 = cbd.Weight(liealg, [2])
        decomp = liealg.multi_fusion([wt2, wt2, wt2], 4)
        dec_wt = cbd.Weight(liealg, [4])
        self.assertTrue(dec_wt in decomp, "Tensor decomp incorrect")
        self.assertEqual(1, decomp[dec_wt], "Tensor decomp incorrect")
        dec_wt = cbd.Weight(liealg, [2])
        self.assertTrue(dec_wt in decomp, "Tensor decomp incorrect")
        self.assertEqual(3, decomp[dec_wt], "Tensor decomp incorrect")
        dec_wt = cbd.Weight(liealg, [0])
        self.assertTrue(dec_wt in decomp, "Tensor decomp incorrect")
        self.assertEqual(1, decomp[dec_wt], "Tensor decomp incorrect")

    def testDegree(self):
        liealg = cbd.TypeALieAlgebra(1)
        wt1 = cbd.Weight(liealg, [1])
        wt2 = cbd.Weight(liealg, [3])
        wt3 = cbd.Weight(liealg, [5])
        #self.assertEqual(1, liealg.degree(wt1,wt2,wt2,wt3, 5), "Degree incorrect")
        #self.assertEqual(0, liealg.degree(wt1, wt2, wt2, wt3, 6), "Degree incorrect")

        liealg = cbd.TypeALieAlgebra(4)
        wt1 = cbd.Weight(liealg, [0,1,0,0])
        wt2 = cbd.Weight(liealg, [0,0,1,0])
        self.assertEqual(2, liealg.degree(wt1, wt1, wt2, wt2, 1), "Degree incorrect")




if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testSL2OrbitIter','Test.testSL2Tensor','Test.testSL2Fusion','testSL2MultiFusion','testDegree']
    unittest.main()