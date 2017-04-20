'''
Created on Nov 2, 2016

@author: mjschust
'''
import unittest
import fusion_prod.cbbundle as cbd

class Test(unittest.TestCase):


    def testSL4Weights(self):
        
        liealg = cbd.TypeALieAlgebra(3)
        wt = tuple([0, 2, 1])
        self.assertItemsEqual([3,3,1,0], liealg._convert_funds_to_epsilons(wt), "Epsilon coords incorrect")
        rt = cbd._Root(liealg, [1, 1, 0])
        self.assertItemsEqual([1,1,-1], rt, "Fundamental coords incorrect")
        self.assertItemsEqual([1,0,-1,0], liealg._convert_funds_to_epsilons(rt), "Epsilon coords incorrect")
        self.assertEqual(2, liealg.killing_form(wt, rt), "Killing product incorrect")
        self.assertItemsEqual([1,0,1], liealg.reflect_to_chamber(rt), "Reflection incorrect")
        rt_list = liealg.get_positive_roots()

        
    def testSL2Weights(self):
        
        liealg = cbd.TypeALieAlgebra(1)
        wt = tuple([1])
        self.assertItemsEqual([1,0], liealg._convert_funds_to_epsilons(wt), "Epsilon coords incorrect")
        rt = cbd._Root(liealg, [1])
        self.assertItemsEqual([2], rt, "Fundamental coords incorrect")
        self.assertItemsEqual([2,0], liealg._convert_funds_to_epsilons(rt), "Epsilon coords incorrect")
        self.assertEqual(1, liealg.killing_form(wt, rt), "Killing product incorrect")
        wt2 = tuple([-1])
        self.assertItemsEqual([1], liealg.reflect_to_chamber(wt), "Reflection incorrect")
        rt_list = liealg.get_positive_roots()

        
        

        

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testSL4Weights']
    unittest.main()