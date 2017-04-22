'''
Created on Nov 2, 2016

@author: mjschust
'''
import unittest
import conformal_blocks.cbbundle as cbd

class Test(unittest.TestCase):


    def test_SL4Weights(self):
        
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
    #import sys;sys.argv = ['', 'Test.test_SL4Weights']
    unittest.main()