'''
Created on Nov 19, 2016

@author: mjschust
'''
import unittest
import fusion_prod.cbd as cbd

class Test(unittest.TestCase):


    def testConformalBlocksRanks(self):
        liealg = cbd.TypeALieAlgebra(3)
        wt1 = cbd._Weight(liealg, [0, 1, 4])
        wt2 = cbd._Weight(liealg, [1, 0, 2])
        cbb = cbd.ConformalBlocksBundle(liealg, [wt1,wt1,wt2,wt2,wt2,wt2], 5)
        self.assertEqual(99, cbb.getRank(), "Rank incorrect")
        
        liealg = cbd.TypeALieAlgebra(1)
        wt1 = cbd._Weight(liealg, [1])
        cbb = cbd.ConformalBlocksBundle(liealg, [wt1,wt1,wt1,wt1,wt1,wt1,wt1,wt1,wt1,wt1], 1)
        self.assertEqual(1, cbb.getRank(), "Rank incorrect")
        
        liealg = cbd.TypeALieAlgebra(2)
        wt1 = cbd._Weight(liealg, [0, 1])
        cbb = cbd.ConformalBlocksBundle(liealg, [wt1,wt1,wt1,wt1,wt1,wt1,wt1,wt1,wt1], 2)
        #print(cbb.getRank())

    def testConfBlocksDivisors(self):
        liealg = cbd.TypeALieAlgebra(1)
        cbb = cbd.SymmetricConformalBlocksBundle(liealg, [1], 4, 1)
        self.assertEqual(1, cbb.getRank(), "Rank incorrect")
        self.assertEqual(1, cbb.getNormalizedDivisorRay()[0], "Divisor incorrect")
        
        cbb = cbd.SymmetricConformalBlocksBundle(liealg, [1], 6, 1)
        self.assertEqual(1, cbb.getRank(), "Rank incorrect")
        ray = cbb.getNormalizedDivisorRay()
        self.assertEqual(2, ray[0], "Divisor incorrect")
        self.assertEqual(1, ray[1], "Divisor incorrect")

        
        cbb = cbd.SymmetricConformalBlocksBundle(liealg, [1], 8, 1)
        self.assertEqual(1, cbb.getRank(), "Rank incorrect")
        ray = cbb.getNormalizedDivisorRay()
        self.assertEqual(3, ray[0], "Divisor incorrect")
        self.assertEqual(2, ray[1], "Divisor incorrect")
        self.assertEqual(4, ray[2], "Divisor incorrect")
        
        cbb = cbd.SymmetricConformalBlocksBundle(liealg, [1], 10, 1)
        self.assertEqual(1, cbb.getRank(), "Rank incorrect")
        ray = cbb.getNormalizedDivisorRay()
        #print(90*40320*cbb.getDivisor()[0],90*40320*cbb.getDivisor()[1],90*40320*cbb.getDivisor()[2],90*40320*cbb.getDivisor()[3])
        #print(cbb.getNormalizedDivisorRay())

    def testFCurve(self):
        liealg = cbd.TypeALieAlgebra(4)
        wt = cbd._Weight(liealg, [0, 1, 0, 0])
        wt2 = cbd._Weight(liealg, [0, 0, 1, 0])
        wt3 = cbd._Weight(liealg, [0, 0, 0, 1])
        cbb = cbd.ConformalBlocksBundle(liealg, [wt, wt2, wt2, wt2, wt3], 1)
        self.assertEqual(1, cbb.getRank(), "Rank incorrect")
        self.assertEqual(0, cbb.intersect_F_curve([[1, 2], [3], [4], [5]]), "F-curve intersection incorrect")
        self.assertEqual(1, cbb.intersect_F_curve([[1], [2, 3], [4], [5]]), "F-curve intersection incorrect")
        self.assertEqual(2, cbb.intersect_F_curve([[1], [2], [3], [4, 5]]), "F-curve intersection incorrect")

        liealg = cbd.TypeALieAlgebra(6)
        wt = cbd._Weight(liealg, [0, 0, 1, 0, 0, 0])
        wt2 = cbd._Weight(liealg, [0, 0, 0, 1, 0, 0])
        wt3 = cbd._Weight(liealg, [0, 0, 0, 0, 0, 1])
        cbb = cbd.ConformalBlocksBundle(liealg, [wt,wt2,wt2,wt2,wt3], 1)
        self.assertEqual(1, cbb.getRank(), "Rank incorrect")
        self.assertEqual(1, cbb.intersect_F_curve([[1],[2,3],[4],[5]]), "F-curve intersection incorrect")
        self.assertEqual(0, cbb.intersect_F_curve([[1, 2], [3], [4], [5]]), "F-curve intersection incorrect")
        self.assertEqual(3, cbb.intersect_F_curve([[1], [2], [3], [4,5]]), "F-curve intersection incorrect")
        

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testConformalBlocksRanks']
    unittest.main()