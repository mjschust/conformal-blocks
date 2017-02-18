'''
Created on Nov 19, 2016

@author: mjschust
'''
import unittest
import fusion_prod.cbd as cbd

class Test(unittest.TestCase):


    def testConformalBlocksRanks(self):
        liealg = cbd.TypeALieAlgebra(3)
        wt1 = tuple([0, 1, 4])
        wt2 = tuple([1, 0, 2])
        cbb = cbd.ConformalBlocksBundle(liealg, [wt1,wt1,wt2,wt2,wt2,wt2], 5)
        self.assertEqual(99, cbb.get_rank(), "Rank incorrect")
        
        liealg = cbd.TypeALieAlgebra(1)
        wt1 = tuple([1])
        cbb = cbd.ConformalBlocksBundle(liealg, [wt1,wt1,wt1,wt1,wt1,wt1,wt1,wt1,wt1,wt1], 1)
        self.assertEqual(1, cbb.get_rank(), "Rank incorrect")
        
        liealg = cbd.TypeALieAlgebra(2)
        wt1 = tuple([0, 1])
        cbb = cbd.ConformalBlocksBundle(liealg, [wt1,wt1,wt1,wt1,wt1,wt1,wt1,wt1,wt1], 2)
        #print(cbb.getRank())

        liealg = cbd.TypeALieAlgebra(3)
        cbb = cbd.SymmetricConformalBlocksBundle(liealg, [1, 0, 2], 8, 1)
        self.assertEqual(0, cbb.get_rank(), "Rank incorrect")

    def testSymmetricDivisors(self):
        liealg = cbd.TypeALieAlgebra(1)
        cbb = cbd.SymmetricConformalBlocksBundle(liealg, [1], 4, 1)
        self.assertEqual(1, cbb.get_rank(), "Rank incorrect")
        self.assertEqual(1, cbb.get_norm_sym_divisor_ray()[0], "Divisor incorrect")
        
        cbb = cbd.SymmetricConformalBlocksBundle(liealg, [1], 6, 1)
        self.assertEqual(1, cbb.get_rank(), "Rank incorrect")
        ray = cbb.get_norm_sym_divisor_ray()
        self.assertEqual(2, ray[0], "Divisor incorrect")
        self.assertEqual(1, ray[1], "Divisor incorrect")

        
        cbb = cbd.SymmetricConformalBlocksBundle(liealg, [1], 8, 1)
        self.assertEqual(1, cbb.get_rank(), "Rank incorrect")
        ray = cbb.get_norm_sym_divisor_ray()
        self.assertEqual(3, ray[0], "Divisor incorrect")
        self.assertEqual(2, ray[1], "Divisor incorrect")
        self.assertEqual(4, ray[2], "Divisor incorrect")
        
        cbb = cbd.SymmetricConformalBlocksBundle(liealg, [1], 10, 1)
        self.assertEqual(1, cbb.get_rank(), "Rank incorrect")
        ray = cbb.get_norm_sym_divisor_ray()

    def testConformalBlocksDivisors(self):
        liealg = cbd.TypeALieAlgebra(1)
        cbb = cbd.ConformalBlocksBundle(liealg, [[1], [1], [1], [1]], 1)
        self.assertEqual(1, cbb.get_rank(), "Rank incorrect")
        self.assertEqual(1, cbb.get_norm_sym_divisor_ray()[0], "Divisor incorrect")

        cbb = cbd.ConformalBlocksBundle(liealg, [[1], [1], [1], [1], [1], [1]], 1)
        self.assertEqual(1, cbb.get_rank(), "Rank incorrect")
        ray = cbb.get_norm_sym_divisor_ray()
        self.assertEqual(2, ray[0], "Divisor incorrect")
        self.assertEqual(1, ray[1], "Divisor incorrect")

        cbb = cbd.ConformalBlocksBundle(liealg, [[1], [1], [1], [1], [1], [1], [1], [1]], 1)
        self.assertEqual(1, cbb.get_rank(), "Rank incorrect")
        ray = cbb.get_norm_sym_divisor_ray()
        self.assertEqual(3, ray[0], "Divisor incorrect")
        self.assertEqual(2, ray[1], "Divisor incorrect")
        self.assertEqual(4, ray[2], "Divisor incorrect")

        liealg = cbd.TypeALieAlgebra(5)
        cbb = cbd.ConformalBlocksBundle(liealg, [[0,1,0,0,0],[0,1,0,0,0],[0,1,0,0,0],[0,1,0,0,0],[0,1,0,0,0],[0,1,0,0,0]], 1)
        self.assertEqual(1, cbb.get_rank(), "Rank incorrect")
        ray = cbb.get_norm_sym_divisor_ray()
        self.assertEqual(1, ray[0], "Divisor incorrect")
        self.assertEqual(3, ray[1], "Divisor incorrect")

        liealg = cbd.TypeALieAlgebra(6)
        cbb = cbd.ConformalBlocksBundle(liealg, [[0, 0, 0, 0, 0, 1], [0, 0, 0, 0, 0, 1], [0, 0, 1, 0, 0, 0], [0, 0, 1, 0, 0, 0],
                                                 [0, 1, 0, 0, 0, 0], [1, 0, 0, 0, 0, 0]], 1)
        self.assertEqual(1, cbb.get_rank(), "Rank incorrect")
        ray = cbb.get_norm_sym_divisor_ray()
        self.assertEqual(4, ray[0], "Divisor incorrect")
        self.assertEqual(5, ray[1], "Divisor incorrect")

        cbb = cbd.ConformalBlocksBundle(liealg,
                                        [[0, 0, 0, 0, 0, 1], [0, 0, 0, 0, 0, 1], [0, 0, 1, 0, 0, 0], [0, 0, 1, 0, 0, 0],
                                         [0, 1, 0, 0, 0, 0], [1, 0, 0, 0, 0, 0]], 2)
        self.assertEqual(7, cbb.get_rank(), "Rank incorrect")
        ray = cbb.get_norm_sym_divisor_ray()
        self.assertEqual(1, ray[0], "Divisor incorrect")
        self.assertEqual(3, ray[1], "Divisor incorrect")

        cbb = cbd.ConformalBlocksBundle(liealg,
                                        [[0, 0, 0, 0, 0, 1], [0, 0, 0, 0, 0, 1], [0, 0, 1, 0, 0, 0], [0, 0, 1, 0, 0, 0],
                                         [0, 1, 0, 0, 0, 0], [1, 0, 0, 0, 0, 0]], 3)
        self.assertEqual(8, cbb.get_rank(), "Rank incorrect")
        ray = cbb.get_norm_sym_divisor_ray()
        self.assertEqual(0, ray[0], "Divisor incorrect")
        self.assertEqual(0, ray[1], "Divisor incorrect")

        liealg = cbd.TypeALieAlgebra(3)
        cbb = cbd.ConformalBlocksBundle(liealg, [[0, 1, 4], [0, 1, 4], [1, 0, 1], [1, 0, 1], [1, 0, 1], [1, 0, 1]], 5)
        self.assertEqual(10, cbb.get_rank(), "Rank incorrect")
        ray = cbb.get_norm_sym_divisor_ray()
        self.assertEqual(8, ray[0], "Divisor incorrect")
        self.assertEqual(9, ray[1], "Divisor incorrect")

    def testFCurve(self):
        liealg = cbd.TypeALieAlgebra(4)
        wt = tuple([0, 1, 0, 0])
        wt2 = tuple([0, 0, 1, 0])
        wt3 = tuple([0, 0, 0, 1])
        cbb = cbd.ConformalBlocksBundle(liealg, [wt, wt2, wt2, wt2, wt3], 1)
        self.assertEqual(1, cbb.get_rank(), "Rank incorrect")
        self.assertEqual(0, cbb.intersect_F_curve([[1, 2], [3], [4], [5]]), "F-curve intersection incorrect")
        self.assertEqual(1, cbb.intersect_F_curve([[1], [2, 3], [4], [5]]), "F-curve intersection incorrect")
        self.assertEqual(2, cbb.intersect_F_curve([[1], [2], [3], [4, 5]]), "F-curve intersection incorrect")

        liealg = cbd.TypeALieAlgebra(6)
        wt = tuple([0, 0, 1, 0, 0, 0])
        wt2 = tuple([0, 0, 0, 1, 0, 0])
        wt3 = tuple([0, 0, 0, 0, 0, 1])
        cbb = cbd.ConformalBlocksBundle(liealg, [wt,wt2,wt2,wt2,wt3], 1)
        self.assertEqual(1, cbb.get_rank(), "Rank incorrect")
        self.assertEqual(1, cbb.intersect_F_curve([[1],[2,3],[4],[5]]), "F-curve intersection incorrect")
        self.assertEqual(0, cbb.intersect_F_curve([[1, 2], [3], [4], [5]]), "F-curve intersection incorrect")
        self.assertEqual(3, cbb.intersect_F_curve([[1], [2], [3], [4,5]]), "F-curve intersection incorrect")

    def testChernClass(self):
        liealg = cbd.TypeALieAlgebra(1)
        cbb = cbd.SymmetricConformalBlocksBundle(liealg, [1], 6, 1)
        self.assertEqual(1, cbb.get_rank(), "Rank incorrect")
        print(cbb.get_sym_chern_class())

        weights = [(1,),(1,),(1,),(1,),(0,)]
        self.assertEqual(1, liealg.get_rank(weights, 1))
        self.assertEqual(0, liealg.get_chern_root2((1,), (1,), (1,), (1,), (0,), 1))

        liealg = cbd.TypeALieAlgebra(3)
        cbb = cbd.SymmetricConformalBlocksBundle(liealg, [1,0,1], 6, 2)
        print(cbb.get_rank())
        print(cbb.get_symmetrized_divisor())
        print(cbb.get_sym_chern_class())
        

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testConformalBlocksRanks']
    unittest.main()