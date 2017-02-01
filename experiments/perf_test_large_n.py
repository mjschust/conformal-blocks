from __future__ import division
import fusion_prod.cbd as cbd
import cProfile

#Tests the performance of rank and divisor calculations for large (>6) number of points
#and small weights.
#Original time: 63 seconds
#After flattening IrrRep: 53 seconds
def experiment():
    rank = 4
    level = 4
    num_points = 9

    liealg = cbd.TypeALieAlgebra(rank, store_fusion=True)
    V = cbd.SymmetricConformalBlocksBundle(liealg, [0,1,1,0], num_points, level)
    print(V.getRank())

if __name__ == '__main__':
    #t0 = time.clock()
    #experiment()
    #print(time.clock() -t0)
    cProfile.run('experiment()', sort='cumtime')