from __future__ import division
import fusion_prod.cbd as cbd
import cProfile, time

#First test
#----------
#Tested the performance of rank calculations for r=4, l=4, n=9, weight=[0,1,1,0]
#and small weights.
#Optimization: 63 seconds -> 10.5 sec (profiler time)
#
#Second test (after above optimizations)
#----------
#Tested the performance of rank and divisor calculations for r=5, l=3, n=9
#Original:17 seconds
#After cdef function: 14.5
#After variable cdefs and optimizations: 12.3
#Removing cdefs: 15.3
#Using multi_fusion instead of factoring: 4.5 sec
#Optimized orbit iterator: 3 sec
#
#Third test (after above optimizations)
#----------
#Tested the performance of rank and divisor calculations for r=5, l=4, n=10
#Original:52 seconds
def experiment():
    rank = 5
    level = 4
    num_points = 10

    liealg = cbd.TypeALieAlgebra(rank, store_fusion=True)
    print("Weight", "Rank", "Divisor", "Cosine")
    for wt in liealg.get_weights(level):
        cbb = cbd.SymmetricConformalBlocksBundle(liealg, wt, num_points, level)
        if cbb.get_rank() == 0: continue
        divisor = cbb.get_symmetrized_divisor()
        print(wt, cbb.get_rank(), divisor)

if __name__ == '__main__':
    t0 = time.clock()
    experiment()
    print(time.clock() -t0)
    #cProfile.run('experiment()', sort='cumtime')