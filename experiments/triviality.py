'''
Created on Nov 21, 2016

@author: mjschust
'''
from __future__ import division
import time
import conformal_blocks.cbbundle as cbd
import math, cProfile

#Computes all non-trivial 4-point conformal blocks divisors of specified Lie rank and level.
def experiment():
    rank = 5
    level = 4

    liealg = cbd.TypeALieAlgebra(rank, store_fusion=True)
    print("Weight", "Rank", "Divisor")
    for wt in liealg.get_weights(level):
        cbb = cbd.SymmetricConformalBlocksBundle(liealg, wt, 4, level)
        if cbb.get_rank() == 0: continue
        #tot_weight = wt.fund_coords[0] + 2*wt.fund_coords[1] + 3*wt.fund_coords[2]
        #if tot_weight <= level: continue
        #tot_weight = 3*wt.fund_coords[0] + 2*wt.fund_coords[1] + wt.fund_coords[2]
        #if tot_weight <= level: continue
        #if level >= tot_weight // (r+1): continue
        divisor = cbb.get_symmetrized_divisor()
        if divisor[0] == 0: continue
        print(wt, cbb.get_rank(), divisor)

if __name__ == '__main__':
    #t0 = time.clock()
    experiment()
    #print(time.clock() -t0)
    #cProfile.run('experiment()', sort='cumtime')
