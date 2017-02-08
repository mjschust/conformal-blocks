'''
Created on Nov 21, 2016

@author: mjschust
'''
from __future__ import division
import fusion_prod.cbd as cbd
import math
import cProfile, time

#Computes the cosine of the angle between two vectors
def vec_cos(v1,v2):
    l1 = math.sqrt(sum([x*x for x in v1]))
    l2 = math.sqrt(sum([x*x for x in v2]))
    if l1 == 0 or l2 == 0: return 0
    
    dot_sum = 0
    for i in range(len(v1)):
        dot_sum += v1[i]*v2[i]
    return dot_sum/(l1*l2)

#Iterates through symmetric conformal blocks divisors of specified Lie rank, level,
#and number of points, and compares them to a given divisor by computing the
#cosine of the angle between them.
def experiment():
    rank = 5
    level = 3
    num_points = 9
    goal = [1, 1, 2]

    liealg = cbd.TypeALieAlgebra(rank, store_fusion=True)
    print("Weight", "Rank", "Divisor", "Cosine")
    for wt in liealg.get_weights(level):
        cbb = cbd.SymmetricConformalBlocksBundle(liealg, wt, num_points, level)
        if cbb.getRank() == 0: continue
        divisor = cbb.get_norm_sym_divisor_ray()
        print(wt, cbb.getRank(), divisor, vec_cos(divisor, goal))

if __name__ == '__main__':
    t0 = time.clock()
    experiment()
    print(time.clock() -t0)
    #cProfile.run('experiment()', sort='cumtime')

        