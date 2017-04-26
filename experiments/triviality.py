'''
Created on Nov 21, 2016

@author: mjschust
'''
from __future__ import division
import conformal_blocks.cbbundle as cbd
import math, cProfile, time
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

#Computes all non-trivial 4-point conformal blocks divisors of specified Lie rank and level.
def experiment():
    rank = 3
    level = 10

    liealg = cbd.TypeALieAlgebra(rank, store_fusion=True, exact=False)
    print("Weight", "Rank", "Divisor")
    trivial_x = []
    trivial_y = []
    trivial_z = []
    nontrivial_x = []
    nontrivial_y = []
    nontrivial_z = []
    for wt in liealg.get_weights(level):
        cbb = cbd.SymmetricConformalBlocksBundle(liealg, wt, 4, level)
        if cbb.get_rank() == 0: continue
        #tot_weight = wt.fund_coords[0] + 2*wt.fund_coords[1] + 3*wt.fund_coords[2]
        #if tot_weight <= level: continue
        #tot_weight = 3*wt.fund_coords[0] + 2*wt.fund_coords[1] + wt.fund_coords[2]
        #if tot_weight <= level: continue
        #if level >= tot_weight // (r+1): continue
        divisor = cbb.get_symmetrized_divisor()
        if divisor[0] == 0:
            trivial_x.append(wt[0])
            trivial_y.append(wt[2])
            trivial_z.append(wt[1])
        else:
            nontrivial_x.append(wt[0])
            nontrivial_y.append(wt[2])
            nontrivial_z.append(wt[1])

        print(wt, cbb.get_rank(), divisor[0])

    # Plot the results
    #fig, ax = plt.subplots()
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(trivial_x, trivial_y, zs=trivial_z, c='black', label="Trivial divisor")
    ax.scatter(nontrivial_x, nontrivial_y, zs=nontrivial_z, c='red', label="Non-trivial divisor")
    ax.set_xlabel('a_1')
    ax.set_ylabel('a_3')
    ax.set_zlabel('a_2')
    ax.legend()
    ax.grid(True)
    plt.show()

if __name__ == '__main__':
    #t0 = time.clock()
    experiment()
    #print(time.clock() -t0)
    #cProfile.run('experiment()', sort='cumtime')
