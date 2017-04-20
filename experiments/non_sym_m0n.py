from __future__ import division
import fusion_prod.cbbundle as cbd
import cProfile, time, random

#First test
#----------
#
def experiment():
    """
    Computes the rank and divisor of conformal block bundles with random weights.
    :return: Null
    """
    rank = 5
    level = 3
    num_points = 10
    tries = 100

    liealg = cbd.TypeALieAlgebra(rank, store_fusion=True, exact=False)
    A_l = liealg.get_weights(level)
    print("Weight", "Rank", "Divisor")
    for i in range(tries):
        weights = [random.choice(A_l) for i in range(num_points)]
        if sum([sum(liealg._convert_funds_to_epsilons(wt)) for wt in weights]) % (rank+1) != 0: continue
        cbb = cbd.ConformalBlocksBundle(liealg, weights, level)
        if cbb.get_rank() > 0:
            divisor = cbb.get_symmetrized_divisor()
            print(weights, cbb.get_rank(), divisor)
        else:
            print(weights, cbb.get_rank(), 0)


if __name__ == '__main__':
    t0 = time.clock()
    experiment()
    print(time.clock() -t0)
    #cProfile.run('experiment()', sort='cumtime')