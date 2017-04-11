from __future__ import division
import fusion_prod.cbd as cbd
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
    level = 2
    num_points = 5
    tries = 100

    liealg = cbd.TypeALieAlgebra(rank, store_fusion=True)
    A_l = liealg.get_weights(level)
    print("Weight", "Rank", "Divisor", "Cosine")
    for i in range(tries):
        weights = [random.choice(A_l) for i in range(num_points)]
        if sum([sum(liealg._convert_funds_to_epsilons(wt)) for wt in weights]) % (rank+1) != 0: continue
        cbb = cbd.ConformalBlocksBundle(liealg, weights, level)
        if cbb.get_rank() > 0:
            divisor = cbb.get_symmetrized_divisor()
            chern_number = liealg.chern_number(*(weights + [level]))
            print(weights, cbb.get_rank(), divisor, chern_number)
        else:
            print(weights, cbb.get_rank(), 0, 0)


if __name__ == '__main__':
    t0 = time.clock()
    experiment()
    print(time.clock() -t0)
    #cProfile.run('experiment()', sort='cumtime')