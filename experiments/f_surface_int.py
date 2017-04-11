import fusion_prod.cbd as cbd
import random, time, cProfile

def experiment():
    """
    Computes the rank and divisor of conformal block bundles with random weights.
    :return: Null
    """
    rank = 1
    level = 2
    num_points = 5
    f_surface1 = [[1], [2], [3], [4], [5]]
    f_surface2 = f_surface1
    tries = 100

    liealg = cbd.TypeALieAlgebra(rank)
    A_l = liealg.get_weights(level)
    print("Weight", "Rank", "Divisor", "Cosine")
    for i in range(tries):
        weights = [random.choice(A_l) for i in range(num_points)]
        if sum([sum(liealg._convert_funds_to_epsilons(wt)) for wt in weights]) % (rank+1) != 0: continue
        cbb = cbd.ConformalBlocksBundle(liealg, weights, level)
        if(cbb.get_rank() > 0):
            divisor = cbb.get_symmetrized_divisor()
            int_num1 = cbb.intersect_F_surface(f_surface1)
            int_num2 = cbb.intersect_F_surface(f_surface2)
            print(weights, cbb.get_rank(), divisor, int_num1, int_num2, liealg.chern_number(*(weights + [level])))

if __name__ == '__main__':
    t0 = time.clock()
    experiment()
    print(time.clock() -t0)
    #cProfile.run('experiment()', sort='cumtime')