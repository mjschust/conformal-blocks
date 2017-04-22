from __future__ import division
import conformal_blocks.cbbundle as cbd
import cProfile, time, random, subprocess

def experiment():
    """
    Generates Macaulay 2 test cases, runs them using Swinarski's program, then outputs new unit tests if successful
    :return: Null
    """
    rank = 3
    level = 2
    num_points = 5
    tries = 10

    liealg = cbd.TypeALieAlgebra(rank)
    A_l = liealg.get_weights(level)
    m2file = open("TestFCurve.m2", "w")
    m2file.write("loadPackage(\"ConformalBlocks\");\n")
    m2file.write("sl_" + str(rank+1) + " = simpleLieAlgebra(\"A\", " + str(rank) + ");\n")
    test_cases = []
    for i in range(tries):
        weights = [random.choice(A_l) for i in range(num_points)]
        cbb = cbd.ConformalBlocksBundle(liealg, weights, level)
        f_curve = random.choice(cbb.get_F_curves())
        test_cases.append((weights, f_curve))

        wt_str = "{"
        for wt in weights:
            if len(wt) == 1:
                wt_str += "{" + str(wt)[1] + "}, "
            else:
                wt_str += "{" + str(wt)[1:-1] + "}, "
        wt_str = wt_str[:-2] + "}"

        f_str = "{"
        for part in f_curve:
            if len(part) == 1:
                f_str += "{" + str(part)[1] + "}, "
            else:
                f_str += "{" + str(part)[1:-1] + "}, "
        f_str = f_str[:-2] + "}"

        m2file.write("V = conformalBlockVectorBundle(sl_" + str(rank+1) + ", " + str(level)  + ", " + wt_str + ", 0);\n")
        m2file.write("if " + str(cbb.intersect_F_curve(f_curve)) +
                     " != FCurveDotConformalBlockDivisor(" + f_str + ", V) then error(\"Bundle " + "(sl_" + str(rank+1) +
                     ", " + str(level)  + ", " + wt_str + ") incorrect intersection with F-curve " + f_str + "\");\n")

    m2file.write("print(\"OK\");\n")
    m2file.close()

    test_out = subprocess.check_output(["M2", "--script", "TestFCurve.m2"])
    if test_out == "OK\n":
        print("liealg = cbd.TypeALieAlgebra(" + str(rank) + ")")
        for case in test_cases:
            weights = case[0]
            f_curve = case[1]
            cbb = cbd.ConformalBlocksBundle(liealg, weights, level)
            print("cbb = cbd.ConformalBlocksBundle(liealg, " + str(weights) + ", " + str(level) + ")")
            print("f_curve = " + str(f_curve))
            print("self.assertEqual(" + str(cbb.intersect_F_curve(f_curve)) +
                  ", cbb.intersect_F_curve(f_curve))")
            print("")
        print("OK")
    else:
        print(test_out)


if __name__ == '__main__':
    experiment()
