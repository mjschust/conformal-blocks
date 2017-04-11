import itertools

vars_str = ""
var_dict = {}
for i in range(1, 5+1):
    for j in range(i+1, 5+1):
        var_str = 'd_' + str(i) + str(j)
        var_dict[var_str] = frozenset((i, j))
        vars_str += var_str + ','

for i in range(1, 5+1):
    for j in range(i+1, 5+1):
        for k in range(j+1, 5+1):
            var_str = 'd_' + str(i) + str(j) + str(k)
            var_dict[var_str] = frozenset((i, j, k))
            vars_str += var_str + ','

vars_str = vars_str[0:-1]

R = PolynomialRing(QQ, 20, names=vars_str)

gen_dict = {}
for mono in R.gens():
    ind = var_dict[str(mono)]
    gen_dict[ind] = mono

comp_dict = {}
A = {1, 2, 3, 4, 5}
for J in gen_dict:
    comp_dict[J] = frozenset(A.difference(J))

ideal_list = []
for J in gen_dict:
    ideal_list.append(gen_dict[J] - gen_dict[comp_dict[J]])

for perm in itertools.permutations([1, 2, 3, 4, 5]):
    i, j, k, l, m = perm[0], perm[1], perm[2], perm[3], perm[4]

    f = gen_dict[frozenset((i, j))] + gen_dict[frozenset((i, j, m))] \
        - gen_dict[frozenset((i, k))] - gen_dict[frozenset((i, k, m))]
    ideal_list.append(f)

    f = gen_dict[frozenset((i, k))] + gen_dict[frozenset((i, k, m))] \
        - gen_dict[frozenset((i, l))] - gen_dict[frozenset((i, l, m))]
    ideal_list.append(f)

for J_1 in gen_dict:
    for J_2 in gen_dict:
        if J_1.issubset(J_2) or J_1.issubset(comp_dict[J_2]) or J_2.issubset(J_1) or J_2.issubset(comp_dict[J_1]):
            continue
        else:
            ideal_list.append(gen_dict[J_1] * gen_dict[J_2])

I = R.ideal(ideal_list)
S = R.quotient_ring(I)
#print(S)
y_1 = 1/2*(gen_dict[frozenset((1,2))] + gen_dict[frozenset((1,3))] + gen_dict[frozenset((1,4))] + \
           gen_dict[frozenset((1,5))]) + 1/6*(gen_dict[frozenset((2,3))] + gen_dict[frozenset((2,4))] + \
           gen_dict[frozenset((2,5))] + gen_dict[frozenset((3,4))] + gen_dict[frozenset((3,5))] + gen_dict[frozenset((4,5))])