# cython: profile=True
from __future__ import division
from collections import defaultdict
import math, fractions, itertools
'''
Created on Nov 10, 2016

@author: mjschust
'''


class ConformalBlocksBundle(object):
    """
    A class representing a conformal blocks vector bundle.
    """

    def __init__(self, liealg, weights, level):
        """
        :param liealg: A SimpleLieAlgebra object.
        :param weights: A list of lists of integers: the weights of the conformal blocks bundle.
        :param level: A positive integer: the level of the conformal blocks bundle.
        """
        self.liealg = liealg
        new_weights = []
        for wt in weights:
            new_weights.append(tuple(wt))
        self.weights = new_weights
        self.level = level
        self._rank = -1

    def getRank(self):
        """
        Computes the rank of the conformal blocks bundle.  The algorithm uses factorization, then
        the fusion product to compute the 3-point ranks.

        :return: An integer: the rank of the bundle.
        """
        if self._rank < 0:
            self._rank = long(round(self._compute_CB_rank(self.weights, self.level)))

        return self._rank

    def _compute_CB_rank(self, weights, level):
        # Find weights with largest and smallest corresponding rep's
        min_dim = max_dim = self.liealg.get_rep_dim(weights[0])
        min_index = max_index = 0
        for i in range(len(weights)):
            dim = self.liealg.get_rep_dim(weights[i])
            if dim < min_dim:
                min_dim = dim
                min_index = i
            if dim > max_dim:
                max_dim = dim
                max_index = i
        # Covers the case when all dimensions are the same
        if min_index == max_index:
            max_index = min_index + 1

        fus_prod = self.liealg.fusion(weights[min_index], weights[max_index], level)
        # indices = min_index, max_index
        # factor_list = [wt for (i, wt) in enumerate(weights) if i not in indices]
        factor_list = []
        for i in range(len(weights)):
            if i != min_index and i != max_index:
                factor_list.append(weights[i])
        multi_fus_prod = self.liealg.multi_fusion(factor_list, level)

        ret_val = 0
        for mu_star in fus_prod:
            mult = fus_prod[mu_star]
            mu = self.liealg.get_dual_weight(mu_star)
            if mu in multi_fus_prod:
                ret_val += mult * multi_fus_prod[mu]

        return ret_val

    #Original version of the above method.  Uses less memory but runs an order of magnitude slower.
    def _alt_compute_CB_rank(self, weights, level):
        # Find weights with largest and smallest corresponding rep's
        min_dim = max_dim = self.liealg.get_rep_dim(weights[0])
        min_index = max_index = 0
        for i in range(len(weights)):
            dim = self.liealg.get_rep_dim(weights[i])
            if dim < min_dim:
                min_dim = dim
                min_index = i
            if dim > max_dim:
                max_dim = dim
                max_index = i
        # Covers the case when all dimensions are the same
        if min_index == max_index:
            max_index = min_index + 1

        fus_prod = self.liealg.fusion(weights[min_index], weights[max_index], level)
        # indices = min_index, max_index
        # factor_list = [wt for (i, wt) in enumerate(weights) if i not in indices]
        factor_list = []
        for i in range(len(weights)):
            if i != min_index and i != max_index:
                factor_list.append(weights[i])

        # Three point case is given by the fusion product
        if len(factor_list) == 1:
            dual_wt3 = self.liealg.get_dual_weight(factor_list[0])
            if dual_wt3 in fus_prod:
                return fus_prod[dual_wt3]
            else:
                return 0

        # If more than three points, factor
        ret_val = 0
        for wt in fus_prod:
            mult = fus_prod[wt]
            if mult > 0:
                ret_val = ret_val + mult * self._alt_compute_CB_rank(factor_list + [wt], level)

        return ret_val

    def get_symmetrized_divisor(self):
        """
        Computes the symmetrized divisor associated to the conformal blocks bundle.

        :return: A list of numbers: the divisor given in the standard basis D_1, D_2,... of
            the symmetric nef cone.
        """

        ret_val = []
        n = len(self.weights)
        weighted_rank = 0
        for wt in self.weights:
            weighted_rank += self.liealg.casimirScalar(wt)
        weighted_rank = self.getRank() * weighted_rank / (n * (n - 1))

        point_indices = [i for i in range(0, n)]
        for i in range(2, n // 2 + 1):
            coord = i * (n - i) * weighted_rank
            sum = 0
            for subset in itertools.combinations(point_indices, i):
                #Could be more efficient here
                wt_list1 = []
                wt_list2 = []
                for j in range(0,n):
                    if j in subset:
                        wt_list1.append(self.weights[j])
                    else:
                        wt_list2.append(self.weights[j])

                prod = self.liealg.multi_fusion(wt_list1, self.level)
                for mu_star in prod.keys():
                    mu = self.liealg.get_dual_weight(mu_star)
                    V1 = ConformalBlocksBundle(self.liealg, wt_list1 + [mu], self.level)
                    V2 = ConformalBlocksBundle(self.liealg, wt_list2 + [mu_star], self.level)
                    sum += self.liealg.casimirScalar(mu) * V1.getRank() * V2.getRank()

            sum = sum*math.factorial(i)*math.factorial(n-i)/math.factorial(n)
            coord = (coord - sum) / (2 * (self.level + self.liealg.dual_coxeter()))
            ret_val.append(coord)

        return ret_val

    def get_norm_sym_divisor_ray(self):
        """
        Computes the symmetrized divisor associated to the conformal blocks bundle and normalizes the
        vector by clearing denominators.

        :return: A list of numbers: the divisor ray given in the standard basis D_1, D_2,... of
            the symmetric nef cone.
        """
        divisor = self.get_symmetrized_divisor()
        n_fact = math.factorial(len(self.weights))
        int_div = [long(round(n_fact * x)) for x in divisor]
        div_gcd = reduce(lambda x, y: fractions.gcd(x, y), int_div)
        if div_gcd > 0:
            return [x // div_gcd for x in int_div]
        else:
            return [x for x in int_div]

    def intersect_F_curve(self, partition):
        """
        Computes the intersection of the divisor associated to this conformal blocks bundle with
        the given F-curve.

        :param partition: A list of 4 lists of integers partitioning the set {1, ..., # points}: the
            F-curve to be intersected.
        :return: An integer: the intersection number.
        """
        ret_val = 0
        wt_list1 = [self.weights[point - 1] for point in partition[0]]
        wt_list2 = [self.weights[point - 1] for point in partition[1]]
        wt_list3 = [self.weights[point - 1] for point in partition[2]]
        wt_list4 = [self.weights[point - 1] for point in partition[3]]

        prod1 = self.liealg.multi_fusion(wt_list1, self.level)
        prod2 = self.liealg.multi_fusion(wt_list2, self.level)
        prod3 = self.liealg.multi_fusion(wt_list3, self.level)
        prod4 = self.liealg.multi_fusion(wt_list4, self.level)

        for wt1 in prod1.keys():
            if prod1[wt1] == 0: continue
            for wt2 in prod2.keys():
                if prod2[wt2] == 0: continue
                for wt3 in prod3.keys():
                    if prod3[wt3] == 0: continue
                    mu_list = [wt1, wt2, wt3]
                    mu_prod = self.liealg.multi_fusion(mu_list, self.level)
                    for wt4 in prod4.keys():
                        if prod4[wt4] == 0: continue
                        if mu_prod[self.liealg.get_dual_weight(wt4)] == 0: continue
                        ret_val += self.liealg.degree(wt1, wt2, wt3, wt4, self.level) * prod1[wt1] * prod2[wt2] * prod3[
                            wt3] * prod4[wt4]

        return ret_val


class SymmetricConformalBlocksBundle(ConformalBlocksBundle):
    """
    A class representing a symmetric conformal blocks vector bundle.
    """

    def __init__(self, liealg, wt, num_points, level):
        """
        :param liealg: A SimpleLieAlgebra object.
        :param wt: A list of integers: the weight of the conformal blocks bundle, repeated at each
            point.
        :param num_points: A positive integer: the number of points of the conformal blocks bundle.
        :param level: A positive integer: the level of the conformal blocks bundle.
        """
        ConformalBlocksBundle.__init__(self, liealg, [wt for i in range(num_points)], level)

    def get_symmetrized_divisor(self):
        """
        Computes the symmetrized divisor associated to the conformal blocks bundle.  Algorithm is
        optimized for the symmetric case.

        :return: A list of numbers: the divisor given in the standard basis D_1, D_2,... of
            the symmetric nef cone.
        """
        ret_val = []
        n = len(self.weights)
        wt = self.weights[0]
        for i in range(2, n // 2 + 1):
            coord = i * (n - i) * self.getRank() * self.liealg.casimirScalar(wt) / (n - 1)
            sum_list = [0]
            self._weightedFactor(wt, wt, 1, i - 1, n - i, sum_list, {})
            coord = (coord - sum_list[0]) / (2 * (self.level + self.liealg.dual_coxeter()))
            ret_val.append(coord)

        return ret_val

    def _weightedFactor(self, wt, wt2, mult, wts_rem, ic, ret_val, rank_dict):
        prod = self.liealg.fusion(wt, wt2, self.level)

        for wt3 in prod.keys():
            if wts_rem > 1:
                self._weightedFactor(wt, wt3, mult * prod[wt3], wts_rem - 1, ic, ret_val, rank_dict)
            else:
                if not wt3 in rank_dict:
                    wt_list = [wt for i in range(ic)]
                    wt_list.append(wt3)
                    cbb = ConformalBlocksBundle(self.liealg, wt_list, self.level)
                    rank_dict[wt3] = cbb.getRank()

                ret_val[0] += self.liealg.casimirScalar(self.liealg.get_dual_weight(wt3)) * mult * prod[wt3] * \
                              rank_dict[wt3]



class IrrRep(object):
    """
    This class represents an irreducible finite dimensional representation of a simple Lie
    algebra.
    """

    def __init__(self, liealg, high_weight):
        """
        :param liealg: A SimpleLieAlgebra object.
        :param high_weight: A list of integers: the highest weight of the representation.
        """

        # if not high_weight.isDominant(): raise ValueError("Highest weight must be dominant")
        self.liealg = liealg
        self.high_weight = tuple(high_weight)
        self._dom_weights = set()
        self._dom_char = {}

    def get_dimension(self):
        """
        Computes the dimension of the representation.  Implements Weyl's dimension formula.

        :return: An integer: the dimension of the representation.
        """
        return self.liealg.get_rep_dim(self.high_weight)


    def get_dominant_character(self):
        """
        Computes the dominant character of the representation.  Implements Freudenthal's recursion
        formula to compute weight multiplicities.  Implementation is heavily influenced by
        the implementation in the LiE system.

        :return: A dictionary where the keys are tuples and the values are positive integers
            corresponding to the multiplicity of the wt space.
        """
        #if self.high_weight in self._dom_char:
        #    return self._dom_char

        pos_roots = self.liealg.get_positive_roots()
        root_level_dict = {}
        for root in pos_roots:
            level = root.get_root_level()
            if level in root_level_dict:
                root_level_dict[level].append(root)
            else:
                root_level_dict[level] = [root]

        level = 0
        weight_level_dict = {}
        weight_level_dict[0] = [self.high_weight]
        self._dom_weights.add(self.high_weight)
        while True:
            done = True
            for key in weight_level_dict.keys():
                if level <= key:
                    done = False
                    break
            if done:
                break

            if not level in weight_level_dict:
                level = level + 1
                continue

            for wt in weight_level_dict[level]:
                for root_lev in root_level_dict.keys():
                    for root in root_level_dict[root_lev]:
                        new_weight = self.liealg._sub_weights(wt, root)
                        if self.liealg.is_dominant(new_weight):
                            if level + root_lev in weight_level_dict:
                                if not new_weight in weight_level_dict[level + root_lev]:
                                    weight_level_dict[level + root_lev].add(new_weight)
                                    self._dom_weights.add(new_weight)
                            else:
                                weight_level_dict[level + root_lev] = set([new_weight])
                                self._dom_weights.add(new_weight)

            level = level + 1

        ret_dict = {}
        ret_dict[self.high_weight] = 1
        sorted_levels = sorted(weight_level_dict.keys())
        for level in sorted_levels:
            for wt in weight_level_dict[level]:
                self._compute_mult(wt, pos_roots)

        return self._dom_char

    def _compute_mult(self, wt, pos_roots):
        '''
        This implements Freudenthal's recursion formula.
        Expects a dominant wt in _dom_weights.
        '''
        if wt in self._dom_char:
            return self._dom_char[wt]
        if wt == self.high_weight:
            self._dom_char[wt] = 1
            return 1

        mult_sum = 0
        for root in pos_roots:
            n = 0
            new_weight = wt
            a = self.liealg.killing_form(wt, root)
            b = self.liealg.killing_form(root, root)
            while True:
                n = n + 1
                new_weight = self.liealg._add_weights(new_weight, root)
                new_dom_weight = self.liealg.reflect_to_chamber(new_weight)
                if not new_dom_weight in self._dom_weights: break

                mult_sum = mult_sum + (a + n * b) * self._compute_mult(new_dom_weight, pos_roots)

        rho = self.liealg.get_rho()
        multiplicity = 2 * mult_sum / (
            self.liealg.length_squared(self.liealg._add_weights(self.high_weight, rho)) - self.liealg.length_squared(
            self.liealg._add_weights(wt, rho)))
        self._dom_char[wt] = multiplicity
        return multiplicity


class SimpleLieAlgebra(object):
    """
    A template class for a simple Lie Algebra.  Objects of this class should be constructed by creating
    an object of the appropriate subclass.  Unimplemented methods must be implemented by
    subclasses of each type.
    """

    def __init__(self, rank, store_fusion=True):
        """
        :param rank: A positive integer: the rank of the Lie algebra.
        :param store_fusion: A boolean: if true the lie algebra will save computed fusion products;
        """
        self.rank = rank
        self.store_fusion = store_fusion
        if store_fusion: self._fusion_dict = {}
        self._pos_roots = []
        self._rep_dim_dict = {}
        self._fte_dict = {}

    def get_rep_dim(self, high_weight):
        """
        Computes the dimension of the representation with given highest weight.
        Implements Weyl's dimension formula.  !!!Currently expects a tuple!!!

        :param high_weight: A list or tuple of non-negative integers of length equal to self.rank:
            represents fundamental weight coordinates of the highest weight of the
            irreducible representation.
        :return: An integer: the dimension of the representation.
        """
        #high_weight = tuple(high_weight)
        if high_weight in self._rep_dim_dict: return self._rep_dim_dict[high_weight]

        lam = high_weight
        rho = self.get_rho()
        pos_roots = self.get_positive_roots()

        numer = 1
        denom = 1
        for root in pos_roots:
            a = self.killing_form(lam, root)
            b = self.killing_form(rho, root)
            numer = numer * (a + b)
            denom = denom * b

        self._rep_dim_dict[high_weight] = numer / denom
        return numer / denom

    def tensor(self, wt1, wt2):
        """
        Computes the tensor product decomposition of the irreducible representations with highest
        weights wt1 and wt2.

        :param wt1: A list or tuple of non-negative integers of length equal to self.rank:
            represents fundamental weight coordinates of the highest weight of the
            irreducible representation.
        :param wt2: Same as wt1.
        :return: A dictionary with keys that are tuples of integers and values that are integers:
            the keys correspond to highest weights and the values are the multiplicities of the
            corresponding representation in the tensor product of wt1 and wt2.
        """

        # Want wt to have larger dimension
        rep1, rep2 = IrrRep(self, wt1), IrrRep(self, wt2)
        if self.get_rep_dim(rep1.high_weight) < self.get_rep_dim(rep2.high_weight):
            rep1, rep2 = rep2, rep1
        wt = rep1.high_weight
        wt2 = rep2.high_weight
        rho = self.get_rho()
        dom_char = rep2.get_dominant_character()
        lam_rho_sum = self._add_weights(wt, rho)
        ret_dict = {}

        # Traverse entire character
        for dom_weight in dom_char.keys():
            for orbit_weight in self.get_orbit_iter(dom_weight):
                new_sum = self._add_weights(lam_rho_sum, orbit_weight)
                new_dom_weight, parity = self.reflect_to_chamber_with_parity(new_sum)
                new_dom_weight = self._sub_weights(new_dom_weight, rho)
                if not self.is_dominant(new_dom_weight): continue
                if new_dom_weight in ret_dict:
                    ret_dict[new_dom_weight] = ret_dict[new_dom_weight] + dom_char[dom_weight] * parity
                else:
                    ret_dict[new_dom_weight] = dom_char[dom_weight] * parity

        return ret_dict

    def fusion(self, wt1, wt2, ell):
        """
        Computes the fusion product decomposition of the irreducible representations with highest
        weights wt1 and wt2 at level ell.

        :param wt1: A list or tuple of non-negative integers of length equal to self.rank:
            represents fundamental weight coordinates of the highest weight of the
            irreducible representation.
        :param wt2: Same as wt1.
        :param ell: A positive integer: corresponds to the level of the fusion product.
        :return:  A dictionary with keys that are tuples of integers and values that are integers:
            the keys correspond to highest weights and the values are the multiplicities of the
            corresponding representation in the fusion product of wt1 and wt2.
        """

        #rep1, rep2 = IrrRep(self, wt1), IrrRep(self, wt2)
        wt1, wt2 = tuple(wt1), tuple(wt2)
        if self.store_fusion and (wt1, wt2, ell) in self._fusion_dict:
            return self._fusion_dict[(wt1, wt2, ell)]
        rep1, rep2 = IrrRep(self, wt1), IrrRep(self, wt2)

        ten_decom = self.tensor(wt1, wt2)
        ret_dict = {}
        rho = self.get_rho()

        for wt in ten_decom.keys():
            if self.get_level(wt) == ell + 1: continue

            wt_rho = self._add_weights(wt, rho)
            new_weight, parity = self.reflect_to_alcove_with_parity(wt_rho, ell + self.rank + 1)
            lev_ell_weight = self._sub_weights(new_weight, rho)
            if not self.is_dominant(lev_ell_weight) or self.get_level(lev_ell_weight) > ell: continue

            if lev_ell_weight in ret_dict:
                ret_dict[lev_ell_weight] = ret_dict[lev_ell_weight] + ten_decom[wt] * parity
            else:
                ret_dict[lev_ell_weight] = ten_decom[wt] * parity

        if self.store_fusion:
            self._fusion_dict[(rep1.high_weight, rep2.high_weight, ell)] = ret_dict

        return ret_dict

    def multi_fusion(self, wts, level):
        """
        Computes the fusion product of a list of representations.

        :param wts: A list of lists of integers: the list of weights.
        :param level: A positive integer: corresponds to the level of the fusion product.
        :return: A dictionary with keys that are tuples of integers and values that are integers:
            the keys correspond to highest weights and the values are the multiplicities of the
            corresponding representation in the fusion product of the weights.
        """
        wt_dict = defaultdict(int)
        rem_wts = list(wts)
        cur_wt = rem_wts.pop()
        wt_dict[cur_wt] = 1
        while (len(rem_wts) > 0):
            new_wt_dict = defaultdict(int)
            cur_wt = rem_wts.pop()
            for wt in wt_dict.keys():
                prod = self.fusion(cur_wt, wt, level)
                for wt2 in prod.keys():
                    new_wt_dict[wt2] += wt_dict[wt] * prod[wt2]

            wt_dict = new_wt_dict

        return wt_dict

    def degree(self, wt1, wt2, wt3, wt4, level):
        """
        Computes the degree of a four-point conformal blocks vector bundle.  Implements Fakhruddin's
        formula.

        :param wt1: A list of integers: a weight of the bundle.
        :param wt2: A list of integers: a weight of the bundle.
        :param wt3: A list of integers: a weight of the bundle.
        :param wt4: A list of integers: a weight of the bundle.
        :param level: A positive integer: the level of the bundle.
        :return: A positive integer: the degree of the bundle.
        """
        cbb = ConformalBlocksBundle(self, [wt1, wt2, wt3, wt4], level)
        ret_val = cbb.getRank() * (
        self.casimirScalar(wt1) + self.casimirScalar(wt2) + self.casimirScalar(wt3) + self.casimirScalar(wt4))

        sum = 0
        prod1 = self.fusion(wt1, wt2, level)
        prod2 = self.fusion(wt3, wt4, level)
        for mu in prod1.keys():
            mu_star = self.get_dual_weight(mu)
            if mu_star in prod2:
                sum += self.casimirScalar(mu_star) * prod1[mu] * prod2[mu_star]
        prod1 = self.fusion(wt1, wt3, level)
        prod2 = self.fusion(wt2, wt4, level)
        for mu in prod1.keys():
            mu_star = self.get_dual_weight(mu)
            if mu_star in prod2:
                sum += self.casimirScalar(mu_star) * prod1[mu] * prod2[mu_star]
        prod1 = self.fusion(wt1, wt4, level)
        prod2 = self.fusion(wt2, wt3, level)
        for mu in prod1.keys():
            mu_star = self.get_dual_weight(mu)
            if mu_star in prod2:
                sum += self.casimirScalar(mu_star) * prod1[mu] * prod2[mu_star]
        ret_val -= sum
        ret_val = ret_val / (2 * (level + self.dual_coxeter()))

        return int(round(ret_val))


    def killing_form(self, wt1, wt2):
        """
        Computes the Killing form product of two weights.

        :param wt1: A list of numbers: a weight.
        :param wt2: A list of numbers: a weight.
        :return: A number: the Killing form product of wt1 and wt2.
        """
        raise NotImplementedError

    def length_squared(self, wt):
        """
        Computes the squared length of a weight using the Killing form of the algebra.

        :param wt: A list of numbers: a weight.
        :return: A positive number: the squared length of wt.
        """
        return self.killing_form(wt, wt)

    def casimirScalar(self, wt):
        """
        Computes the Casimir scalar of a weight.

        :param wt: A list of numbers: a weight.
        :return: A number: the Casimir scalar of wt.
        """
        twoRho = tuple([2 for i in range(self.rank)])
        wt2 = self._add_weights(wt, twoRho)
        return self.killing_form(wt, wt2)

    def dual_coxeter(self):
        """
        Computes the dual coxeter number of the Lie algebra.

        :return: An integer: the dual coxeter number.
        """
        raise NotImplementedError

    def is_dominant(self, wt):
        """
        Checks if the weight is dominant.

        :param wt: A list of numbers; a weight of the Lie algebra.
        :return: A boolean.
        """
        for coord in wt:
            if coord < 0:
                return False

        return True

    def get_level(self, wt):
        """
        Computes the level of the weight, which is defined as the Killing product of wt and
        the highest root of the algebra.

        :param wt:  A list of numbers: a weight.
        :return: A number: the level of the weight.
        """
        raise NotImplementedError

    def get_dual_weight(self, wt):
        """
        Computes the highest weight of the contragredient representation associated to wt.

        :param wt: A list of positive integers: a dominant integral weight.
        :return: A list of positive integers: a dominant integral weight.
        """
        raise NotImplementedError

    def get_rho(self):
        """
        Computes one-half the sum of the positive roots of the algebra, which is usually denoted
        by rho.

        :return: A list of positive integers: the weight rho.
        """
        ret_coords = []

        for i in range(self.rank):
            ret_coords.append(1)

        return tuple(ret_coords)

    def get_positive_roots(self):
        """
        Computes a list of all positive roots of the Lie algebra.

        :return: A list of weights: the positive roots of the algebra.
        """
        raise NotImplementedError

    def get_weights(self, level):
        """
        Computes a list of all weights of level less than or equal to the given level.

        :param level: A positive integer.
        :return: A list of weights with level less than level.
        """
        return [tuple(coords) for coords in self._get_weights(level, self.rank)]

    def _get_weights(self, level, rank):
        ret_list = []
        if rank == 1:
            for i in range(level + 1):
                ret_list.append([i])
        else:
            r_minus_one_list = self._get_weights(level, rank - 1)
            for coord in r_minus_one_list:
                for i in range(level - reduce(lambda x, y: x + y, coord) + 1):
                    ret_list.append(coord + [i])

        return ret_list

    def reflect_to_chamber(self, wt):
        """
        Reflects a weight into the dominant chamber of the the Lie algebra.

        :param wt: A list of numbers: a weight.
        :return: A list of non-negative numbers: a dominant weight.
        """
        raise NotImplementedError

    def reflect_to_chamber_with_parity(self, wt):
        """
        Reflects a weight into the dominant chamber of the the Lie algebra, keeping track of the
        parity of the reflections.

        :param wt: A list of numbers: a weight.
        :return: A list of non-negative numbers: a dominant weight; and a number equal to +1 or -1.
        """
        raise NotImplementedError

    def reflect_to_alcove_with_parity(self, wt, ell):
        """
        Reflects a weight into the level ell fundamental chamber of the the Lie algebra, keeping track of the
        parity of the reflections.

        :param wt: A list of numbers: a weight.
        :param ell: A positive integer: the level.
        :return: A list of non-negative numbers: a dominant weight; and a number equal to +1 or -1.
        """
        raise NotImplementedError

    def get_orbit_iter(self, wt):
        """
        Returns an iterable object that iterates through the Weyl group orbit of the given weight.

        :param wt: A list of numbers: a weight.
        :return: An iterator object.
        """
        raise NotImplementedError

    def _convert_funds_to_epsilons(self, coords):
        '''
        '''
        raise NotImplementedError

    def _convert_epsilons_to_funds(self, coords):
        '''
        '''
        raise NotImplementedError

    def _convert_funds_to_roots(self, coords):
        '''
        '''
        raise NotImplementedError

    def _convert_roots_to_funds(self, coords):
        '''
        '''
        raise NotImplementedError

    def _add_weights(self, wt1, wt2):
        '''
        Adds two weights and returns the sum as a new weight object
        '''
        ret_coords = []

        for i in range(len(wt1)):
            ret_coords.append(wt1[i] + wt2[i])

        return tuple(ret_coords)

    def _sub_weights(self, wt1, wt2):
        '''
        Subtracts wt1 from wt2 the difference as a new weight object
        '''
        ret_coords = []

        for i in range(len(wt1)):
            ret_coords.append(wt1[i] - wt2[i])

        return tuple(ret_coords)



class TypeALieAlgebra(SimpleLieAlgebra):
    """
    A type A Lie algebra.
    """

    def killing_form(self, wt1, wt2):
        ret_val = 0
        ep_coords1 = self._convert_funds_to_epsilons(wt1)
        ep_coords2 = self._convert_funds_to_epsilons(wt2)

        for i in range(self.rank + 1):
            ret_val = ret_val + ep_coords1[i] * ep_coords2[i]

        n = sum(ep_coords1) * sum(ep_coords2)

        ret_val = ret_val - n / (self.rank + 1)

        return ret_val

    def dual_coxeter(self):
        return self.rank + 1

    def get_level(self, wt):
        return sum(wt)

    def get_dual_weight(self, wt):
        return tuple(wt[::-1])

    def get_positive_roots(self):
        if len(self._pos_roots) > 0: return self._pos_roots

        ret_list = []
        coords = []
        for i in range(self.rank):
            coords.append(0)

        for i in range(len(coords)):
            for j in range(i, len(coords)):
                coords[j] = 1
                ret_list.append(_Root(self, list(coords)))

            for j in range(i, len(coords)):
                coords[j] = 0

        self._pos_roots = ret_list
        return ret_list

    def reflect_to_chamber(self, wt):
        ret_coords = self._insertsort(self._convert_funds_to_epsilons(wt))

        for i in range(len(ret_coords)):
            ret_coords[i] = ret_coords[i] - ret_coords[-1]

        return tuple(self._convert_epsilons_to_funds(ret_coords))

    def _insertsort(self, coords):
        ret_list = list(coords)

        for i in range(1, len(ret_list)):
            j = i
            while j > 0 and ret_list[j - 1] < ret_list[j]:
                ret_list[j - 1], ret_list[j] = ret_list[j], ret_list[j - 1]
                j = j - 1

        return ret_list

    def reflect_to_chamber_with_parity(self, wt):
        ret_coords, parity = self._insertsort_parity(self._convert_funds_to_epsilons(wt))

        for i in range(len(ret_coords)):
            ret_coords[i] = ret_coords[i] - ret_coords[-1]

        return tuple(self._convert_epsilons_to_funds(ret_coords)), parity

    def _insertsort_parity(self, coords):
        ret_list = list(coords)

        parity = 1
        for i in range(1, len(ret_list)):
            j = i
            while j > 0 and ret_list[j - 1] < ret_list[j]:
                ret_list[j - 1], ret_list[j] = ret_list[j], ret_list[j - 1]
                parity = parity * -1
                j = j - 1

        return ret_list, parity

    def reflect_to_alcove_with_parity(self, wt, ell):
        ret_coords, parity = self._insertsort_parity(self._convert_funds_to_epsilons(wt))
        ret_coords = [x - ret_coords[-1] for x in ret_coords]

        while (ret_coords[0] > ell):
            ret_coords[-1] = ret_coords[0] - ell
            ret_coords[0] = ell
            ret_coords, fin_parity = self._insertsort_parity(ret_coords)
            ret_coords = [x - ret_coords[-1] for x in ret_coords]
            parity = parity * -1 * fin_parity

        return tuple(self._convert_epsilons_to_funds(ret_coords)), parity

    def get_orbit_iter(self, wt):
        return self._TypeAOrbitIterator(self, wt)

    def _convert_funds_to_epsilons(self, coords):
        if coords in self._fte_dict: return self._fte_dict[coords]

        ret_coords = [0]
        part = 0
        for i in reversed(range(len(coords))):
            part += coords[i]
            ret_coords.insert(0, part)

        self._fte_dict[coords] = ret_coords
        return ret_coords

    def _convert_epsilons_to_funds(self, coords):

        ret_coords = []
        for i in range(len(coords) - 1):
            ret_coords.append(coords[i] - coords[i + 1])

        return ret_coords

    def _convert_roots_to_funds(self, coords):
        if len(coords) == 1: return [2 * coords[0]]

        ret_coords = []
        ret_coords.append(2 * coords[0] - coords[1])
        for i in range(1, len(coords) - 1):
            ret_coords.append(2 * coords[i] - coords[i + 1] - coords[i - 1])

        ret_coords.append(2 * coords[-1] - coords[-2])
        return ret_coords

    class _TypeAOrbitIterator(object):
        '''
        Optimized iterator object that traverses the Weyl group orbit of a given weight.

        Attributes:
            no public attributes
        '''
        def __init__(self, liealg, wt):
            ep_coords = liealg._convert_funds_to_epsilons(liealg.reflect_to_chamber(wt))

            #Construct list of items and multiplicities
            self._item_list = [ep_coords[0]]
            cur_item = ep_coords[0]
            rem_list = [0]
            for item in ep_coords:
                if item < cur_item:
                    self._item_list.append(item)
                    rem_list.append(1)
                    cur_item = item
                else:
                    rem_list[-1] += 1

            #Contruct matrix of remaining items, and initial index list
            index_list = []
            rem_mat = [list(rem_list)]
            for i in range(len(ep_coords)):
                j = 0
                while rem_mat[i][j] == 0:
                    j += 1
                index_list.append(j)
                rem_mat.append(list(rem_mat[i]))
                rem_mat[i+1][j] -= 1

            self._index_list = index_list
            self._rem_mat = rem_mat
            self.done = False
            self.liealg = liealg

        def __iter__(self):
            return self

        def next(self):
            if self.done: raise StopIteration()
            r = len(self._index_list)
            num_items = len(self._item_list)

            #Construct new weight
            ep_coords = []
            for index in self._index_list:
                ep_coords.append(self._item_list[index])
            ret_val = tuple(self.liealg._convert_epsilons_to_funds(ep_coords))

            #Find index to increment
            i = r-2
            j = 0
            while i >= 0:
                j = self._index_list[i] + 1
                while j < num_items:
                    if self._rem_mat[i][j] > 0: break
                    j += 1
                if j < num_items: break
                i -= 1

            #If we're finished, return the last weight
            if i < 0:
                self.done = True
                return ret_val

            #Increment indices
            self._index_list[i] = j
            self._rem_mat[i+1] = list(self._rem_mat[i])
            self._rem_mat[i+1][j] -= 1
            i += 1
            while i < r:
                j = 0
                while self._rem_mat[i][j] == 0:
                    j += 1
                self._index_list[i] = j
                self._rem_mat[i + 1] = list(self._rem_mat[i])
                self._rem_mat[i + 1][j] -= 1
                i += 1

            return ret_val


class _Root(tuple):
    """
    This class represents an element of the root lattice of the given simple Lie algebra.
    """

    def __new__(cls, liealg, coords):
        """

        :param liealg:
        :param coords:
        :return:
        """
        return super(_Root, cls).__new__(cls, liealg._convert_roots_to_funds(coords))

    def __init__(self, liealg, coords):
        """
        :param liealg: A SimpleLieAlgebra object.
        :param coords: A list of integers: coordinates for the weight in terms of the basis of simple roots.
        """
        self.liealg = liealg
        self.root_coords = coords

    def get_root_level(self):
        ret_val = 0
        for coord in self.root_coords:
            ret_val = ret_val + coord

        return ret_val
