from __future__ import division
from collections import defaultdict
import math, fractions
# cython: profile=True
'''
Created on Nov 10, 2016

@author: mjschust
'''


class ConformalBlocksBundle(object):
    '''
    A class representing a conformal blocks vector bundle.

    Attributes:
        liealg: simple lie algebra
        weights: list of weights
        level: positive integer
    '''

    def __init__(self, liealg, weights, level):
        '''
        Constructor
        '''
        self.liealg = liealg
        new_weights = []
        for wt in weights:
            new_weights.append(Weight(liealg, wt))
        self.weights = new_weights
        self.level = level
        self._rank = -1

    def getRank(self):
        '''
        Computes the rank of the conformal blocks bundle
        '''
        if self._rank < 0:
            self._rank = self._compute_CB_rank(self.weights, self.level)

        return self._rank

    def _compute_CB_rank(self, weights, level):
        # Find weights with largest and smallest corresponding rep's
        rep = IrrRep(self.liealg, weights[0])
        min_dim = max_dim = rep.get_dimension()
        min_index = max_index = 0
        for i in range(len(weights)):
            rep = IrrRep(self.liealg, weights[i])
            dim = rep.get_dimension()
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
        indices = min_index, max_index
        factor_list = [wt for (i, wt) in enumerate(weights) if i not in indices]

        # Three point case is given by the fusion product
        if len(factor_list) == 1:
            dual_wt3 = self.liealg.get_dual_weight(factor_list[0])
            if dual_wt3 in fus_prod:
                return fus_prod[dual_wt3]
            else:
                return 0

        # If more than three points, factor
        ret_val = 0
        for wt in fus_prod.keys():
            if fus_prod[wt] > 0:
                ret_val = ret_val + fus_prod[wt] * self._compute_CB_rank(factor_list + [wt], level)

        return ret_val

    def intersect_F_curve(self, partition):
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
    '''
    A class representing a conformal blocks vector bundle.

    Attributes:
        liealg: simple lie algebra
        weights: list of weights
        level: positive integer
    '''

    def __init__(self, liealg, wt, num_points, level):
        '''
        Constructor
        '''
        ConformalBlocksBundle.__init__(self, liealg, [wt for i in range(num_points)], level)

    def getDivisor(self):
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

    def getNormalizedDivisorRay(self):
        divisor = self.getDivisor()
        n_fact = math.factorial(len(self.weights))
        int_div = [long(n_fact * x) for x in divisor]
        div_gcd = reduce(lambda x, y: fractions.gcd(x, y), int_div)
        if div_gcd > 0:
            return [x // div_gcd for x in int_div]
        else:
            return [x for x in int_div]

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
    '''
    This class represents an irreducible finite dimensional representation of a simple Lie
    algebra.

    Attributes:
    liealg: object of type SimpleLieAlgebra; the simple Lie algebra acting on the representation
    high_weight: highest weight of the representation
    '''

    def __init__(self, liealg, high_weight):
        '''
        Constructor for the class.

        Parameter:
        high_weigtht: dominant integral weight
        '''

        # if not high_weight.isDominant(): raise ValueError("Highest weight must be dominant")
        self.liealg = liealg
        self.high_weight = Weight(liealg, high_weight)
        self._dom_weights = set()
        self._dom_char = {}

    def get_dimension(self):
        if self.high_weight in self.liealg._rep_dim_dict: return self.liealg._rep_dim_dict[self.high_weight]

        lam = self.high_weight
        rho = self.liealg.get_rho()
        pos_roots = self.liealg.get_positive_roots()

        numer = 1
        denom = 1
        for root in pos_roots:
            a = self.liealg.killing_form(lam, root)
            b = self.liealg.killing_form(rho, root)
            numer = numer * (a + b)
            denom = denom * b

        self.liealg._rep_dim_dict[self.high_weight] = numer / denom
        return numer / denom

    def get_dominant_character(self):
        '''
        Returns the dominant character of the representation.  Implementation is heavily influenced by
        the implementation in the LiE system.

        Return value:
        A dictionary where the keys are Weight objects and the values are positive integers
        corresponding to the multiplicity of the wt space.
        '''
        if self.high_weight in self._dom_char:
            return self._dom_char

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
                        new_weight = self.liealg.sub_weights(wt, root)
                        if new_weight.isDominant():
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
                new_weight = self.liealg.add_weights(new_weight, root)
                new_dom_weight = self.liealg.reflect_to_chamber(new_weight)
                if not new_dom_weight in self._dom_weights: break

                mult_sum = mult_sum + (a + n * b) * self._compute_mult(new_dom_weight, pos_roots)

        rho = self.liealg.get_rho()
        multiplicity = 2 * mult_sum / (
        self.liealg.length_squared(self.liealg.add_weights(self.high_weight, rho)) - self.liealg.length_squared(
            self.liealg.add_weights(wt, rho)))
        self._dom_char[wt] = multiplicity
        return multiplicity


class SimpleLieAlgebra(object):
    '''
    A template class for a simple Lie Algebra.  Unimplemented methods must be implemented by
    subclasses of each type.  Objects of this class should be constructed by creating
    an object of the appropriate subclass.

    Attributes:
        rank: a positive integer representing the Lie rank of the algebra.
        store_fusion: if true the lie algebra will save computed fusion products
    '''

    def __init__(self, rank, store_fusion=False):
        '''
        Constructor
        '''
        self.rank = rank
        self.store_fusion = store_fusion
        if store_fusion: self._fusion_dict = {}
        self._pos_roots = []
        self._rep_dim_dict = {}

    def convert_funds_to_epsilons(self, coords):
        '''
        '''
        raise NotImplementedError

    def convert_epsilons_to_funds(self, coords):
        '''
        '''
        raise NotImplementedError

    def convert_funds_to_roots(self, coords):
        '''
        '''
        raise NotImplementedError

    def convert_roots_to_funds(self, coords):
        '''
        '''
        raise NotImplementedError

    def killing_form(self, wt1, wt2):
        '''
        '''
        raise NotImplementedError

    def reflect_to_chamber(self, wt):
        '''
        '''
        raise NotImplementedError

    def reflect_to_chamber_with_parity(self, wt):
        '''
        '''
        raise NotImplementedError

    def reflect_to_alcove_with_parity(self, wt, ell):
        '''
        '''
        raise NotImplementedError

    def get_positive_roots(self):
        '''
        Returns a list of all positive roots of the Lie algebra.
        '''
        raise NotImplementedError

    def get_orbit_iter(self, wt):
        '''
        Returns an iterable object that iterates through the Weyl group orbit of the given weight.

        Parameters:
            weight: weight whose orbit we want to traverse

        Returns:
            iterator object
        '''
        raise NotImplementedError

    def get_dual_weight(self, wt):
        '''
        '''
        return NotImplementedError

    def dual_coxeter(self):
        '''
        '''
        return NotImplementedError

    def get_level(self, wt):
        '''
        '''
        return NotImplementedError

    def length_squared(self, wt):
        return self.killing_form(wt, wt)

    def add_weights(self, wt1, wt2):
        '''
        Adds two weights and returns the sum as a new weight object
        '''
        ret_coords = []

        for i in range(len(wt1)):
            ret_coords.append(wt1[i] + wt2[i])

        return Weight(self, ret_coords)

    def sub_weights(self, wt1, wt2):
        '''
        Subtracts wt1 from wt2 the difference as a new weight object
        '''
        ret_coords = []

        for i in range(len(wt1)):
            ret_coords.append(wt1[i] - wt2[i])

        return Weight(self, ret_coords)

    def get_rho(self):
        ret_coords = []

        for i in range(self.rank):
            ret_coords.append(1)

        return Weight(self, ret_coords)

    def casimirScalar(self, wt):
        twoRho = Weight(self, [2 for i in range(self.rank)])
        wt2 = self.add_weights(wt, twoRho)
        return self.killing_form(wt, wt2)

    def get_weights(self, level):
        return [Weight(self, coords) for coords in self._get_weights(level, self.rank)]

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

    def tensor(self, wt1, wt2):
        """
        Computes the tensor product decomposition of the irreducible representations with highest
        weights wt1 and wt2.

        :param wt1: A list or tuple of non-negative integers of length equal to self.rank;
            represents fundamental weight coordinates of the highest weight of the
            irreducible representation.
        :param wt2: Same as wt1.
        :return: A dictionary with keys that are tuples of integers and values that are integers;
            the keys correspond to highest weights and the values are the multiplicities of the
            corresponding representation in the tensor product of wt1 and wt2.
        """

        # Want wt to have larger dimension
        rep1, rep2 = IrrRep(self, wt1), IrrRep(self, wt2)
        if rep1.get_dimension() < rep2.get_dimension():
            rep1, rep2 = rep2, rep1
        wt = rep1.high_weight
        wt2 = rep2.high_weight
        rho = self.get_rho()
        dom_char = rep2.get_dominant_character()
        lam_rho_sum = self.add_weights(wt, rho)
        ret_dict = {}

        # Traverse entire character
        for dom_weight in dom_char.keys():
            for orbit_weight in self.get_orbit_iter(dom_weight):
                new_sum = self.add_weights(lam_rho_sum, orbit_weight)
                new_dom_weight, parity = self.reflect_to_chamber_with_parity(new_sum)
                new_dom_weight = self.sub_weights(new_dom_weight, rho)
                if not new_dom_weight.isDominant(): continue
                if new_dom_weight in ret_dict:
                    ret_dict[new_dom_weight] = ret_dict[new_dom_weight] + dom_char[dom_weight] * parity
                else:
                    ret_dict[new_dom_weight] = dom_char[dom_weight] * parity

        return ret_dict

    def fusion(self, wt1, wt2, ell):
        """
        Computes the fusion product decomposition of the irreducible representations with highest
        weights wt1 and wt2 at level ell.

        :param wt1: A list or tuple of non-negative integers of length equal to self.rank;
            represents fundamental weight coordinates of the highest weight of the
            irreducible representation.
        :param wt2: Same as wt1.
        :param ell: A positive integer; corresponds to the level of the fusion product.
        :return:  A dictionary with keys that are tuples of integers and values that are integers;
            the keys correspond to highest weights and the values are the multiplicities of the
            corresponding representation in the fusion product of wt1 and wt2.
        """

        rep1, rep2 = IrrRep(self, wt1), IrrRep(self, wt2)
        if self.store_fusion and (rep1.high_weight, rep2.high_weight, ell) in self._fusion_dict:
            return self._fusion_dict[(rep1.high_weight, rep2.high_weight, ell)]

        ten_decom = self.tensor(wt1, wt2)
        ret_dict = {}
        rho = self.get_rho()

        for wt in ten_decom.keys():
            if self.get_level(wt) == ell + 1: continue

            wt_rho = self.add_weights(wt, rho)
            new_weight, parity = self.reflect_to_alcove_with_parity(wt_rho, ell + self.rank + 1)
            lev_ell_weight = self.sub_weights(new_weight, rho)
            if not lev_ell_weight.isDominant() or self.get_level(lev_ell_weight) > ell: continue

            if lev_ell_weight in ret_dict:
                ret_dict[lev_ell_weight] = ret_dict[lev_ell_weight] + ten_decom[wt] * parity
            else:
                ret_dict[lev_ell_weight] = ten_decom[wt] * parity

        if self.store_fusion:
            self._fusion_dict[(rep1.high_weight, rep2.high_weight, ell)] = ret_dict

        return ret_dict

    def multi_fusion(self, wts, level):
        '''
        Computes
        :param wts:
        :param level:
        :return:
        '''
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


class TypeALieAlgebra(SimpleLieAlgebra):
    '''
    A type A Lie algebra.

    Inherited attributes:
        rank: a positive integer representing the Lie rank of the algebra.
    '''

    def convert_funds_to_epsilons(self, coords):
        '''
        '''
        # if len(coords) == 1: return list(coords)

        ret_coords = [0]
        part = 0
        for i in reversed(range(len(coords))):
            part += coords[i]
            ret_coords.insert(0, part)

        return ret_coords

    def convert_epsilons_to_funds(self, coords):
        '''
        '''
        # if len(coords) == 2: return [coords[0]- coords[1]]

        ret_coords = []
        for i in range(len(coords) - 1):
            ret_coords.append(coords[i] - coords[i + 1])

        return ret_coords

    def convert_roots_to_funds(self, coords):
        '''
        '''
        if len(coords) == 1: return [2 * coords[0]]

        ret_coords = []
        ret_coords.append(2 * coords[0] - coords[1])
        for i in range(1, len(coords) - 1):
            ret_coords.append(2 * coords[i] - coords[i + 1] - coords[i - 1])

        ret_coords.append(2 * coords[-1] - coords[-2])
        return ret_coords

    def killing_form(self, wt1, wt2):
        '''
        '''
        ret_val = 0

        for i in range(self.rank + 1):
            ret_val = ret_val + wt1.epsilon_coords[i] * wt2.epsilon_coords[i]

        n = sum(wt1.epsilon_coords) * sum(wt2.epsilon_coords)

        ret_val = ret_val - n / (self.rank + 1)

        return ret_val

    def reflect_to_chamber(self, wt):
        '''
        '''
        ret_coords = self._insertsort(wt.epsilon_coords)

        for i in range(len(ret_coords)):
            ret_coords[i] = ret_coords[i] - ret_coords[-1]

        return Weight(self, self.convert_epsilons_to_funds(ret_coords))

    def _insertsort(self, coords):
        ret_list = list(coords)

        for i in range(1, len(ret_list)):
            j = i
            while j > 0 and ret_list[j - 1] < ret_list[j]:
                ret_list[j - 1], ret_list[j] = ret_list[j], ret_list[j - 1]
                j = j - 1

        return ret_list

    def reflect_to_chamber_with_parity(self, wt):
        '''
        '''
        ret_coords, parity = self._insertsort_parity(wt.epsilon_coords)

        for i in range(len(ret_coords)):
            ret_coords[i] = ret_coords[i] - ret_coords[-1]

        return Weight(self, self.convert_epsilons_to_funds(ret_coords)), parity

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
        ret_coords, parity = self._insertsort_parity(wt.epsilon_coords)
        ret_coords = [x - ret_coords[-1] for x in ret_coords]

        while (ret_coords[0] > ell):
            ret_coords[-1] = ret_coords[0] - ell
            ret_coords[0] = ell
            ret_coords, fin_parity = self._insertsort_parity(ret_coords)
            ret_coords = [x - ret_coords[-1] for x in ret_coords]
            parity = parity * -1 * fin_parity

        return Weight(self, self.convert_epsilons_to_funds(ret_coords)), parity

    def get_positive_roots(self):
        '''
        Returns a list of all positive roots of the Lie algebra.
        '''
        if len(self._pos_roots) > 0: return self._pos_roots

        ret_list = []
        coords = []
        for i in range(self.rank):
            coords.append(0)

        for i in range(len(coords)):
            for j in range(i, len(coords)):
                coords[j] = 1
                ret_list.append(Root(self, list(coords)))

            for j in range(i, len(coords)):
                coords[j] = 0

        self._pos_roots = ret_list
        return ret_list

    def get_dual_weight(self, wt):
        return Weight(self, wt[::-1])

    def dual_coxeter(self):
        return self.rank + 1

    def get_level(self, wt):
        '''
        '''
        return wt.epsilon_coords[0]

    def get_orbit_iter(self, wt):
        '''
        Returns an iterable object that iterates through the Weyl group orbit of the given weight.

        Parameters:
            weight: weight whose orbit we want to traverse

        Returns:
            iterator object
        '''
        return self._TypeAOrbitIterator(wt)

    class _TypeAOrbitIterator(object):
        '''
        Iterator object that traverses the Weyl group orbit of a given weight

        Attributes:
            no public attributes
        '''

        def __init__(self, wt):
            self._init_weight = wt.liealg.reflect_to_chamber(wt)
            oms = _OrderedMultiSet(self._init_weight.epsilon_coords)
            self._oms_list = []
            self._index_list = []
            for i in self._init_weight.epsilon_coords:
                self._oms_list.append(oms)
                self._index_list.append(0)
                oms = oms.remove(0)
            self.done = False
            self.liealg = wt.liealg

        def __iter__(self):
            return self

        def next(self):
            if self.done: raise StopIteration()

            new_ep_coords = []
            for i in range(len(self._oms_list)):
                new_ep_coords.append(self._oms_list[i].get(self._index_list[i]))

            # Increment indices
            i = len(self._index_list) - 1
            while i >= 0:
                if self._index_list[i] < self._oms_list[i].num_unique_elements() - 1:
                    break
                i = i - 1
            if i < 0:
                self.done = True
                return Weight(self.liealg, self.liealg.convert_epsilons_to_funds(new_ep_coords))

            self._index_list[i] = self._index_list[i] + 1
            self._oms_list[i + 1] = self._oms_list[i].remove(self._index_list[i])
            for j in range(i + 1, len(self._index_list)):
                if j + 1 < len(self._index_list):
                    self._oms_list[j + 1] = self._oms_list[j].remove(0)
                self._index_list[j] = 0

            return Weight(self.liealg, self.liealg.convert_epsilons_to_funds(new_ep_coords))


class _OrderedMultiSet(object):
    '''
    A simple immutable multiset structure.
    '''

    def __init__(self, initial_list):
        self._multi_set = {}

        for el in initial_list:
            if el in self._multi_set:
                self._multi_set[el] = self._multi_set[el] + 1
            else:
                self._multi_set[el] = 1
        self._key_list = sorted(self._multi_set.keys())

    def num_unique_elements(self):
        return len(self._key_list)

    def get(self, i):
        return self._key_list[i]

    def remove(self, ind):
        ret_list = []
        for el in self._key_list:
            if el == self.get(ind):
                for i in range(self._multi_set[el] - 1):
                    ret_list.append(el)
            else:
                for i in range(self._multi_set[el]):
                    ret_list.append(el)

        return _OrderedMultiSet(ret_list)


class Weight(list):
    '''
    This class represents a weight of a given simple Lie algebra.  The class contains methods that
    are type-independent.

    Attributes:
        liealg: a SimpleLieAlgebra object, the algebra associated to the weight
        epsilon_coords: coordinates of the weight in terms of epsilon coordinates

    Methods:

    '''

    def __init__(self, liealg, coords):
        '''
        Constructor for this class.

        Parameters:
            liealg: a SimpleLieAlgebra object
            coords: coordinates for the weight in terms of the basis of fundamental weights
        '''
        list.__init__(self, coords)
        self.liealg = liealg
        self.epsilon_coords = liealg.convert_funds_to_epsilons(coords)

    def __hash__(self):
        '''
        Hash function for the weight.  Only takes into account the fundamental coordinates, not the
        Lie algebra.
        '''
        return hash(tuple(self))

    def isDominant(self):
        '''
        Checks if the weight is dominant
        '''

        for coord in self:
            if coord < 0:
                return False

        return True


class Root(Weight):
    '''
    This class represents an element of the root lattice of the given simple Lie algebra.

    Attributes:
        root_coords: coordinates in terms of the basis of simple roots of the Lie algebra

    Inherited Attributes:
        liealg: a SimpleLieAlgebra object, the algebra associated to the weight
        fund_coords: coordinates of the weight in terms of the basis of fundamental weights
        epsilon_coords: coordinates of the weight in terms of epsilon coordinates
    '''

    def __init__(self, liealg, coords):
        '''
        Constructor for this class.

        Parameters:
            liealg: a SimpleLieAlgebra object
            coords: coordinates for the weight in terms of the basis of simple roots
        '''
        Weight.__init__(self, liealg, liealg.convert_roots_to_funds(coords))
        self.root_coords = coords
        self.epsilon_coords = liealg.convert_funds_to_epsilons(self)

    def get_root_level(self):
        ret_val = 0
        for coord in self.root_coords:
            ret_val = ret_val + coord

        return ret_val