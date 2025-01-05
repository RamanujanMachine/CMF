from __future__ import annotations
from typing import Union, List, Dict, Set, Tuple
from functools import cached_property
import copy
import itertools
from tqdm import tqdm

import sympy as sp
from sympy.abc import n

from ramanujantools import Matrix, Limit, GenericPolynomial


def trim_trailing_zeros(sequence: List[int]) -> List[int]:
    ending = len(sequence)
    while sequence[ending - 1] == 0:
        ending -= 1
    return sequence[0:ending]


class LinearRecurrence:
    r"""
    Represents a linear recurrence of the form
    $a_0(n) p(n + 1) = \sum_{i=1}^{N}a_i(n) p(n + 1 - i)$

    Note that the beginning index can be decided later,
    therefore this class actually represents all recurrences of the form
    $a_0(n + s) p(n + 1) = \sum_{i=1}^{N}a_i(n + s) p(n + 1 - i)$
    for all integers $s$.
    """

    def __init__(self, recurrence: Union[Matrix, List[sp.Expr]]):
        r"""
        Construct the recurrence.

        The recurrence argument can be one of two types:
            1. A list of the coefficients of the recurrence [a_0(n), ..., a_N(n)]
            2. A matrix which is companionized and used as the recurrence sequence
        """
        if isinstance(recurrence, Matrix):
            recurrence_matrix = recurrence.as_companion()
            col = recurrence_matrix.col(-1)
            lead = col.denominator_lcm
            coeffs = [sp.simplify(p * lead) for p in reversed(col)]
            relation = [lead] + coeffs
        else:
            relation = recurrence
        if len(relation) == 0:
            raise ValueError("Attempted to construct an empty recurrence!")
        relation = [sp.factor(sp.simplify(p)) for p in relation]
        self.relation = trim_trailing_zeros(relation)

    def __eq__(self, other: Matrix) -> bool:
        """
        Returns True iff two requrences are identical (even up to gcd).
        """
        return self.relation == other.relation

    def __mul__(self, scalar: int) -> LinearRecurrence:
        return LinearRecurrence([p * scalar for p in self.relation])

    def __rmul__(self, scalar: int) -> LinearRecurrence:
        return self * scalar

    def __truediv__(self, scalar: int) -> LinearRecurrence:
        return LinearRecurrence([p / scalar for p in self.relation])

    def __floordiv__(self, scalar: int) -> LinearRecurrence:
        return LinearRecurrence([p / scalar for p in self.relation])

    def __repr__(self) -> str:
        return f"LinearRecurrence({self.relation})"

    def __str__(self) -> str:
        p = sp.Function("p")(n)
        lhs = f"{(self.relation[0])*p.subs({n: n+1})}"
        rhs = " + ".join(
            [f"{c*p.subs({n: n-i})}" for i, c in enumerate(self.relation[1:])]
        )
        return f"{lhs} = {rhs}"

    @cached_property
    def gcd(self) -> sp.Expr:
        """
        Returns the GCD of all recurrence coefficients
        """
        return sp.gcd([r.as_numer_denom()[0] for r in self.relation])

    @cached_property
    def denominator_lcm(self) -> sp.Expr:
        return sp.lcm([r.as_numer_denom()[1] for r in self.relation])

    @cached_property
    def recurrence_matrix(self) -> Matrix:
        """
        Returns the companion form recurrence matrix corresponding to the recurrence
        """
        denominator = sp.simplify(self.relation[0])
        column = [c / denominator for c in self.relation[1:]]
        return Matrix.companion_form(list(reversed(column)))

    def depth(self) -> int:
        """
        Returns the depth of the recurrence
        """
        return len(self.relation) - 1

    def degrees(self) -> List[int]:
        """
        Returns a list of the degrees of all coefficients
        """
        return [sp.Poly(p, n).degree() for p in self.relation]

    def subs(self, substitutions: Dict) -> LinearRecurrence:
        """
        Substitutes symbols in the recurrence.
        """
        return LinearRecurrence([p.subs(substitutions) for p in self.relation])

    def free_symbols(self) -> Set[sp.Symbol]:
        """
        Returns all free symbols of the recurrence (including `n`)
        """
        return set.union(*[p.free_symbols for p in self.relation])

    def parameters(self) -> Set[sp.Symbol]:
        """
        Returns all symbolic parameters of the recurrence (excluding `n`)
        """
        return self.free_symbols() - self.axes()

    def normalize(self) -> LinearRecurrence:
        """
        Normalizes the recurrence
        """
        relation = [p * self.denominator_lcm / self.gcd for p in self.relation]
        return LinearRecurrence(relation).simplify()

    def simplify(self) -> LinearRecurrence:
        """
        Simplifies the coefficients of the recurrence
        """
        relation = [p * self.denominator_lcm / self.gcd for p in self.relation]
        return LinearRecurrence([sp.factor(p.simplify()) for p in relation])

    def limit(self, iterations: Union[int, List[int]], start=1) -> Limit:
        r"""
        Returns the Limit matrix of the recursion up to a certain depth
        """
        return self.recurrence_matrix.limit({n: 1}, iterations, {n: start})

    def compose(self, composition: sp.Expr) -> LinearRecurrence:
        r"""
        Composes the recurrence into a higher depth recurrence.

        Given a recurrence
        $a_0(n) p(n + 1) = \sum_{i=1}^{N}a_i(n) p(n + 1 - i)$

        We can derive that
        $a_0(n) p(n) = \sum_{i=1}^{N}a_i(n - 1) p(n - i)$

        Selecting a composition polynomial $d$, we can subtract the latter from the recursion d(n) times.

        Example:
            >>> s = LinearRecurrence([sp.Function("a")(n), sp.Function("b")(n), sp.Function("c")(n)])
            >>> s
            LinearRecurrence([a(n), b(n), c(n)])
            >>> s.compose(sp.Function("d")(n))
            LinearRecurrence([a(n), -a(n - 1)*d(n) + b(n), b(n - 1)*d(n) + c(n), c(n - 1)*d(n)])
        """
        relation = self.relation
        modification = (
            [0]
            + [-composition * self.relation[0].subs({n: n - 1})]
            + [composition * c.subs({n: n - 1}) for c in self.relation[1:]]
        )
        relation = [
            sum(d) for d in itertools.zip_longest(relation, modification, fillvalue=0)
        ]
        return LinearRecurrence(relation)

    @staticmethod
    def all_divisors(p: sp.Poly) -> List[sp.Poly]:
        r"""
        Returns all divisors of polynomial `p`.
        Assumes p is polynomial in `n`.
        """
        p = sp.Poly(p, n)
        content, factors_list = p.factor_list()
        factors = []
        for factor, order in factors_list:
            factors.append([factor**d for d in range(order + 1)])
        if len(content.free_symbols) == 0:
            for root, order in content.factors().items():
                factors.append([root**d for d in range(order + 1)])
        combinations = itertools.product(*factors)
        divisors = []
        for combination in combinations:
            divisors.append(sp.prod(combination))
        return divisors

    def possible_decompositions(self) -> List[sp.Poly]:
        r"""
        Returns all candidates for a composition polynomial $d(n)$
        that could have been used to compose a lesser depth recursion into this one.
        """
        return LinearRecurrence.all_divisors(self.relation[-1])

    def decompose_poly(self, decomposition: sp.Poly) -> LinearRecurrence:
        r"""
        Attempts to decompose this recursion using a decomposition polynomial.
        In case of success, returns a recurrence such that `recurrence.compose(decomposition) == self`

        If `self` contains parameters (other than n), will attempt to find a matching substitution.
        In case of success, the returned recurrence will be substituted with the solution.
        i.e, for that solution, `recurrence.compose(decomposition) == self.subs(solution)`
        """
        decomposition = sp.Poly(decomposition, n).as_expr()
        decomposed = [self.relation[0]]
        for index in range(1, len(self.relation) - 1):
            next = (-1 if index == 1 else 1) * decomposed[index - 1].subs({n: n - 1})
            decomposed.append(sp.simplify(self.relation[index] - next * decomposition))
        expected_tail = decomposed[-1].subs({n: n - 1}) * decomposition
        if sp.simplify(expected_tail - self.relation[-1]) == 0:
            return LinearRecurrence(decomposed)
        if len(self.free_symbols()) > 0:
            solutions = sp.solve(sp.Poly(expected_tail - self.relation[-1], n).coeffs())
            assert len(solutions) < 2
            for solution in solutions:
                decomposed = [d.subs(solution) for d in decomposed]
                return LinearRecurrence(decomposed)
        return None

    def decompose(self, inflation_degree=0) -> Tuple[LinearRecurrence, sp.Poly]:
        r"""
        Attempts to decompose this recursion by enumerating over all possible
        decomposition candidates and attempting to decompose using them.

        If `inflation_degree` is not zero, will inflate `self` by a generic polynomial of that degree.
        The solver will then attempt to find a solution for that inflation that manages to decompose the recurrence.
        """
        results = []
        if inflation_degree > 0:
            inflation, _ = GenericPolynomial.of_degree(inflation_degree, "c", n)
            recurrence = self.inflate(inflation)
        elif inflation_degree < 0:
            inflation, _ = GenericPolynomial.of_degree(-inflation_degree, "c", n)
            recurrence = self.deflate(inflation)
        else:
            recurrence = self
        for decomposition in tqdm(self.possible_decompositions()):
            decomposed = recurrence.decompose_poly(decomposition)
            if decomposed is not None:
                results.append((decomposed, decomposition))
        return results

    def inflate(self, c: sp.Expr) -> LinearRecurrence:
        r"""
        Inflates the recurrence by a polynomial c.

        The inflated recurrence satisfies
        $a_0(n) p(n + 1) = \sum_{i=1}^{N}(\prod_{j=0}^{i-1}c(n-i))a_i(n) p(n + 1 - i)$

        The inflated recurrence converges to the same limit up to different initial values.
        """
        c = sp.simplify(c).as_expr()
        current = c
        relation = copy.deepcopy(self.relation)
        for i in range(1, len(self.relation)):
            relation[i] *= current
            current *= c.subs({n: n - i})
        return LinearRecurrence(relation)

    def deflate(self, c: sp.Expr) -> LinearRecurrence:
        r"""
        Deflates the recurrence by a polynomial c.

        Equivalent to `self.inflate(1 / c)`.
        """
        recurrence = self.inflate(1 / c)
        return recurrence

    def apteshuffle(self) -> LinearRecurrence:
        if self.depth() != 3:
            raise ValueError(
                "apteshuffle is only supported for depth 3 recursions for now"
            )
        c, b, a = self.recurrence_matrix.col(-1)
        return (
            LinearRecurrence([1, -b, -c * a.subs({n: n - 1}), c * c.subs({n: n - 1})])
            .simplify()
            .subs({n: n + 1})
        )

    def eigenvals(self):
        return self.recurrence_matrix.eigenvals(poincare=True)

    def errors(self):
        return self.recurrence_matrix.errors()

    def gcd_slope(self, depth=20):
        return self.recurrence_matrix.gcd_slope(depth)

    def kamidelta(self, depth=20):
        return self.recurrence_matrix.kamidelta(depth)
