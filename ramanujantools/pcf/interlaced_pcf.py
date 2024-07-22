import sympy as sp

from ramanujantools import Matrix
from ramanujantools.pcf import PCF


class InterlacedPCF(PCF):
    r"""
    Represents Interlaced Polynomial Continued Fractions. These are PCFs constructed by partial
    numerators and denominators that alternate between a few different polynomials.
    An Interlaced PCF is not a PCF, but can be converted into one.

    The polynomial matrix representation $M_n$ is given by the product of the
    matrices representing each numerator-denominator pair.
    For example, consider a continued fraction with three polynomials appearing as the partial
    numerator - $b_1, b_2, b_3$ - and two polynomials appearing as the partial denominator - $a_1, a_2$:
    The Interlaced PCF is represented by
    $M_n = \prod_{i=1}^6 \begin{pmatrix} 0 & b_i(n) \cr 1 & a_i(n) \end{pmatrix}$,
    where $6 = lcm(3, 2)$ is the period of the continued fraction and the sequences of polynomials are repeated
    to match the period of the continued fraction.
    """

    def __init__(self, a0, a_dict, b_dict):
        """
        Initializes an Interlaced PCF with a number `a0` and `a_dict`, and `b_dict`, dictionaries of polynomials.
        For example: `inter_pcf = InterlacedPCF(3, {1: 3}, {1: -n, 2: n})`
        Note: the order of the polynomials in the dictionaries is the order
        in which they will appear in the continued fraction.
        """
        self.a0 = a0
        """The free constant of the continued fraction."""

        self.a_dict = {f'a_{i+1}': val for (i, val) in enumerate(list(a_dict.values()))}
        """A dictionary of polynomials representing `a_n`, the partial denominator of the continued fraction."""

        self.b_dict = {f'b_{i+1}': val for (i, val) in enumerate(list(b_dict.values()))}
        """A dictionary of polynomials representing `b_n`, the partial numerator of the continued fraction."""

        self.period = sp.lcm(len(a_dict), len(b_dict))
        """The lcm of the lengths of `a_dict` and `b_dict`."""

        self.a_list = list(a_dict.values()) * (self.period // len(a_dict))
        """The partial denominators of the continued fraction repeated to match the period."""

        self.b_list = list(b_dict.values()) * (self.period // len(b_dict))
        """The partial numerators of the continued fraction repeated to match the period."""

        M = Matrix.eye(2)
        for a, b in zip(self.a_list, self.b_list):
            M = M * Matrix([[0, b], [1, a]])

        self.M_mat = M
        """The polynomial matrix representing the interlaced continued fraction."""

    def __repr__(self):
        return f"InterlacedPCF({self.a0}, {self.a_dict}, {self.b_dict})"

    def __str__(self):
        return self.__repr__()

    def M(self):
        r"""
        Returns the matrix representation of the interlaced continued fraction.

        $M = \prod\limits_{i=1}^p \begin{pmatrix} 0, b_i(n) \cr 1, a_i(n) \end{pmatrix}$
        where $p$ is the period of the Interlaced PCF.
        """
        return self.M_mat

    def A(self):
        r"""
        Returns the matrix that represents the $a_0$ part:

        $A = \begin{pmatrix} 1, a_0 \cr 0, 1 \end{pmatrix}$
        """
        return Matrix([[1, self.a0], [0, 1]])

    def pcf(self):
        """
        Returns the PCF representation of the interlaced continued fraction.
        Note: an Interlaced PCF is not a PCF, but can be converted into one.
        """
        return self.M_mat.as_pcf().pcf
