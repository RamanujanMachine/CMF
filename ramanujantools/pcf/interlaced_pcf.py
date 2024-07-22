import sympy as sp
from sympy.abc import n

from typing import Dict, List, Collection
from multimethod import multimethod

from ramanujantools import Matrix, Limit
from ramanujantools.pcf import PCF


class InterlacedPCF(PCF):
    r"""
    Represents Interlaced Polynomial Continued Fractions. These are PCFs constructed by partial 
    numerators and denominators that alternate between a few different polynomials.
        
    The polynomial matrix representation $M_n$ is given by the product of the 
    matrices representing each numerator-denominator pair.
    For example, consider a continued fraction with three sequences appearing as the partial 
    numerator - $b_1, b_2, b_3$ - and two sequences appearing as the partial denominator 
    - $a_1, a_2$:
    The CF's matrix is represented by $M_n = \prod_{i=1}^6 \begin{pmatrix} 0 & b_i(n) \cr 1 & a_i(n) \end{pmatrix}$.
    Where $6 = lcm(3, 2)$ is the period of the continued fraction and the sequences are repeated 
    to match the period of the continued fraction.
    """

    def __init__(self, a0, a_dict, b_dict):
        """
        Initializes an Interlaced PCF with
        - a0: The free constant of the continued fraction.
        - a_dict: A dictionary of polynomials representing `a_n`, the partial denominator of the continued fraction.
        - b_dict: A dictionary of polynomials representing `b_n`, the partial numerator of the continued fraction.
        Note: the order of the polynomials in the dictionaries is the order in which they will appear in the continued fraction.
        """
        self.a0 = a0
        self.a_dict = {f'a_{i+1}': val for (i, val) in enumerate(list(a_dict.values()))}
        self.b_dict = {f'b_{i+1}': val for (i, val) in enumerate(list(b_dict.values()))}
        self.period = sp.lcm(len(a_dict), len(b_dict))

        self.a_list = list(a_dict.values()) * (self.period // len(a_dict))
        self.b_list = list(b_dict.values()) * (self.period // len(b_dict))

        M = Matrix.eye(2)
        for a, b in zip(self.a_list, self.b_list):
            M = M * Matrix([[0, b], [1, a]])
        self.M_mat = M

    def __repr__(self):
        return f"InterlacedPCF({self.a0}, {self.a_dict}, {self.b_dict})"
    
    def __str__(self):
        return self.__repr__()
    
    def M(self):
        return self.M_mat
    
    def A(self):
        return Matrix([[1, self.a0], [0, 1]])
    
    def pcf(self):
        return self.M_mat.as_pcf().pcf

