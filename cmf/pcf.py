import sympy as sp
from sympy.abc import n

from cmf import Matrix


def is_pcf(M: Matrix):
    return M[0, 0] == 1 and M[1, 0] == 0


class PCF:
    @staticmethod
    def A_matrix(M: Matrix):
        a_0 = M[1, 1].subs(n, 0)
        return Matrix([[1, a_0], [0, 1]])

    def __init__(self, M: Matrix):
        self.m_M = M
        if not is_pcf(self.M):
            raise ValueError("The given matrix M is not of a PCF form!")
        self.m_A = PCF.A_matrix(self.m_M)

    def __init__(self, an, bn):
        self.m_M = Matrix([[0, bn], [1, an]])
        self.m_A = PCF.A_matrix(self.m_M)

    def subs(self, substitutions):
        return PCF(self.m_M.subs(substitutions))

    def walk(self, iterations, start=1) -> Matrix:
        return self.m_M.walk([1], iterations, start)

    def limit(self, depth, start=[1], vector=Matrix([[0], [1]])) -> sp.Float:
        return (self.m_A * self.walk(depth, start)).limit(vector)
