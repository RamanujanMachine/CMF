from sympy.abc import n

from ramanujantools.pcf import PCF, InterlacedPCF
from ramanujantools import Matrix


def test_period():
    a0 = 1
    a_dict = {'a1': 1, 'a2': -n, 'a3': 3 * n**3 - 5}
    b_dict = {'b1': n, 'b2': 2 * n**2}
    pcf = InterlacedPCF(a0, a_dict, b_dict)
    assert 6 == pcf.period


def test_ab_lists():
    a0 = 1
    a_dict = {'a1': 1, 'a2': -n, 'a3': 3 * n**3 - 5}
    b_dict = {'b1': n, 'b2': 2 * n**2}
    pcf = InterlacedPCF(a0, a_dict, b_dict)
    assert [1, -n, 3 * n**3 - 5, 1, -n, 3 * n**3 - 5] == pcf.a_list
    assert [n, 2 * n**2, n, 2 * n**2, n, 2 * n**2] == pcf.b_list


def test_M():
    pcf = InterlacedPCF(3, {1: 3}, {1: -n, 2: n})
    assert Matrix([[-n, -3 * n], [3, n + 9]]) == pcf.M()


def test_pcf():
    pcf = InterlacedPCF(3, {1: 3}, {1: -n, 2: n}).pcf()
    assert PCF(8, n**2) == pcf