from sympy.abc import n

from ramanujantools.pcf import PCF, InterlacedPCF
from ramanujantools import Matrix


def test_period():
    a0 = 1
    a_dict = {1: 1, 2: 1, 3: 1}
    b_dict = {1: 1, 2: 1}
    pcf = InterlacedPCF(a0, a_dict, b_dict)
    assert pcf.period == 6


def test_a_list():
    a0 = 1
    a_dict = {1: 1, 2: 1, 3: 1}
    b_dict = {1: 1, 2: 1}
    pcf = InterlacedPCF(a0, a_dict, b_dict)
    assert pcf.a_list == [1, 1, 1, 1, 1, 1]


def test_M():
    m = InterlacedPCF(3, {1: 3}, {1: -n, 2: n}).M()
    assert Matrix([[-n, -3 * n], [3, n + 9]]) == m


def test_pcf():
    pcf = InterlacedPCF(3, {1: 3}, {1: -n, 2: n}).pcf()
    assert PCF(8, n**2) == pcf