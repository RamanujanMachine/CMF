from pytest import approx
from sympy.abc import c, n
import mpmath as mp

from ramanujantools import Matrix
from ramanujantools.pcf import PCF
from ramanujantools.pcf.pcf import is_deflatable


def test_repr():
    pcf = PCF(1 + n, 3 - n)
    assert pcf == eval(repr(pcf))


def test_degree():
    pcf = PCF(1 + n - n**2, 3 - n**9)
    assert (2, 9) == pcf.degree()


def test_singular_points():
    a_n = n + 1
    b_n = (n + 1) * (n - 17) * (n + 59) * (n - 102)
    pcf = PCF(a_n, b_n)
    assert pcf.singular_points() == [{n: 17}, {n: 102}]


def test_limit_as_float():
    pcf = PCF(5 + 10 * n, 1 - 9 * n**2)
    expected = (4 ** (1 / 3) + 1) / (4 ** (1 / 3) - 1)
    assert expected == approx(pcf.limit(100).as_float(), 1e-4)


def test_walk_list():
    iterations = [1, 2, 3, 17, 29, 53, 99]
    pcf = PCF(5 + 10 * n, 1 - 9 * n**2)
    assert pcf.walk(iterations) == [pcf.walk(i) for i in iterations]


def test_walk_start():
    iterations = [1, 2, 3, 17, 29, 53, 99]
    p = PCF(n + 7, 3 * n**2 - 1)
    expected = p.walk(sum(iterations))
    actual = Matrix.eye(2)
    for i in range(len(iterations)):
        actual *= p.walk(iterations[i], start=sum(iterations[0:i]))
    assert expected == actual


def test_inflate_constant():
    c = 3
    a_n = n + 4
    b_n = n**2
    pcf = PCF(a_n, b_n)
    assert PCF(c * a_n, c**2 * b_n) == pcf.inflate(c)


def test_inflate_symbol():
    a_n = n + 4
    b_n = n**2
    pcf = PCF(a_n, b_n)
    assert PCF(c * a_n, c**2 * b_n) == pcf.inflate(c)


def test_inflate_poly():
    a_n = n + 4
    b_n = n**2
    c_n = n**7 + 5 * n - 3
    pcf = PCF(a_n, b_n)
    assert PCF(c_n * a_n, c_n.subs({n: n - 1}) * c_n * b_n) == pcf.inflate(c_n)


def test_inflate_matrix_equivalent():
    a_n = n + 4
    b_n = n**2
    c_n = n**7 + 5 * n - 3
    pcf = PCF(a_n, b_n)
    assert pcf.M().inflate(c_n) == pcf.inflate(c_n).M()


def test_deflate_constant():
    a_n = n + 4
    b_n = n**2
    pcf = PCF(c * a_n, c**2 * b_n)
    assert PCF(a_n, b_n) == pcf.deflate(c)


def test_deflate_poly():
    a_n = n + 4
    b_n = n**2
    c_n = n**7 + 5 * n - 3
    pcf = PCF(c_n * a_n, c_n.subs({n: n - 1}) * c_n * b_n)
    assert PCF(a_n, b_n) == pcf.deflate(c_n)


def test_deflate_matrix_equivalent():
    a_n = n + 4
    b_n = n**2
    c_n = n**7 + 5 * n - 3
    pcf = PCF(c_n * a_n, c_n.subs(n, n - 1) * c_n * b_n)
    assert pcf.M().deflate(c_n) == pcf.deflate(c_n).M()


def test_is_deflatable():
    a_factors = {n**2 + 1: 1}
    b_factors = {n**2 - 2*n + 2: 1, n**2 + 1: 1}
    assert is_deflatable(a_factors, b_factors, n**2 + 1)


def test_deflate_all():
    c_n = c**2 * (7 * n - 13 * c) * (2 * n - 5) ** 4 * (3 * n + 11)
    a_n = n + c
    b_n = n**2 - c * n
    pcf1 = PCF(c_n * a_n, c_n.subs({n: n - 1}) * c_n * b_n)
    pcf2 = PCF(n**2 + 1, n**4 - 2*n**3 + 3*n**2 - 2*n + 2)
    pcf3 = PCF(9*n**4 + 72*n**3 + 201*n**2 + 228*n + 85,
              -18*n**8 - 225*n**7 - 1131*n**6 - 2904*n**5 \
                - 3932*n**4 - 2400*n**3 + 66*n**2 + 769*n + 255)
    assert PCF(a_n, b_n) == pcf1.deflate_all()
    assert PCF(1, 1) == pcf2.deflate_all()
    assert PCF(3*n**2 + 9*n + 5, -2*n**4 - 9*n**3 - 9*n**2 + n + 3) \
        == pcf3.deflate_all()


def test_blind_delta():
    pcf = PCF(34 * n**3 + 51 * n**2 + 27 * n + 5, -(n**6))
    depth = 2000
    delta = pcf.delta(depth)
    assert delta > 0.08


def test_precision_e():
    pcf = PCF(n, n)
    assert pcf.limit(2**10 + 1).precision() == 2642


def test_precision_phi():
    pcf = PCF(1, 1)
    assert pcf.limit(2**10 + 1).precision() == 427


def test_delta_sequence_agrees_with_delta():
    pcf = PCF(2 * n + 1, n**2)
    depth = 50
    limit = 4 / mp.pi

    actual_deltas = pcf.delta_sequence(depth, limit)
    expected_deltas = []
    for dep in range(1, depth + 1):
        expected_deltas.append(pcf.delta(dep, limit))

    assert expected_deltas == actual_deltas


def test_blind_delta_sequence_agrees_with_blind_delta():
    pcf = PCF(2 * n + 1, n**2)
    depth = 50
    limit = pcf.limit(2 * depth).as_float()

    actual_values = pcf.delta_sequence(depth)
    expected_deltas = []
    for dep in range(1, depth + 1):
        expected_deltas.append(pcf.delta(dep, limit))

    assert expected_deltas == actual_values
