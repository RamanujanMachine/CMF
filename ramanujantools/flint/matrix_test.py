import sympy as sp
from sympy.abc import n

from ramanujantools import Matrix
from ramanujantools.flint import FlintMatrix


def flintify(matrix: Matrix) -> FlintMatrix:
    return FlintMatrix.from_sympy(matrix)


def test_factor():
    matrix = Matrix(
        [
            [1, n**2 + n, n**2 - n + 5],
            [3 * n + 9, n**2 - 1, 1 / (n + 1)],
            [n**2 + 7, n - 2, (n + 1) * (n - 3)],
        ],
    )

    assert matrix.applyfunc(sp.factor) == flintify(matrix).factor()


def test_mul():
    m1 = Matrix([[0, n**2], [1, 1 / n]])
    m2 = Matrix([[3, n - 2], [n, 5]])

    assert flintify(m1 * m2) == flintify(m1) * flintify(m2)
    assert flintify(m2 * m1) == flintify(m2) * flintify(m1)


def test_mul_scalar():
    m = Matrix([[0, n**2], [1, 1 / n]])

    assert flintify(m * 17) == flintify(m) * 17
    assert flintify(m * 17) == 17 * flintify(m)


def test_walk():
    matrix = Matrix(
        [
            [1, n**2 + n, n**2 - n + 5],
            [3 * n + 9, n**2 - 1, 1 / (n + 1)],
            [n**2 + 7, n - 2, (n + 1) * (n - 3)],
        ],
    )

    expected = (matrix * matrix.subs({n: n + 1}) * matrix.subs({n: n + 2})).factor()

    assert expected == flintify(matrix).walk({n: 1}, 3, {n: n}).factor()
