import pytest
from sympy.abc import x, y, n

from ramanujan import GenericPolynomial
from ramanujan.pcf import PCF
from ramanujan.cmf.ffbar import solve_ffbar, from_pcf
from .solver import solve

from ramanujan.cmf.known_cmfs import cmf1


@pytest.mark.parametrize("deg", [1, 2])
def test_solver_full_poly(deg: int):
    f, _ = GenericPolynomial.of_combined_degree(deg=deg, var_name="c", variables=[x, y])
    fbar, _ = GenericPolynomial.of_combined_degree(
        deg=deg, var_name="d", variables=[x, y]
    )
    assert len(solve_ffbar(f, fbar)) > 0


def test_from_ffbar_gauss():
    pcf = PCF(2 * n + 1, -(n**2))
    solutions = from_pcf(pcf)
    assert len(solutions) == 1
    cmf = solutions[0]
    assert pcf.M() == cmf.M(x).subs(x, n)


def test_ffbar_from_cmf1():
    cmf = cmf1()
    pcf = cmf.as_pcf({x: 1, y: 0}).pcf
    solutions = from_pcf(pcf)
    for ffbar in solutions:
        assert (
            len(solve([ffbar.M(x) - cmf.M(x), ffbar.M(y) - cmf.M(y)])) > 0
        )  # Making sure that our CMF is equivalent up to symbol names
