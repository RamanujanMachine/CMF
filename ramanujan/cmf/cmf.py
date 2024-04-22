from __future__ import annotations
import itertools
from multimethod import multimethod
from typing import Collection, Dict, List, Set, Union

import sympy as sp
from sympy.abc import n

from ramanujan import Matrix, simplify, zero
from ramanujan.pcf import PCFFromMatrix


class CMF:
    r"""
    Represents a Conservative Matrix Field (CMF).

    A CMF is defined by a set of axes and their relevant matrices that satisfy the perservation quality:
    for every two axes x and y, $Mx(x, y) \cdot My(x+1, y) = My(x, y) \cdot Mx(x, y+1)$
    """

    def __init__(self, matrices: Dict[sp.Symbol, Matrix]):
        """
        Initializes a CMF with `Mx` and `My` matrices
        """

        self.matrices = matrices

        assert (
            n not in self.matrices.keys()
        ), "Do not use symbol n as an axis, it's reserved for PCF conversions"

        for x, y in itertools.combinations(self.matrices.keys(), 2):
            Mx = self.matrices[x]
            My = self.matrices[y]
            Mxy = simplify(Mx * My.subs({x: x + 1}))
            Myx = simplify(My * Mx.subs({y: y + 1}))
            if simplify(Mxy - Myx) != Matrix([[0, 0], [0, 0]]):
                raise ValueError(f"M{x} and M{y} matrices are not conserving!")

    def __eq__(self, other) -> bool:
        return self.matrices == other.matrices

    def __repr__(self) -> str:
        return f"CMF({self.matrices})"

    def M(self, axis: sp.Symbol) -> Matrix:
        """
        Returns the axis matrix for a given axis.
        """
        return self.matrices[axis]

    def axes(self) -> Set[sp.Symbol]:
        """
        Returns the symbols of all axes of the CMF.
        """
        return set(self.matrices.keys())

    def axis_vector(self, axis: sp.Symbol) -> Dict[sp.Symbol, int]:
        """
        Given a CMF axis symbol `axis`,
        Returns the vector which is 1 in that axis and 0 in all other axes.
        """
        return {i: int(i == axis) for i in self.axes()}

    def parameters(self) -> Set[sp.Symbol]:
        """
        Returns all (non-axis) symbolic parameters of the CMF.
        """
        return self.free_symbols() - self.axes()

    def free_symbols(self) -> Set[sp.Symbol]:
        """
        Returns all symbolic variables of the CMF, both axes and parameters.
        """
        return set.union(
            *list(map(lambda matrix: matrix.free_symbols, self.matrices.values()))
        )

    def subs(self, *args, **kwargs) -> CMF:
        """Returns a new CMF with substituted Mx and My."""
        return CMF(
            matrices={
                symbol: matrix.subs(*args, **kwargs)
                for symbol, matrix in self.matrices.items()
            }
        )

    def simplify(self):
        """Returns a new CMF with simplified Mx and My"""
        return CMF(
            matrices={symbol: simplify(matrix) for symbol, matrix in self.matrices}
        )

    def default_origin(self):
        """
        Returns the default origin value, which is 1 for every axis.
        """
        return {axis: 1 for axis in self.axes()}

    def trajectory_matrix(self, trajectory: dict, start: dict = None) -> Matrix:
        """

        Returns the corresponding matrix for walking in trajectory.
        If `start` is given, the new matrix will be reduced to a single variable `n`.
        Args:
            trajectory: a dict containing the amount of steps in each direction.
            start: a dict representing the starting point of the multiplication.
        Returns:
            A matrix that represents a single step in the desired trajectory
        """
        assert (
            self.axes() == trajectory.keys()
        ), f"Trajectory axes {trajectory.keys()} does not match CMF axes {self.axes()}"

        if start:
            assert (
                self.axes() == start.keys()
            ), f"Start axes {start.keys()} does not match CMF axes {self.axes()}"

        position = {axis: axis for axis in self.axes()}
        m = sp.eye(2)
        for axis in self.axes():
            print(axis)
            m *= self.M(axis).walk(self.axis_vector(axis), trajectory[axis], position)
            position[axis] += trajectory[axis]
        if start:
            m = CMF.substitute_trajectory(m, trajectory, start)
        return simplify(m)

    def as_pcf(self, trajectory) -> PCFFromMatrix:
        """
        Returns the PCF equivalent to the CMF in a certain trajectory, up to a mobius transform.
        If `start` is None, will use the default origin value.
        """
        return self.trajectory_matrix(trajectory, self.default_origin()).as_pcf()

    @staticmethod
    def substitute_trajectory(
        trajectory_matrix: Matrix, trajectory: dict, start: dict
    ) -> Matrix:
        """
        Returns trajectory_matrix reduced to a single variable `n`.

        This transformation is possible only when the starting point is known.
        Each incrementation of the variable `n` represents a full step in `trajectory`.
        """
        from sympy.abc import n

        assert (
            trajectory.keys() == start.keys()
        ), f"Key mismatch between trajectory ({trajectory.keys()}) and start ({start.keys()})"

        def sub(i):
            return start[i] + (n - 1) * trajectory[i]

        return trajectory_matrix.subs([(axis, sub(axis)) for axis in trajectory.keys()])

    @multimethod
    def walk(
        self,
        trajectory: dict,
        iterations: Collection[int],
        start: Union[dict, type(None)] = None,
    ) -> List[Matrix]:
        r"""
        Returns a list of trajectorial walk multiplication matrices in the desired depths.

        The walk operation is defined as $\prod_{i=0}^{n-1}M(s_0 + i \cdot t_0, ..., s_k + i \cdot t_k)$,
        where `M=trajectory_matrix(trajectory, start)`, and `n / size(trajectory)` (L1 size - total amount of steps)

        Args:
            trajectory: a dict containing the amount of steps in each direction.
            iterations: the amount of multiplications to perform, either an integer or a list.
            start: a dict representing the starting point of the multiplication, `default_origin` by default.
        Returns:
            The walk multiplication as defined above.
            If `iterations` is a list, returns a list of matrices.
        """
        trajectory_matrix = self.trajectory_matrix(
            trajectory, start or self.default_origin()
        )
        actual_iterations = [i // sum(trajectory.values()) for i in iterations]
        return trajectory_matrix.walk({n: 1}, actual_iterations, {n: 1})

    @multimethod
    def walk(  # noqa: F811
        self,
        trajectory: dict,
        iterations: int,
        start: Union[dict, type(None)] = None,
    ) -> Matrix:
        return self.walk(trajectory, [iterations], start)[0]

    @multimethod
    def limit(
        self,
        trajectory: dict,
        iterations: Collection[int],
        start: Union[dict, type(None)] = None,
        vector: Matrix = zero(),
    ) -> List[Matrix]:
        """
        Returns the convergence limit of walking in a certain trajectory.

        This is essentially the same as `self.walk(trajectory, iterations, start) * vector`

        Args:
            trajectory: a dict containing the amount of steps in each direction.
            iterations: the amount of multiplications to perform, either an integer or a list.
            start: a dict representing the starting point of the multiplication, `default_origin` by default.
            vector: The final vector to multiply the matrix by (the zero vector by default)
        Returns:
            The limit of the walk multiplication as defined above.
            If `iterations` is a list, returns a list of matrices.
        """
        m_walk = self.walk(trajectory, iterations, start)
        return [mat * vector for mat in m_walk]

    @multimethod
    def limit(  # noqa: F811
        self,
        trajectory: dict,
        iterations: int,
        start: Union[dict, type(None)] = None,
        vector: Matrix = zero(),
    ) -> Matrix:
        return self.limit(trajectory, [iterations], start, vector)[0]
