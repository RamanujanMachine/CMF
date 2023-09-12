import sympy as sp
from sympy.abc import n, x, y

from ramanujan import Matrix, Vector, simplify
from ramanujan.pcf import PCFFromMatrix


class CMF:
    r"""
    Represents a Conservative Matrix Field (CMF).

    A CMF is defined by two matrices $Mx, My$ that satisfy the perservation quality:
    $Mx(x, y) \cdot My(x+1, y) = My(x, y) \cdot Mx(x, y+1)$
    """

    def __init__(self, Mx: Matrix, My: Matrix):
        """
        Initializes a CMF with `Mx` and `My` matrices
        """

        self.Mx = Mx
        """The Mx matrix of the CMF"""
        self.My = My
        """The My matrix of the CMF"""

        Mxy = simplify(self.Mx @ self.My({x: x + 1}))
        Myx = simplify(self.My @ self.Mx({y: y + 1}))
        if simplify(Mxy - Myx) != Matrix([[0, 0], [0, 0]]):
            raise ValueError("The given Mx and My matrices are not conserving!")

    def __repr__(self):
        return f"CMF({self.Mx}, {self.My})"

    def subs(self, *args, **kwrags):
        """Returns a new CMF with substituted Mx and My."""
        return CMF(self.Mx.subs(*args, **kwrags), self.My.subs(*args, **kwrags))

    def simplify(self):
        """Returns a new CMF with simplified Mx and My"""
        return CMF(simplify(self.Mx), simplify(self.My))

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
        
        # Take `trajectory[x]` steps in the x direction
        m = self.Mx.walk({x: 1, y: 0}, trajectory[x], {x: x, y: y})
        # Take `trajectory[y]` steps in the y direction.
        m = m @ self.My.walk({x: 0, y: 1}, trajectory[y], {x: x + trajectory[x], y: y})

        if start is not None:
            m = CMF.substitute_trajectory(m, trajectory, start)
        return simplify(m)

    def as_pcf(self, trajectory, start: dict = {x: 1, y: 1}) -> PCFFromMatrix:
        """
        Returns the PCF equivalent to the CMF in a certain trajectory, up to a mobius transform.
        """
        return self.trajectory_matrix(trajectory, start).as_pcf()

    @staticmethod
    def substitute_trajectory(
        trajectory_matrix: Matrix, trajectory: dict, start: dict
    ) -> Matrix:
        """
        Returns trajectory_matrix reduced to a single variable `n`.

        This transformation is possible only when the starting point is known.
        Each incrementation of the variable `n` represents a full step in `trajectory`.
        """
        return trajectory_matrix.subs({
            x: start[x] + (n-1)*trajectory[x],
            y: start[y] + (n-1)*trajectory[y],
            })

    def walk(
        self, trajectory: dict, iterations: int, start: dict = {x: 1, y: 1}
    ) -> Matrix:
        r"""
        Returns the multiplication matrix of walking in a certain trajectory.

        The walk operation is defined as $\prod_{i=0}^{n-1}M(s_0 + i \cdot t_0, ..., s_k + i \cdot t_k)$,
        where `M=trajectory_matrix(trajectory, start)`, and `n / size(trajectory)` (L1 size - total amount of steps)

        Args:
            trajectory: a dict containing the amount of steps in each direction.
            iterations: the amount of multiplications to perform
            start: a dict representing the starting point of the multiplication.
        Returns:
            the walk multiplication as defined above.
        """
        trajectory_matrix = self.trajectory_matrix(trajectory, start)
        return trajectory_matrix.walk(
            {n: 1},
            iterations // sum(trajectory.values()),
            {n: 1},
        )

    def limit(
        self,
        trajectory: dict,
        iterations: int,
        start: dict = {x: 1, y: 1},
        vector: Vector = Vector.zero(),
    ) -> sp.Float:
        """
        Returns the convergence limit of walking in a certain trajectory.

        This is essentially the same as `self.walk(trajectory, iterations, start).limit(vector)`

        Args:
            trajectory: a dict containing the amount of steps in each direction.
            iterations: the amount of multiplications to perform
            start: a dict representing the starting point of the multiplication.
            vector: The final vector to multiply the matrix by (the zero vector by default)
        Returns:
            the walk multiplication as defined above.
        """
        return self.walk(trajectory, iterations, start).limit(vector)
