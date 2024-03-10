import numpy as np
from quat import Quat


class GrpElemO3:

    def __init__(self, q, inv):
        if (isinstance(q, Quat)):
            self.q = q
        elif (isinstance(q, np.ndarray) and q.size == 4):
            self.q = Quat(q[0], q[1], q[2], q[3])
        else:
            raise TypeError(f'Invalid quaternion parameter {q}')

        if (inv == 1 or inv == -1):
            self.inv = inv
        else:
            raise TypeError(f'Invalid inversion parameter {inv}')

    def __str__(self):
        return f'{self.q} {self.inv}'

    def __mul__(self, rhs):
        if isinstance(rhs, GrpElemO3):
            # multiply two O(3) group elements
            return GrpElemO3(self.q * rhs.q, self.inv * rhs.inv)

        elif isinstance(rhs, np.ndarray) and rhs.size == 3:
            # rotate a vector by this group element
            v = self.q * rhs * self.inv
            return v

        else:
            raise TypeError(f'Can\'t multiply GrpElemO3 by {rhs}')

    def Identity():
        return GrpElemO3(Quat(1.0, 0.0, 0.0, 0.0), 1)
