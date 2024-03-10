import numpy as np


class Quat:
    def __init__(self, w, x, y, z):
        self.w = w
        self.x = x
        self.y = y
        self.z = z

    def __str__(self):
        return f'{self.w} {self.x} {self.y} {self.z}'

    def __mul__(self, rhs):
        if isinstance(rhs, Quat):
            # multiply two quaternions
            w = self.w * rhs.w - self.x * rhs.x - self.y * rhs.y - self.z * rhs.z
            x = self.w * rhs.x + self.x * rhs.w + self.y * rhs.z - self.z * rhs.y
            y = self.w * rhs.y - self.x * rhs.z + self.y * rhs.w + self.z * rhs.x
            z = self.w * rhs.z + self.x * rhs.y - self.y * rhs.x + self.z * rhs.w
            return Quat(w, x, y, z)

        elif isinstance(rhs, np.ndarray) and rhs.size == 4:
            # multiply a 4-vector by a quaternion
            qv = self * Quat(rhs[3], rhs[0], rhs[1], rhs[2])
            return np.array([qv.w, qv.x, qv.y, qv.z])

        elif isinstance(rhs, np.ndarray) and rhs.size == 3:
            # multiply a 3-vector by a quaternion
            qv = self * Quat(0, rhs[0], rhs[1], rhs[2]) * self.Inverse()
            return np.array([qv.x, qv.y, qv.z])

        elif isinstance(rhs, float) or isinstance(rhs, int):
            # multiply by an int or float a quaternion
            return Quat(self.w * rhs, self.x * rhs, self.y * rhs, self.z * rhs)

        else:
            raise TypeError(f'Can\'t multiply Quat by {rhs}')

    def Inverse(self):
        return Quat(self.w, -self.x, -self.y, -self.z)
