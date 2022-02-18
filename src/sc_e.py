import numpy as np

from constants import *
from utils import stround


def l(a2f, r, t):
    """
    a2f kernel
    :param a2f: a2f dictionary of the System object
    :param r: matrix index difference
    :param t: temperature in K
    :return: kernel integral
    """
    freqs = a2f['freqs, THz']
    a2f = a2f['a2f (gamma), THz']
    dw = freqs[1] - freqs[0]
    return np.sum(2 * a2f * freqs * dw / (np.square(freqs) + np.square(r * t * k_K_THz * 2 * np.pi)))


class Superconducting(object):
    """
    Superconducting class describes
    """
    a2f = dict()
    dim = int()
    t = list()
    k = list()
    b = list()

    def __init__(self, a2f):
        self.a2f = a2f

    def get_tc_e(self, mu, t_min=0.1, dt=0.001, dim=24):
        self.dim = dim
        t = t_min
        b = 100
        a = 10
        while b > 0:
            s = np.arange(1, dim + 1)
            for y in range(1, dim + 1):
                s1 = 0
                for r in range(1, y + 1):
                    s1 = s1 + l(self.a2f, r, t)
                s[y - 1] = s1
            d = 0
            k = np.zeros((dim, dim), dtype=np.float32)
            for n in range(1, dim + 1):
                for m in range(1, dim + 1):
                    if n == m:
                        if m == 1:
                            d = 2 * (m - 1) + 1 + l(self.a2f, 0, t)
                        else:
                            d = 2 * (m - 1) + 1 + l(self.a2f, 0, t) + 2 * s[m - 1]
                    else:
                        d = 0
                    k[n - 1, m - 1] = l(self.a2f, m - n, t) + l(self.a2f, m + n - 1, t) - 2 * mu - d
            w, _ = np.linalg.eig(k)
            b = np.max(w)
            self.t.append(t)
            self.b.append(b)
            dt = a * b
            if dt > 10:
                dt = 10
            if dt < 0.001:
                dt = 0.001
            t = t + dt
        print(f'Eliashberg Tc = {stround(t - dt)}+-{stround(dt)} K')
        return t
