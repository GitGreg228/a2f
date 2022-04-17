import numpy as np

from pymatgen.core.composition import Composition
from tc import hc, dctc, delta

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
    Superconducting class is used for calculation of Eliashberg Tc.
    It stores a2f (taken as input) + Kmn matrix, Kmn matrix maximum eigenvalue on each T.
    """
    a2f = dict()
    dim = int()
    tc_e = float()
    t = list()
    k = list()
    b = list()
    gamma = float()
    hc = float()

    def __init__(self, a2f):
        self.a2f = a2f

    def get_tc_e(self, mu, t_min=1, dim=24):
        """
        :param mu:  Coulomb pseudopotential
        :param t_min: staring temperature in K
        :param dim: Dimension of the matrix
        :return:
        """
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
            self.k.append(k)
            dt = a * b
            if dt > 10:
                # a = 0.1 * a
                dt = 10
            if dt < 0.001:
                dt = 0.001
                # a = 10 * a
            t = t + dt
        self.tc_e = t - dt
        print(f'Eliashberg Tc = {stround(t - dt)}+-{stround(dt)} K')
        return t

    def get_all(self, direct, nef):
        tc = self.tc_e
        _lambda = direct['lambda (gamma)'][-1]
        print(f'Lambda {"%.3f" % round(_lambda, 3)}')
        wlog = direct['wlog (gamma), K'][-1]
        print(f'wlog {"%.3f" % round(wlog, 3)} K')
        weight = 1 / N_A
        nef = nef / (k_Ry_J * weight)
        self.gamma = 2 / 3 * (np.pi * k_B)**2 * nef * (1 + _lambda) * 1000
        print(f'Sommerfeld gamma\t {"%.3f" % round(self.gamma, 3)} mJ/(mol K2)')
        self.hc = hc(tc, wlog, self.gamma)
        #print(self.hc)
        self.delta = delta(tc, wlog)
        print(f'Delta\t {"%.3f" % round(self.delta, 3)} meV')
        self.dctc = dctc(tc, wlog, self.gamma)
        print(f'DeltaC/Tc\t {"%.3f" % round(self.dctc, 3)} mJ/(mol K2)')
