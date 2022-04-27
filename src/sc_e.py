import numpy as np

from pymatgen.core.composition import Composition
from tc import hc, dctc, delta

from constants import *
from utils import stround, parse_formula, floatround, format_e


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
    delta = float()
    gamma = float()
    dctc = float()
    hc = float()
    result = dict()

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

    def get_all(self, system, nef, structure):
        direct = system.direct
        a2f = system.a2f
        self.result['formula'] = parse_formula(structure)
        self.result['Smoothing, THz'] = system.smoothing
        self.result['Resolution, pts'] = system.resolution
        self.result['mu'] = system.mu
        self.result['volume, A^3'] = floatround(structure.volume)
        self.result['density, g/cm^3'] = floatround(structure.density)
        tc = self.tc_e
        self.result['lambda'] = {
            'direct': {
                'lambda': floatround(direct['lambda (lambda)'][-1]),
                'gamma': floatround(direct['lambda (gamma)'][-1])
            },
            'a2f': {
                'lambda': floatround(a2f['lambda (lambda)'][-1]),
                'gamma': floatround(a2f['lambda (gamma)'][-1])
            },
        }
        self.result['wlog, K'] = {
            'direct': {
                'lambda': int(direct['wlog (lambda), K'][-1]),
                'gamma': int(direct['wlog (gamma), K'][-1])
            },
            'a2f': {
                'lambda': int(a2f['wlog (lambda), K'][-1]),
                'gamma': int(a2f['wlog (gamma), K'][-1])
            },
        }
        self.result['w2, K'] = {
            'direct': {
                'lambda': int(direct['w2 (lambda), K'][-1]),
                'gamma': int(direct['w2 (gamma), K'][-1])
            },
            'a2f': {
                'lambda': int(a2f['w2 (lambda), K'][-1]),
                'gamma': int(a2f['w2 (gamma), K'][-1])
            },
        }
        self.result['Tc, K'] = {
            'McMillan': {
                'direct': {
                    'lambda': floatround(direct['Tc_McM (lambda), K'][-1]),
                    'gamma': floatround(direct['Tc_McM (gamma), K'][-1])
                },
                'a2f': {
                    'lambda': floatround(a2f['Tc_McM (lambda), K'][-1]),
                    'gamma': floatround(a2f['Tc_McM (gamma), K'][-1])
                },
            },
            'Allen-Dynes': {
                'direct': {
                    'lambda': floatround(direct['Tc_AD (lambda), K'][-1]),
                    'gamma': floatround(direct['Tc_AD (gamma), K'][-1])
                },
                'a2f': {
                    'lambda': floatround(a2f['Tc_AD (lambda), K'][-1]),
                    'gamma': floatround(a2f['Tc_AD (gamma), K'][-1])
                },
            },
            'Eliashberg': floatround(self.tc_e)
        }
        _lambda = direct['lambda (gamma)'][-1]
        wlog = direct['wlog (gamma), K'][-1]
        self.result['nef'] = {
            'per Unit cell': {
                'states/spin/Ry': format_e(nef),
                'states/spin/eV': format_e(nef / k_Ry_eV),
                'states/spin/J': format_e(nef / k_Ry_J),
                'states/spin/K': format_e(nef / k_Ry_K),
                'states/spin mJ/K^2' : format_e(nef * 1000 * k_B / k_Ry_K),
            },
            'per A^3': {
                'states/spin/Ry': format_e(nef / structure.volume),
                'states/spin/eV': format_e(nef / k_Ry_eV / structure.volume),
                'states/spin/J': format_e(nef / k_Ry_J / structure.volume),
                'states/spin/K': format_e(nef / k_Ry_K / structure.volume),
                'states/spin mJ/K^2': format_e(nef * 1000 * k_B / k_Ry_K / structure.volume),
            },
            'per g': {
                'states/spin/Ry': format_e(nef / structure.volume / structure.density / k_A_cm ** 3),
                'states/spin/eV': format_e(nef / k_Ry_eV / structure.volume / structure.density / k_A_cm ** 3),
                'states/spin/J': format_e(nef / k_Ry_J / structure.volume / structure.density / k_A_cm ** 3),
                'states/spin/K': format_e(nef / k_Ry_K / structure.volume / structure.density / k_A_cm ** 3),
                'states/spin mJ/K^2': format_e(nef * 1000 * k_B / k_Ry_K / structure.volume / structure.density / k_A_cm ** 3),
            }
        }
        gamma = 2 / 3 * (np.pi * k_B) ** 2 * nef * (1 + _lambda) / k_Ry_J  # J/Unit cell/K^2
        self.gamma = gamma
        self.result['Sommerfeld gamma'] = {
            'J/Unit cell/K^2': format_e(gamma),
            'J/A^3/K^2': format_e(gamma / structure.volume),
            'J/cm^3/K^2': format_e(gamma / structure.volume / k_A_cm ** 3),
            'Erg/cm^3/K^2': format_e(gamma * k_J_Erg / structure.volume / k_A_cm ** 3),
            'J/mol/K^2': format_e(gamma * N_A),
            'J/g/K^2': format_e(gamma / structure.volume / structure.density / k_A_cm ** 3),
        }
        self.dctc = dctc(tc, wlog, gamma)
        self.result['DeltaC/Tc'] = {
            'J/Unit cell/K^2': format_e(self.dctc),
            'J/A^3/K^2': format_e(self.dctc / structure.volume),
            'J/cm^3/K^2': format_e(self.dctc / structure.volume / k_A_cm ** 3),
            'J/mol/K^2': format_e(self.dctc * N_A),
            'J/g/K^2': format_e(self.dctc / structure.volume / structure.density / k_A_cm ** 3),
        }
        self.delta = delta(tc, wlog)
        self.result['Delta'] = {
            'J': format_e(self.delta),
            'meV': format_e(self.delta * k_J_meV)
        }
        self.hc = hc(tc, wlog, gamma * k_J_Erg / structure.volume / k_A_cm ** 3)
        self.result['Hc'] = {
            'G': format_e(self.hc),
            'T': format_e(self.hc * k_G_T)
        }
        return self.result
