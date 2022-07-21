import numpy as np

from scipy import integrate
from pymatgen.core.composition import Composition
from tc import hc, dctc, delta, beta

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
    I = 2 * a2f * freqs / (np.square(freqs) + np.square(r * t * k_K_THz * 2 * np.pi))
    return integrate.simps(I, x=freqs)


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
    beta = float()
    dctc = float()
    hc = float()
    result = dict()

    def __init__(self, a2f):
        self.a2f = a2f

    def get_tc_e(self, mu, t_min=0, dim=24):
        """
        :param mu:  Coulomb pseudopotential
        :param t_min: staring temperature in K
        :param dim: Dimension of the matrix
        :return:
        """
        self.dim = dim
        t = t_min
        dt = 10
        while dt > 0.001:
            b = 100
            while b > 0:
                s = np.zeros(dim + 1, dtype=np.float32)
                # print(s)
                for y in range(1, dim + 2):
                    s1 = 0
                    for r in range(1, y + 1):
                        # print(y, r)
                        s1 = s1 + l(self.a2f, r, t)
                    s[y - 1] = s1
                # print(s)
                d = 0
                k = np.zeros((dim + 1, dim + 1), dtype=np.float32)
                for n in range(1, dim + 2):
                    for m in range(1, dim + 2):
                        if n == m:
                            if m == 1:
                                d = 2 * (m - 1) + 1 + l(self.a2f, 0, t)
                            else:
                                d = 2 * (m - 1) + 1 + l(self.a2f, 0, t) + 2 * s[m - 2]
                        else:
                            d = 0
                        k[n - 1, m - 1] = l(self.a2f, m - n, t) + l(self.a2f, m + n - 1, t) - 2 * mu - d
                w, _ = np.linalg.eig(k)
                b = np.max(w)
                self.t.append(t)
                self.b.append(b)
                self.k.append(k)
                if b > 0:
                    t = t + dt
            if t > 0:
                t = t - dt
            if dt > 0.001:
                dt = dt * 0.1
            if dt < 0.001:
                dt = 0.001
        self.tc_e = t - dt
        if self.tc_e <= 0:
            self.tc_e = 0.001
        print(f'Eliashberg Tc = {stround(self.tc_e)}+-{stround(dt)} K')
        return t

    def get_all(self, system, nef, structure):
        direct = system.direct
        a2f = system.a2f
        formula, fu = parse_formula(structure, get_gcd=True)
        self.result['formula'] = formula
        self.result["# of chemical formula units (f.u.)"] = fu
        self.result['Smoothing, THz'] = system.smoothing
        self.result['Resolution, pts'] = system.resolution
        self.result['mu'] = system.mu
        self.result['volume, A^3'] = floatround(structure.volume)
        self.result['volume of chemical formula unit, A^3'] = floatround(structure.volume / fu)
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
        gamma = 2 / 3 * (np.pi * k_B) ** 2 * nef * (1 + _lambda) / k_Ry_J  # J/Unit cell/K^2
        self.gamma = gamma
        self.result['Sommerfeld gamma'] = {
            'mJ/mol/K^2': round(1000 * gamma * N_A / fu, 3),
            'mJ/g/K^2': round(1000 * gamma / structure.volume / structure.density / k_A_cm ** 3, 5),
            'mJ/cm^3/K^2': round(1000 * gamma / structure.volume / k_A_cm ** 3, 4),
            # 'Erg/cm^3/K^2': round(10000000 * gamma / structure.volume / k_A_cm ** 3, 4),
        }
        self.dctc = dctc(tc, wlog, gamma)
        self.result['DeltaC/Tc'] = {
            'mJ/mol/K^2': round(1000 * self.dctc * N_A / fu, 3),
            'mJ/g/K^2': round(1000 * self.dctc / structure.volume / structure.density / k_A_cm ** 3, 3),
            'mJ/cm^3/K^2': round(1000 * self.dctc / structure.volume / k_A_cm ** 3, 4)
        }
        self.delta = delta(tc, wlog)
        self.result['Delta'] = {
            'meV': round(self.delta * k_J_meV, 3),
            '2Delta/kBTc': round(2 * self.delta / (k_B * tc), 3),
            'J': format_e(self.delta)
        }
        self.hc = hc(tc, wlog, 2 * gamma * N_A / fu)
        self.result['Hc2, T'] = round(self.hc, 3),
        self.beta = beta(_lambda, system.mu)
        self.result['beta (McMillan isotope coefficient)'] = round(self.beta, 5)
        self.result['DOS'] = {
            'per Unit cell': {
                'states/spin/Ry': round(nef, 4),
                'states/spin/eV': round(nef / k_Ry_eV, 4),
            },
            f'per chemical formula unit ({formula})': {
                'states/spin/Ry': round(nef / fu, 4),
                'states/spin/eV': round(nef / k_Ry_eV / fu, 4),
            },
            'per A^3': {
                'states/spin/Ry': round(nef / structure.volume, 4),
                'states/spin/eV': round(nef / k_Ry_eV / structure.volume, 4),
            },
            f'per mole of {formula}': {
                'states/spin/Ry': format_e(nef / fu * N_A),
                'states/spin/eV': format_e(nef / k_Ry_eV / fu * N_A),
            },
            'per g': {
                'states/spin/Ry': format_e(nef / structure.volume / structure.density / k_A_cm ** 3),
                'states/spin/eV': format_e(nef / k_Ry_eV / structure.volume / structure.density / k_A_cm ** 3),
            }
        }
        print(f"lambda: {'%.3f' % round(direct['lambda (gamma)'][-1], 3)}")
        print(f"wlog: {round(direct['wlog (gamma), K'][-1])} K")
        print(f"w2: {round(direct['w2 (gamma), K'][-1])} K")
        print(f"McM Tc: {round(direct['Tc_McM (gamma), K'][-1], 1)} K")
        print(f"AD Tc: {round(direct['Tc_AD (gamma), K'][-1], 1)} K")
        print(f"E Tc: {round(tc, 1)} K")
        print(f"DOS: {'%.3f' % round(self.result['DOS']['per A^3']['states/spin/Ry'], 3)} states/spin/Ry/A^3")
        print(f"Sommerfeld gamma: {'%.3f' % round(self.result['Sommerfeld gamma']['mJ/cm^3/K^2'], 3)} mJ/cm^3/K^2")
        print(f"DeltaC/Tc: {'%.3f' % round(self.result['DeltaC/Tc']['mJ/cm^3/K^2'], 3)} mJ/cm^3/K^2")
        print(f"Delta: {'%.1f' % round(self.result['Delta']['meV'], 1)} meV")
        print(f"2Delta/kBTc: {'%.2f' % round(self.result['Delta']['2Delta/kBTc'], 2)}")
        print(f"Hc2: {'%.1f' % round(self.result['Hc2, T'][0], 1)} T")
        return self.result
