import numpy as np

from constants import *


def tc_mcm_formula(wlog, _lambda, mu):
    """
    McMillan Tc
    :param wlog: logarithmic average frequency
    :param _lambda: electron-phonon coupling parameter
    :param mu: Coulomb pseudopotential
    :return: McMillan Tc
    """
    return wlog / 1.2 * np.exp(-1.04 * (1 + _lambda) / (_lambda * (1 - 0.62 * mu) - mu))


def tc_ad_formula(wlog, w2, _lambda, mu):
    """
    Allen-Dynes Tc
    :param wlog: logarithmic average frequency
    :param w2: mean square frequency
    :param _lambda: electron-phonon coupling parameter
    :param mu: Coulomb pseudopotential
    :return:
    """
    f1 = np.cbrt(1 + np.power(_lambda / (2.46 * (1 + 3.8 * mu)), 1.5))
    f2 = 1 - (np.power(_lambda, 2) * (1 - w2[-1] / wlog[-1])) / (
            np.power(_lambda, 2) + 3.312 * np.power(1 + 6.3 * mu, 2))
    return f1 * f2 * tc_mcm_formula(wlog, _lambda, mu)


def tc_mcm(d, mu):
    """
    McMillan Tc for given System object
    :param d: System.a2f or System.direct
    :param mu: Coulomb pseudopotential
    :return: dictionary with McMillan Tc computed using lambdas and gammas
    """
    wlog_lambda, wlog_gamma = d['wlog (lambda), K'], d['wlog (gamma), K']
    lambda_lambda, lambda_gamma = d['lambda (lambda)'][-1], d['lambda (gamma)'][-1]
    tc_lambda, tc_gamma = map(lambda x, y: tc_mcm_formula(x, y, mu), [wlog_lambda, wlog_gamma],
                              [lambda_lambda, lambda_gamma])
    return {
        'Tc_McM (lambda), K': tc_lambda,
        'Tc_McM (gamma), K': tc_gamma,
    }


def tc_ad(d, mu):
    """
    Allen-Dynes Tc for given System object
    :param d: System.a2f or System.direct
    :param mu: Coulomb pseudopotential
    :return: dictionary with McMillan Tc computed using lambdas and gammas
    """
    wlog_lambda, wlog_gamma = d['wlog (lambda), K'], d['wlog (gamma), K']
    w2_lambda, w2_gamma = d['w2 (lambda), K'], d['w2 (gamma), K']
    lambda_lambda, lambda_gamma = d['lambda (lambda)'][-1], d['lambda (gamma)'][-1]
    tc_lambda, tc_gamma = map(lambda x, y, z: tc_ad_formula(x, y, z, mu), [wlog_lambda, wlog_gamma],
                              [w2_lambda, w2_gamma], [lambda_lambda, lambda_gamma])
    return {
        'Tc_AD (lambda), K': tc_lambda,
        'Tc_AD (gamma), K': tc_gamma,
    }


def hc(tc, wlog, gamma):
    return tc * np.sqrt(gamma / (0.168 * (1 - 12.2 * (tc / wlog) ** 2 * np.log(wlog / (3 * tc)))))


def dctc(tc, wlog, gamma):
    return gamma * 1.43 * (1 + 53 * (tc / wlog) ** 2 * np.log(wlog / (3 * tc)))


def delta(tc, wlog):
    return 0.5 * k_B * k_J_meV * tc * 3.53 * (1 + 12.5 * (tc / wlog) ** 2 * np.log(wlog / (2 * tc)))
