import numpy as np


def tc_mcm_formula(wlog, _lambda, mu):
    # assert isinstance(_lambda, np.float32)
    out = wlog / 1.2 * np.exp(-1.04 * (1 + _lambda) / (_lambda * (1 - 0.62 * mu) - mu))
    return out


def tc_ad_formula(wlog, w2, _lambda, mu):
    f1 = np.cbrt(1 + np.power(_lambda / (2.46 * (1 + 3.8 * mu)), 1.5))
    f2 = 1 - (np.power(_lambda, 2) * (1 - w2[-1] / wlog[-1])) / (np.power(_lambda, 2) + 3.312 * np.power(1 + 6.3 * mu, 2))
    out = f1 * f2 * tc_mcm_formula(wlog, _lambda, mu)
    return out


def tc_mcm(d, mu):
    wlog_lambda, wlog_gamma = d['wlog (lambda), K'], d['wlog (gamma), K']
    lambda_lambda, lambda_gamma = d['lambda (lambda)'][-1], d['lambda (gamma)'][-1]
    tc_lambda, tc_gamma = map(lambda x, y: tc_mcm_formula(x, y, mu), [wlog_lambda, wlog_gamma],
                              [lambda_lambda, lambda_gamma])
    out_dict = {
        'Tc_McM (lambda), K': tc_lambda,
        'Tc_McM (gamma), K': tc_gamma,
    }
    return out_dict


def tc_ad(d, mu):
    wlog_lambda, wlog_gamma = d['wlog (lambda), K'], d['wlog (gamma), K']
    w2_lambda, w2_gamma = d['w2 (lambda), K'], d['w2 (gamma), K']
    lambda_lambda, lambda_gamma = d['lambda (lambda)'][-1], d['lambda (gamma)'][-1]
    tc_lambda, tc_gamma = map(lambda x, y, z: tc_ad_formula(x, y, z, mu), [wlog_lambda, wlog_gamma],
                              [w2_lambda, w2_gamma], [lambda_lambda, lambda_gamma])
    out_dict = {
        'Tc_AD (lambda), K': tc_lambda,
        'Tc_AD (gamma), K': tc_gamma,
    }
    return out_dict
