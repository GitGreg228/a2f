import numpy as np
from scipy import integrate

from constants import *


def lambdas_a2f(a2f, dtype=np.float16):
    """
    Computes lambda from given a2f.
    """
    freqs, a2f_lambda, a2f_gamma = a2f['freqs, THz'], a2f['a2f (lambda), THz'], a2f['a2f (gamma), THz']
    dw = freqs[1] - freqs[0]
    lambda_lambda = np.zeros(freqs.shape, dtype=dtype)
    lambda_gamma = np.zeros(freqs.shape, dtype=dtype)
    for i in range(len(freqs)):
        # lambda_lambda[i] = 2 * integrate.simps(a2f_lambda[:i+1] / freqs[:i+1])
        # lambda_gamma[i] = 2 * integrate.simps(a2f_gamma[:i + 1] / freqs[:i + 1])
        lambda_lambda[i] = 2 * np.sum(a2f_lambda[:i + 1] / freqs[:i + 1])
        lambda_gamma[i] = 2 * np.sum(a2f_gamma[:i + 1] / freqs[:i + 1])
    out_dict = {
        'lambda (lambda)': lambda_lambda,
        'lambda (gamma)': lambda_gamma,
    }
    return out_dict


def wlogs_a2f(a2f, dtype=np.float16):
    """
    Computes wlog from given a2f.
    """
    freqs, a2f_lambda, a2f_gamma = a2f['freqs, THz'], a2f['a2f (lambda), THz'], a2f['a2f (gamma), THz']
    lambda_lambda, lambda_gamma = a2f['lambda (lambda)'][-1], a2f['lambda (gamma)'][-1]
    wlog_lambda = wlog_gamma = np.zeros(freqs.shape, dtype=dtype)
    for i in range(len(freqs)):
        # wlog_lambda[i] = np.exp(2 / lambda_lambda * integrate.simps(a2f_lambda[:i+1] * np.log(freqs[:i+1] / k_K_THz) / freqs[:i+1]))
        # wlog_gamma[i] = np.exp(2 / lambda_gamma * integrate.simps(a2f_gamma[:i + 1] * np.log(freqs[:i+1] / k_K_THz) / freqs[:i+1]))
        wlog_lambda[i] = np.exp(2 / lambda_lambda * np.sum(a2f_lambda[:i + 1] * np.log(freqs[:i + 1] / k_K_THz) / freqs[:i + 1]))
        wlog_gamma[i] = np.exp(2 / lambda_gamma * np.sum(a2f_gamma[:i + 1] * np.log(freqs[:i + 1] / k_K_THz) / freqs[:i + 1]))
    out_dict = {
        'wlog (lambda), K': wlog_lambda,
        'wlog (gamma), K': wlog_gamma,
    }
    return out_dict


def w2s_a2f(a2f, dtype=np.float16):
    """
    Computes wlog from given a2f.
    """
    freqs, a2f_lambda, a2f_gamma = a2f['freqs, THz'], a2f['a2f (lambda), THz'], a2f['a2f (gamma), THz']
    lambda_lambda, lambda_gamma = a2f['lambda (lambda)'][-1], a2f['lambda (gamma)'][-1]
    w2_lambda = w2_gamma = np.zeros(freqs.shape, dtype=dtype)
    for i in range(len(freqs)):
        # wlog_lambda[i] = np.exp(2 / lambda_lambda * integrate.simps(a2f_lambda[:i+1] * np.log(freqs[:i+1] / k_K_THz) / freqs[:i+1]))
        # wlog_gamma[i] = np.exp(2 / lambda_gamma * integrate.simps(a2f_gamma[:i + 1] * np.log(freqs[:i+1] / k_K_THz) / freqs[:i+1]))
        w2_lambda[i] = np.sqrt(2 / lambda_lambda * np.sum(a2f_lambda[:i + 1] * freqs[:i + 1] / k_K_THz ** 2))
        w2_gamma[i] = np.sqrt(2 / lambda_gamma * np.sum(a2f_gamma[:i + 1] * freqs[:i + 1] / k_K_THz ** 2))
    out_dict = {
        'w2 (lambda), K': w2_lambda,
        'w2 (gamma), K': w2_gamma,
    }
    return out_dict
