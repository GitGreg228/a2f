import numpy as np
from scipy import integrate

from constants import *


def lambdas_a2f(system, dtype=np.float32):
    """
    Computes lambda from given a2f.
    """
    a2f = system.a2f
    freqs, a2f_lambda, a2f_gamma = a2f['freqs, THz'], a2f['a2f (lambda), THz'], a2f['a2f (gamma), THz']
    lambda_lambda = np.zeros(freqs.shape, dtype=dtype)
    lambda_gamma = np.zeros(freqs.shape, dtype=dtype)
    for i in range(1, len(freqs)):
        if system.integration == 'simps':
            lambda_lambda[i] = 2 * integrate.simps(a2f_lambda[:i + 1] / freqs[:i + 1], x=freqs[:i + 1])
            lambda_gamma[i] = 2 * integrate.simps(a2f_gamma[:i + 1] / freqs[:i + 1], x=freqs[:i + 1])
        elif system.integration == 'sum':
            dw = freqs[2] - freqs[1]
            lambda_lambda[i] = 2 * np.sum(a2f_lambda[:i + 1] / freqs[:i + 1] * dw)
            lambda_gamma[i] = 2 * np.sum(a2f_gamma[:i + 1] / freqs[:i + 1] * dw)
    return {
        'lambda (lambda)': lambda_lambda,
        'lambda (gamma)': lambda_gamma,
    }


def wlogs_a2f(system, dtype=np.float16):
    """
    Computes wlog from given a2f.
    """
    a2f = system.a2f
    freqs, a2f_lambda, a2f_gamma = a2f['freqs, THz'], a2f['a2f (lambda), THz'], a2f['a2f (gamma), THz']

    lambda_lambda, lambda_gamma = a2f['lambda (lambda)'][-1], a2f['lambda (gamma)'][-1]
    wlog_lambda = wlog_gamma = np.zeros(freqs.shape, dtype=dtype)
    for i in range(1, len(freqs)):
        if system.integration == 'simps':
            wlog_lambda[i] = np.exp(2 / lambda_lambda * integrate.simps(a2f_lambda[:i + 1] * np.log(freqs[:i + 1] / k_K_THz) / freqs[:i + 1], x=freqs[:i + 1]))
            wlog_gamma[i] = np.exp(2 / lambda_gamma * integrate.simps(a2f_gamma[:i + 1] * np.log(freqs[:i + 1] / k_K_THz) / freqs[:i + 1], x=freqs[:i + 1]))
        elif system.integration == 'sum':
            dw = freqs[1] - freqs[0]
            wlog_lambda[i] = np.exp(
                2 / lambda_lambda * np.sum(a2f_lambda[:i + 1] * np.log(freqs[:i + 1] / k_K_THz) / freqs[:i + 1]) * dw)
            wlog_gamma[i] = np.exp(
                2 / lambda_gamma * np.sum(a2f_gamma[:i + 1] * np.log(freqs[:i + 1] / k_K_THz) / freqs[:i + 1]) * dw)
    return {
        'wlog (lambda), K': wlog_lambda,
        'wlog (gamma), K': wlog_gamma,
    }


def w2s_a2f(system, dtype=np.float16):
    """
    Computes wlog from given a2f.
    """
    a2f = system.a2f
    freqs, a2f_lambda, a2f_gamma = a2f['freqs, THz'], a2f['a2f (lambda), THz'], a2f['a2f (gamma), THz']
    dw = freqs[1] - freqs[0]
    lambda_lambda, lambda_gamma = a2f['lambda (lambda)'][-1], a2f['lambda (gamma)'][-1]
    w2_lambda = w2_gamma = np.zeros(freqs.shape, dtype=dtype)
    for i in range(len(freqs)):
        if system.integration == 'simps':
            w2_lambda[i] = np.sqrt(
                2 / lambda_lambda * integrate.simps(a2f_lambda[:i + 1] * freqs[:i + 1] / k_K_THz ** 2, x=freqs[:i + 1]))
            w2_gamma[i] = np.sqrt(
                2 / lambda_gamma * integrate.simps(a2f_gamma[:i + 1] * freqs[:i + 1] / k_K_THz ** 2, x=freqs[:i + 1]))
        elif system.integration == 'sum':
            dw = freqs[1] - freqs[0]
            w2_lambda[i] = np.sqrt(2 / lambda_lambda * np.sum(a2f_lambda[:i + 1] * freqs[:i + 1] / k_K_THz ** 2) * dw)
            w2_gamma[i] = np.sqrt(2 / lambda_gamma * np.sum(a2f_gamma[:i + 1] * freqs[:i + 1] / k_K_THz ** 2) * dw)
    return {
        'w2 (lambda), K': w2_lambda,
        'w2 (gamma), K': w2_gamma,
    }
