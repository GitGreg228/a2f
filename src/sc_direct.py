import numpy as np

from constants import *
from utils import stairs


def a2f_direct(qpoints, smoothing, dtype=np.float16):
    """
    Computes Eliashberg functions a2F directly from lambdas and gammas without any integration.
    """
    freqs, a2f_lambda, a2f_gamma = list(), list(), list()
    for q_point in qpoints:
        dyn_elph = q_point['DynElph']
        lambdas = dyn_elph.lambdas
        gammas = dyn_elph.gammas
        nef = dyn_elph.DOS
        weight = q_point['weight']
        q_freqs = np.sqrt(np.abs(dyn_elph.sqr_freqs))*np.sign(dyn_elph.sqr_freqs) * k_Ry_THz
        smooth = np.exp(-np.square(smoothing / q_freqs))
        q_a2f_lambda = smooth * lambdas * q_freqs / 2 * weight
        q_a2f_gamma = smooth * gammas / q_freqs / nef / k_Ry_GHz / np.pi / 2 * weight * k_Ry_THz ** 2
        a2f_lambda.append(q_a2f_lambda)
        a2f_gamma.append(q_a2f_gamma)
        freqs.append(q_freqs)
    freqs, a2f_lambda, a2f_gamma = map(lambda x: np.array(x).flatten(), (freqs, a2f_lambda, a2f_gamma))
    idxs = freqs.argsort()
    freqs = freqs[idxs]
    a2f_lambda = a2f_lambda[idxs]
    a2f_gamma = a2f_gamma[idxs]
    a2f_lambda[a2f_lambda <= 0] = 0
    a2f_gamma[a2f_gamma <= 0] = 0
    out_dict = {
        'freqs, THz': freqs,
        'a2f (lambda), THz': a2f_lambda,
        'a2f (gamma), THz': a2f_gamma
    }
    return out_dict


def lambdas_direct(qpoints, smoothing, dtype=np.float16):
    """
    Computes lambdas directly from lambdas and gammas without any integration.
    WARNING: Avoid using dtype=np.float32 or better, sometimes it leads to computational errors.
    """
    freqs, lambdas_lambda, lambdas_gamma = list(), list(), list()
    for q_point in qpoints:
        dyn_elph = q_point['DynElph']
        lambdas = dyn_elph.lambdas
        gammas = dyn_elph.gammas
        nef = dyn_elph.DOS
        weight = q_point['weight']
        q_freqs = np.sqrt(np.abs(dyn_elph.sqr_freqs)) * np.sign(dyn_elph.sqr_freqs)
        smooth = np.exp(-np.square(smoothing / (q_freqs * k_Ry_THz)))
        q_lambdas_lambda = smooth * lambdas * weight
        q_lambdas_gamma = smooth * gammas * weight / np.square(q_freqs) / nef / k_Ry_GHz / np.pi
        lambdas_lambda.append(q_lambdas_lambda)
        lambdas_gamma.append(q_lambdas_gamma)
        freqs.append(q_freqs)
    freqs, lambdas_lambda, lambdas_gamma = map(lambda x: np.array(x).flatten(), (freqs, lambdas_lambda, lambdas_gamma))
    idxs = freqs.argsort()
    # freqs = freqs[idxs]
    lambdas_lambda = lambdas_lambda[idxs]
    lambdas_gamma = lambdas_gamma[idxs]
    lambda_lambda_stairs = stairs(lambdas_lambda)
    lambda_gamma_stairs = stairs(lambdas_gamma)
    out_dict = {
        'lambdas (lambda)': lambdas_lambda,
        'lambdas (gamma)': lambdas_gamma,
        'lambda (lambda)': lambda_lambda_stairs,
        'lambda (gamma)': lambda_gamma_stairs
    }
    return out_dict


def wlogs_direct(direct):
    """
    Computes wlogs directly from lambdas and gammas without any integration.
    First you need to calculate lambdas.
    """
    freqs, lambdas_lambda, lambdas_gamma = direct['freqs, THz'], direct['lambdas (lambda)'], direct['lambdas (gamma)']
    lambda_lambda = np.sum(lambdas_lambda)
    lambda_gamma = np.sum(lambdas_gamma)
    freqs_K = freqs[freqs > 0] / k_Ry_THz * k_Ry_K
    zeros = np.zeros(freqs[freqs <= 0].shape)
    wlogs_lambda = lambdas_lambda[freqs > 0] * np.log(freqs_K) / lambda_lambda
    wlogs_gamma = lambdas_gamma[freqs > 0] * np.log(freqs_K) / lambda_gamma
    wlogs_lambda[wlogs_lambda <= 0] = 0
    wlogs_gamma[wlogs_gamma <= 0] = 0
    wlogs_lambda_stairs = np.concatenate((np.zeros(wlogs_lambda[wlogs_lambda <= 0].shape),
                                          np.exp(stairs(wlogs_lambda[wlogs_lambda > 0]))))
    wlogs_gamma_stairs = np.concatenate((np.zeros(wlogs_gamma[wlogs_gamma <= 0].shape),
                                         np.exp(stairs(wlogs_gamma[wlogs_gamma > 0]))))
    wlogs_lambda, wlogs_gamma, wlogs_lambda_stairs, wlogs_gamma_stairs = \
        map(lambda x: np.concatenate((zeros, x)), (wlogs_lambda, wlogs_gamma, wlogs_lambda_stairs, wlogs_gamma_stairs))
    out_dict = {
        'wlogs (lambda), K': wlogs_lambda,
        'wlogs (gamma), K': wlogs_gamma,
        'wlog (lambda), K': wlogs_lambda_stairs,
        'wlog (gamma), K': wlogs_gamma_stairs,
    }
    return out_dict


def w2s_direct(direct):
    """
    Computes w2s directly from lambdas and gammas without any integration.
    First you need to calculate lambdas.
    """
    freqs, lambdas_lambda, lambdas_gamma = direct['freqs, THz'], direct['lambdas (lambda)'], direct['lambdas (gamma)']
    lambda_lambda = np.sum(lambdas_lambda)
    lambda_gamma = np.sum(lambdas_gamma)
    freqs_K = freqs / k_Ry_THz * k_Ry_K
    w2s_lambda = lambdas_lambda * np.square(freqs_K) / lambda_lambda
    w2s_gamma = lambdas_gamma * np.square(freqs_K) / lambda_gamma
    w2s_lambda_stairs = np.sqrt(stairs(w2s_lambda, dtype=np.longdouble))
    w2s_gamma_stairs = np.sqrt(stairs(w2s_gamma, dtype=np.longdouble))
    out_dict = {
        'w2s (lambda), K': w2s_lambda,
        'w2s (gamma), K': w2s_gamma,
        'w2 (lambda), K': w2s_lambda_stairs,
        'w2 (gamma), K': w2s_gamma_stairs,
    }
    return out_dict
