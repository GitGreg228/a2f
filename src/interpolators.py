import numpy as np
from scipy.interpolate import interp1d
from scipy.ndimage.filters import gaussian_filter

from utils import stairs


def pre_interpolation(x, y):
    """
    Returns y and x corresponding to positive x
    """
    # x, idxs = np.unique(x, return_index=True)
    # print(x, idxs)
    # _x = np.insert(x[x > 0], 0, 0)
    # _y = np.insert(y[x > 0], 0, 0)
    _x = x[x > 0]
    _y = y[x > 0]
    return _x, _y


def interp_manager(name):
    """
    Desides which function: (lambda, wlog, w2) interpolation will be used
    :param name: function name ('lambda', 'wlog' or 'w2')
    :return:
    """
    if name == 'lambda':
        return interp_lambda
    else:
        return 1


def interp_lambda(x, y, r, sigma):
    """
    Interpolation using for a2f calculation by the formula:
    a2f(w) = w * d(lambda(w))/dw
    :param x: frequencies
    :param y: a2f
    :param r: resolution of desired output a2f
    :param sigma: gaussian filter parameter
    :return: new frequencies, new a2f
    """
    x, y = pre_interpolation(x, y)
    x_new = np.linspace(np.min(x), np.max(x), r)
    dx = x_new[1] - x_new[0]
    lambdas = stairs(y / x)
    f = interp1d(x, lambdas, kind='linear')
    lambdas = f(x_new)
    a2w = np.diff(lambdas)
    a2w = np.insert(a2w, 0, 0)
    a2w = a2w * x_new / dx
    y_new = gaussian_filter(a2w, sigma=sigma)
    return x_new, y_new
