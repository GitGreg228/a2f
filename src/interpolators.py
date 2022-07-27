import numpy as np
from scipy.interpolate import interp1d
from scipy import integrate
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
    _y = y[x > 0]
    _x = x[x > 0]
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


def interp_lambda(x, direct_lambda, direct_lambda_smooth, r, sigma, smoothing):
    """
    Interpolation using for a2f calculation by the formula:
    a2f(w) = w * d(lambda(w))/dw
    """
    x_smooth = x
    x, y = pre_interpolation(x, direct_lambda)
    x_new = np.linspace(np.min((0.001, np.min(x))), np.max(x), r)
    # x_new = np.linspace(0, np.max(x), r)
    dx = x_new[1] - x_new[0]
    lambdas = np.insert(y, 0, 0) / 2
    x = np.insert(x, 0, 0)
    f = interp1d(x, lambdas, kind='linear')
    lambdas = f(x_new)
    a2w = np.diff(lambdas)
    a2w = np.insert(a2w, 0, 0)
    a2w = a2w * x_new / dx
    y_new = gaussian_filter(a2w, sigma=sigma)
    smooth = np.exp(-np.square(smoothing / x_new))
    y_new = y_new * smooth
    diff = 2 * integrate.simps(y_new / x_new, x=x_new) / direct_lambda_smooth[-1]
    y_new = y_new / diff
    return x_new, y_new
