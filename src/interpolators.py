import numpy as np

from scipy.interpolate import interp1d
from scipy.ndimage.filters import gaussian_filter

from utils import stairs


def pre_interpolation(x, y):
    """
    This function removes duplicates in x array
    and sums values from y array that correspond
    to removed duplicates.
    """
    # x, idxs = np.unique(x, return_index=True)
    # print(x, idxs)
    # _x = np.insert(x[x > 0], 0, 0)
    # _y = np.insert(y[x > 0], 0, 0)
    _x = x[x > 0]
    _y = y[x > 0]
    return _x, _y


def interp_manager(name):
    if name == 'lambda':
        return interp_lambda
    else:
        return 1


def interp_lambda(x, y, r, sigma):
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
    while np.max(y_new) > np.max(y) and sigma < 3:
        sigma = sigma + 0.1
        y_new = gaussian_filter(a2w, sigma=sigma)
    return x_new, y_new

