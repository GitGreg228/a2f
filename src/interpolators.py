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
    #x, idxs = np.unique(x, return_index=True)
    #print(x, idxs)
    return x[x > 0], y[x > 0]


def interp_manager(name):
    if name == 'lambda':
        return interp_lambda
    else:
        return 1


def interp_lambda(x, y, r):
    x, y = pre_interpolation(x, y)
    x_new = np.linspace(np.min(x), np.max(x), r)
    dx = x_new[1] - x_new[0]
    lambdas = stairs(y / x)
    f = interp1d(x, lambdas, kind='linear')
    lambdas = f(x_new)
    a2w = np.diff(lambdas)/np.diff(x_new)
    a2w = np.insert(a2w, 0, 0)
    a2w = gaussian_filter(a2w, sigma=4)
    y_new = a2w * x_new * dx
    return x_new, y_new

