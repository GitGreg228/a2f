import pandas as pd
import os

import numpy as np


def compare_q_points(q_point1, q_point2):
    """
    Assuming that there are same q-points in Dyn.elphs and in ph.outs,
    we can extract q-point weights from ph.out files.
    However, q-points coordinates have different numbers of
    numbers after floating point in that files.
    This function computes if q-points are same considering
    first 5 numbers after floating point.
    :param q_point1:
    :param q_point2:
    :return:
    """
    result = True
    assert (len(q_point1) == len(q_point2))
    l = len(q_point1)
    for i in range(l):
        result = result * (round(q_point1[i],5) == round(q_point2[i],5))
    return bool(result)


def mkdirs(path, dirname):
    dir_name = os.path.join(path, dirname)
    if not os.path.isdir(dir_name):
        os.mkdir(dir_name)


def save_dict(system, path, float_format='%.12f'):
    smoothing = system.smoothing
    if smoothing == int(smoothing):
        smoothing = int(smoothing)
    smoothing = str(smoothing).replace('.', ',')
    if system.direct:
        name = os.path.join(path, 'results', f'direct_s{smoothing}.csv')
        df = pd.DataFrame.from_dict(system.direct)
        df.to_csv(name, index=False, header=True, sep='\t', float_format=float_format)
    if system.a2f:
        resolution, interp_name = system.resolution, system.interp_name
        name = os.path.join(path, 'results', f'a2f_s{smoothing}_r{resolution}_{interp_name}.csv')
        df = pd.DataFrame.from_dict(system.a2f)
        df.to_csv(name, index=False, header=True, sep='\t', float_format=float_format)


def stairs(arr, dtype=np.float16):
    """
    This function turns an array of values into stairs.
    For example, [1, 2, 3, 2] will be
    turned into [1, 3, 6, 8].
    :param arr:
    :param dtype:
    :return:
    """
    stairs = np.empty(arr.shape, dtype=dtype)
    for i in range(len(arr)):
        stairs[i] = np.sum(arr[:i + 1])
    return stairs

