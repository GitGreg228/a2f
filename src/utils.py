import json
import os
import re
from functools import reduce
from math import gcd

import numpy as np
import pandas as pd
from pymatgen.io.cif import CifWriter
from pymatgen.io.vasp.inputs import Poscar
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer


def stround(a, d=3):
    a = float(a)
    return f'%.{str(d)}f' % round(a, d)


def floatround(a, d=3):
    a = float(a)
    return round(a, d)


def format_e(n):
    a = '%E' % n
    return a.split('E')[0].rstrip('0').rstrip('.') + 'E' + a.split('E')[1]


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
        result = result * (round(q_point1[i], 5) == round(q_point2[i], 5))
    return bool(result)


def mkdirs(path, dirname):
    dir_name = os.path.join(path, dirname)
    if not os.path.isdir(dir_name):
        os.mkdir(dir_name)


def save_dict(system, path, float_format='%.8f', b=0):
    smoothing = system.smoothing
    smoothing = str(smoothing).replace('.', ',')
    if system.direct:
        if b > 0:
            name = os.path.join(path, f'direct_s{smoothing}_b{str(b).replace(".", ",")}.csv')
        else:
            name = os.path.join(path, f'direct_s{smoothing}.csv')
        df = pd.DataFrame.from_dict(system.direct)
        df.to_csv(name, float_format=float_format, index=False)
    if system.a2f:
        resolution, sigma = system.resolution, system.sigma
        if b > 0:
            name = os.path.join(path, f'a2f_s{smoothing}_r{resolution}_g{sigma}_b{str(b).replace(".", ",")}.csv')
        else:
            name = os.path.join(path, f'a2f_s{smoothing}_r{resolution}_g{sigma}.csv')
        df = pd.DataFrame.from_dict(system.a2f)
        df.to_csv(name, float_format=float_format, index=False)


def parse_formula(structure, get_gcd=False):
    formula = structure.formula
    # split = formula.split()
    split = re.split(r'(\d+)', formula)
    numeric = list()
    alpha = list()
    for word in split:
        if word.isnumeric():
            numeric.append(int(word))
        else:
            alpha.append(word)
    if numeric:
        GCD = reduce(gcd, numeric)
    else:
        GCD = 1
    result = str()
    for i in range(len(alpha)):
        result = result + alpha[i]
        if not i == len(numeric):
            if not numeric[i] // GCD == 1:
                result = result + str(numeric[i] // GCD)
    if get_gcd:
        return result.replace(' ', ''), GCD
    else:
        return result.replace(' ', '')


def save_structure(struct, tol, path):
    formula, fu = parse_formula(struct, get_gcd=True)
    try:
        analyzer = SpacegroupAnalyzer(struct, symprec=0.2)
        num = analyzer.get_space_group_number()
        sym = analyzer.get_space_group_symbol()
        symm = analyzer.get_symmetrized_structure()
        CifWriter(symm, symprec=tol).write_file(os.path.join(path, 'results', f'{formula}_{str(num)}.cif'))
    except TypeError:
        print('Unable to symmetrize the structure')
        pass
    if int(fu) == 1:
        print(f'\nThe structure under consideration is {sym}-{formula}.')
        Poscar(struct).write_file(os.path.join(path, 'results', f'{formula}.vasp'))
    else:
        print(f'\nThe structure under consideration is {fu}x of {sym}-{formula}.')
        Poscar(struct).write_file(os.path.join(path, 'results', f'{fu}x{formula}.vasp'))


def stairs(arr, dtype=np.float16):
    """
    This function turns an array of values into stairs.
    For example, [1, 2, 3, 2] will be
    turned into [1, 3, 6, 8].
    :param arr:
    :param dtype:
    :return:
    """
    stairs = np.zeros(arr.shape, dtype=dtype)
    for i in range(len(arr)):
        stairs[i] = np.sum(arr[:i])
    return stairs


def print_direct(system):
    _lambda = stround(system.direct["lambda (gamma)"][-1])
    wlog = round(system.direct["wlog (gamma), K"][-1])
    w2 = round(system.direct["w2 (gamma), K"][-1])
    print(f'With the smoothing of {system.smoothing} THz, direct lambda = {_lambda}, '
          f'direct wlog = {wlog} K, direct w2 = {w2} K.')


def print_a2f(system):
    resolution = system.resolution
    sigma = system.sigma
    _lambda = stround(system.a2f["lambda (gamma)"][-1])
    wlog = round(system.a2f["wlog (gamma), K"][-1])
    w2 = round(system.a2f["w2 (gamma), K"][-1])
    print(f'After calculating alpha2f on resolution r = {resolution} '
          f'and applying gaussian filter with sigma = {sigma}, '
          f'interp lambda = {_lambda}, interp wlog = {wlog} K, interp w2 = {w2} K.')


def print_tc(system):
    mu = system.mu
    tc_mcm_direct = stround(system.direct['Tc_McM (gamma), K'][-1])
    tc_ad_direct = stround(system.direct['Tc_AD (gamma), K'][-1])
    tc_mcm_a2f = stround(system.a2f['Tc_McM (gamma), K'][-1])
    tc_ad_a2f = stround(system.a2f['Tc_AD (gamma), K'][-1])
    print(f'Mu = {mu}. Calculated McMillan Tc = {tc_mcm_direct} ({tc_mcm_a2f}) K, '
          f'calculated Allen-Dynes Tc = {tc_ad_direct} ({tc_ad_a2f}) K')


def save_result(result, path):
    with open(os.path.join(path, 'result.json'), 'w', encoding='utf-8') as f:
        json.dump(result, f, ensure_ascii=False, sort_keys=False, indent=4)


def update_summary(summary, result, broadening, formula, s, mu, path):
    keys = ['Broadening, Ry',
            'lambda (direct)',
            'lambda (interp)',
            'wlog, K (direct)',
            'wlog, K (interp)',
            'w2, K (direct)',
            'w2, K (interp)',
            'Tc McM, K (direct)',
            'Tc McM, K (interp)',
            'Tc AD, K (direct)',
            'Tc AD, K (interp)',
            'Tc E, K',
            'Sommerfeld gamma, mJ/cm^3/K^2',
            'DOS, states/spin/Ry/A^3',
            'DeltaC/Tc, mJ/cm^3/K^2',
            'Delta, meV',
            'Hc2, T',
            'beta']
    if not summary:
        summary = {key: list() for key in keys}
    summary['Broadening, Ry'].append(broadening)
    summary['lambda (direct)'].append(result['lambda']['direct']['gamma'])
    summary['lambda (interp)'].append(result['lambda']['a2f']['gamma'])
    summary['wlog, K (direct)'].append(result['wlog, K']['direct']['gamma']['K'])
    summary['wlog, K (interp)'].append(result['wlog, K']['a2f']['gamma']['K'])
    summary['w2, K (direct)'].append(result['w2, K']['direct']['gamma']['K'])
    summary['w2, K (interp)'].append(result['w2, K']['a2f']['gamma']['K'])
    summary['Tc McM, K (direct)'].append(result['Tc, K']['McMillan']['direct']['gamma'])
    summary['Tc McM, K (interp)'].append(result['Tc, K']['McMillan']['a2f']['gamma'])
    summary['Tc AD, K (direct)'].append(result['Tc, K']['Allen-Dynes']['direct']['gamma'])
    summary['Tc AD, K (interp)'].append(result['Tc, K']['Allen-Dynes']['a2f']['gamma'])
    summary['Tc E, K'].append(result['Tc, K']['Eliashberg'])
    summary['Sommerfeld gamma, mJ/cm^3/K^2'].append(result['Sommerfeld gamma']['mJ/cm^3/K^2'])
    summary['DOS, states/spin/Ry/A^3'].append(result['DOS']['per A^3']['states/spin/Ry'])
    summary['DeltaC/Tc, mJ/cm^3/K^2'].append(result['DeltaC/Tc']['mJ/cm^3/K^2'])
    summary['Delta, meV'].append(result['Delta']['meV'])
    summary['Hc2, T'].append(result['Hc2, T'][0])
    summary['beta'].append(result['beta (McMillan isotope coefficient)'])
    df = pd.DataFrame.from_dict(summary)
    s = str(s).replace('.', ',')
    mu = str(mu).replace('.', ',')
    df.to_csv(os.path.join(path, f'summary_{formula}_s{s}_mu{mu}.csv'), index=False)
    return summary
