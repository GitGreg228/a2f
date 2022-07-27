import argparse

import matplotlib.pyplot as plt
import numpy as np

from plotters import plot_system, plot_article_view
from qe_outputs import Folder, PhOuts, DynElph, Dyn
from sc_e import Superconducting
from system import System
from utils import save_dict, mkdirs, save_structure, parse_formula, print_direct, print_a2f, print_tc, save_result


def main():
    """
    The manager of all processes that parses input and create output using classes and functions from other files.
    :return:
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('-p', type=str, default='.', help='Path to the directory')
    parser.add_argument('-s', type=float, default=3, help='Smoothing in THz')
    parser.add_argument('-r', type=int, default=1000, help='Interpolation resolution. 0 => r = len(freqs[freqs > 0])')
    parser.add_argument('-g', type=float, default=5, help='gaussian filter sigma')
    parser.add_argument('--mu', type=float, default=0.1, help='Coulomb pseudopotential')
    parser.add_argument('--tol', type=float, default=0.2, help='Structure tolerance')
    parser.add_argument('--int', type=str, choices=['simps', 'sum'], default='sum', help='Integration method')
    args = parser.parse_args()
    mkdirs(args.p, 'results')

    folder = Folder(args.p)
    dyn_elphs = [DynElph(path) for path in folder.dyn_elph_paths]
    try:
        ph_outs = PhOuts(folder.ph_outs_paths)
        structure = ph_outs.struc()
        save_structure(structure, args.tol, args.p)
    except FileNotFoundError or KeyError:
        print('Found no ph.out files. The structure will not be read and written.')
        structure = 'Unknown'
        pass
    try:
        dyns = [Dyn(path) for path in folder.dyn_paths]
        q_points = [dyn.q_point for dyn in dyns]
        weights = [dyn.weight for dyn in dyns]
        weights_sum = sum(weights)
        weights = [weight / weights_sum for weight in weights]
        weights = list(zip(q_points, weights))
    except TypeError or FileNotFoundError:
        print('Found no Dyn files, trying to extract weight info from ph.out file(s)')
        weights = ph_outs.weights()
    system = System(dyn_elphs, weights)

    if args.s == int(args.s):
        args.s = int(args.s)
    if args.g == int(args.g):
        args.g = int(args.g)

    direct = system.get_direct(args.s)
    print_direct(system)
    a2f = system.get_a2f(args.r, args.g, args.int)
    print_a2f(system)
    system.get_tc(args.mu)
    print_tc(system)

    save_dict(system, args.p)
    if isinstance(structure, str):
        plot_system(system, 'Unknown', args.p)
        plot_article_view(system, 'Unknown', args.p)
    else:
        plot_system(system, parse_formula(structure), args.p)
        plot_article_view(system, parse_formula(structure), args.p)

    sc = Superconducting(a2f)
    sc.get_tc_e(args.mu)
    nef = dyn_elphs[0].dos()
    result = sc.get_all(system, nef, structure)
    save_result(result, args.p)


if __name__ == '__main__':
    main()
