import argparse

from plotters import plot_system
from qe_ouputs import Folder, PhOuts, DynElph
from sc_e import Superconducting
from system import System
from utils import save_dict, mkdirs, save_structure, parse_formula, print_direct, print_a2f, print_tc


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-p', type=str, default='.', help='Path to the directory')
    parser.add_argument('-s', type=float, default=3, help='Smoothing in THz')
    parser.add_argument('-r', type=int, default=0, help='Interpolation resolution. 0 => r = len(freqs[freqs > 0])')
    parser.add_argument('-g', type=float, default=1, help='gaussian filter sigma')
    parser.add_argument('--mu', type=float, default=0.1, help='mu')
    parser.add_argument('--tol', type=float, default=0.2, help='Structure tolerance')
    args = parser.parse_args()

    mkdirs(args.p, 'results')

    folder = Folder(args.p)
    dyn_elphs = [DynElph(path) for path in folder.dyn_elph_paths]
    ph_outs = PhOuts(folder.ph_outs_paths)
    structure = ph_outs.struc()
    save_structure(structure, args.tol, args.p)

    weights = ph_outs.weights()
    system = System(dyn_elphs, weights)

    if args.s == int(args.s):
        args.s = int(args.s)
    if args.g == int(args.g):
        args.g = int(args.g)

    _ = system.get_direct(args.s)
    print_direct(system)
    _ = system.get_a2f(args.r, args.g)
    print_a2f(system)
    system.get_tc(args.mu)
    print_tc(system)

    sc = Superconducting(system.a2f)
    sc.get_tc_e(args.mu)
    save_dict(system, args.p)
    plot_system(system, parse_formula(structure), args.p)


if __name__ == '__main__':
    main()
