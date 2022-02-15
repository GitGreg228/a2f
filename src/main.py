import argparse
import os

import matplotlib.pyplot as plt

from qe_ouputs import Folder, PhOuts, DynElph
from system import System
from utils import save_dict, mkdirs, save_structure
from plotters import plot_system


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-p', type=str, help='Path to the file')
    parser.add_argument('-s', type=float, default=3.0, help='Smoothing in THz')
    parser.add_argument('-r', type=int, default=500, help='Interpolation resolution')
    parser.add_argument('-i', type=str, default='lambda', help='Interpolation type')
    parser.add_argument('--mu', type=float, default=0.1, help='mu')
    parser.add_argument('--tol', type=float, default=0.2, help='Structure tolerance')
    args = parser.parse_args()

    folder = Folder(args.p)
    dyn_elphs = [DynElph(path) for path in folder.dyn_elph_paths]
    ph_outs = PhOuts(folder.ph_outs_paths)
    weights = ph_outs.weights
    structure = ph_outs.struc()
    system = System(dyn_elphs, weights)

    mkdirs(args.p, 'results')

    #print(structure.get_primitive_structure(tolerance=0.2))

    save_structure(structure, args.tol, args.p)

    #_ = system.get_direct(args.s)
    #_ = system.get_a2f(args.s, args.r, args.i)
    #system.get_tc(args.mu)
    #save_dict(system, args.p)

    #plot_system(system, args.p)


if __name__ == '__main__':
    main()
