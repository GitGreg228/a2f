import pandas as pd
import os

import numpy as np


class Folder(object):
    """
    Parses a folder with QE ph.x calculations,
    then gives paths to different output files
    and also the prefix (formula)
    """
    __path = ''
    dyn_elph_paths = list()

    def __init__(self, path):
        self.__path = path
        self.dyn_elph_paths = self.dyn_elphs()
        self.ph_outs_paths = self.ph_outs()

    def dyn_elphs(self):
        dyn_elphs_paths = list(filter(lambda x: 'dyn' in x and 'elph' in x,
                                      os.listdir(self.__path)))
        full_dyn_elphs_paths = [os.path.join(self.__path, dyn_elphs_path) for dyn_elphs_path in dyn_elphs_paths]
        if full_dyn_elphs_paths:
            print(f'Found dyn_elphs in {", ".join(full_dyn_elphs_paths)}')
            return full_dyn_elphs_paths
        else:
            print(f'Error: Unable to detect *dyn*.elph* files in {self.__path}')
            raise FileNotFoundError

    def ph_outs(self):
        ph_outs_paths = list(filter(lambda x: '.ph' in x and 'out' in x,
                                    os.listdir(self.__path)))
        full_ph_outs_paths = [os.path.join(self.__path, ph_outs_path) for ph_outs_path in ph_outs_paths]
        if full_ph_outs_paths:
            print(f'Found ph.outs in {", ".join(full_ph_outs_paths)}')
            return full_ph_outs_paths
        else:
            print('Warning: Unable to detech ph.out files in ', self.__path)
            return [self.__path]

    def a2f_dats(self, *p):
        a2f_dats_paths = list(filter(lambda x: 'a2F.dat' in x, os.listdir(self.__path)))
        if a2f_dats_paths:
            if len(a2f_dats_paths) == 1:
                return a2f_dats_paths[0]
            else:
                print('Warning: Detected several a2F.dat files in ', self.__path)
                print('Specify which to parse')
                return a2f_dats_paths[p]
        else:
            print('Warning: Unable to detech a2F.dat files in ', self.__path)
            return [self.__path]

    def formula(self):
        formula = os.path.basename(self.dyn_elphs()[0]).split(".")[0]
        return formula


class DynElph(object):
    """
    Parses a single *dyn*.elph* file,
    returning all it contains:
    q-point, lambdas, gammas and squared frequencies
    """
    lines = list()
    q_point = tuple()
    sqr_freqs = list()
    dos_line = str()
    E_F = float()
    DOS = float()
    lambda_gamma_lines = list()
    lambdas = list()
    gammas = list()

    def __init__(self, path):
        with open(path) as read_obj:
            lines = read_obj.readlines()
            read_obj.close()
        self.lines = lines
        self.dos_line = next(line for line in lines if 'DOS' in line).split()
        self.lambda_gamma_lines = list(filter(lambda x: 'lambda' in x and 'gamma' in x, lines))
        self.q_point()
        self.sqr_freqs()
        self.e_fermi()
        self.dos()
        self.lambdas()
        self.gammas()

    def q_point(self):
        self.q_point = tuple(float(x) for x in self.lines[0].split()[0:3])
        return self.q_point

    def sqr_freqs(self):
        lines = self.lines
        sqr_freqs = list()
        idx = lines.index(next(line for line in lines if 'method' in line))
        for line in lines[1:idx]:
            for freq in line.strip().split():
                sqr_freqs.append(float(freq))
        self.sqr_freqs = np.array(sqr_freqs)
        return self.sqr_freqs

    def e_fermi(self):
        self.E_F = float(self.dos_line[7])
        return self.E_F

    def dos(self):
        self.DOS = float(self.dos_line[2])
        return self.DOS

    def lambdas(self):
        lambdas = list()
        for line in self.lambda_gamma_lines:
            lambdas.append(float(line.split()[2]))
        self.lambdas = np.array(lambdas)
        return self.lambdas

    def gammas(self):
        gammas = list()
        for line in self.lambda_gamma_lines:
            gammas.append(float(line.split()[4]))
        self.gammas = np.array(gammas)
        return self.gammas

    def as_dict(self):
        dyn_elph_dict = {
            'q_point': self.q_point(),
            'sqr_freqs': self.sqr_freqs(),
            'e_fermi': self.e_fermi(),
            'dos': self.dos(),
            'lambdas': self.lambdas(),
            'gammas': self.gammas()
        }
        return dyn_elph_dict


class PhOuts(object):
    """
    Parses given list of ph.out files,
    returns weights (and soon will be parsing frequencies).
    Works both with multiple ph.out files
    and single ph.out file.
    """
    __lines = list()
    weights = list()
    __paths = list()

    def __init__(self, paths):
        self.__paths = paths
        if len(paths) > 1 or (len(paths) == 1 and not os.path.isdir(paths[0])):
            for path in paths:
                with open(path, 'r') as f:
                    _lines = f.readlines()
                self.__lines.append(_lines)
        self.weights()

    def weights(self):
        """
        output format: (q-points, weights)
        :return:
        """
        weights = list()
        q_points = list()
        if self.__lines:
            for _lines in self.__lines:
                weight_lines = list(filter(lambda x: 'number of k points=' in x and 'method' in x, _lines))
                occupancies = [i for i, x in enumerate(_lines) if 'Number of q in the star' in x]
                for occupancy in occupancies:
                    q_point = tuple([float(x) for x in _lines[occupancy+2].strip().split()[-3:]])
                    weight = int(_lines[occupancy].split()[-1])
                    print(f'q = ({", ".join(["%.3f" % round(_q, 3) for _q in q_point])}) '
                          f'with number of q in the star {weight}')
                    q_points.append(q_point)
                    weights.append(weight)
        else:
            print('Assuming no-symmetry calculation, consequently setting uniform weights')
            weights = [1] * len(Folder(self.__paths[0]).dyn_elphs())
        weights_sum = sum(weights)
        weights = [weight / weights_sum for weight in weights]
        self.weights = list(zip(q_points, weights))
        return self.weights


