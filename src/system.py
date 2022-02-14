from qe_ouputs import DynElph
from utils import compare_q_points
from sc_direct import a2f_direct, lambdas_direct, wlogs_direct, w2s_direct
from sc_a2f import lambdas_a2f, wlogs_a2f, w2s_a2f
from interpolators import interp_manager
from tc import tc_mcm, tc_ad


class System(object):
    """
    Calculates a2F, Tc and other parameters
    of superconducting state basing on files
    in given weights and dyn.elphs.
    """
    q_points = list()
    smoothing = float()
    resolution = int()
    interp_name = str()
    direct = dict()
    a2f = dict()
    mu = float()

    def __init__(self, dyn_elphs, weights):
        for dyn_elph in dyn_elphs:
            q_point = dyn_elph.q_point
            for weight in weights:
                if compare_q_points(q_point, weight[0]):
                    q_point_dict = {
                        'DynElph': dyn_elph,
                        'weight': weight[1]
                    }
                    self.q_points.append(q_point_dict)

    def get_direct(self, smoothing):
        self.smoothing = smoothing
        self.direct.update(a2f_direct(self.q_points, smoothing))
        freqs = self.direct['freqs, THz']
        self.direct.update(lambdas_direct(self.q_points, smoothing))
        self.direct.update(wlogs_direct(self.direct))
        self.direct.update(w2s_direct(self.direct))
        return self.direct

    def get_a2f(self, smoothing, resolution, interp_name):
        self.smoothing, self.resolution, self.interp_name = smoothing, resolution, interp_name
        interpolator = interp_manager(interp_name)
        freqs, a2f_lambda = interpolator(self.direct['freqs, THz'], self.direct['a2f (lambda), THz'], resolution)
        _, a2f_gamma = interpolator(self.direct['freqs, THz'], self.direct['a2f (gamma), THz'], resolution)
        self.a2f['freqs, THz'], self.a2f['a2f (lambda), THz'], self.a2f['a2f (gamma), THz'] = \
            freqs, a2f_lambda, a2f_gamma
        self.a2f.update(lambdas_a2f(self.a2f))
        self.a2f.update(wlogs_a2f(self.a2f))
        self.a2f.update(w2s_a2f(self.a2f))
        return self.a2f

    def get_tc(self, mu):
        self.mu = mu
        self.direct.update(tc_mcm(self.direct, mu))
        self.direct.update(tc_ad(self.direct, mu))
        self.a2f.update(tc_mcm(self.a2f, mu))
        self.a2f.update(tc_ad(self.a2f, mu))
