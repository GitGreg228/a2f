from interpolators import interp_lambda
from sc_a2f import lambdas_a2f, wlogs_a2f, w2s_a2f
from sc_direct import a2f_direct, lambdas_direct, wlogs_direct, w2s_direct
from tc import tc_mcm, tc_ad
from utils import compare_q_points


class System(object):
    """
    This class stores all numbers from *.dyn*.elph* files and weights. It also calculates a2f and other parameters:
    lambda, wlog, w2 and Tc (McMillan + Allen-Dynes)
    """
    q_points = list()
    smoothing = float()
    resolution = int()
    integration = str()
    sigma = float()
    direct = dict()
    a2f = dict()
    mu = float()

    def __init__(self, dyn_elphs, weights):
        """
        :param dyn_elphs: list of DynElph objects (the class is defined in qe_outputs)
        :param weights: the list of weights taken from ph.out
        """
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
        self.direct.update(lambdas_direct(self.q_points, smoothing))
        self.direct.update(wlogs_direct(self.direct))
        self.direct.update(w2s_direct(self.direct))
        return self.direct

    def get_a2f(self, resolution, sigma, integration):
        self.resolution, self.sigma, self.integration = resolution, sigma, integration
        if not self.resolution:
            self.resolution = len(self.direct['freqs, THz'][self.direct['freqs, THz'] > 0])
        interpolator = interp_lambda
        freqs, a2f_lambda = interpolator(self.direct['freqs, THz'], self.direct['lambda (lambda, not smooth)'],
                                         self.direct['lambda (lambda)'], self.resolution, sigma, self.smoothing)
        _, a2f_gamma = interpolator(self.direct['freqs, THz'], self.direct['lambda (gamma, not smooth)'],
                                         self.direct['lambda (gamma)'], self.resolution, sigma, self.smoothing)
        self.a2f['freqs, THz'], self.a2f['a2f (lambda), THz'], self.a2f['a2f (gamma), THz'] = \
            freqs, a2f_lambda, a2f_gamma
        self.a2f.update(lambdas_a2f(self))
        self.a2f.update(wlogs_a2f(self))
        self.a2f.update(w2s_a2f(self))
        return self.a2f

    def get_tc(self, mu):
        self.mu = mu
        self.direct.update(tc_mcm(self.direct, mu))
        self.direct.update(tc_ad(self.direct, mu))
        self.a2f.update(tc_mcm(self.a2f, mu))
        self.a2f.update(tc_ad(self.a2f, mu))
