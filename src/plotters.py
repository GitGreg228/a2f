import os

import matplotlib.pyplot as plt

from scipy.ndimage.filters import gaussian_filter


def plot_entity(system, xname, yname, ax, label='', gf=False, fp=True):
    x = system[xname]
    y = system[yname]
    if gf:
        y = gaussian_filter(y, sigma=3)
    _label = label
    if fp:
        if y[-1] <= 10:
            _label = _label + f' ({"%.3f" % round(y[-1], 3)})'
        else:
            _label = _label + f' ({str(round(y[-1]))})'
    ax.plot(x[x > 0], y[x > 0], '--', label=_label)


def plot_entities(system, xname, yname, ax, label='', gf=False, direct=True, a2f=True, fp=True):
    if direct:
        plot_entity(system.direct, xname, yname, ax, label='direct' + label, gf=gf, fp=fp)
    if a2f:
        plot_entity(system.a2f, xname, yname, ax, label='interp' + label, gf=gf, fp=fp)
    ax.legend()


def plot_system(system, path):
    plt.rcParams['lines.markersize'] = 1
    fig, axs = plt.subplots(3, 2, figsize=(12, 16))
    plot_entities(system, 'freqs, THz', 'a2f (lambda), THz', axs[0, 0], label=r' $\alpha^2f$ ($\lambda$)', gf=False, direct=False, fp=False)
    plot_entities(system, 'freqs, THz', 'a2f (gamma), THz', axs[0, 0], label=r' $\alpha^2f$ ($\gamma$)', gf=False, direct=False, fp=False)
    plot_entities(system, 'freqs, THz', 'lambda (lambda)', axs[0, 1], label=r' $\lambda$ ($\lambda$)')
    plot_entities(system, 'freqs, THz', 'lambda (gamma)', axs[0, 1], label=r' $\lambda$ ($\gamma$)')
    plot_entities(system, 'freqs, THz', 'wlog (lambda), K', axs[1, 0], label=r' $\omega_{\mathrm{log}}$ ($\lambda$)')
    plot_entities(system, 'freqs, THz', 'wlog (gamma), K', axs[1, 0], label=r' $\omega_{\mathrm{log}}$ ($\gamma$)')
    plot_entities(system, 'freqs, THz', 'w2 (lambda), K', axs[1, 1], label=r' $\omega^2$ ($\lambda$)')
    plot_entities(system, 'freqs, THz', 'w2 (gamma), K', axs[1, 1], label=r' $\omega^2$ ($\gamma$)')
    plot_entities(system, 'freqs, THz', 'Tc_McM (lambda), K', axs[2, 0], label=r' $T_{\mathrm{C}}$ ($\lambda$)')
    plot_entities(system, 'freqs, THz', 'Tc_McM (gamma), K', axs[2, 0], label=r' $T_{\mathrm{C}}$ ($\gamma$)')
    plot_entities(system, 'freqs, THz', 'Tc_AD (lambda), K', axs[2, 1], label=r' $T_{\mathrm{C}}$ ($\lambda$)')
    plot_entities(system, 'freqs, THz', 'Tc_AD (gamma), K', axs[2, 1], label=r' $T_{\mathrm{C}}$ ($\gamma$)')

    smoothing, resolution, interp_name = system.smoothing, system.resolution, system.interp_name
    plt.suptitle(f'Smoothing: {smoothing} THz, resolution = {resolution}, interpolation: {interp_name}')
    plt.savefig(os.path.join(path, 'results', f'plot_s{smoothing}_r{resolution}_{interp_name}.pdf'))
    plt.show()
