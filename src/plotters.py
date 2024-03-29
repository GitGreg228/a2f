import os

import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axisartist.parasite_axes import HostAxes, ParasiteAxes

from utils import stround


def plot_entity(system, xname, yname, ax, label='', fp=True, lw=2.5, ms=2.5, zorder=1):
    x = system[xname]
    y = system[yname]
    _label = label
    if fp:
        if y[-1] <= 10:
            _label = _label + f' ({stround(y[-1])})'
        else:
            _label = _label + f' ({str(round(y[-1]))})'
    style = 'o-'
    if lw == 0:
        style = 'o'
    if ms == 0:
        style = '-'
    if 'gamma' in yname:
        style = style.replace('o', 'x')
        style = style.replace('-', ':')
        zorder = zorder + 1
    ax.plot(x[x > 0], y[x > 0], style, linewidth=lw, markersize=ms, label=_label, zorder=zorder)


def fill_entity(system, xname, yname, ax, label='', linestyle='-', fc='lightgray', ec='gray'):
    x = system[xname]
    y = system[yname]
    ax.fill_between(x, y, linestyle=linestyle, label=label, facecolor=fc, edgecolor=ec, alpha=0.9)
    ax.legend()


def plot_spectra(system, xname, yname, ax, label='', marker='o'):
    x = system[xname]
    y = system[yname]
    # ax.vlines(x[y > 0], 0, y[y > 0], color=color, linestyle=linestyle, linewidth=1, label=label)
    ax.scatter(x[y > 0], y[y > 0], label=label, s=2.5, marker=marker)
    ax.legend()


def plot_entities(system, xname, yname, ax, label='', direct=True, a2f=True, fp=True):
    if direct:
        plot_entity(system.direct, xname, yname, ax, label='direct' + label, fp=fp, ms=4, lw=0, zorder=3)
    if a2f:
        plot_entity(system.a2f, xname, yname, ax, label='interp' + label, fp=fp, ms=0, lw=1.5, zorder=1)
    ax.legend()


def plot_system(system, formula, path, b=0):
    plt.close('all')
    plt.rcParams['lines.markersize'] = 2.5
    plt.rcParams['lines.linewidth'] = 2.5
    plt.rcParams['font.size'] = 12
    fig, axs = plt.subplots(3, 2, figsize=(12, 16))

    for ax in axs.flatten():
        ax.set_xlabel('$\omega$, THz')

    axs[0, 0].set_ylim(0, 1.05 * max(system.a2f['a2f (gamma), THz']))
    axs[0, 0].set_title(r'Eliashberg function $\alpha^2f(\omega)$')
    axs[0, 0].set_ylabel(r'$\alpha^2f$')
    axs[0, 1].set_title(r'EPC coefficient $\lambda$')
    axs[0, 1].set_ylabel(r'$\lambda$')
    axs[1, 0].set_title(r'Logarithmic average frequency $\omega_{\mathrm{log}}$')
    axs[1, 0].set_ylabel(r'$\omega_{\mathrm{log}}$, K')
    axs[1, 1].set_title(r'Mean square frequency $\omega^2$')
    axs[1, 1].set_ylabel(r'$\omega^2$, K')
    axs[2, 0].set_title(r'McMillan $T_{\mathrm{C}}$')
    axs[2, 0].set_ylabel(r'$T_{\mathrm{C}}$, K')
    axs[2, 1].set_title(r'Allen-Dynes $T_{\mathrm{C}}$')
    axs[2, 1].set_ylabel(r'$T_{\mathrm{C}}$, K')

    fill_entity(system.a2f, 'freqs, THz', 'a2f (lambda), THz', axs[0, 0], label=r' interp ($\lambda$)', fc='aliceblue',
                ec='deepskyblue')
    fill_entity(system.a2f, 'freqs, THz', 'a2f (gamma), THz', axs[0, 0], label=r' interp ($\gamma$)', fc='papayawhip',
                ec='lightsalmon')
    plot_spectra(system.direct, 'freqs, THz', 'a2f (lambda), THz', axs[0, 0], label='direct ($\lambda$)', marker='x')
    plot_spectra(system.direct, 'freqs, THz', 'a2f (gamma), THz', axs[0, 0], label='direct ($\gamma$)', marker='o')
    plot_entities(system, 'freqs, THz', 'lambda (lambda)', axs[0, 1], label=r'($\lambda$)')
    plot_entities(system, 'freqs, THz', 'lambda (gamma)', axs[0, 1], label=r'($\gamma$)')
    plot_entities(system, 'freqs, THz', 'wlog (lambda), K', axs[1, 0], label=r'($\lambda$)')
    plot_entities(system, 'freqs, THz', 'wlog (gamma), K', axs[1, 0], label=r'($\gamma$)')
    plot_entities(system, 'freqs, THz', 'w2 (lambda), K', axs[1, 1], label=r'($\lambda$)')
    plot_entities(system, 'freqs, THz', 'w2 (gamma), K', axs[1, 1], label=r'($\gamma$)')
    plot_entities(system, 'freqs, THz', 'Tc_McM (lambda), K', axs[2, 0], label=r'($\lambda$)')
    plot_entities(system, 'freqs, THz', 'Tc_McM (gamma), K', axs[2, 0], label=r'($\gamma$)')
    plot_entities(system, 'freqs, THz', 'Tc_AD (lambda), K', axs[2, 1], label=r'($\lambda$)')
    plot_entities(system, 'freqs, THz', 'Tc_AD (gamma), K', axs[2, 1], label=r'($\gamma$)')

    smoothing, resolution, sigma = system.smoothing, system.resolution, system.sigma

    info = f'Smoothing: {smoothing} THz, resolution = {resolution}, $\sigma$ = {sigma}, $\mu$=' + f'{system.mu}'
    if b > 0:
        info = info + f', broadening = {str(b)} Ry'
    fig.suptitle(f'{formula}\n{info}')

    fig.tight_layout(rect=[0, 0, 1, 0.95], h_pad=2.5)

    if b > 0:
        name = f'plot_s{smoothing}_r{resolution}_g{sigma}_b{str(b).replace(".", ",")}.pdf'
    else:
        name = f'plot_s{smoothing}_r{resolution}_g{sigma}.pdf'
    plt.savefig(os.path.join(path, name))
    # plt.show()


def plot_article_view(system, formula, path):
    plt.close('all')
    plt.rcParams['lines.markersize'] = 3
    plt.rcParams['lines.linewidth'] = 3
    plt.rcParams['font.size'] = 20

    fig = plt.figure(figsize=(12, 6))

    a2f = system.a2f
    x = a2f['freqs, THz']

    ax1 = HostAxes(fig, [0.15, 0.1, 0.65, 0.8])
    ax1.axis["right"].set_visible(False)
    ax1.set_xlim(0, max(x))
    # ax1.set_xticks(np.linspace(0, 7 * np.round(np.max(x) / 7), 8))
    ax2 = ParasiteAxes(ax1, sharex=ax1)
    ax3 = ParasiteAxes(ax1, sharex=ax1)
    ax4 = ParasiteAxes(ax1, sharex=ax1)
    ax1.parasites.append(ax2)
    ax1.parasites.append(ax3)
    ax1.parasites.append(ax4)
    fig.add_axes(ax1)

    y = a2f['a2f (gamma), THz']
    ax1.fill_between(x, y, fc='lightgrey', ec='black')
    ax1.set_ylim(0, 1.05 * max(y))
    ax1.set_ylabel(r'$\alpha^2f$')
    ax1.set_xlabel('$\omega$, THz')

    y = a2f['lambda (gamma)']
    new_axisline = ax2.get_grid_helper().new_fixed_axis
    ax2.axis["right2"] = new_axisline(loc="right", axes=ax2)
    ax2.set_ylim(0, 1.05 * max(y))
    ax2.plot(x, y, color='steelblue')
    ax2.set_ylabel('$\lambda$', color='steelblue', labelpad=1)

    y = a2f['wlog (gamma), K']
    new_axisline = ax3.get_grid_helper().new_fixed_axis
    ax3.axis["right2"] = new_axisline(loc="right", axes=ax3, offset=(90, 0))
    ax3.set_ylim(0, 1.1 * max(y))
    ax3.set_yticks(np.linspace(0, 5 * np.round(np.max(y) / 5), 6))
    ax3.plot(x, y, linestyle='-', color='darkorange')
    ax3.set_ylabel(r'$\omega_{\mathrm{log}}$, K', color='darkorange', labelpad=1)

    y = a2f['Tc_AD (gamma), K']
    new_axisline = ax4.get_grid_helper().new_fixed_axis
    ax4.axis["right2"] = new_axisline(loc="right", axes=ax4, offset=(180, 0))
    ax4.set_ylim(0, 1.15 * max(y))
    ax4.set_yticks(np.linspace(0, 5 * np.round(np.max(y) / 5), 6))
    ax4.plot(x, y, linestyle='-', color='darkgreen')
    ax4.set_ylabel(r'Allen-Dynes $T_{\mathrm{C}}$, K', color='darkgreen', labelpad=1)

    smoothing, resolution, sigma = system.smoothing, system.resolution, system.sigma
    fig.savefig(os.path.join(path, f'plot_article_{formula}.pdf'), bbox_inches='tight')

    # plt.show()


def plot_summary(summary, formula, smoothing, sigma, resolution, mu, path):
    plt.close('all')
    plt.rcParams['lines.markersize'] = 3
    plt.rcParams['lines.linewidth'] = 2
    plt.rcParams['font.size'] = 15

    fig = plt.figure(figsize=(12, 6))

    ax1 = HostAxes(fig, [0.15, 0.1, 0.65, 0.8])
    ax1.axis["right"].set_visible(False)
    # ax1.set_xlim(0, np.max(summary['Broadening, Ry']))
    # ax1.set_xticks(np.linspace(0, 7 * np.round(np.max(x) / 7), 8))
    ax2 = ParasiteAxes(ax1, sharex=ax1)
    ax1.parasites.append(ax2)
    fig.add_axes(ax1)
    new_axisline = ax2.get_grid_helper().new_fixed_axis
    ax2.axis["right2"] = new_axisline(loc="right", axes=ax2)

    keys = sorted([key for key in summary.keys() if key.startswith('Tc')])
    for key in keys:
        x = summary['Broadening, Ry']
        y = summary[key]
        fancy_key = key.replace('Tc', r'$T_{\mathrm{C}}$').replace('McM', '(McMillan)').replace('AD',
                                                                                                '(Allen-Dynes)').replace(
            'E', ('Eliashberg'))
        ax1.plot(x, y, label=fancy_key)
    # keys = sorted([key for key in summary.keys() if key.startswith('lambda')])
    keys = ['lambda (direct)']
    for key in keys:
        x = summary['Broadening, Ry']
        y = summary[key]
        fancy_key = key.replace('lambda', r'$\lambda$')
        ax2.plot(x, y, 'o', markersize=4, label=fancy_key)
    ax1.legend()
    # ax2.legend()
    ax1.set_xticks(summary['Broadening, Ry'])
    ax1.set_xticklabels(list([str(int(label*1000)) for label in summary['Broadening, Ry']]), fontsize=10)
    ax1.set_xlabel('Broadening, mRy')
    ax1.set_ylabel(r'$T_{\mathrm{C}}$')
    ax2.set_ylabel('$\lambda$', color='k')
    info = f'Smoothing: {smoothing} THz, resolution = {resolution}, $\sigma$ = {sigma}, $\mu$=' + f'{mu}'
    ax1.set_title(f'Convergence of {formula}\n{info}')
    s = str(smoothing).replace('.', ',')
    # r = str(resolution).replace('.', ',')
    mu = str(mu).replace('.', ',')
    fig.savefig(os.path.join(path, f'summary_{formula}_s{s}_mu{mu}.pdf'), bbox_inches='tight')
