import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import ticker as mticker
from matplotlib import cm
from matplotlib import colors
from matplotlib.ticker import ScalarFormatter, NullFormatter
import os.path as path
from os import makedirs
from Synthesis.units import *
from Synthesis.post.simulation import simulation


def line_n_planets(pop, a_up_lim = 30, thresholds=None):
    if thresholds is None:
        thresholds = [0.1, 0.4, 1.6]

    thresholds = sorted(thresholds)
    numbers = {key: [] for key in thresholds}
    occurences = {key: [] for key in thresholds}
    systems = []

    for sim in pop.SIMS.values():
        zipped = list(zip(sim.snaps[sim.N_snaps - 1].satellites['M'].values * M_S / M_E,
                          sim.snaps[sim.N_snaps - 1].satellites['a'].values * R_S / au))
        zipped = [(m, a) for (m, a) in zipped if a <= a_up_lim]
        systems.append(zipped)
    number_of_systems = {}
    for keys, lists in numbers.items():
        count = 0
        for i, item in enumerate(systems):
            filtered = [(m, a) for (m, a) in item if
                        m >= keys]
            leng = len(filtered)
            lists += [leng]
            systems[i] = filtered
            if filtered:
                count += 1
            else:
                print(filtered)
            print(count)
        number_of_systems[keys] = count

    max_number = max(numbers[thresholds[0]])
    for keys, value in numbers.items():
        for i in range(max_number + 1):
            occurences[keys] += [value.count(i)]

    plt.rcParams.update({'figure.autolayout': True})
    plt.style.use('seaborn-paper')
    plt.rcParams.update({'font.size': pop.fontsize})
    plt.rcParams.update({"legend.title_fontsize": pop.legend_fontsize})

    fig, ax = plt.subplots(figsize=pop.figsize)
    for keys, value in occurences.items():
        ax.plot(range(max_number + 1), np.array(value)/ len(systems), linestyle='dashed', linewidth=0.7,
                markersize=6,
                marker='o', label=f'Threshold {keys}' + r' $\mathrm{M_{\oplus}}$')
        ax.fill_between(range(max_number + 1),np.array(value)/ len(systems), alpha=0.3)
    ax.set(xlabel='Number of Planets', ylabel=r'Normalized Occurence', xticks=range(max_number + 1))
    ax.legend(loc='best')
    if pop.plot_config == 'presentation':
        ax.set(title=r'Number of Planets Occurences')
    save_name = 'line_planet_number_occurences'
    if a_up_lim < 30:
        save_name += '_lim'
    fig.savefig(path.join(pop.PLOT, save_name + '.png'), transparent=False, dpi=pop.dpi,
                bbox_inches="tight")
    plt.close(fig)


def a_wm_distribution(pop):
    Masses = []
    A = []
    WM = []
    SWM = []
    for sim in pop.SIMS.values():
        Masses.extend(sim.snaps[sim.N_snaps - 1].satellites['M'].values * M_S / M_E)
        A.extend(sim.snaps[sim.N_snaps - 1].satellites['a'].values * R_S / au)
        WM.extend(sim.snaps[sim.N_snaps - 1].satellites['WM'].values * M_S / M_E)
        SWM.extend(sim.snaps[sim.N_snaps - 1].satellites['SWM'].values * M_S / M_E)

    cmap = 'cividis_r'
    cmin = min(Masses)
    cmax = max(Masses)

    norm = colors.LogNorm(cmin, cmax)

    fig, ax = plt.subplots(ncols=1)
    fig.set_size_inches(15.5, 10.5)
    ax.scatter(A, Masses)
    # fig.colorbar(cm.ScalarMappable(cmap=cmap, norm=norm), orientation='vertical', label="Mass in Earth Masses",ax=ax)
    ax.set_xlabel('Distane from Star in au')
    ax.set_ylabel('Mass in Earth Masses')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.legend()
    fig.suptitle('Distance Mass Distribution')
    fig.savefig(path.join(pop.PLOT, 'a_mass_distribution.png'))
    plt.close(fig)

def final_wmf_radial_distribution(pop, cumulative=False, density=False, thresholds=[10e-12]):
    Masses = []
    SWM = []
    WM = []

    for sim in pop.SIMS.values():
        Masses.extend(sim.snaps[sim.N_snaps - 1].satellites['M'].values * M_S / M_E)
        SWM.extend(sim.snaps[sim.N_snaps - 1].satellites['SWM'].values * M_S / M_E)
        WM.extend(sim.snaps[sim.N_snaps - 1].satellites['WM'].values * M_S / M_E)

    WMF = [wm / m for (m, wm) in zip(Masses, WM)]
    SWMF = [swm / m for (m, swm) in zip(Masses, SWM)]

    data = []
    if not density:
        thresholds = [thresh * M_S / M_E for thresh in thresholds]
    for thresh in thresholds:
        data.append([(m, wmf, swmf) for (m, wmf, swmf) in zip(Masses, WMF, SWMF) if m > thresh])

    N_bins = 30
    bins = [np.linspace(min(WMF), max(WMF), N_bins),
            np.linspace(min(SWMF), max(SWMF), N_bins)]

    plt.rcParams.update({'figure.autolayout': True})
    plt.style.use('seaborn-paper')
    plt.rcParams.update({'font.size': pop.fontsize})

    fig, ax = plt.subplots(ncols=2)
    fig.set_size_inches(pop.figsize)

    ax[0].set_title('Water Mass Fraction Distribution')
    ax[1].set_title('Solid Water Mass Fraction Distribution')

    if not density:
        for dat, thresh in zip(data, thresholds):
            ax[0].hist([x[1] for x in dat], bins=bins[0], density=density, cumulative=cumulative,
                       label=str(thresh))
            ax[1].hist([x[2] for x in dat], bins=bins[1], density=density, cumulative=cumulative,
                       label=str(thresh))
    else:
        ax[0].hist(WMF, bins=bins[0], density=density, cumulative=cumulative, stacked=True)
        ax[1].hist(SWMF, bins=bins[1], density=density, cumulative=cumulative, stacked=True)

    ax[0].set_xlabel('Water Mass Fraction')
    ax[1].set_xlabel('Solid Water Mass Fraction')

    if density == False:
        ax[0].set_ylabel('Counts')
        ax[1].set_ylabel('Counts')
    else:
        ax[0].set_ylabel('Probability')
        ax[1].set_ylabel('Probability')

    if pop.plot_config == 'presentation':
        ax[0].set(title=r'Final WMF Distribution')
        ax[1].set(title=r'Final Solid WMF Distribution')
    fig.savefig(path.join(pop.PLOT, 'final_wm_distribution' + str(cumulative) + str(density) + '.png'))
    plt.close(fig)