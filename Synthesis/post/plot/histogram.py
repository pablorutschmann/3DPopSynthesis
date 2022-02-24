import matplotlib.pyplot as plt
import os.path as path

import numpy as np

from Synthesis.units import *


# import Synthesis.post.histogram as hist

def histogram_mass(pop, m_low_lim=0, a_up_lim=30):
    Masses = []
    Orb_Dist = []
    for sim in pop.SIMS.values():
        Masses += list(sim.snaps[sim.N_snaps - 1].satellites['M'].values * M_S / M_E)
        Orb_Dist += list(sim.snaps[sim.N_snaps - 1].satellites['a'].values * R_S / au)

    data = zip(Masses, Orb_Dist)
    data = [item for item in data if item[0] >= m_low_lim and item[1] <= a_up_lim]

    Masses, Orb_Dist = zip(*data)

    plt.rcParams.update({'figure.autolayout': True})
    plt.style.use('seaborn-paper')
    plt.rcParams.update({'font.size': pop.fontsize})

    N_bins = 15
    bins = 10 ** np.linspace(np.log10(min(Masses)), np.log10(max(Masses)), N_bins)
    fig, ax = plt.subplots(figsize=pop.figsize)
    # ax.hist(Masses, bins=bins)
    values, base, _ = plt.hist(Masses, bins=bins, rwidth=0.95)
    ax.axvline(1, color='red', linewidth=1)
    ax.axvline(M_M / M_E, color='red', linewidth=1)
    ax.axvline(M_V / M_E, color='red', linewidth=1)
    ax.axvline(M_ME / M_E, color='red', linewidth=1)
    ax_bis = ax.twinx()
    # values, base = np.histogram(Masses, bins = 10 ** np.linspace(np.log10(min(Masses)), np.log10(max(Masses)), N_bins))

    values = np.append(0, values)
    ax_bis.plot(base, np.cumsum(values) / np.cumsum(values)[-1], color='black', linestyle='dashed', markersize=0.1)
    ax.set(xlabel=r'Mass [$M_E$]', ylabel=r'Counts')
    ax_bis.set(ylabel='Cumulative Distribution')
    ax.set_xscale('log')
    ax_bis.set_xscale('log')
    if pop.plot_config == 'presentation':
        ax.set(title=r'Histrogram of Planet Masses')
    save_name = 'histogram_mass'
    if a_up_lim < 30 and m_low_lim > 0:
        save_name += '_lim'
    fig.savefig(path.join(pop.PLOT, save_name + '.png'), transparent=False, dpi=pop.dpi, bbox_inches="tight")
    plt.close(fig)


def histogram_weighted_mass(pop, m_low_lim=0, a_up_lim=30):
    Masses = []
    Orb_Dist = []

    for sim in pop.SIMS.values():
        Masses += list(sim.snaps[sim.N_snaps - 1].satellites['M'].values * M_S / M_E)
        Orb_Dist += list(sim.snaps[sim.N_snaps - 1].satellites['a'].values * R_S / au)

    data = zip(Masses, Orb_Dist)

    data = [item for item in data if item[0] >= m_low_lim and item[1] <= a_up_lim]

    Masses, Orb_Dist = zip(*data)

    plt.rcParams.update({'figure.autolayout': True})
    plt.style.use('seaborn-paper')
    plt.rcParams.update({'font.size': pop.fontsize})

    N_bins = 15
    bins = 10 ** np.linspace(np.log10(min(Masses)), np.log10(max(Masses)), N_bins)
    fig, ax = plt.subplots(figsize=pop.figsize)
    # ax.hist(Masses, bins=bins)
    values, base, _ = plt.hist(Masses, bins=bins, rwidth=0.95, weights=Masses / np.sum(Masses))
    ax.axvline(1, color='red', linewidth=1)
    ax.axvline(M_M / M_E, color='red', linewidth=1)
    ax.axvline(M_V / M_E, color='red', linewidth=1)
    ax.axvline(M_ME / M_E, color='red', linewidth=1)
    ax_bis = ax.twinx()
    # values, base = np.histogram(Masses, bins = 10 ** np.linspace(np.log10(min(Masses)), np.log10(max(Masses)), N_bins))

    values = np.append(0, values)
    ax_bis.plot(base, np.cumsum(values) / np.cumsum(values)[-1], color='black', linestyle='dashed', markersize=0.1)
    ax.set(xlabel=r'Mass [$M_E$]', ylabel=r'Weighed Counts')
    ax_bis.set(ylabel='Cumulative Distribution')
    ax.set_xscale('log')
    ax_bis.set_xscale('log')
    if pop.plot_config == 'presentation':
        ax.set(title=r'Weighted Histrogram of Planet Masses')
    save_name = 'histogram_weighted_mass'
    if a_up_lim < 30 and m_low_lim > 0:
        save_name += '_lim'
    fig.savefig(path.join(pop.PLOT, save_name + '.png'), transparent=False, dpi=pop.dpi, bbox_inches="tight")
    plt.close(fig)


def histogram_weighted_mass_nonlog(pop, m_low_lim=0, a_up_lim=30):
    Masses = []
    Orb_Dist = []

    for sim in pop.SIMS.values():
        Masses += list(sim.snaps[sim.N_snaps - 1].satellites['M'].values * M_S / M_E)
        Orb_Dist += list(sim.snaps[sim.N_snaps - 1].satellites['a'].values * R_S / au)

    data = zip(Masses, Orb_Dist)

    data = [item for item in data if item[0] >= m_low_lim and item[1] <= a_up_lim]

    Masses, Orb_Dist = zip(*data)

    plt.rcParams.update({'figure.autolayout': True})
    plt.style.use('seaborn-paper')
    plt.rcParams.update({'font.size': pop.fontsize})

    fig, ax = plt.subplots(figsize=pop.figsize)
    # ax.hist(Masses, bins=bins)
    values, base, _ = plt.hist(Orb_Dist, bins=15, rwidth=0.95, weights=Masses / np.sum(Masses))
    ax.axvline(1, color='red', linewidth=1)
    ax.axvline(M_M / M_E, color='red', linewidth=1)
    ax.axvline(M_V / M_E, color='red', linewidth=1)
    ax.axvline(M_ME / M_E, color='red', linewidth=1)
    ax_bis = ax.twinx()
    # values, base = np.histogram(Masses, bins = 10 ** np.linspace(np.log10(min(Masses)), np.log10(max(Masses)), N_bins))

    values = np.append(0, values)
    ax_bis.plot(base, np.cumsum(values) / np.cumsum(values)[-1], color='black', linestyle='dashed', markersize=0.1)
    ax.set(xlabel=r'M [$M_E$]', ylabel=r'Weighed Counts')
    ax_bis.set(ylabel='Cumulative Distribution')
    # ax.set_xscale('log')
    # ax_bis.set_xscale('log')
    if pop.plot_config == 'presentation':
        ax.set(title=r'Weighted Histrogram of Planet Masses')
    save_name = 'histogram_weighted_orb_dist_nonlog'
    if a_up_lim < 30 and m_low_lim > 0:
        save_name += '_lim'
    fig.savefig(path.join(pop.PLOT, save_name + '.png'), transparent=False, dpi=pop.dpi, bbox_inches="tight")
    plt.close(fig)


def histogram_a(pop, m_low_lim=0, a_up_lim=30):
    Masses = []
    Orb_Dist = []

    for sim in pop.SIMS.values():
        Masses += list(sim.snaps[sim.N_snaps - 1].satellites['M'].values * M_S / M_E)
        Orb_Dist += list(sim.snaps[sim.N_snaps - 1].satellites['a'].values * R_S / au)

    data = zip(Masses, Orb_Dist)

    data = [item for item in data if item[0] >= m_low_lim and item[1] <= a_up_lim]

    Masses, Orb_Dist = zip(*data)

    plt.rcParams.update({'figure.autolayout': True})
    plt.style.use('seaborn-paper')
    plt.rcParams.update({'font.size': pop.fontsize})

    N_bins = 15
    bins = 10 ** np.linspace(np.log10(min(Orb_Dist)), np.log10(max(Orb_Dist)), N_bins)
    fig, ax = plt.subplots(figsize=pop.figsize)
    # ax.hist(Masses, bins=bins)
    values, base, _ = plt.hist(Orb_Dist, bins=bins, rwidth=0.95)
    ax.axvline(1, color='red', linewidth=1)
    ax.axvline(0.387, color='red', linewidth=1)
    ax.axvline(0.732, color='red', linewidth=1)
    ax.axvline(1.52, color='red', linewidth=1)
    ax_bis = ax.twinx()
    # values, base = np.histogram(Masses, bins = 10 ** np.linspace(np.log10(min(Masses)), np.log10(max(Masses)), N_bins))

    values = np.append(0, values)
    ax_bis.plot(base, np.cumsum(values) / np.cumsum(values)[-1], color='black', linestyle='dashed', markersize=0.1)
    ax.set(xlabel=r'Orbital distance [$\mathrm{au}$]', ylabel=r'Counts')
    ax_bis.set(ylabel='Cumulative Distribution')
    ax.set_xscale('log')
    ax_bis.set_xscale('log')
    if pop.plot_config == 'presentation':
        ax.set(title=r'Histrogram of Planet Orbital Distances')
    save_name = 'histogram_a'
    if a_up_lim < 30 and m_low_lim > 0:
        save_name += '_lim'
    fig.savefig(path.join(pop.PLOT, save_name + '.png'), transparent=False, dpi=pop.dpi, bbox_inches="tight")
    plt.close(fig)


def histogram_weighted_a(self, m_low_lim=0, a_up_lim=30):
    Masses = []
    Orb_Dist = []

    for sim in self.SIMS.values():
        Masses += list(sim.snaps[sim.N_snaps - 1].satellites['M'].values * M_S / M_E)
        Orb_Dist += list(sim.snaps[sim.N_snaps - 1].satellites['a'].values * R_S / au)

    data = zip(Masses, Orb_Dist)

    data = [item for item in data if item[0] >= m_low_lim and item[1] <= a_up_lim]

    Masses, Orb_Dist = zip(*data)

    plt.rcParams.update({'figure.autolayout': True})
    plt.style.use('seaborn-paper')
    plt.rcParams.update({'font.size': self.fontsize})

    N_bins = 15
    bins = 10 ** np.linspace(np.log10(min(Orb_Dist)), np.log10(max(Orb_Dist)), N_bins)
    fig, ax = plt.subplots(figsize=self.figsize)
    # ax.hist(Masses, bins=bins)
    values, base, _ = plt.hist(Orb_Dist, bins=bins, rwidth=0.95, density=True, weights=Masses / np.sum(Masses))
    ax.axvline(1, color='red', linewidth=1)
    ax.axvline(0.387, color='red', linewidth=1)
    ax.axvline(0.732, color='red', linewidth=1)
    ax.axvline(1.52, color='red', linewidth=1)
    ax_bis = ax.twinx()
    # values, base = np.histogram(Masses, bins = 10 ** np.linspace(np.log10(min(Masses)), np.log10(max(Masses)), N_bins))

    values = np.append(0, values)
    ax_bis.plot(base, np.cumsum(values) / np.cumsum(values)[-1], color='black', linestyle='dashed', markersize=0.1)
    ax.set(xlabel=r'Orbital distance [$\mathrm{au}$]', ylabel=r'Weighted Counts')
    ax_bis.set(ylabel='Cumulative Distribution')
    ax.set_xscale('log')
    ax_bis.set_xscale('log')
    if self.plot_config == 'presentation':
        ax.set(title=r'Weighted Histogram of Planet Orbital Distances')
    save_name = 'histogram_weighted_a'
    if a_up_lim < 30 and m_low_lim > 0:
        save_name += '_lim'
    fig.savefig(path.join(self.PLOT, save_name + '.png'), transparent=False, dpi=self.dpi, bbox_inches="tight")
    plt.close(fig)


def histogram_a_twmf(pop, m_low_lim=0, a_up_lim=30):
    Masses = []
    Orb_Dist = []
    WM = []
    SWM = []
    System = []

    for key, sim in pop.SIMS.items():
        Masses += list(sim.snaps[sim.N_snaps - 1].satellites['M'].values * M_S / M_E)
        Orb_Dist += list(sim.snaps[sim.N_snaps - 1].satellites['a'].values * R_S / au)
        WM += list(sim.snaps[sim.N_snaps - 1].satellites['WM'].values * M_S / M_E)
        SWM += list(sim.snaps[sim.N_snaps - 1].satellites['SWM'].values * M_S / M_E)
        System += [key for i in sim.snaps[sim.N_snaps - 1].satellites['M'].values]
    data = zip(Masses, Orb_Dist, WM, SWM, System)

    data = [item for item in data if item[0] >= m_low_lim and item[1] <= a_up_lim]

    Masses, Orb_Dist, WM, SWM, System = zip(*data)
    Masses = np.array(Masses)
    WM = np.array(WM)
    SWM = np.array(SWM)

    WMF = WM / Masses
    SWMF = SWM / Masses
    TWMF = WMF + SWMF

    plt.rcParams.update({'figure.autolayout': True})
    plt.style.use('seaborn-paper')
    plt.rcParams.update({'font.size': pop.fontsize})

    N_bins = 15
    bins = 10 ** np.linspace(np.log10(min(Orb_Dist)), np.log10(max(Orb_Dist)), N_bins)
    fig, ax = plt.subplots(figsize=pop.figsize)
    # ax.hist(Masses, bins=bins)
    values, base, _ = plt.hist(Orb_Dist, bins=bins, rwidth=0.95, density=True, weights=TWMF / np.sum(TWMF))
    ax.axvline(1, color='red', linewidth=1)
    ax.axvline(0.387, color='red', linewidth=1)
    ax.axvline(0.732, color='red', linewidth=1)
    ax.axvline(1.52, color='red', linewidth=1)
    ax_bis = ax.twinx()
    # values, base = np.histogram(Masses, bins = 10 ** np.linspace(np.log10(min(Masses)), np.log10(max(Masses)), N_bins))

    values = np.append(0, values)
    ax_bis.plot(base, np.cumsum(values) / np.cumsum(values)[-1], color='black', linestyle='dashed', markersize=0.1)
    ax.set(xlabel=r'a [$\mathrm{au}$]', ylabel=r'Weighted Counts')
    ax_bis.set(ylabel='Cumulative Distribution')
    ax.set_xscale('log')
    ax_bis.set_xscale('log')
    if pop.plot_config == 'presentation':
        ax.set(title=r'TWMF Weighted Histrogram of Terrestrial Planets Orbital Distances')
    save_name = 'histogram_a_twmf'
    if a_up_lim < 30 and m_low_lim > 0:
        save_name += '_lim'
    fig.savefig(path.join(pop.PLOT, save_name + '.png'), transparent=False, dpi=pop.dpi, bbox_inches="tight")
    plt.close(fig)


def histogram_a_wmf(pop, m_low_lim=0, a_up_lim=30):
    Masses = []
    Orb_Dist = []
    WM = []
    SWM = []
    System = []

    for key, sim in pop.SIMS.items():
        Masses += list(sim.snaps[sim.N_snaps - 1].satellites['M'].values * M_S / M_E)
        Orb_Dist += list(sim.snaps[sim.N_snaps - 1].satellites['a'].values * R_S / au)
        WM += list(sim.snaps[sim.N_snaps - 1].satellites['WM'].values * M_S / M_E)
        SWM += list(sim.snaps[sim.N_snaps - 1].satellites['SWM'].values * M_S / M_E)
        System += [key for i in sim.snaps[sim.N_snaps - 1].satellites['M'].values]
    data = zip(Masses, Orb_Dist, WM, SWM, System)

    data = [item for item in data if item[0] >= m_low_lim and item[1] <= a_up_lim]

    Masses, Orb_Dist, WM, SWM, System = zip(*data)
    Masses = np.array(Masses)
    WM = np.array(WM)
    SWM = np.array(SWM)

    WMF = WM / Masses
    SWMF = SWM / Masses
    TWMF = WMF + SWMF

    plt.rcParams.update({'figure.autolayout': True})
    plt.style.use('seaborn-paper')
    plt.rcParams.update({'font.size': pop.fontsize})

    N_bins = 15
    bins = 10 ** np.linspace(np.log10(min(Orb_Dist)), np.log10(max(Orb_Dist)), N_bins)
    fig, ax = plt.subplots(figsize=pop.figsize)
    # ax.hist(Masses, bins=bins)
    values, base, _ = plt.hist(Orb_Dist, bins=bins, rwidth=0.95, density=True, weights=WMF)
    ax.axvline(1, color='red', linewidth=1)
    ax.axvline(0.387, color='red', linewidth=1)
    ax.axvline(0.732, color='red', linewidth=1)
    ax.axvline(1.52, color='red', linewidth=1)
    ax_bis = ax.twinx()
    # values, base = np.histogram(Masses, bins = 10 ** np.linspace(np.log10(min(Masses)), np.log10(max(Masses)), N_bins))

    values = np.append(0, values)
    ax_bis.plot(base, np.cumsum(values) / np.cumsum(values)[-1], color='black', linestyle='dashed', markersize=0.1)
    ax.set(xlabel=r'a [$\mathrm{au}$]', ylabel=r'Weighted Counts')
    ax_bis.set(ylabel='Cumulative Distribution')
    ax.set_xscale('log')
    ax_bis.set_xscale('log')
    if pop.plot_config == 'presentation':
        ax.set(title=r'WMF Weighted Histrogram of Terrestrial Planets Orbital Distances')
    save_name = 'histogram_a_wmf'
    if a_up_lim < 30 and m_low_lim > 0:
        save_name += '_lim'
    fig.savefig(path.join(pop.PLOT, save_name + '.png'), transparent=False, dpi=pop.dpi, bbox_inches="tight")
    plt.close(fig)


def histogram_wmf(pop, m_low_lim=0, a_up_lim=30):
    Masses = []
    Orb_Dist = []
    WM = []
    SWM = []

    for sim in pop.SIMS.values():
        Masses += list(sim.snaps[sim.N_snaps - 1].satellites['M'].values * M_S / M_E)
        Orb_Dist += list(sim.snaps[sim.N_snaps - 1].satellites['a'].values * R_S / au)
        WM += list(sim.snaps[sim.N_snaps - 1].satellites['WM'].values * M_S / M_E)
        SWM += list(sim.snaps[sim.N_snaps - 1].satellites['SWM'].values * M_S / M_E)

    data = zip(Masses, Orb_Dist, WM, SWM)

    data = [item for item in data if item[0] >= m_low_lim and item[1] <= a_up_lim]

    Masses, Orb_Dist, WM, SWM = zip(*data)
    Masses = np.array(Masses)
    WM = np.array(WM)
    SWM = np.array(SWM)

    WMF = WM / Masses
    SWMF = SWM / Masses

    plt.rcParams.update({'figure.autolayout': True})
    plt.style.use('seaborn-paper')
    plt.rcParams.update({'font.size': pop.fontsize})

    # N_bins = 15
    # bins = 10 ** np.linspace(np.log10(min(Orb_Dist)), np.log10(max(Orb_Dist)), N_bins)
    fig, ax = plt.subplots(figsize=pop.figsize)
    # ax.hist(Masses, bins=bins)
    # values, base, _ = plt.hist(Orb_Dist, bins=bins, rwidth=0.95)
    values, base, _ = plt.hist(WMF / Masses, bins=15, rwidth=0.95, alpha=0.5)
    values, base, _ = plt.hist(SWMF / Masses, bins=15, rwidth=0.95, alpha=0.5)
    # ax.axvline(1, color='red', linewidth=1)
    # ax.axvline(0.387, color='red', linewidth=1)
    # ax.axvline(0.732, color='red', linewidth=1)
    # ax.axvline(1.52, color='red', linewidth=1)
    # ax_bis = ax.twinx()
    # values, base = np.histogram(Masses, bins = 10 ** np.linspace(np.log10(min(Masses)), np.log10(max(Masses)), N_bins))

    # values = np.append(0, values)
    # ax_bis.plot(base, np.cumsum(values) / np.cumsum(values)[-1], color='black', linestyle='dashed', markersize=0.1)
    ax.set(xlabel=r'Water Mass Fraction', ylabel=r'Counts')
    # ax_bis.set(ylabel='Cumulative Distribution')
    # ax.set_xscale('log')
    # ax_bis.set_xscale('log')
    if pop.plot_config == 'presentation':
        ax.set(title=r'Histrogram of Terrestrial Planets Orbital Distances')
    save_name = 'histogram_wmf'
    if a_up_lim < 30 and m_low_lim > 0:
        save_name += '_lim'
    fig.savefig(path.join(pop.PLOT, save_name + '.png'), transparent=False, dpi=pop.dpi, bbox_inches="tight")
    fig.savefig(path.join(pop.PLOT, 'histogram_wmf.png'), transparent=False, dpi=pop.dpi, bbox_inches="tight")
    plt.close(fig)



def histogram_totalmass(pop, a_up_lim=30, m_low_lim=0):
    TotalMasses = []

    for sim in pop.SIMS.values():
        TotalMasses.append(np.sum(
            [item for item in list(sim.snaps[sim.N_snaps - 1].satellites['M'].values * M_S / M_E) if
             item >= m_low_lim]))

    plt.rcParams.update({'figure.autolayout': True})
    plt.style.use('seaborn-paper')
    plt.rcParams.update({'font.size': pop.fontsize})

    N_bins = 15
    bins = 10 ** np.linspace(np.log10(min(TotalMasses)), np.log10(max(TotalMasses)), N_bins)
    fig, ax = plt.subplots(figsize=pop.figsize)
    # ax.hist(Masses, bins=bins)
    values, base, _ = plt.hist(TotalMasses, bins=bins, rwidth=0.95)
    # ax.axvline(1, color='red', linewidth=1)

    ax_bis = ax.twinx()
    # values, base = np.histogram(Masses, bins = 10 ** np.linspace(np.log10(min(Masses)), np.log10(max(Masses)), N_bins))

    values = np.append(0, values)
    ax_bis.plot(base, np.cumsum(values) / np.cumsum(values)[-1], color='black', linestyle='dashed', markersize=0.1)
    ax.set(xlabel=r'Total Mass [$\mathrm{M_{\oplus}}$]', ylabel=r'Counts')
    ax_bis.set(ylabel='Cumulative Distribution')
    ax.set_xscale('log')
    ax_bis.set_xscale('log')
    if pop.plot_config == 'presentation':
        ax.set(title=r'Histrogram of Total Masses')
    save_name = 'histogram_totalmass'
    if a_up_lim < 30 and m_low_lim > 0:
        save_name += '_lim'
    fig.savefig(path.join(pop.PLOT, save_name + '.png'), transparent=False, dpi=pop.dpi, bbox_inches="tight")
    plt.close(fig)


def histogram_totalmass_thresh(pop, a_up_lim=30, thresholds=None):
    if thresholds is None:
        thresholds = [0.0, 0.1, 1.0, 10]

    thresholds = sorted(thresholds)

    systems = []
    for sim in pop.SIMS.values():
        zipped = zip(sim.snaps[sim.N_snaps - 1].satellites['M'].values * M_S / M_E,
                     sim.snaps[sim.N_snaps - 1].satellites['a'].values * R_S / au)
        zipped = [(m, a) for (m, a) in zipped if a <= a_up_lim]
        masses, orb_dist = zip(*zipped)
        systems.append(np.array(masses))

    total_masses = {key: [] for key in thresholds}
    for keys, lists in total_masses.items():
        for i, item in enumerate(systems):
            filtered = item[np.where(item > keys)]
            lists.append(np.sum(filtered))

    solar_sys = [m / M_E for (m, a) in terrestrial if a / au <= a_up_lim]
    total_mass_solar_sys = np.sum(solar_sys)

    plt.rcParams.update({'figure.autolayout': True})
    plt.style.use('seaborn-paper')
    plt.rcParams.update({'font.size': pop.fontsize})
    plt.rcParams.update({"legend.title_fontsize": pop.legend_fontsize})

    N_bins = 15
    bins = 10 ** np.linspace(np.log10(min(total_masses[thresholds[0]])), np.log10(max(total_masses[thresholds[0]])),
                             N_bins)
    fig, ax = plt.subplots(figsize=pop.figsize)
    for key, lists in total_masses.items():
        label = f'{key}' + r' $\mathrm{M_{\oplus}}$'
        ax.hist(lists, bins=bins, rwidth=0.95, label=label)

    ax.axvline(total_mass_solar_sys, color='red', linewidth=1)

    ax.set(xlabel=r'Total Mass [$\mathrm{M_{\oplus}}$]', ylabel=r'Counts')
    ax.set_xscale('log')
    ax.legend(title='Thresholds', frameon=False)
    #plt.legend(frameon=False)
    if pop.plot_config == 'presentation':
        ax.set(title=r'Histrogram of Total Masses')
    save_name = 'histogram_totalmass_thresh'
    if a_up_lim < 30:
        save_name += '_lim'
    fig.savefig(path.join(pop.PLOT, save_name + '.png'), transparent=False, dpi=pop.dpi, bbox_inches="tight")
    plt.close(fig)

def histogram_a_thresh(pop, a_up_lim=30, thresholds=None):
    if thresholds is None:
        thresholds = [0.0, 0.5, 1.0]

    thresholds = sorted(thresholds)

    Masses = []
    Orb_Dist = []
    System = []
    for key, sim in pop.SIMS.items():
        Masses += list(sim.snaps[sim.N_snaps - 1].satellites['M'].values * M_S / M_E)
        Orb_Dist += list(sim.snaps[sim.N_snaps - 1].satellites['a'].values * R_S / au)
        System += [key for i in sim.snaps[sim.N_snaps - 1].satellites['M'].values]

    data = zip(Masses, Orb_Dist, System)
    data = [item for item in data if item[1] <= a_up_lim]

    Orb_Dict = {}
    Weights = {}
    for thresh in thresholds:
            Orb_Dict[thresh] = np.array([item[1] for item in data if item[0] >= thresh])
            print(np.array([item[-1] for item in data if item[0] >= thresh]))
            Weights[thresh] = np.size(np.unique(np.array([item[-1] for item in data if item[0] >= thresh])))
    solar_sys = [m / M_E for (m, a) in terrestrial if a / au <= a_up_lim]
    total_mass_solar_sys = np.sum(solar_sys)

    plt.rcParams.update({'figure.autolayout': True})
    plt.style.use('seaborn-paper')
    plt.rcParams.update({'font.size': pop.fontsize})
    plt.rcParams.update({"legend.title_fontsize": pop.legend_fontsize})

    N_bins = 15
    bins = 10 ** np.linspace(np.log10(min(Orb_Dict[thresholds[0]])), np.log10(max(Orb_Dict[thresholds[0]])),
                             N_bins)
    fig, ax = plt.subplots(figsize=pop.figsize)
    stacked = []
    for thr, arr in Orb_Dict.items():
        label = f'{thr}' + r' $\mathrm{M_{\oplus}}$'
        stacked.append(arr)
        ax.hist(arr, bins=bins, rwidth=0.95, label=label, alpha=1, weights=np.full_like(arr,1/Weights[thr]))
    #ax.hist(stacked,bins=bins, rwidth=0.95, label=Orb_Dict.keys(), alpha=0.4)
    ax.axvline(1, color='red', linewidth=1)
    ax.axvline(0.387, color='red', linewidth=1)
    ax.axvline(0.732, color='red', linewidth=1)
    ax.axvline(1.52, color='red', linewidth=1)

    ax.set(xlabel=r'Orbital Distance [$\mathrm{au}$]', ylabel=r'Counts')
    ax.set_xscale('log')
    ax.legend(title='Thresholds', frameon=False)
    plt.legend(frameon=False)
    if pop.plot_config == 'presentation':
        ax.set(title=r'Histrogram of Total Masses')
    save_name = 'histogram_a_thresh'
    if a_up_lim < 30:
        save_name += '_lim'
    fig.savefig(path.join(pop.PLOT, save_name + '.png'), transparent=False, dpi=pop.dpi, bbox_inches="tight")
    plt.close(fig)