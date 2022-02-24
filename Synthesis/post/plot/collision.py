import matplotlib.pyplot as plt
import os.path as path

import numpy as np

from Synthesis.units import *


def collisions_time(pop, t_low_lim=0):
    times = []
    for sim in pop.SIMS.values():
        times.extend(sim.collisions.index)
    times = np.array(times)
    times = times[np.where(times > t_low_lim)]
    print(times)
    times /= 1e6
    plt.rcParams.update({'figure.autolayout': True})
    plt.style.use('seaborn-paper')
    plt.rcParams.update({'font.size': pop.fontsize})
    N_bins = 15
    fig, ax = plt.subplots(figsize=pop.figsize)
    values, base, _ = plt.hist(times, bins=N_bins, rwidth=0.95, log=True)
    ax_bis = ax.twinx()
    values = np.append(0, values)
    ax_bis.plot(base, np.cumsum(values) / np.cumsum(values)[-1], color='black', linestyle='dashed', markersize=0.1)
    ax.set(xlabel=r'Time [$\mathrm{Myr}$]', ylabel=r'Counts')
    ax_bis.set(ylabel='Cumulative Distribution')
    if pop.plot_config == 'presentation':
        ax.set(title=r'Histrogram of Collision Times')
    save_name = 'histogram_col_time'
    if t_low_lim > 0:
        save_name += '_lim'
    fig.savefig(path.join(pop.PLOT, save_name + '.png'), transparent=False, dpi=pop.dpi, bbox_inches="tight")
    plt.close(fig)


def collisions_weighted_time(pop, t_low_lim = 0):
    times = []
    primary_masses = []
    secondary_masses = []
    final_masses = []
    wm2 = []
    for sim in pop.SIMS.values():
        times.extend(sim.collisions.index)
        primary_masses.extend(sim.collisions['mass1'] * M_S / M_E)
        secondary_masses.extend(sim.collisions['mass2'] * M_S / M_E)
        final_masses.extend(sim.collisions['mass'] * M_S / M_E)
        wm2.extend(sim.collisions['wm2'] * M_S / M_E)

    zipped = zip(times, primary_masses, secondary_masses, final_masses, wm2)
    filtered = [item for item in zipped if item[0] > t_low_lim]

    times, primary_masses, secondary_masses, final_masses, wm2 = zip(*filtered)

    zipped = zip(primary_masses,secondary_masses)
    colliders = [np.abs(tup[0] - tup[1]) for tup in zipped]
    colliders = np.array(colliders) / np.array(secondary_masses)
    weights = colliders / np.sum(colliders)

    times = np.array(times) / 1e6
    plt.rcParams.update({'figure.autolayout': True})
    plt.style.use('seaborn-paper')
    plt.rcParams.update({'font.size': pop.fontsize})
    N_bins = 15
    fig, ax = plt.subplots(figsize=pop.figsize)
    values, base, _ = plt.hist(times, bins=N_bins, rwidth=0.95, density=True)
    ax_bis = ax.twinx()
    values = np.append(0, values)
    ax_bis.plot(base, np.cumsum(values) / np.cumsum(values)[-1], color='black', linestyle='dashed', markersize=0.1)
    ax.set(xlabel=r'Time [$\mathrm{Myr}$]', ylabel=r'Weighted counts')
    ax_bis.set(ylabel='Cumulative Distribution')
    if pop.plot_config == 'presentation':
        ax.set(title=r'Histrogram of Collision Times')
    save_name = 'histogram_weighted_col_time'
    if t_low_lim > 0:
        save_name += '_lim'
    fig.savefig(path.join(pop.PLOT, save_name + '.png'), transparent=False, dpi=pop.dpi, bbox_inches="tight")
    plt.close(fig)


def collisions_mass(pop):
    times = []
    primary_masses = []
    secondary_masses = []
    final_masses = []
    wm2 = []
    for sim in pop.SIMS.values():
        times.extend(sim.collisions.index)
        primary_masses.extend(sim.collisions['mass1'] * M_S / M_E)
        secondary_masses.extend(sim.collisions['mass2'] * M_S / M_E)
        final_masses.extend(sim.collisions['mass'] * M_S / M_E)
        wm2.extend(sim.collisions['wm2'] * M_S / M_E)

    secondary_masses = np.array(secondary_masses)
    primary_masses = np.array(primary_masses)
    final_masses = np.array(final_masses)
    lists = [primary_masses,secondary_masses]

    plt.rcParams.update({'figure.autolayout': True})
    plt.style.use('seaborn-paper')
    plt.rcParams.update({'font.size': pop.fontsize})

    N_bins = 15
    bins = 10 ** np.linspace(np.log10(min(secondary_masses)), np.log10(max(primary_masses)), N_bins)
    fig, ax = plt.subplots(figsize=pop.figsize)
    #values, base, _ = plt.hist(stacked, bins=bins, rwidth=0.95, stacked=True)
    ax.hist(lists, bins=bins, rwidth=0.95, alpha=1, stacked=True)
    #ax.hist(primary_masses, bins=bins, rwidth=0.95, alpha=0.5)
    #ax.hist(secondary_masses, bins=bins, rwidth=0.95, alpha=0.5)
    #plt.hist(data2, bins=bins, rwidth=0.95, alpha=0.5, label='Secondary Masses')
    #ax_bis = ax.twinx()
    #values = np.append(0, values)
    #ax_bis.plot(base, np.cumsum(values) / np.cumsum(values)[-1], color='black', linestyle='dashed', markersize=0.1)
    ax.set(xlabel=r'Mass [$\mathrm{M_{\oplus}}$]', ylabel=r'Counts')
    #ax_bis.set(ylabel='Cumulative Distribution')
    ax.set_xscale('log')
    #ax_bis.set_xscale('log')
    if pop.plot_config == 'presentation':
        ax.set(title=r'Histrogram of Collision Times')
    save_name = 'histogram_col_mass'
    fig.savefig(path.join(pop.PLOT, save_name + '.png'), transparent=False, dpi=pop.dpi, bbox_inches="tight")
    plt.close(fig)


def collisions_orb_dist(pop):
    xs = []
    ys = []
    zs = []
    for sim in pop.SIMS.values():
        xs.extend(sim.collisions['x'].values)
        ys.extend(sim.collisions['y'].values)
        zs.extend(sim.collisions['z'].values)

    xs = np.array(xs)
    ys = np.array(ys)
    zs = np.array(zs)

    orb_dist = np.sqrt(np.power(xs, 2) + np.power(ys, 2) + np.power(zs, 2)) * R_S / au
    plt.rcParams.update({'figure.autolayout': True})
    plt.style.use('seaborn-paper')
    plt.rcParams.update({'font.size': pop.fontsize})
    N_bins = 15
    bins = 10 ** np.linspace(np.log10(min(orb_dist)), np.log10(max(orb_dist)), N_bins)
    fig, ax = plt.subplots(figsize=pop.figsize)
    values, base, _ = plt.hist(orb_dist, bins=bins, rwidth=0.95)
    ax_bis = ax.twinx()
    values = np.append(0, values)
    ax_bis.plot(base, np.cumsum(values) / np.cumsum(values)[-1], color='black', linestyle='dashed', markersize=0.1)
    ax.set(xlabel=r'Orbital distance [$\mathrm{au}$]', ylabel=r'Counts')
    ax_bis.set(ylabel='Cumulative Distribution')
    ax.set_xscale('log')
    ax_bis.set_xscale('log')
    if pop.plot_config == 'presentation':
        ax.set(title=r'Histrogram of Collision Places')
    save_name = 'histogram_col_orb_dist'
    fig.savefig(path.join(pop.PLOT, save_name + '.png'), transparent=False, dpi=pop.dpi, bbox_inches="tight")
    plt.close(fig)


def collisions_engulfed(pop):
    lost_mass = []
    for id, sim in pop.SIMS.items():
        masses = sim.lost_satellites['mass'].values
        cols = sim.lost_satellites['collision'].values
        filter_null = cols == 0.0
        summed = np.sum(masses[filter_null])
        lost_mass.append(summed)
    lost_mass = np.array(lost_mass) * M_S / M_J
    filter_null = lost_mass > 0
    filter_not_null = lost_mass == 0
    non_lost_mass = lost_mass[filter_not_null]
    lost_mass = lost_mass[filter_null]
    print(f'Number of System with no mass loss: {len(non_lost_mass)}')

    plt.rcParams.update({'figure.autolayout': True})
    plt.style.use('seaborn-paper')
    plt.rcParams.update({'font.size': pop.fontsize})

    N_bins = 15
    bins = 10 ** np.linspace(np.log10(min(lost_mass)), np.log10(max(lost_mass)), N_bins)
    fig, ax = plt.subplots(figsize=pop.figsize)
    values, base, _ = plt.hist(lost_mass, bins=bins, rwidth=0.95)
    ax_bis = ax.twinx()
    values = np.append(0, values)
    ax_bis.plot(base, np.cumsum(values) / np.cumsum(values)[-1], color='black', linestyle='dashed', markersize=0.1)
    ax.set(xlabel=r'Mass [$\mathrm{M_J}$]', ylabel=r'Counts')
    ax_bis.set(ylabel='Cumulative Distribution')
    ax.set_xscale('log')
    ax_bis.set_xscale('log')
    if pop.plot_config == 'presentation':
        ax.set(title=r'Histrogram of Lost Mass')
    save_name = 'histogram_lost_mass'
    fig.savefig(path.join(pop.PLOT, save_name + '.png'), transparent=False, dpi=pop.dpi, bbox_inches="tight")
    plt.close(fig)
