import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib import colors
import os.path as path
from Synthesis.units import *


def scatter_parameters(pop):
    TotalMasses = []
    SigmaCoeffs = []
    Reference = []

    for sim in pop.SIMS.values():
        TotalMasses.append(sim.Total_Mass)
        SigmaCoeffs.append(sim.Sigma_Exponent)
        Reference.append(sim.Sigma_Norm * pow(R_S / au, sim.Sigma_Exponent) * (M_S / R_S2))

    plt.rcParams.update({'figure.autolayout': True})
    plt.style.use('seaborn-paper')
    plt.rcParams.update({'font.size': pop.fontsize})

    cmap = pop.cmap_standart
    cmin = min(Reference)
    cmax = max(Reference)

    norm = colors.LogNorm(cmin, cmax)

    fig, ax = plt.subplots(figsize=pop.figsize)
    ax.scatter(SigmaCoeffs, TotalMasses, c=Reference, cmap=cmap, norm=norm, s=12)
    x_labels = ax.get_xticklabels()
    plt.setp(x_labels, horizontalalignment='center')
    ax.set(xlabel='Surface Density Power Law Exponent', ylabel=r'Total Disk Mass [$M_{\odot}$]', xticks=SigmaCoeffs)
    ax2 = ax.twinx()
    mn, mx = ax.get_ylim()
    ax2.set_ylim(M_S / M_J * mn, M_S / M_J * mx)
    ax2.set_ylabel('Total Disk Mass [$M_{J}$]')
    fig.colorbar(cm.ScalarMappable(cmap=cmap, norm=norm), orientation='vertical',
                 label=r'Reference Value at $1 \mathrm{au}$ [$\mathrm{g}\mathrm{cm}^{-2}$]', ax=ax2, pad=0.12)
    # ax.set_yscale('log')
    if pop.plot_config == 'presentation':
        ax.set(title=r'Synthesis Parameters')

    fig.savefig(path.join(pop.PLOT, 'scatter_parameters.png'), transparent=False, dpi=pop.dpi,
                bbox_inches="tight")
    plt.close(fig)

def scatter_parameters_numbers(pop, m_low_lim=0, a_up_lim=30):
    TotalMasses = []
    SigmaCoeffs = []
    Reference = []
    Masses = []
    Orb_Dist = []
    Numbers = []
    Systems = []

    for sim in pop.SIMS.values():
        TotalMasses.append(sim.Total_Mass)
        SigmaCoeffs.append(sim.Sigma_Exponent)
        Masses = list(sim.snaps[sim.N_snaps - 1].satellites['M'].values * M_S / M_E)
        Orb_Dist = list(sim.snaps[sim.N_snaps - 1].satellites['a'].values * R_S / au)
        system = zip(Masses, Orb_Dist)
        filtered = [item for item in system if item[0] >= m_low_lim and item[1] <= a_up_lim]
        Numbers.append(len(filtered))
    print(Numbers)
    plt.rcParams.update({'figure.autolayout': True})
    plt.style.use('seaborn-paper')
    plt.rcParams.update({'font.size': pop.fontsize})
    cmap = pop.cmap_standart
    cmin = min(Numbers)
    cmax = max(Numbers)

    norm = colors.Normalize(cmin, cmax)

    fig, ax = plt.subplots(figsize=pop.figsize)
    ax.scatter(SigmaCoeffs, TotalMasses, c=Numbers, cmap=cmap, norm=norm, s=12)
    x_labels = ax.get_xticklabels()
    plt.setp(x_labels, horizontalalignment='center')
    ax.set(xlabel='Surface Density Power Law Exponent', ylabel=r'Total Disk Mass [$M_{\odot}$]', xticks=SigmaCoeffs)
    ax2 = ax.twinx()
    mn, mx = ax.get_ylim()
    ax2.set_ylim(M_S / M_J * mn, M_S / M_J * mx)
    ax2.set_ylabel('Total Disk Mass [$M_{J}$]')
    fig.colorbar(cm.ScalarMappable(cmap=cmap, norm=norm), orientation='vertical',
                 label=r'Number of Planets', ax=ax2, pad=0.12)
    # ax.set_yscale('log')
    if pop.plot_config == 'presentation':
        ax.set(title=r'Synthesis Parameters')

    fig.savefig(path.join(pop.PLOT, 'scatter_parameters_numbers.png'), transparent=False, dpi=pop.dpi,
                bbox_inches="tight")
    plt.close(fig)


def scatter_ecc_inc(pop, m_low_lim=0, a_up_lim=30):
    Masses = []
    Orb_Dist = []
    Ecc = []
    Inc = []

    for sim in pop.SIMS.values():
        Masses += list(sim.snaps[sim.N_snaps - 1].satellites['M'].values * M_S / M_E)
        Orb_Dist += list(sim.snaps[sim.N_snaps - 1].satellites['a'].values * R_S / au)
        Ecc += list(sim.snaps[sim.N_snaps - 1].satellites['e'].values)
        Inc += list(sim.snaps[sim.N_snaps - 1].satellites['i'].values)

    data = zip(Masses, Orb_Dist, Ecc, Inc)

    data = [item for item in data if item[0] >= m_low_lim and item[1] <= a_up_lim]

    Masses, Orb_Dist, Ecc, Inc = zip(*data)

    print(f'Number of Object: {len(Masses)}')

    plt.rcParams.update({'figure.autolayout': True})
    plt.style.use('seaborn-paper')
    plt.rcParams.update({'font.size': pop.fontsize})
    plt.rcParams.update({"legend.title_fontsize": pop.legend_fontsize})

    cmap = pop.cmap_standart
    cmin = min(Orb_Dist)
    cmax = max(Orb_Dist)

    norm = colors.LogNorm(cmin, cmax)

    fig, ax = plt.subplots(figsize=pop.figsize)
    ax.scatter(Ecc, Inc, c=Orb_Dist, cmap=cmap, norm=norm, s=3)
    x_labels = ax.get_xticklabels()
    plt.setp(x_labels, horizontalalignment='center')
    ax.set(xlabel='Eccentricity', ylabel=r'$\sin(\mathrm{inclination})$')
    fig.colorbar(cm.ScalarMappable(cmap=cmap, norm=norm), orientation='vertical',
                 label=r'Orbital Distance [$\mathrm{au}$]',
                 ax=ax)
    ax.set_xscale('log')
    ax.set_yscale('log')
    if pop.plot_config == 'presentation':
        ax.set(title=r'Eccentricity and Inclination')
    save_name = 'scatter_ecc_inc'
    if a_up_lim < 30 and m_low_lim > 0:
        save_name += '_lim'
    fig.savefig(path.join(pop.PLOT, save_name + '.png'), transparent=False, dpi=pop.dpi, bbox_inches="tight")
    plt.close(fig)

def scatter_a_mass(pop, m_low_lim=0, a_up_lim=30):
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

    data = [(m,a,wm/m,swm/m) for (m,a,wm,swm) in data if m >= m_low_lim and a <= a_up_lim]

    Masses, Orb_Dist, WMF, SWMF = zip(*data)

    TWMF = np.array(WMF) + np.array(SWMF)

    print(f'Number of Object: {len(Masses)}')

    plt.rcParams.update({'figure.autolayout': True})
    plt.style.use('seaborn-paper')
    plt.rcParams.update({'font.size': pop.fontsize})
    plt.rcParams.update({"legend.title_fontsize": pop.legend_fontsize})

    cmap = pop.cmap_standart
    cmin = min(TWMF)
    cmax = max(TWMF)

    norm = colors.Normalize(cmin, cmax)

    fig, ax = plt.subplots(figsize=pop.figsize)
    ax.scatter(Orb_Dist, Masses, c=TWMF, cmap=cmap, norm=norm, s=3)
    x_labels = ax.get_xticklabels()
    plt.setp(x_labels, horizontalalignment='center')
    ax.set(xlabel=r'Orbital Distance [$\mathrm{au}$]', ylabel=r'Mass [$\mathrm{M_{\oplus}}$]')
    fig.colorbar(cm.ScalarMappable(cmap=cmap, norm=norm), orientation='vertical',
                 label=r'Total WMF',
                 ax=ax)
    ax.set_xscale('log')
    ax.set_yscale('log')
    if pop.plot_config == 'presentation':
        ax.set(title=r'Eccentricity and Inclination')
    save_name = 'scatter_a_mass'
    if a_up_lim < 30 and m_low_lim > 0:
        save_name += '_lim'
    fig.savefig(path.join(pop.PLOT, save_name + '.png'), transparent=False, dpi=pop.dpi, bbox_inches="tight")
    plt.close(fig)


def scatter_radial_twmf(pop, m_low_lim=0, a_up_lim=30):
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

    WMF = np.array(WM) / np.array(Masses)
    SWMF = np.array(WM) / np.array(Masses)
    TWMF = WMF + SWMF

    data = zip(Masses, Orb_Dist, TWMF, System)

    data = [item for item in data if item[0] >= m_low_lim and item[1] <= a_up_lim]

    earth_analogs = [item for item in data if item[0] >= 0.101 and item[1] <= a_up_lim and item[2] > 0.001]
    print(earth_analogs)

    Masses, Orb_Dist, TWMF, System = zip(*data)

    ms, obs, twmf, system = zip(*earth_analogs)

    plt.rcParams.update({'figure.autolayout': True})
    plt.style.use('seaborn-paper')
    plt.rcParams.update({'font.size': pop.fontsize})
    plt.rcParams.update({"legend.title_fontsize": pop.legend_fontsize})

    cmap = pop.cmap_standart
    cmin = min(TWMF)
    cmax = max(TWMF)

    norm = colors.Normalize(cmin, cmax)

    fig, ax = plt.subplots(figsize=pop.figsize)
    ax.scatter(Orb_Dist, Masses, c=TWMF, cmap=cmap, norm=norm, s=5, alpha=0.1)
    # ax.scatter(obs, ms, c=twmf, cmap=cmap, norm=norm, s=10)
    ax.axvline(1, color='black', linewidth=0.7, linestyle='--')
    ax.axhline(1, color='black', linewidth=0.7, linestyle='--')
    x_labels = ax.get_xticklabels()
    plt.setp(x_labels, horizontalalignment='center', )

    ax.set(xlabel='Radial Distance [$\mathrm{au}$]', ylabel=r'Mass [$\mathrm{M_{\oplus}}$]')
    fig.colorbar(cm.ScalarMappable(cmap=cmap, norm=norm), orientation='vertical', label=r'Total WMF', ax=ax)
    ax.set_xscale('log')
    ax.set_yscale('log')
    if pop.plot_config == 'presentation':
        ax.set(title=r'Total WMF Radial Distribution')
        save_name = 'scatter_radial_twmf'
    if a_up_lim < 30 and m_low_lim > 0:
        save_name += '_lim'
    fig.savefig(path.join(pop.PLOT, save_name + '.png'), transparent=False, dpi=pop.dpi, bbox_inches="tight")
    fig.savefig(path.join(pop.PLOT, 'scatter_radial_twmf.png'), transparent=False, dpi=pop.dpi,
                bbox_inches="tight")
    plt.close(fig)

    for earth in earth_analogs:
        print(f'System: {earth[-1]}')
        print(f'Mass: {earth[0]}')
        print(f'Orb Dist: {earth[1]}')
        print(f'TWMF: {earth[2]}')
        print(f'Exponent: {pop.SIMS[earth[-1]].Sigma_Exponent}')
        print(f'Disk Mass: {pop.SIMS[earth[-1]].Total_Mass * M_S / M_J}')
        print(" ")
