import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm
from matplotlib import colors
from matplotlib import patches
import os.path as path
from Synthesis.units import *
from tqdm import tqdm
from scipy.integrate import quad

def Power_Law(x, a, b):
    return a * np.power(x, b)

def scatter_parameters(pop):
    TotalMasses = []
    SigmaCoeffs = []
    Reference = []

    for sim in pop.SIMS.values():
        TotalMasses.append(sim.Total_Mass)
        SigmaCoeffs.append(sim.Sigma_Exponent)
        print(sim.Sigma_Exponent)
        print(sim.Sigma_Norm * (R_S / au)**sim.Sigma_Exponent / denstos)
        Reference.append(sim.Sigma_Norm / (R_S / au)**sim.Sigma_Exponent / denstos * pow(au/R_S, sim.Sigma_Exponent) / denstos)


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
    ax.set(xlabel='Surface Density Power Law Exponent', ylabel=r'Total Mass [$M_{\odot}$]', xticks=SigmaCoeffs)
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
    Means = []
    Systems = []

    for id,sim in pop.SIMS.items():
        TotalMasses.append(sim.Total_Mass)
        SigmaCoeffs.append(sim.Sigma_Exponent)
        Masses = list(sim.snaps[sim.N_snaps - 1].satellites['M'].values * M_S / M_E)
        Orb_Dist = list(sim.snaps[sim.N_snaps - 1].satellites['a'].values * R_S / au)
        system = zip(Masses, Orb_Dist)
        filtered = [item for item in system if item[0] >= m_low_lim and item[1] <= a_up_lim]
        # mean = np.max([item[0] for item in filtered])/np.sum([item[0] for item in filtered])
        # Means.append(mean)
        Numbers.append(len(filtered))
    #Means = np.array(Means) / np.sum(Means)
    print(Numbers)
    Numbers = np.array(Numbers)
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

def scatter_parameters_lost_mass(pop, m_low_lim=0, a_up_lim=30):
    TotalMasses = []
    SigmaCoeffs = []
    lost_mass = []
    numbers = []
    Reference = []
    TM = []
    for id, sim in pop.SIMS.items():
        TotalMasses.append(sim.Total_Mass)
        SigmaCoeffs.append(sim.Sigma_Exponent)
        TM.append(np.sum([item for item in list(sim.snaps[sim.N_snaps - 1].satellites['M'].values * M_S / M_E)]))
        Reference.append(sim.Sigma_Norm / (R_S / au)**sim.Sigma_Exponent / denstos * pow(au/R_S, sim.Sigma_Exponent) / denstos)
        masses = sim.lost_satellites['mass'].values * M_S / M_J
        cols = sim.lost_satellites['collision'].values
        filter_null = cols == 0.0
        filtered = masses[filter_null]
        summed = np.sum(filtered)
        numbers.append(len(filtered))
        # summed = np.sum(filtered)
        lost_mass.append(summed)

    lost_mass = np.array(lost_mass)

    #Means = np.array(Means) / np.sum(Means)
    # print(Numbers / np.array(Means))
    Numbers = np.array(SigmaCoeffs)
    plt.rcParams.update({'figure.autolayout': True})
    plt.style.use('seaborn-paper')
    plt.rcParams.update({'font.size': pop.fontsize})

    # arr = np.unique(SigmaCoeffs)
    cmap = plt.get_cmap(pop.cmap_standart,len(SigmaCoeffs))

    norm = colors.BoundaryNorm(np.linspace(-1.625, -0.375, len(np.unique(SigmaCoeffs))+1, endpoint=True), cmap.N)



    fig, ax = plt.subplots(figsize=pop.figsize)
    ax.scatter(Reference, lost_mass, c=SigmaCoeffs, cmap=cmap, norm=norm, s=12)
    x_labels = ax.get_xticklabels()
    plt.setp(x_labels, horizontalalignment='center')
    ax.set(ylabel='Total Lost Mass [$\mathrm{M_J}$]', xlabel=r'Reference Value at 1 $\mathrm{au}$ [$\mathrm{g cm^{-2}}$]')
    # ax2 = ax.twinx()
    # mn, mx = ax.get_ylim()
    # ax2.set_ylim(M_S / M_J * mn, M_S / M_J * mx)
    # ax2.set_ylabel('Total Disk Mass [$M_{J}$]')
    fig.colorbar(cm.ScalarMappable(cmap=cmap,norm=norm), orientation='vertical',
                 label=r'Power-Law Exponent', ax=ax, ticks=np.unique(SigmaCoeffs))
    ax.set_yscale('log')
    ax.set_xscale('log')
    if pop.plot_config == 'presentation':
        ax.set(title=r'Synthesis Parameters')

    fig.savefig(path.join(pop.PLOT, 'scatter_reference_lost_mass.png'), transparent=False, dpi=pop.dpi,
                bbox_inches="tight")
    plt.close(fig)



def scatter_parameters_AMD(pop, m_low_lim=0, a_up_lim=30):
    TotalMasses = []
    SigmaCoeffs = []
    Reference = []
    Masses = []
    Orb_Dist = []
    Numbers = []
    Means = []
    Systems = []
    AMDS = []

    for sim in tqdm(pop.SIMS.values()):
        TotalMasses.append(sim.Total_Mass)
        SigmaCoeffs.append(sim.Sigma_Exponent)
        AMD, N = sim.get_AMD(m_low_lim, a_up_lim)

        AMDS.append(AMD)

    Means = np.array(Means) / np.sum(Means)
    print(Numbers * Means)
    Numbers = np.array(AMDS)
    plt.rcParams.update({'figure.autolayout': True})
    plt.style.use('seaborn-paper')
    plt.rcParams.update({'font.size': pop.fontsize})
    cmap = pop.cmap_standart
    cmin = min(Numbers)
    cmax = max(Numbers)

    norm = colors.LogNorm(cmin, cmax)

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
                 label=r'AMD', ax=ax2, pad=0.12)
    # ax.set_yscale('log')
    if pop.plot_config == 'presentation':
        ax.set(title=r'Synthesis Parameters')

    fig.savefig(path.join(pop.PLOT, 'scatter_parameters_amd.png'), transparent=False, dpi=pop.dpi,
                bbox_inches="tight")
    plt.close(fig)


def scatter_parameters_RMC(pop, m_low_lim=0, a_up_lim=30):
    TotalMasses = []
    SigmaCoeffs = []
    Reference = []
    Masses = []
    Orb_Dist = []
    Numbers = []
    Means = []
    Systems = []
    RMCS = []

    for sim in tqdm(pop.SIMS.values()):
        TotalMasses.append(sim.Total_Mass)
        SigmaCoeffs.append(sim.Sigma_Exponent)
        RMC, N = sim.get_RMC(m_low_lim, a_up_lim)

        RMCS.append(RMC)

    Means = np.array(Means) / np.sum(Means)
    print(Numbers * Means)
    Numbers = np.array(RMCS)
    plt.rcParams.update({'figure.autolayout': True})
    plt.style.use('seaborn-paper')
    plt.rcParams.update({'font.size': pop.fontsize})
    cmap = pop.cmap_standart
    cmin = min(RMCS)
    cmax = max(RMCS)

    norm = colors.Normalize(cmin, cmax)

    fig, ax = plt.subplots(figsize=pop.figsize)
    ax.scatter(SigmaCoeffs, TotalMasses, c=RMCS, cmap=cmap, norm=norm, s=12)
    x_labels = ax.get_xticklabels()
    plt.setp(x_labels, horizontalalignment='center')
    ax.set(xlabel='Surface Density Power Law Exponent', ylabel=r'Total Disk Mass [$M_{\odot}$]', xticks=SigmaCoeffs)
    ax2 = ax.twinx()
    mn, mx = ax.get_ylim()
    ax2.set_ylim(M_S / M_J * mn, M_S / M_J * mx)
    ax2.set_ylabel('Total Disk Mass [$M_{J}$]')
    fig.colorbar(cm.ScalarMappable(cmap=cmap, norm=norm), orientation='vertical',
                 label=r'RMC', ax=ax2, pad=0.12)
    # ax.set_yscale('log')
    if pop.plot_config == 'presentation':
        ax.set(title=r'Synthesis Parameters')

    fig.savefig(path.join(pop.PLOT, 'scatter_parameters_rmc_nonlog.png'), transparent=False, dpi=pop.dpi,
                bbox_inches="tight")
    plt.close(fig)



def scatter_collision_number(pop, m_low_lim=0, a_up_lim=30):
    TotalMasses = []
    SigmaCoeffs = []
    times = []

    for sim in tqdm(pop.SIMS.values()):
        TotalMasses.append(sim.Total_Mass)
        SigmaCoeffs.append(sim.Sigma_Exponent)
        times.append(len(sim.collisions.index))


    plt.rcParams.update({'figure.autolayout': True})
    plt.style.use('seaborn-paper')
    plt.rcParams.update({'font.size': pop.fontsize})
    cmap = pop.cmap_standart
    cmin = min(times)
    cmax = max(times)

    norm = colors.Normalize(cmin, cmax)

    fig, ax = plt.subplots(figsize=pop.figsize)
    ax.scatter(SigmaCoeffs, TotalMasses, c=times, cmap=cmap, norm=norm, s=12)
    x_labels = ax.get_xticklabels()
    plt.setp(x_labels, horizontalalignment='center')
    ax.set(xlabel='Surface Density Power Law Exponent', ylabel=r'Total Disk Mass [$M_{\odot}$]', xticks=SigmaCoeffs)
    ax2 = ax.twinx()
    mn, mx = ax.get_ylim()
    ax2.set_ylim(M_S / M_J * mn, M_S / M_J * mx)
    ax2.set_ylabel('Total Disk Mass [$M_{J}$]')
    fig.colorbar(cm.ScalarMappable(cmap=cmap, norm=norm), orientation='vertical',
                 label=r'Number of Collisions', ax=ax2, pad=0.12)
    # ax.set_yscale('log')
    if pop.plot_config == 'presentation':
        ax.set(title=r'Synthesis Parameters')

    fig.savefig(path.join(pop.PLOT, 'scatter_parameters_collision_number.png'), transparent=False, dpi=pop.dpi,
                bbox_inches="tight")
    plt.close(fig)


def scatter_ecc_inc(pop, m_low_lim=0, a_up_lim=30):
    Masses = []
    Orb_Dist = []
    Ecc = []
    Inc = []
    Types = []

    for sim in pop.SIMS.values():
        Masses += list(sim.snaps[sim.N_snaps - 1].satellites['M'].values * M_S / M_E)
        Orb_Dist += list(sim.snaps[sim.N_snaps - 1].satellites['a'].values * R_S / au)
        Ecc += list(sim.snaps[sim.N_snaps - 1].satellites['e'].values)
        Inc += list(sim.snaps[sim.N_snaps - 1].satellites['i'].values)
        Types += list(sim.snaps[sim.N_snaps - 1].satellites['Type'].values)

    data = zip(Masses, Orb_Dist, Ecc, Inc, Types)

    data = [item for item in data if item[0] >= m_low_lim and item[1] <= a_up_lim]

    Masses, Orb_Dist, Ecc, Inc, Types = zip(*data)

    number_of_no_accretion = len([item for item in data if np.abs(0.01-item[0])/item[0] < 0.01 and item[-1] == 1])


    print(f'Number of Object: {len(Masses)}')
    print(f'Number of Embryos with no significant accretion: {number_of_no_accretion}, {number_of_no_accretion/len(Masses)}')

    plt.rcParams.update({'figure.autolayout': True})
    plt.style.use('seaborn-paper')
    plt.rcParams.update({'font.size': pop.fontsize})
    plt.rcParams.update({"legend.title_fontsize": pop.legend_fontsize})

    cmap = pop.cmap_standart
    cmin = min(Orb_Dist)
    cmax = max(Orb_Dist)

    norm = colors.LogNorm(cmin, cmax)

    fig, ax = plt.subplots(figsize=pop.figsize)
    ax.scatter(Ecc, np.sin(np.array(Inc)/360 * 2 * np.pi), c=Orb_Dist, cmap=cmap, norm=norm, s=3)
    # ax.scatter(Ecc, np.sin(np.array(Inc)), c=Orb_Dist, cmap=cmap, norm=norm, s=3)
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

    data = [(m, a, wm / m, swm / m) for (m, a, wm, swm) in data if m >= m_low_lim and a <= a_up_lim]

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
    ax.scatter(Orb_Dist, Masses, c=TWMF, cmap=cmap, norm=norm, s=2, alpha=1)
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
    Ecc = []
    System = []

    for key, sim in pop.SIMS.items():
        Masses += list(sim.snaps[sim.N_snaps - 1].satellites['M'].values * M_S / M_E)
        Ecc += list(sim.snaps[sim.N_snaps - 1].satellites['e'].values * M_S / M_E)
        Orb_Dist += list(sim.snaps[sim.N_snaps - 1].satellites['a'].values * R_S / au)
        WM += list(sim.snaps[sim.N_snaps - 1].satellites['WM'].values * M_S / M_E)
        SWM += list(sim.snaps[sim.N_snaps - 1].satellites['SWM'].values * M_S / M_E)
        System += [key for i in sim.snaps[sim.N_snaps - 1].satellites['M'].values]

    WMF = np.array(WM) / np.array(Masses)
    SWMF = np.array(SWM) / np.array(Masses)
    TWMF = WMF + SWMF
    total_number = len(Masses)
    print(f'Total Number of planets: {total_number}')
    data = zip(Masses, Orb_Dist, WMF, SWMF, TWMF, Ecc, System)

    data = [item for item in data if item[0] >= 0.3 and item[0] <= 3]
    mass_lim_number = len(data)
    print(f'Number of planets in mass limit: {mass_lim_number}, {mass_lim_number/total_number}')

    data_copy = data.copy()
    data_wmf = [item for item in data_copy if item[2] > 0.0]
    Masses, Orb_Dist, WMF, SWMF, TWMF, Ecc, System = zip(*data_wmf)
    n_ea_ml_nz_wmf = len(Masses)
    print(f'Number of planets in mass limit with nonzero liquid watermass fraction: {n_ea_ml_nz_wmf}, {n_ea_ml_nz_wmf/mass_lim_number}')

    plt.rcParams.update({'figure.autolayout': True})
    plt.style.use('seaborn-paper')
    plt.rcParams.update({'font.size': pop.fontsize})

    N_bins = 15
    bins = 10 ** np.linspace(np.log10(min(WMF)), np.log10(max(WMF)), N_bins)
    fig, ax = plt.subplots(figsize=pop.figsize)
    # ax.hist(Masses, bins=bins)
    # values, base, _ = plt.hist(Orb_Dist, bins=bins, rwidth=0.95)
    ax.hist(WMF, bins=bins, rwidth=0.95)
    ax.axvline(OE/M_E, color='red', linewidth=1)
    ax.set(xlabel=r'Mass Fraction', ylabel=r'Counts')
    ax.set_xscale('log')
    # ax.set_yscale('log')
    if pop.plot_config == 'presentation':
        ax.set(title=r'Histrogram of Terrestrial Planets Orbital Distances')
    save_name = 'histogram_earth_analogs_wmf'
    if a_up_lim < 30 and m_low_lim > 0:
        save_name += '_lim'
    fig.savefig(path.join(pop.PLOT, save_name + '.png'), transparent=False, dpi=pop.dpi, bbox_inches="tight")
    plt.close(fig)

    data_copy = data.copy()
    data_wmf_lim = [item for item in data_copy if item[2] > 0.0 and item[3] > 0.00075]
    Masses, Orb_Dist, WMF, SWMF, TWMF, Ecc, System = zip(*data_wmf_lim)
    # n_ea_ml_nz_wmf = len(Masses)
    # print(f'Number of planets in mass limit with nonzero liquid watermass fraction: {n_ea_ml_nz_wmf}, {n_ea_ml_nz_wmf/mass_lim_number}')

    plt.rcParams.update({'figure.autolayout': True})
    plt.style.use('seaborn-paper')
    plt.rcParams.update({'font.size': pop.fontsize})

    N_bins = 15
    bins = 10 ** np.linspace(np.log10(min(WMF)), np.log10(max(WMF)), N_bins)
    fig, ax = plt.subplots(figsize=pop.figsize)
    # ax.hist(Masses, bins=bins)
    # values, base, _ = plt.hist(Orb_Dist, bins=bins, rwidth=0.95)
    ax.hist(WMF, bins=bins, rwidth=0.95)
    ax.axvline(OE/M_E, color='red', linewidth=1)
    ax.set(xlabel=r'Mass Fraction', ylabel=r'Counts')
    ax.set_xscale('log')
    # ax.set_yscale('log')
    if pop.plot_config == 'presentation':
        ax.set(title=r'Histrogram of Terrestrial Planets Orbital Distances')
    save_name = 'histogram_earth_analogs_wmf_lim'
    if a_up_lim < 30 and m_low_lim > 0:
        save_name += '_lim'
    fig.savefig(path.join(pop.PLOT, save_name + '.png'), transparent=False, dpi=pop.dpi, bbox_inches="tight")
    plt.close(fig)

    data_copy = data.copy()
    data_swmf_lim = [item for item in data_copy if item[3] > 0.0 and item[2] > 0.00025]
    Masses, Orb_Dist, WMF, SWMF, TWMF, Ecc, System = zip(*data_swmf_lim)
    # n_ea_ml_nz_wmf = len(Masses)
    # print(f'Number of planets in mass limit with nonzero liquid watermass fraction: {n_ea_ml_nz_wmf}, {n_ea_ml_nz_wmf/mass_lim_number}')

    plt.rcParams.update({'figure.autolayout': True})
    plt.style.use('seaborn-paper')
    plt.rcParams.update({'font.size': pop.fontsize})

    N_bins = 15
    bins = 10 ** np.linspace(np.log10(min(WMF)), np.log10(max(WMF)), N_bins)
    fig, ax = plt.subplots(figsize=pop.figsize)
    # ax.hist(Masses, bins=bins)
    # values, base, _ = plt.hist(Orb_Dist, bins=bins, rwidth=0.95)
    ax.hist(WMF, bins=bins, rwidth=0.95)
    ax.axvline(OE/M_E, color='red', linewidth=1)
    ax.set(xlabel=r'Mass Fraction', ylabel=r'Counts')
    ax.set_xscale('log')
    # ax.set_yscale('log')
    if pop.plot_config == 'presentation':
        ax.set(title=r'Histrogram of Terrestrial Planets Orbital Distances')
    save_name = 'histogram_earth_analogs_swmf_lim'
    if a_up_lim < 30 and m_low_lim > 0:
        save_name += '_lim'
    fig.savefig(path.join(pop.PLOT, save_name + '.png'), transparent=False, dpi=pop.dpi, bbox_inches="tight")
    plt.close(fig)

    data_copy = data.copy()
    data_swmf = [item for item in data_copy if item[3] > 0.0]
    Masses, Orb_Dist, WMF, SWMF, TWMF, Ecc, System = zip(*data_swmf)
    n_ea_ml_nz_swmf = len(Masses)
    print(f'Number of planets in mass limit with nonzero hydrated solids watermass fraction: {n_ea_ml_nz_swmf}, {n_ea_ml_nz_swmf/mass_lim_number}')

    plt.rcParams.update({'figure.autolayout': True})
    plt.style.use('seaborn-paper')
    plt.rcParams.update({'font.size': pop.fontsize})

    N_bins = 15
    bins = 10 ** np.linspace(np.log10(min(SWMF)), np.log10(max(SWMF)), N_bins)
    fig, ax = plt.subplots(figsize=pop.figsize)
    # ax.hist(Masses, bins=bins)
    # values, base, _ = plt.hist(Orb_Dist, bins=bins, rwidth=0.95)
    ax.hist(WMF, bins=bins, rwidth=0.95)
    ax.axvline(3 * OE/M_E, color='red', linewidth=1)
    ax.set(xlabel=r'Mass Fraction', ylabel=r'Counts')
    ax.set_xscale('log')
    # ax.set_yscale('log')
    if pop.plot_config == 'presentation':
        ax.set(title=r'Histrogram of Terrestrial Planets Orbital Distances')
    save_name = 'histogram_earth_analogs_swmf'
    if a_up_lim < 30 and m_low_lim > 0:
        save_name += '_lim'
    fig.savefig(path.join(pop.PLOT, save_name + '.png'), transparent=False, dpi=pop.dpi, bbox_inches="tight")
    plt.close(fig)

    data_copy = data.copy()
    data_twmf = [item for item in data_copy if item[2] > 0.0 and item[3] > 0.0]
    Masses, Orb_Dist, WMF, SWMF, TWMF, Ecc, System = zip(*data_twmf)
    n_ea_ml_nz_twmf = len(Masses)
    print(f'Number of planets in mass limit with nonzero wmf and swmf: {n_ea_ml_nz_twmf}, {n_ea_ml_nz_twmf/mass_lim_number}')

    plt.rcParams.update({'figure.autolayout': True})
    plt.style.use('seaborn-paper')
    plt.rcParams.update({'font.size': pop.fontsize})
    ratios = np.array(WMF)/np.array(SWMF)
    N_bins = 15
    bins = 10 ** np.linspace(np.log10(min(ratios)), np.log10(max(ratios)), N_bins)
    fig, ax = plt.subplots(figsize=pop.figsize)
    # ax.hist(Masses, bins=bins)
    # values, base, _ = plt.hist(Orb_Dist, bins=bins, rwidth=0.95)
    ax.hist(ratios, bins=bins, rwidth=0.95)
    ax.axvline(1/3, color='red', linewidth=1)
    ax.set(xlabel=r'Ratio', ylabel=r'Counts')
    ax.set_xscale('log')
    ax.set_yscale('log')
    if pop.plot_config == 'presentation':
        ax.set(title=r'Histrogram of Terrestrial Planets Orbital Distances')
    save_name = 'histogram_earth_analogs_twmf'
    if a_up_lim < 30 and m_low_lim > 0:
        save_name += '_lim'
    fig.savefig(path.join(pop.PLOT, save_name + '.png'), transparent=False, dpi=pop.dpi, bbox_inches="tight")
    plt.close(fig)

    data = [item for item in data if item[4] > 0]
    non_zero_wm = data.copy()
    non_zero_wmf_number = len(data)
    print(f'Number of planets in mass limit with non zero TWMF: {non_zero_wmf_number}, {non_zero_wmf_number/mass_lim_number} ({non_zero_wmf_number/total_number})')



    earth_analogs = [item for item in data if item[0] >= 0.101 and item[1] <= a_up_lim and item[2] > 0.001]
    #print(earth_analogs)

    Masses, Orb_Dist, WMF, SWMF, TWMF, Ecc, System = zip(*data)



    plt.rcParams.update({'figure.autolayout': True})
    plt.style.use('seaborn-paper')
    plt.rcParams.update({'font.size': pop.fontsize})
    plt.rcParams.update({"legend.title_fontsize": pop.legend_fontsize})

    cmap = pop.cmap_standart
    TWMF = TWMF
    cmin = min(TWMF)
    cmax = max(TWMF)

    norm = colors.LogNorm(cmin, cmax)

    fig, ax = plt.subplots(figsize=pop.figsize)
    ax.scatter(Orb_Dist, Masses, c=TWMF, cmap=cmap, norm=norm, s=7, alpha=1)
    # ax.scatter(obs, ms, c=twmf, cmap=cmap, norm=norm, s=10)
    ax.axvline(1, color='black', linewidth=0.7, linestyle='--')
    ax.axhline(1, color='black', linewidth=0.7, linestyle='--')
    x_labels = ax.get_xticklabels()
    plt.setp(x_labels, horizontalalignment='center', )

    ax.set(xlabel='Orbital Distance [$\mathrm{au}$]', ylabel=r'Mass [$\mathrm{M_{\oplus}}$]')
    fig.colorbar(cm.ScalarMappable(cmap=cmap, norm=norm), orientation='vertical', label=r'Total WMF', ax=ax)
    ax.set_xscale('log')
    ax.set_yscale('log')
    if pop.plot_config == 'presentation':
        ax.set(title=r'Total WMF Radial Distribution')
    save_name = 'scatter_radial_twmf'
    if a_up_lim < 30 and m_low_lim > 0:
        save_name += '_lim'
    fig.savefig(path.join(pop.PLOT, save_name + '.png'), transparent=False, dpi=pop.dpi, bbox_inches="tight")
    plt.close(fig)

    def scatter_pie(earth_analogs):

        plt.rcParams.update({'figure.autolayout': True})
        plt.style.use('seaborn-paper')
        plt.rcParams.update({'font.size': pop.fontsize})
        plt.rcParams.update({"legend.title_fontsize": pop.legend_fontsize})

        fig, ax = plt.subplots(figsize=pop.figsize)

        colors = ['red', 'blue']
        labels = ['Hydrated Silica', 'Water/Ice']

        red_patch = patches.Patch(color='red', label='Hydrated Silica')
        blue_patch = patches.Patch(color='blue', label='Water/Ice')
        handles = [red_patch, blue_patch]

        Masses, Orb_Dist, WMF, SWMF, TWMF, Ecc, System = zip(*earth_analogs)

        mean_mass = np.min(Masses)
        mass_scaling = mean_mass / 90000
        mass_scaling = 0.000000001

        def pie_1d(r1, r2):
            # calculate the points of the first pie marker
            # these are just the origin (0, 0) + some (cos, sin) points on a circle
            x1 = np.cos(2 * np.pi * np.linspace(0, r1))
            y1 = np.sin(2 * np.pi * np.linspace(0, r1))
            xy1 = np.row_stack([[0, 0], np.column_stack([x1, y1])])
            s1 = np.abs(xy1).max()

            x2 = np.cos(2 * np.pi * np.linspace(r1, 1))
            y2 = np.sin(2 * np.pi * np.linspace(r1, 1))
            xy2 = np.row_stack([[0, 0], np.column_stack([x2, y2])])
            s2 = np.abs(xy2).max()

            # x3 = np.cos(2 * np.pi * np.linspace(r2, 1))
            # y3 = np.sin(2 * np.pi * np.linspace(r2, 1))
            # xy3 = np.row_stack([[0, 0], np.column_stack([x3, y3])])
            # s3 = np.abs(xy3).max()

            return xy1, s1, xy2, s2#, xy3, s3

        # cale the masses to the marker sizes
        # def NormalizeData(m):
        #     return (np.log10(m) - np.log10(np.min(TWMF))) / (np.log10(np.max(TWMF)) - np.log10(np.min(TWMF)))

        def NormalizeData(m):
            return (np.log10(m) - np.log10(np.min(Masses))) / (np.log10(np.max(Masses)) - np.log10(np.min(Masses)))

        # def NormalizeData(m):
        #     return (m - (np.min(TWMF))) / ((np.max(TWMF)) - (np.min(TWMF)))

        # def NormalizeData(m):
        #     return (m - (np.min(Masses))) / ((np.max(Masses)) - (np.min(Masses)))

        earth_point = (1,1,0.00025,0.00075,0.001,0,0)
        def plot_one(row,earth=False):
            WMF_ratio = row[2]/row[4]
            SWMF_Ratio = 1
            #xy1, s1, xy2, s2, xy3, s3 = pie_1d(WMF_ratio, SWMF_ratio)
            xy1, s1, xy2, s2 = pie_1d(WMF_ratio, 1)

            scale = NormalizeData(row[0]) * 50
            if earth == True:
                ax.scatter(row[1], row[4], s=s2 * scale * 2, facecolor='green')
            ax.scatter(row[1], row[4], marker=xy1, s=s1 * scale , facecolor='blue')
            ax.scatter(row[1], row[4], marker=xy2, s=s2 * scale, facecolor='red')

            #ax.scatter(row[1], row[6], marker=xy3, s=s3 * scale , facecolor='red')

        for index, row in enumerate(earth_analogs):
            plot_one(row)
        plot_one(earth_point,True)
        #ax.set_ylim(-1 * min(self.satellites['e']), 1.1 * max(self.satellites['e']))
        ax.set_xlabel(r'Orbital Distance [$\mathrm{au}$]')
        ax.set_ylabel('Total Water Mass Fractions')
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.legend(handles=handles, title='Components')
        fig.savefig(path.join(pop.PLOT, 'scatter_ratios.png'), transparent=False, dpi=pop.dpi, bbox_inches="tight")
        plt.close(fig)

    #filter twmf close to earth
    data = [item for item in data if item[2] >= 0.00025 and item[3] >= 0.00075]
    systems_id = [sys[-1] for sys in data]
    print(f'Number of systems with Earth Candidate {len(np.unique(systems_id))}, {len(np.unique(systems_id))/pop.NSIMS} ')
    scatter_pie(non_zero_wm)

    earth_analogs2 = data.copy()

    wmf_sim_number = len(data)
    #print(data)
    print(f'Number of planets in mass limit and WMF above 0.00025 and SWMF above 0.00075: {len(data)}, {wmf_sim_number/mass_lim_number}  ({wmf_sim_number/total_number})')
    # for earth in data:
    #     print(f'System: {earth[-1]}')
    #     print(f'Mass: {earth[0]}')
    #     print(f'Orb Dist: {earth[1]}')
    #     print(f'WMF: {earth[2]}')
    #     print(f'SWMF: {earth[2]}')
    #     print(f'TWMF: {earth[2]}')
    #     print(f'Exponent: {pop.SIMS[earth[-1]].Sigma_Exponent}')
    #     print(f'Disk Mass: {pop.SIMS[earth[-1]].Total_Mass * M_S / M_J}')
    #     print(" ")

    ms, obs, wmf, swmf, twmf, ecc, system = zip(*earth_analogs2)

    plt.rcParams.update({'figure.autolayout': True})
    plt.style.use('seaborn-paper')
    plt.rcParams.update({'font.size': pop.fontsize})
    plt.rcParams.update({"legend.title_fontsize": pop.legend_fontsize})

    cmap = pop.cmap_standart
    cmin = min(twmf)
    cmax = max(twmf)

    norm = colors.Normalize(cmin, cmax)

    fig, ax = plt.subplots(figsize=pop.figsize)
    ax.scatter(wmf, swmf, c=twmf, cmap=cmap, norm=norm, s=2, alpha=1)
    x_labels = ax.get_xticklabels()
    plt.setp(x_labels, horizontalalignment='center')
    ax.set(xlabel=r'Water Mass Fraction]', ylabel=r'Solids Water Mass Fraction')
    ax.axvline(0.00025, color='black', linewidth=0.7, linestyle='--')
    ax.axhline(0.00075, color='black', linewidth=0.7, linestyle='--')
    fig.colorbar(cm.ScalarMappable(cmap=cmap, norm=norm), orientation='vertical',
                 label=r'Total WMF',
                 ax=ax)
    ax.set_xscale('log')
    ax.set_yscale('log')
    if pop.plot_config == 'presentation':
        ax.set(title=r'Eccentricity and Inclination')
    save_name = 'scatter_wmf_swmf'
    if a_up_lim < 30 and m_low_lim > 0:
        save_name += '_lim'
    fig.savefig(path.join(pop.PLOT, save_name + '.png'), transparent=False, dpi=pop.dpi, bbox_inches="tight")
    plt.close(fig)


    #filter roughly earth mass already roughly in the right positions
    data = [item for item in data if item[1] <= 2]
    earth_like_number = len(data)
    print(f'Number of planets in mass limit and WMF above 0.00025 and SWMF above 0.00075 at correct positions: {earth_like_number}, {earth_like_number/wmf_sim_number}  ({earth_like_number/total_number})')

    SE = []
    TM = []
    RE = []
    for id in np.unique([sys[-1] for sys in earth_analogs2]):
        SE.append(pop.SIMS[id].Sigma_Exponent)
        TM.append(pop.SIMS[id].Total_Mass)
        RE.append(pop.SIMS[id].Sigma_Norm / (R_S / au)**pop.SIMS[id].Sigma_Exponent / denstos * pow(au/R_S, pop.SIMS[id].Sigma_Exponent) / denstos)

    SE = np.array(SE)
    TM = np.array(TM)

    cmap = pop.cmap_standart
    cmin = min(RE)
    cmax = max(RE)

    norm = colors.LogNorm(cmin, cmax)
    fig, ax = plt.subplots(figsize=pop.figsize)
    ax.scatter(SE, TM, c=RE, cmap=cmap, norm=norm, s=12)
    x_labels = ax.get_xticklabels()
    plt.setp(x_labels, horizontalalignment='center')
    ax.set(xlabel='Surface Density Power Law Exponent', ylabel=r'Total Disk Mass [$M_{\odot}$]', xticks=SE)
    ax2 = ax.twinx()
    mn, mx = ax.get_ylim()
    ax2.set_ylim(M_S / M_J * mn, M_S / M_J * mx)
    ax2.set_ylabel('Total Disk Mass [$M_{J}$]')
    fig.colorbar(cm.ScalarMappable(cmap=cmap, norm=norm), orientation='vertical',
                 label=r'Reference Value at $1 \mathrm{au}$ [$\mathrm{g}\mathrm{cm}^{-2}$]', ax=ax2, pad=0.12)
    # ax.set_yscale('log')
    if pop.plot_config == 'presentation':
        ax.set(title=r'Synthesis Parameters')

    fig.savefig(path.join(pop.PLOT, 'scatter_parameters_earth_analogs.png'), transparent=False, dpi=pop.dpi,
                bbox_inches="tight")
    plt.close(fig)






