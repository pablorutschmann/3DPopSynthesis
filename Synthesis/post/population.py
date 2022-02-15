import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib import colors
import os.path as path
from os import makedirs
from ..units import *
from .simulation import simulation
from .metric import *

figsize = (8, 6)
fontsize = 60
dpi = 600
cmap_standart = 'viridis_r'


class population:
    def __init__(self, PATH, N_sims):
        self.MAIN = PATH
        self.NAME = path.basename(self.MAIN)
        self.PLOT = path.join(self.MAIN, 'plots')
        self.NSIMS = N_sims
        self.SIMS = {}
        if not path.exists(self.PLOT):
            makedirs(self.PLOT)

        for i in range(1, self.NSIMS + 1):
            self.SIMS[i] = simulation(self.system_path(i))
        print(self.SIMS)

    def scatter_parameters(self):
        TotalMasses = []
        SigmaCoeffs = []
        Reference = []

        for sim in self.SIMS.values():
            TotalMasses.append(sim.Total_Mass)
            SigmaCoeffs.append(sim.Sigma_Exponent)
            Reference.append(sim.Sigma_Norm * pow(R_S / au, sim.Sigma_Exponent) * (M_S / R_S2))

        plt.rcParams.update({'figure.autolayout': True})
        plt.style.use('seaborn-paper')
        plt.rcParams.update({'font.size': fontsize})

        cmap = cmap_standart
        cmin = min(Reference)
        cmax = max(Reference)

        norm = colors.LogNorm(cmin, cmax)

        fig, ax = plt.subplots(figsize=figsize)
        ax.scatter(SigmaCoeffs, TotalMasses, c=Reference, cmap=cmap, norm=norm, s=12)
        x_labels = ax.get_xticklabels()
        plt.setp(x_labels, horizontalalignment='center')
        ax.set(xlabel='Surface Density Power Law Exponent', ylabel=r'Total Disk Mass [$M_{\odot}$]', xticks=SigmaCoeffs)
        ax2 = ax.twinx()
        mn, mx = ax.get_ylim()
        ax2.set_ylim(M_S / M_J * mn, M_S / M_J * mx)
        ax2.set_ylabel('Total Disk Mass [$M_{J}$]')
        fig.colorbar(cm.ScalarMappable(cmap=cmap, norm=norm), orientation='vertical',
                     label=r'Reference Value at $1 \mathrm{au}$ [$\mathrm{g}\mathrm{cm}^{-2}$]', ax=ax2, pad=0.09)
        # ax.set_yscale('log')
        fig.savefig(path.join(self.PLOT, 'scatter_parameters.png'), transparent=False, dpi=dpi, bbox_inches="tight")
        plt.close(fig)

    def scatter_ecc_inc(self, m_low_lim=0, a_up_lim=30):
        Masses = []
        Orb_Dist = []
        Ecc = []
        Inc = []

        for sim in self.SIMS.values():
            Masses += list(sim.snaps[sim.N_snaps - 1].satellites['M'].values * M_S / M_E)
            Orb_Dist += list(sim.snaps[sim.N_snaps - 1].satellites['a'].values * R_S / au)
            Ecc += list(sim.snaps[sim.N_snaps - 1].satellites['e'].values)
            Inc += list(sim.snaps[sim.N_snaps - 1].satellites['i'].values)

        data = zip(Masses, Orb_Dist, Ecc, Inc)

        data = [item for item in data if item[0] >= m_low_lim and item[1] <= a_up_lim]

        Masses, Orb_Dist, Ecc, Inc = zip(*data)

        plt.rcParams.update({'figure.autolayout': True})
        plt.style.use('seaborn-paper')
        plt.rcParams.update({'font.size': fontsize})

        cmap = cmap_standart
        cmin = min(Masses)
        cmax = max(Masses)

        norm = colors.LogNorm(cmin, cmax)

        fig, ax = plt.subplots(figsize=figsize)
        ax.scatter(Ecc, np.sin(Inc), c=Orb_Dist, cmap=cmap, norm=norm, s=12)
        x_labels = ax.get_xticklabels()
        plt.setp(x_labels, horizontalalignment='center')
        ax.set(xlabel='Eccentricity', ylabel=r'$\sin(\mathrm{inclination})$')
        fig.colorbar(cm.ScalarMappable(cmap=cmap, norm=norm), orientation='vertical', label=r'Mass [$\mathrm{M_E}$]',
                     ax=ax)
        ax.set_xscale('log')
        ax.set_yscale('log')
        fig.savefig(path.join(self.PLOT, 'scatter_ecc_inc.png'), transparent=False, dpi=dpi, bbox_inches="tight")
        plt.close(fig)

    def scatter_AMD_RMC(self, m_low_lim, a_up_lim):
        RMC_sol = 89.9
        AMD_sol = 0.0018

        AMDS = []
        RMCS = []
        NS = []

        for sys in self.SIMS.values():
            AMD, N = sys.get_AMD(m_low_lim, a_up_lim)
            AMDS.append(AMD)
            RMC, N = sys.get_RMC(m_low_lim, a_up_lim)
            RMCS.append(RMC)
            NS.append(N)

        cmap = cmap_standart
        cmin = min(NS)
        cmax = max(NS)

        norm = colors.Normalize(cmin, cmax)

        plt.rcParams.update({'figure.autolayout': True})
        plt.style.use('seaborn-paper')
        plt.rcParams.update({'font.size': fontsize})

        fig, ax = plt.subplots(figsize=figsize)
        ax.scatter(RMCS, AMDS, c=NS, cmap=cmap, norm=norm, s=12)
        ax.scatter(RMC_sol, AMD_sol, c='red')
        x_labels = ax.get_xticklabels()
        plt.setp(x_labels, horizontalalignment='center')
        ax.set(xlabel='RMC', ylabel=r'AMC')
        fig.colorbar(cm.ScalarMappable(cmap=cmap, norm=norm), orientation='vertical', label=r'Number of Satellites',
                     ax=ax)
        ax.set_xscale('log')
        ax.set_yscale('log')
        fig.savefig(path.join(self.PLOT, 'scatter_AMD_RMC.png'), transparent=False, dpi=dpi, bbox_inches="tight")
        plt.close(fig)

    def histogram_mass(self, m_low_lim=0, a_up_lim=30):
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
        plt.rcParams.update({'font.size': fontsize})

        N_bins = 15
        bins = 10 ** np.linspace(np.log10(min(Masses)), np.log10(max(Masses)), N_bins)
        fig, ax = plt.subplots(figsize=figsize)
        # ax.hist(Masses, bins=bins)
        values, base, _ = plt.hist(Masses, bins=bins, rwidth=0.95)
        ax_bis = ax.twinx()
        # values, base = np.histogram(Masses, bins = 10 ** np.linspace(np.log10(min(Masses)), np.log10(max(Masses)), N_bins))

        values = np.append(0, values)
        ax_bis.plot(base, np.cumsum(values) / np.cumsum(values)[-1], color='black', linestyle='dashed', markersize=0.1)
        ax.set(xlabel=r'M [$M_E$]', ylabel=r'Counts')
        ax_bis.set(ylabel='Cumulative Distribution')
        ax.set_xscale('log')
        ax_bis.set_xscale('log')
        fig.savefig(path.join(self.PLOT, 'histogram_mass.png'), transparent=False, dpi=dpi, bbox_inches="tight")
        plt.close(fig)

    def histogram_weighted_mass(self, m_low_lim=0, a_up_lim=30):
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
        plt.rcParams.update({'font.size': fontsize})

        N_bins = 15
        bins = 10 ** np.linspace(np.log10(min(Masses)), np.log10(max(Masses)), N_bins)
        fig, ax = plt.subplots(figsize=figsize)
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
        ax.set(xlabel=r'M [$M_E$]', ylabel=r'Weighed Counts')
        ax_bis.set(ylabel='Cumulative Distribution')
        ax.set_xscale('log')
        ax_bis.set_xscale('log')
        fig.savefig(path.join(self.PLOT, 'histogram_weighted_mass.png'), transparent=False, dpi=dpi,
                    bbox_inches="tight")
        plt.close(fig)

    def histogram_a(self, m_low_lim=0, a_up_lim=30):
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
        plt.rcParams.update({'font.size': fontsize})

        N_bins = 15
        bins = 10 ** np.linspace(np.log10(min(Orb_Dist)), np.log10(max(Orb_Dist)), N_bins)
        fig, ax = plt.subplots(figsize=figsize)
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
        ax.set(xlabel=r'a [$\mathrm{au}$]', ylabel=r'Counts')
        ax_bis.set(ylabel='Cumulative Distribution')
        ax.set_xscale('log')
        ax_bis.set_xscale('log')
        fig.savefig(path.join(self.PLOT, 'histogram_a.png'), transparent=False, dpi=dpi, bbox_inches="tight")
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
        plt.rcParams.update({'font.size': fontsize})

        N_bins = 15
        bins = 10 ** np.linspace(np.log10(min(Orb_Dist)), np.log10(max(Orb_Dist)), N_bins)
        fig, ax = plt.subplots(figsize=figsize)
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
        ax.set(xlabel=r'a [$\mathrm{au}$]', ylabel=r'Counts')
        ax_bis.set(ylabel='Cumulative Distribution')
        ax.set_xscale('log')
        ax_bis.set_xscale('log')
        fig.savefig(path.join(self.PLOT, 'histogram_weighted_a.png'), transparent=False, dpi=dpi, bbox_inches="tight")
        plt.close(fig)

    def histogram_totalmass(self, m_low_lim=0):
        TotalMasses = []

        for sim in self.SIMS.values():
            TotalMasses.append(np.sum(
                [item for item in list(sim.snaps[sim.N_snaps - 1].satellites['M'].values * M_S / M_E) if
                 item >= m_low_lim]))

        plt.rcParams.update({'figure.autolayout': True})
        plt.style.use('seaborn-paper')
        plt.rcParams.update({'font.size': fontsize})

        N_bins = 15
        bins = 10 ** np.linspace(np.log10(min(TotalMasses)), np.log10(max(TotalMasses)), N_bins)
        fig, ax = plt.subplots(figsize=figsize)
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
        fig.savefig(path.join(self.PLOT, 'histogram_totalmass.png'), transparent=False, dpi=dpi, bbox_inches="tight")
        plt.close(fig)

    def distances(self, m_low_lim, a_up_lim):
        systems = []
        total_masses = []
        sigma_exponents = []
        m_low_lim = m_low_lim
        a_up_lim = a_up_lim

        for sim in self.SIMS.values():
            total_masses.append(sim.Total_Mass)
            sigma_exponents.append(sim.Sigma_Exponent)

            zipped = zip(sim.snaps[sim.N_snaps - 1].satellites['M'].values,
                         sim.snaps[sim.N_snaps - 1].satellites['a'].values)

            # Remove under threshold Masses
            filtered = [(m * M_S / M_E, a * R_S / au) for (m, a) in zipped if
                        m >= m_low_lim * M_E / M_S and a <= a_up_lim * au / R_S]
            systems.append(filtered)

        func = lambda x: distance(x, terrestrial)
        distances = list(map(func, systems))

        plt.rcParams.update({'figure.autolayout': True})
        plt.style.use('seaborn-paper')
        plt.rcParams.update({'font.size': fontsize})

        fig, ax = plt.subplots(ncols=1)
        fig.set_size_inches(15.5, 10.5)
        ax.scatter(sigma_exponents, distances)

        ax.set_xlabel('Number of initial Planetesimals')

        ax.set_ylabel('Distance')
        ax.legend()
        fig.suptitle('Distance Metric to Solar System')
        fig.savefig(path.join(self.PLOT, 'distances_exponent.png'))
        plt.close(fig)

        fig, ax = plt.subplots(ncols=1)
        fig.set_size_inches(15.5, 10.5)
        ax.scatter(total_masses, distances)

        ax.set_xlabel('Number of initial Planetesimals')

        ax.set_ylabel('Distance')
        ax.legend()
        fig.suptitle('Distance Metric to Solar System')
        fig.savefig(path.join(self.PLOT, 'distances_mass.png'))
        plt.close(fig)

        cmap = 'viridis_r'
        cmin = min(distances)
        cmax = max(distances)

        norm = colors.Normalize(cmin, cmax)

        fig, ax = plt.subplots(figsize=figsize)
        ax.scatter(sigma_exponents, total_masses, c=distances, cmap=cmap, norm=norm, s=12)
        x_labels = ax.get_xticklabels(list(set(sigma_exponents)))
        plt.setp(x_labels, horizontalalignment='center')
        ax.set(xlabel='Sigma Exponent', ylabel=r'Total Mass', xticks=sigma_exponents)
        fig.colorbar(cm.ScalarMappable(cmap=cmap, norm=norm), orientation='vertical', label=r'Distance',
                     ax=ax)
        ax.set_xscale('log')
        ax.set_yscale('log')
        fig.savefig(path.join(self.PLOT, 'scatter_distances.png'), transparent=False, dpi=dpi, bbox_inches="tight")
        plt.close(fig)

    def line_n_planets(self, a_up_lim, thresholds=None):
        if thresholds is None:
            thresholds = [0.05, 0.1, 0.15]

        thresholds = sorted(thresholds)
        numbers = {key: [] for key in thresholds}
        occurences = {key: [] for key in thresholds}
        systems = []

        for sim in self.SIMS.values():
            zipped = list(zip(sim.snaps[sim.N_snaps - 1].satellites['M'].values * M_S / M_E,
                              sim.snaps[sim.N_snaps - 1].satellites['a'].values * R_S / au))
            zipped = [(m, a) for (m, a) in zipped if a <= a_up_lim]
            systems.append(zipped)
        number_of_systems = len(systems)
        for keys, lists in numbers.items():
            for i, item in enumerate(systems):
                filtered = [(m, a) for (m, a) in item if
                            m >= keys]
                leng = len(filtered)
                lists += [leng]
                systems[i] = filtered

        max_number = max(numbers[thresholds[0]])
        for key, value in numbers.items():
            for i in range(max_number + 1):
                occurences[key] += [value.count(i)]

        plt.rcParams.update({'figure.autolayout': True})
        plt.style.use('seaborn-paper')
        plt.rcParams.update({'font.size': fontsize})

        fig, ax = plt.subplots(figsize=figsize)
        for key, value in occurences.items():
            ax.plot(range(max_number + 1), np.array(value) / number_of_systems, linestyle='dashed', linewidth=0.5, markersize=4,
                    marker='o', label=f'Threshold {key}' + r' $\mathrm{M_{\oplus}}$')
        ax.set(xlabel='Number of Planets', ylabel=r'Normalized Occurence', xticks=range(max_number + 1))
        ax.legend(loc='best')
        fig.savefig(path.join(self.PLOT, 'line_planet_number_occurences.png'), transparent=False, dpi=dpi,
                    bbox_inches="tight")
        plt.close(fig)

    def a_wm_distribution(self):
        Masses = []
        A = []
        WM = []
        SWM = []
        for sim in self.SIMS.values():
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
        fig.savefig(path.join(self.PLOT, 'a_mass_distribution.png'))
        plt.close(fig)

    def final_wmf_radial_distribution(self, cumulative=False, density=False, thresholds=[10e-12]):
        Masses = []
        SWM = []
        WM = []

        for sim in self.SIMS.values():
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

        fig, ax = plt.subplots(ncols=2)
        fig.set_size_inches(15.5, 10.5)

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

        fig.suptitle('Final Water Mass Fraction Distribution')
        fig.savefig(path.join(self.PLOT, 'final_wm_distribution' + str(cumulative) + str(density) + '.png'))
        plt.close(fig)

    def system_path(self, i):
        sys = "system_" + str(i)
        return path.join(self.MAIN, sys)

    def system_output_path(self, i):
        return path.join(self.system_path(i), 'outputs')

    def system_output_path(self, i):
        return path.join(self.system_path(i), 'inputs')

    def __str__(self):
        return f"""
      Synthesis Run:
          Run Name: {self.NAME}
          Main Path: {self.MAIN}
          Number Of Simulations: {self.NSIMS}
      """


if __name__ == "__main__":
    TEST = population("test", 9)
