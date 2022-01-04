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


class population:
    def __init__(self, PATH, N_sims):
        self.MAIN = PATH
        self.NAME = path.basename(self.MAIN)
        self.PLOT = path.join(self.MAIN, 'plots')
        self.NSIMS = N_sims
        self.SIMS = {}
        if not path.exists(self.PLOT):
            makedirs(self.PLOT)

        for i in range(1, self.NSIMS):
            self.SIMS[i] = simulation(self.system_path(i))
        print(self.SIMS)

    def total_mass_distribution(self, cumulative=False, density=False, thresholds=[10e-12]):
        Masses = []

        for sim in self.SIMS.values():
            Masses.append(np.sum(sim.snaps[sim.N_snaps - 1].satellites['M'].values))

        data = []
        if density == False:
            thresholds = [M_S / M_E * thresh for thresh in thresholds]
            for thresh in thresholds:
                data.append([x for x in Masses if x > thresh])

        N_bins = 30
        bins = np.logspace(np.log10(min(Masses)), np.log10(max(Masses)), N_bins)

        print(len(Masses))
        fig, ax = plt.subplots(ncols=1)
        fig.set_size_inches(15.5, 10.5)
        if density == False:
            for dat, thresh in zip(data, thresholds):
                ax.hist(dat, bins=bins, cumulative=cumulative, label=str(thresh))
        else:
            ax.hist(Masses, bins=bins, density=True, cumulative=cumulative, stacked=True)
        ax.vlines(np.mean(Masses), 0, 700, color='r', linestyle='dashed', linewidth=2, label='Mean Mass')
        ax.vlines(M_M / M_E, 0, 700, color='r', linestyle='dashed', linewidth=2, label='Mars Mass')
        ax.set_xlabel('Mass in Earth Masses')
        if density == False:
            ax.set_ylabel('Counts')
        else:
            ax.set_ylabel('Probability')
        ax.set_xscale('log')
        ax.legend()
        fig.suptitle('Final Mass Distribution')
        fig.savefig(path.join(self.PLOT, 'final_mass_distribution' + str(cumulative) + str(density) + '.png'))
        plt.close(fig)

    def a_mass_distribution(self):
        Masses = []
        A = []
        WM = []
        SWM = []
        for sim in self.SIMS.values():
            Masses.extend(sim.snaps[sim.N_snaps - 1].satellites['M'].values * M_S / M_E)
            A.extend(sim.snaps[sim.N_snaps - 1].satellites['a'].values * R_S / au)
            WM.extend(sim.snaps[sim.N_snaps - 1].satellites['WM'].values * M_S / M_E)
            SWM.extend(sim.snaps[sim.N_snaps - 1].satellites['SWM'].values * M_S / M_E)

        Masses = np.asarray(Masses)
        WM = np.asarray(WM)
        SWM = np.asarray(SWM)
        total_wmf = (WM + SWM) / Masses

        cmap = 'cividis_r'
        cmin = min(total_wmf)
        cmax = max(total_wmf)

        norm = colors.Normalize(cmin, cmax)

        fig, ax = plt.subplots(ncols=1)
        fig.set_size_inches(15.5, 10.5)
        ax.scatter(A, Masses, c=total_wmf, cmap=cmap, norm=norm)
        fig.colorbar(cm.ScalarMappable(cmap=cmap, norm=norm), orientation='vertical', label="Water Mass Fraction",
                     ax=ax)
        ax.set_xlabel('Distane from Star in au')
        ax.set_ylabel('Mass in Earth Masses')
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.legend()
        fig.suptitle('Distance Mass Distribution')
        fig.savefig(path.join(self.PLOT, 'a_mass_distribution.png'))
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

    def final_mass_distribution(self, cumulative=False, density=False, thresholds=[10e-12]):
        Masses = []

        for sim in self.SIMS.values():
            Masses.extend(sim.snaps[sim.N_snaps - 1].satellites['M'].values * M_S / M_E)
        data = []
        if density == False:
            thresholds = [M_S / M_E * thresh for thresh in thresholds]
            for thresh in thresholds:
                data.append([x for x in Masses if x > thresh])

        N_bins = 30
        bins = np.logspace(np.log10(min(Masses)), np.log10(max(Masses)), N_bins)

        print(len(Masses))
        fig, ax = plt.subplots(ncols=1)
        fig.set_size_inches(15.5, 10.5)
        if density == False:
            for dat, thresh in zip(data, thresholds):
                ax.hist(dat, bins=bins, cumulative=cumulative, label=str(thresh))
        else:
            ax.hist(Masses, bins=bins, density=True, cumulative=cumulative, stacked=True)
        ax.vlines(np.mean(Masses), 0, 700, color='r', linestyle='dashed', linewidth=2, label='Mean Mass')
        ax.vlines(M_M / M_E, 0, 700, color='r', linestyle='dashed', linewidth=2, label='Mars Mass')
        ax.set_xlabel('Mass in Earth Masses')
        if density == False:
            ax.set_ylabel('Counts')
        else:
            ax.set_ylabel('Probability')
        ax.set_xscale('log')
        ax.legend()
        fig.suptitle('Final Mass Distribution')
        fig.savefig(path.join(self.PLOT, 'final_mass_distribution' + str(cumulative) + str(density) + '.png'))
        plt.close(fig)

    def final_radius_distribution(self, cumulative=False, density=False, thresholds=[10e-12]):
        Masses = []
        A = []

        for sim in self.SIMS.values():
            Masses.extend(sim.snaps[sim.N_snaps - 1].satellites['M'].values * M_S / M_E)
            A.extend(sim.snaps[sim.N_snaps - 1].satellites['a'].values * R_S / au)
        data = []
        if not density:
            thresholds = [thresh * M_S / M_E for thresh in thresholds]
            for thresh in thresholds:
                data.append([(m, a) for (m, a) in zip(Masses, A) if m > thresh])

        N_bins = 30
        bins = np.logspace(np.log10(min(A)), np.log10(max(A)), N_bins)
        print(len(Masses))
        fig, ax = plt.subplots(ncols=1)
        fig.set_size_inches(15.5, 10.5)
        if not density:
            for dat, thresh in zip(data, thresholds):
                ax.hist([x[1] for x in dat], bins=bins, density=density, cumulative=cumulative,
                        label=str(thresh))
        else:
            ax.hist(A, bins=bins, density=density, cumulative=cumulative, stacked=True)

        ax.set_xlabel('Distance from Star in AU')
        ax.set_xscale('log')

        if not density:
            ax.set_ylabel('Counts')
        else:
            ax.set_ylabel('Probability')
        ax.legend()
        fig.suptitle('Final Radius Distribution')
        fig.savefig(path.join(self.PLOT, 'final_radius_distribution' + str(cumulative) + str(density) + '.png'))
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

    def distances(self, threshold=(0.05 * M_E / M_S)):
        systems = []
        parameters = []
        thresh = threshold * M_S / M_E

        for sim in self.SIMS.values():
            zipped = zip(sim.snaps[sim.N_snaps - 1].satellites['M'].values * M_S / M_E,
                         sim.snaps[sim.N_snaps - 1].satellites['a'].values * R_S / au)

            # Remove under threshold Masses
            filtered = [(m, a) for (m, a) in zipped if m > thresh]
            systems.append(filtered)

            parameters.append(sim.N_planetesimals)

        fun = lambda x: distance(x, terrestrial)
        distances = list(map(fun, systems))
        print(distances)

        fig, ax = plt.subplots(ncols=1)
        fig.set_size_inches(15.5, 10.5)
        ax.scatter(parameters, distances)

        ax.set_xlabel('Number of initial Planetesimals')

        ax.set_ylabel('Distance')
        ax.legend()
        fig.suptitle('Distance Metric to Solar System')
        fig.savefig(path.join(self.PLOT, 'distances.png'))
        plt.close(fig)

    def ecc_inc_distribution(self):
        Masses = []
        Inc = []
        Ecc = []

        for sim in self.SIMS.values():
            Masses.extend(sim.snaps[sim.N_snaps - 1].satellites['M'].values * M_S / M_E)
            Inc.extend(sim.snaps[sim.N_snaps - 1].satellites['i'].values)
            Ecc.extend(sim.snaps[sim.N_snaps - 1].satellites['e'].values)

        cmap = 'cividis_r'
        cmin = min(Masses)
        cmax = max(Masses)

        norm = colors.LogNorm(cmin, cmax)

        fig, ax = plt.subplots(ncols=1)
        fig.set_size_inches(15.5, 10.5)
        ax.scatter(Ecc, np.sin(Inc), c=Masses, cmap=cmap, norm=norm)
        fig.colorbar(cm.ScalarMappable(cmap=cmap, norm=norm), orientation='vertical', label="Mass in Earth Masses",
                     ax=ax)
        ax.set_xlabel('Eccentricity')
        ax.set_ylabel('sin(Inclination)')
        ax.legend()
        fig.suptitle('Ecc and Inc Distribution')
        fig.savefig(path.join(self.PLOT, 'inc_ecc_distribution.png'))
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
