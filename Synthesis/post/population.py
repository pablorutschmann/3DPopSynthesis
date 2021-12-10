import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os.path as path
from os import makedirs
from ..units import *
from .simulation import simulation
import pickle

class population:
    def __init__(self, PATH, N_sims):
        self.MAIN = PATH
        self.PLOT = path.join(self.MAIN, 'plots')
        self.NSIMS = N_sims
        self.SIMS = {}
        if not path.exists(self.PLOT):
            makedirs(self.PLOT)

        for i in range(1,self.NSIMS):
            self.SIMS[i] = simulation(self.system_path(i))
        print(self.SIMS)


    def final_mass_distribution(self):
        Masses = []

        for sim in self.SIMS.values():
            Masses.extend(sim.snaps[sim.N_snaps-1].satellites['M'].values * M_S / M_E)

        print(len(Masses))
        fig, ax = plt.subplots(ncols=1)
        fig.set_size_inches(15.5, 10.5)
        ax.hist(Masses, bins=100, density=False, range=(10e-6,max(Masses)))
        fig.suptitle('Final Mass Distribution')
        fig.savefig(path.join(self.PLOT, 'final_mass_distribution.png'))
        plt.close(fig)

    def system_path(self,i):
        sys = "system_" + str(i)
        return path.join(self.MAIN,sys)

    def system_output_path(self,i):
        return path.join(self.system_path(i),'outputs')

    def system_output_path(self,i):
        return path.join(self.system_path(i),'inputs')


if __name__ == "__main__":
    TEST = population("test", 9)