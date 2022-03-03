import matplotlib.pyplot as plt
import os.path as path
from os import makedirs
from .simulation import simulation
from tqdm import tqdm

class population:
    def __init__(self, PATH, N_sims):
        self.MAIN = PATH
        self.NAME = path.basename(self.MAIN)
        self.PLOT = path.join(self.MAIN, 'plots')
        self.NSIMS = N_sims
        self.SIMS = {}
        if not path.exists(self.PLOT):
            makedirs(self.PLOT)

        for i in tqdm(range(1, self.NSIMS + 1)):
            self.SIMS[i] = simulation(self.system_path(i))

        print(self.SIMS)

        self.plot_config = 'paper'
        self.figsize = (8, 6)
        self.fontsize = 60
        self.legend_fontsize = 9
        self.dpi = 600
        self.cmap_standart = 'viridis_r'

        print(self)

    def switch_plot_config(self, config):
        if config == 'paper':
            self.plot_config = 'paper'
            self.figsize = (8, 6)
            self.fontsize = 80
            self.legend_fontsize = 9
            self.dpi = 600
            plt.rcParams.update({"legend.title_fontsize": self.legend_fontsize})
            plt.rcParams.update({'font.size': self.fontsize})
            print(f'Plot configuration changed to {config}.')

        elif config == 'presentation':
            self.plot_config = 'presentation'
            self.figsize = (8, 7)
            self.fontsize = 130
            self.legend_fontsize = 11
            self.dpi = 600
            plt.rcParams.update({"legend.title_fontsize": self.legend_fontsize})
            plt.rcParams.update({'font.size': self.fontsize})
            print(f'Plot configuration changed to {config}.')
        else:
            print(f'Option {config} not found.')

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
