import numpy as np
import pandas as pd
import os.path
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib import colors
from .units import *


class snapshot:
    def __init__(self, path):
        self.path = path

        self.plot_path = os.path.join(get_parent(get_parent(self.path)), 'plots')

        # import satellites.txt
        path = self.path + "/satellites.txt"
        df = pd.read_csv(path, sep="	", index_col="#ID", dtype='float64')
        df.sort_index(inplace=True)
        df.index = df.index.astype(int)
        self.satellites = df

        # import parameters.txt
        path = self.path + "/parameters.txt"
        dic = {}
        file = open(path)
        for line in file:
            key, value = line.split()
            dic[key] = value
        self.time = dic['Time']
        self.updatetime = dic['UpdateTime']
        self.index = int(dic['SaveIndex'])
        self.timestopformation = dic['TimeStopFormation']
        self.globaldt = dic['GlobalDt']

        # import disk.txt
        path = self.path + "/disk.txt"
        df = pd.read_csv(path, sep="	", index_col="#index", dtype='float64')
        df.sort_index(inplace=True)
        df.index = df.index.astype(int)
        self.disk = df

    def plot_satellites(self):
        fig, ax = plt.subplots(ncols=1)
        fig.set_size_inches(15.5, 10.5)

        cmap = 'turbo_r'
        WMF_satellites = self.satellites['WM']/self.satellites['M']
        cmin = min(WMF_satellites)
        cmax = max(WMF_satellites)
        norm = colors.Normalize(cmin, cmax)
        mean_mass = np.min(mtome * self.satellites['M']) / 50000

        ax.scatter(rtoau * self.satellites['a'], self.satellites['e'], c=WMF_satellites, cmap=cmap, norm=norm,
                   s=self.satellites['M'] / mean_mass, alpha=1)
        # ax1.scatter(init_sma, init_e, c=init_sma, cmap='cividis', s=init_s * scaling, alpha=1)
        fig.colorbar(cm.ScalarMappable(cmap=cmap, norm=norm), orientation='vertical', label="Water Mass Fraction",
                     ax=ax)
        ax.set_ylim(-1 * min(self.satellites['e']), 1.1 * max(self.satellites['e']))
        ax.set_xlabel('Semi Mayor Axis in AU', fontsize=15)
        ax.set_ylabel('Eccentricity', fontsize=15)
        fig.suptitle('Satellites at time ' + str(self.time) + ' Myrs')
        fig.savefig(self.plot_path + '/satellites' + str(self.index).zfill(4) + '.png')
        print('Plot saved at: ' + self.plot_path + '/satellites' + str(self.index).zfill(4) + '.png')

    def fig_disk(self, fig, ax, field):
        ax.plot(rtoau * self.disk['r'], self.disk[field], label=self.time)

        return fig, ax


# Helper Functions

def get_parent(path):
    return os.path.normpath(os.path.join(path, os.pardir))


if __name__ == "__main__":
    test = snapshot('/Users/prut/CLionProjects/3DPopSynthesis/Runs/fulltest/outputs/Snapshot_0009')
    # test.plot()

    test.plot_satellites()
