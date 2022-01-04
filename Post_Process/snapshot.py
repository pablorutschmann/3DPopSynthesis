import numpy as np
import pandas as pd
import os.path
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib import colors
from matplotlib import patches
from Post_Process.units import *


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

        self.Mass = np.sum(self.satellites['M']) * M_S
        #pi (R^2 - r^2)
        self.Area = np.pi * ( max(self.satellites['a'])**2 - min(self.satellites['a'])**2) * R_S2

        self.Sigma = self.Mass / self.Area

    def plot_satellites(self):
        fig, ax = plt.subplots(ncols=1)
        fig.set_size_inches(15.5, 10.5)

        cmap = 'turbo_r'
        WMF_satellites = self.satellites['WM'] / self.satellites['M']
        cmin = min(WMF_satellites)
        cmax = max(WMF_satellites)
        norm = colors.Normalize(cmin, cmax)

        # Scale the masses to the marker sizes
        interval_min = 40
        interval_max = 320
        ss = (self.satellites['M'] - np.min(self.satellites['M'])) / (np.max(self.satellites['M']) - np.min(self.satellites['M'])) * (interval_max - interval_min) + interval_min

        ax.scatter(rtoau * self.satellites['a'], self.satellites['e'], c=WMF_satellites, cmap=cmap, norm=norm, s=ss, alpha=1)
        fig.colorbar(cm.ScalarMappable(cmap=cmap, norm=norm), orientation='vertical', label="Water Mass Fraction",ax=ax)
        ax.set_ylim(-1 * min(self.satellites['e']), 1.1 * max(self.satellites['e']))
        ax.set_xlabel('Semi Mayor Axis in AU', fontsize=15)
        ax.set_ylabel('Eccentricity', fontsize=15)
        fig.suptitle('Satellites at time ' + str(self.time) + ' Myrs')
        fig.savefig(self.plot_path + '/satellites' + str(self.index).zfill(4) + '.png')
        plt.close(fig)
        print('Plot saved at: ' + self.plot_path + '/satellites' + str(self.index).zfill(4) + '.png')

    def plot_histogram(self):
        fig, ax = plt.subplots()
        ax.set_xlabel('Radius in au')
        ax.set_ylabel('Counts')
        ax.set_title('Distribution of Satellites')

        ax.hist(self.satellites['a']*R_S/au, bins=50)
        savepath = os.path.join(self.plot_path, 'histogram.png')
        fig.savefig(savepath)
        plt.close(fig)

    def plot_satellites_ratio(self):
        fig, ax = plt.subplots(ncols=1)
        fig.set_size_inches(15.5, 10.5)

        colors = ['red', 'green', 'blue']
        labels = ['Solids', 'Hydrated Silica', 'Water/Ice']

        red_patch = patches.Patch(color='red', label='Solids')
        green_patch = patches.Patch(color='green', label='Hydrated Silica')
        blue_patch = patches.Patch(color='blue', label='Water/Ice')
        handles = [red_patch, green_patch, blue_patch]

        mean_mass = np.min(mtome * self.satellites['M'])
        mass_scaling = mean_mass / 90000
        mass_scaling = 0.000000001

        def pie_1d(r1, r2):
            # calculate the points of the first pie marker
            # these are just the origin (0, 0) + some (cos, sin) points on a circle
            x1 = np.cos(2 * np.pi * np.linspace(0, r1))
            y1 = np.sin(2 * np.pi * np.linspace(0, r1))
            xy1 = np.row_stack([[0, 0], np.column_stack([x1, y1])])
            s1 = np.abs(xy1).max()

            x2 = np.cos(2 * np.pi * np.linspace(r1, r2))
            y2 = np.sin(2 * np.pi * np.linspace(r1, r2))
            xy2 = np.row_stack([[0, 0], np.column_stack([x2, y2])])
            s2 = np.abs(xy2).max()

            x3 = np.cos(2 * np.pi * np.linspace(r2, 1))
            y3 = np.sin(2 * np.pi * np.linspace(r2, 1))
            xy3 = np.row_stack([[0, 0], np.column_stack([x3, y3])])
            s3 = np.abs(xy3).max()

            return xy1, s1, xy2, s2, xy3, s3

        # cale the masses to the marker sizes
        def point_scaling(row):
            interval_min = 40
            interval_max = 320
            ss = (row - np.min(self.satellites['M'])) / (np.max(self.satellites['M']) - np.min(self.satellites['M'])) * (interval_max - interval_min) + interval_min
            return ss

        def plot_one(row):
            WMF_ratio = row['WM'] / row['M']
            SWMF_ratio = row['SWM'] / row['M'] + WMF_ratio
            xy1, s1, xy2, s2, xy3, s3 = pie_1d(WMF_ratio, SWMF_ratio)
            scale = point_scaling(row['M'])

            ax.scatter(rtoau * row['a'], row['e'], marker=xy1, s=s1 ** 2 * scale, facecolor='blue')
            ax.scatter(rtoau * row['a'], row['e'], marker=xy2, s=s2 ** 2 * scale, facecolor='green')
            ax.scatter(rtoau * row['a'], row['e'], marker=xy3, s=s3 ** 2 * scale, facecolor='red')

        for index, row in self.satellites.iterrows():
            plot_one(row)

        ax.set_ylim(-1 * min(self.satellites['e']), 1.1 * max(self.satellites['e']))
        ax.set_xlabel('Semi Mayor Axis in AU', fontsize=15)
        ax.set_ylabel('Eccentricity', fontsize=15)
        ax.legend(handles=handles, title='Components')
        fig.suptitle('Satellites at time ' + str(self.time) + ' Myrs')
        fig.savefig(self.plot_path + '/satellites_ratios_' + str(self.index).zfill(4) + '.png')
        plt.close(fig)
        print('Plot saved at: ' + self.plot_path + '/satellites_ratios_' + str(self.index).zfill(4) + '.png')

    def fig_disk(self, fig, ax, field):
        ax.plot(rtoau * self.disk['r'], self.disk[field], label=self.time)

        return fig, ax

# Helper Functions

def get_parent(path):
    return os.path.normpath(os.path.join(path, os.pardir))


if __name__ == "__main__":
    test = snapshot('/Users/prut/CLionProjects/3DPopSynthesis/Runs/testsmf/outputs/Snapshot_0009')
    # # test.plot()
    #
    test.plot_satellites_ratio()
