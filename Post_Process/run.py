# run.py
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os.path
from collections import OrderedDict
import Post_Process as PP
from Post_Process.units import *

class run:
    def __init__(self, path, types = 'embryo'):
        self.path = path
        self.output_path = os.path.join(self.path, 'outputs')
        self.input_path = os.path.join(self.path, 'inputs')
        self.plot_path = os.path.join(self.path, 'plots')
        if not os.path.exists(self.plot_path):
            os.makedirs(self.plot_path)

        self.labels = {}
        self.labels['SigmaGas'] = ['Gas Surface Density', 'Sigma in [units]']
        self.labels['SigmaDust'] = ['Dust Surface Density', 'Sigma in [units]']
        self.labels['Temp'] = ['Temperature', 'Kelvin']

        # import collisions.txt
        path = os.path.join(self.output_path, 'collisions.txt')
        df = pd.read_csv(path, sep="	", index_col="#time", dtype='float64')
        df.sort_index(inplace=True)
        df.index = df.index.astype(float)
        self.collisions = df

        # import lost_satellites.txt
        path = os.path.join(self.output_path, 'lost_satellites.txt')
        df = pd.read_csv(path, sep="	", index_col="#ID", dtype='float64')
        df.index = df.index.astype(int)
        self.lost_satellites = df

        # import lost_satellites.txt
        path = os.path.join(self.output_path, 'satellite_list.txt')
        df = pd.read_csv(path, sep="	", index_col="#ID", dtype='float64')
        df.sort_index(inplace=True)
        df.index = df.index.astype(int)
        self.satellite_list = df

        # create dictionary of snapshots
        snaps_unordered = {}
        # iterate through list of Snapshot directory name and create class instance Snapshot and save it to dictionary
        for snap_path in [f.path for f in os.scandir(self.output_path) if (f.is_dir() and 'Snapshot' in f.path)]:
            shot = PP.snapshot(snap_path)
            snaps_unordered[shot.index] = shot
        self.snaps = OrderedDict(sorted(snaps_unordered.items(), key=lambda t: t[0]))

        # transform data to satellite specific
        self.satellites = {}
        self.embryos = {}
        self.planetesimals = {}
        self.get_satellites(types)


    def get_satellites(self, types = 'embryo'):
        for index in self.satellite_list.index:
            if types == 'embryo' and self.satellite_list.loc[index, 'Type'] == 1:
                self.embryos[index] = PP.satellite(index, self)
                continue
            elif types == 'planetesimal' and self.satellite_list.loc[index, 'Type'] == 0:
                self.planetesimals[index] = PP.satellite(index, self)
                continue
            elif types == 'all':
                if self.satellite_list.loc[index, 'Type'] == 0:
                    self.planetesimals[index] = PP.satellite(index, self)
                elif self.satellite_list.loc[index, 'Type'] == 1:
                    self.embryos[index] = PP.satellite(index, self)

        self.satellites.update(self.embryos)
        self.satellites.update(self.planetesimals)

    def plot_snapshots(self, ids=None):
        if ids == None:
            ids = [list(self.snaps.keys())[0], list(self.snaps.keys())[-1]]
        for i in ids:
            self.snaps[i].plot_satellites()
            self.snaps[i].plot_satellites_ratio()

    def plot_disk_evol(self, field, log=True, N=10, keys=None):
        if keys == None:
            keys = list(self.snaps.keys())

        if N < len(keys):
            idx = np.round(np.linspace(0, len(keys) - 1, N)).astype(int)
            keys = [keys[x] for x in idx]

        title, units = self.labels[field]

        fig, ax = plt.subplots(ncols=1)
        fig.set_size_inches(15.5, 10.5)

        for key in keys:
            fig, ax = self.snaps[key].fig_disk(fig, ax, field)
        if field == 'Temp':
            xmin = rtoau * self.snaps[keys[0]].disk['r'].min()
            xmax = rtoau * self.snaps[keys[0]].disk['r'].max()
            ax.hlines(170, xmin, xmax, linestyles='dotted')
        ax.set_xlabel('Radius in AU', fontsize=15)
        ax.set_ylabel(units, fontsize=15)
        if log:
            ax.set_xscale('log')
            ax.set_yscale('log')
        ax.legend(title='Time [Myrs]')
        fig.suptitle(title + ' Evolution')
        filename = os.path.join(self.plot_path, title.replace(" ", "_") + '.png')
        fig.savefig(filename)
        plt.close(fig)
        print('Plot saved at: ' + os.path.join(self.plot_path, filename))

    def plot_disk_evol_all(self, N=10):
        self.plot_disk_evol(field="SigmaGas", N=N)
        self.plot_disk_evol(field="SigmaDust", N=N)
        self.plot_disk_evol(field="Temp", N=N)

    def plot_accretion(self, types = 'embryo'):
        fig, ax = plt.subplots()
        ax.set_xlabel('time in years')
        ax.set_ylabel('mass in [Units]')
        ax.set_title('Mass Evolution of remaining satellites')

        for item in self.satellites.values():
            if types == 'embryo' and item.type == 1:
                fig, ax = item.fig_accretion(fig, ax)
                continue
            elif types == 'planetesimal' and item.type == 0:
                fig, ax = item.fig_accretion(fig, ax)
                continue
            else:
                fig, ax = item.fig_accretion(fig, ax)

        filename = 'accretion_' + types + '.png'
        savepath = os.path.join(self.plot_path, filename)
        fig.savefig(savepath)
        plt.close(fig)

    def plot_wm(self):
        for item in self.embryos.values():
            item.wm_time(self.plot_path)



if __name__ == "__main__":
    test = run('/Users/prut/CLionProjects/3DPopSynthesis/Runs/testembryo')
    # test.satellites.
    print(test.embryos[7].acc)
    print(test.embryos[7].data['a'])
    print(rtoau * 1547.68)

    for key, item in test.snaps.items():
        rad = rtoau * item.IceLineRadius
        print(rad)
