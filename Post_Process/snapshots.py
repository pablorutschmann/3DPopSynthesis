import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib import colors
from pathlib import Path

class sat_data:
    def __init__(self, dir):
        self.dir = dir
        self.path = self.dir + 'outputs/'

        Path(self.dir + 'plots/' ).mkdir(parents=False, exist_ok=True)

        #get number of snapshots
        N_snaps = 0
        self.numbers = []
        for file in os.listdir(self.path):
                if "Snapshot" in file:
                    number = int(file[-4:])
                    self.numbers.append(number)
                    N_snaps += 1
        self.n_snaps = N_snaps
        self.numbers.sort()
        self.raw_data = {}
        for i in self.numbers:
            path = self.path + "Snapshot_" + str(i).zfill(4) + "/satellites.txt"
            df = pd.read_csv(path, sep="	", index_col="#ID", dtype='float64')
            df.sort_index(inplace=True)
            df.index = df.index.astype(int)
            self.raw_data[str(i)] = df

    def plot_remaining(self):
        last_snap = self.raw_data[str(self.numbers[-1])]
        first_snap = self.raw_data[str(self.numbers[0])]

        all_init_masses = first_snap['M']
        all_init_sma = first_snap['a']
        all_init_e = first_snap['e']

        #ids of the remaining satellites
        id_rem = list(last_snap.index.values)
        id_rem = [id for id in id_rem if id < 21]

        #get data of remaining satellites
        data_rem = {}
        for snap,df in self.raw_data.items():
            data_rem[snap] = df.loc[id_rem, :]

        #plot initial conditions and final product
        fig, (ax1, ax2) = plt.subplots(ncols=2)
        fig.set_size_inches(18.5, 10.5)
        final_masses = data_rem[str(self.numbers[-1])]['M']
        final_sma = data_rem[str(self.numbers[-1])]['a']
        final_e = data_rem[str(self.numbers[-1])]['e']
        final_s = final_masses/max(all_init_masses)

        init_masses = data_rem[str(self.numbers[0])]['M']
        init_sma = data_rem[str(self.numbers[0])]['a']
        init_e = data_rem[str(self.numbers[0])]['e']

        init_s = init_masses/max(final_masses)
        all_init_s = all_init_masses/max(all_init_masses)
        init_sma = (init_sma / max(init_sma))
        scaling = 5

        cmap = cm.cividis
        bounds = init_sma
        norm = colors.Normalize(vmin=min(all_init_sma), vmax=max(all_init_sma))


        ax1.scatter(all_init_sma, all_init_e, c=all_init_sma, cmap='cividis', s=all_init_s*scaling, alpha=1)
        #ax1.scatter(init_sma, init_e, c=init_sma, cmap='cividis', s=init_s * scaling, alpha=1)
        ax2.scatter(final_sma,final_e, c=init_sma, cmap='cividis', s=final_s*scaling, alpha=1)
        fig.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap),
                     orientation='vertical',
                     label="starting sma", ax=[ax1,ax2])
        # ax1.set_ylim(-0.05,0.4)
        # ax2.set_ylim(-0.05, 0.4)
        ax1.set_xlabel('Semi Mayor Axis in J_r', fontsize=15)
        ax1.set_ylabel('Eccentricity', fontsize=15)
        ax2.set_xlabel('Semi Mayor Axis in J_r', fontsize=15)
        ax2.set_ylabel('Eccentricity', fontsize=15)
        fig.suptitle('Final Satellites and their starting Conditions')

        #plt.show()
        fig.savefig(self.dir + 'plots/' + 'sma.png')


    def accretion(self):
        path = self.path + "/collisions.txt"
        coll = pd.read_csv(path, sep="	", index_col="#time", dtype='float64')
        #print(coll.columns)
        #print(coll[['ID1','mass1','ID2','mass2','ID']])

        last_snap = self.raw_data[str(self.n_snaps - 1)]

        # ids of the remaining satellites
        id_rem = list(last_snap.index.values)
        max_time = 0
        fig, ax = plt.subplots()
        for rem in id_rem:
            acc = coll.loc[coll['ID'] == rem]
            times = [0,]
            times += list(acc.index.values)
            times.append(70000)
            masses=[]
            if acc.iloc[0]['ID1'] == rem:
                print(acc.iloc[0]['mass1'])
                masses.append(acc.iloc[0]['mass1'])
            elif acc.iloc[0]['ID2'] == rem:
                masses.append(acc['mass2'].iloc[0])
            masses += list(acc['mass'].values)
            masses.append(masses[-1])
            print(masses)



            ax.set_xlabel('time in years')
            ax.set_ylabel('mass in M_j')
            ax.set_title('Mass Evolution of remaining satellites')
            ax.step(times, masses, where='post')

        #plt.show()
        fig.savefig(self.dir + 'plots/' + 'accretion.png')






