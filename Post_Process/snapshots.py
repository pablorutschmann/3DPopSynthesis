import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib import colors
from pathlib import Path
import astropy.constants as astroconst
from astropy import units as u

# AU in cm
au = astroconst.au.decompose(u.cgs.bases).value

# Solar Radius in grams
R_S = astroconst.R_sun.decompose(u.cgs.bases).value

# Solar Radius squared
R_S2 = R_S ** 2

# Solar Mass in grams
M_S = astroconst.M_sun.decompose(u.cgs.bases).value

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
        all_init_sma = first_snap['a'] * R_S / au
        all_init_e = first_snap['e']
        all_init_wmf = first_snap['WMF']

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
        final_sma = data_rem[str(self.numbers[-1])]['a'] * R_S / au
        final_e = data_rem[str(self.numbers[-1])]['e']
        final_wmf = data_rem[str(self.numbers[-1])]['WMF']
        final_s = final_masses/max(all_init_masses)

        init_masses = data_rem[str(self.numbers[0])]['M']
        init_sma = data_rem[str(self.numbers[0])]['a']* R_S / au
        init_e = data_rem[str(self.numbers[0])]['e']

        init_s = init_masses/max(final_masses)
        all_init_s = all_init_masses/max(all_init_masses)
        init_sma = (init_sma / max(init_sma)) * R_S / au
        scaling = 20

        #norm = plt.Normalize(0,0.5)

        #cmap = colors.LinearSegmentedColormap.from_list("", ["red","yellow","blue"], N=1000)
        cmap = 'turbo_r'

        cmin = min([min(all_init_wmf),min(final_wmf)])
        cmax = max([max(all_init_wmf),max(final_wmf)])
        print(cmin)
        print(final_wmf)
        norm = colors.Normalize(cmin,cmax)

        ax1.scatter(all_init_sma, all_init_e, c=all_init_wmf, cmap=cmap, norm=norm, s=all_init_s*scaling, alpha=1)
        #ax1.scatter(init_sma, init_e, c=init_sma, cmap='cividis', s=init_s * scaling, alpha=1)
        ax2.scatter(final_sma,final_e, c=final_wmf, cmap=cmap, norm=norm, s=final_s*scaling, alpha=1)
        fig.colorbar(cm.ScalarMappable(cmap=cmap,norm=norm), orientation='vertical', label="Water Mass Fraction", ax=[ax1,ax2])
        ax1.set_ylim(-0.05,0.4)
        ax2.set_ylim(-0.05, 0.4)
        ax1.set_xlabel('Semi Mayor Axis in AU', fontsize=15)
        ax1.set_ylabel('Eccentricity', fontsize=15)
        ax2.set_xlabel('Semi Mayor Axis in AU', fontsize=15)
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

class disk_data:
    def __init__(self, dir):
        self.dir = dir
        self.path = self.dir + 'outputs/'

        Path(self.dir + 'plots/' ).mkdir(parents=False, exist_ok=True)

        labelnames = ['Radius', 'Delta Radius', 'Surface Density Gas', 'Surface Density Dust', 'Initial Surface Density Dust', 'Temperature', 'Annulus Area','Keplerian Velocity', 'Surface Density Exponent', 'Temperature Exponent', 'Opacity']
        units = ['AU', 'AU', r'$M_s R_s^{-2}$', r'$M_s R_s^{-2}$', r'$M_s R_s^{-2}$', r'$Kelvin$', r'$R_s^{2}$', 'Velocity', 'unitless', 'unitless', 'unitless']
        names = ['r', 'dr', 'SigmaGas', 'SigmaDust', 'SigmaDustBar', 'Temp', 'Area','OmegaK', 'SigmaExponent', 'TExponent', 'Opacity']
        self.labels = {}
        self.labels['unit'] = units
        self.labels['label'] = labelnames
        self.labels['name'] = names

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
        self.times = []
        for i in self.numbers:
            path = self.path + "Snapshot_" + str(i).zfill(4)
            diskfile = path + "/disk.txt"
            df = pd.read_csv(diskfile, sep="	", index_col="#index", dtype='float64')
            df.sort_index(inplace=True)
            df.index = df.index.astype(int)
            d = {}
            parafile = path + "/parameters.txt"
            with open(parafile) as f:
                for line in f:
                    (key, val) = line.split()
                    d[key] = val
            self.times.append(d['Time'])
            self.raw_data[str(i)] = df


    def temp(self,r):
        return (5780 * (1 / r)**0.5 * 0.7520883995742579)


    def plot_evol(self,name):

        index = self.labels['name'].index(name)
        name = self.labels['name'][index]
        unit = self.labels['unit'][index]
        label = self.labels['label'][index]

        fig, ax = plt.subplots()

        rs = self.raw_data[str(self.numbers[0])]['r'] * R_S / au

        data = []
        for key, item in self.raw_data.items():
            data.append(item[name])
        #ax.plot(rs,self.temp(rs))
        for i, item in enumerate(data):
            ax.plot(rs,item, label=str(self.times[i]))
        ax.set_xlabel('Radius in AU', fontsize=15)
        ax.set_ylabel(unit, fontsize=15)
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.legend()
        fig.suptitle(label + 'Profile Evolution')

        #plt.show()
        fig.savefig(self.dir + 'plots/' + label + '_profile_evolution_' + '.png')



