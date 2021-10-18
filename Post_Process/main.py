import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.animation import FuncAnimation
from snapshots.py import sat_data

test = snapshots.sat_data('system_0000/outputs')

sat = pd.read_csv("/Users/prut/Desktop/3DPop/system_0000/outputs/Snapshot_0100/satellites.txt", sep="	")

# print(sat.head())
# print(sat.wos)

sns.relplot(x="a", y="i", size="M", data=sat)
#plt.show()

snapshots = []
for i in range(101):
    path = "/Users/prut/Desktop/3DPop/system_0000/outputs/Snapshot_" + str(i).zfill(4) + "/satellites.txt"
    snap = pd.read_csv(path, sep="	", index_col="#ID")
    snap.sort_index(inplace=True)
    snapshots.append(snap)
print(snapshots[1].head)

# col = ['Mass','sma']
#
# sat_num = 1
# dic = {'Mass': [],
#        'sma': [],
#        'ecc': []}
# for i in range(101):
#     cur = snapshots[i]
#     cur = cur.loc[cur['#ID'] == sat_num]
#     #print(cur)
#     #id = snapshots[i].[sat_num,'#ID']
#     mass = cur.iat[0,1]
#     #print(mass)
#     sma = cur.iat[0,9]
#     #print(sma)
#     dic['Mass'].append(mass)
#     dic['sma'].append(sma)
#
# evol = pd.DataFrame.from_dict(dic)
# evol.index.name = 'index'
#
# sns.relplot(x="index", y="sma", size="Mass", data=evol)
# plt.show()

class importer:
    def __init__(self, dir, N):
        self.path = dir
        self.n_snaps = N
        self.raw_data = {}
        for i in range(self.n_snaps):
            path = self.path + "Snapshot_" + str(i).zfill(4) + "/satellites.txt"
            df = pd.read_csv(path, sep="	", index_col="#ID", dtype='float64')
            df.sort_index(inplace=True)
            self.raw_data[str(i)] = df
        self.n_sat = len(self.raw_data['0'])
        self.sats = {}
        for i in range(self.n_sat):
            df = pd.DataFrame(columns=self.raw_data['0'].columns,index=[int(x) for x in self.raw_data.keys()], dtype='float64')
            for key,value in self.raw_data.items():
                print(int(key))
                df.loc[int(key)] = value.loc[i+1]
            self.sats[str(i+1)] = df

test = importer("/Users/prut/Desktop/3DPop/system_0000/outputs/",101)

print(test.sats)

evol = test.sats['1']
print(evol.head)
dat = evol.iloc[1]
sns.relplot(x="a", y="e", size="M", data=dat)
plt.show()

