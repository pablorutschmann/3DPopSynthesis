import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib import colors
from matplotlib import patches
import Synthesis.post as post
from os import getcwd
import os.path as path

CWD = getcwd()

PATH = path.join(CWD, 'RUNTIME')
HISTORY = path.join(PATH, 'history.txt')
N_SIMS = 99

RUNTIME_POP = post.population(PATH, N_SIMS)


def get_runtimes():
    for i, sim in RUNTIME_POP.SIMS.items():
        query = 'system_{index}/inputs finished at'.format(index=i)
        with open(HISTORY, "r") as hist:
            for line in hist:
                if query in line:
                    # Month,Day,Time,Year = line.split()[-4:]
                    # Hours, Minutes, Seconds
                    runtime = float(next(hist).split()[1])
                    print(i, type(runtime))
                    sim.Runtime = runtime


get_runtimes()

columns = ['id', 'NPlanetesimals', 'Total Mass', 'Sigma Exponent', 'Runtime']
data = []

for i, sim in RUNTIME_POP.SIMS.items():
    row = [i, sim.N_Planetesimals, sim.Total_Mass, sim.Sigma_Exponent, sim.Runtime]
    data.append(row)

df = pd.DataFrame(columns=columns, data=data)
df.set_index(keys='id', inplace=True)

# print(df['Total Mass'])
# ax = df.plot.scatter(x='NPlanetesimals', y='Runtime', title= "Runtimes for different Number of Planetesimals")
# plot.show()

PLOT = path.join(PATH, 'runtimes_sigma.png')

fig, ax = plt.subplots(ncols=1)
fig.set_size_inches(15.5, 10.5)

cmap = 'turbo_r'
cmin = min(df["Total Mass"])
cmax = max(df["Total Mass"])
norm = colors.Normalize(cmin, cmax)

# extract all colors from the .jet map
cmaplist = ['red', 'green', 'blue', 'orange']

# create the new map
cmap = colors.LinearSegmentedColormap.from_list(
    'Custom cmap', cmaplist)

# define the bins and normalize
print(np.min(df['Sigma Exponent']))
print(np.max(df['Sigma Exponent']))
bounds = np.linspace(np.min(df['Sigma Exponent'].values), np.max(df['Sigma Exponent'].values), 5)

print(bounds)
norm = colors.BoundaryNorm(bounds, cmap.N)

# Scale the masses to the marker sizes
interval_min = 40
interval_max = 320
# ss = (df["Sigma Exponent"] - np.min(df["Sigma Exponent"])) / (np.max(df["Sigma Exponent"]) - np.min(df["Sigma Exponent"])) * (interval_max - interval_min) + interval_min

color_list = ['red', 'green', 'blue', 'orange']

ax.scatter(df["NPlanetesimals"], df["Runtime"], c=df['Sigma Exponent'], cmap=cmap, norm=norm)
fig.colorbar(cm.ScalarMappable(cmap=cmap, norm=norm), orientation='vertical', label="Sigma Exponent", ax=ax)
ax.set_xlabel('Number of Planetesimals', fontsize=15)
ax.set_ylabel('Runtime in hours', fontsize=15)
fig.suptitle('runtimes for different Parameters')
fig.savefig(PLOT)
plt.close(fig)
print('Plot saved!')





PLOT = path.join(PATH, 'runtimes_TM.png')

fig, ax = plt.subplots(ncols=1)
fig.set_size_inches(15.5, 10.5)

# extract all colors from the .jet map
cmaplist = ['red', 'green', 'blue', 'yellow', 'orange']

# create the new map
cmap = colors.LinearSegmentedColormap.from_list('Custom cmap', cmaplist)

# define the bins and normalize
bounds = np.linspace(np.min(df['Total Mass'].values), np.max(df['Total Mass'].values), 5)

print(bounds)
norm = colors.BoundaryNorm(bounds, cmap.N)

# Scale the masses to the marker sizes
interval_min = 40
interval_max = 320
# ss = (df["Sigma Exponent"] - np.min(df["Sigma Exponent"])) / (np.max(df["Sigma Exponent"]) - np.min(df["Sigma Exponent"])) * (interval_max - interval_min) + interval_min

color_list = ['red', 'green', 'blue', 'orange']

ax.scatter(df["NPlanetesimals"], df["Runtime"], c=df['Total Mass'], cmap=cmap, norm=norm)
fig.colorbar(cm.ScalarMappable(cmap=cmap, norm=norm), orientation='vertical', label="Total Disk Mass", ax=ax)
ax.set_xlabel('Number of Planetesimals', fontsize=15)
ax.set_ylabel('Runtime in hours', fontsize=15)
fig.suptitle('runtimes for different Parameters')
fig.savefig(PLOT)
plt.close(fig)
print('Plot saved!')
