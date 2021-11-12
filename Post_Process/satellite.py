import pandas as pd
import matplotlib.pyplot as plt
from Post_Process.units import *


class satellite:
    def __init__(self, ID, run):
        # Check if valid ID
        if ID in run.satellite_list.index:
            self.ID = ID
        else:
            raise Exception('Satellite ' + str(ID) + ' does not exist!')\


        self.type = run.satellite_list.loc[ID, 'Type']

        self.init_time = run.satellite_list.loc[ID, 'init_time']
        self.init_temp = run.satellite_list.loc[ID, 'init_temp']

        # check whether it was lost
        if ID in run.lost_satellites.index:
            self.lost = True
        else:
            self.lost = False

        # get data from each snapshot
        columns = run.snaps[list(run.snaps)[0]].satellites.columns
        data = {k: [] for k in columns}
        data['time'] = []
        for key, item in run.snaps.items():
            if self.ID in item.satellites.index:
                for col in columns:
                    data[col].append(item.satellites.loc[ID, col])
                data['time'].append(float(item.time))
        self.data = pd.DataFrame.from_dict(data)
        self.data.set_index('time', inplace=True)

        # Get Accretion Data
        self.acc = self.data[['M', 'WM']]
        self.acc = self.acc.rename(columns={"M": "mass", "WM": "wm"})

        for time, row in run.collisions.iterrows():
            if row['ID'] == self.ID:
                # print(row)
                new_row = row[['mass', 'wm']]
                new_row.index.name = 'time'
                self.acc = self.acc.append(new_row)
        self.acc.sort_index(inplace=True)

    def fig_accretion(self, fig, ax):
        ax.step(list(self.acc.index), mtome * self.acc['mass'], where='post')
        return fig, ax

    def wm_time(self, path):
        fig, ax = plt.subplots(ncols=1)
        fig.set_size_inches(15.5, 10.5)

        ax.scatter(self.data.index, self.data['WM']/self.data['M'])

        ax.set_xlabel('Time in Years', fontsize=15)
        ax.set_ylabel('Water Mass', fontsize=15)
        fig.suptitle('Water Mass Time Evolution of Satellite ' + str(self.ID))
        fig.savefig(path + '/satellite' + str(self.ID) + '_wm.png')
        plt.close(fig)
        print('Plot saved at: ' + path + '/satellite' + str(self.ID) + '_wm.png')


if __name__ == "__main__":
    print('hello')
