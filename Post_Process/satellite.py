import pandas as pd

class satellite:
    def __init__(self, ID, run):
        # Check if valid ID
        if ID in run.satellite_list.index:
            self.ID = ID
        else:
            raise Exception('Satellite ' + str(ID) + ' does not exist!')

        self.init_time = run.satellite_list.loc[ID,'init_time']
        self.init_temp = run.satellite_list.loc[ID,'init_temp']

        # check whether it was lost
        if ID in run.lost_satellites.index:
            self.lost = True
        else:
            self.lost = False

        # get data from each snapshot
        columns = run.snaps[list(run.snaps)[0]].satellites.columns
        data = {k:[] for k in columns}
        data['time'] = []
        for key,item in run.snaps.items():
            if self.ID in item.satellites.index:
                for col in columns:
                    data[col].append(item.satellites.loc[ID,col])
                data['time'].append(item.time)
        self.data = pd.DataFrame.from_dict(data)
        self.data.set_index('time',inplace=True)

        # Get Accretion Data



if __name__ == "__main__":
    print('hello')