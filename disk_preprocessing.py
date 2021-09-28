import os, sys
import numpy as np
from scipy.optimize import curve_fit
import astropy.constants as const
from os import listdir
from os.path import isfile, join
from math import floor
from tqdm import tqdm
import itertools
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit



au = (const.au.value * 100)

R_S = const.R_sun.value * 100
print(au/R_S)
def import_ppd(path):
    columns = ['x [cm]', 'y [cm]', 'z [cm]', 'dens [g/cm3]', 'temp [K]', 'velo_x [cm/s]', 'velo_y [cm/s]', 'velo_z [cm/s]', 'opac[cm^2/g]']
    disk = pd.read_csv(path, delim_whitespace=True, names=columns, index_col=False)
    # disk['x [cm]'] /= au
    # disk['y [cm]'] /= au
    # disk['z [cm]'] /= au
    radius, phi = cart2pol(disk['x [cm]'],disk['y [cm]'])
    disk['r [AU]'] = radius
    disk['phi [rad]'] = phi + np.pi

    return disk

def cart2pol(x, y):
    radius = np.sqrt(x**2 + y**2)
    phi = np.arctan2(y, x)
    return(radius, phi)


#dust = import_ppd('ppd/lev0_dust.dat')
gas = import_ppd('ppd/lev0_gas.dat')
#print(dust.head(5))
r_resh = np.reshape(gas['r [AU]'].values, [680,215,20])
print(r_resh[1,:,1]-r_resh[2,:,2])
# print()
# #dust = dust.sort_values(by=['r [AU]'])
# print(len(dust['phi [rad]'].unique()))
# #print(gas.head(5))
#
# angles = np.linspace(0,2*np.pi,680)
# #print(angles)
#
# grid = np.around(dust['phi [rad]'].unique(),4)



def integrate_1D(df,plot=False):
    N_f = 5
    N_b = 5

    M_sun = const.M_sun.value * 1000
    R_sun2 = (const.R_sun.value * 100)**2

    def truncate(arr,N_front,N_back):
        trunc = arr[N_front:-N_back]
        return trunc

    # density
    dens = np.reshape(df['dens [g/cm3]'].values, [680,215,20])
    sigma2d = np.sum(dens, axis=2)
    sigma1d = truncate(np.mean(sigma2d, axis=(0)) / M_sun * R_sun2,N_f,N_b)

    #opacity
    opac = np.reshape(df['opac[cm^2/g]'].values, [680,215,20])
    opac1d = truncate(np.mean(opac,axis=(0,2)),N_f,N_b)

    #temperature
    temp = np.reshape(df['temp [K]'].values, [680,215,20])
    temp1d = truncate(np.mean(temp,axis=(0,2)),N_f,N_b)

    rs = np.reshape(df['r [AU]'].values, [680,215,20])
    rs = truncate(rs[0,:,0],N_f,N_b)/R_S

    def func(x, a, b):
        return a * x + b

    popt, pcov = curve_fit(func, np.log(rs), np.log(sigma1d))

    x = np.linspace(0.5,5,500)* const.au.value / const.R_sun.value
    if plot == True:

        fig, ax = plt.subplots()

        ax.set_xlabel('radius in AU')
        ax.set_ylabel('surface density in g / cm^2')
        ax.set_title('gas surface density profile')
        #ax.loglog(r_resh[0,:,0]/au, sigma1d)
        ax.plot(np.log(rs), np.log(sigma1d), 'ko', label="Original Data")
        ax.plot(np.log(x), func(np.log(x), *popt), 'r-', label="Fitted Curve")
        plt.show()

    return sigma1d,opac1d,temp1d,x

def write_disk(df,dim=1):
    if dim == 1:
        sigma,_,temp = integrate_1D(df)
    elif dim == 2:
        sigma,_,temp = integrate_2D(df)
    else:
        raise Exception('Not Possible Intergration')

    rs = np.reshape(df['r [AU]'].values, [680,215,20])
    rs = rs[0,:,0]/R_S

    out = {}
    out['#r [R_s]'] = rs
    out['sigma [M_S/R_S^2]'] = sigma
    out['T [K]'] = temp

    df_out = pd.DataFrame.from_dict(out)
    print(df_out['sigma [M_S/R_S^2]'])

    df_out.to_csv(r'disk215.txt', index=None, sep=' ', mode='a')




a = integrate_1D(gas,True)

#a = write_disk(gas)


#def one_diskmodel(df):



#def two_diskmodel(df):