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

#Defining Constants

# AU in cm
au = (const.au.value * 100)

# Solar Radius in cm
R_S = const.R_sun.value * 100

# Solar Radius squared
R_sun2 = (const.R_sun.value * 100)**2

# Solar Mass in grams
M_sun = const.M_sun.value * 1000

# Earth Masses
M_E = const.M_earth

# Earth Radius
R_E = const.R_earth

print(1.5 * 10**(-21) * M_sun / R_sun2 * 202799192448.87842)

# Defining Parameters

phi_cells = 680
r_cells = 215
th_cells = 20

# Truncate Cells
N_f = 7
N_b = 6


# New Disk Radii Limits
R_min = 0.5
R_max = 5
N = 1000

#Helper Functions
def cart2pol(x, y):
    radius = np.sqrt(x**2 + y**2)
    phi = np.arctan2(y, x)
    return(radius, phi)


# Import Function, Saves Disk to Dataframe in cgs units
def import_ppd(path):
    # Initial Columns
    columns = ['x [cm]', 'y [cm]', 'z [cm]', 'dens [g/cm3]', 'temp [K]', 'velo_x [cm/s]', 'velo_y [cm/s]', 'velo_z [cm/s]', 'opac[cm^2/g]']

    # Read the file
    disk = pd.read_csv(path, delim_whitespace=True, names=columns, index_col=False)

    # Convert to cylindrcal coordinates
    radius, phi = cart2pol(disk['x [cm]'],disk['y [cm]'])
    disk['r [AU]'] = radius
    disk['phi [rad]'] = phi + np.pi

    return disk


def integrate_1D(name,space,df,plot=False):

    def truncate(arr,N_front,N_back):
        trunc = arr[N_front:-N_back]
        return trunc

    # density
    dens = np.reshape(df['dens [g/cm3]'].values, [phi_cells,r_cells,th_cells])
    #height one cell 202799192448.87842


    #Integrate in z direction
    sigma2d = np.sum(dens, axis=2) * 203599192448.87842
    #Average azimutal direction
    sigma1d = truncate(np.mean(sigma2d, axis=(0)),N_f,N_b)
    #convert to Solar Radius and Mass units
    sigma1d = sigma1d / M_sun * R_sun2

    # opacity
    opac = np.reshape(df['opac[cm^2/g]'].values, [phi_cells,r_cells,th_cells])
    # average in vertical and azimutal directions
    opac1d = truncate(np.mean(opac,axis=(0,2)),N_f,N_b)
    #convert to Solar Radius and Mass units
    opac1d = opac1d * M_sun / R_sun2

    # temperature
    temp = np.reshape(df['temp [K]'].values, [phi_cells,r_cells,th_cells])
    # average in vertical and azimutal directions
    temp1d = truncate(np.mean(temp,axis=(0,2)),N_f,N_b)

    # Radii of data values
    rs = np.reshape(df['r [AU]'].values, [phi_cells,r_cells,th_cells])
    print(rs.min()/au)
    print(rs.max()/au)
    rs = truncate(rs[0,:,0],N_f,N_b)/R_S


    r_all = np.linspace(R_max,rs.max()*R_S / au, N) * au / R_S

    # Extrapolation linspace or logspace of radii in R_S
    if space == "lin":
        x = np.linspace(R_min,R_max,N) * au / R_S
    elif space == 'log':
        x = np.logspace(np.power(R_min,(1/np.e)), np.power(R_max,(1/np.e)), N, base=np.e) * au / R_S

    def extrapolate(type,x_init,y_init):

        def func(x, a, b):
            return a * x + b

        # Fit the linear function to log log data
        popt, pcov = curve_fit(func, np.log(x_init), np.log(y_init))
        if type == 'Surface Density' and name == 'Gas':
            print(popt)

        # Retransformed function and parameters
        def func_exp(x,a,b):
            B = np.exp(b)
            return B * x**a

        if plot == True:

            def raymond(r):
                y = 4000 * au / (r* R_S)
                return y / M_sun * R_sun2

            fig, ax = plt.subplots()
            unit = ' in M_sun / R_sun^2'
            ax.set_xlabel('radius in Solar Radii')
            ax.set_ylabel(type + unit)
            ax.set_title(type + ' profile')
            #ax.loglog(r_resh[0,:,0]/au, sigma1d)
            ax.plot(x_init, y_init, 'ko', label="Original Data")
            ax.plot(x, func_exp(x,*popt), 'r-', label="Disk Model")
            ax.plot(np.concatenate((x,r_all)), raymond(np.concatenate((x,r_all))), label="Raymond")
            ax.plot(r_all,func_exp(r_all,*popt), label="Fitted Curve")
            plt.legend()
            plt.savefig('Disk_Profiles/'  + type + '_' + space + '_' + str(N) + '.png')

        return np.exp(func(np.log(x), *popt))

    sigma_extra = extrapolate('Surface Density', rs, sigma1d)
    opac_extra = extrapolate('Opacity', rs, opac1d)
    temp_extra = extrapolate('Temperature', rs, temp1d)

    return x, sigma_extra, opac_extra, temp_extra

def write_disk(name,space,df,dim=1):

    if dim == 1:
        x,sigma,_,temp = integrate_1D(name,space, df, True)
    elif dim == 2:
        x,sigma,_,temp = integrate_2D(df)
    else:
        raise Exception('Not Possible Intergration')


    out = {}
    out['#r [R_s]'] = x
    out['sigma [M_S/R_S^2]'] = sigma
    out['T [K]'] = temp

    df_out = pd.DataFrame.from_dict(out)
    if space == 'log':
        outpath = 'Disk_Profiles/disk_log_' + str(N) + '.txt'
    else:
        outpath = 'Disk_Profiles/disk_' + str(N) + '.txt'

    df_out.to_csv(outpath, index=None, sep=' ', mode='w')

    return df_out



dust = import_ppd('ppd/lev0_dust.dat')
gas = import_ppd('ppd/lev0_gas.dat')

#a = integrate_1D(gas,True)
#b = integrate_1D(dust,True)
space = "lin"

dust_out = write_disk('Dust', space, import_ppd('ppd/lev0_dust.dat'))
gas_out = write_disk('Gas', space, import_ppd('ppd/lev0_gas.dat'))