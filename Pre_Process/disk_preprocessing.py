import numpy as np
import astropy.constants as astroconst
from astropy import units as u
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from math import log

# Defining Constants

# AU in cm
au = astroconst.au.decompose(u.cgs.bases).value

# Jupiter Radius in cm
R_J = astroconst.R_jup.decompose(u.cgs.bases).value

# Jupiter Radius squared
R_J2 = (astroconst.R_jup.decompose(u.cgs.bases).value) ** 2

# Jupiter Mass in grams
M_J = astroconst.M_jup.decompose(u.cgs.bases).value

# Solar Mass in grams
R_S = astroconst.R_sun.decompose(u.cgs.bases).value

# Solar Mass in grams
M_S = astroconst.M_sun.decompose(u.cgs.bases).value

# Earth Masses
M_E = astroconst.M_earth.decompose(u.cgs.bases).value

# Earth Radius
R_E = astroconst.R_earth.decompose(u.cgs.bases).value

# Gravitational Constant
G = astroconst.G.decompose(u.cgs.bases).value

# Seconds in a year
year = 31536000

# print(1.5 * 10 ** (-21) * M_J / R_J2 * 202799192448.87842)

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
space ='log'

print(R_min *  au / R_S)

print(R_max *  au / R_S)

# Helper Functions
def cart2pol(x, y):
    radius = np.sqrt(x ** 2 + y ** 2)
    phi = np.arctan2(y, x)
    return (radius, phi)


# Import Function, Saves Disk to Dataframe in cgs units
def import_ppd(path):
    # Initial Columns
    columns = ['x [cm]', 'y [cm]', 'z [cm]', 'dens [g/cm3]', 'temp [K]', 'velo_x [cm/s]', 'velo_y [cm/s]',
               'velo_z [cm/s]', 'opac[cm^2/g]']

    # Read the file
    disk = pd.read_csv(path, delim_whitespace=True, names=columns, index_col=False, dtype=np.float64)

    # Convert to cylindrcal coordinates
    radius, phi = cart2pol(disk['x [cm]'], disk['y [cm]'])
    disk['r [AU]'] = radius
    disk['phi [rad]'] = phi + np.pi
    print('Imported disk from: ' + path)
    return disk


def integrate_1D(name, space, df, plot=False):
    def truncate(arr, N_front, N_back):
        trunc = arr[N_front:-N_back]
        return trunc

    # density
    dens = np.reshape(df['dens [g/cm3]'].values, [phi_cells, r_cells, th_cells])
    # height one cell 202799192448.87842

    # Integrate in z direction
    sigma2d = np.sum(dens, axis=2) * 203599192448.87842
    # Average azimutal direction
    sigma1d = truncate(np.mean(sigma2d, axis=(0)), N_f, N_b)
    # convert to Solar Radius and Mass units
    sigma1d = sigma1d / M_J * R_J2

    # opacity
    opac = np.reshape(df['opac[cm^2/g]'].values, [phi_cells, r_cells, th_cells])
    # average in vertical and azimutal directions
    opac1d = truncate(np.mean(opac, axis=(0, 2)), N_f, N_b)
    # convert to Solar Radius and Mass units
    opac1d = opac1d * M_J / R_J2

    # temperature
    temp = np.reshape(df['temp [K]'].values, [phi_cells, r_cells, th_cells])
    # average in vertical and azimutal directions
    temp1d = truncate(np.mean(temp, axis=(0, 2)), N_f, N_b)

    # Radii of data values
    rs = np.reshape(df['r [AU]'].values, [phi_cells, r_cells, th_cells])
    rs = truncate(rs[0, :, 0], N_f, N_b) / R_J

    r_all = np.linspace(R_max, rs.max() * R_J / au, N) * au / R_J

    # Extrapolation linspace or logspace of radii in R_J, include one point more for N cells
    if space == "lin":
        x = np.linspace(R_min, R_max, N + 1) * au / R_J
    elif space == 'log':
        x = np.logspace(1, log(R_max, R_min), N + 1, base=R_min) * au / R_J

    # calculate delta_x
    dx = np.zeros(N)
    for i in range(N):
        dx[i] = x[i + 1] - x[i]

    # remove last element
    x = x[:-1]

    # Calculate Keplerian Veloctiy



    def kepl_velo(r):
        # convert r to cgs
        r = r * R_J
        v_kep = np.sqrt(G * M_S / r ** 3)
        # convert back to Jupiter Units

        return v_kep / year

    vx = kepl_velo(x)

    # Calculate Area of Annulus
    A = np.zeros(N)
    for i in range(N):
        A[i] = np.pi * ((x[i] + dx[i]) ** 2 - x[i] ** 2)

    def extrapolate(type, x_init, y_init):

        def func(x, a, b):
            return a * x + b

        # Fit the linear function to log log data
        popt, pcov = curve_fit(func, np.log(x_init), np.log(y_init))

        pow_coeff = popt[0]

        # Retransformed function and parameters
        def func_exp(x, a, b):
            B = np.exp(b)
            return B * x ** a

        if plot == True:
            def raymond(r):
                y = 4000 * (au / (r * R_J))
                return y / M_J * R_J2

            def binkert(r):
                y = 80 * (au / (r * R_J))**0.5
                return y / M_J * R_J2

            def temp_ronco (r):
                y = 280 * (au / (r * R_J))
                return y

            fig, ax = plt.subplots()

            ax.set_xlabel('radius in Jupiter Radii')

            ax.set_title(type + ' profile')
            # ax.loglog(r_resh[0,:,0]/au, sigma1d)
            ax.plot(x_init, y_init, 'ko', label="Original Data")
            ax.plot(x, func_exp(x, *popt), 'r-', label="Disk Model")
            if type == 'Surface Density' and name == 'gas':
                unit = ' in M_J / R_jup^2'
                ax.plot(np.concatenate((x, r_all)), binkert(np.concatenate((x, r_all))), label="Raymond")
            if type == 'Surface Density' and name == 'dust':
                unit = ' in M_J / R_jup^2'
                ax.plot(np.concatenate((x, r_all)), binkert(np.concatenate((x, r_all)))*0.01, label="Raymond")
            if type == 'Temperature':
                unit = ' Kelvin'
                ax.plot(np.concatenate((x, r_all)), temp_ronco(np.concatenate((x, r_all))), label="Ronco")
            else:
                unit = ''
            ax.set_ylabel(type + unit)
            ax.plot(r_all, func_exp(r_all, *popt), label="Fitted Curve")
            plt.legend()
            plt.savefig('Disk_Plots/' + name + '_'  + type + '_' + space + '_' + str(N) + '.png')

        return np.exp(func(np.log(x), *popt)), pow_coeff

    sigma_extra, pow_coeff_sigma = extrapolate('Surface Density', rs, sigma1d)
    opac_extra, pow_coeff_opac = extrapolate('Opacity', rs, opac1d)
    temp_extra, pow_coeff_temp = extrapolate('Temperature', rs, temp1d)

    pow_coeffs = [pow_coeff_sigma, pow_coeff_opac, pow_coeff_temp]

    print("Inter/Extrapolated values for " + name)
    return x, dx, sigma_extra, opac_extra, temp_extra, A, vx, pow_coeffs


def write_disk(space, dim=1):
    gas = import_ppd('../ppd/lev0_gas.dat')
    dust = import_ppd('../ppd/lev0_dust.dat')

    if dim == 1:
        x, dx, sigma, opac, temp, Area, v_kepl, pow_coeffs = integrate_1D('gas', space, gas, True)
        _, _, sigma_dust, _, temp_dust, _, _, pow_coeffs_dust = integrate_1D('dust', space, dust, True)
    elif dim == 2:
        x, dx, sigma, _, temp, Area, v_kepl, pow_coeffs = integrate_2D(dust)
    else:
        raise Exception('Not Possible Intergration')

    # disk profile
    # index 0-(N-1), radius R_J, delta radius R_J, sigma gas, sigma dust M_J R_J2, sigma dust bar, temperature, area of annulus, keplerian orbital velocity omega, sigma exponent, temp exponent, opacity

    # sigma dust for now with gas to dust ratio of 0.01
    sigma_dust = sigma * 0.01

    out = {}
    out['r [R_J]'] = x
    out['dr [R_J]'] = dx
    out['sigma gas [M_J/R_J^2]'] = sigma
    out['sigma dust [M_J/R_J^2]'] = sigma_dust
    out['sigma dustbar [M_J/R_J^2]'] = sigma_dust
    out['T [K]'] = temp
    out['Area [R_J^2]'] = Area
    out['Keplerian Velocity [R_J / yr]'] = v_kepl
    out['Power Coefficient Density'] = np.full(N, pow_coeffs[0])
    out['Power Coefficient Temperature'] = np.full(N, pow_coeffs[-1])
    out['Gas Opacity []'] = opac

    df_out = pd.DataFrame.from_dict(out)
    if space == 'log':
        outpath = 'disks/disk_log_' + str(N) + '.txt'
    else:
        outpath = 'disks/disk_' + str(N) + '.txt'

    df_out.to_csv(outpath, header=False, sep='\t', mode='w')

    print("Written disk file")
    return df_out


# a = integrate_1D(gas,True)
# b = integrate_1D(dust,True)


disk_out = write_disk(space)

print('Solar Mass in Jupiter Masses: ' + str(M_S/M_J))

print('Solar Radius in Jupiter Radii: ' + str(R_S/R_J))

print(0.5 * au / R_J )
print(5 * au / R_J)