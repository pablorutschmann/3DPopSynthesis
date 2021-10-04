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

class disk:
    def __init__(self, path):
        # Parameters
        self.inpath = path
        self.spacing = 'None'
        self.N_front = -1
        self.N_back = -1
        self.R_min = -1 # Minimum Radius in AU for extrapolation
        self.R_max = -1 # Maximum Radius in AU for extrapolation
        self.N = -1 # Number of cells in extrapolation axis

        # Shape of Data, [phi_cells, r_cells, th_cells]
        self.phi_cells = 680
        self.r_cells = 215
        self.th_cells = 20

        self.tot_mass_gas = -1
        self.tot_mass_dust = -1
        self.tot_mass = -1

        self.gas_profiles = {}
        self.gas_profiles_extra = {}
        self.dust_profiles = {}
        self.dust_profiles_extra = {}
        self.out = {}

        #import data
        self.gas = self.import_ppd('lev0_gas.dat')
        self.dust = self.import_ppd('lev0_dust.dat')

    # Import Function, Saves Disk to Dataframe in cgs units
    def import_ppd(self, path):
        # Initial Columns
        columns = ['x [cm]', 'y [cm]', 'z [cm]', 'dens [g/cm3]', 'temp [K]', 'velo_x [cm/s]', 'velo_y [cm/s]',
                   'velo_z [cm/s]', 'opac[cm^2/g]']

        # Read the file
        disk = pd.read_csv(self.inpath + path, delim_whitespace=True, names=columns, index_col=False)

        # Convert to cylindrcal coordinates
        radius, phi = self.cart2pol(disk['x [cm]'], disk['y [cm]'])
        disk['r [AU]'] = radius
        disk['phi [rad]'] = phi + np.pi
        print('Imported disk from: ' + path)
        return disk

    def integrate_1D(self, phase):

        if phase == 'gas':
            df = self.gas
        if phase == 'dust':
            df = self.dust

        def truncate(arr, N_front,N_back):
            trunc = arr[N_front:-N_back]
            return trunc

        # density
        dens = np.reshape(df['dens [g/cm3]'].values, [self.phi_cells, self.r_cells, self.th_cells])
        # height one cell 202799192448.87842
        # Integrate in z direction
        sigma2d = np.sum(dens, axis=2) * 203599192448.87842
        # Average azimutal direction
        sigma1d = truncate(np.sum(sigma2d, axis=(0)), self.N_front, self.N_back)
        # convert to Solar Radius and Mass units
        sigma1d = sigma1d / M_J * R_J2

        self.save_profile(phase, False, 'Surface Density [M_J / R_J^2]', sigma1d)

        # opacity
        opac = np.reshape(df['opac[cm^2/g]'].values, [self.phi_cells, self.r_cells, self.th_cells])
        # average in vertical and azimutal directions
        opac1d = truncate(np.mean(opac, axis=(0, 2)), self.N_front, self.N_back)
        # convert to Solar Radius and Mass units
        opac1d = opac1d * M_J / R_J2

        self.save_profile(phase, False, 'Opacity [M_J / R_J^2]', opac1d)

        # temperature
        temp = np.reshape(df['temp [K]'].values, [self.phi_cells, self.r_cells, self.th_cells])
        # average in vertical and azimutal directions
        temp1d = truncate(np.mean(temp, axis=(0, 2)), self.N_front, self.N_back)

        self.save_profile(phase, False, 'Temperature [K]', temp1d)

        # Radii of data values
        rs = np.reshape(df['r [AU]'].values, [self.phi_cells, self.r_cells, self.th_cells])
        rs = truncate(rs[0, :, 0], self.N_front, self.N_back) / R_J

        self.save_profile(phase, False, 'Radius [R_J]', rs)


    def extrapolate(self, spacing, R_min, R_max, N):
        self.spacing = spacing
        self.R_min = R_min
        self.R_max = R_max
        self.N = N

        # Extrapolation linspace or logspace of radii in R_J, include one point more for N cells
        if self.spacing == 'lin':
            x = np.linspace(self.R_min, self.R_max, self.N + 1) * au / R_J
        elif self.spacing == 'log':
            x = np.logspace(1, log(self.R_max, self.R_min), self.N + 1, base=self.R_min) * au / R_J


        r_all = np.linspace(R_max, self.gas_profiles['Radius [R_J]'].max() * R_J / au, N) * au / R_J

        # calculate delta_x
        dx = np.zeros(N)
        for i in range(N):
            dx[i] = x[i + 1] - x[i]

        # remove last element
        x = x[:-1]

        self.out['r [R_J]'] = x
        self.out['dr [R_J]'] = dx

        # Calculate Keplerian Veloctiy
        def kepl_velo(r):
            # convert r to cgs
            r = r * R_J
            v_kep = np.sqrt(G * M_S / r ** 3)
            # convert back to Jupiter Units

            return v_kep / year

        vx = kepl_velo(x)

        self.out['Keplerian Velocity [R_J / yr]'] = vx

        # Calculate Area of Annulus
        A = np.zeros(N)
        for i in range(N):
            A[i] = np.pi * ((x[i] + dx[i]) ** 2 - x[i] ** 2)

        self.out['Area [R_J^2]'] = A

        for phase in ['gas','dust']:
            if phase == 'gas':
                profiles = self.gas_profiles
            if phase == 'dust':
                profiles = self.dust_profiles

            rs = profiles['Radius [R_J]']

            sigma_extra, pow_coeff_sigma = self.get_extrapolate(rs, profiles['Surface Density [M_J / R_J^2]'])
            self.save_profile(phase, True, 'Surface Density [M_J / R_J^2]', sigma_extra)
            self.save_profile(phase, True, 'Power Coefficient Sigma', pow_coeff_sigma)

            opac_extra, pow_coeff_opac = self.get_extrapolate(rs, profiles['Opacity [M_J / R_J^2]'])
            self.save_profile(phase, True, 'Opacity []', opac_extra)
            self.save_profile(phase, True, 'Power Coefficient Opac', pow_coeff_opac)

            temp_extra, pow_coeff_temp = self.get_extrapolate(rs, profiles['Temperature [K]'])
            self.save_profile(phase, True, 'Temperature [K]', temp_extra)
            self.save_profile(phase, True, 'Power Coefficient Temp', pow_coeff_temp)

        self.out['sigma gas [M_J/R_J^2]'] = self.gas_profiles_extra['Surface Density [M_J / R_J^2]']
        self.out['sigma dust [M_J/R_J^2]'] = self.dust_profiles_extra['Surface Density [M_J / R_J^2]']
        self.out['sigma dustbar [M_J/R_J^2]'] = self.dust_profiles_extra['Surface Density [M_J / R_J^2]']
        self.out['T [K]'] = self.gas_profiles_extra['Temperature [K]']
        self.out['Power Coefficient Density'] = np.full(N, self.gas_profiles_extra['Power Coefficient Sigma'])
        self.out['Power Coefficient Temperature'] = np.full(N, self.gas_profiles_extra['Power Coefficient Temp'])
        self.out['Gas Opacity []'] = self.gas_profiles_extra['Opacity []']

        print("Inter/Extrapolated Data!")


    def write_disk(self):
        out = self.out.copy()
        print(len(out))
        print(type(out))
        for x,y in out.items():
            print(x, y.shape)
        print(out.values)
        # disk profile
        # index 0-(N-1), radius R_J, delta radius R_J, sigma gas, sigma dust M_J R_J2, sigma dust bar, temperature, area of annulus, keplerian orbital velocity omega, sigma exponent, temp exponent, opacity

        # sigma dust for now with gas to dust ratio of 0.01
        sigma = self.raymond(self.out['r [R_J]'])
        sigma_dust = sigma * 0.01

        out['sigma gas [M_J/R_J^2]'] = sigma
        out['sigma dust [M_J/R_J^2]'] = sigma_dust
        out['sigma dustbar [M_J/R_J^2]'] = sigma_dust
        out['Power Coefficient Density'] = np.full(self.N, 0)
        out['Power Coefficient Temperature'] = np.full(self.N, 0)
        out['Gas Opacity []'] = np.full(self.N, 0)


        for key in ['r [R_J]', 'dr [R_J]', 'sigma gas [M_J/R_J^2]', 'sigma dust [M_J/R_J^2]', 'sigma dustbar [M_J/R_J^2]', 'T [K]', 'Area [R_J^2]', 'Keplerian Velocity [R_J / yr]', 'Power Coefficient Density', 'Power Coefficient Temperature', 'Gas Opacity []']:
            #self.out = self.out.pop[key]
            out[key] = out.pop(key)

        df_out = pd.DataFrame.from_dict(out)
        if self.spacing == 'log':
            outpath = 'disks/disk_log_' + str(self.N) + '.txt'
        else:
            outpath = 'disks/disk_' + str(self.N) + '.txt'

        self.tot_mass_gas = np.dot(out['sigma gas [M_J/R_J^2]'],out['Area [R_J^2]'])
        self.tot_mass_dust = np.dot(out['sigma dust [M_J/R_J^2]'],out['Area [R_J^2]'])
        self.tot_mass = self.tot_mass_gas + self.tot_mass_dust

        print('Total Mass of Gas: ' + str(self.tot_mass_gas))
        print('Total Mass of Dust: ' + str(self.tot_mass_dust))
        print('Total Mass: ' + str(self.tot_mass))

        df_out.to_csv(outpath, header=True, sep='\t', mode='w')

        print("Written disk file!")


    # Helper Functions
    def cart2pol(self, x, y):
        radius = np.sqrt(x ** 2 + y ** 2)
        phi = np.arctan2(y, x)
        return (radius, phi)

    def func(self, x, a, b):
        return a * x + b

    # Retransformed function and parameters
    def func_exp(self,x, a, b):
        B = np.exp(b)
        return B * x ** a

    def raymond(self, r):
        y = 4000 * au / (r * R_J)
        return y / M_J * R_J2

    def temp_ronco(self, r):
        y = 280 * (au / (r * R_J))**(0.5)
        return y

    def get_extrapolate(self, x_init, y_init):
        # Fit the linear function to log log data
        popt, pcov = curve_fit(self.func, np.log(x_init), np.log(y_init))

        pow_coeff = popt[0]

        return self.func_exp(self.out['r [R_J]'], *popt), pow_coeff

    def save_profile(self, type, extra, key, arr):
        if extra:
            if type == 'gas':
                self.gas_profiles_extra[key] = arr
            if type == 'dust':
                self.dust_profiles_extra[key] = arr
        if not extra:
            if type == 'gas':
                self.gas_profiles[key] = arr
            if type == 'dust':
                self.dust_profiles[key] = arr

