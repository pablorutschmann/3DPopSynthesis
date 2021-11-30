import numpy as np
import astropy.constants as astroconst
from astropy import units as u
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.integrate import tplquad
from math import log
import os.path

from scipy.optimize import fsolve

# Defining Constants
# AU in cm
au = astroconst.au.decompose(u.cgs.bases).value

# Jupiter Radius in cm
#R_J = astroconst.R_jup.decompose(u.cgs.bases).value
R_J = astroconst.R_jup.decompose(u.cgs.bases).value

# Jupiter Radius squared
R_J2 = R_J ** 2

# Jupiter Mass in grams
#M_J = astroconst.M_jup.decompose(u.cgs.bases).value
M_J = astroconst.M_jup.decompose(u.cgs.bases).value

# Solar Radius in grams
R_S = astroconst.R_sun.decompose(u.cgs.bases).value

# Solar Radius squared
R_S2 = R_S ** 2

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

rho = 5.513 / M_S * R_S**3

G_S = (G / R_S**3) * M_S * year**2


if __name__ == "__main__":
    #Effective Temperature of Sun
    T_S = 5780

    T_J = 130

    SMA_J = 5.2
    #
    print('alpha')
    alpha = T_J/T_S * (SMA_J/(R_S/au))**0.5
    print(alpha)

    #inputs in au
    def temp(r):
        return (T_S * (R_S/au / r)**0.5 * alpha) - 170

    def temp(r):
        return (5780 * (1 / r)**0.5 * 0.7520883995742579)

    def temp_ronco(r):
        y = 280 * (au / r)**(0.5) - 170
        return y

    r_ice = fsolve(temp_ronco,3)
    print("ICE LIne Radius")
    print(r_ice/au)

    temps = [temp(x) for x in np.linspace(0.5,30,100)]

    print(temps)



    m1,m2 = 2.80226e-07, 2.80226e-07
    wmf1, wmf2 = 0.1, 0.1
    wm1 = m1 * wmf1
    wm2 = m2 * wmf2

    wmf = (wm1 + wm2)/(m1+m2)
    print("WMF: ")
    print(wmf)

    radius = 100 * 1000 * 100 #cm

    vol = 4/3 * np.pi * radius**3

    mass = 5.513 * vol

    print(mass / M_S)

    print(0.01 * M_E /M_S)

    print(5.1e-07 * M_S / M_E)

    # radi = (mass * 3 / 4 / np.pi / 2)**(1/3) / 100 / 1000
    #
    # print(radi)
    #
    # mass = 2 * vol / M_E
    #
    # print(mass)

    total_plaenetsima_mass = 3.0e-08 * 100;
    print("TPM")
    print(total_plaenetsima_mass * M_S / M_E)

    # 100 ME / Myrs
    # = 100 ME / (1000000 Yrs)
    # = (100 ME / 1000000) / Yrs
    print(100 * M_E / M_S / 1000000)

    print(year * 9.47834e-07)

    kb = 1.380649e-23 * 1000 * 100 * 100 / M_S / R_S2 * 31536000

    kb = 1.380649e-23 * 1000 * 100 * 100
    print(kb)

    mh20 = 18.02

    Pascal = 1000 / 100 / M_S * R_S * 31536000**2
    print(Pascal)

    wm = 1.0e-8 * 0.2

    T = 203.15

    def P_Vap(T):
        A = -2663.5
        B = 12.537
        P_PA = np.power(10, A/T + B)
        return P_PA

    print(P_Vap(T))

    def Vap_Rate(T):
        return np.sqrt(mh20/ ( 2 * np.pi * kb * T)) * P_Vap(T)/10


    print(Vap_Rate(T))

    a1 = 0.12 / M_S * R_S2
    print(a1)
    def dmdt(T):
        # g / cm**2 / orbital period
        return a1/np.sqrt(T) * np.exp(-1865/T)


    print(dmdt(T)*10e6)
    # print(1.1061960808000799E-22 * 10e6)
    #
    #
    # print(300 / 1000000)
    # print(1.4849409446705092E-17 * R_S / year)
    # print(15 * 100 / R_S * year)
    #
    # def kepl_velo(r):
    #     # convert r to cgs
    #     v_kep = np.sqrt(G * M_S / r ** 3)
    #     # convert back to Solar Units
    #     return v_kep
    #
    # print(kepl_velo(au))
    #
    # print(300 * 0.0001)

    a1 = 0.12 / M_S * R_S2

    print(G * year**2)
    print(float(np.format_float_positional((G / R_J**3) * M_J * year**2, precision=4, unique=False, fractional=False,trim='k')))
    rj = 69911 * 1000 * 100
    print(float(np.format_float_positional((G / rj**3) * M_J * year**2, precision=4, unique=False, fractional=False,trim='k')))


    def gauss(n=11,sigma=1):
        r = range(-int(n/2),int(n/2)+1)
        return [1 / (sigma * np.sqrt(2*np.pi)) * np.exp(-float(x)**2/(2*sigma**2)) for x in r]



    def sf(x,prec=3):
        x = float(np.format_float_positional(x, precision=prec, unique=False, fractional=False,trim='k'))
        return x
    vec_sf = np.vectorize(sf)

    print(list(vec_sf(gauss(30,5))))
    print(np.sum(vec_sf(gauss(30,5))))

class disk:
    def __init__(self, path, pick = False):
        #self.wmf = np.vectorize(self.wmf_1d)
        #self.swmf = np.vectorize(self.swmf_1d)
        self.WMF = np.vectorize(self.wmf)

        # Parameters
        self.inpath = path
        self.spacing = 'None'
        self.N_front = -1
        self.N_back = -1
        self.R_min = -1 # Minimum Radius in AU for extrapolation
        self.R_max = -1 # Maximum Radius in AU for extrapolation
        self.N = -1 # Number of cells in extrapolation axis
        self.from_pickle = pick

        # Shape of Data, [phi_cells, r_cells, th_cells]
        self.phi_cells = 680
        self.r_cells = 215
        self.th_cells = 20
        self.dth = 0.006461
        self.dphi = self.vec_sf(2 * np.pi / 680,prec=6)
        self.dr = 7.18092128e+11

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
        if os.path.isfile(self.inpath + path + '.pickle') and self.from_pickle == True:
            disk = pd.read_pickle(self.inpath + path + '.pickle')
        elif self.from_pickle == False or not(os.path.isfile(self.inpath + path + '.pickle')):
            # Initial Columns
            columns = ['x [cm]', 'y [cm]', 'z [cm]', 'dens [g/cm3]', 'temp [K]', 'velo_x [cm/s]', 'velo_y [cm/s]',
                       'velo_z [cm/s]', 'opac[cm^2/g]']

            # Read the file
            disk = pd.read_csv(self.inpath + path, delim_whitespace=True, names=columns, index_col=False, dtype=np.float64)

            # Convert to cylindrcal coordinates
            #radius, phi = self.cart2pol(disk['x [cm]'], disk['y [cm]'])
            radius, phi, theta = self.cart2sph(disk['x [cm]'], disk['y [cm]'], disk['z [cm]'])
            disk['r [AU]'] = self.vec_sf(radius)
            disk['phi [rad]'] = self.vec_sf(self.vec_sf(phi) + self.vec_sf(np.pi),prec=8)
            disk['theta [rad]'] = self.vec_sf(theta)

            #radius, phi, theta = self.vec_sf(self.cart2sph(self.vec_sf(disk['x [cm]']), self.vec_sf(disk['y [cm]']), self.vec_sf(disk['z [cm]'])),prec=8)
            # disk['r [AU]'] = radius
            # disk['phi [rad]'] = phi
            # disk['theta [rad]'] = theta
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
        # Radii of data values
        rs = np.reshape(df['r [AU]'].values, [self.phi_cells, self.r_cells, self.th_cells])
        phis = np.reshape(df['phi [rad]'].values, [self.phi_cells, self.r_cells, self.th_cells])
        thetas = np.reshape(df['theta [rad]'].values, [self.phi_cells, self.r_cells, self.th_cells])

        def check_unique(arr, type):
            print('Unique ' + type )
            flag = False
            unique = np.unique(arr)
            print(len(unique))

        check_unique(rs, 'Radius')
        ass = []
        unique = np.unique(rs)
        for i in range(self.r_cells-1):
            ass.append(unique[i+1]/unique[i])

        a = np.mean(ass)
        check_unique(phis, 'Phi')
        check_unique(thetas, 'Theta')

        unique = unique.tolist()
        #logarithmically extend by 1 on both sides for distance measure
        unique.append(unique[-1]*a)
        unique.insert(0,unique[0]/a)
        # Calculate Delta Radius
        dx = [(a-1) * x for x in unique]

        # density
        dens = np.reshape(df['dens [g/cm3]'].values, [self.phi_cells, self.r_cells, self.th_cells])
        # Integrate in z direction

        def cell_height(arr):
            dx = []
            for i in range(len(arr)-1):
                dx.append(arr[i+1]-arr[i])
            heights = []
            heights.append(dx[0])
            for i in range(0,len(dx)-1):
                heights.append(dx[i]/2 + dx[i+1]/2)
            heights.append(dx[-1])
            return np.array((heights))

        #calulate cell Volume
        Vols = []
        # for i in range(len(rs[0,:,0])):
        #     for th in thetas[0,0,:]:
        #         r_tuple = (dx[i],dx[i+1])
        #         Vols.append(self.cell_vol(rs[0,i,0],r_tuple,self.dphi,th,self.dth))
        #
        # Vols = np.asarray(Vols)
        # Vols = np.reshape(Vols, (self.th_cells,self.r_cells))
        # np.save('ppd/Vols', Vols)
        Vols = np.load('/Users/prut/CLionProjects/3DPopSynthesis/ppd/Vols.npy')

        #integrate vertically
        dens_mean = np.sum(dens, axis=(0))
        annulus_masses = np.einsum('jk,kj->j',dens_mean, Vols)

        #sigma1d = np.sum(sigma2d, axis=(0))
        #sigma2d = np.sum(dens, axis=2) * 202799192448.87842
        # Average azimutal direction

        ann = []
        for r in rs[0,:,0]:
            ann.append(np.pi * ((r + self.dr/2)**2 - (r - self.dr/2)**2))

        #Surface Density
        sigma1d = annulus_masses/ann

        #Mirror Surface Density
        sigma1d *= 2

        total_mass = np.dot(sigma1d,ann) / M_S
        print('Total mass: ' + str(total_mass))

        # Convert to CODEUNTIS
        sigma1d = sigma1d / M_S * R_S2

        # Truncate front and back artifacts
        sigma1d = truncate(sigma1d,self.N_front, self.N_back)

        self.save_profile(phase, False, 'Surface Density [M_S / R_S^2]', sigma1d)

        # opacity
        opac = np.reshape(df['opac[cm^2/g]'].values, [self.phi_cells, self.r_cells, self.th_cells])
        # average in vertical and azimutal directions
        opac1d = truncate(np.mean(opac, axis=(0, 2)), self.N_front, self.N_back)
        # convert to Solar Radius and Mass units
        opac1d = opac1d * M_S / R_S2

        self.save_profile(phase, False, 'Opacity [M_S / R_S^2]', opac1d)

        # temperature
        temp = np.reshape(df['temp [K]'].values, [self.phi_cells, self.r_cells, self.th_cells])
        # average in vertical and azimutal directions
        temp1d = truncate(np.mean(temp, axis=(0, 2)), self.N_front, self.N_back)

        self.save_profile(phase, False, 'Temperature [K]', temp1d)

        rs = truncate(rs[0, :, 0], self.N_front, self.N_back) / R_S

        self.save_profile(phase, False, 'Radius [R_S]', rs)


    def extrapolate(self, spacing, R_min, R_max, N):
        self.spacing = spacing
        self.R_min = R_min
        self.R_max = R_max
        self.N = N

        # Extrapolation linspace or logspace of radii in R_S, include one point more for N cells
        if self.spacing == 'lin':
            x = np.linspace(self.R_min, self.R_max, self.N + 1) * au / R_S
        elif self.spacing == 'log':
            a = np.power(self.R_max/self.R_min, 1/self.N)
            x = np.geomspace(self.R_min,self.R_max,self.N) * au / R_S
            dx = [(a-1) * item for item in x]
            dx.append((a-1)*x[-1]*a)
            dx.insert(0,(a-1)*x/a)



        r_all = np.linspace(R_max, self.gas_profiles['Radius [R_S]'].max() * R_S / au, N) * au / R_S

        # calculate delta_x
        # dx = np.zeros(N)
        # for i in range(N):
        #     dx[i] = x[i + 1] - x[i]
        #
        # # remove last element
        # x = x[:-1]

        self.out['r [R_S]'] = x
        dx.pop(0)
        dx.pop(-1)
        self.out['dr [R_S]'] = np.array(dx)

        # Get Water Mass Fractions
        self.out['WMF []'],self.out['SWMF []'] = self.WMF(x)

        # Calculate Keplerian Veloctiy
        def kepl_velo(r):
            # convert r to cgs
            r = r * R_S
            v_kep = np.sqrt(G * M_S / r ** 3)
            # convert back to Solar Units
            return v_kep * year

        vx = kepl_velo(x)

        self.out['Keplerian Velocity [R_S / yr]'] = vx

        # Calculate Area of Annulus
        A = np.zeros(N)
        for i in range(N):
            A[i] = np.pi * ((x[i] + dx[i]/2)**2 - (x[i] - dx[i]/2)** 2)

        self.out['Area [R_S^2]'] = A

        for phase in ['gas','dust']:
            if phase == 'gas':
                profiles = self.gas_profiles
            if phase == 'dust':
                profiles = self.dust_profiles

            rs = profiles['Radius [R_S]']

            sigma_extra, popt_sigma = self.get_extrapolate(rs, profiles['Surface Density [M_S / R_S^2]'])
            self.save_profile(phase, True, 'Surface Density [M_S / R_S^2]', sigma_extra)
            self.save_profile(phase, True, 'Power Coefficient Sigma', popt_sigma)

            opac_extra, popt_opac = self.get_extrapolate(rs, profiles['Opacity [M_S / R_S^2]'])
            self.save_profile(phase, True, 'Opacity []', opac_extra)
            self.save_profile(phase, True, 'Power Coefficient Opac', popt_opac)

            temp_extra, popt_temp = self.get_extrapolate(rs, profiles['Temperature [K]'])
            self.save_profile(phase, True, 'Temperature [K]', temp_extra)
            self.save_profile(phase, True, 'Power Coefficient Temp', popt_temp)

        self.out['sigma gas [M_S/R_S^2]'] = self.gas_profiles_extra['Surface Density [M_S / R_S^2]']
        self.out['sigma dust [M_S/R_S^2]'] = self.dust_profiles_extra['Surface Density [M_S / R_S^2]']
        self.out['sigma dustbar [M_S/R_S^2]'] = self.dust_profiles_extra['Surface Density [M_S / R_S^2]']
        self.out['T [K]'] = self.gas_profiles_extra['Temperature [K]']
        self.out['Power Coefficient Density'] = np.full(N, self.gas_profiles_extra['Power Coefficient Sigma'][0])
        self.out['Power Coefficient Temperature'] = np.full(N, self.gas_profiles_extra['Power Coefficient Temp'][0])
        self.out['Gas Opacity []'] = self.gas_profiles_extra['Opacity []']

        print("Inter/Extrapolated Data!")


    def write_disk(self):
        out = self.out.copy()
        # print(len(out))
        # print(type(out))
        # for x,y in out.items():
        #     print(x, type(y))


        # disk profile
        # index 0-(N-1), radius R_S, delta radius R_S, sigma gas, sigma dust M_S R_S2, sigma dust bar, temperature, area of annulus, keplerian orbital velocity omega, sigma exponent, temp exponent, opacity, wmf , swmf

        # sigma dust for now with gas to dust ratio of 0.01
        sigma = self.raymond(self.out['r [R_S]'])
        sigma_dust = sigma * 0.01

        # out['sigma gas [M_S/R_S^2]'] = sigma
        # out['sigma dust [M_S/R_S^2]'] = sigma_dust
        # out['sigma dustbar [M_S/R_S^2]'] = sigma_dust
        out['Power Coefficient Density'] = np.full(self.N, -1)
        #out['Power Coefficient Temperature'] = np.full(self.N, -0.5)
        out['Gas Opacity []'] = np.full(self.N, 0)

        #Headwind Factor
        out["Headwind Factor []"] = self.eta_peb(self.out['r [R_S]'])


        for key in ['r [R_S]', 'dr [R_S]', 'sigma gas [M_S/R_S^2]', 'sigma dust [M_S/R_S^2]', 'sigma dustbar [M_S/R_S^2]', 'T [K]', 'Area [R_S^2]', 'Keplerian Velocity [R_S / yr]', 'Power Coefficient Density', 'Power Coefficient Temperature', 'Gas Opacity []', 'WMF []', 'SWMF []', 'Headwind Factor []']:
            #self.out = self.out.pop[key]
            out[key] = out.pop(key)

        df_out = pd.DataFrame.from_dict(out)
        df_out.index.name = 'index'
        if self.spacing == 'log':
            outpath = 'disks/disk_log_' + str(self.N) + '.txt'
        else:
            outpath = 'disks/disk_' + str(self.N) + '.txt'

        self.tot_mass_gas = np.dot(out['sigma gas [M_S/R_S^2]'],out['Area [R_S^2]'])
        self.tot_mass_dust = np.dot(out['sigma dust [M_S/R_S^2]'],out['Area [R_S^2]'])
        self.tot_mass = self.tot_mass_gas + self.tot_mass_dust

        print('Total Mass of Gas: ' + str(self.tot_mass_gas*M_S/M_J))
        print('Total Mass of Dust: ' + str(self.tot_mass_dust*M_S/M_J))
        print('Total Mass: ' + str(self.tot_mass*M_S/M_J))
        print(df_out.columns)
        df_out.to_csv(outpath, header=False, sep='\t', mode='w')

        print("Written disk file!")


    # Helper Functions
    def cart2pol(self, x, y):
        radius = np.sqrt(x ** 2 + y ** 2)
        phi = np.arctan2(y, x)
        return (radius, phi)

    def cart2sph(self, x, y, z):
        radius = np.sqrt(x ** 2 + y ** 2 + z ** 2)
        theta = np.arctan(z / np.sqrt(x ** 2 + y ** 2))
        #phi = np.arctan2(self.vec_sf(y, prec = 5),self.vec_sf(x, prec = 5))
        phi = np.arctan2(y,x)
        for i in range(len(x)):
            if x[i] <= 0:
                #theta[i] -= 2 * theta[i]
                a=1
        # for i in range(len(y)):
        #         phi[i] += np.pi

        return (radius, phi , theta)

    def cell_vol(self, r, dr , dphi, theta, dtheta):
        def theta_trafo(th):
            return np.pi/2 - th
        theta = theta_trafo(theta)
        theta1 = theta - dtheta/2
        theta2 = theta + dtheta/2
        dr1, dr2 = dr
        r1 = r - dr1/2
        r2 = r + dr2/2

        phi1 = 0
        phi2 = dphi

        def integrant(phi,th,r):
            return r**2 * np.sin(th)
        V, dV = tplquad(integrant, r1, r2, lambda x: theta1, lambda x: theta2,lambda x, y: phi1, lambda x, y: phi2)
        return V

    # Significant Digits Rounding Function
    def sf(x,prec=13):
        x = float(np.format_float_positional(x, precision=prec, unique=False, fractional=False,trim='k'))
        return x
    vec_sf = np.vectorize(sf)

    def func(self, x, a, b):
        return a * x + b

    # Retransformed function and parameters
    def func_exp(self,x, a, b):
        B = np.exp(b)
        return B * x ** a

    def raymond(self, r):
        y = 4000 * au / (r * R_S)
        return y / M_S * R_S2

    def binkert(self, r):
        y = 80 * (au / (r * R_S))**0.5
        return y / M_S * R_S2

    def temp_ronco(self, r):
        y = 280 * (au / (r * R_S))**(0.5)
        return y

    def eta_peb(self,r):
        r_au = r / au * R_S
        y = 1.5e-3 * r**(0.5)
        return y

    def get_extrapolate(self, x_init, y_init):
        # Fit the linear function to log log data
        popt, pcov = curve_fit(self.func, np.log(x_init), np.log(y_init))

        return self.func_exp(self.out['r [R_S]'], *popt), popt

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

    def mirror_theta(self, arr):
        flipped = np.flip(arr, axis=0)
        mirrored = np.concatenate((arr, flipped), axis=0)
        return mirrored

    def wmf_1d(self,r):
        r_au = r / au * R_S
        if r_au < 2:
            wmf = 0.0
        elif 2 <= r_au < 2.5:
            wmf = 0.2
        elif 2.5 <= r_au:
            wmf = 0.4

        return wmf

    def swmf_1d(self,r):
        r_au = r / au * R_S
        # if r_au < 2:
        #     wmf = 0.0
        # elif 2 <= r_au < 2.5:
        #     wmf = 0.1
        # elif 2.5 <= r_au:
        #     wmf = 0.5
        swmf = 0.1

        return swmf

    def wmf(self, r):
        # From Morbidelli et al, 2012
        IceLine = 2.71280277
        r_au = r / au * R_S

        if r_au < IceLine:
            wmf = 0.0
            if r_au < 1.5:
                swmf = 0.0
            elif 1.5 <= r_au < 2.0:
                swmf = 0.01
            else:
                swmf = 0.1

        elif IceLine <= r_au < IceLine + 0.3 :
            wmf = 0.2
            swmf = 0.2

        elif IceLine + 0.3 <= r_au:
            wmf = 0.4
            swmf = 0.2

        return 0.2, 0.2

    WMF = np.vectorize(wmf)










