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
