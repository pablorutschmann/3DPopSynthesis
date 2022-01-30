import astropy.constants as astroconst
from astropy import units as u
import numpy as np

# Defining Constants
# AU in cm
au = astroconst.au.decompose(u.cgs.bases).value

# Jupiter Radius in cm
# R_J = astroconst.R_jup.decompose(u.cgs.bases).value
R_J = astroconst.R_jup.decompose(u.cgs.bases).value

# Jupiter Radius squared
R_J2 = R_J ** 2

# Jupiter Mass in grams
# M_J = astroconst.M_jup.decompose(u.cgs.bases).value
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

# Mars Mass
M_M = 0.107 * M_E

# Gravitational Constant
G = astroconst.G.decompose(u.cgs.bases).value

# Seconds in a year
year = 31536000

# Volume Density
rho = 5.513 / M_S * R_S ** 3

# AU to Solar Radius Conversion Factor
autor = au / R_S

# Solar Radius to AU Conversion Factor
rtoau = R_S / au

# Solar Mass to Earth Mass Conversion Factor
mtome = M_S / M_E

# Solar Mass to Jupiter Mass Conversion Factor
mtomj = M_S / M_J

# Densities in cgs to Solar Units
denstos = R_S2 / M_S

# Earth Ocean in Earth Masses
OE = 2.4E-4 * M_E

# Solar Mass to Earths Ocean
mtooe = M_S / OE

print(7.63764e-07 * M_S / M_E)

print(65.0 * R_S / au)
print(6540 * R_S / au)

M_ceres = 9.1e20 * 1000

initm = 3.0035e-08
initm2 = 3.0035e-11
print("Ceres")
print(initm * M_S / M_ceres)
print(initm2 * M_S / M_ceres)

print(initm * M_S / M_E)
print(initm2 * M_S / M_E)

r = pow(3 * initm * M_S / 4 / np.pi / 5.513,1/3)
print(r/R_E)
print(r / 100 / 1000)

r = pow(3 * initm2 * M_S / 4 / np.pi / 5.513,1/3)
print(r/R_E)
print(r / 100 / 1000)
