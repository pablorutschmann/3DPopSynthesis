import astropy.constants as astroconst
from astropy import units as u

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

