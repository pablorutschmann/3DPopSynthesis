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

# Gravitational Constant
G = astroconst.G.decompose(u.cgs.bases).value

# Seconds in a year
year = 31536000

# Volume Density
rho = 5.513 / M_S * R_S ** 3

# Solar Radius to AU conversion Factor
rtoau = R_S / au

# Solar Mass to Earth Mass Conversion Factor
mtome = M_S / M_E

# Solar Mass to Jupiter Mass Conversion Factor
mtomj = M_S / M_J
