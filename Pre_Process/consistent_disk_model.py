import numpy as np
import astropy.constants as astroconst
from astropy import units as u
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

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

rho = 5.513 / M_S * R_S ** 3

G_S = (G / R_S ** 3) * M_S * year ** 2

R_min = 0.5
R_max = 30
N = 1000

# getting the radius and diffrerences

a = np.power(R_max / R_min, 1 / N)
x = np.geomspace(R_min, R_max, N)
dx = [(a - 1) * item for item in x]
dx.append((a - 1) * x[-1] * a)
dx.insert(0, (a - 1) * x / a)

dx.pop(0)
dx.pop(-1)

x = np.array(x)
dx = np.array(dx)

# Area

A = np.zeros(N)
for i in range(N):
    A[i] = np.pi * ((x[i] + dx[i] / 2) ** 2 - (x[i] - dx[i] / 2) ** 2)
A = np.array(A)


# Density

def raymond(r):
    # coeff = -1.0
    coeff = np.random.uniform(-1.5, -0.5)
    y = 4000 * (r) ** coeff
    return y


def binkert(r):
    coeff = -0.5
    y = 800 * (r) ** coeff
    return y

def hayashi(r):
    coeff = -1.5
    y = 1700 * (r) ** coeff
    return y


# fig, ax = plt.subplots()
# ax.set_xlabel('Radius in AU')
# ax.plot(x, binkert(x), label="Binkert")
# ax.plot(x, raymond(x), label="Raymond")
# ax.plot(x, hayashi(x), label="Hayashi")
# plt.legend()
# ax.set_xscale("log", base=10)
# ax.set_yscale("log", base=10)
#plt.show()


gas_density = binkert(x)

dust_density = gas_density * 0.01

x *= au
dx *= au
A *= au ** 2


# planetesimals

def kepl_velo(r):
    # convert r to cgs
    r *= au
    v_kep = np.sqrt(G * M_S / r)
    # convert back to Solar Units
    return v_kep


# pebbles
PebbleFlux = 300 * M_E / (1e6 * year)
stokes_number = 0.1


def Eta(r):
    y = 1.5e-3 * r ** (0.5)
    return y


def v_r(r):
    return 2 * (stokes_number / (stokes_number ** 2 + 1)) * Eta(r) * kepl_velo(r)

def Pebble_Density(r):
    return PebbleFlux / (2 * np.pi * r * v_r(r))

pebble_density = Pebble_Density(x)

print(gas_density[0:-1:len(gas_density) // 20])
print(dust_density[0:-1:len(dust_density) // 20])
print(pebble_density[0:-1:len(pebble_density) // 20])

TM_gas = np.dot(A, gas_density)

TM_dust = np.dot(A, dust_density)

TM_pebble = np.dot(A, pebble_density)

print(TM_gas)
print(TM_dust)
print(TM_pebble)
print("total mass")
TM_disk = TM_pebble + TM_dust + TM_gas

print(TM_disk / M_S)

# Planetesimals

N_planetesimals = 100

rho = 0.93357 * M_S / R_S ** 3

initial_mass = 3.0e-08 * M_S

log_rmin = np.log10(R_min)

log_rmax = np.log10(R_max)

print(10**log_rmin)

radii = 10 ** np.random.uniform(log_rmin, log_rmax, N_planetesimals)
radii = np.sort(radii)
print(radii[::10])

# Planetesimal Mass Distribution
# Radius of Planetesiomals in cm according to Coleman 2021

def Mass(radius):
    return 4 / 3 * np.pi * (radius)**3 * rho

# at log(1) au, radius log(30) km

# at log(10) au, radius log(100) km

dx = np.log10(10) - np.log10(1)
dy = np.log10(90) - np.log10(20)
dxdy = dx/dy
a0 = 100 * 1000 * 100 * np.power(10,np.log10(30) - dx/dy)
a1 = dx/dy
print("a0: ",a0)
print("a1: ",a1)

# initial Planetesimal Radius in cm
def Radius(r):
    # a0:  6548731.842129445
    # a1:  1.6609640474436815
    return  a0 * r**a1


b0 = a0 / R_S
b1 = a1 * au / R_S
print("b0: ",b0)
print("b1: ",b1)

def Radius_S(r_s):
    # b1 = a1 * AU/R_S
    # b0 = a0 /R_S
    return  b0 * r**b1





fig, ax = plt.subplots()
ax.set_xlabel('Distance in AU')
ax.set_ylabel('d')
ax.scatter(radii,Mass(Radius(radii))/M_E, label="Initial Masses of Planetesimals")
ax.set_xscale("log", base=10)
ax.set_yscale("log", base=10)
plt.legend()
#plt.show()
print(np.sum(Mass(Radius(radii)))/M_S)
print(Radius(1)/100/1000)
print(Mass(Radius(30))/M_J)





TM_planetesimals = initial_mass * np.sum(radii)


print(TM_planetesimals / M_S)

# Planetesimal Surface Density

dr = []
for i,item in enumerate(radii[:-1]):
    delta = radii[i+1] - item
    dr.append(delta)

dr.insert(0, radii[0]-R_min)


areas = np.zeros(N_planetesimals)
for i in range(N_planetesimals-1):
    areas[i] = np.pi * ((radii[i+1]) ** 2 - (radii[i]) ** 2)
areas = np.array(areas)

print("Area")
print(np.sum(A))
print(np.sum(areas))

surface_dens = initial_mass * areas

averagesurface = N_planetesimals * initial_mass / np.sum(A)

# ax.hlines(y=averagesurface, xmin=R_min,xmax=R_max, label="Planetesimals")
# plt.show()


def pla_dens(rad, a0, a1, a2):
    coeff = a1
    y = a2 +  a0 * rad**coeff
    return y

popt,_ = curve_fit(pla_dens, radii, surface_dens)

print(popt)










