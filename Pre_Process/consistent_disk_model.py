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
N = 100

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


# Density Profiles
def raymond(r):
    coeff = -1.0
    #coeff = np.random.uniform(-1.5, -0.5)
    y = 4000 * r ** coeff
    return y


def binkert(r):
    coeff = -0.5
    y = 100 * r ** coeff
    return y

def hayashi(r):
    coeff = -1.5
    y = 1700 * r ** coeff
    return y

gas_density = binkert(x)

dust_density = gas_density * 0.01

planetesimal_density = gas_density * 0.001


TM_gas = np.dot(A * au ** 2, gas_density)

TM_dust = np.dot(A * au ** 2, dust_density)

print("Total Gas Mass: ", TM_gas / M_J)
print("Total Dust Mass: ", TM_dust / M_J)

TM_disk = TM_dust + TM_gas
print("Total Mass: ", TM_disk / M_J)

# Planetesimals

N_planetesimals = 100

# rho = 0.93357 * M_S / R_S ** 3
rho = 0.93357
# initial_mass = 3.0e-08 * M_S

initial_radius = 100  # km

initial_mass = 4 / 3 * np.pi * (1000 * 100 * initial_radius) ** 3 * rho

print("initial Mass: ", initial_mass / M_E)

log_rmin = np.log10(R_min)
log_rmax = np.log10(R_max)


radii_log = 10 ** np.random.uniform(log_rmin, log_rmax, N_planetesimals)
radii = np.random.uniform(R_min, R_max, N_planetesimals)

radii = np.sort(radii)
radii_log = np.sort(radii_log)

dx = np.log10(10) - np.log10(1)
dy = np.log10(90) - np.log10(20)
dxdy = dx / dy
a0 = 100 * 1000 * 100 * np.power(10, np.log10(20) - dx / dy)
a1 = dx / dy
print("a0: ", a0)
print("a1: ", a1)


# initial Planetesimal Radius in cm
def Radius(r):
    # a0:  6548731.842129445
    # a1:  1.6609640474436815
    return a0 * r ** a1




def get_sd(arr, rdist = False):
    disk_ids = []

    for item in arr:
        index = np.argmin(np.abs(x - item))
        disk_ids.append(index)

    surface_density = np.zeros(N)
    for i in range(N):
        sum = 0
        for j in range(len(arr)):
            if disk_ids[j] == i:
                if not(rdist):
                    sum += initial_mass
                else:
                    sum += 4 / 3 * np.pi * (Radius(x[i])) ** 3 * rho
        surface_density[i] = sum / (A[i] * au ** 2)
    return surface_density

surface_density = get_sd(radii)
surface_density_log = get_sd(radii_log)
surface_density_self = get_sd(radii_self)

def SD(r):
    mean_sd = N_planetesimals * initial_mass / (np.sum(A) * au ** 2)
    return mean_sd


mean_sd = N_planetesimals * initial_mass / (np.sum(A) * au ** 2)

fig, ax = plt.subplots()
ax.set_xlabel('Distance in AU')
ax.set_ylabel('Planetesimal Surface Density in CGS')
ax.set_xlim(0.5,30)
ax.plot(x, binkert(x), label="Binkert")
# ax.plot(x, hayashi(x), label="Hayashi")
# ax.plot(x, raymond(x), label="Raymond")
ax.plot(x, 0.01 * binkert(x), label="Binkert Dust")
# ax.plot(x, np.full(shape=N, fill_value=mean_sd), label="Average Planetesimal Surface Density")
ax.scatter(x, surface_density, label="Planetesimal Surface Density, lin")
ax.scatter(x, surface_density_log, label="Planetesimal Surface Density, log")
ax.scatter(x, surface_density_self, label="Planetesimal Surface Density, self")
ax.set_xscale("log", base=10)
ax.set_yscale("log", base=10)
plt.legend()
plt.savefig("planetesimal_SD.png")

d0 = 80
d1 = -1


def Planetesimal_SurfaceDensity(r):
    return d0 * r ** d1 * r ** 2


import scipy.integrate as integrate

norm, _ = integrate.quad(Planetesimal_SurfaceDensity, R_min, R_max)

import scipy.stats as st


class PDF(st.rv_continuous):
    def _pdf(self, x):
        return Planetesimal_SurfaceDensity(x) / norm  # Normalized over its range, in this case [0,1]


# def PDF(r):
#     return Planetesimal_SurfaceDensity(r) / norm

my_cv = PDF(a=R_min, b=R_max, name='PDF')
sample = my_cv.rvs(size=N_planetesimals)

# unity, _ = integrate.quad(PDF, R_min, R_max)
# print("Unity Test: ", unity)

fig, ax = plt.subplots()
# ax.set_xlabel('Distance in AU')
# ax.set_ylabel('Probability')
ax.hist(my_cv.rvs(size=N_planetesimals * 100), bins=100, label="PDF")
# ax.set_xscale("log", base=10)
# ax.set_yscale("log", base=10)
plt.legend()
plt.savefig("PDF.png")

disk_ids = []

for item in sample:
    index = np.argmin(np.abs(x - item))
    disk_ids.append(index)

surface_density = np.zeros(N)

for i in range(N):
    sum = 0
    for j in range(N_planetesimals):
        if disk_ids[j] == i:
            sum += initial_mass
    surface_density[i] = sum / A[i]

# Convert to CGS
x *= au
dx *= au
A *= au ** 2


# Planetesimal Mass Distribution
# Radius of Planetesiomals in cm according to Coleman 2021

def Mass(radius):
    return 4 / 3 * np.pi * (radius) ** 3 * rho


# at log(1) au, radius log(20) km

# at log(10) au, radius log(100) km

dx = np.log10(10) - np.log10(1)
dy = np.log10(90) - np.log10(20)
dxdy = dx / dy
a0 = 100 * 1000 * 100 * np.power(10, np.log10(20) - dx / dy)
a1 = dx / dy
print("a0: ", a0)
print("a1: ", a1)


# initial Planetesimal Radius in cm
def Radius(r):
    # a0:  6548731.842129445
    # a1:  1.6609640474436815
    return a0 * r ** a1


b0 = a0 * (R_S / au) ** a1 / R_S
print("b0: ", b0)
print(6.1929286044739439e-15 * M_S / M_E * 2 * 10e6)


def Radius_S(r_s):
    # b0 = a0 * (R_S /au)**a1
    # b0 = 2.2745423088987758e-08
    return b0 * r_s ** a1


print(Radius(1) / R_S / R_E)
print(Radius_S(215.03215567054764) * R_S / R_E)

fig, ax = plt.subplots()
ax.set_xlabel('Distance in AU')
ax.set_ylabel('Initial Mass in Erath Masses')
ax.plot(radii, Mass(Radius(radii)) / M_E, label="Initial Masses of Planetesimals")
ax.set_xscale("log", base=10)
ax.set_yscale("log", base=10)
plt.legend()
plt.savefig("Embryos.png")

print(np.sum(Mass(Radius(radii))) / M_S)
print(Radius(1) / 100 / 1000)
print(Mass(Radius(30)))

TM_planetesimals = np.sum(initial_mass * N_planetesimals)

print("Total Mass of Planetesimals: ", TM_planetesimals / M_J)

# Planetesimal Surface Density

dr = []
for i, item in enumerate(radii[:-1]):
    delta = radii[i + 1] - item
    dr.append(delta)

dr.insert(0, radii[0] - R_min)

areas = np.zeros(N_planetesimals)
for i in range(N_planetesimals - 1):
    areas[i] = np.pi * ((radii[i + 1]) ** 2 - (radii[i]) ** 2)
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
    y = a2 + a0 * rad ** coeff
    return y


popt, _ = curve_fit(pla_dens, radii, surface_dens)

print(popt)
