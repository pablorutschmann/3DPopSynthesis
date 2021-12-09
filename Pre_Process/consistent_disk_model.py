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


# Planetesimals

N_planetesimals = 100

# rho = 0.93357 * M_S / R_S ** 3
rho = 0.93357
# initial_mass = 3.0e-08 * M_S

initial_radius = 200/R_E * R_E *1000 *100 # km

initial_mass = 4 / 3 * np.pi * (initial_radius) ** 3 * rho

mass = 1.6e-11 * M_S
initR = np.abs(np.power(mass * 3 / 4 / np.pi / rho,1/3))
print(mass/M_E)
print(initR/1000/100)

print(0.8 *au/R_S)

print("initial Mass in M_E: ", initial_mass / M_E)
print("initial Mass in M_S: ", initial_mass / M_S)
print("initial Radius in km: ", initial_radius/1000/100)
print("initial Radius in R_E: ", initial_radius/R_E)

TM_Planetesimals = N_planetesimals * initial_mass

print("Total Gas Mass: ", TM_gas / M_J)
print("Total Dust Mass: ", TM_dust / M_J)
print("Total Planetesimals Mass: ", TM_Planetesimals / M_J)

TM_disk = TM_dust + TM_gas + TM_Planetesimals
print("Total Mass: ", TM_disk / M_J)

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


def get_sd(arr):
    disk_ids = []

    for item in arr:
        index = np.argmin(np.abs(x - item))
        disk_ids.append(index)

    surface_density = np.zeros(N)
    for i in range(N):
        sum = 0
        for j in range(len(arr)):
            if disk_ids[j] == i:
                sum += initial_mass
                surface_density[i] = sum / (A[i] * au ** 2)
    return surface_density

def get_cum_sd(arr):
    disk_ids = []

    for item in arr:
        index = np.argmin(np.abs(x - item))
        disk_ids.append(index)

    surface_density = np.zeros(N)
    Area = 0
    Mass = 0
    for i in range(N):
        Area += A[i] * au ** 2
        for j in range(len(arr)):
            if disk_ids[j] == i:
                Mass += initial_mass
        surface_density[i] = Mass/Area
    return surface_density

surface_density = get_sd(radii)
surface_density_log = get_sd(radii_log)

mean_sd = N_planetesimals * initial_mass / (np.sum(A) * au ** 2)

grad = np.gradient(get_cum_sd(radii))
grad = np.gradient(get_cum_sd(radii_log))


fig, ax = plt.subplots()
ax.set_xlabel('Distance in AU')
ax.set_ylabel('Planetesimal Surface Density in CGS')
ax.set_xlim(0.5,30)
ax.plot(x, binkert(x), label="Binkert")
# ax.plot(x, hayashi(x), label="Hayashi")
# ax.plot(x, raymond(x), label="Raymond")
ax.plot(x, 0.01 * binkert(x), label="Binkert Dust")
# ax.plot(x, np.full(shape=N, fill_value=mean_sd), label="Average Planetesimal Surface Density")
#ax.scatter(x, surface_density, label="Planetesimal Surface Density, lin")
#ax.scatter(x, surface_density_log, label="Planetesimal Surface Density, log")
ax.plot(x, get_cum_sd(radii), label="Cumulative Planetesimal Surface Density, lin")
ax.plot(x, get_cum_sd(radii_log), label="Cumulative Planetesimal Surface Density, log")
#ax.hlines(y=mean_sd,xmin=R_min,xmax=R_max)
ax.set_xscale("log", base=10)
ax.set_yscale("log", base=10)
plt.legend()
plt.savefig("planetesimal_SD.png")
plt.close()


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
plt.close()

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

TM_planetesimals = np.sum(initial_mass * N_planetesimals)

print("Total Mass of Planetesimals: ", TM_planetesimals / M_J)