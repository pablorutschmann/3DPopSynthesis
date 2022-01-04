from ..units import *
import numpy as np
from scipy.integrate import dblquad


# From Alibert 2017

def weight(M):
    return np.log10(M*M_E/M_J) + 1

def f_p(M, A, MP, AP):
    # Funtion for each planet
    # Variables: Mass and Semi Major Axis
    # Parameter: Planet in form of Pandas Row
    sigma_m = 1
    sigma_a = 0.2

    return np.exp(-np.power((np.log10(M) - np.log10(MP)) / (sigma_m), 2) - np.power((np.log10(A) - np.log10(AP)) / (sigma_a), 2)) * weight(MP)


def PSI_S(M, A, System):
    # Variables: Mass and Semi Major Axis
    # Parameters: List of Tuples containing (MP, AP)

    sum = 0.0

    for (MP, AP) in System:
        sum += f_p(M, A, MP, AP)

    return sum


def distance(s1, s2):
    def f(M, A):
        return np.power(PSI_S(M, A, s1) - PSI_S(M, A, s2), 2) / M / A

    integral, err = dblquad(f, 0.1, 30, 0.1, 50)

    return np.sqrt(integral)


# Solar System Planets (Mass, Orbital Distance) in Earth Units

# Terrestrial Planets
mercury = (3.3010e26 / M_E, 5.7909227e10 / au)
venus = (4.8673e27 / M_E,1.0820948e11 / au)
earth = (1, 1)
mars = (6.4169e26 / M_E , 2.2794382e11 / au)

# Gas and Ice Giant Planets
jupiter = (M_J / M_E, 7.7834082e11 / au)
saturn = (5.6832e29 / M_E, 1.4266664e12 / au)
uranus = (8.6810e28 / M_E, 2.8706582e12 / au)
neptune = (1.0241e29, 4.4983964e12 / au)

terrestrial = [mercury, venus, earth, mars]

solar_system = [mercury, venus, earth, mars, jupiter, saturn, uranus, neptune]

