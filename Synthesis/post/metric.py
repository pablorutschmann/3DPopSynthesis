from Synthesis.units import *
import numpy as np
from scipy.integrate import dblquad
from scipy.integrate import tplquad


# From Alibert 2017

def weight(M,A):
    return np.log10(M) - np.log(M_ME/M_E)
    #return M
    # return 1 / (1 + np.exp(-0.5 * ((np.log10(M)-np.log10(0.4818825))/1)**2))
    # return 1 / (1 + np.exp(-0.5 * ((np.log10(M)-np.log10(0.4818825))/1)))
    # return 1 / (1 + np.exp(- (M-1)))
    # return 1 / (1 + np.exp(- (M-1))) * 1 / (1 + np.exp(- (-A+3)/0.7))
    # return 1 / (1 + np.exp(- (-M+3)))
    # return np.exp(-0.5 * ((np.log10(M)-np.log10(0.4818825))/10)**2)
    # return 1

def f_p(M, A, MP, AP):
    # Function for each planet
    # Variables: Mass and Semi Major Axis
    # Parameter: Planet in form of Pandas Row
    sigma_m = 1/0.3
    sigma_a = 1/0.3
    out = np.exp(-0.5 * np.power((np.log10(M) - np.log10(MP)) / (sigma_m), 2) - 0.5 * np.power(
        (np.log10(A) - np.log10(AP)) / (sigma_a), 2)) * weight(MP,AP)
    return out


def PSI_S(M, A, System):
    # Variables: Mass and Semi Major Axis
    # Parameters: List of Tuples containing (MP, AP)

    sum = 0.0

    for (MP, AP) in System:
        sum += f_p(M, A, MP, AP)

    return sum


def distance(s1, s2):
    tm1 = np.sum([item[0] for item in s1])
    tm2 = np.sum([item[0] for item in s2])
    def f(M, A):
        return np.power(PSI_S(M, A, s1)/tm1 - PSI_S(M, A, s2)/tm2, 2) / M / A

    integral, err = dblquad(f, 0.001, 100, 0.1, 50)

    return np.sqrt(integral)

average_mass = np.mean([m for (m,a) in terrestrial])
print(average_mass)
#
# print(PSI_S(1,1,terrestrial))
#
# print(PSI_S(10,2,terrestrial))
#
#
# print(1)
# print(distance([(0.3,0.3)],terrestrial))
#
# print(distance(terrestrial + [(0.3,0.3)],terrestrial))
#
# print(distance(terrestrial + [(10,2)],terrestrial))
#
# print(distance([(0.3,0.3), (0.8,1.2)],terrestrial))
#
# print(distance([(1,1)],terrestrial))
#
# print(distance([(10,2)],terrestrial))
#
# print(distance([(60,3)],terrestrial))
#
# print(distance([(10,3)],terrestrial))
#
# print(distance([(1,1),(10,3)],terrestrial))
#
# print(distance([(0.3,1)],terrestrial))


def weight_wmf(m,wmf):
    # return (1 / (1 + np.exp(wmf/0.001)) + 0.5) * (np.log10(m) - np.log(M_ME/M_E))
    # return (1 - np.exp(- 0.5 * (x-0.001)**2 / 0.1**2)) * (np.log10(m) - np.log(M_ME/M_E))
    #return np.log10(m) - np.log(M_ME/M_E)
    return 1


def earth_f_p(M, A, MP, AP, WMFP):
    # Function for each planet
    # Variables: Mass and Semi Major Axis
    # Parameter: Planet in form of Pandas Row
    sigma_m = 1 / 0.3
    sigma_a = 1 / 0.3

    return np.exp(-np.power((np.log10(M) - np.log10(MP)) / (2 * sigma_m), 2) - np.power(
        (np.log10(A) - np.log10(AP)) / (2 * sigma_a), 2))


def earth_PSI_S(M, A, System):
    # Variables: Mass and Semi Major Axis
    # Parameters: List of Tuples containing (MP, AP)

    sum = 0.0

    for (MP, AP, WMFP) in System:
        sum += earth_f_p(M, A, MP, AP, WMFP)

    return sum


def earth_distance(s1, s2):
    def f(M, A):
        return np.power(earth_PSI_S(M, A, s1) - earth_PSI_S(M, A, s2), 2) / M / A

    integral, err = dblquad(f, 0.001, 100, 0.1, 50)

    return np.sqrt(integral)

print(earth_distance([(1,0.8,0.1)],[(1,1,0.001)]))

print(earth_distance([(0.8,1,0.002)],[(1,1,0.001)]))

print(earth_distance([(1.8,1,0.001)],[(1,1,0.001)]))

