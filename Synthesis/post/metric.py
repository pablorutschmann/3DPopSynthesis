from Synthesis.units import *
import numpy as np
from scipy.integrate import dblquad
from scipy.integrate import tplquad


# From Alibert 2017

def weight(M):
    return np.log10(M * M_E / M_J) + 1


def f_p(M, A, MP, AP):
    # Function for each planet
    # Variables: Mass and Semi Major Axis
    # Parameter: Planet in form of Pandas Row
    sigma_m = 1 / 0.3
    sigma_a = 1 / 0.3

    return np.exp(-np.power((np.log10(M) - np.log10(MP)) / (2 * sigma_m), 2) - np.power(
        (np.log10(A) - np.log10(AP)) / (2 * sigma_a), 2))
    # * weight(MP)


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

    integral, err = dblquad(f, 0.01, 5, 0.1, 20)

    return np.sqrt(integral)


def earth_f_p(M, A, WMF, MP, AP, WMFP):
    # Function for each planet
    # Variables: Mass and Semi Major Axis
    # Parameter: Planet in form of Pandas Row
    sigma_m = 1 / 0.1
    sigma_a = 1 / 0.1
    sigma_wmf = 1 / 0.02

    return np.exp(-np.power((np.log10(M) - np.log10(MP)) / (2 * sigma_m), 2) - np.power(
        (np.log10(A) - np.log10(AP)) / (2 * sigma_a), 2) - np.power(
        (WMF - WMFP) / (2 * sigma_wmf), 2))
    # * weight(MP)


def earth_PSI_S(M, A, WMF, System):
    # Variables: Mass and Semi Major Axis
    # Parameters: List of Tuples containing (MP, AP)

    sum = 0.0

    for (MP, AP, WMFP) in System:
        sum += earth_f_p(M, A, WMF, MP, AP, WMFP)

    return sum


def earth_distance(s1, s2):
    def f(M, A, WMF):
        return np.power(earth_PSI_S(M, A, WMF, s1) - earth_PSI_S(M, A, WMF, s2), 2) / M / A

    integral, err = tplquad(f, 0.01, 5, 0.1, 20, 0, 1)

    return np.sqrt(integral)

