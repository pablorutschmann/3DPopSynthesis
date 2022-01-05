import numpy as np
#from .units import *

# Radius always in AU
# Use exactly these function names: WMF, Surface_Density, Temperature, Eta

def wmf(r):
    # From Morbidelli et al, 2012
    IceLine = 2.71280277
    if r < IceLine:
        wmf = 0.0
        if r < 1.5:
            swmf = 0.0
        elif 1.5 <= r < 2.0:
            swmf = 0.01
        else:
            swmf = 0.1

    elif IceLine <= r < IceLine + 0.3 :
        wmf = 0.2
        swmf = 0.2

    elif IceLine + 0.3 <= r:
        wmf = 0.4
        swmf = 0.2

    return wmf, swmf

WMF = np.vectorize(wmf)

def raymond(r,coeff):
    # coeff = -1.0
    coeff = np.random.uniform(-1.5,-0.5)
    y = 4000 * (r)**coeff
    return y, coeff

def binkert(r, coeff):
    norm = 100
    y = norm * r**coeff
    return y, coeff, norm

Surface_Density = binkert

def Temperature(r):
    coeff = -0.5
    y = 280 * (r)**coeff
    return y, coeff

def Eta(r):
    y = 1.5e-3 * r**(0.5)
    return y

if __name__ == "__main__":
    # aus = np.linspace(0.5,30,5)
    # print(aus)
    # for i in range(5):
    #     sd, coeff = Surface_Density(aus)
    #     print(sd)
    #     print(coeff)

    print(Temperature(2.71280277))



