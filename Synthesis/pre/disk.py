import numpy as np
from ..units import *
import pandas as pd
import os.path as path
from ..functions import *

class disk:
    def __init__(self):
        # Parameters
        self.spacing = 'None'
        self.R_min = -1 # Minimum Radius in AU for extrapolation
        self.R_max = -1 # Maximum Radius in AU for extrapolation
        self.N = -1 # Number of cells in extrapolation axis

        self.tot_mass_gas = -1
        self.tot_mass_dust = -1
        self.tot_mass = -1

        self.out = {}

    def prepare(self, spacing, R_min, R_max, N):
        self.spacing = spacing
        self.R_min = R_min
        self.R_max = R_max
        self.N = N

        # Creation of linspace or logspace of radii in R_S, include one point more for N cells
        if self.spacing == 'lin':
            x = np.linspace(self.R_min, self.R_max, self.N + 1) * autor
        elif self.spacing == 'log':
            a = np.power(self.R_max/self.R_min, 1/self.N)
            x = np.geomspace(self.R_min,self.R_max,self.N) * autor
            dx = [(a-1) * item for item in x]
            dx.append((a-1)*x[-1]*a)
            dx.insert(0,(a-1)*x/a)

        # Chache Radii in AU for Sampling (conversion before file writing)
        self.out['r [R_S]'] = x * rtoau
        dx.pop(0)
        dx.pop(-1)
        self.out['dr [R_S]'] = np.array(dx)

        # Get Keplerian Velocity
        vx = kepl_velo(x)
        self.out['Keplerian Velocity [R_S / yr]'] = vx

        # Get Area of Annulus
        A = Areas(x,dx)
        self.out['Area [R_S^2]'] = A

        # Get Opacity
        self.out['Gas Opacity []'] = np.zeros(N)

    def sample(self, INPUT):
        # Get Surface Densities and Power Coefficients
        Sigma_Gas, Sigma_Coeff = Surface_Density(self.out['r [R_S]'])
        # Unit Conversion from CGS to Solar Units
        Sigma_Gas *= denstos
        Sigma_Dust = 0.01 * Sigma_Gas
        self.out['sigma gas [M_S/R_S^2]'] = Sigma_Gas
        self.out['sigma dust [M_S/R_S^2]'] = Sigma_Dust
        self.out['sigma dustbar [M_S/R_S^2]'] = Sigma_Dust
        self.out['Power Coefficient Density'] = np.full(self.N, Sigma_Coeff)

        # Get Temperature
        self.out['T [K]'], Temp_Coeff = Temperature(self.out['r [R_S]'])
        self.out['Power Coefficient Temperature'] = np.full(self.N, Temp_Coeff)

        # Get Water Mass Fractions
        self.out['WMF []'],self.out['SWMF []'] = WMF(self.out['r [R_S]'])

        #Get Headwind Factor
        self.out["Headwind Factor []"] = Eta(self.out['r [R_S]'])

        # Write out the disk file for that Sample in the inputs folder
        self.write_disk(INPUT)

    def write_disk(self,INPUT):

        for key in ['r [R_S]', 'dr [R_S]', 'sigma gas [M_S/R_S^2]', 'sigma dust [M_S/R_S^2]', 'sigma dustbar [M_S/R_S^2]', 'T [K]', 'Area [R_S^2]', 'Keplerian Velocity [R_S / yr]', 'Power Coefficient Density', 'Power Coefficient Temperature', 'Gas Opacity []', 'WMF []', 'SWMF []', 'Headwind Factor []']:
            #self.out = self.out.pop[key]
            self.out[key] = self.out.pop(key)

        df_out = pd.DataFrame.from_dict(self.out)
        df_out.index.name = 'index'
        #Convert Radius back to Solar Units
        df_out['r [R_S]'] *= autor

        SAVEPATH = path.join(INPUT,"disk.txt")
        df_out.to_csv(SAVEPATH, header=False, sep='\t', mode='w')


# HELPER FUCNTIONS

# Calculate Keplerian Veloctiy
def kepl_velo(r):
    # convert r to cgs
    r = r * R_S
    v_kep = np.sqrt(G * M_S / r ** 3)
    # convert back to Solar Units
    return v_kep * year

# Calculate Area of Annulus
def Areas(x,dx):
    A = np.zeros(len(dx))
    for i in range(len(dx)):
        A[i] = np.pi * ((x[i] + dx[i]/2)**2 - (x[i] - dx[i]/2)** 2)
    return A


if __name__ == "__main__":
    foo = 1