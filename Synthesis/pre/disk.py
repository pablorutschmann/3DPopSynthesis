import numpy as np
from ..units import *
import pandas as pd
import os.path as path
from ..functions import *
from scipy.integrate import quad


class disk:
    def __init__(self):
        # Parameters
        self.spacing = 'None'
        self.R_min = -1  # Minimum Radius in AU for extrapolation
        self.R_max = -1  # Maximum Radius in AU for extrapolation
        self.N = -1  # Number of cells in extrapolation axis

        self.tot_mass_gas = -1
        self.tot_mass_dust = -1
        self.tot_mass = -1

        self.DustGasRatio = 0.01

        self.out = {}

    def prepare(self, TM_min, TM_max, R_min, R_max, N, spacing, Sigma_min, Sigma_max, N_Sigma):
        self.TM_min = TM_min
        self.TM_max = TM_max
        self.R_min = R_min
        self.R_max = R_max
        self.N = N
        self.spacing = spacing
        self.Sigma_min = Sigma_min
        self.Sigma_max = Sigma_max
        self.N_Sigma = N_Sigma
        self.Sigmas = np.round(np.linspace(self.Sigma_min,self.Sigma_max, self.N_Sigma),2)

        # Creation of linspace or logspace of radii in R_S, include one point more for N cells
        if self.spacing == 'lin':
            x = np.linspace(self.R_min, self.R_max, self.N + 1) * autor
        elif self.spacing == 'log':
            a = np.power(self.R_max / self.R_min, 1 / self.N)
            x = np.geomspace(self.R_min, self.R_max, self.N) * autor
            dx = [(a - 1) * item for item in x]
            dx.append((a - 1) * x[-1] * a)
            dx.insert(0, (a - 1) * x / a)

        # Chache Radii in AU for Sampling (conversion before file writing)
        self.out['r [R_S]'] = x * rtoau
        dx.pop(0)
        dx.pop(-1)
        self.out['dr [R_S]'] = np.array(dx)

        # Get Keplerian Velocity
        vx = kepl_velo(x)
        self.out['Keplerian Velocity [R_S / yr]'] = vx

        # Get Area of Annulus
        A = Areas(x, dx)
        self.out['Area [R_S^2]'] = A

        # Get Opacity
        self.out['Gas Opacity []'] = np.zeros(N)

    def sample(self, INPUT):
        # Get Surface Densities and Power Coefficients with constant total mass
        Sigma_Coeff = np.random.choice(self.Sigmas)
        integral, _ = quad(lambda x: Power_Law(x, 1, Sigma_Coeff) * x, self.R_min * au / R_S,
                           self.R_max * au / R_S)
        TotalMass = np.random.uniform(self.TM_min, self.TM_max)
        Sigma_Norm = TotalMass / (2 * np.pi * integral)

        Sigma_Gas = Power_Law(self.out['r [R_S]'] * au / R_S, Sigma_Norm, Sigma_Coeff)
        # Unit Conversion from CGS to Solar Units
        Sigma_Gas = Sigma_Gas
        Sigma_Dust = self.DustGasRatio * Sigma_Gas
        self.out['sigma gas [M_S/R_S^2]'] = Sigma_Gas
        self.out['sigma dust [M_S/R_S^2]'] = Sigma_Dust
        self.out['sigma dustbar [M_S/R_S^2]'] = Sigma_Dust

        # Get Temperature
        self.out['T [K]'], Temp_Coeff = Temperature(self.out['r [R_S]'])

        # Get Water Mass Fractions
        self.out['WMF []'], self.out['SWMF []'] = WMF(self.out['r [R_S]'])

        # Get Headwind Factor
        self.out["Headwind Factor []"] = Eta(self.out['r [R_S]'])

        # Write out the disk file for that Sample in the inputs folder
        self.write_disk(INPUT)

        return  TotalMass, Sigma_Coeff, Sigma_Norm * denstos * (R_S / au) ** Sigma_Coeff, Temp_Coeff

    def write_disk(self, INPUT):

        for key in ['r [R_S]', 'dr [R_S]', 'sigma gas [M_S/R_S^2]', 'sigma dust [M_S/R_S^2]',
                    'sigma dustbar [M_S/R_S^2]', 'T [K]', 'Area [R_S^2]', 'Keplerian Velocity [R_S / yr]',
                    'Gas Opacity []', 'WMF []', 'SWMF []',
                    'Headwind Factor []']:
            # self.out = self.out.pop[key]
            self.out[key] = self.out.pop(key)

        df_out = pd.DataFrame.from_dict(self.out)
        df_out.index.name = 'index'
        # Convert Radius back to Solar Units
        df_out['r [R_S]'] *= autor

        SAVEPATH = path.join(INPUT, "disk.txt")
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
def Areas(x, dx):
    A = np.zeros(len(dx))
    for i in range(len(dx)):
        A[i] = np.pi * ((x[i] + dx[i] / 2) ** 2 - (x[i] - dx[i] / 2) ** 2)
    return A


def Power_Law(x, a, b):
    return a * np.power(x, b)


if __name__ == "__main__":
    foo = 1
