import numpy as np
import matplotlib.pyplot as plt
import astropy.constants as astroconst
from astropy import units as u

# Solar Radius in grams
R_S = astroconst.R_sun.decompose(u.cgs.bases).value

# AU in cm
au = astroconst.au.decompose(u.cgs.bases).value

class plotter:
    def __init__(self, disk):
        self.disk = disk

    def plot_density(self,scale = 'log'):

        x_init = self.disk.gas_profiles['Radius [R_S]']
        x_out = self.disk.out['r [R_S]']
        x_all = np.linspace(x_out.min(),x_init.max(),1000)

        unit = ' in M_S / R_S^2'
        type = ' Surface_Density'

        for phase in ['Gas','Dust']:
            fig, ax = plt.subplots()
            ax.set_xlabel('Radius in AU')
            ax.set_title(type + ' profile')
            # ax.loglog(r_resh[0,:,0]/au, sigma1d)
            if phase == 'Gas':
                profile = self.disk.gas_profiles['Surface Density [M_S / R_S^2]']
                coeffs = self.disk.gas_profiles_extra['Power Coefficient Sigma']
                print(coeffs)
                ax.plot(x_all  * R_S / au, self.disk.binkert(x_all), label="Binkert")
                ax.plot(x_all  * R_S / au, self.disk.raymond(x_all), label="Raymond")
            if phase == 'Dust':
                profile = self.disk.dust_profiles['Surface Density [M_S / R_S^2]']
                coeffs = self.disk.dust_profiles_extra['Power Coefficient Sigma']
                ax.plot(x_all  * R_S / au, self.disk.binkert(x_all)*0.01, label="Binkert")
                ax.plot(x_all  * R_S / au, self.disk.raymond(x_all)*0.01, label="Raymond")

            #ax.plot(x_init * R_S / au, profile, 'ko', label="Original Data")
            #ax.plot(x_out * R_S / au, self.disk.func_exp(x_out, *coeffs), 'r-', label="Disk Model")
            ax.set_ylabel(phase + ' '  + type + unit)
            #ax.plot(x_init * R_S / au, self.disk.func_exp(x_init, *coeffs), label="Fitted Curve")
            plt.legend()
            if scale == 'log':
                ax.set_xscale("log", base=10)
                ax.set_yscale("log", base=10)
                plt.savefig('Disk_Plots/log_' + phase + '_'  + type + '_' + self.disk.spacing + '_' + str(self.disk.N) + '.png')
            elif scale == 'lin':
                plt.savefig('Disk_Plots/' + phase + '_'  + type + '_' + self.disk.spacing + '_' + str(self.disk.N) + '.png')


