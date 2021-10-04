import disk

mydisk = disk.disk('ppd/')

mydisk.gas.to_pickle('ppd/lev0_gas.dat.pickle')
mydisk.dust.to_pickle('ppd/lev0_dust.dat.pickle')

mydisk.N_front = 7
mydisk.N_back = 6

mydisk.integrate_1D('gas')
mydisk.integrate_1D('dust')

spacing = 'log'
R_min = 0.5 # Minimum Radius in AU for extrapolation
R_max = 5 # Maximum Radius in AU for extrapolation
N = 1000 # Number of cells in extrapolation axis

mydisk.extrapolate(spacing, R_min, R_max, N)

mydisk.write_disk()