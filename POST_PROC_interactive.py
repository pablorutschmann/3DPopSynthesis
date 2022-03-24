import Synthesis.post as post
from Synthesis.post.plot import *
from Synthesis.units import *
from importlib import reload

POP = post.population('SynthesisRuns/combined', 242)
POP.switch_plot_config('paper')

print(POP)

# Limits for Terrestrial Planets
m_low_lim = 1 * M_ME / M_E
a_up_lim = 2

# m_low_lim = 0
# a_up_lim = 30

thresholds = [0.05,0.1,0.15]


# New and Improved


## GLOBAL STATISTICS
scatter.scatter_parameters(POP)
scatter.scatter_ecc_inc(POP)
scatter.scatter_a_mass(POP)
histogram.histogram_mass(POP)
histogram.histogram_weighted_mass(POP)
histogram.histogram_a(POP)
histogram.histogram_weighted_a(POP)

misc.line_n_planets(POP)
histogram.histogram_totalmass(POP)

## COLLISIONS
collision.collisions_orb_dist(POP)
collision.collisions_time(POP)
collision.collisions_mass(POP)
collision.collisions_engulfed(POP)


## TERRESTRIAL REGION STATISTICS
histogram.histogram_mass(POP, m_low_lim, a_up_lim)
histogram.histogram_weighted_mass(POP, m_low_lim, a_up_lim)
histogram.histogram_a(POP, m_low_lim, a_up_lim)
histogram.histogram_weighted_a(POP, m_low_lim, a_up_lim)
histogram.histogram_totalmass(POP, m_low_lim)
histogram.histogram_totalmass_thresh(a_up_lim)
misc.line_n_planets(a_up_lim,thresholds)



#distance.scatter_AMD_RMC(POP, m_low_lim, a_up_lim)




# old and improved
#distance.distances(POP, m_low_lim,a_up_lim)


# old
histogram.histogram_wmf(POP, m_low_lim,a_up_lim)
scatter.scatter_radial_twmf(POP, m_low_lim,a_up_lim)

histogram.histogram_a_twmf(POP, m_low_lim,a_up_lim)
histogram.histogram_a_wmf(POP, m_low_lim,a_up_lim)

histogram.histogram_weighted_mass_nonlog(POP)


#final_wmf_radial_distribution(thresholds=[0,1.0e-11])
# final_wmf_radial_distribution(cumulative=True, thresholds=[0,1.0e-11])
#final_wmf_radial_distribution(cumulative=False, density=True, thresholds=[1.0e-11,9.0e-11])
# final_wmf_radial_distribution(cumulative=True, density=True, thresholds=[1.0e-11,9.0e-11])
