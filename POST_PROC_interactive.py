import Synthesis.post as post
from Synthesis.post.plot import *
from Synthesis.units import *
from importlib import reload

POP = post.population('SynthesisRuns/combined',145)
POPULATION.switch_plot_config('paper')

print(POPULATION)

# Limits for Terrestrial Planets
m_low_lim = 1 * M_ME / M_E
a_up_lim = 2

# m_low_lim = 0
# a_up_lim = 30

thresholds = [0.05,0.1,0.15]


# New and Improved


## GLOBAL STATISTICS
scatter.scatter_parameters(POPULATION)
scatter.scatter_ecc_inc(POPULATION)
histogram.histogram_mass(POPULATION)
histogram.histogram_weighted_mass(POPULATION)
histogram.histogram_a(POPULATION)
histogram.histogram_weighted_a(POPULATION)

misc.line_n_planets(POPULATION)
histogram.histogram_totalmass(POPULATION)

## COLLISIONS
collision.collisions_orb_dist(POPULATION)
collision.collisions_time(POPULATION)
collision.collisions_mass(POPULATION)
collision.collisions_engulfed(POPULATION)


## TERRESTRIAL REGION STATISTICS
histogram.histogram_mass(POPULATION, m_low_lim, a_up_lim)
histogram.histogram_weighted_mass(POPULATION, m_low_lim, a_up_lim)
histogram.histogram_a(POPULATION, m_low_lim, a_up_lim)
histogram.histogram_weighted_a(POPULATION, m_low_lim, a_up_lim)
histogram.histogram_totalmass(POPULATION, m_low_lim)
histogram.histogram_totalmass_thresh(a_up_lim)
misc.line_n_planets(a_up_lim,thresholds)



#distance.scatter_AMD_RMC(POPULATION, m_low_lim, a_up_lim)




# old and improved
#distance.distances(POPULATION, m_low_lim,a_up_lim)


# old
histogram.histogram_wmf(POPULATION, m_low_lim,a_up_lim)
scatter.scatter_radial_twmf(POPULATION, m_low_lim,a_up_lim)

histogram.histogram_a_twmf(POPULATION, m_low_lim,a_up_lim)
histogram.histogram_a_wmf(POPULATION, m_low_lim,a_up_lim)

histogram.histogram_weighted_mass_nonlog(POPULATION)


#final_wmf_radial_distribution(thresholds=[0,1.0e-11])
# final_wmf_radial_distribution(cumulative=True, thresholds=[0,1.0e-11])
#final_wmf_radial_distribution(cumulative=False, density=True, thresholds=[1.0e-11,9.0e-11])
# final_wmf_radial_distribution(cumulative=True, density=True, thresholds=[1.0e-11,9.0e-11])
