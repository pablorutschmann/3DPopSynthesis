import Synthesis.post as post
import Synthesis.post.plot as plt
from Synthesis.units import *



POPULATION = post.population('SynthesisRuns/combined',3)
POPULATION.switch_plot_config('paper')

print(TEST)

# Limits for Terrestrial Planets
m_low_lim = 1 * M_ME / M_E
a_up_lim = 2

# m_low_lim = 0
# a_up_lim = 30

thresholds = [0.05,0.1,0.15]


# New and Improved


## GLOBAL STATISTICS
plt.scatter.scatter_parameters(POPULATION)
plt.scatter.scatter_ecc_inc(POPULATION)
plt.histogram.histogram_mass(POPULATION)
plt.histogram.histogram_weighted_mass(POPULATION)
plt.histogram.histogram_a(POPULATION)
plt.histogram.histogram_weighted_a(POPULATION)

plt.misc.line_n_planets(POPULATION)
plt.histogram.histogram_totalmass(POPULATION)

## COLLISIONS
plt.collision.collisions_orb_dist(POPULATION)
plt.collision.collisions_time(POPULATION)
plt.collision.collisions_mass(POPULATION)
plt.collision.collisions_engulfed(POPULATION)


## TERRESTRIAL REGION STATISTICS
plt.histogram.histogram_mass(POPULATION, m_low_lim, a_up_lim)
plt.histogram.histogram_weighted_mass(POPULATION, m_low_lim, a_up_lim)
plt.histogram.histogram_a(POPULATION, m_low_lim, a_up_lim)
plt.histogram.histogram_weighted_a(POPULATION, m_low_lim, a_up_lim)
plt.histogram.histogram_totalmass(POPULATION, m_low_lim)
plt.histogram.histogram_totalmass_thresh(a_up_lim)
plt.misc.line_n_planets(a_up_lim,thresholds)



#plt.distance.scatter_AMD_RMC(POPULATION, m_low_lim, a_up_lim)




# old and improved
#plt.distance.distances(POPULATION, m_low_lim,a_up_lim)


# old
plt.histogram.histogram_wmf(POPULATION, m_low_lim,a_up_lim)
plt.scatter.scatter_radial_twmf(POPULATION, m_low_lim,a_up_lim)

plt.histogram.histogram_a_twmf(POPULATION, m_low_lim,a_up_lim)
plt.histogram.histogram_a_wmf(POPULATION, m_low_lim,a_up_lim)

plt.histogram.histogram_weighted_mass_nonlog(POPULATION)


#final_wmf_radial_distribution(thresholds=[0,1.0e-11])
# final_wmf_radial_distribution(cumulative=True, thresholds=[0,1.0e-11])
#final_wmf_radial_distribution(cumulative=False, density=True, thresholds=[1.0e-11,9.0e-11])
# final_wmf_radial_distribution(cumulative=True, density=True, thresholds=[1.0e-11,9.0e-11])
