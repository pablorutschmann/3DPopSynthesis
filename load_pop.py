import Synthesis.post as post
from Synthesis.units import *
from importlib import reload

from Synthesis.post.plot import *

POPULATION = post.population('SynthesisRuns/combined',145)
POPULATION.switch_plot_config('paper')

# Limits for Terrestrial Planets
m_low_lim = 1 * M_ME / M_E
a_up_lim = 2

# m_low_lim = 0
# a_up_lim = 30

thresholds = [0.05,0.1,0.15]


## GLOBAL STATISTICS
scatter.scatter_parameters(POPULATION)
scatter.scatter_ecc_inc(POPULATION)
histogram.histogram_mass(POPULATION)
histogram.histogram_weighted_mass(POPULATION)
histogram.histogram_a(POPULATION)
histogram.histogram_weighted_a(POPULATION)

misc.line_n_planets(POPULATION,a_up_lim=30)
histogram.histogram_totalmass(POPULATION)

## COLLISIONS
collision.collisions_orb_dist(POPULATION)
collision.collisions_time(POPULATION, t_low_lim=100000)
collision.collisions_mass(POPULATION)
collision.collisions_engulfed(POPULATION)


## TERRESTRIAL REGION STATISTICS
histogram.histogram_mass(POPULATION, m_low_lim, a_up_lim)
histogram.histogram_weighted_mass(POPULATION, m_low_lim, a_up_lim)
histogram.histogram_a(POPULATION, m_low_lim, a_up_lim)
histogram.histogram_weighted_a(POPULATION, m_low_lim, a_up_lim)
histogram.histogram_totalmass(POPULATION, m_low_lim)
histogram.histogram_totalmass_thresh(POPULATION, a_up_lim)
misc.line_n_planets(POPULATION, a_up_lim, thresholds)
