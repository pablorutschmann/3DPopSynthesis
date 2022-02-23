import Synthesis.post as post
from Synthesis.units import *
from importlib import reload

from Synthesis.post.plot import *

TEST = post.population('SynthesisRuns/combined',3)
TEST.switch_plot_config('paper')

# Limits for Terrestrial Planets
m_low_lim = 1 * M_ME / M_E
a_up_lim = 2

# m_low_lim = 0
# a_up_lim = 30

thresholds = [0.05,0.1,0.15]