import Synthesis.post as post
from Synthesis.post.metric import distance
from Synthesis.units import *
from tqdm import tqdm
from os.path import join

POP = post.population('SynthesisRuns/combined',242)
POP.switch_plot_config('paper')

print(POP)

# Limits for Terrestrial Planets
m_low_lim = 1 * M_ME / M_E
a_up_lim = 2

# m_low_lim = 0
# a_up_lim = 30

thresholds = [0.05,0.1,0.15]

#tsne(POP,m_low_lim,a_up_lim)

def get_dist(pop,mlow,aup):
    systems = []
    total_masses = []
    sigma_exponents = []


    for id, sim in tqdm(pop.SIMS.items()):
        if id > 250:
            break
        total_masses.append(sim.Total_Mass)
        sigma_exponents.append(sim.Sigma_Exponent)

        zipped = zip(sim.snaps[sim.N_snaps - 1].satellites['M'].values,
                     sim.snaps[sim.N_snaps - 1].satellites['a'].values)

        # Remove under threshold Masses
        filtered = [(m * M_S / M_E, a * R_S / au) for (m, a) in zipped if
                    m >= mlow * M_E / M_S and a <= aup * au / R_S]
        systems.append(filtered)
    num=len(systems)
    distances = []
    # print(systems)
    # func = lambda x: distance(x, terrestrial)
    distances = np.array(list(map(lambda x :distance(x,terrestrial),systems)))

    np.save(join(pop.PLOT,f'dist_terr_log_M_{num}'), distances)
    print('saved to Pickle')
    print(distances)


get_dist(POP,m_low_lim,a_up_lim)
