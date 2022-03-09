import Synthesis.post as post
from Synthesis.post.metric import distance
from Synthesis.units import *
from tqdm import tqdm
from os.path import join

POP = post.population('SynthesisRuns/combined', 242)
POP.switch_plot_config('paper')

print(POP)

# Limits for Terrestrial Planets
m_low_lim = 1 * M_ME / M_E
a_up_lim = 2

# m_low_lim = 0
# a_up_lim = 30

thresholds = [0.05, 0.1, 0.15]


# tsne(POP,m_low_lim,a_up_lim)

def get_dist(pop, mlow, aup):
    systems = []
    sys_matrix = []
    total_masses = []
    sigma_exponents = []
    TotalMasses = []
    SigmaCoeffs = []
    Reference = []
    RMCS = []
    AMDS = []
    Numbers = []

    for id, sim in tqdm(pop.SIMS.items()):
        if id > 250:
            break
        TotalMasses.append(sim.Total_Mass)
        SigmaCoeffs.append(sim.Sigma_Exponent)

    for id, sim in tqdm(pop.SIMS.items()):
        if id > 250:
            break
        total_masses.append(sim.Total_Mass)
        sigma_exponents.append(sim.Sigma_Exponent)

        zipped = zip(sim.snaps[sim.N_snaps - 1].satellites['M'].values * M_S / M_E,
                     sim.snaps[sim.N_snaps - 1].satellites['a'].values * R_S / au)

        # Remove under threshold Masses
        filtered = [(m, a) for (m, a) in zipped if
                    m >= m_low_lim and a <= a_up_lim]
        # Numbers.append(len(filtered))
        # Sorted = sorted(filtered, key=lambda x: x[0])[:4]
        # Sorted = sorted(Sorted,key=lambda x: x[1])
        Numbers.append(len(filtered))
        cum_mass = np.sum([m for (m, a) in filtered])
        data_point = []
        for item in filtered:
            # print([*item])
            data_point.extend([*item])
        systems.append(data_point)
        # Reference.append(sim.Sigma_Norm * pow(R_S / (0.3 * au), sim.Sigma_Exponent) * (M_S / R_S2))
        system = [(m, a) for (m, a) in filtered]
        # print(system)
        sys_matrix.append(system)
        # Reference.append(distance(system, terrestrial))
        Reference.append(sim.Sigma_Norm / (R_S / au) ** sim.Sigma_Exponent / denstos * pow(au / R_S, sim.Sigma_Exponent) / denstos)
        # Reference.append(cum_mass)
    Numbers = Reference

    terr = [(m, a) for (m, a) in terrestrial]
    sys_matrix.append(terr)
    num = len(sys_matrix)

    # def get_distances(i):
    #     dist =  np.zeros(num)
    #     for j in range(i,num):
    #         if i == j:
    #             dist[j] = 0.0
    #         else:
    #             dist[j] = distance(sys_matrix[i],sys_matrix[j])
    #     return dist
    #
    #
    #     for j in range(i,num):
    #         if i == j:
    #             dist[j] = 0.0
    #         else:
    #             dist[j] = distance(sys_matrix[i],sys_matrix[j])
    #     return dist

    def get_distances(n, arr):
        out = np.zeros(n)
        sys = arr[0]
        others = arr[1:]

        def dist_arr(other):
            return distance(sys, other)

        dist_func = np.vectorize(dist_arr)

        distances_array = dist_func(others)

        # out = np.pad(distances_array, (0,0) , 'constant', constant_values=(n-len(others), 0))
        out[-len(distances_array):] = distances_array
        print(len(others))

        return out

    # rows = np.array([np.array(sys_matrix[i:], dtype='object') for i in range(num-1)], dtype='object')
    # print(rows)
    #
    # upper = np.array(list(map(lambda x :get_distances(arr=x,n=num),rows)))
    # print(upper)
    # upper = np.vstack([upper,np.zeros(num)])
    #
    # upper = np.reshape(upper,(num,num))
    # matrix_metric = upper + upper.T

    rows = np.array([np.array(sys_matrix[i:], dtype='object') for i in range(num - 1)], dtype='object')
    print(rows)

    upper = np.array(list(map(lambda x: get_distances(arr=x, n=num), rows)))
    print(upper)
    upper = np.vstack([upper, np.zeros(num)])

    upper = np.reshape(upper, (num, num))
    matrix_metric = upper + upper.T

    np.save(join(pop.PLOT, f'metric_matrix_log_M_{num}'), matrix_metric)
    print('saved to Pickle')
    print(matrix_metric)


get_dist(POP, m_low_lim, a_up_lim)
