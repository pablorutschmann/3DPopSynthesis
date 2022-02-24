import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib import colors
import os.path as path
from Synthesis.post.metric import *
from tqdm import tqdm


def scatter_AMD_RMC(pop, m_low_lim = 0, a_up_lim = 30):
    RMC_sol = 89.9
    AMD_sol = 0.0018

    AMDS = []
    RMCS = []
    NS = []

    for id, sys in tqdm(pop.SIMS.items()):
        if id > 200:
            break
        else:
            AMD, N = sys.get_AMD(m_low_lim, a_up_lim)
            AMDS.append(AMD)
            RMC, N = sys.get_RMC(m_low_lim, a_up_lim)
            RMCS.append(RMC)
            NS.append(N)

    AMDS = np.array((AMDS))
    RMCS = np.array(RMCS)
    sigma_AMD = np.std(AMDS)
    sigma_RMC = np.std(RMCS)
    print(f'{sigma_AMD=}')
    print(f'{sigma_RMC=}')

    def dist_to_solar(rmc, amd):
        return np.sqrt((rmc - RMC_sol) ** 2 / sigma_RMC ** 2 + (amd - AMD_sol) ** 2 / sigma_AMD ** 2)

    distances = dist_to_solar(RMCS, AMDS)
    print(distances)

    cmap = pop.cmap_standart
    cmin = min(NS)
    cmax = max(NS)

    norm = colors.Normalize(cmin, cmax)

    plt.rcParams.update({'figure.autolayout': True})
    plt.style.use('seaborn-paper')
    plt.rcParams.update({'font.size': pop.fontsize})

    fig, ax = plt.subplots(figsize=pop.figsize)
    ax.scatter(RMCS / sigma_RMC, AMDS / sigma_AMD, c=NS, cmap=cmap, norm=norm, s=12)
    ax.scatter(RMC_sol / sigma_RMC, AMD_sol / sigma_AMD, c='red')
    ax.scatter(RMCS[np.argmin(distances)] / sigma_RMC, AMDS[np.argmin(distances)] / sigma_AMD,
               c=NS[np.argmin(distances)])
    x_labels = ax.get_xticklabels()
    plt.setp(x_labels, horizontalalignment='center')
    ax.set(xlabel='RMC', ylabel=r'AMC')
    fig.colorbar(cm.ScalarMappable(cmap=cmap, norm=norm), orientation='vertical', label=r'Number of Satellites',
                 ax=ax)
    ax.set_xscale('log')
    ax.set_yscale('log')
    if pop.plot_config == 'presentation':
        ax.set(title=r'Radial Mass Concentration and Angular Momentum Deficit')

    save_name = 'scatter_AMD_RMC'
    if a_up_lim < 30 and m_low_lim > 0:
        save_name += '_lim'
    fig.savefig(path.join(pop.PLOT, save_name + '.png'), transparent=False, dpi=pop.dpi, bbox_inches="tight")
    plt.close(fig)


def distances(pop, m_low_lim, a_up_lim):
    systems = []
    total_masses = []
    sigma_exponents = []
    m_low_lim = m_low_lim
    a_up_lim = a_up_lim

    for sim in pop.SIMS.values():
        total_masses.append(sim.Total_Mass)
        sigma_exponents.append(sim.Sigma_Exponent)

        zipped = zip(sim.snaps[sim.N_snaps - 1].satellites['M'].values,
                     sim.snaps[sim.N_snaps - 1].satellites['a'].values)

        # Remove under threshold Masses
        filtered = [(m * M_S / M_E, a * R_S / au) for (m, a) in zipped if
                    m >= m_low_lim * M_E / M_S and a <= a_up_lim * au / R_S]
        systems.append(filtered)

    func = lambda x: distance(x, terrestrial)
    distances = list(map(func, systems))

    plt.rcParams.update({'figure.autolayout': True})
    plt.style.use('seaborn-paper')
    plt.rcParams.update({'font.size': pop.fontsize})
    plt.rcParams.update({"legend.title_fontsize": pop.legend_fontsize})

    fig, ax = plt.subplots(ncols=1)
    fig.set_size_inches(15.5, 10.5)
    ax.scatter(sigma_exponents, distances)

    ax.set_xlabel('Number of initial Planetesimals')

    ax.set_ylabel('Distance')
    ax.legend()
    fig.suptitle('Distance Metric to Solar System')
    fig.savefig(path.join(pop.PLOT, 'distances_exponent.png'))
    plt.close(fig)

    fig, ax = plt.subplots(ncols=1)
    fig.set_size_inches(15.5, 10.5)
    ax.scatter(total_masses, distances)

    ax.set_xlabel('Number of initial Planetesimals')

    ax.set_ylabel('Distance')
    ax.legend()
    fig.suptitle('Distance Metric to Solar System')
    fig.savefig(path.join(pop.PLOT, 'distances_mass.png'))
    plt.close(fig)

    cmap = 'viridis_r'
    cmin = min(distances)
    cmax = max(distances)

    norm = colors.Normalize(cmin, cmax)

    fig, ax = plt.subplots(figsize=pop.figsize)
    ax.scatter(sigma_exponents, total_masses, c=distances, cmap=cmap, norm=norm, s=12)
    x_labels = ax.get_xticklabels(list(set(sigma_exponents)))
    plt.setp(x_labels, horizontalalignment='center')
    ax.set(xlabel='Sigma Exponent', ylabel=r'Total Mass', xticks=sigma_exponents)
    fig.colorbar(cm.ScalarMappable(cmap=cmap, norm=norm), orientation='vertical', label=r'Distance',
                 ax=ax)
    ax.set_xscale('log')
    ax.set_yscale('log')
    if pop.plot_config == 'presentation':
        ax.set(title=r'Distances between Systems')
    save_name = 'scatter_distances'
    if a_up_lim < 30 and m_low_lim > 0:
        save_name += '_lim'
    fig.savefig(path.join(pop.PLOT, save_name + '.png'), transparent=False, dpi=pop.dpi, bbox_inches="tight")
    fig.savefig(path.join(pop.PLOT, 'scatter_distances.png'), transparent=False, dpi=pop.dpi, bbox_inches="tight")
    plt.close(fig)
