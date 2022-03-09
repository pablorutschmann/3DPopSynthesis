import os.path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm
from matplotlib import colors
import os.path as path
from Synthesis.post.metric import *
from Synthesis.post.metric import earth_distance
from tqdm import tqdm
from importlib import reload

from sklearn.manifold import TSNE

RMC_sol = 89.9
AMD_sol = 0.0018


def histogram_AMD(pop, m_low_lim=0, a_up_lim=30):
    AMDS = []
    NS = []

    for id, sys in tqdm(pop.SIMS.items()):
        if id > 250:
            break
        else:
            AMD, N = sys.get_AMD(m_low_lim, a_up_lim)
            AMDS.append(AMD)
            NS.append(N)

    AMDS = np.array((AMDS))

    plt.rcParams.update({'figure.autolayout': True})
    plt.style.use('seaborn-paper')
    plt.rcParams.update({'font.size': pop.fontsize})

    N_bins = 15
    bins = 10 ** np.linspace(np.log10(min(AMDS)), np.log10(max(AMDS)), N_bins)
    fig, ax = plt.subplots(figsize=pop.figsize)
    # ax.hist(Masses, bins=bins)
    values, base, _ = plt.hist(AMDS, bins=bins, rwidth=0.95)
    print(f'Bin edges: {base[np.argwhere(base >= AMD_sol)[0,0]]}, {base[np.argwhere(base >= AMD_sol)[0,0] + 1]}')
    print(f'Number of systems in that bin: {values[np.argwhere(base >= AMD_sol)[0,0]-1]}')
    print(f'Fraction: {values[np.argwhere(base >= AMD_sol)[0,0]-1]/np.sum(values)}')
    ax.axvline(AMD_sol, color='red', linewidth=1)
    # ax_bis = ax.twinx()
    # values, base = np.histogram(Masses, bins = 10 ** np.linspace(np.log10(min(Masses)), np.log10(max(Masses)), N_bins))

    # values = np.append(0, values)
    # ax_bis.plot(base, np.cumsum(values) / np.cumsum(values)[-1], color='black', linestyle='dashed', markersize=0.1)
    ax.set(xlabel=r'AMD', ylabel=r'Counts')
    # ax_bis.set(ylabel='Cumulative Distribution')
    ax.set_xscale('log')
    # ax_bis.set_xscale('log')
    if pop.plot_config == 'presentation':
        ax.set(title=r'Histrogram of AMD')
    save_name = 'histogram_AMD'
    if a_up_lim < 30 and m_low_lim > 0:
        save_name += '_lim'
    fig.savefig(path.join(pop.PLOT, save_name + '.png'), transparent=False, dpi=pop.dpi, bbox_inches="tight")
    plt.close(fig)


def histogram_RMC(pop, m_low_lim=0, a_up_lim=30):
    RMCS = []
    NS = []

    for id, sys in tqdm(pop.SIMS.items()):
        if id > 250:
            break
        else:
            RMC, N = sys.get_RMC(m_low_lim, a_up_lim)
            RMCS.append(RMC)
            NS.append(N)

    RMCS = np.array((RMCS))

    plt.rcParams.update({'figure.autolayout': True})
    plt.style.use('seaborn-paper')
    plt.rcParams.update({'font.size': pop.fontsize})

    N_bins = 15
    bins = 10 ** np.linspace(np.log10(min(RMCS)), np.log10(max(RMCS)), N_bins)
    fig, ax = plt.subplots(figsize=pop.figsize)
    # ax.hist(Masses, bins=bins)
    values, base, _ = plt.hist(RMCS, bins=bins, rwidth=0.95)
    print(f'Bin edges: {base[np.argwhere(base >= AMD_sol)[0,0]]}, {base[np.argwhere(base >= RMC_sol)[0,0] + 1]}')
    print(f'Number of systems in that bin: {values[np.argwhere(base >= RMC_sol)[0,0]-1]}')
    print(f'Fraction: {values[np.argwhere(base >= RMC_sol)[0,0]-1]/np.sum(values)}')
    ax.axvline(RMC_sol, color='red', linewidth=1)
    # ax_bis = ax.twinx()
    # values, base = np.histogram(Masses, bins = 10 ** np.linspace(np.log10(min(Masses)), np.log10(max(Masses)), N_bins))

    # values = np.append(0, values)
    # ax_bis.plot(base, np.cumsum(values) / np.cumsum(values)[-1], color='black', linestyle='dashed', markersize=0.1)
    ax.set(xlabel=r'RMC', ylabel=r'Counts')
    # ax_bis.set(ylabel='Cumulative Distribution')
    ax.set_xscale('log')
    # ax_bis.set_xscale('log')
    if pop.plot_config == 'presentation':
        ax.set(title=r'Histrogram of RMC')
    save_name = 'histogram_RMC'
    if a_up_lim < 30 and m_low_lim > 0:
        save_name += '_lim'
    fig.savefig(path.join(pop.PLOT, save_name + '.png'), transparent=False, dpi=pop.dpi, bbox_inches="tight")
    plt.close(fig)


def scatter_AMD_RMC(pop, m_low_lim=0, a_up_lim=30):
    AMDS = []
    RMCS = []
    NS = []
    SIGMAS = []
    TOTALMASSES = []

    for id, sys in tqdm(pop.SIMS.items()):
        if id > 250:
            break
        else:
            AMD, N = sys.get_AMD(m_low_lim, a_up_lim)
            AMDS.append(AMD)
            RMC, N = sys.get_RMC(m_low_lim, a_up_lim)
            RMCS.append(RMC)
            NS.append(N)
            SIGMAS.append(sys.Sigma_Exponent)
            TOTALMASSES.append(sys.Total_Mass * M_S / M_J)

    AMDS = np.array((AMDS))
    RMCS = np.array(RMCS)
    SIGMAS = np.array(SIGMAS)
    TOTALMASSES = np.array(TOTALMASSES)
    sigma_AMD = np.std(AMDS)
    sigma_RMC = np.std(RMCS)
    print(f'{sigma_AMD=}')
    print(f'{sigma_RMC=}')

    def dist_to_solar(rmc, amd):
        return np.sqrt((rmc - RMC_sol) ** 2 / sigma_RMC ** 2 + (amd - AMD_sol) ** 2 / sigma_AMD ** 2)

    distances = dist_to_solar(RMCS, AMDS)
    print(distances)

    cmap = pop.cmap_standart
    cmin = min(distances)
    cmax = max(distances)

    norm = colors.Normalize(cmin, cmax)

    plt.rcParams.update({'figure.autolayout': True})
    plt.style.use('seaborn-paper')
    plt.rcParams.update({'font.size': pop.fontsize})

    fig, ax = plt.subplots(figsize=pop.figsize)
    #ax.scatter(AMDS/sigma_AMD, RMCS/sigma_RMC, c=distances, cmap=cmap, norm=norm, s=12)
    #ax.scatter(SIGMAS, TOTALMASSES, c=distances[np.argmin(distances)])
    ax.scatter(SIGMAS, TOTALMASSES, c=distances, cmap=cmap, norm=norm, s=12)
    x_labels = ax.get_xticklabels()
    plt.setp(x_labels, horizontalalignment='center')
    ax.set(xlabel='Power Law Exponent', ylabel=r'Total Disk Mass [$M_J$]')
    fig.colorbar(cm.ScalarMappable(cmap=cmap, norm=norm), orientation='vertical', label=r'Distance',
                 ax=ax)
    # ax.set_xscale('log')
    ax.set_yscale('log')
    if pop.plot_config == 'presentation':
        ax.set(title=r'Radial Mass Concentration and Angular Momentum Deficit')

    save_name = 'scatter_para_dist_AMD_RMC'
    if a_up_lim < 30 and m_low_lim > 0:
        save_name += '_lim'
    fig.savefig(path.join(pop.PLOT, save_name + '.png'), transparent=False, dpi=pop.dpi, bbox_inches="tight")
    plt.close(fig)


def distances(pop, m_low_lim=0, a_up_lim=30):
    systems = []
    total_masses = []
    sigma_exponents = []
    numbers = []
    reference = []


    for id, sim in pop.SIMS.items():
        if id > 250:
            break
        total_masses.append(sim.Total_Mass)
        sigma_exponents.append(sim.Sigma_Exponent)
        reference.append(sim.Sigma_Norm / (R_S / au)**sim.Sigma_Exponent / denstos * pow(au/R_S, sim.Sigma_Exponent) / denstos)

        zipped = zip(sim.snaps[sim.N_snaps - 1].satellites['M'].values,
                     sim.snaps[sim.N_snaps - 1].satellites['a'].values)

        # Remove under threshold Masses
        filtered = [(m * M_S / M_E, a * R_S / au) for (m, a) in zipped if
                    m >= m_low_lim * M_E / M_S and a <= a_up_lim * au / R_S]
        systems.append(filtered)
        numbers.append(len(filtered))
    distances = []
    num=len(systems)
    # print(systems)
    # func = lambda x: distance(x, terrestrial)
    # distances = list(map(func, systems))
    try:
        distances = np.load(file=path.join(pop.PLOT,f'dist_terr_log_M_{num}.npy'))

    except:
        for sys in tqdm(systems):
            distances.append(distance(sys, terrestrial))
        np.save(arr=distances,file=path.join(pop.PLOT,f'dist_terr_{num}.npy'))
    print(distances)
    plt.rcParams.update({'figure.autolayout': True})
    plt.style.use('seaborn-paper')
    plt.rcParams.update({'font.size': pop.fontsize})
    plt.rcParams.update({"legend.title_fontsize": pop.legend_fontsize})

    cmap = pop.cmap_standart
    cmin = min(distances)
    cmax = max(distances)

    norm = colors.LogNorm(cmin, cmax)

    fig, ax = plt.subplots(figsize=pop.figsize)
    ax.scatter(sigma_exponents, total_masses, c=distances, cmap=cmap, norm=norm, s=12)
    x_labels = ax.get_xticklabels(list(set(sigma_exponents)))
    plt.setp(x_labels, horizontalalignment='center')
    ax.set(xlabel='Sigma Exponent', ylabel=r'Total Disk Mass [$M_{\odot}$]', xticks=sigma_exponents)
    ax2 = ax.twinx()
    mn, mx = ax.get_ylim()
    ax2.set_ylim(M_S / M_J * mn, M_S / M_J * mx)
    ax2.set_ylabel('Total Disk Mass [$M_{J}$]')
    fig.colorbar(cm.ScalarMappable(cmap=cmap, norm=norm), orientation='vertical', label=r'Distance',
                 ax=ax2, pad=0.12)

    if pop.plot_config == 'presentation':
        ax.set(title=r'Distances between Systems')
    save_name = 'scatter_distances_log_M'
    if a_up_lim < 30 and m_low_lim > 0:
        save_name += '_lim'
    fig.savefig(path.join(pop.PLOT, save_name + '.png'), transparent=False, dpi=pop.dpi, bbox_inches="tight")
    #fig.savefig(path.join(pop.PLOT, 'scatter_distances.png'), transparent=False, dpi=pop.dpi, bbox_inches="tight")
    plt.close(fig)

    plt.rcParams.update({'figure.autolayout': True})
    plt.style.use('seaborn-paper')
    plt.rcParams.update({'font.size': pop.fontsize})
    plt.rcParams.update({"legend.title_fontsize": pop.legend_fontsize})
    cmap = pop.cmap_standart
    cmin = min(numbers)
    cmax = max(numbers)

    norm = colors.Normalize(cmin, cmax)

    fig, ax = plt.subplots(figsize=pop.figsize)
    ax.scatter(reference, distances, c=numbers, norm=norm, s=12)
    # x_labels = ax.get_xticklabels(list(set(reference)))
    # plt.setp(x_labels, horizontalalignment='center')
    ax.set(xlabel='Reference Value [$\mathrm{g cm^{-2}}$]', ylabel=r'Distance')
    fig.colorbar(cm.ScalarMappable(cmap=cmap, norm=norm), orientation='vertical', label=r'Number of Planets',
                 ax=ax)
    ax.set_yscale('log')
    ax.set_xscale('log')
    if pop.plot_config == 'presentation':
        ax.set(title=r'Distances between Systems')
    save_name = 'scatter_distances_log_M_reference'
    if a_up_lim < 30 and m_low_lim > 0:
        save_name += '_lim'
    fig.savefig(path.join(pop.PLOT, save_name + '.png'), transparent=False, dpi=pop.dpi, bbox_inches="tight")
    #fig.savefig(path.join(pop.PLOT, 'scatter_distances.png'), transparent=False, dpi=pop.dpi, bbox_inches="tight")
    plt.close(fig)

def distances_reference(pop, m_low_lim=0, a_up_lim=30):
    systems = []
    total_masses = []
    sigma_exponents = []
    reference = []
    numbers = []


    for id, sim in tqdm(pop.SIMS.items()):
        if id > 250:
            break
        total_masses.append(sim.Total_Mass)
        sigma_exponents.append(sim.Sigma_Exponent)
        reference.append(sim.Sigma_Norm / (R_S / au)**sim.Sigma_Exponent / denstos * pow(au/R_S, sim.Sigma_Exponent) / denstos)

        zipped = zip(sim.snaps[sim.N_snaps - 1].satellites['M'].values * M_S / M_E,
                     sim.snaps[sim.N_snaps - 1].satellites['a'].values  * R_S / au)

        # Remove under threshold Masses
        filtered = [(m, a) for (m, a) in zipped if
                    m >= m_low_lim and a <= a_up_lim]

        systems.append(filtered)
        numbers.append(len(filtered))
    print(numbers)
    # print(systems)
    # func = lambda x: distance(x, terrestrial)
    # distances = list(map(func, systems))
    num=len(systems)
    try:
        distances = np.load(file=path.join(pop.PLOT,f'dist_terr_log_M_{num}.npy'))
    except:
        distances = []
        for sys in tqdm(systems):
            distances.append(distance(sys, terrestrial))
        np.save(arr=distances,file=path.join(pop.PLOT,f'dist_terr_{num}'))
    print(distances)

    plt.rcParams.update({'figure.autolayout': True})
    plt.style.use('seaborn-paper')
    plt.rcParams.update({'font.size': pop.fontsize})
    plt.rcParams.update({"legend.title_fontsize": pop.legend_fontsize})

    cmap = pop.cmap_standart
    cmin = min(numbers)
    cmax = max(numbers)

    norm = colors.Normalize(cmin, cmax)

    fig, ax = plt.subplots(figsize=pop.figsize)
    ax.scatter(reference, distances, c=numbers, norm=norm, s=12)
    # x_labels = ax.get_xticklabels(list(set(reference)))
    # plt.setp(x_labels, horizontalalignment='center')
    ax.set(xlabel='Reference Value [$\mathrm{g cm^{-2}}$]', ylabel=r'Distance')
    fig.colorbar(cm.ScalarMappable(cmap=cmap, norm=norm), orientation='vertical', label=r'Number of Planets',
                 ax=ax)
    # ax.set_yscale('log')
    ax.set_xscale('log')
    if pop.plot_config == 'presentation':
        ax.set(title=r'Distances between Systems')
    save_name = 'scatter_distances_log_M_reference'
    if a_up_lim < 30 and m_low_lim > 0:
        save_name += '_lim'
    fig.savefig(path.join(pop.PLOT, save_name + '.png'), transparent=False, dpi=pop.dpi, bbox_inches="tight")
    #fig.savefig(path.join(pop.PLOT, 'scatter_distances.png'), transparent=False, dpi=pop.dpi, bbox_inches="tight")
    plt.close(fig)


def tsne(pop, m_low_lim=0, a_up_lim=30):
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
        if id > 3:
            break
        TotalMasses.append(sim.Total_Mass)
        SigmaCoeffs.append(sim.Sigma_Exponent)

    for id, sim in tqdm(pop.SIMS.items()):
        if id > 3:
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
        #print(system)
        sys_matrix.append(system)
        #Reference.append(distance(system, terrestrial))
        Reference.append(sim.Sigma_Norm / (R_S / au)**sim.Sigma_Exponent / denstos * pow(au/R_S, sim.Sigma_Exponent) / denstos)
        # Reference.append(cum_mass)
    Numbers = Reference

    terr = [(m/M_E, a/au) for (m,a) in terrestrial]
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


    def get_distances(n,arr):
        out = np.zeros(n)
        sys = arr[0]
        others = arr[1:]
        def dist_arr(other):
            return distance(sys,other)

        dist_func = np.vectorize(dist_arr)

        distances_array = dist_func(others)

        #out = np.pad(distances_array, (0,0) , 'constant', constant_values=(n-len(others), 0))
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


    try:
        matrix_metric = np.load(path.join(pop.PLOT,f'metric_matrix_log_M_{num}.npy'))
        print('Loaded from pickle')
    except:
        rows = np.array([np.array(sys_matrix[i:], dtype='object') for i in range(num-1)], dtype='object')
        print(rows)

        upper = np.array(list(map(lambda x :get_distances(arr=x,n=num),rows)))
        print(upper)
        upper = np.vstack([upper,np.zeros(num)])

        upper = np.reshape(upper,(num,num))
        matrix_metric = upper + upper.T

        np.save(path.join(pop.PLOT,f'metric_matrix_log_M_{num}'), matrix_metric)
        print('saved to Pickle')
        print(matrix_metric)

    # AMDS = np.array((AMDS))
    # RMCS = np.array(RMCS)
    # sigma_AMD = np.std(AMDS)
    # sigma_RMC = np.std(R  CS)
    # print(f'{sigma_AMD=}')
    # print(f'{sigma_RMC=}')
    #
    # def dist_to_solar(rmc, amd):
    #     return np.sqrt((rmc - RMC_sol) ** 2 / sigma_RMC ** 2 + (amd - AMD_sol) ** 2 / sigma_AMD ** 2)
    #
    # distances = dist_to_solar(RMCS, AMDS)
    # print(f'{distances=}')

    max_len = 0
    for sys in systems:
        if len(sys) > max_len:
            max_len = len(sys)

    data = []
    for sys in systems:
        sys = np.pad(sys, pad_width=(0, max_len - len(sys)), mode="constant", constant_values=0.0)
        data.append(sys)

    terrestrial_data_point = []
    for item in terrestrial:
        terrestrial_data_point.extend([item[0] / M_E, item[1] / au])

    terrestrial_data_point = np.pad(terrestrial_data_point, pad_width=(0, max_len - len(terrestrial_data_point)),
                                    mode="constant", constant_values=0.0)
    # print(data)

    data = np.vstack(data)
    # print(data)
    print(f'{data.shape=}')

    all_data = np.vstack([data, terrestrial_data_point])

    model = TSNE(2, perplexity=15, metric='precomputed')

    X_embedded = model.fit_transform(matrix_metric)

    systems_embedded = X_embedded[:-1, :]
    terrestrial_embedded = X_embedded[-1, :]

    def dist_to_earth_embed(point):
        return np.sqrt((point[0] - terrestrial_embedded[0])**2 + (point[1] - terrestrial_embedded[1])**2)

    dist_embed = np.apply_along_axis(arr=systems_embedded, func1d=dist_to_earth_embed, axis=1)
    Numbers = dist_embed
    print(dist_embed)
    print(f'{systems_embedded.shape=}')

    plt.rcParams.update({'figure.autolayout': True})
    plt.style.use('seaborn-paper')
    plt.rcParams.update({'font.size': pop.fontsize})
    plt.rcParams.update({"legend.title_fontsize": pop.legend_fontsize})

    Numbers = Reference
    cmap = pop.cmap_standart
    cmin = min(Numbers)
    cmax = max(Numbers)

    norm = colors.Normalize(cmin, cmax)

    fig, ax = plt.subplots(figsize=pop.figsize)
    ax.scatter(systems_embedded[:, 0], systems_embedded[:, 1], c=Numbers, cmap=cmap, norm=norm, s=12)
    ax.scatter(terrestrial_embedded[0], terrestrial_embedded[1], c='red')
    x_labels = ax.get_xticklabels()
    plt.setp(x_labels, horizontalalignment='center')
    ax.set_xlabel('Axis 1')
    ax.set_ylabel('Axis 2')
    fig.colorbar(cm.ScalarMappable(cmap=cmap, norm=norm), orientation='vertical', label=r'AMD,RMC Distance',
                 ax=ax)
    # ax.set_xscale('log')
    # ax.set_yscale('log')
    if pop.plot_config == 'presentation':
        ax.set(title=r'Radial Mass Concentration and Angular Momentum Deficit')

    save_name = 'scatter_tsne'
    if a_up_lim < 30 and m_low_lim > 0:
        save_name += '_lim'
    fig.savefig(path.join(pop.PLOT, save_name + '.png'), transparent=False, dpi=pop.dpi, bbox_inches="tight")
    plt.close(fig)

    plt.rcParams.update({'figure.autolayout': True})
    plt.style.use('seaborn-paper')
    plt.rcParams.update({'font.size': pop.fontsize})

    cmap = pop.cmap_standart
    cmin = min(dist_embed)
    cmax = max(dist_embed)

    norm = colors.LogNorm(cmin, cmax)

    fig, ax = plt.subplots(figsize=pop.figsize)
    ax.scatter(SigmaCoeffs, TotalMasses, c=dist_embed, cmap=cmap, norm=norm, s=12)
    x_labels = ax.get_xticklabels()
    plt.setp(x_labels, horizontalalignment='center')
    ax.set(xlabel='Surface Density Power Law Exponent', ylabel=r'Total Disk Mass [$M_{\odot}$]', xticks=SigmaCoeffs)
    ax2 = ax.twinx()
    mn, mx = ax.get_ylim()
    ax2.set_ylim(M_S / M_J * mn, M_S / M_J * mx)
    ax2.set_ylabel('Total Disk Mass [$M_{J}$]')
    fig.colorbar(cm.ScalarMappable(cmap=cmap, norm=norm), orientation='vertical',
                 label=r'Reference Value at $1 \mathrm{au}$ [$\mathrm{g}\mathrm{cm}^{-2}$]', ax=ax2, pad=0.12)
    # ax.set_yscale('log')
    if pop.plot_config == 'presentation':
        ax.set(title=r'Synthesis Parameters')

    fig.savefig(path.join(pop.PLOT, 'scatter_parameters_embed_dist.png'), transparent=False, dpi=pop.dpi,
                bbox_inches="tight")
    plt.close(fig)

    #
    # fig, ax = plt.subplots(ncols=1)
    # fig.set_size_inches(pop.figsize)
    # ax.scatter(systems_embedded[:, 0], systems_embedded[:, 1], cmap=cmap)
    # ax.scatter(terrestrial_embedded[0], terrestrial_embedded[1], color='r')
    #
    # plt.close(fig)


def distance_earth(pop, m_low_lim=0, a_up_lim=30):
    systems = []
    total_masses = []
    sigma_exponents = []

    for id, sim in tqdm(pop.SIMS.items()):
        if id > 50:
            break
        total_masses.append(sim.Total_Mass)
        sigma_exponents.append(sim.Sigma_Exponent)

        zipped = zip(sim.snaps[sim.N_snaps - 1].satellites['M'].values,
                     sim.snaps[sim.N_snaps - 1].satellites['a'].values)

        # Remove under threshold Masses
        filtered = [(m * M_S / M_E, a * R_S / au) for (m, a) in zipped if
                    m >= m_low_lim * M_E / M_S and a <= a_up_lim * au / R_S]
        systems.append(filtered)
    distances = []
    # print(systems)
    # func = lambda x: distance(x, terrestrial)
    # distances = list(map(func, systems))
    def get_min_dist(sys):
        distances = []
        for planet in enumerate(sys):
            print(planet[1])
            distances.append(earth_distance([planet[1]],[(1,1)]))

    for sys in tqdm(systems):
        distances.append(distance(sys, [(M_E,au)]))

    plt.rcParams.update({'figure.autolayout': True})
    plt.style.use('seaborn-paper')
    plt.rcParams.update({'font.size': pop.fontsize})
    plt.rcParams.update({"legend.title_fontsize": pop.legend_fontsize})

    cmap = 'viridis_r'
    cmin = min(distances)
    cmax = max(distances)

    norm = colors.Normalize(cmin, cmax)

    fig, ax = plt.subplots(figsize=pop.figsize)
    ax.scatter(sigma_exponents, total_masses, c=distances, cmap=cmap, norm=norm, s=12)
    x_labels = ax.get_xticklabels(list(set(sigma_exponents)))
    plt.setp(x_labels, horizontalalignment='center')
    ax.set(xlabel='Sigma Exponent', ylabel=r'Total Disk Mass [$M_{\odot}]', xticks=sigma_exponents)
    fig.colorbar(cm.ScalarMappable(cmap=cmap, norm=norm), orientation='vertical', label=r'Distance',
                 ax=ax)
    #ax.set_yscale('log')
    if pop.plot_config == 'presentation':
        ax.set(title=r'Distances between Systems')
    save_name = 'scatter_metric_earth'
    if a_up_lim < 30 and m_low_lim > 0:
        save_name += '_lim'
    fig.savefig(path.join(pop.PLOT, save_name + '.png'), transparent=False, dpi=pop.dpi, bbox_inches="tight")

    plt.close(fig)

def distance_earth2(pop, m_low_lim=0, a_up_lim=30):
    systems = []
    total_masses = []
    sigma_exponents = []
    reference = []


    for id, sim in tqdm(pop.SIMS.items()):
        if id > 242:
            break
        total_masses.append(sim.Total_Mass)
        sigma_exponents.append(sim.Sigma_Exponent)
        reference.append(sim.Sigma_Norm / (R_S / au)**sim.Sigma_Exponent / denstos * pow(au/R_S, sim.Sigma_Exponent) / denstos)

        WM = np.array(sim.snaps[sim.N_snaps - 1].satellites['WM'].values)
        SWM = np.array(sim.snaps[sim.N_snaps - 1].satellites['SWM'].values)
        Masses = np.array(sim.snaps[sim.N_snaps - 1].satellites['M'].values)
        Orb_Dist = np.array(sim.snaps[sim.N_snaps - 1].satellites['a'].values)
        TWMF = (WM + SWM) / Masses

        zipped = zip(Masses, Orb_Dist, TWMF)

        # Remove under threshold Masses
        # filtered = [(m * M_S / M_E, a * R_S / au, wmf) for (m, a, wmf) in zipped if
        #             m >=  m_low_lim * M_E / M_S and a <= a_up_lim * au / R_S and wmf > 0]
        filtered = [(m * M_S / M_E, a * R_S / au, wmf) for (m, a, wmf) in zipped if
                    m >=  0.1 * M_E / M_S and m <= 5 * M_E / M_S and a <= a_up_lim * au / R_S and wmf > 0]
        systems.append(filtered)
    distances = []
    num = len(systems)
    # print(systems)
    # func = lambda x: distance(x, terrestrial)
    # distances = list(map(func, systems))
    #print(earth_distance([(0.6,1.2,0.002)], [earth_wmf]))
    # def get_min_dist(sys):
    #     distances = []
    #     for planet in enumerate(sys):
    #         print(planet[1])
    #         distances.append(earth_distance([planet[1]],[earth_wmf]))
    #
    #     return np.min(distances), sys[np.argmin(distances)]
    # planets = []

        # dist, planet = get_min_dist(sys)
        # distances.append(dist)
        # planets.append(planet)

    try:
        distances = np.load(path.join(pop.PLOT,f'dist_earth_min_{num}.npy'))
        distances_all = np.load(path.join(pop.PLOT,f'dist_earth_{num}.npy'))
        planets = np.load(path.join(pop.PLOT,f'planets_earth_min_{num}.npy'))
        print('Loaded from pickle')
    except:
        distances = []
        for sys in tqdm(systems):
            distances.append(earth_distance(sys, [earth_wmf]))

        np.save(path.join(pop.PLOT,f'dist_earth_{num}'), distances)
        print('saved to Pickle')
        print(distances)

    plt.rcParams.update({'figure.autolayout': True})
    plt.style.use('seaborn-paper')
    plt.rcParams.update({'font.size': pop.fontsize})
    plt.rcParams.update({"legend.title_fontsize": pop.legend_fontsize})

    # for i,item in enumerate(distances):
    #     if item < 0:
    #         distances[i] = max(distances)*0.1


    cmap = pop.cmap_standart
    cmin = min(distances_all)
    cmax = max(distances_all)

    norm = colors.Normalize(cmin, cmax)

    fig, ax = plt.subplots(figsize=pop.figsize)
    ax.scatter(sigma_exponents, total_masses, c=distances_all, cmap=cmap, norm=norm, s=12)
    x_labels = ax.get_xticklabels(list(set(sigma_exponents)))
    plt.setp(x_labels, horizontalalignment='center')
    ax2 = ax.twinx()
    mn, mx = ax.get_ylim()
    ax2.set_ylim(M_S / M_J * mn, M_S / M_J * mx)
    ax2.set_ylabel('Total Disk Mass [$M_{J}$]')
    ax.set(xlabel='Sigma Exponent', ylabel=r'Total Disk Mass [$M_{\odot}$]', xticks=sigma_exponents)
    fig.colorbar(cm.ScalarMappable(cmap=cmap, norm=norm), orientation='vertical', label=r'Minimal Distance',
                 ax=ax2, pad=0.12)
    #ax.set_yscale('log')
    if pop.plot_config == 'presentation':
        ax.set(title=r'Distances between Systems')
    save_name = 'scatter_metric_earth'
    if a_up_lim < 30 and m_low_lim > 0:
        save_name += '_lim'
    fig.savefig(path.join(pop.PLOT, save_name + '.png'), transparent=False, dpi=pop.dpi, bbox_inches="tight")
    plt.close(fig)


    wmfs = [w for (m,a,w) in planets]

    cmap = pop.cmap_standart
    cmin = min(wmfs)
    cmax = max(wmfs)

    norm = colors.Normalize(cmin, cmax)

    fig, ax = plt.subplots(figsize=pop.figsize)
    ax.scatter(reference, distances, c=wmfs, norm=norm, s=12)
    # x_labels = ax.get_xticklabels(list(set(reference)))
    # plt.setp(x_labels, horizontalalignment='center')
    ax.set(xlabel='Reference Value [$\mathrm{g cm^{-2}}$]', ylabel=r'Minimal Distance')
    fig.colorbar(cm.ScalarMappable(cmap=cmap, norm=norm), orientation='vertical', label=r'Total Water Mass Fraction',
                 ax=ax)
    #ax.set_yscale('log')
    ax.set_xscale('log')
    if pop.plot_config == 'presentation':
        ax.set(title=r'Distances between Systems')
    save_name = 'scatter_distances_earth_reference'
    if a_up_lim < 30 and m_low_lim > 0:
        save_name += '_lim'
    fig.savefig(path.join(pop.PLOT, save_name + '.png'), transparent=False, dpi=pop.dpi, bbox_inches="tight")
    #fig.savefig(path.join(pop.PLOT, 'scatter_distances.png'), transparent=False, dpi=pop.dpi, bbox_inches="tight")
    plt.close(fig)


