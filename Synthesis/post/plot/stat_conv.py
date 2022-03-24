import numpy as np
import random
from Synthesis.units import *
import os.path as path
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from tqdm import tqdm

def func(x, a, b):
    return a * np.exp(-b * x)


def statistical_convergence(POP):
    sample_sizes = np.arange(10,POP.NSIMS+1,1)
    print(sample_sizes)
    def sample():
        mass_std = []
        mass_mean = []
        orb_dist_std = []
        orb_dist_mean = []
        tm_mean = []
        tm_std = []

        for sub in sample_sizes:
            keys = POP.SIMS.keys()
            ids = random.sample(keys, sub)
            masses = []
            orb_dist = []
            tot_masses = []
            for id in ids:
                mass = list(POP.SIMS[id].snaps[POP.SIMS[id].N_snaps - 1].satellites['M'].values * M_S / M_E)
                masses += mass
                tot_masses.append(np.sum(mass))
                orb_dist += list(POP.SIMS[id].snaps[POP.SIMS[id].N_snaps - 1].satellites['a'].values * R_S / au)
            mass_std.append(np.std(masses))
            mass_mean.append((np.mean(masses)))

            orb_dist_std.append(np.std(orb_dist))
            orb_dist_mean.append((np.mean(orb_dist)))

            tm_std.append(np.std(tot_masses))
            tm_mean.append((np.mean(tot_masses)))
        return tm_mean, tm_std

    tm_mean, tm_std = sample()
    # print(f'STD Mass: {mass_std}')
    # print(f'MEAN Mass: {mass_mean}')
    # print(f'STD Orb Dist: {orb_dist_std}')
    # print(f'MEAN Orb Dist: {orb_dist_mean}')

    def get_ratios(arr):
        out = [1]
        for i,item in enumerate(arr[:-1]):
            out.append((item - arr[i+1]) / item)
        return np.abs(np.array(out))

    def movingaverage(interval, window_size):
        window= np.ones(int(window_size))/float(window_size)
        return np.convolve(interval, window, 'same')

    def converges_at(a, b, limit=0.001):
        return - np.log(limit / a) / b

    plt.rcParams.update({'figure.autolayout': True})
    plt.style.use('seaborn-paper')
    plt.rcParams.update({'font.size': POP.fontsize})
    plt.rcParams.update({"legend.title_fontsize": POP.legend_fontsize})

    fig, ax = plt.subplots(figsize=POP.figsize)

    popt_mean, pcov_mean = curve_fit(func, sample_sizes, get_ratios(tm_mean), bounds=(0, [3., 0.5]))
    print(func(242,*popt_mean))
    print(f'Total Mass converges at {converges_at(*popt_mean)} samples.')
    ax.plot(sample_sizes, get_ratios(tm_mean), linewidth=0.7, markersize=2, label=r'$\Delta_M$ Mean Total Mass')
    ax.plot(sample_sizes, movingaverage(get_ratios(tm_mean), 30), 'r', label='Moving Average')

    ax.plot(sample_sizes, func(sample_sizes, *popt_mean), 'g--', label='Exponential Decay Fit')

    # ax.plot(sample_sizes, get_ratios(tm_std), linestyle='dashed', linewidth=0.7,
    #         markersize=6,
    #         marker='o', label='STD Total Mass')
    # ax.plot(sample_sizes, movingaverage(get_ratios(tm_std), 20), label='Moving Average')
    # ax.set_yscale('log')
    ax.set(xlabel='Sample Size', ylabel=r'Rate of Change')
    ax.legend(loc='best')
    save_name = 'stat_conv'
    fig.savefig(path.join(POP.PLOT, save_name + '.png'), transparent=False, dpi=POP.dpi,
                bbox_inches="tight")
    plt.close(fig)

    N_Samples = 300
    mean_conv_numbers = []
    std_conv_numbers = []
    tft_rates = []

    for i in tqdm(range(N_Samples)):
        tm_mean, tm_std = sample()
        popt_mean, pcov_mean = curve_fit(func, sample_sizes, get_ratios(tm_mean), bounds=(0, [1.5, 0.5]))
        popt_std, pcov_std = curve_fit(func, sample_sizes, get_ratios(tm_std), bounds=(0, [1.5, 0.5]))
        mean_conv_numbers.append(converges_at(*popt_mean))
        tft_rates.append(func(242,*popt_mean))
        std_conv_numbers.append(converges_at(*popt_std))

    mean_conv_mean = np.mean(mean_conv_numbers)
    mean_conv_std = np.std(mean_conv_numbers)
    print(f'Mean Convergences at {mean_conv_mean} +- {mean_conv_std / np.sqrt(N_Samples)} samples')

    std_conv_mean = np.mean(std_conv_numbers)
    std_conv_std = np.std(std_conv_numbers)
    print(f'Std Convergences at {std_conv_mean} +- {std_conv_std / np.sqrt(N_Samples)} samples')

    print(f'Reached Convergence Rate {np.mean(tft_rates)} +- {np.std(tft_rates)}')