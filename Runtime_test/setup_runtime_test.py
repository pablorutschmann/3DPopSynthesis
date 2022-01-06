import numpy as np
import sys
import os.path as path
from os import getcwd
from os import makedirs
from os import system
from units import *
from create_options import write_option_file
from disk import disk

sigmas = np.round(np.linspace(-1.5, -0.5, 5), 2)
total_masses = np.round(np.linspace(10, 100, 4), 2) * M_J / M_S

planetesimal_number = np.array([100, 200, 500, 800, 1000])

combinations = np.array(np.meshgrid(planetesimal_number, total_masses, sigmas)).T.reshape(-1,3)

NAME = 'RUNTIME'
N_SIMS = len(combinations[:,0])
N_EMBRYO = 10
RUNTIME = 24
EVOTIME = 10000

# Current Working Directory
CWD = getcwd()

# Path of the Executable
EXECUTABLE = path.join(CWD, 'CodeBase/3DPopSyn')

# Directory of the Synthesis Run
RUN = path.join(CWD, NAME)
makedirs(RUN, exist_ok=True)

# LOG File Path
LOG = path.join(RUN, "log")

# History File Path
HISTORY = path.join(RUN, 'history.txt')

def dir_structure(index):
    SYSTEM = path.join(RUN, 'system_' + str(index))
    INPUT = path.join(SYSTEM, 'inputs')
    OUTPUT = path.join(SYSTEM, 'outputs')

    makedirs(INPUT, exist_ok=True)
    makedirs(OUTPUT, exist_ok=True)
    # Create LOG File
    with open(LOG, 'w') as log:
        pass

    return INPUT, OUTPUT

# Create Disk Object
disk_object = disk()
#Disk Parameters
R_min = 0.3
R_max = 30
N = 1000
spacing = "log"

disk_object.prepare(combinations, R_min, R_max, N, spacing)

for i in range(1, N_SIMS + 1):
    INPUT, OUTPUT = dir_structure(i)
    # Create DiskFile
    SIGMA_COEFF, SIGMA_NORM, TEMP_COEFF, N_PLANETESIMALS, TOTAL_MASS = disk_object.sample(INPUT,i)
    # Create Options File
    write_option_file(INPUT, RUNTIME, EVOTIME, SIGMA_COEFF, SIGMA_NORM, TEMP_COEFF, N_EMBRYO, N_PLANETESIMALS, TOTAL_MASS)

# Creating History File
with open(HISTORY, 'w') as history:
    pass

#  INPUT AND OUTPUT WITH $LSB_JOBINDEX
JI_SYSTEM = path.join(RUN, 'system_\$LSB_JOBINDEX')
JI_INPUT = path.join(JI_SYSTEM, 'inputs')
JI_OUTPUT = path.join(JI_SYSTEM, 'outputs')
JI_LOG = path.join(RUN, 'log')

# Command for Euler Job Array
command = 'bsub -J "{name}[1-{N}]%100" -n 1 -r -W {runtime}:00 -oo {log} "{exe} {input} {output} {history}"'.format(
    name=NAME,
    N=str(N_SIMS),
    runtime=RUNTIME,
    log=JI_LOG,
    exe=EXECUTABLE,
    input=JI_INPUT,
    output=JI_OUTPUT,
    history=HISTORY)
system(command)
print(command)

# /Users/prut/CLionProjects/3DPopSynthesis/Runtime_test/RUNTIME_TEST/system_1/inputs /Users/prut/CLionProjects/3DPopSynthesis/Runtime_test/RUNTIME_TEST/system_1/outputs /Users/prut/CLionProjects/3DPopSynthesis/Runtime_test/RUNTIME_TEST/history.txt