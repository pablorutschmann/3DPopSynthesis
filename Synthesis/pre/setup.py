import sys
import os.path as path
from os import getcwd
from os import makedirs
from os import system
from ..units import *
from .create_options import write_option_file
from .disk import disk


def setup(NAME, N_sims, Runtime, Evotime):
    N_SIMS = int(N_sims)
    RUNTIME = int(Runtime)
    EVOTIME = float(Evotime)

    # Current Working Directory
    CWD = getcwd()

    # Path of the Executable
    EXECUTABLE = path.join(CWD, 'CodeBase/3DPopSyn')

    # Directory of the Synthesis Run
    RUN = path.join(CWD, 'SynthesisRuns', NAME)

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

    #Synthesis Parameters
    #Total Disk Mass Range
    TM_min = 10 * M_J / M_S
    TM_max = 100 * M_J / M_S

    # Surface Density Coefficient Range
    Sigma_min = -1.5
    Sigma_max = -0.5
    # Number of Bins for Sigma
    N_Sigma = 5

    # Number of Embryos and Planetesimals
    N_EMBRYO = 10
    N_PLANETESIMAL = 100

    disk_object.prepare(TM_min, TM_max, R_min, R_max, N, spacing, Sigma_min, Sigma_max, N_Sigma)

    for i in range(1, N_SIMS + 1):
        INPUT, OUTPUT = dir_structure(i)
        # Create DiskFile
        SIGMA_COEFF, SIGMA_NORM, TEMP_COEFF = disk_object.sample(INPUT)
        # Create Options File
        write_option_file(INPUT, RUNTIME, EVOTIME, SIGMA_COEFF, SIGMA_NORM, TEMP_COEFF, N_EMBRYO, N_PLANETESIMAL)

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
    return command


if __name__ == "__main__":
    Name = sys.argv[1]
    N = sys.argv[2]
    RUNTIME = sys.argv[3]
    EVOTIME = sys.argv[4]

    # Name = 'Test'
    # N = 5
    cmd = setup(Name, N, RUNTIME, EVOTIME)
    print(cmd)
