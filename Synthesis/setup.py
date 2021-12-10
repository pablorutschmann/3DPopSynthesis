import sys
import os.path as path
from os import getcwd
from os import makedirs
from os import system
from create_options import write_option_file
import disk as dsk


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
    LOG = path.join(RUN,"log")

    # History File Path
    HISTORY = path.join(RUN,'history.txt')


    def dir_structure(index):
        SYSTEM = path.join(RUN, 'system_' + str(index))
        INPUT = path.join(SYSTEM, 'inputs')
        OUTPUT = path.join(SYSTEM, 'outputs')
        makedirs(INPUT, exist_ok=True)
        makedirs(OUTPUT, exist_ok=True)
        # LOG File Path
        LOG = path.join(RUN,"log")
        # Create LOG File
        with open(LOG, 'w') as log:
            pass

        return INPUT, OUTPUT

    # Create Disk Object
    disk = dsk.disk()
    spacing = "log"
    R_min = 0.5
    R_max = 30
    N = 1000
    disk.prepare(spacing, R_min, R_max, N)

    for i in range(1,N_SIMS+1):
        INPUT, OUTPUT = dir_structure(i)
        # Create DiskFile
        disk.sample(INPUT)
        # Create Options File
        write_option_file(INPUT, RUNTIME, EVOTIME)

    # Creating History File
    with open(HISTORY, 'w') as history:
        pass

    #  INPUT AND OUTPUT WITH $LSB_JOBINDEX
    JI_SYSTEM = path.join(RUN, 'system_\$LSB_JOBINDEX')
    JI_INPUT = path.join(JI_SYSTEM, 'inputs')
    JI_OUTPUT = path.join(JI_SYSTEM, 'outputs')
    JI_LOG =  path.join(RUN, 'log')

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
