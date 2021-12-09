import sys
import os.path as path
from os import getcwd
from os import makedirs
from os import system
from create_options import write_option_file
import disk as dsk


def setup(NAME, N_sims):
    N_SIMS = int(N_sims)

    # Current Working Directory
    CWD = getcwd()

    # Path of the Executable
    EXECUTABLE = path.join(CWD, 'CodeBase/3DPopSyn')

    # Directory of the Synthesis Run
    RUN = path.join(CWD, 'SynthesisRuns', NAME)

    # History File Path
    HISTORY = path.join(RUN,'history.txt')


    def dir_structure(index):
        SYSTEM = path.join(RUN, 'system_' + str(index))
        INPUT = path.join(SYSTEM, 'inputs')
        OUTPUT = path.join(SYSTEM, 'outputs')
        makedirs(INPUT, exist_ok=True)
        makedirs(OUTPUT, exist_ok=True)
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
        write_option_file(INPUT)

    # Creating History File
    with open(HISTORY, 'w') as history:
        pass

    #  INPUT AND OUTPUT WITH $LSB_JOBINDEX
    JI_SYSTEM = path.join(RUN, 'system_\$LSB_JOBINDEX')
    JI_INPUT = path.join(JI_SYSTEM, 'inputs')
    JI_OUTPUT = path.join(JI_SYSTEM, 'outputs')
    MAXTIME=120
    command = 'bsub -J "{name}[1-{N}]%100" -n 1 -r -W {maxtime}:00 -o {system}/log "{exe} {input} {output} {history}"'.format(
        name=NAME,
        N=str(N_SIMS),
        maxtime=MAXTIME,
        system=JI_SYSTEM,
        exe=EXECUTABLE,
        input=JI_INPUT,
        output=JI_OUTPUT,
        history=HISTORY)
    system(command)
    return command


if __name__ == "__main__":
    Name = sys.argv[1]
    N = sys.argv[2]
    #
    # Name = 'Test'
    # N = 5
    cmd = setup(Name, N)
    print(cmd)
