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
    EXECUTABLE = path.join(CWD, '3DPopSynth')

    # Directory 0f the Synthesis Run
    RUN = path.join(CWD, 'SynthesisRuns', NAME)

    print(RUN)

    def dir_structure(index):
        SYSTEM = path.join(RUN, 'system_' + str(index))
        INPUT = path.join(SYSTEM, 'inputs/')
        OUTPUT = path.join(SYSTEM, 'outputs/')
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

    for i in range(N_SIMS):
        INPUT, OUTPUT = dir_structure(i)
        # Create DiskFile
        disk.sample(INPUT)
        # Create Options File
        write_option_file(INPUT)

    #  INPUT AND OUTPUT WITH $LSB_JOBINDEX
    JI_SYSTEM = path.join(RUN, 'system_\$LSB_JOBINDEX')
    JI_INPUT = path.join(JI_SYSTEM, 'inputs/')
    JI_OUTPUT = path.join(JI_SYSTEM, 'outputs/')
    JI_HISTORY = path.join(JI_SYSTEM, 'history.txt')

    command = 'bsub -J "{name}[1-{N}]%100" -r -W {maxtime}:00 -o {run}/log "{exe} {input} {output} {history}"'.format(
        name=NAME,
        N=str(N_SIMS),
        maxtime=MAXTIME,
        run=RUN,
        exe=EXECUTABLE,
        input=JI_INPUT,
        output=JI_OUTPUT,
        history=JI_HISTORY)
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
