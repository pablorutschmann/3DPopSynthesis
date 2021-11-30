import sys
import os.path as path
from os import getcwd
from os import makedirs

def setup(NAME, N_SIMS):
    runname = NAME

    N_sims = int(N_SIMS)

    CWD = getcwd()
    PARENT = path.dirname(CWD)

    RUN = path.join(PARENT, 'SynthesisRuns', runname)

    print(RUN)

    def dir_structure(index):
        SYSTEM = path.join(RUN,'system_' + str(index).zfill(4))

        INPUT = path.join(SYSTEM, 'inputs')
        OUTPUT = path.join(SYSTEM, 'outputs')
        makedirs(INPUT)
        makedirs(OUTPUT)

    for i in range(N_sims):
        dir_structure(i)
        #Create DiskFile

        #Create Options File

    return RUN

if __name__ == "__main__":
    Name = sys.argv[1]
    N = sys.argv[2]
    Path = setup(Name, N)






