import sys
import os.path as path
from os import getcwd
from os import makedirs
from options import write_option_file

def setup(NAME, N_SIMS):
    runname = NAME

    N_sims = int(N_SIMS)

    CWD = getcwd()

    RUN = path.join(CWD, 'SynthesisRuns', runname)

    print(RUN)

    def dir_structure(index):
        SYSTEM = path.join(RUN,'system_' + str(index).zfill(4))

        INPUT = path.join(SYSTEM, 'inputs/')
        OUTPUT = path.join(SYSTEM, 'outputs/')
        makedirs(INPUT, exist_ok=True)
        makedirs(OUTPUT, exist_ok=True)
        return INPUT, OUTPUT

    for i in range(N_sims):
        INPUT, OUTPUT = dir_structure(i)
        #Create DiskFile

        #Create Options File
        write_option_file(INPUT)

        #Make Command in commands file
        

    return RUN

if __name__ == "__main__":
    # Name = sys.argv[1]
    # N = sys.argv[2]

    Name = 'Test'
    N = 5
    Path = setup(Name, N)






