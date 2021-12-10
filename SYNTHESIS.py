import sys

import Synthesis.pre as pre

if __name__ == "__main__":
    Name = sys.argv[1]
    N = sys.argv[2]
    RUNTIME = sys.argv[3]
    EVOTIME = sys.argv[4]

    # Name = 'Test'
    # N = 5
    cmd = pre.setup(Name, N, RUNTIME, EVOTIME)
    print(cmd)