import sys

from Synthesis.pre.setup import setup

if __name__ == "__main__":
    Name = sys.argv[1]
    N = sys.argv[2]
    RUNTIME = sys.argv[3]
    EVOTIME = sys.argv[4]

    # Name = 'Test'
    # N = 5
    cmd = setup(Name, N, RUNTIME, EVOTIME)
    print(cmd)