import sys
import os.path as path
from os import getcwd
from setup import setup

if __name__ == "__main__":
    Name = sys.argv[1]
    N = sys.argv[2]
    Path = setup(Name, N)




# python Synthesis/test.py | python Synthesis/setup.py