import snapshots
import sys
print(sys.argv[1])

path = sys.argv[1]

test = snapshots.sat_data('/Users/prut/CLionProjects/3DPopSynthesis/Runs/' + path + '/')

test.plot_remaining()

#test.accretion()