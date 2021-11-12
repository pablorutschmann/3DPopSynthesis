import Post_Process as post
import sys
import os.path

if __name__ == "__main__":

    print(sys.argv[1])

    pathname = sys.argv[1]

    pathname = os.path.join('/Users/prut/CLionProjects/3DPopSynthesis/Runs/',pathname)

    current_run = post.run(pathname)

    current_run.plot_snapshots()

    current_run.plot_disk_evol_all()

    current_run.plot_accretion()

    current_run.plot_wm()

    #current_run.plot_wm()

