import Post_Process as post
import sys
from os.path import join
from os import getcwd

if __name__ == "__main__":
    print(sys.argv[1])

    pathname = sys.argv[1]

    CWD = getcwd()

    Runs = join(CWD, 'Runs')

    pathname = join(Runs, pathname)

    current_run = post.run(pathname)

    current_run.Print_Sigma(0)

    current_run.plot_snapshots()

    current_run.plot_disk_evol_all()

    current_run.plot_accretion()


    for i in range(1,len(current_run.satellites),len(current_run.satellites)//10):
        current_run.satellites[i].plot_evo(current_run.plot_path)
