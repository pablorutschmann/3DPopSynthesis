import sys
from os.path import join
from os.path import isdir
import os
from os import getcwd
import shutil

runs_name = 'run'
N_cores = 48
N_snapshots = 100
final_snapshot = f'Snapshot_0{N_snapshots}'
CWD = getcwd()

pop_path = join(CWD,"SynthesisRuns")

comb = join(pop_path,'combined')
os.mkdir(comb)
global_counter = 1

for run in [ run for run in sorted(os.listdir(pop_path)) if runs_name in run]:
    sys_path = join(pop_path,run)
    n_run = int(run.replace(runs_name,''))
    print(sys_path)
    for i in range(1,N_cores):
        check_path = join(sys_path,f'system_{i}/outputs',final_snapshot)
        print(f'{check_path=}')
        if isdir(check_path):
            original = join(sys_path,f'system_{i}')
            new_system_n = global_counter
            target = join(comb,f'system_{new_system_n}')
            print(f'{original=}')
            print(f'{target=}')
            shutil.copytree(original,target)
            global_counter += 1
        else:
            continue






