module load gcc/4.8.2 python/3.6.0

cd 3DPopSynthesis/CodeBase

bsub -J compile -W 00:01 -oo ../job_history.txt "make CPPFLAGS=-o3"