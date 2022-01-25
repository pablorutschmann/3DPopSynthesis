rsync -Pav /Users/prut/CLionProjects/3DPopSynthesis/CodeBase/  rupablo@euler.ethz.ch:/cluster/home/rupablo/3DPopSynthesis/CodeBase/

rsync -Pav /Users/prut/CLionProjects/3DPopSynthesis/disks/  rupablo@euler.ethz.ch:/cluster/home/rupablo/3DPopSynthesis/disks/

rsync -Pav /Users/prut/CLionProjects/3DPopSynthesis/Post_Process/  rupablo@euler.ethz.ch:/cluster/home/rupablo/3DPopSynthesis/Post_Process/

rsync -Pav /Users/prut/CLionProjects/3DPopSynthesis/POST_PROCESS.py  rupablo@euler.ethz.ch:/cluster/home/rupablo/3DPopSynthesis/POST_PROCESS.py

rsync -Pav /Users/prut/CLionProjects/3DPopSynthesis/run.sh  rupablo@euler.ethz.ch:/cluster/home/rupablo/3DPopSynthesis/run.sh

rsync -Pav /Users/prut/CLionProjects/3DPopSynthesis/restart.sh  rupablo@euler.ethz.ch:/cluster/home/rupablo/3DPopSynthesis/restart.sh



ssh euler

cd 3DPopSynthesis/CodeBase

bsub -o ../job_history.txt "make CPPFLAGS=-O0"

rsync -Pav /Users/prut/CLionProjects/3DPopSynthesis/system_0000/*  rupablo@euler.ethz.ch:/cluster/home/rupablo/3DPopSynthesis/

rm -r $HOME/3DPopSynthesis/Runs/test/outputs/*

rsync -Pav /Users/prut/CLionProjects/3DPopSynthesis/system_0000/options.txt  rupablo@euler.ethz.ch:/cluster/home/rupablo/3DPopSynthesis/Runs/test/inputs/

bsub -oo $HOME/3DPopSynthesis/Runs/test/log CodeBase/3DPopSyn $HOME/3DPopSynthesis/Runs/test/inputs $HOME/3DPopSynthesis/Runs/test/outputs /home/rupablo/3DPopSynthesis/history.txt

bsub -W 23:50 -oo $HOME/3DPopSynthesis/Runs/test/log $HOME/3DPopSynthesis/CodeBase/3DPopSyn $HOME/3DPopSynthesis/Runs/long/inputs $HOME/3DPopSynthesis/Runs/long/outputs /home/rupablo/3DPopSynthesis/history.txt

#Restart
bsub -W 23:59 -o $HOME/3DPopSynthesis/Runs/long/log $HOME/3DPopSynthesis/CodeBase/3DPopSyn $HOME/3DPopSynthesis/Runs/long/inputs $HOME/3DPopSynthesis/Runs/long/outputs /home/rupablo/3DPopSynthesis/history.txt

bsub -J job1 -W 23:59 -o $HOME/3DPopSynthesis/Runs/long/log $HOME/3DPopSynthesis/CodeBase/3DPopSyn $HOME/3DPopSynthesis/Runs/long/inputs $HOME/3DPopSynthesis/Runs/long/outputs /home/rupablo/3DPopSynthesis/history.txt
bsub -J job2 -w "done(job1)" -W 23:59 -o $HOME/3DPopSynthesis/Runs/long/log $HOME/3DPopSynthesis/CodeBase/3DPopSyn $HOME/3DPopSynthesis/Runs/long/inputs $HOME/3DPopSynthesis/Runs/long/outputs /home/rupablo/3DPopSynthesis/history.txt
bsub -J job3 -w "done(job2)" -W 23:59 -o $HOME/3DPopSynthesis/Runs/long/log $HOME/3DPopSynthesis/CodeBase/3DPopSyn $HOME/3DPopSynthesis/Runs/long/inputs $HOME/3DPopSynthesis/Runs/long/outputs /home/rupablo/3DPopSynthesis/history.txt

rm -r $PWD/Runs/debug_test/outputs ;mkdir $PWD/Runs/debug_test/outputs

cp disks/disk_log_1000.txt Runs/debug_test/inputs/disk.txt

rsync -Pav rupablo@euler.ethz.ch:/cluster/home/rupablo/3DPopSynthesis/SynthesisRuns/initsynth /Users/prut/CLionProjects/3DPopSynthesis/SynthesisRuns/

python SYNTHESIS.py thresh 500 120 1e6

rsync -Pav rupablo@euler.ethz.ch:/cluster/home/rupablo/3DPopSynthesis/RUNTIME /Users/prut/CLionProjects/3DPopSynthesis/RUNTIME

rsync -Pav rupablo@euler.ethz.ch:/cluster/home/rupablo/3DPopSynthesis/Runs /Users/prut/CLionProjects/3DPopSynthesis/MigrationRuns

python SYNTHESIS.py final2 48 120 1000000
/cluster/home/rupablo/3DPopSynthesis


bsub -J "final2[1-48]" -n 1 -r -W 120:00 -oo /cluster/home/rupablo/3DPopSynthesis/SynthesisRuns/final2/log "/cluster/home/rupablo/3DPopSynthesis/CodeBase/3DPopSyn /cluster/home/rupablo/3DPopSynthesis/SynthesisRuns/final2/system_\$LSB_JOBINDEX/inputs /cluster/home/rupablo/3DPopSynthesis/SynthesisRuns/final2/system_\$LSB_JOBINDEX/outputs /cluster/home/rupablo/3DPopSynthesis/SynthesisRuns/final2/history.txt"