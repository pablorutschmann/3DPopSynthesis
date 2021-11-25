rsync -Pav /Users/prut/CLionProjects/3DPopSynthesis/CodeBase/  rupablo@euler.ethz.ch:/cluster/home/rupablo/3DPopSynthesis/CodeBase/

rsync -Pav /Users/prut/CLionProjects/3DPopSynthesis/disks/  rupablo@euler.ethz.ch:/cluster/home/rupablo/3DPopSynthesis/disks/

rsync -Pav /Users/prut/CLionProjects/3DPopSynthesis/Post_Process/  rupablo@euler.ethz.ch:/cluster/home/rupablo/3DPopSynthesis/Post_Process/

rsync -Pav /Users/prut/CLionProjects/3DPopSynthesis/POST_PROCESS.py  rupablo@euler.ethz.ch:/cluster/home/rupablo/3DPopSynthesis/POST_PROCESS.py

rsync -Pav /Users/prut/CLionProjects/3DPopSynthesis/run.sh  rupablo@euler.ethz.ch:/cluster/home/rupablo/3DPopSynthesis/run.sh

rsync -Pav /Users/prut/CLionProjects/3DPopSynthesis/restart.sh  rupablo@euler.ethz.ch:/cluster/home/rupablo/3DPopSynthesis/restart.sh



ssh euler

cd 3DPopSynthesis/CodeBase

bsub -o job_history.txt make


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