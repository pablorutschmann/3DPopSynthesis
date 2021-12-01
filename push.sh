#!/bin/zsh

COMP=${1:-0}

rsync -Pav /Users/prut/CLionProjects/3DPopSynthesis/CodeBase/  rupablo@euler.ethz.ch:/cluster/home/rupablo/3DPopSynthesis/CodeBase/

rsync -Pav /Users/prut/CLionProjects/3DPopSynthesis/Synthesis/  rupablo@euler.ethz.ch:/cluster/home/rupablo/3DPopSynthesis/Synthesis/

rsync -Pav /Users/prut/CLionProjects/3DPopSynthesis/disks/  rupablo@euler.ethz.ch:/cluster/home/rupablo/3DPopSynthesis/disks/

rsync -Pav /Users/prut/CLionProjects/3DPopSynthesis/system_0000/options.txt  rupablo@euler.ethz.ch:/cluster/home/rupablo/3DPopSynthesis/options.txt

rsync -Pav /Users/prut/CLionProjects/3DPopSynthesis/Post_Process/  rupablo@euler.ethz.ch:/cluster/home/rupablo/3DPopSynthesis/Post_Process/

rsync -Pav /Users/prut/CLionProjects/3DPopSynthesis/POST_PROCESS.py  rupablo@euler.ethz.ch:/cluster/home/rupablo/3DPopSynthesis/POST_PROCESS.py

rsync -Pav /Users/prut/CLionProjects/3DPopSynthesis/run_euler.sh  rupablo@euler.ethz.ch:/cluster/home/rupablo/3DPopSynthesis/run_euler.sh

rsync -Pav /Users/prut/CLionProjects/3DPopSynthesis/restart.sh  rupablo@euler.ethz.ch:/cluster/home/rupablo/3DPopSynthesis/restart.sh

rsync -Pav /Users/prut/CLionProjects/3DPopSynthesis/compile.sh  rupablo@euler.ethz.ch:/cluster/home/rupablo/compile.sh

if [[ $COMP == 0 ]] ; then
  echo "Not Compiling!"
else
  ssh euler 'bash compile.sh'
fi
