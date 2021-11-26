#!/bin/zsh

name=$1

#paths
input=$PWD/Runs/$name/inputs
output=$PWD/Runs/$name/outputs

echo $input
echo $output
# Start Simulation
Executable/3DPopSyn $input $output $PWD/system_0000/history.txt

# Star Post-Processing
echo "Starting Plotting!"

python3 POST_PROCESS.py $name

echo "Done!"
