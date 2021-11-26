#!/bin/bash

name=$1

# Create Directory Structure

# Remove Run directory if it exists
[[ -d $HOME/3DPopSynthesis/Runs/$name ]] && rm -r $HOME/3DPopSynthesis/Runs/$name

# Set Path Variables
run=$HOME/3DPopSynthesis/Runs/$name
input=$run/inputs
output=$run/outputs

mkdir -p $run

mkdir -p $output
mkdir -p $input

# copy disk to inputs
cp $HOME/3DPopSynthesis/disks/disk_log_1000.txt $input/disk.txt

# copy options to inputs
cp $HOME/3DPopSynthesis/options.txt $input/options.txt

# Start Simulation
echo $input
echo $output

bsub -J $name -r -W 72:00 -o $run/log $HOME/3DPopSynthesis/CodeBase/3DPopSyn $input $output $HOME/3DPopSynthesis/history.txt
bsub -J ${name}1 -w "done($name)" -r -W 72:00 -o $run/log $HOME/3DPopSynthesis/CodeBase/3DPopSyn $input $output $HOME/3DPopSynthesis/history.txt

# Star Post-Processing
echo "Starting Plotting!"

bsub -J ${name}_POST -w "done(${name}1)" -W 10:00 "python POST_PROCESS.py $name"

echo "Done!"
