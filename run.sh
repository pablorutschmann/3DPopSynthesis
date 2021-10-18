name=$1

# move disk to inputs

cp disks/disk_log_1000.txt system_0000/inputs/disk.txt

# Empty outputs folder

rm -rf system_0000/outputs/
mkdir system_0000/outputs

cp -r restart system_0000/outputs/

# Start Simulation

Executable/3DPopSyn system_0000/inputs system_0000/outputs system_0000/history.txt

cp -r system_0000 Runs/$name

python3 Post_Process/testing.py $name
