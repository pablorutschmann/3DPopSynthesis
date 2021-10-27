name=$1

# move disk to inputs

cp disks/disk_log_1000.txt system_0000/inputs/disk.txt

# Empty outputs folder

rm -rf system_0000/outputs/
mkdir system_0000/outputs

# cp -r restart system_0000/outputs/

# Start Simulation

Executable/3DPopSyn system_0000/inputs system_0000/outputs history.txt

[[ -d Runs/$name ]] && rm -r Runs/$name

cp -r system_0000 Runs/$name

echo "Starting Plotting!"

python3 POST_PROCESS.py $name

echo "Done!"
