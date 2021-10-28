name=$1

# Create Directory Sructure

# Remove Run directory if it exists
[[ -d Runs/$name ]] && rm -r Runs/$name

#[[ -d Runs/testnew ]] && rm -r Runs/testnew
#mkdir -p Runs/testnew/
#
#mkdir -p Runs/testnew/outputs/
#mkdir -p Runs/testnew/inputs/



mkdir -p Runs/$name/

mkdir -p Runs/$name/outputs/
mkdir -p Runs/$name/inputs/

# copy disk to inputs
cp disks/disk_log_1000.txt Runs/$name/inputs/disk.txt

# copy options to inputs
cp system_0000/options.txt Runs/$name/inputs/options.txt

# Start Simulation
Executable/3DPopSyn Runs/$name/inputs Runs/$name/outputs system_0000/history.txt

# Star Post-Processing
echo "Starting Plotting!"

python3 POST_PROCESS.py $name

echo "Done!"
