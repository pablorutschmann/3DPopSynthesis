name=$1

# Create Directory Structure

# Remove Run directory if it exists
[[ -d $PWD/Runs/$name ]] && rm -r $PWD/Runs/$name

# Set Path Variables
run=$PWD/Runs/$name
input=$run/inputs
output=$run/outputs

mkdir -p $run

mkdir -p $output
mkdir -p $input

# copy disk to inputs
cp $PWD/disks/disk_log_1000.txt $input/disk.txt

# copy options to inputs
cp $PWD/system_0000/options.txt $input/options.txt

# Start Simulation
echo $input
echo $output

$PWD/3DPopSyn $input $output $PWD/system_0000/history.txt

# Star Post-Processing
echo "Starting Plotting!"

python3 POST_PROCESS.py $name

echo "Done!"
