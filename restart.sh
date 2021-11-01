name=$1

# Start Simulation
Executable/3DPopSyn Runs/$name/inputs Runs/$name/outputs system_0000/history.txt

# Star Post-Processing
echo "Starting Plotting!"

python3 POST_PROCESS.py $name

echo "Done!"
