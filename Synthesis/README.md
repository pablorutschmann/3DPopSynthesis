# Synthesis Framework Module

This module allows for automation of population synthesis with the simulation code stored in `CodeBase`. It consists of
two sub-modules, one which handles the preprocessing and prepares the appropriate directories and files and one that
handles the data analysis of the simulation data.

## Requirenments

- Python3
- numpy
- pandas
- astropy

## Sub-Modules

### `pre`

This module handles the generation of the input directories of the entire population. It creates the following directory
structure:

* [Synthesis_Runs](./Synthesis_Runs)
    * [name](./Synthesis_Runs/name)
        * [history.txt](./Synthesis_Runs/name/history.txt)
        * [system_1](./Synthesis_Runs/name/system_1)
            * [inputs](./Synthesis_Runs/name/system_1/inputs)
                * [options.txt](./Synthesis_Runs/name/system_1//inputs/options.txt)
                * [disk.txt](./Synthesis_Runs/name/system_1//inputs/disk.txt)
            * [outputs](./Synthesis_Runs/name/system_1/outputs)
        * [system_2](./Synthesis_Runs/name/system_2)
        * [system_3](./Synthesis_Runs/name/system_3)
        * [system_4](./Synthesis_Runs/name/system_4)
        * ...

Options that are the same for all simulation and are not specifiable at execution, can be set in the
file `Synthesis/pre/options_template.txt`. In the file `Synthesis/functions.py` you change the functions used for the
disk profiles, like the power-law, temperature profile and water distribution. To change the population synthesis
parameters, you have to go into the file `setup.py` and change the relevant values.

#### Running pre

The python script `SYNTHESIS.py` acts as an interface for the preprocessing module. It takes 4 arguments:

- `Name`: The name of the synthesis run
- `N`: The number of simulations in the population
- `RUNTIME`: The real time in hours which should be allocated for an individual simulation
- `Evotime`: Evolution time of the simulations in years

An example command would look as follows:

```bash
python SYNTHESIS.py name 100 120 1000
```

This creates 100 simulations in the folder `name`, which are limited to `120` hours runtime and evolve up to `1000`
years. For now this only works for the ETH euler cluster. It submits the simulations as a job
array ([Euler: Job Array](https://scicomp.ethz.ch/wiki/Job_arrays)). This way, 48 simulations can be run in parallel.
For more simulations, the other are started as soon as another finishes, such that there are always 48 simulations
running. Since the euler cluster is restricted to only 120 hours per job, heavy simulations would exceed this, are
stopped and need to be restarted. The easiest way to run more than 48 simulations is to only submit 48 simulations at a
time, but do this severall times. This allows for more controlled restarting by just restarting all jobs. The ones
that have already finished, finish immediately and the next in queue is started while the other continue the simulation.

### `post`

in-progress

