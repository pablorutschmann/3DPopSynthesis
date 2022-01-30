# Planet Formation Simulation Code for 3D-Population Synthesis

Based on Cilibrasi et al. [2021]. (doi:  10.1093/mnras/stab1179)

## Preparation

### Directory Structure

Before running a simulation, a directory structure has to be prepared and the appropriate input files have to be placed.
To run a simulation in the folder `system_0000` the following structure must be in place:

* [history.txt](./history.txt)
* [system_0000](./system_0000)
    * [inputs](./inputs)
        * [options.txt](./inputs/options.txt)
        * [disk.txt](./inputs/disk.txt)
    * [outputs](./outputs)

The `history.txt` can be placed anywhere, but the path has to be specified when executing(more detail later). It is used
to keep track of all simulation runs and stores starting, end and runtime. In the `input` directory the `options.txt`
and `disk.txt files have to be saved. More on the structure of these files in the next section. The output directory
should be empty.

### Input Files

#### `options.txt` File

The `options.txt` file is a space seperated list of keys and values. An example file can be found in the repository.
Here we list all options with a short description:

- `MP`: The mass of the central star in Jupiter units.
- `RP`: The radius of the central star in Jupiter units.
- `DustToGas`: The dust-to-gas ratio of the disk model.
- `NEmbryos`: The number of planetary embryos in the simulation.
- `Nplanetesimals`: The number of smaller planetesimals in the simulation.
- `Spacing`: The spacing between embryos in Hill radii.
- `SigmaExponent`: The exponent of the power law used for the disk model (always negative).
  `SigmaNorm`: The normalization constant of the power law used for the disk model, in units of central star.
- `TempExponent`: The temperature exponent of the power law used for the disk model (always negative).
- `R_min`: Minimum spawning orbital distance of embryos and planetesimals.
- `R_max`: Maximum spawning orbital distance of embryos and planetesimals.
- `EmbryoInitMass`: Initial mass of the embryos in central star units.
- `InitMass`: Initial mass of the planetesimals in central star units.
- `Rho`: The density of the solids in central star units.
- `FeedRadius`: The extent of the feeding zone in Hill radii.
- `ThresholdMass`: The mass above which bodies are treated with full close-encounter treatment.
- `MaxTime`: The virtual simulation time in years.
- `SaveInterval`: The virtual simulation time interval between saving snapshots. E.g. for `MaxTime = 100`
  and `SaveInteval = 10`, `100/10 = 10` snapshots will be created.
- `MaxRunTime`: The maximum real running time in hours.
- `RotationFraction`: Maximum fraction of the orbital period a body can move within a timestep.
- `DiskPrecision`: Factor determining the disk evolution timestep by `dt = DiskPrecision * DispersionTime`.
- `DtMax`: Maximum timestep in years.
- `DisperionTime`: The dispersion timescale in years.
- `CoolingTime`: The colling timescale in years.
- `SublimationFactor`: Scaling factor adjusting the sublimation rate.
- `AccretionCoeff`: Scaling factor adjusting the accretion rate.
- `PebbleFlux`: The inward pebble flux in central stars units per year.
- `StokesNumber`: The dimensionless Stokes number of the pebbles.
- `RCavity`: The radius of the cavity around the central star in central star radii.
- `RDistruction`: The minimum orbital radius before a body is accreted and destroyed by the central star, in central
  star radii.
- `MigrationType`: The type of migration used in the code:
    - `0`: Tanaka
    - `1`: Dangelo
    - `2`: Paardekooper
    - `3`: Jimenez and Lin
    - `4`: Jimenez
- `Migration`: Turning the migration on/off (1/0).
- `NBody`: Turning the N-body interactions on/off (1/0).
- `EccDamping`: Turning the eccentricity damping on/off (1/0).
- `IncDamping`: Turning the inclination damping on/off (1/0).`
- `Accretion`: Turning the accretion on/off (1/0).`
- `PebbleAccretion`: Turning the pebble accretion on/off (1/0).`
- `TraceWM`: Turning the tracing of the water mass on/off (1/0).`
- `Sublimation`: Turning the sublimation on/off (1/0).`
- `Cooling`: Turning the disk cooling on/off (1/0).`
- `GasDispersion`: Turning the gas disk dispersion on/off (1/0).`
- `DustDispersion`: Turning the dust disk dispersion on/off (1/0).`
- `PebbleDispersion`: Turning the dispersion of the pebble flux on/off (1/0).`
- `PebbleFiltering`: Turning the pebble filtering by large bodies on/off (1/0).`
- `Refilling`: Turning the refilling of the disk on/off (1/0).`
- `Snapshots`: Turning snapshots on/off (1/0).`
- `SeedFormation`: Turning the formation of seeds on/off (1/0).`

#### `disk.txt` File

The disk file needs to consist of 1000 cells. Each cell corresponds to a row in the disk file. The columns have to be
ordered as follows:

- `index`: index of the cells starting from 0.
- `radius`: Midpoint orbital distance of cell in central star radii.
- `dradius`: Annular width of the cell in central star radii.
- `sigmagas`: Gas density value in central star units.
- `sigmadust`: Dust density value in central star units.
- `sigmadustbar`: Backup dust density value in central star units.
- `temperature`: Temperature value in Kelvin.
- `area`: Area of the annulus in central star units.
- `keplerianvelocity` : Keplerian velocity at `radius` in central star radii per year.
- `gasopacity`: Gas opacity of the disk. This should be filled with zeros, since the code calculates the opacities.
- `wmf`: Water mass fraction bodies inherit when initialized at `radius`.
- `swmf`: Solid water mass fraction bodies inherit when initialized at `radius`.
- `headwindfactor`: Headwind factor. Only used for pebble accretion. Fill with zeros in case of no pebble accretion.

The values within a row should be `space` or `tab` separated. Column names must not be included in the input disk file.

## Installation and Compiling

The code was developed in C++ 11 using CMake version 3.20 for ease of debugging and compilation. In the repository you
can find the `CMakeLists.txt`, and we highly recommend using it. Simple compilation can also be done using a normal
Makefile. Once the repository is pulled into the desired folder. While being in the `3DPopSynthesis` directory, it can
be compiled using the following command.

```bash
make build
```

This calls the file `Makefile` and creates a folder `build_outputs`, where all intermediate build files are stored. The
executable `3DPopSynis created in the `3DPopSynthesis` directory. Now the simulation is ready to be started.

## Running a single Simulation

Once the directories and input files are set up, the compiled programm can be run. It takes 3 arguments:

- `in`: The path to the `inputs` directory.
- `out`: The Path to the `outputs` directory.
- 'hist': The path to the `history.txt` file.

To start the simulation run:

```bash
./3DPopSyn in out hist
```

If you and the example folder `system_0000` are in the same directory as the compiled Program, the command to start the
simulation would look as follows:

```bash
./3DPopSyn system_0000/inputs system_0000/outputs history.txt
```

If the code has been run but crashed or ended before finishing, it can be restarted by entering the same command again.
The code automatically detects the last state and continues from there.

### Output Files

Once the simulation is completed you can see several files and directories in the `output` folder. The three files
keeping track of the satellites and their collisions are: `satellite_list.txt`, `collisions.txt`
and `lost_satellites.txt`. For every snapshot a directory with the files `satellites.txt`, `disk.txt`
and `parameters.txt` is created. Finally, a folder called `restart` is used as backup in the case the simulation
crashes. `restart` looks exactly like a snapshot but without the column names in the files. All output files consist of
rows for each entry. The values per row are `space` or `tab` separated.

#### `satellite_list.txt`

- `#ID`: Index of the body
- `init_time`: Time of initialization
- `Type`: Type of object (1: Embryo, 0: Planetesimal)
- `mass`: Mass
- `wm`: Water mass
- `swm`: Solid water mass
- `r2d`: 2-dimensional distance to central star (sqrt(x^2 + y^2))
- `theta`: Azimuthal angle
- `x`,`y`,`z`: x,y,z cartesian coordinates
- `init_temp`: Local Temperature at creation

#### `lost_satellites.txt`

- `#ID`: Index of the body
- `time`: Time of being lost
- `type`: Type of object (1: Embryo, 0: Planetesimal)
- `mass`: Mass
- `wm`: Water mass
- `swm`: Solid water mass
- `r`: 3-dimensional distance to central star
- `x`,`y`,`z`: x,y,z cartesian coordinates
- `xv`,`yv`,`zv`: x,y,z cartesian components of the velocity
- `a`: Semi-major axis
- `ecc`: Eccentricity
- `inc`: Inclination
- `formation_time`: Time for it to reach size of interest
- `collision`: reason of destruction:
    - `0`: engulfed by central star
    - `1`: collision with other body
    - `2`: ejected from system (eccentricity >= 1)
- `collision_index`: Index of the collision, otherwise `-1`

#### `collision_list.txt`

- `time`: Time of the collision
- `ID1`: Index of the body 1
- `type1`: Type of body 1 (1: Embryo, 0: Planetesimal)
- `mass1`: Mass of body 1
- `wm1`: Water mass of body 1
- `swm1`: Solid water mass
- `x1`,`y1`,`z1`: x,y,z cartesian coordinates of body 1
- `xv1`,`yv1`,`zv1`: x,y,z cartesian components of the velocity of body 1
- `ID2`: Index of the body 2
- `type2`: Type of body 2 (1: Embryo, 0: Planetesimal)
- `mass2`: Mass of body 2
- `wm2`: Water mass of body 2
- `swm2`: Solid water mass
- `x2`,`y2`,`z2`: x,y,z cartesian coordinates of body 2
- `xv2`,`yv2`,`zv2`: x,y,z cartesian components of the velocity of body 2
-
    - `ID`: Index of the resulting body
- `type`: Type of the resulting body  (1: Embryo, 0: Planetesimal)
- `mass`: Mass of the resulting body
- `wm`: Water mass of the resulting body
- `swm`: Solid water mass of the resulting body
- `x`,`y`,`z`: x,y,z cartesian coordinates of the resulting body
- `xv`,`yv`,`zv`: x,y,z cartesian components of the velocity of the resulting body

#### `Snapshot_*/satellites.txt`

- `#ID`: Index of the body
- `Type`: Type of object (1: Embryo, 0: Planetesimal)
- `M`: Mass
- `WM`: Water mass
- `SWM`: Solid water mass
- `x`,`y`,`z`: x,y,z cartesian coordinates
- `xv`,`yv`,`zv`: x,y,z cartesian components of the velocity
- `a`: Semi-major axis
- `e`: Eccentricity
- `i`: Inclination
- `N`: Current timestep index
- `dt`: Current timestep
- `init_time`: Time of initialization
- `formation_time`: Time for it to reach size of interest
- `P`: Gap opening parameter

#### `Snapshot_*/disk.txt`

The same structure as the input disk file.

#### `Snapshot_*/parameters.txt`

- `Time`: Time within the simulation
- `UpdateTime`: Last time at which the disk was updated
- `SaveIndex`: Index of the snapshot
- `TimeStopFormation`: Time when new object are not created anymore
- `GlobalDt`: Global timestep which is multiplied with `N` from `satellites.txt` to get the indivual timestep.
- `IceLineId`: Index of the Iceline within the disk file.
- `IceLineRadius`: Orbital distance of the Iceline
- `PebbleFlux`: Current pebble flux

## Code Overview

The code is structured into modules which define different parts and handle specific tasks. The file `main.cpp` serves
as an interface to the  `EvolutionModel` class defined in the `EvolutionModule.*` files. It sets up a model, starts it
and counts the runtime. In the `EvolutionModule.cpp` all other modules are combined to calculate all relevant
quantities. A `DiskModule::DiskModel` is set up containing all relevant information for the evolution of the disk. It
includes functions to calculate disk properties and evolve the disk in time. A list of `SatelliteModule::SatelliteModel`
is created in for teh `EvolutionModel`. Again, the `SatelliteModel` contains all relevant parameters an individual
satellite. The dynamical evolution is however handled in with the `NBodyModule`. With pointers pointing on the
coordinates and velocities it changes theses with time. When ever the satellites interact with the disk, the calcuations
are done within the `EvolutionModule`.







