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
can find the `CMakeLists.txt`, and we highly recommend using it. Simple compilation can also be done using the normal
Makefile. Once the repository is pulled into the desired folder it can be compiled using the following command.

```cmake

```


