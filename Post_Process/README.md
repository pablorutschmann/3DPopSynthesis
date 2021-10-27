# 3DPopSynthesis Post-Process Python Package

Rutschmann Pablo, September 2021

rupablo@student.ethz.ch

This package allows for simple and automated post-processing of 3DPopSynthesis Code outputs. The package can import the output data from the directory structure and save it in Python data structures. Once imported, the package offers several functionality to generate plots and figures automatically.


##Build status

The current build 

## Installation

### Requirements
* Python 3.8.8 or later
* numpy
* pandas
* matplotlib
* scipy
* astropy

The package was developed in Python 3.8.8. Download the folder and move it to your working directory or Package Library.
You can import the package into your python scripts using the following command.
```python
import Post_Process as post
```
## Features

### `run` Class
#### Initialization
For every run of the Simulation Code 3DPopSynthesis, you should create a class instance. This can be done by calling:
```python
my_run  = post.run(path)
```
where path is the location of the `system_0000` directory. This import all the data from the run. It imports the `collisions.txt`, `lost_satellites.txt` and `satellite_list.txt` and saves them as Pandas Dataframes. Then it creates a dictionary of the snapshots, where the values are instances of the Class `snapshot` (more below). Finally, the data is transformed to have satellite specific information. Again a dictionary of `satellite` Class instances is created.

#### Plotting
* Plotting Snapshots:
```python
my_run.plot_snapshots(ids=None)
```

This generates the `snapshot.plot_satellites()` figures for a list of Snapshot ID's. The default is generating for all Snapshots.

* Plotting Disk Evolution
```python
my_run.plot_disk_evol(field, log=True, N=10, keys=None):
```
Calling this function creates figure showing the Time Evolution of the `field` Disk Profile. `N` is the number of snapshots to use for the time evolution representation. With the `keys` argument you can input a list of keys form the Snapshot dictionary to plot.

```python
my_run.plot_disk_evol_all(N=10):
```
This funtion automaticall generates the figure for the most relevant fields: Gas Surface Density `SigmaGas`, Dust Surface Density `SigmaDust` and Temperature `Temp`.

* Plotting Accretion
```python
my_run.plot_accretion()
```
To generate a figure show the accretion and collisions of the satellites this function can be called. It uses the `satelite` Class to get satellite specific data.

### `snapshot` Class

If you choose to only look at specific Snapshots, you can access the attributes and methods of the `snapshot` Class. 

#### Initialization
A `snapshot` Class instance can be initialized by calling:
```python
my_snapshot = post.snapshot(path)
```
where path is the directory of the Snapshot directory. This imports `satelllites.txt` and `disk.txt` into Pandas Dataframes while saving the parameters from the `parameters.txt` file to Class variables.

#### Plotting
With:
```python
my_snapshot.plot_satellites()
```
you can generate a figure which plots the Semi Major Axis and Eccentricities of the current Satellites. The size of the points correspond to the mass, while the color indicates the Water Mass Fraction.

### `Satellite` Subclass
To initialize a `Satellite` Class instance you need to pass the `ID` of the Satellite and a `run` Class instance, like follows:
```python
my_satellite = post.satellite(ID,my_run)
```
This transform the unspecific data from `my_run` to only `ID` satellite relevant information.

## License
MIT © Pablo Rutschmann, 01.09.2021