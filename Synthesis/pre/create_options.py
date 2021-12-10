from random import randint
import os.path as path
import pathlib


def write_option_file(PATH, RUNTIME, EVOTIME):

    options = {}
    # stream = pkg_resources.resource_stream(__name__, 'options_template.txt')
    with open(str(pathlib.Path(__file__).resolve().parent)+"/options_template.txt", 'r') as stream:
        for line in stream:
            (key, val) = line.split()
            options[key] = val

    # Declaring Options
    #   Units: Mass [Solar Mass], Distance [Solar Radius], Time [year]

    # Maximum Runtime
    options["MaxRunTime"] = 0.99 * RUNTIME

    # Maximum Evolution Time
    options["MaxTime"] = EVOTIME

    # Save Interval
    options["SaveInterval"] = EVOTIME // 10

    # Number of Embryos
    options["NEmbryos"] = 0

    # Number of Planetesiamls
    NPlanetesimals = randint(300,500)
    options["NPlanetesimals"] = NPlanetesimals

    # Initial Mass of Embryos
    options["EmbryoInitMass"] = 3.0e-07

    # Density of Embryos
    options["EmbryoRho"] = 0.93357

    # Spacing between Embryos in Hill Radii
    options['Spacing'] = 10

    # Initial Mass of Planetesimals
    options["InitMass"] = 1.6e-11

    # Density of Planetesimals
    options["Rho"] = 0.93357

    # Disperion Timescale
    options["DispersionTime"] = 500000.0

    # Colling Timescale
    options["CoolingTime"] = 100000.0

    # PebbleFlux
    options["PebbleFlux"] = 1e-9

    # Stokes Number of the Pebbles
    options["StokesNumber"] = 0.1

    #write dictionary to txt file
    max_len = max(len(l) for l in options.keys())
    with open(path.join(PATH,"options.txt"), 'w+') as f:
        for key, value in options.items():
            f.write('{}{}\n'.format(key.ljust(max_len+3), value))

if __name__ == "__main__":
    write('here')
