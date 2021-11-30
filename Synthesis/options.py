def write(path):

    options = {}
    with open("options_template.txt") as f:
        for line in f:
            (key, val) = line.split()
            options[key] = val
    print(options)

    # Number of Embryos
    NEmbryos = np.random.uniform(low = 7, high = 14)
    options["NEmbryos"] = NEmbryos

    # Number of Planetesiamls
    NPlanetesimals = 100
    options["NPlanetesimals"] = NPlanetesimals

    # Initial Mass of Embryos
    EmbryoInitMass =




    with open('options.txt', 'w') as f:
        f.write("{0}\t{1}\n".format("NEmbryos", NEmbryos)








if __name__ == "__main__":
    write('here')
