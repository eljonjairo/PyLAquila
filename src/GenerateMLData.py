#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#
# Save each simulation data with its src fields and DG synthetic wave fields
# John Diaz july 2024

from Classes.LAquilaData import Data
from os import scandir
from tqdm import tqdm


if __name__ == '__main__':

    src_folder = "Outputs/Fields/"  # Folder with fields for each simulation
    data_folder = "Outputs/SyntheticVelo_1.0Hz/"  # Folder with DG synthetic velo
    output_folder = "Outputs/MLData_1.0Hz/"  # Output folder for simulation objects

    # Set the stations dict with name and trim time in secs
    # stations = {'AQK': 24, 'AQU': 24, 'AQV': 24, 'AQA': 24, 'AQG': 24,
    #             'GSA': 25, 'MTR': 26, 'FMG': 28, 'ANT': 28, 'AVZ': 29,
    #             'CSO1': 30, 'LSS': 30, 'SUL': 31}

    stations = {'AQK': 0, 'AQU': 0, 'AQV': 0, 'AQA': 0, 'AQG': 0,
                'GSA': 0, 'MTR': 0, 'FMG': 0, 'ANT': 0, 'AVZ': 0,
                'CSO1': 0, 'LSS': 0, 'SUL': 0}
    sims = [f for f in scandir(src_folder) if f.is_file()]
    n_sim = len(sims)

    n_saved = 0
    ini = 0
    for sim in tqdm(sims[ini:]):

        name = sim.name.split('.')[0]
        print()
        print(f" Processing simulation: {name}")
        data = Data(name)
        data.set_dg_data(stations, data_folder)
        data.extract_src_fields(sim.path)
        data.plot_simulation()
        print()
        save = input(" Save the data object (q) ")
        if save == 'q':
            data.save_simulation_data(output_folder)
            n_saved += 1
        else:
            # remove(sim.path)
            print(f" {sim.name} not saved")

    print(f" {n_saved} models.")