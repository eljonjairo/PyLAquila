#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#
# Save each simulation data with its src fields and DG synthetic wave fields
# John Diaz july 2024

from Classes.LAquilaData import Data
from os import scandir
from tqdm import tqdm
import time


if __name__ == '__main__':

    src_folder = "Outputs/Fields/"  # Folder with fields for each simulation
    data_folder = "Outputs/SyntheticVelo_1.0Hz/"  # Folder with DG synthetic velo
    output_folder = "Outputs/MLData_1.0Hz/"  # Output folder for simulation objects
    low_freq = 1.0  # low pass freq
    high_freq = 0.01  # high pass freq

    # Set the stations dict with name and trim time in secs
    # stations = {'AQK': 24, 'AQU': 24, 'AQV': 24, 'AQA': 24, 'AQG': 24,
    #             'GSA': 25, 'MTR': 26, 'FMG': 28, 'ANT': 28, 'AVZ': 29,
    #             'CSO1': 30, 'LSS': 30, 'SUL': 31}

    stations = {'AQK': 0, 'AQU': 0, 'AQV': 0, 'AQA': 0, 'AQG': 0,
                'GSA': 0, 'MTR': 0, 'FMG': 0, 'ANT': 0, 'AVZ': 0,
                'CSO1': 0, 'LSS': 0, 'SUL': 0}

    # Scan src_folder for files with src distributions for each available simulation
    sims = [f for f in scandir(src_folder) if f.is_file()]
    n_sim = len(sims)

    ini = 0
    for sim in tqdm(sims[ini:]):

        name = sim.name.split('.')[0]
        print()
        print(f" Processing simulation: {name}")
        data = Data(name)  # Create a Data object with the name of the simulation
        # low_freq and high_freq are optional
        # data.set_dg_data(stations, data_folder)
        data.set_dg_data(stations, data_folder, low_freq=low_freq)
        # data.set_dg_data(stations, data_folder, high_freq=high_freq)
        # data.set_dg_data(stations, data_folder, low_freq=low_freq,
        # high_freq=high_freq)
        data.extract_src_fields(sim.path)
        data.plot_simulation()
        print()
        # save = input(" Save the data object (q) ").strip()
        data.save_simulation_data(output_folder)

        # time.sleep(1)



