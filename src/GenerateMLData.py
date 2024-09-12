#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#
# Generate pickle files with dg and src data fir each simulation
#
# John Diaz December 2022


import numpy as np
from os import scandir
from pathlib import Path
from Classes.LAquilaData import Data
from tqdm import tqdm

if __name__ == '__main__':
    # DG folder
    DGFolder = "../DGrun_1.0Hz_alpha_beta/"

    # Folder with fields for each simulation
    src_folder = "/home/jon/PhD/codes/LAquila_2023/LAquila/3DFaults/3DFaults_mat/LAquilaCirella2012bdhF500m"

    # Output Folder
    output_folder = "Outputs/MLData_alpha_beta_1.0Hz/"  # Output folder for simulation objects

    # Stations name in the same order as GEODG3D.acquisition DG file.
    stations = ["AQK", "AQU", "AQV", "AQA", "AQG", "GSA", "MTR", "FMG", "ANT",
                "AVZ", "CSO1", "LSS", "SUL"]

    # 1.0 Hz simulations
    nt = 816
    dt = 0.0490204

    # Cutoff freqs
    freqmin = 0.02  # Cirella
    freqmax = 1.0

    print("  ")
    print(" START PROGRAM ")
    print("  ")

    nstats = len(stations)
    ini = 0

    # Start looking for the simulations folders
    subfolders = [f.path for f in scandir(DGFolder) if f.is_dir()]

    for sub_f in tqdm(subfolders[ini:1]):
        file = Path(sub_f)
        name = sub_f.split('/')[2]
        file_vx = file / "VX_1"
        if file_vx.exists():
            file_vy = file / "VY_1"
            file_vz = file / "VZ_1"
            print("Loading DG velocity synthetic files:")
            print(f"  {file_vx}")
            print(f"  {file_vy}")
            print(f"  {file_vz}")

            # Create a Data object with the name of the simulation
            data = Data(name, nstats, nt, dt, freqmin, freqmax)

            # Download DG velo files and convert from m/s to in cm/s
            velo_syn_x = np.reshape(np.fromfile(file_vx, dtype=np.float32),
                                    (nstats, nt), order='F')*100
            velo_syn_y = np.reshape(np.fromfile(file_vy, dtype=np.float32),
                                    (nstats, nt), order='F')*100
            velo_syn_z = np.reshape(np.fromfile(file_vz, dtype=np.float32),
                                    (nstats, nt), order='F')*100

            # Set velo traces
            data.set_velo_acc(velo_syn_x, velo_syn_y, velo_syn_z, stations)
            data.set_velo_spec()
            data.set_arias()
            data.set_src(src_folder)
            data.plot_data()
            data.save_simulation_data(output_folder)

        else:
            print(f" file {file_vx} not found")
