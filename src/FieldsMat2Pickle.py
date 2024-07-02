#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#
# Save each test fields (slip, rise time, rupture time) into pickle file
#
# John Diaz july 2024

from scipy.io import loadmat
from os import scandir
import pickle


class Source:
     def __init__(self, name):
        self.name = name
        self.SLIP = None  # Slip distribution
        self.RISE = None  # Rise Time distribution
        self.RUPT_T = None  # Rupture Time distribution
        self.RUPT_V = None  # Rupture Velocity distribution

if __name__ == '__main__':
    # Input folder
    input_folder = "/home/jon/PhD/codes/LAquila_2023/LAquila/Source/SourceGen/srate_files/Sources_mat/"

    # Output folder
    output_folder = "Outputs/Fields/"

    files = [f.path for f in scandir(input_folder) if f.is_file()]

    for file in files:
        f_name = file.split('/')[-1]
        f_model = "LAquila_ID" + f_name.split('_')[-1].split('.')[0]
        source = Source(f_model)
        print(f" Loading file: {f_name}")
        input_source = loadmat(file)
        source.SLIP = input_source['SourceParam']['slip'][-1][0]
        source.RISE = input_source['SourceParam']['riseT'][-1][0]
        source.RUPT_T = input_source['SourceParam']['rupTime'][-1][0]
        source.RUPT_V = input_source['SourceParam']['Vr'][-1][0]

        out_file = output_folder + source.name + ".pickle"
        object_file = open(out_file, 'wb')
        pickle.dump(source, object_file)
        object_file.close()
        print(f" Fault object saved in: {out_file} ")
